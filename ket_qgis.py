# -*- coding: utf-8 -*-
"""
KET for QGIS - Knickzone Extraction Tool (single-file edition)
==============================================================

DEM ひとつから、遷急区間 (knickzone) をワンステップで抽出します。
ArcGIS + TauDEM 版 KET (hdtopography/KET) の処理を numpy + GDAL だけで
書き直したもので、追加インストールは不要です (どちらも QGIS に同梱)。

Use it in two ways:

1. QGIS: Processing Toolbox > Scripts (python icon) > "Add Script to Toolbox..."
   and choose this file. The tool appears as  Scripts > KET > Knickzone
   Extraction (KET).
2. Command line (Python with numpy + GDAL, e.g. the OSGeo4W shell):
       python ket_qgis.py dem.tif -o out_folder

Processing steps:
    1. Depression filling (priority-flood; replaces TauDEM PitRemove / ArcGIS Fill)
    2. D8 flow direction, flow accumulation, upstream main stems
    3. Stream network = cells whose contributing area >= threshold
    4. Multi-scale stream gradient Gd for d = from, from+by, ..., to
       direction='both' (default): e1 = d/2 upstream along the main stem
           (largest contributing area; if equal, the lower cell), stopping at
           the ridge; e2 = d/2 downstream  -> Gd = (e1 - e2) / distance
           (Hayakawa & Oguchi 2006, 2009)
       direction='down': e3 = d downstream -> Gd = (e0 - e3) / distance
           (KET v1 / TauDEM SlopeAveDown; Zahra et al. 2017)
    5. Rd = -(least-squares slope of Gd against d)   [1/m; >0 = locally steep]
    6. Knickzone = connected stream cells with Rd > threshold; length, height
       (elevation range) and other attributes

References
    Hayakawa, Y.S., Oguchi, T. (2006) Geomorphology 78, 90-106.
    Hayakawa, Y.S., Oguchi, T. (2009) Geomorphology 111, 27-37.
    Zahra, T., Paudel, U., Hayakawa, Y.S., Oguchi, T. (2017) Open Geosciences 9, 73-88.
"""

import os
import csv
import heapq
import math
import time

import numpy as np
from osgeo import gdal, ogr, osr

# NOTE: gdal.UseExceptions() is only switched on for the command line (see _main),
# so that the global GDAL state of a running QGIS session is left untouched.

__version__ = "2.0.0"


# ---------------------------------------------------------------------------
#  Core (numpy + GDAL only; no QGIS needed)
# ---------------------------------------------------------------------------

class KETError(Exception):
    pass


def _noop_log(msg):
    print(msg)


def _noop_progress(pct):
    pass


def _never_cancel():
    return False


def _traditional_axis_order(srs):
    """GDAL >= 3 only; older GDAL (2.x) has no axis-mapping API and needs nothing."""
    try:
        srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    except (AttributeError, TypeError):
        pass


def read_dem(path, log=_noop_log):
    """Open a DEM. If it is in a geographic CRS (lat/lon, e.g. GSI 基盤地図 DEM),
    it is reprojected on the fly to the matching UTM zone so that distances are
    in metres."""
    ds = gdal.Open(path)
    if ds is None:
        raise KETError("Cannot open DEM: %s" % path)
    srs = osr.SpatialReference()
    wkt = ds.GetProjection()
    if not wkt:
        raise KETError("The DEM has no coordinate reference system. "
                       "Assign one (projected, in metres) and try again.")
    srs.ImportFromWkt(wkt)
    _traditional_axis_order(srs)

    if srs.IsGeographic():
        gt = ds.GetGeoTransform()
        cx = gt[0] + gt[1] * ds.RasterXSize / 2.0
        cy = gt[3] + gt[5] * ds.RasterYSize / 2.0
        zone = int((cx + 180.0) // 6.0) + 1
        epsg = (32600 if cy >= 0 else 32700) + zone
        res = abs(gt[5]) * 111320.0  # approx. metres per degree of latitude
        res = float("%.2g" % res)
        log("DEM is in a geographic CRS -> reprojecting to UTM zone %d "
            "(EPSG:%d), cell size %.2f m" % (zone, epsg, res))
        nd = ds.GetRasterBand(1).GetNoDataValue()
        ds = gdal.Warp("", ds, format="MEM", dstSRS="EPSG:%d" % epsg,
                       xRes=res, yRes=res, resampleAlg="bilinear",
                       srcNodata=nd, dstNodata=-9999.0, outputType=gdal.GDT_Float64)
        if ds is None:
            raise KETError("Reprojection to UTM failed. Please reproject the DEM to a "
                           "projected CRS (metres) yourself and try again.")
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(epsg)
        _traditional_axis_order(srs)
    elif srs.GetLinearUnits() <= 0 or abs(srs.GetLinearUnits() - 1.0) > 1e-6:
        log("Warning: CRS linear unit is not metre (factor %.4f). Distances and "
            "areas are reported in CRS units x factor." % srs.GetLinearUnits())

    gt = ds.GetGeoTransform()
    if abs(gt[2]) > 0 or abs(gt[4]) > 0:
        raise KETError("Rotated rasters are not supported. Please warp the DEM first.")
    band = ds.GetRasterBand(1)
    z = band.ReadAsArray().astype(np.float64)
    nd = band.GetNoDataValue()
    valid = np.isfinite(z)
    if nd is not None:
        valid &= z != nd
    valid &= z > -1e4  # catch -9999 / -3.4e38 style no-data without metadata
    return z, valid, gt, srs.ExportToWkt()


def fill_depressions(Zp, Vp, ncp, eps=1e-6, progress=_noop_progress,
                     is_canceled=_never_cancel):
    """Priority-flood + epsilon (Barnes et al., 2014) on a padded, flattened DEM.
    Returns the filled DEM and the indices of the seed (edge) cells."""
    offs = (-ncp - 1, -ncp, -ncp + 1, -1, 1, ncp - 1, ncp, ncp + 1)
    nrp = Vp.size // ncp
    V2 = Vp.reshape(nrp, ncp)
    # seeds: valid cells touching an invalid cell (incl. the padding)
    inv = ~V2
    touch = np.zeros_like(V2)
    touch[1:-1, 1:-1] = (inv[:-2, :-2] | inv[:-2, 1:-1] | inv[:-2, 2:] |
                         inv[1:-1, :-2] | inv[1:-1, 2:] |
                         inv[2:, :-2] | inv[2:, 1:-1] | inv[2:, 2:])
    seeds = np.flatnonzero(touch.ravel() & Vp)

    W = Zp.tolist()
    closed = (~Vp).tolist()
    heap = [(W[i], i) for i in seeds.tolist()]
    for i in seeds.tolist():
        closed[i] = True
    heapq.heapify(heap)
    pop, push = heapq.heappop, heapq.heappush
    total = int(Vp.sum())
    done = 0
    step = max(total // 50, 1)
    while heap:
        e, c = pop(heap)
        done += 1
        if done % step == 0:
            progress(100.0 * done / total)
            if is_canceled():
                raise KETError("Canceled")
        for o in offs:
            n = c + o
            if closed[n]:
                continue
            closed[n] = True
            zn = W[n]
            if zn <= e:
                zn = e + eps
                W[n] = zn
            push(heap, (zn, n))
    return np.array(W, dtype=np.float64), seeds


def d8_receivers(F, Vp, ncp, dx, dy):
    """Steepest-descent (D8) receiver of each cell (-1 = outlet) and the
    horizontal length of that step."""
    diag = math.hypot(dx, dy)
    offs = (-ncp - 1, -ncp, -ncp + 1, -1, 1, ncp - 1, ncp, ncp + 1)
    dists = (diag, dy, diag, dx, dx, diag, dy, diag)
    idx = np.flatnonzero(Vp)
    best = np.zeros(idx.size)
    rec = np.full(idx.size, -1, dtype=np.int64)
    step = np.zeros(idx.size)
    Fc = F[idx]
    for o, d in zip(offs, dists):
        nb = idx + o
        ok = Vp[nb]
        s = np.where(ok, (Fc - np.where(ok, F[nb], 0.0)) / d, -np.inf)
        better = s > best
        best[better] = s[better]
        rec[better] = nb[better]
        step[better] = d
    R = np.full(Vp.size, -1, dtype=np.int64)
    L = np.zeros(Vp.size)
    R[idx] = rec
    L[idx] = step
    return R, L


def flow_accumulation(F, R, Vp):
    """Number of cells draining through each cell (including itself)."""
    idx = np.flatnonzero(Vp)
    order = idx[np.argsort(-F[idx], kind="stable")]  # upstream -> downstream
    acc = Vp.astype(np.float64)
    accl = acc.tolist()
    Rl = R.tolist()
    for i in order.tolist():
        r = Rl[i]
        if r >= 0:
            accl[r] += accl[i]
    return np.array(accl)


def upstream_mainstem(F, R, acc, Vp):
    """Main-stem upstream pointer U (-1 = no donor, i.e. ridge / channel head).
    Among the cells draining into a cell, the one with the largest contributing
    area is chosen; if the areas are equal, the lower one (gentler gradient)."""
    idx = np.flatnonzero(Vp & (R >= 0))
    rec = R[idx]
    order = np.lexsort((F[idx], -acc[idx], rec))  # primary: receiver
    rec_s = rec[order]
    first = np.ones(rec_s.size, dtype=bool)
    first[1:] = rec_s[1:] != rec_s[:-1]
    U = np.full(Vp.size, -1, dtype=np.int64)
    U[rec_s[first]] = idx[order][first]
    return U


def flow_lengths(F, R, U, L, Vp):
    """Flow distance to the outlet (along D8) and to the head of the main stem
    (along U, i.e. up to the ridge) for every cell."""
    idx = np.flatnonzero(Vp)
    order = idx[np.argsort(F[idx], kind="stable")]  # low -> high
    Rl, Ul, Ll = R.tolist(), U.tolist(), L.tolist()
    dout = [0.0] * Vp.size
    for i in order.tolist():  # receivers are always lower
        r = Rl[i]
        if r >= 0:
            dout[i] = dout[r] + Ll[i]
    dhead = [0.0] * Vp.size
    for i in reversed(order.tolist()):  # main-stem donors are always higher
        u = Ul[i]
        if u >= 0:
            dhead[i] = dhead[u] + Ll[u]
    return np.array(dout), np.array(dhead)


def walk_endpoints(F, P, Ls, cells, targets, is_canceled=_never_cancel):
    """Follow pointer P (R = downstream, U = upstream main stem) from each start
    cell. For every target distance t (ascending), stop at the first cell where
    the travelled distance reaches t, and return its elevation E and the
    travelled distance D. If the path ends first (outlet / DEM edge downstream,
    ridge upstream), the last cell and the distance travelled so far are used.
    Ls[c] = horizontal length of the step c -> P[c].
    Returns E, D with shape (len(cells), len(targets))."""
    targets = np.asarray(targets, dtype=np.float64)
    nt = targets.size
    E = np.empty((cells.size, nt))
    D = np.empty((cells.size, nt))
    row = np.arange(cells.size)
    cur = cells.copy()
    dist = np.zeros(cells.size)
    kdone = np.zeros(cells.size, dtype=np.int64)

    def record(sel, upto):
        r_s, k0, e_s, d_s = row[sel], kdone[sel], F[cur[sel]], dist[sel]
        k1 = upto
        for j in range(int((k1 - k0).max())):
            k = k0 + j
            m = k < k1
            E[r_s[m], k[m]] = e_s[m]
            D[r_s[m], k[m]] = d_s[m]

    while row.size:
        if is_canceled():
            raise KETError("Canceled")
        nxt = P[cur]
        end = nxt < 0
        if end.any():
            record(end, np.full(int(end.sum()), nt))
            keep = ~end
            row, cur, dist, kdone, nxt = row[keep], cur[keep], dist[keep], kdone[keep], nxt[keep]
            if row.size == 0:
                break
        dist = dist + Ls[cur]
        cur = nxt
        knew = np.searchsorted(targets, dist, side="right")
        crossed = knew > kdone
        if crossed.any():
            record(crossed, knew[crossed])
            kdone = np.maximum(kdone, knew)
        keep = kdone < nt
        row, cur, dist, kdone = row[keep], cur[keep], dist[keep], kdone[keep]
    return E, D


def multiscale_gradient(F, R, U, L, Lup, cells, dlist, direction="both",
                        is_canceled=_never_cancel):
    """Gd for each start cell and each d.
    direction='both' : Gd = (e1 - e2) / (distance e1..e2), e1 = d/2 upstream on
                       the main stem, e2 = d/2 downstream  (SAKE; Hayakawa &
                       Oguchi 2006, 2009)
    direction='down' : Gd = (e0 - e3) / (distance), e3 = d downstream
                       (KET v1 / TauDEM SlopeAveDown; Zahra et al. 2017)"""
    dlist = np.asarray(dlist, dtype=np.float64)
    if direction == "down":
        Ed, Dd = walk_endpoints(F, R, L, cells, dlist, is_canceled)
        num = F[cells][:, None] - Ed
        den = Dd
    elif direction == "both":
        half = dlist / 2.0
        Ed, Dd = walk_endpoints(F, R, L, cells, half, is_canceled)
        Eu, Du = walk_endpoints(F, U, Lup, cells, half, is_canceled)
        num = Eu - Ed
        den = Du + Dd
    else:
        raise KETError("direction must be 'both' or 'down'")
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)


def relative_steepness(G, by):
    """Rd = -a, where a is the least-squares slope of Gd against d
    (Hayakawa & Oguchi 2006, 2009; Zahra et al. 2017, Eq. 4). Positive Rd =
    gradient decreases with d = locally steep. For equally spaced d:
    a = (12*sum(n*G_n) - 6*(N+1)*sum(G_n)) / (by * N*(N+1)*(N-1))
    NOTE: KET v1 (ArcGIS) saved +a as 'Finalraster', i.e. v1 Finalraster = -Rd."""
    N = G.shape[1]
    if N < 2:
        raise KETError("At least two d values are needed (check from/to/by).")
    n = np.arange(1, N + 1, dtype=np.float64)
    a = (12.0 * (G * n).sum(axis=1) - 6.0 * (N + 1) * G.sum(axis=1)) / \
        (by * N * (N + 1) * (N - 1))
    return -a


def _chains(members, R):
    """Split a set of cells (a tree along D8 paths) into polylines.
    Each chain starts at a cell with no upstream member and runs downstream
    until it leaves the set or joins an already-drawn chain (the junction
    cell is included so lines connect)."""
    mem = set(members)
    has_donor = set()
    for c in members:
        r = R[c]
        if r in mem:
            has_donor.add(r)
    visited = set()
    out = []
    for h in members:
        if h in has_donor:
            continue
        line = [h]
        visited.add(h)
        c = R[h]
        while c in mem:
            line.append(c)
            if c in visited:
                break
            visited.add(c)
            c = R[c]
        out.append(line)
    return out


def _union_find_groups(cells, R):
    parent = {c: c for c in cells}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for c in cells:
        r = R[c]
        if r in parent:
            ra, rb = find(c), find(r)
            if ra != rb:
                parent[ra] = rb
    groups = {}
    for c in cells:
        groups.setdefault(find(c), []).append(c)
    return list(groups.values())


def _ogr_driver(path):
    ext = os.path.splitext(path)[1].lower()
    return {".shp": "ESRI Shapefile", ".geojson": "GeoJSON",
            ".json": "GeoJSON"}.get(ext, "GPKG")


def _write_lines(path, wkt, layername, fields, records):
    """records: list of (list_of_coord_lists, dict_of_attributes)"""
    drv = ogr.GetDriverByName(_ogr_driver(path))
    if os.path.exists(path):
        drv.DeleteDataSource(path)
    ds = drv.CreateDataSource(path)
    if ds is None:
        raise KETError("Cannot create output file: %s" % path)
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    lyr = ds.CreateLayer(layername, srs, ogr.wkbMultiLineString)
    for name, ftype in fields:
        lyr.CreateField(ogr.FieldDefn(name, ftype))
    defn = lyr.GetLayerDefn()
    lyr.StartTransaction()
    for lines, attrs in records:
        g = ogr.Geometry(ogr.wkbMultiLineString)
        for coords in lines:
            ls = ogr.Geometry(ogr.wkbLineString)
            for x, y in coords:
                ls.AddPoint_2D(x, y)
            g.AddGeometry(ls)
        f = ogr.Feature(defn)
        f.SetGeometry(g)
        for k, v in attrs.items():
            f.SetField(k, v)
        lyr.CreateFeature(f)
    lyr.CommitTransaction()
    ds = None


def _write_raster(path, arr, gt, wkt, nodata=-9999.0):
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(path, arr.shape[1], arr.shape[0], 1, gdal.GDT_Float32,
                    options=["COMPRESS=DEFLATE", "TILED=YES"])
    if ds is None:
        raise KETError("Cannot create output raster: %s" % path)
    ds.SetGeoTransform(gt)
    ds.SetProjection(wkt)
    b = ds.GetRasterBand(1)
    out = np.where(np.isfinite(arr), arr, nodata).astype(np.float32)
    b.WriteArray(out)
    b.SetNoDataValue(nodata)
    ds = None


def run_ket(dem_path, out_knickzones, out_streams=None, out_rd=None, out_csv=None,
            stream_area_km2=0.1, direction="both",
            d_from=0.0, d_to=0.0, d_by=0.0,
            rd_threshold=0.0, auto_sd=1.0, min_length=20.0, min_height=0.0,
            rd_all_cells=False,
            log=_noop_log, progress=_noop_progress, is_canceled=_never_cancel):
    """Run the whole KET workflow.
    direction : 'both' (d/2 upstream along the main stem + d/2 downstream) or
                'down' (d downstream; KET v1 behaviour)
    d_from, d_to, d_by : metres; 0 = automatic (6, 36 and 3 cell sizes, i.e. the
                ratio of Hayakawa & Oguchi's 320-1760 m by 160 m on a 50-m DEM)
    rd_threshold : Rd (1/m) above which stream cells are knickzones;
                0 = automatic = mean + auto_sd * SD of Rd on the stream network
                (Zahra et al. 2017; Hayakawa & Oguchi 2006 used the mean SD)
    Returns a summary dict."""
    t0 = time.time()
    if direction not in ("both", "down"):
        raise KETError("direction must be 'both' or 'down'")
    z, valid, gt, wkt = read_dem(dem_path, log)
    nr, nc = z.shape
    dx, dy = abs(gt[1]), abs(gt[5])
    cell_area = dx * dy
    nvalid = int(valid.sum())
    log("DEM: %d x %d cells (%d valid), cell size %.3g x %.3g m" % (nc, nr, nvalid, dx, dy))
    if nvalid < 100:
        raise KETError("The DEM has too few valid cells.")
    if nvalid > 30e6:
        log("Warning: large DEM (%.1f M cells). This may take a long time and a lot of "
            "memory; consider clipping or resampling." % (nvalid / 1e6))

    # ---- scales d
    cs = max(dx, dy)
    auto = []
    if d_from <= 0:
        d_from = 6 * cs
        auto.append("from")
    if d_to <= 0:
        d_to = 36 * cs
        auto.append("to")
    if d_by <= 0:
        d_by = 3 * cs
        auto.append("by")
    if d_to <= d_from:
        raise KETError("'d to' must be larger than 'd from'.")
    nd = int(math.floor((d_to - d_from) / d_by + 1e-9)) + 1
    dlist = d_from + d_by * np.arange(nd)
    if nd < 2:
        raise KETError("Only one d value: decrease 'd by' or widen the d range.")
    log("Search direction: %s" % ("upstream (main stem) d/2 + downstream d/2"
                                  if direction == "both" else "downstream d (KET v1)"))
    log("Scales d: %g to %g m by %g m (%d values)%s" %
        (dlist[0], dlist[-1], d_by, nd, " [auto: %s]" % ", ".join(auto) if auto else ""))

    # ---- padded flat arrays
    ncp = nc + 2
    Vp = np.zeros((nr + 2, ncp), dtype=bool)
    Vp[1:-1, 1:-1] = valid
    Zp = np.zeros((nr + 2, ncp))
    Zp[1:-1, 1:-1] = np.where(valid, z, 0.0)
    Vp = Vp.ravel()
    Zp = Zp.ravel()

    log("[1/5] Filling depressions ...")
    F, _ = fill_depressions(Zp, Vp, ncp,
                            progress=lambda p: progress(0.35 * p),
                            is_canceled=is_canceled)
    log("[2/5] D8 flow direction, accumulation and main stems ...")
    R, L = d8_receivers(F, Vp, ncp, dx, dy)
    acc = flow_accumulation(F, R, Vp)
    U = upstream_mainstem(F, R, acc, Vp)
    Lup = np.where(U >= 0, L[np.maximum(U, 0)], 0.0)
    dout, dhead = flow_lengths(F, R, U, L, Vp)
    progress(50)

    area_cells = stream_area_km2 * 1e6 / cell_area
    stream = Vp & (acc >= area_cells)
    scells = np.flatnonzero(stream)
    log("[3/5] Stream network: contributing area >= %g km2 (%d cells) -> %d stream cells"
        % (stream_area_km2, int(math.ceil(area_cells)), scells.size))
    if scells.size < 10:
        raise KETError("Almost no stream cells were found. Make the stream threshold "
                       "(km2) smaller than the area of your DEM (%.3g km2)."
                       % (nvalid * cell_area / 1e6))

    log("[4/5] Multi-scale gradients Gd and relative steepness Rd ...")
    targets = np.flatnonzero(Vp) if rd_all_cells else scells
    Rd = np.full(Vp.size, np.nan)
    chunk = 500_000
    for i0 in range(0, targets.size, chunk):
        part = targets[i0:i0 + chunk]
        G = multiscale_gradient(F, R, U, L, Lup, part, dlist, direction, is_canceled)
        Rd[part] = relative_steepness(G, d_by)
        progress(50 + 35 * min(1.0, (i0 + chunk) / targets.size))

    # ---- threshold
    rds = Rd[scells]
    rds = rds[np.isfinite(rds)]
    if rds.size == 0:
        raise KETError("Rd could not be computed on the stream network "
                       "(the d range may be longer than the streams).")
    mu, sd = float(rds.mean()), float(rds.std())
    log("Rd on streams (1/m): mean %.3e, SD %.3e, median %.3e, 95%% %.3e, max %.3e"
        % (mu, sd, np.median(rds), np.percentile(rds, 95), rds.max()))
    if rd_threshold == 0:
        rd_threshold = mu + auto_sd * sd
        log("Rd threshold (auto = mean + %g SD): %.3e" % (auto_sd, rd_threshold))
    else:
        log("Rd threshold: %.3e" % rd_threshold)

    log("[5/5] Extracting knickzones (Rd > threshold) ...")
    kz = stream & np.isfinite(Rd) & (Rd > rd_threshold)
    kzcells = np.flatnonzero(kz).tolist()
    Rl = R.tolist()
    groups = _union_find_groups(kzcells, Rl)

    def xy(p):
        r, c = divmod(p, ncp)
        r -= 1
        c -= 1
        return (gt[0] + (c + 0.5) * gt[1], gt[3] + (r + 0.5) * gt[5])

    km2 = cell_area / 1e6
    records = []
    rows_csv = []
    kid = 0
    n_short = n_low = 0
    for g in groups:
        gs = set(g)
        length = sum(L[c] for c in g if Rl[c] in gs)
        if length <= min_length:
            n_short += 1
            continue
        zs = F[g]
        height = float(zs.max() - zs.min())
        if height < min_height:
            n_low += 1
            continue
        top = g[int(np.argmax(zs))]
        bot = g[int(np.argmin(zs))]
        rdv = Rd[g]
        kid += 1
        attrs = {
            "kz_id": kid,
            "length_m": round(float(length), 3),
            "height_m": round(height, 3),
            "gradient": round(height / length, 5) if length > 0 else None,
            "n_cells": len(g),
            "rd_mean": float(np.mean(rdv)),
            "rd_max": float(np.max(rdv)),
            "z_top": round(float(zs.max()), 3),
            "z_bottom": round(float(zs.min()), 3),
            "z_mean": round(float(zs.mean()), 3),
            "z_std": round(float(zs.std()), 3),
            "area_km2": round(float(acc[bot] * km2), 5),
            "dist_head_m": round(float(dhead[top]), 2),
            "dist_outlet_m": round(float(dout[bot]), 2),
            "x_top": xy(top)[0], "y_top": xy(top)[1],
            "x_bottom": xy(bot)[0], "y_bottom": xy(bot)[1],
        }
        lines = [[xy(p) for p in ch] for ch in _chains(g, Rl) if len(ch) > 1]
        records.append((lines, attrs))
        rows_csv.append(attrs)

    fields = [("kz_id", ogr.OFTInteger), ("length_m", ogr.OFTReal), ("height_m", ogr.OFTReal),
              ("gradient", ogr.OFTReal), ("n_cells", ogr.OFTInteger), ("rd_mean", ogr.OFTReal),
              ("rd_max", ogr.OFTReal), ("z_top", ogr.OFTReal), ("z_bottom", ogr.OFTReal),
              ("z_mean", ogr.OFTReal), ("z_std", ogr.OFTReal), ("area_km2", ogr.OFTReal),
              ("dist_head_m", ogr.OFTReal), ("dist_outlet_m", ogr.OFTReal),
              ("x_top", ogr.OFTReal), ("y_top", ogr.OFTReal),
              ("x_bottom", ogr.OFTReal), ("y_bottom", ogr.OFTReal)]
    _write_lines(out_knickzones, wkt, "knickzones", fields, records)
    log("Knickzones: %d written to %s (removed: %d with length <= %g m, %d with height < %g m)"
        % (kid, out_knickzones, n_short, min_length, n_low, min_height))

    if out_csv:
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=[f[0] for f in fields])
            w.writeheader()
            w.writerows(rows_csv)

    if out_streams:
        srec = []
        for sid, ch in enumerate(_chains(scells.tolist(), Rl), start=1):
            if len(ch) < 2:
                continue
            srec.append(([[xy(p) for p in ch]],
                         {"str_id": sid, "area_km2": round(float(acc[ch[-1]] * km2), 5)}))
        _write_lines(out_streams, wkt, "streams",
                     [("str_id", ogr.OFTInteger), ("area_km2", ogr.OFTReal)], srec)

    if out_rd:
        _write_raster(out_rd, Rd.reshape(nr + 2, ncp)[1:-1, 1:-1], gt, wkt)

    progress(100)
    dt = time.time() - t0
    log("Done in %.1f s" % dt)
    return {"knickzones": kid, "rd_threshold": rd_threshold, "d": dlist.tolist(),
            "direction": direction, "stream_cells": int(scells.size), "seconds": dt}


# ---------------------------------------------------------------------------
#  QGIS Processing wrapper (only defined when running inside QGIS)
# ---------------------------------------------------------------------------

try:
    from qgis.PyQt.QtCore import QCoreApplication
    from qgis.core import (QgsProcessingAlgorithm, QgsProcessingParameterRasterLayer,
                           QgsProcessingParameterNumber, QgsProcessingParameterBoolean,
                           QgsProcessingParameterEnum,
                           QgsProcessingParameterVectorDestination,
                           QgsProcessingParameterRasterDestination,
                           QgsProcessingParameterFileDestination,
                           QgsProcessingParameterDefinition,
                           QgsProcessingLayerPostProcessorInterface,
                           QgsProcessingException, QgsLineSymbol, QgsProcessing)
    _HAS_QGIS = True
except ImportError:
    _HAS_QGIS = False

if _HAS_QGIS:
    # enum names changed in QGIS 3.36; support both old and new versions
    try:
        from qgis.core import Qgis
        _DOUBLE = Qgis.ProcessingNumberParameterType.Double
        _LINE = Qgis.ProcessingSourceType.VectorLine
        _ADV = Qgis.ProcessingParameterFlag.Advanced
    except (ImportError, AttributeError):
        _DOUBLE = QgsProcessingParameterNumber.Double
        _LINE = QgsProcessing.TypeVectorLine
        _ADV = QgsProcessingParameterDefinition.FlagAdvanced

    class _LineStyler(QgsProcessingLayerPostProcessorInterface):
        _keep = []  # keep Python references alive until QGIS uses them

        def __init__(self, color, width):
            super().__init__()
            self.color, self.width = color, width

        def postProcessLayer(self, layer, context, feedback):
            try:
                sym = QgsLineSymbol.createSimple({"color": self.color, "width": str(self.width),
                                                  "capstyle": "round"})
                layer.renderer().setSymbol(sym)
                layer.triggerRepaint()
            except Exception:
                pass

        @classmethod
        def make(cls, color, width):
            s = cls(color, width)
            cls._keep.append(s)
            return s

    class KnickzoneExtractionAlgorithm(QgsProcessingAlgorithm):
        DEM = "DEM"
        AREA = "STREAM_AREA_KM2"
        DIR = "DIRECTION"
        DFROM = "D_FROM"
        DTO = "D_TO"
        DBY = "D_BY"
        TH = "RD_THRESHOLD"
        KSD = "AUTO_SD"
        MINLEN = "MIN_LENGTH"
        MINH = "MIN_HEIGHT"
        ALL = "RD_ALL_CELLS"
        OUT_KZ = "KNICKZONES"
        OUT_STR = "STREAMS"
        OUT_RD = "RD"
        OUT_CSV = "TABLE"
        DIRS = ("both", "down")

        def tr(self, s):
            return QCoreApplication.translate("KET", s)

        def createInstance(self):
            return KnickzoneExtractionAlgorithm()

        def name(self):
            return "ket_knickzone_extraction"

        def displayName(self):
            return self.tr("Knickzone Extraction (KET)")

        def group(self):
            return self.tr("KET")

        def groupId(self):
            return "ket"

        def shortHelpString(self):
            return self.tr(
                "<p><b>Knickzone Extraction Tool (KET) v%s</b><br>"
                "DEM を入れて実行するだけで遷急区間 (knickzone) を抽出します。"
                "Just choose a DEM and press Run.</p>"
                "<p><b>DEM</b>: any GDAL raster. Projected CRS in metres is best; "
                "lat/lon DEMs (e.g. GSI 基盤地図) are reprojected to UTM automatically.</p>"
                "<p><b>Stream threshold</b>: minimum contributing area (km²) of streams.</p>"
                "<p><b>Search direction</b>: <i>upstream + downstream</i> takes d/2 up the main "
                "stem (largest contributing area; if equal, the lower cell) and d/2 downstream, "
                "as in Hayakawa &amp; Oguchi (2006, 2009). Upstream search stops at the ridge. "
                "<i>Downstream only</i> takes d downstream from the cell (KET v1).</p>"
                "<p><b>d from / to / by</b>: measurement lengths (m). 0 = automatic: 6, 36 and 3 "
                "cell sizes (= 320–1760 m by 160 m on a 50-m DEM).</p>"
                "<p><b>Rd</b> = −(slope of Gd against d), unit 1/m; positive = locally steep. "
                "Stream cells with Rd above the threshold are knickzones. 0 = automatic "
                "(mean + 1 SD of Rd on the streams). Published values: 1.42×10⁻⁵ "
                "(Hayakawa &amp; Oguchi 2006, 50-m DEM), 1.92×10⁻⁴ (Zahra et al. 2017, 10-m DEM).</p>"
                "<p><b>Minimum length / height</b>: smaller knickzones are removed.</p>"
                "<p>Outputs: knickzone lines (length, height, gradient, Rd, drainage area, distance "
                "from head and to outlet), streams, Rd raster, CSV table.</p>"
                "<p>Hayakawa &amp; Oguchi (2006) Geomorphology 78: 90–106; (2009) Geomorphology "
                "111: 27–37; Zahra et al. (2017) Open Geosciences 9: 73–88.</p>" % __version__)

        def initAlgorithm(self, config=None):
            self.addParameter(QgsProcessingParameterRasterLayer(self.DEM, self.tr("DEM")))
            self.addParameter(QgsProcessingParameterNumber(
                self.AREA, self.tr("Stream threshold: contributing area (km²)"),
                _DOUBLE, 0.1, minValue=0.0))
            self.addParameter(QgsProcessingParameterEnum(
                self.DIR, self.tr("Search direction"),
                options=[self.tr("Upstream (main stem) d/2 + downstream d/2  [SAKE]"),
                         self.tr("Downstream only, d  [KET v1]")],
                defaultValue=0))
            self.addParameter(QgsProcessingParameterNumber(
                self.TH, self.tr("Rd threshold, 1/m (0 = automatic)"), _DOUBLE, 0.0))
            self.addParameter(QgsProcessingParameterNumber(
                self.MINLEN, self.tr("Minimum knickzone length (m)"), _DOUBLE, 20.0, minValue=0.0))
            self.addParameter(QgsProcessingParameterNumber(
                self.MINH, self.tr("Minimum knickzone height (m)"), _DOUBLE, 0.0, minValue=0.0))
            adv = [
                QgsProcessingParameterNumber(self.DFROM, self.tr("d from (m, 0 = 6 cells)"),
                                             _DOUBLE, 0.0, minValue=0.0),
                QgsProcessingParameterNumber(self.DTO, self.tr("d to (m, 0 = 36 cells)"),
                                             _DOUBLE, 0.0, minValue=0.0),
                QgsProcessingParameterNumber(self.DBY, self.tr("d by (m, 0 = 3 cells)"),
                                             _DOUBLE, 0.0, minValue=0.0),
                QgsProcessingParameterNumber(self.KSD, self.tr("Automatic threshold: mean + k × SD, k ="),
                                             _DOUBLE, 1.0, minValue=0.0),
                QgsProcessingParameterBoolean(self.ALL, self.tr("Compute Rd for every DEM cell (slow)"),
                                              False),
            ]
            for p in adv:
                p.setFlags(p.flags() | _ADV)
                self.addParameter(p)
            self.addParameter(QgsProcessingParameterVectorDestination(
                self.OUT_KZ, self.tr("Knickzones"), _LINE))
            self.addParameter(QgsProcessingParameterVectorDestination(
                self.OUT_STR, self.tr("Streams"), _LINE, optional=True))
            self.addParameter(QgsProcessingParameterRasterDestination(
                self.OUT_RD, self.tr("Rd (relative steepness) raster"), optional=True))
            self.addParameter(QgsProcessingParameterFileDestination(
                self.OUT_CSV, self.tr("Knickzone table (CSV)"), "CSV files (*.csv)",
                optional=True, createByDefault=False))

        def processAlgorithm(self, parameters, context, feedback):
            dem = self.parameterAsRasterLayer(parameters, self.DEM, context)
            if dem is None:
                raise QgsProcessingException("Invalid DEM")
            out_kz = self.parameterAsOutputLayer(parameters, self.OUT_KZ, context)
            out_str = self.parameterAsOutputLayer(parameters, self.OUT_STR, context) or None
            out_rd = self.parameterAsOutputLayer(parameters, self.OUT_RD, context) or None
            out_csv = self.parameterAsFileOutput(parameters, self.OUT_CSV, context) or None
            dval = lambda k: self.parameterAsDouble(parameters, k, context)
            try:
                run_ket(dem.source(), out_kz, out_str, out_rd, out_csv,
                        stream_area_km2=dval(self.AREA),
                        direction=self.DIRS[self.parameterAsEnum(parameters, self.DIR, context)],
                        d_from=dval(self.DFROM), d_to=dval(self.DTO), d_by=dval(self.DBY),
                        rd_threshold=dval(self.TH), auto_sd=dval(self.KSD),
                        min_length=dval(self.MINLEN), min_height=dval(self.MINH),
                        rd_all_cells=self.parameterAsBool(parameters, self.ALL, context),
                        log=feedback.pushInfo, progress=feedback.setProgress,
                        is_canceled=feedback.isCanceled)
            except KETError as e:
                raise QgsProcessingException(str(e))

            for path, color, width, name in ((out_str, "#3a7bd5", 0.4, "Streams"),
                                             (out_kz, "#F71C54", 1.2, "Knickzones")):
                if path and context.willLoadLayerOnCompletion(path):
                    det = context.layerToLoadOnCompletionDetails(path)
                    det.name = name
                    det.setPostProcessor(_LineStyler.make(color, width))
            res = {self.OUT_KZ: out_kz}
            if out_str:
                res[self.OUT_STR] = out_str
            if out_rd:
                res[self.OUT_RD] = out_rd
            if out_csv:
                res[self.OUT_CSV] = out_csv
            return res


# ---------------------------------------------------------------------------
#  Command line
# ---------------------------------------------------------------------------

def _main():
    import argparse
    gdal.UseExceptions()
    ogr.UseExceptions()
    osr.UseExceptions()
    ap = argparse.ArgumentParser(description="KET %s - Knickzone Extraction Tool (numpy/GDAL)"
                                 % __version__)
    ap.add_argument("dem", help="input DEM (any GDAL raster)")
    ap.add_argument("-o", "--outdir", default="ket_output", help="output folder")
    ap.add_argument("--area", type=float, default=0.1, help="stream threshold, km2 (default 0.1)")
    ap.add_argument("--direction", choices=("both", "down"), default="both",
                    help="both = d/2 up (main stem) + d/2 down [default]; down = d downstream (KET v1)")
    ap.add_argument("--dfrom", type=float, default=0, help="d from, m (0 = 6 cells)")
    ap.add_argument("--dto", type=float, default=0, help="d to, m (0 = 36 cells)")
    ap.add_argument("--dby", type=float, default=0, help="d by, m (0 = 3 cells)")
    ap.add_argument("--th", type=float, default=0, help="Rd threshold, 1/m (0 = mean + k SD)")
    ap.add_argument("--k", type=float, default=1.0, help="k for the automatic threshold")
    ap.add_argument("--minlen", type=float, default=20, help="minimum knickzone length, m")
    ap.add_argument("--minheight", type=float, default=0, help="minimum knickzone height, m")
    ap.add_argument("--all", action="store_true", help="compute Rd for all cells")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    j = lambda f: os.path.join(a.outdir, f)
    last = [-1]

    def prog(p):
        if int(p) // 10 != last[0]:
            last[0] = int(p) // 10
            print("  %3d %%" % (last[0] * 10), flush=True)

    try:
        run_ket(a.dem, j("knickzones.gpkg"), j("streams.gpkg"), j("rd.tif"), j("knickzones.csv"),
                stream_area_km2=a.area, direction=a.direction,
                d_from=a.dfrom, d_to=a.dto, d_by=a.dby,
                rd_threshold=a.th, auto_sd=a.k, min_length=a.minlen, min_height=a.minheight,
                rd_all_cells=a.all, progress=prog)
    except KETError as e:
        raise SystemExit("ERROR: %s" % e)


if __name__ == "__main__":
    _main()
