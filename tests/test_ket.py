"""Consistency tests for ket_qgis (run:  python tests/test_ket.py).
Builds a synthetic DEM and checks the vectorised routines against
straightforward cell-by-cell implementations."""
import os, sys, tempfile
import numpy as np
from osgeo import gdal, osr

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ket_qgis as K  # noqa: E402

gdal.UseExceptions()


def synthetic_dem(path, nr=300, nc=450, cs=10.0, seed=0):
    """Valley draining to x=0 with a 60-m step (knickzone) at x = 2000-2300 m."""
    rng = np.random.default_rng(seed)
    x = (np.arange(nc) + 0.5) * cs
    y = (np.arange(nr) + 0.5) * cs
    X, Y = np.meshgrid(x, y)
    prof = 0.02 * X + 40 * (1 - np.exp(-X / 2000)) + np.clip((X - 2000) / 300, 0, 1) * 60
    centre = nr * cs / 2 + 150 * np.sin(X / 700)
    valley = 0.25 * np.abs(Y - centre) + 0.0002 * (Y - centre) ** 2
    kx = np.fft.fftfreq(nc)[None, :]
    ky = np.fft.fftfreq(nr)[:, None]
    k = np.hypot(kx, ky)
    k[0, 0] = 1
    n = np.real(np.fft.ifft2(np.fft.fft2(rng.normal(size=(nr, nc))) / k ** 1.6))
    z = 500 + prof + valley + n / n.std() * 4
    ds = gdal.GetDriverByName("GTiff").Create(path, nc, nr, 1, gdal.GDT_Float32)
    ds.SetGeoTransform((400000, cs, 0, 4800000, 0, -cs))
    s = osr.SpatialReference()
    s.ImportFromEPSG(32654)
    ds.SetProjection(s.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(z.astype(np.float32))
    ds = None
    return path


def prepare(path):
    z, valid, gt, wkt = K.read_dem(path, lambda m: None)
    nr, nc = z.shape
    ncp = nc + 2
    Vp = np.zeros((nr + 2, ncp), bool)
    Vp[1:-1, 1:-1] = valid
    Zp = np.zeros((nr + 2, ncp))
    Zp[1:-1, 1:-1] = z
    Vp, Zp = Vp.ravel(), Zp.ravel()
    F, _ = K.fill_depressions(Zp, Vp, ncp)
    R, L = K.d8_receivers(F, Vp, ncp, abs(gt[1]), abs(gt[5]))
    acc = K.flow_accumulation(F, R, Vp)
    U = K.upstream_mainstem(F, R, acc, Vp)
    Lup = np.where(U >= 0, L[np.maximum(U, 0)], 0.0)
    return F, Vp, R, L, acc, U, Lup


def brute_endpoint(F, P, step, c, t):
    cur, dist = c, 0.0
    while dist < t and P[cur] >= 0:
        dist += step(cur)
        cur = P[cur]
    return F[cur], dist


def test_all():
    tmp = tempfile.mkdtemp()
    dem = synthetic_dem(os.path.join(tmp, "synth.tif"))
    F, Vp, R, L, acc, U, Lup = prepare(dem)
    idx = np.flatnonzero(Vp)

    # filled DEM drains: every receiver is strictly lower
    r = R[idx]
    ok = r >= 0
    assert (F[r[ok]] < F[idx[ok]]).all()
    assert acc.max() == Vp.sum()

    # main stem: U is a donor with maximal acc; ties -> lowest elevation
    donors = {}
    for c in idx[ok]:
        donors.setdefault(int(R[c]), []).append(int(c))
    for rc, ds in donors.items():
        best = sorted(ds, key=lambda d: (-acc[d], F[d]))[0]
        assert U[rc] == best, rc
    assert all(U[c] == -1 for c in idx if int(c) not in donors)
    n_tie = sum(1 for ds in donors.values()
                if len(ds) > 1 and sorted(acc[ds])[-1] == sorted(acc[ds])[-2])
    print("main-stem pointer OK (%d receivers with tied donor areas)" % n_tie)

    d = np.arange(60.0, 361.0, 30.0)
    cells = np.random.default_rng(1).choice(idx, 1500, replace=False)
    for direction in ("down", "both"):
        G = K.multiscale_gradient(F, R, U, L, Lup, cells, d, direction)
        Gb = np.full_like(G, np.nan)
        for i, c in enumerate(cells):
            for k, dd in enumerate(d):
                if direction == "down":
                    e2, D2 = brute_endpoint(F, R, lambda q: L[q], c, dd)
                    e1, D1 = F[c], 0.0
                else:
                    e2, D2 = brute_endpoint(F, R, lambda q: L[q], c, dd / 2)
                    e1, D1 = brute_endpoint(F, U, lambda q: L[U[q]], c, dd / 2)
                if D1 + D2 > 0:
                    Gb[i, k] = (e1 - e2) / (D1 + D2)
        assert np.array_equal(np.isnan(G), np.isnan(Gb))
        assert np.nanmax(np.abs(G - Gb)) < 1e-12
        Rd = K.relative_steepness(G, 30.0)
        Rreg = np.array([-np.polyfit(d, g, 1)[0] if np.isfinite(g).all() else np.nan for g in G])
        assert np.nanmax(np.abs(Rd - Rreg)) < 1e-12
        print("Gd and Rd (%s) match the cell-by-cell computation" % direction)

    # end-to-end: the knickzone must be found around x = 2000-2300 m
    for direction in ("both", "down"):
        out = os.path.join(tmp, "kz_%s.gpkg" % direction)
        res = K.run_ket(dem, out, out_csv=out + ".csv", stream_area_km2=0.3,
                        direction=direction, d_from=60, d_to=720, d_by=60,
                        rd_threshold=1e-4, log=lambda m: None)
        import csv
        rows = list(csv.DictReader(open(out + ".csv")))
        xs = [(float(r["x_top"]) - 400000, float(r["x_bottom"]) - 400000, float(r["height_m"]))
              for r in rows]
        trunk = [t for t, r in zip(xs, rows) if float(r["area_km2"]) > 5]
        assert len(trunk) == 1, trunk
        top, bot, h = trunk[0]
        print("%s: %d knickzones; trunk knickzone x = %.0f -> %.0f m, height %.1f m "
              "(true step: 2300 -> 2000 m, 60 m)" % (direction, res["knickzones"], top, bot, h))
        assert 1950 < bot < 2150 and 2150 < top < 2400 and h > 40


if __name__ == "__main__":
    test_all()
    print("all tests passed")
