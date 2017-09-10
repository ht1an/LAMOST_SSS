def footprint(cen_ra, cen_dec, n=1000):
    import numpy as np
    if cen_dec < -85 or cen_dec > 85:
        print("primary version")
        print("cen_dec should be in range [-85.0, 85.0]")
        return
    r = 2.5
    d1 = cen_dec - r
    d2 = cen_dec + r
    d = np.linspace(d1, d2, n)
    cos_r = np.cos(np.deg2rad(r))
    cos_d = np.cos(np.deg2rad(d))
    cos_d0 = np.cos(np.deg2rad(cen_dec))
    sin_d = np.sin(np.deg2rad(d))
    sin_d0 = np.sin(np.deg2rad(cen_dec))
    cos_dra = (cos_r - sin_d * sin_d0) / (cos_d * cos_d0)
    da = np.arccos(cos_dra)
    da = np.rad2deg(da)
    a1 = cen_ra - da
    a2 = cen_ra + da
    x = np.hstack([a1, a2[::-1][1:]])
    y1 = d
    y2 = d[::-1][1:]
    y = np.hstack([y1, y2])
    return x, y
