def choose_center(ra, dec, ftycho):
    import numpy as np
    from astropy.io import fits
    import gcirc
    # read catalog from tycho-catalog
    h = fits.open(ftycho)
    t = h[1].data

    ra0  = t['RAmdeg']
    dec0 = t['DEmdeg']
    HIP  = t['HIP']
    idx  = np.isfinite(ra0) * np.isfinite(dec0) * (HIP > 0)
    t    = t[idx]
    ra0  = ra0[idx]
    dec0 = dec0[idx]
    HIP  = HIP[idx]

    ra_min  = np.min(ra)
    ra_max  = np.max(ra)
    dec_min = np.min(dec)
    dec_max = np.max(dec)
    idx = (ra0 >= ra_min) * (ra0 <= ra_max) * (dec0 >= dec_min) * (dec0 <= dec_max)
    ts = t[idx]
    ra0  = ra0[idx]
    dec0 = dec0[idx]
    HIP  = HIP[idx]
    n = len(ts)

    max_radius = 2.5
    mmax = 0
    k = 0

    for i in range(n):
        cen_ra  = ra0[i]
        cen_dec = dec0[i]
        d = gcirc.gcirc(ra, dec, cen_ra, cen_dec, u=1)
        idx = (d <= max_radius)
        m = np.sum(idx)
        if (mmax < m):
            mmax = m
            k = i
            id0 = idx

    return (ra0[k], dec0[k], id0, HIP[k], mmax)

def choose_center2(ra, dec, ftycho, k):
    import numpy as np
    from astropy.io import fits
    import gcirc
    # read catalog from tycho-catalog
    h = fits.open(ftycho)
    t = h[1].data

    ra0  = t['RAmdeg']
    dec0 = t['DEmdeg']
    HIP  = t['HIP']
    idx  = np.isfinite(ra0) * np.isfinite(dec0) * (HIP > 0)
    t    = t[idx]
    ra0  = ra0[idx]
    dec0 = dec0[idx]
    HIP  = HIP[idx]

    ra_min  = np.min(ra)
    ra_max  = np.max(ra)
    dec_min = np.min(dec)
    dec_max = np.max(dec)
    idx = (ra0 >= ra_min) * (ra0 <= ra_max) * (dec0 >= dec_min) * (dec0 <= dec_max)
    ts = t[idx]
    ra0  = ra0[idx]
    dec0 = dec0[idx]
    HIP  = HIP[idx]
    n = len(ts)

    max_radius = 2.5
    mmax = 0
    idarray = np.empty(n * len(ra), dtype=bool).reshape(n, len(ra))
    m = np.zeros(n, dtype=int)
    l = np.arange(n)

    for i in l:
        cen_ra  = ra0[i]
        cen_dec = dec0[i]
        d = gcirc.gcirc(ra, dec, cen_ra, cen_dec, u=1)
        idx = (d <= max_radius)
        idarray[i] = idx
        m[i] = np.sum(idx)

    idx = (m > 0)
    l = l[idx]
    m = m[idx]
    idarray = idarray[idx]

    idx = np.argsort(-m)
    l = l[idx]
    m = m[idx]
    idarray = idarray[idx]

    result = [None] * k
    for i in range(k):
        j = l[i]
        result[i] = [ra0[j], dec0[j], idarray[i], HIP[j], m[i]]

    return result
