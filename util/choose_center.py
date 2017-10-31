def choose_center(ra, dec, ftycho, k=1, max_radius=2.5):
    # assume the range of 'ra' is less than 20-30 degrees,
    # the number of bright stars (center stars) from tycho-catalog is small 
    # in the region of the input catalog, 
    # so looping over all candidates is adopted.
    # the radius of the focal plane is 2.5 degrees by default (LAMOST):
    # max_radius=2.5
    # it should not be too large, otherwise the efficiency may be influenced.
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
    # Test the case if the input catalog distribute across ra = 0
    if ra_max - ra_min > 180:
        ra_min_pos = 0
        ra_max_pos = np.max(ra[ra < 180])
        ra_min_neg = np.min(ra[ra > 180])
        ra_max_neg = 360
        idx = (ra0 >= ra_min_neg) + (ra0 <= ra_max_pos)
        idx = idx * (dec0 >= dec_min) * (dec0 <= dec_max)
    else:
        idx = (ra0 >= ra_min) * (ra0 <= ra_max) * (dec0 >= dec_min) * (dec0 <= dec_max)
    # ugly python indent style
    ts = t[idx]
    ra0  = ra0[idx]
    dec0 = dec0[idx]
    HIP  = HIP[idx]

    # n is the number of bright star (center star) candidates
    n = len(ts)
    # index array, each element store the index that covered by the chosen center star
    idarray = np.empty(n * len(ra), dtype=bool).reshape(n, len(ra))
    # number array, star number that covered by the chosen center star
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

    # return "k" bright stars as center stars
    result = [None] * k
    for i in range(k):
        j = l[i]
        result[i] = [ra0[j], dec0[j], idarray[i], HIP[j], m[i]]

    return result
