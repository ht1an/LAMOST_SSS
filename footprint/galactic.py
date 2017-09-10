def galactic(ra, dec):
    import numpy as np
    gp_ra = 192.859508
    gp_dec = 27.128336
    lp = 122.932
    
    sin_gp_dec = np.sin(np.deg2rad(gp_dec))
    cos_gp_dec = np.cos(np.deg2rad(gp_dec))

    ra0 = np.deg2rad(ra - gp_ra)
    sin_dec = np.sin(np.deg2rad(dec))
    # cos_dec = np.cos(np.deg2rad(dec))
    cos_dec = np.sqrt(1.0 - sin_dec * sin_dec)
    cos_ra0 = np.cos(ra0)
    sin_ra0 = np.sin(ra0)

    sin_b = sin_gp_dec * sin_dec + cos_gp_dec * cos_dec * np.cos(ra0)
    b = np.rad2deg(np.arcsin(sin_b))
    cos_b = np.sqrt(1.0 - sin_b * sin_b)

    sin_l0 = cos_dec * np.sin(ra0) / cos_b
    cos_l0 = (sin_dec - sin_gp_dec * sin_b) / (cos_gp_dec * cos_b)
    l0 = np.rad2deg(np.arccos(cos_l0))

#    idx = (sin_l0 < 0.0)
#    l0[idx] = 360.0 - l0[idx]
#    l = lp - l0
#    idx2 = (l < 0.0)
#    l[idx2] = l[idx2] + 360.0

    l0 = np.where(sin_l0 < 0.0, 360.0 - l0, l0)
    l  = lp - l0
    l  = np.where(l < 0.0, l + 360.0, l)
#    if (sin_l0 < 0):
#        l0 = 360.0 - l0
#    l = lp - l0
#    if (l < 0):
#        l = l + 360.0

    return (l, b)
