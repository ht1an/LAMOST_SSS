import numpy as np
# gp_ra, gp_dec, lp: given by Reid & Brunthaler 2004
# http://adsabs.harvard.edu/abs/2004ApJ...616..872R
gp_ra = 192.859508
gp_dec = 27.128336
lp = 122.932
sin_gp_dec = np.sin(np.deg2rad(gp_dec))
cos_gp_dec = np.cos(np.deg2rad(gp_dec))
def galactic1(ra, dec):
#    import numpy as np
#    gp_ra = 192.859508
#    gp_dec = 27.128336
#    lp = 122.932
#    sin_gp_dec = np.sin(np.deg2rad(gp_dec))
#    cos_gp_dec = np.cos(np.deg2rad(gp_dec))

    d_ra = np.deg2rad(ra - gp_ra)
    sin_dec = np.sin(np.deg2rad(dec))
    # cos_dec = np.cos(np.deg2rad(dec))
    cos_dec = np.sqrt(1.0 - sin_dec * sin_dec)
    cos_dra = np.cos(d_ra)
    sin_dra = np.sin(d_ra)

    sin_b = sin_gp_dec * sin_dec + cos_gp_dec * cos_dec * np.cos(d_ra)
    b = np.rad2deg(np.arcsin(sin_b))
    cos_b = np.sqrt(1.0 - sin_b * sin_b)

    sin_l0 = cos_dec * np.sin(d_ra) / cos_b
    cos_l0 = (sin_dec - sin_gp_dec * sin_b) / (cos_gp_dec * cos_b)
    l0 = np.rad2deg(np.arctan2(sin_l0, cos_l0))
    l  = lp - l0
    l  = np.where(l < 0.0, l + 360.0, l)
    # for type consistency
    l  = l + 0.0

    return (l, b)

def galactic2(l, b):
    sin_b = np.sin(np.deg2rad(b))
    cos_b = np.sqrt(1.0 - sin_b * sin_b)
    d_l = np.deg2rad(lp - l)
    sin_dec = sin_b * sin_gp_dec + cos_b * cos_gp_dec * np.cos(d_l)
    dec = np.rad2deg(np.arcsin(sin_dec))
    cos_dec = np.sqrt(1.0 - sin_dec * sin_dec)
    sin_ra0 = cos_b * np.sin(d_l) / cos_dec
    cos_ra0 = (sin_b - sin_gp_dec * sin_dec) / (cos_gp_dec * cos_dec)
    ra0 = np.rad2deg(np.arctan2(sin_ra0, cos_ra0))
    ra  = ra0 + gp_ra
    ra  = np.where(ra >= 360.0, ra - 360.0, ra)
    # for type consistency
    ra  = ra + 0.0

    return (ra, dec)

def galactic(lon, lat, j):
    if j == 1:
        return galactic1(lon, lat)
    elif j == 2:
        return galactic2(lon, lat)
