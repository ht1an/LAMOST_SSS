# calculate the circle trajectory 
# when given the center point (lon1, lat1) and the distance r
def footprint(lon1, lat1, n=1000, r=2.5, u=1):
    # u = 0, in radian, otherwise in degree
    # rotate the coordinate so that point (lon1, lat1) 
    # locates at (0, 0) in new coordinate
    # R2(-lat1) * R3(lon1)
    import numpy as np
    if u != 0:
        r = np.deg2rad(r)
        x, y = footprint(np.deg2rad(lon1), np.deg2rad(lat1), n, r, u=0)
        return np.rad2deg(x), np.rad2deg(y)

    # in coord1, (lon1, lat1), (lon2, lat2)
    # in coord2, (0, 0), (a, b)
    # hav(x) = sin(0.5*x)^2
    # haversine formula behave better for small r
    b = np.linspace(-r, r, n)
    sin_b = np.sin(b)
    # cos_b = np.cos(b)
    cos_b = np.sqrt(1.0 - sin_b * sin_b)
    hav_r = np.sin(0.5 * r)
    hav_r = hav_r * hav_r
    hav_b = np.sin(0.5 * b)
    hav_b = hav_b * hav_b
    hav_a = (hav_r - hav_b) / cos_b
    cos_a = 1.0 - 2 * hav_a
    sin_a = np.sqrt(1.0 - cos_a * cos_a) # use only the positive part

    # Transform back to coord1
    sin_lat1 = np.sin(lat1)
    cos_lat1 = np.sqrt(1.0 - sin_lat1 * sin_lat1)
    sin_lat2 = sin_lat1 * cos_a * cos_b + cos_lat1 * sin_b
    lat2 = np.arcsin(sin_lat2)
    cos_lat2 = np.sqrt(1.0 - sin_lat2 * sin_lat2)
    sin_dlon_cos_lat2 = sin_a * cos_b
    cos_dlon_cos_lat2 = cos_lat1 * cos_a * cos_b - sin_lat1 * sin_b
    dlon = np.arctan2(sin_dlon_cos_lat2, cos_dlon_cos_lat2)

    # positive part of sin_a
    lon2_pos = lon1 + dlon
    # negative part of sin_a, reverse order
    lon2_neg = lon1 - dlon[::-1][1:]

    x = np.hstack([lon2_pos, lon2_neg])
    y = np.hstack([lat2, lat2[::-1][1:]])

    x = np.where(x < 0.0, x + 2*np.pi, x)
    x = np.where(x > 2*np.pi, x - 2*np.pi, x)
    return x, y
