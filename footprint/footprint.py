def footprint(lon1, lat1, n=1000, r=2.5, u=1):
    # u = 0, in radian, otherwise in degree
    import numpy as np
    if u != 0:
        r = np.deg2rad(r)
        x, y = footprint(np.deg2rad(lon1), np.deg2rad(lat1), n, r, u=0)
        return np.rad2deg(x), np.rad2deg(y)
    if lat1 > np.deg2rad(85) or lat1 < np.deg2rad(-85):
        print('too high to calculate')
        return
    lat2_l = lat1 - r
    lat2_u = lat1 + r
    lat2 = np.linspace(lat2_l, lat2_u, n)
    delta_lat_2 = 0.5 * np.linspace(-r, r, n)
    sin_dlat_2 = np.sin(delta_lat_2)
    cos_lat1 = np.cos(lat1)
    cos_lat2 = np.cos(lat2)
    sin_r_2 = np.sin(0.5*r)
    #sin_dlon_2 = (sin_r_2 * sin_r_2 - sin_dlat_2 * sin_dlat_2) / (cos_lat1 * cos_lat2)
    sin_dlon_2 = (sin_r_2 + sin_dlat_2) * (sin_r_2 - sin_dlat_2) / (cos_lat1 * cos_lat2)
    dlon = 2.0 * np.arcsin(np.sqrt(sin_dlon_2))
    lon2_l = lon1 - dlon
    lon2_u = lon1 + dlon
    x = np.hstack([lon2_l, lon2_u[::-1][1:]])
    y = np.hstack([lat2, lat2[::-1][1:]])
    x = np.where(x < 0.0, x + 2*np.pi, x)
    x = np.where(x > 2*np.pi, x - 2*np.pi, x)
    return x, y
