# Calculate the Great-circle distance
# returned "dist" in radian
# see:
# https://en.wikipedia.org/wiki/Great-circle_distance
def gcirc(lon1, lat1, lon2, lat2, u=0):
    # u = 0, in radian, otherwise in degree
    import numpy as np
    if u != 0:
        return gcirc(np.deg2rad(lon1), np.deg2rad(lat1), np.deg2rad(lon2), np.deg2rad(lat2), 0)

    delta_lon_2 = 0.5 * (lon1 - lon2)
    delta_lat_2 = 0.5 * (lat1 - lat2)
    sin_dlon_2 = np.sin(delta_lon_2)
    sin_dlat_2 = np.sin(delta_lat_2)
    sin_dist = np.sqrt(sin_dlat_2 * sin_dlat_2 + 
            np.cos(lat1) * np.cos(lat2) * sin_dlon_2 * sin_dlon_2)
    dist = 2.0 * np.arcsin(sin_dist)
    return dist
