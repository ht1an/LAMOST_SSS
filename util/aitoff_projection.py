# aitoff projection
# see:
# https://en.wikipedia.org/wiki/Aitoff_projection
def aitoff_projection(theta, phi):
    import numpy as np
    # theta, phi in radian
    theta = theta - np.pi
    cos_phi = np.cos(phi)
    denom = np.sqrt(1 + cos_phi * np.cos(theta/2))
    x = 180 * cos_phi * np.sin(theta/2) / denom
    x = x + 180
    y = 90 * np.sin(phi) / denom
    return x,y
