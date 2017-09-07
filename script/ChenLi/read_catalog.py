def read_catalog(filename, dtype=None):
    import numpy as np

    dt = dtype
    if dtype is not None:
        dt = dtype
    else:
        dt = np.dtype([('RADEG', float), ('DECDEG', float), ('MAG0', float), ('PRI', int), ('OBJID', 'S19')])

    x = np.genfromtxt(filename, dtype=dt, delimiter=',', names=True)
    return x
