def read_target(filename):
    import numpy as np

    dt = np.dtype([
        ('CATNAME' , 'S5'),
        ('X'       , float),
        ('Y'       , float),
        ('LAMBDA'  , float),
        ('GAMMA'   , float),
        ('OBJID'   , 'S19'),
        ('RA'      , float),
        ('DEC'     , float),
        ('MAG0'    , float),
        ('PLATE'   , int  ),
        ('FIBER'   , int  ),
        ('PACK'    , int  ),
        ('OBJTYPE' , 'S4' ),
        ('CAT0'    , 'S10'),
        ('MAGTYPE0', 'S5' ),
        ('MAG1'    , float),
        ('MAG2'    , float),
        ('MAG3'    , float),
        ('MAG4'    , float),
        ('MAG5'    , float),
        ('MAG6'    , float),
        ('MAG7'    , float),
        ('PRI'     , int  ),
        ('CAT1'    , 'S10')
        ])

    x = np.genfromtxt(filename, skip_header=2, dtype=dt, usecols=range(24))
    return x
