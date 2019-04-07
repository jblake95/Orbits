"""
Finding sensible ranges for GEO orbital params
"""

import argparse as ap

def argParse():
    """
    Argument parser settings
    
    Parameters
    ----------
    None
    
    Returns
    -------
    args : array-like
        Array of command line arguments
    """
    parser = ap.ArgumentParser()
    
    parser.add_argument('cat_path',
                        help='path to file containing GEO tle catalog',
                        type=str)
    
    return parser.parse_args()

if __name__ == "__main__":
    
    args = argParse()
    
    with open(args.cat_path) as f:
        lines = f.readlines()
    lines = [i.strip() for i in lines]

    i = 0
    lines1 = []
    lines2 = []
    while i < len(lines):
        lines1.append(lines[i+1])
        lines2.append(lines[i+2])
        i += 3
    
    mmdot = []
    for line1 in lines1:
        mmdot.append(float(line1[33:43]))
    
    incl = []
    raan = []
    e = []
    argp = []
    anom = []
    motion = []
    for line2 in lines2:
        incl.append(float(line2[8:16]))
        raan.append(float(line2[17:25]))
        e.append(float(line2[26:33]))
        argp.append(float(line2[34:42]))
        anom.append(float(line2[43:51]))
        motion.append(float(line2[52:63]))
    print(sorted(incl))
    print('Ranges for GEO params:\n'
          '{} < mmdot < {}\n'
          '{} < incl < {}\n'
          '{} < raan < {}\n'
          '{} < e < {}\n'
          '{} < argp < {}\n'
          '{} < anom < {}\n'
          '{} < motion < {}'.format(str(min(mmdot)),
                                    str(max(mmdot)),
                                    str(min(incl)),
                                    str(max(incl)),
                                    str(min(raan)),
                                    str(max(raan)),
                                    str(min(e)),
                                    str(max(e)),
                                    str(min(argp)),
                                    str(max(argp)),
                                    str(min(anom)),
                                    str(max(anom)),
                                    str(min(motion)),
                                    str(max(motion))))
    
