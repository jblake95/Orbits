"""
Applying D. Foreman-Mackey's emcee to orbit determination
"""

import argparse as ap
from astropy.table import Table

ECCENTRICITY = 0. # fix orbits as circular

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
    
    parser.add_argument('ephem_path',
                        help='path to file containing ephemeris info',
                        type=str)
    
    return parser.parse_args()

def formTLE():
	"""
	Form a TLE from the given parameters
	"""
	return None

def propagateTLE(tle, epoch):
	"""
	Propagate a TLE to the given epoch
	
	Parameters
	----------
	tle : TLE object
	    Object containing two line ephemeris information
	epoch : datetime object
	    datetime object [utc] with replaced tzinfo as utc
	"""
	return tle.radec(epoch)

def lnprior(theta):
	"""
	"""
	return lnprior

def lnlike():
	"""
	"""
	
	return lnlike

if __name__ == "__main__":
	
	args = argParse()
	
	
