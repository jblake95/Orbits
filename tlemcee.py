"""
Applying D. Foreman-Mackey's emcee to orbit determination
"""

import argparse as ap
from astropy.table import Table

"""
NOTES
-----
Things we'll need to feed to the tle:
(ranges based on geo_investigator snippet - see snippets)

(1) Line 1:
epoch in form yyddd.dddddddd where d's are fractional year day
-0.000004<mmdot<0.000002
NB: this appears in tle divided by 2
NB: mmdot2 and drag term seemingly do not affect ephemeris calc 

(2) Line 2
0<incl<60
0<raan<360
Fix e=0000000 to enforce circular orbit
0<argp<360
0<anom<360
0.99<mm<1.01 
"""

# TLE templates - non-essential elements set to zero
TLE1 = '1 00000U 000000   {:2}{:12.8f} {:10}  00000-0  00000+0 0 00000'
TLE2 = '2 00000 {:8.4f} {:8.4f} 0000000 {:8.4f} {:8.4f} {:11.8f}000000'

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

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
    
    parser.add_argument('config_path',
                        help='path to file containing config info',
                        type=str)
    
    return parser.parse_args()

def readConfig(config_path):
	"""
	Extract info from config file
	
	Parameters
	----------
	config_path : str
	    Path to file containing config info
	
	Returns
	-------
	config : astropy Table object
	    Table containing config info
	"""
	try:
		config = Table.read(config_path, format='csv')
	except FileNotFoundError:
		print('Config file not found...')
		quit()
	
	return config

def formTLE(year, yday, mmdot, mmdot2, incl, raan, ):
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
	Log prior function - used to ensure walkers are exploring the range 
	defined by the priors
	
	Parameters
	----------
	theta : array-like
	    current set of parameters from MCMC
	
	lnprior : float
        log prior of current sample (0.0 | -np.inf)
	"""
	
	
	return lnprior

def lnlike(theta, t, ra, raerr, dec, decerr):
	"""
	Log likelihood for a given subset of data
	"""
	
	##TODO: Will need to split process into ra over time and dec over
	##      time by the looks of it, then recombine here
	
	return lnlike

if __name__ == "__main__":
	
	args = argParse()
	
	