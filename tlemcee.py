"""
Applying D. Foreman-Mackey's emcee to orbit determination
"""

from tle import DummyTLE
import json
import numpy as np
import argparse as ap
from collections import OrderedDict
from astropy.table import Table
import emcee
import corner
import matplotlib.pyplot as plt

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
TLE2 = '2 00000 {:8.4f} {:8.4f} {:07d} {:8.4f} {:8.4f} {:11.8f}000000'

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
    
    parser.add_argument('directory',
                        help='directory containing config file',
                        type=str)
    
    parser.add_argument('--threads',
                        help='number of threads to run',
                        type=int,
                        default=1)
    
    return parser.parse_args()

def readConfig(directory):
	"""
	Extract info from config file
	
	Parameters
	----------
	directory : str
	    Directory containing config file
	
	Returns
	-------
	config : astropy Table object
	    Table containing config info
	"""	
	# load config file
	try:
		with open('{}/config.json'.format(directory), 'r') as cfg:
			config = json.load(cfg)
	except FileNotFoundError:
		print('Config file not found...')
		quit()
	
	# load ephemeris info
	if 'ephem_path' in config.keys():
		try:
			ephem = Table.read(config['ephem_path'], format='csv')
		except FileNotFoundError:
			print('Ephemeris file not found...')
			quit()
	else:
		print('Invalid config: please supply ephem_path...')
		quit()
	
	return OrderedDict(config), ephem

def formatTLE(year, yday, mmdot, mmdot2, incl, raan, ):
	"""
	Form a TLE from the given parameters
	"""
	tle1 = TLE1.format()
	
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
	mmdot, incl, raan, argp, anom, mm = theta
	
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
	
	# read config info
	config, ephem = readConfig(args.directory)
	
	# set up initial seed list and weights
	initial = [config['uniform'][key]['init'] for key in config['uniform'].keys()]
	weights = [config['uniform'][key]['wght'] for key in config['uniform'].keys()]
	n_priors = len(initial)
	
	# obtain data arrays
	x_t = np.array(ephem['utc'])
	y_ra = np.array(ephem['ra'])
	yerr_ra = np.array(ephem['raerr'])
	y_dec = np.array(ephem['dec'])
	yerr_dec = np.array(ephem['decerr'])
	
	####################################################################
	###################### set up the sampler ##########################
	####################################################################
	n_dim = len(initial)
	# recommended nwalkers is 4*n_parameters
    # more walkers can help find the global minima, so optional scaling
    n_walkers = 4*len(initial) * config['walker_scaling']
	
    # set up the starting positions
    pos = [initial 
           + weights*np.random.randn(n_dim) for i in range(n_walkers)]
    
    sampler = emcee.EnsembleSampler(n_walkers, n_dim, lnprob,
                                    args=(config, n_priors,
                                          x_lc, y_lc, yerr_lc,
                                          x_rv, y_rv, yerr_rv),
                                    threads=args.threads)
	
