"""
Gauss method of preliminary orbit determination
"""

import argparse as ap
import numpy as np
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import Latitude, Longitude

G = 6.6740831e-11     # SI units
M_EARTH = 5.9722e+24  # Mass of Earth [kg]
R_EARTH = 6378000.0   # Radius of Earth [m]
F = 0.003353          # flattening for Earth

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
    
    parser.add_argument('loc_path',
                        help='path to location config file',
                        type=str)
    
    parser.add_argument('--diagnostics',
                        help='include sanity checks?',
                        action='store_true')
    
    return parser.parse_args()

def parseInput(args):
	"""
	Obtain necessary data/ info from user input
	
	Parameters
	----------
	args : argParse output
	    Arguments provided by user
	
	Returns
	-------
	ephem_tab : astropy Table object
	    Table of ephemeris data
	latitude : astropy Latitude object
	    Latitude of observer
	altitude : float
	    Altitude of observer [m]
	"""
	# load ephemeris data
	try:
		ephem_tab = Table.read(args.ephem_path, format='csv')
	except FileNotFoundError:
		print('No ephemeris file found...')
		quit()
	
	# parse location config info
	try:
		loc_tab = Table.read(args.loc_path, format='csv')
	except FileNotFoundError:
		print('No location config file found...')
		quit()
	
	latitude = Latitude(loc_tab['latitude'], u.deg)
	altitude = loc_tab['altitude']
	
	return ephem_tab, latitude, altitude

def selectObs(ephem_tab, n1=None, n2=None, n3=None):
	"""
	Select the desired three observations to use in the analysis
	
	Parameters
	----------
	ephem_tab : astropy Table object
	    Table of ephemeris data
	n1, n2, n3 : int, optional
	    Indices specifying desired observations to use - if not given,
	    start-middle-end will be used
	    Default = None
	
	Returns
	-------
	obs1, obs2, obs3 : astropy Table objects
	    Three observations to take forward, updated with correct format
	"""
	if not n1 < n2 < n3:
		print('Observation indices must be in ascending order...')
		quit()
	
	if n1 is None or n2 is None or n3 is None:
		n1 = 0
		n3 = len(ephem) - 1
		n2 = round((n3 - n1) / 2)
	
	obs1 = ephem_tab[n1]
	obs2 = ephem_tab[n2]
	obs3 = ephem_tab[n3]
	
	return obs1, obs2, obs3

def positionVector(phi, lst, h):
	"""
	Position vector of an observer at a given time
	
	Parameters
	----------
	phi : float
	    Geodetic latitude of observer's location - angle between
	    the equatorial and normal planes [rad]
	lst : float
	    Local sidereal time for observation [rad]
	h : float
	    Altitude of observer [m]
	
	Returns
	-------
	r : array-like
	    Position vector of observer for given time
	"""
	r_x = (R_EARTH / np.sqrt(1 - (2*F - F**2)*np.sin(phi)**2) + h) \
	      *np.cos(phi)*np.cos(lst)
	r_y = (R_EARTH / np.sqrt(1 - (2*F - F**2)*np.sin(phi)**2) + h) \
	      *np.cos(phi)*np.sin(lst)
	r_z = (R_EARTH*(1 - F)**2 / np.sqrt(1 - (2*F - F**2)*np.sin(phi)**2) + h) \
	      *np.sin(phi)
	
	return np.array([r_x, r_y, r_z])

def cosineVector(alpha, delta):
	"""
	Direction cosine vector for an orbiting body
	
	Parameters
	----------
	alpha, delta : float
	    Topocentric right ascension & declination [rad]
	"""
	rho_hat_x = np.cos(delta)*np.cos(alpha)
	rho_hat_y = np.cos(delta)*np.sin(alpha)
	rho_hat_z = np.sin(delta)
	
	return np.array([rho_hat_x, rho_hat_y, rho_hat_z])

if __name__ == "__main__":
	
	args = argParse()
	
	ephem, phi, h = parseInput(args) # from user input
	
	# select 
