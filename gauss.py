"""
Gauss method of preliminary orbit determination
"""

import argparse as ap
import numpy as np
from astropy.table import Table
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
	ephem : astropy Table object
	    Table of ephemeris data
	phi : astropy Latitude object
	    Latitude of observer
	h : float
	    Altitude of observer
	"""
	# load ephemeris data
	try:
		ephem = Table.read(args.ephem_path, format='csv')
	except FileNotFoundError:
		print('No ephemeris file found...')
		quit()
	
	# parse config info
	try:
		loc = Table.read(args.loc_path, format='csv')
	except FileNotFoundError:
		print('No location config file found...')
		quit()
	
	return ephem, Latitude(loc['latitude'], u.deg), loc['altitude']

def positionVector(phi, lst, h):
	"""
	Determine the position vector of an observer at a given time
	
	Parameters
	----------
	phi : float
	    Geodetic latitude of observer's location - angle between
	    the equatorial and normal planes
	lst : float
	    Local sidereal time for observation [deg]
	h : float
	    Altitude of observer [m]
	
	Returns
	-------
	r : array-like
	    Position vector of observer for given time
	"""
	r_x = (R_EARTH / np.sqrt(1 - (2*F - F**2)*np.sin(phi)**2) + h)*np.cos(phi)*np.cos(lst)
	r_y = (R_EARTH / np.sqrt(1 - (2*F - F**2)*np.sin(phi)**2) + h)*np.cos(phi)*np.sin(lst)
	r_z = (R_EARTH*(1 - F)**2 / np.sqrt(1 - (2*F - F**2)*np.sin(phi)**2) + h)*np.sin(phi)
	
	return np.array([r_x, r_y, r_z])

if __name__ == "__main__":
	
	r = positionVector(40, 44.506, 1000)
	print(r)
