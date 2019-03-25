"""
Gauss method of preliminary orbit determination
"""

import argparse as ap
import numpy as np
from datetime import datetime, timedelta
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import Latitude, Longitude

G = 6.6740831e-11     # SI units
M_EARTH = 5.9722e+24  # Mass of Earth [kg]
R_EARTH = 6378000.    # Radius of Earth [m]
F = 0.003353          # flattening for Earth
MU = 398600.          # gravitational parameter [km^3/s^2]

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

class Observation:
	"""
	Observation of an orbiting body at a given time
	"""
	def __init__(self, obs_tab):
		"""
		Initiate Observation object
		"""
		self.ra = Longitude(obs_tab['ra'], u.deg)
		self.dec = Latitude(obs_tab['dec'], u.deg)
		self.raerr = Longitude(obs_tab['raerr'], u.deg)
		self.decerr = Latitude(obs_tab['decerr'], u.deg)
		self.utc = datetime.strptime(obs_tab['utc'], 
		                             '%Y-%m-%dT%H:%M:%S.%f')

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
	if n1 is None or n2 is None or n3 is None:
		n1 = 0
		n3 = len(ephem_tab) - 1
		n2 = round((n3 - n1) / 2)
	else:
		if not n1 < n2 < n3:
			print('Observation indices must be in ascending order...')
			quit()
	
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
	
	Returns
	-------
	rho_hat : array-like
	    Direction cosine vector for the orbiting body
	"""
	rho_hat_x = np.cos(delta)*np.cos(alpha)
	rho_hat_y = np.cos(delta)*np.sin(alpha)
	rho_hat_z = np.sin(delta)
	
	return np.array([rho_hat_x, rho_hat_y, rho_hat_z])

def performAlgorithm(args, obs_idx=[None, None, None]):
	"""
	Carry out the Gauss method of preliminary orbit determination
	
	Parameters
	----------
	args : argParse output
	    Arguments provided by user 
	obs_idx : array-like, optional
	    Indices corresponding to observations to feed into algorithm
	    (format: [n1, n2, n3] where ni specifies index i)
	    default = None, will automatically select start-middle-end
	
	Returns
	-------
	vel_vec : array-like
	    Velocity vectors for the three observations
	"""
	ephem, phi, h = parseInput(args) # from user input
	
	# select which observations to use
	obs1, obs2, obs3 = selectObs(ephem,
	                             n1=obs_idx[0],
	                             n2=obs_idx[1],
	                             n3=obs_idx[2])
	obs1 = Observation(obs1)
	obs2 = Observation(obs2)
	obs3 = Observation(obs3)
	
	# position and cosine vectors
	## TODO: get lst from frame headers
	## for now use example from Curtis textbook
	r_1 = np.array([3489.8, 3430.2, 4078.5])
	r_2 = np.array([3460.1, 3460.1, 4078.5])
	r_3 = np.array([3429.9, 3490.1, 4078.5])
	
	rho_hat_1 = np.array([0.71643, 0.68074, -0.15270])
	rho_hat_2 = np.array([0.56897, 0.79531, -0.20917])
	rho_hat_3 = np.array([0.41841, 0.87007, -0.26059])
	
	# step 1 - time intervals
	tau_1 = -118.10
	tau_3 = 119.47
	tau = 237.58
	
	# step 2 - rho_hat cross products
	p_1 = np.cross(rho_hat_2, rho_hat_3)
	p_2 = np.cross(rho_hat_1, rho_hat_3)
	p_3 = np.cross(rho_hat_1, rho_hat_2)
	
	# step 3 - rho_hat scalar triple product
	d_0 = np.dot(rho_hat_1, p_1)
	
	# step 4 - compute scalar quantities
	d_11 = np.dot(r_1, p_1)
	d_12 = np.dot(r_1, p_2)
	d_13 = np.dot(r_1, p_3)
	d_21 = np.dot(r_2, p_1)
	d_22 = np.dot(r_2, p_2)
	d_23 = np.dot(r_2, p_3)
	d_31 = np.dot(r_3, p_1)
	d_32 = np.dot(r_3, p_2)
	d_33 = np.dot(r_3, p_3)
	
	# step 5 - calculate scalar position coefficients
	A = (1 / d_0)*(-d_12*(tau_3 / tau) + d_22 + d_32*(tau_1 / tau))
	B = (1 / (6*d_0))*(d_12*(tau_3**2 - tau**2)*tau_3 / tau +
	                   d_32*(tau**2 - tau_1**2)*tau_1 / tau)
	E = np.dot(r_2, rho_hat_2)
	
	# step 6 - squared scalar distance of obs2
	r2_2 = np.dot(r_2, r_2)
	
	# step 7 - coefficients for scalar distance polynomial
	a = -(A**2 + 2*A*E + r2_2)
	b = -2*MU*B*(A + E)
	c = -(MU**2)*(B**2)
	
	print(a, b, c, MU)
	
	return np.array([])

if __name__ == "__main__":
	
	args = argParse()
	
	vel = performAlgorithm(args)
	
	
