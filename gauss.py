"""
Gauss method of preliminary orbit determination
"""

import argparse as ap
import numpy as np
from datetime import datetime, timedelta
from scipy.optimize import newton
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import Latitude, Longitude, Angle
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector

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

def toggleSelector(event):
    """
    Toggle function compatible with the RectangleSelector widget 
    in matplotlib.widgets
    
    Parameters
    ----------
    event : Mouse event
        Mouse event upon matplotlib artist
    """
    print('Key pressed.')
    if event.key in ['Q', 'q'] and toggleSelector.RS.active:
        print('RectangleSelector deactivated.')
        toggleSelector.RS.set_active(False)
    if event.key in ['A', 'a'] and not toggleSelector.RS.active:
        print('RectangleSelector activated.')
        toggleSelector.RS.set_active(True)

def onSelect(eclick, erelease):
    """
    Callback function compatible with the RectangleSelector widget
    in matplotlib.widgets
    
    Parameters
    ----------
    eclick, erelease : Mouse events
        Press and release events corresponding to the placement of
        a rectangle 
    """
    global zero
    zero = (eclick.xdata + erelease.xdata) / 2
    
    print('Selected: x = {}'.format(str(zero)))
    print('The button you used was: {}'.format(str(eclick.button))) 

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

def poly8(x, a, b, c):
    """
    Form of the eighth degree polynomial for use in algorithm 
    
    Parameters
    ----------
    x : float
        Input value
    a, b, c : float
        Polynomial coefficients
    
    Returns
    -------
    y : float
        Output value
    """
    return x**8 + a*x**6 + b*x**3 + c

def poly8prime(x, a, b, c):
    """
    Derivative of the eight degree polynomial for use in algorithm
    
    Parameters
    ----------
    x : float
        Input value
    a, b, c : float
        Polynomial coefficients
    
    Returns
    -------
    y : float
        Output value
    """
    return 8*x**7 + 6*a*x**5 + 3*b*x**2

def orbElementsAlgorithm(r_vec, v_vec):
    """
    Obtain the orbital elements from the state vector
    
    Parameters
    ----------
    r_vec : array-like
        Position vector for the orbiting body [km]
    v_vec : array-like
        Velocity vector for the orbiting body [km/s]
    
    Returns
    -------
    orb_elements : array-like
        Orbital elements for the orbiting body in following order
        a     - semimajor axis [km]
        e     - eccentricity
        i     - inclination [deg]
        Omega - right ascension of the ascending node [deg]
        omega - argument of perigee [deg]
        theta - true anomaly [deg]
    """
    # step 1 - distance
    r = np.sqrt(np.dot(r_vec, r_vec))
    
    # step 2 - speed
    v = np.sqrt(np.dot(v_vec, v_vec))
    
    # step 3 - radial velocity
    v_r = np.dot(v_vec, r_vec) / r
    
    # step 4 - specific angular momentum
    h_vec = np.cross(r_vec, v_vec)
    
    # step 5 - magnitude of specific angular momentum
    h = np.sqrt(np.dot(h_vec, h_vec))
    
    # step 6 - inclination [rad]
    i = Angle(np.arccos(h_vec[2] / h), 
              u.rad)
    
    # step 7 - node line vector
    k_hat = (0, 0, 1)
    N_vec = np.cross(k_hat, h_vec)
    
    # step 8 - magnitude of node line vector
    N = np.sqrt(np.dot(N_vec, N_vec))
    
    # step 9 - right ascension of the ascending node [rad]
    if N_vec[1] >= 0:
        Omega = Angle(np.arccos(N_vec[0] / N), 
                      u.rad)
    else:
        Omega = Angle(2*np.pi - np.arccos(N_vec[0] / N), 
                      u.rad)
    
    # step 10 - eccentricity vector
    e_vec = (1 / MU)*((v**2 - (MU / r))*r_vec - r*v_r*v_vec)
    
    # step 11 - eccentricity
    e = np.sqrt(np.dot(e_vec, e_vec))
    
    # step 12 - argument of perigee
    if e_vec[2] >= 0:
        omega = Angle(np.arccos(np.dot(N_vec, e_vec) / (N*e)), 
                      u.rad)
    else:
        omega = Angle(2*np.pi - np.arccos(np.dot(N_vec, e_vec) / (N*e)), 
                      u.rad)
    
    # step 13 - true anomaly
    if v_r >= 0:
        theta = Angle(np.arccos(np.dot(e_vec, r_vec) / (e*r)), 
                      u.rad)
    else:
        theta = Angle(2*np.pi - np.arccos(np.dot(e_vec, r_vec) / (e*r)), 
                      u.rad)
    
    # additional extras - perigee & apogee radii, semimajor axis, period
    r_p = (h**2 / MU)*(1 / (1 + e*np.cos(0.)))
    r_a = (h**2 / MU)*(1 / (1 + e*np.cos(np.pi)))
    
    a = (1 / 2)*(r_p + r_a)
    T = (2*np.pi / np.sqrt(MU))*a**(3 / 2)
    
    print('--------------------\n'
          'Orbital information:\n'
          '--------------------\n'
          'a[km]      = {}\n'
          'e          = {}\n'
          'i[deg]     = {}\n'
          'Omega[deg] = {}\n'
          'omega[deg] = {}\n'
          'theta[deg] = {}\n'
          'T[hrs]     = {}\n'
          '--------------------'.format(str(a),
                                        str(e),
                                        str(i.deg),
                                        str(Omega.deg),
                                        str(omega.deg),
                                        str(theta.deg),
                                        str(T / 3600.)))
    
    return np.array([a, e, i.deg, Omega.deg, omega.deg, theta.deg])

def stumpffC(z):
    """
    Form of the Stumpff function C(z) in terms of circular and 
    hyperbolic trig functions
    
    Parameters
    ----------
    z : float
        Input value, calculated from the reciprocal semimajor axis and
        universal anomaly as z_i = alpha*chi_i**2
    
    Returns
    -------
    c : float
        Output value
    """
    if z > 0:
        c = (1 - np.cos(np.sqrt(z))) / z
    elif z < 0:
        c = (np.cosh(np.sqrt(-z)) - 1) / -z
    else:
        c = 1 / 2
    
    return c

def stumpffS(z):
    """
    Form of the Stumpff function S(z) in terms of circular and 
    hyperbolic trig functions
    
    Parameters
    ----------
    z : float
        Input value, calculated from the reciprocal semimajor axis and
        universal anomaly as z_i = alpha*chi_i**2
    
    Returns
    -------
    s : float
        Output value
    """
    if z > 0:
        s = (np.sqrt(z) - np.sin(np.sqrt(z))) / np.sqrt(z)**3
    elif z < 0:
        s = (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / np.sqrt(-z)**3
    else:
        s = 1 / 6
    
    return s

def solveUniversalKepler(delta_t, r_0, v_r0, alpha, tolerance=1e-6):
    """
    Solve the universal Kepler's equation for the universal anomaly
    
    Parameters
    ----------
    delta_t : float
        Time interval
    r_0, v_r0 : float
        Orbital radius and radial velocity at time t = t_0
    alpha : float
        Reciprocal of semimajor axis
    tolerance : float
        Error tolerance defining algorithmic success
    
    Returns
    -------
    chi : float
        Universal anomaly
    """
    # step 1 - initial estimate of chi_i
    chi = np.sqrt(MU)*abs(alpha)*delta_t
    
    count = 0
    ratio = np.inf
    while abs(ratio) > tolerance:
        if count > 10:
            print('Looped 10 times. Try increasing tolerance...')
            quit()
        
        # step 4 - if ratio (see below) exceeds tolerance, update chi
        if ratio is not np.inf:
            chi -= ratio
        
        # step 2 - calculate f(chi) & f'(chi)
        z = alpha*chi**2
        f_chi = ((r_0*v_r0 / np.sqrt(MU))*chi**2*stumpffC(z) +
                (1 - alpha*r_0)*chi**3*stumpffS(z) + r_0*chi -
                np.sqrt(MU)*delta_t)
        f_chi_prime = ((r_0*v_r0 / np.sqrt(MU))*chi*(1 - z*stumpffS(z)) +
                      (1 - alpha*r_0)*chi**2*stumpffC(z) + r_0)
        
        # step 3 - calculate ratio
        ratio = f_chi / f_chi_prime
        
        count += 1
    
    # step 5 - if ratio is less than tolerance, accept chi as solution
    return chi

def round_sig(x, sig=5):
    """
    Round a given number to the requested number of significant figures
    
    Parameters
    ----------
    x : float
        Number to be rounded
    sig : int
        Number of significant figures
    
    Returns
    -------
    y : float
        Rounded number
    """
    if x == 0.:
        return 0.
    else:
        return round(x, sig - int(np.floor(np.log10(abs(x)))) - 1)

def improveOrbit(r_vec, v_vec, tau_1, tau_3, p, tolerance=1e-6, sig=5):
    """
    Iterative improvement of the orbit determined with Gauss method
    
    Parameters
    ----------
    r_vec, v_vec : array-like
        The state vectors obtained using the Gauss method
    tau_1, tau_3 : float
        Time intervals from initial stages of the Gauss method
    p : dict
        Dictionary of parameters carried forward from the Gauss method
    tolerance : float, optional
        Error tolerance for universal Kepler equation solver
        Default = 1e-6
    sig : int, optional
        Number of significant figures for improvement algorithm success
        Default = 5
    
    Returns
    -------
    orb_elements : array-like
        Orbital elements for the orbiting body in following order
        a     - semimajor axis [km]
        e     - eccentricity
        i     - inclination [deg]
        Omega - right ascension of the ascending node [deg]
        omega - argument of perigee [deg]
        theta - true anomaly [deg]
        improved using 'exact' values of the Lagrange coefficients
    """
    print('r_vec: ', r_vec)
    print('v_vec: ', v_vec)
    print('tau_1: {} tau_3: {}'.format(str(tau_1),
                                       str(tau_3)))
    i = 0
    rho_1_list = [round_sig(p['rho_1'], sig=sig)]
    rho_2_list = [round_sig(p['rho_2'], sig=sig)]
    rho_3_list = [round_sig(p['rho_3'], sig=sig)]
    while True:
        # step 1 - distance and speed
        r = np.sqrt(np.dot(r_vec, r_vec))
        v = np.sqrt(np.dot(v_vec, v_vec))
        print('r: {} v: {}'.format(str(r),
                                   str(v)))
        # step 2 - reciprocal of semimajor axis
        alpha = 2 / r - v**2 / MU
        print('alpha: {}'.format(str(alpha)))
        # step 3 - radial component of v_vec
        v_r = np.dot(v_vec, r_vec) / r
        print('v_r: {}'.format(str(v_r)))
        # step 4 - solve universal Kepler's equation for chi_i
        chi_1 = solveUniversalKepler(tau_1, 
                                     r, 
                                     v_r, 
                                     alpha,
                                     tolerance=tolerance)
        chi_3 = solveUniversalKepler(tau_3, 
                                     r, 
                                     v_r, 
                                     alpha,
                                     tolerance=tolerance)
        print('chi_1: {} chi_3: {}'.format(str(chi_1),
                                           str(chi_3)))
        # step 5 - use chi_i to determine new Lagrange coefficients
        z_1 = alpha*chi_1**2
        z_3 = alpha*chi_3**2
        print('z_1: {} z_3: {}'.format(str(z_1),
                                       str(z_3)))
        f_1 = 1 - (chi_1**2 / r)*stumpffC(z_1)
        g_1 = tau_1 - (1 / np.sqrt(MU))*chi_1**3*stumpffS(z_1)
        f_3 = 1 - (chi_3**2 / r)*stumpffC(z_3)
        g_3 = tau_3 - (1 / np.sqrt(MU))*chi_3**3*stumpffS(z_3)
        print('f_1: {} f_3: {}'.format(str(f_1),
                                       str(f_3)))
        print('g_1: {} g_3: {}'.format(str(g_1),
                                       str(g_3)))
        # step 6 - determine c_i
        c_1 = g_3 / (f_1*g_3 - f_3*g_1)
        c_3 = -g_1 / (f_1*g_3 - f_3*g_1)
        print('c_1: {} c_3: {}'.format(str(c_1),
                                       str(c_3)))
        # step 7 - update values of rho_i
        rho_1 = (1 / p['d_0'])*(-p['d_11'] + (1 / c_1)*p['d_21'] -
                                (c_3 / c_1)*p['d_31'])
        rho_2 = (1 / p['d_0'])*(-c_1*p['d_12'] + p['d_22'] - 
                                c_3*p['d_32'])
        rho_3 = (1 / p['d_0'])*(-(c_1 / c_3)*p['d_13'] + 
                                (1 / c_3)*p['d_23'] - p['d_33'])
        print('rho_1: {} rho_2: {} rho_3: {}'.format(str(rho_1),
                                                     str(rho_2),
                                                     str(rho_3)))
        rho_1_list.append(round_sig(rho_1, sig=sig))
        rho_2_list.append(round_sig(rho_2, sig=sig))
        rho_3_list.append(round_sig(rho_3, sig=sig))
        
        # step 8 - update geocentric position vectors r_i
        r_1 = p['R_1'] + rho_1*p['rho_hat_1']
        r_2 = p['R_2'] + rho_2*p['rho_hat_2']
        r_3 = p['R_3'] + rho_3*p['rho_hat_3']
        print('r_1: ', r_1)
        print('r_2: ', r_2)
        print('r_3: ', r_3)
        # step 9 - update velocity vector v
        v_2 = (1 / (f_1*g_3 - f_3*g_1))*(-f_3*r_1 + f_1*r_3)
        print('v_2: ', v_2)
        # step 10 - check if rho_i are the 'same' to within sig figs
        if (rho_1_list[i] == rho_1_list[i-1] and 
            rho_2_list[i] == rho_2_list[i-1] and
            rho_3_list[i] == rho_3_list[i-1]):
            break
        elif i > 100:
            print('Failing to converge: try a more generous sig...')
            quit()
        else:
            r_vec, v_vec = r_2, v_2 # update state vector and repeat
        
        """
        print('Step {}: {} {} {}'.format(str(i + 1),
                                         str(rho_1),
                                         str(rho_2),
                                         str(rho_3)))
        """
        print('r_vec: ', r_vec)
        print('v_vec: ', v_vec)
        i += 1
        input('enter.')
    
    orb_elements = 0 # dummy 
    return orb_elements

def gaussAlgorithm(args, obs_idx=[None, None, None], improve=False):
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
    improve : bool, optional
        Toggle to perform iterative improvement of the determined orbit
        default = False
    
    Returns
    -------
    vel_vec : array-like
        Velocity vectors for the three observations
    """
    global zero
    
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
    ## for now use example from Curtis textbook, algorithm 5.5
    R_1 = np.array([3489.8, 3430.2, 4078.5])
    R_2 = np.array([3460.1, 3460.1, 4078.5])
    R_3 = np.array([3429.9, 3490.1, 4078.5])
    
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
    d_11 = np.dot(R_1, p_1)
    d_12 = np.dot(R_1, p_2)
    d_13 = np.dot(R_1, p_3)
    d_21 = np.dot(R_2, p_1)
    d_22 = np.dot(R_2, p_2)
    d_23 = np.dot(R_2, p_3)
    d_31 = np.dot(R_3, p_1)
    d_32 = np.dot(R_3, p_2)
    d_33 = np.dot(R_3, p_3)
    
    # step 5 - calculate scalar position coefficients
    A = (1 / d_0)*(-d_12*(tau_3 / tau) + d_22 + d_32*(tau_1 / tau))
    B = (1 / (6*d_0))*(d_12*(tau_3**2 - tau**2)*tau_3 / tau +
                       d_32*(tau**2 - tau_1**2)*tau_1 / tau)
    E = np.dot(R_2, rho_hat_2)
    
    # step 6 - squared scalar distance for obs2
    R2_2 = np.dot(R_2, R_2)
    
    # step 7 - coefficients of scalar distance polynomial for obs2
    a = -(A**2 + 2*A*E + R2_2)
    b = -2*MU*B*(A + E)
    c = -(MU**2)*(B**2)
    
    # step 8 - find root of scalar distance polynomial for obs2
    roots = np.roots([1,0,a,0,0,b,0,0,c])
    physical_roots = []
    for root in roots:
        if not np.iscomplex(root):
            if root >= 0:
                physical_roots.append(np.real(root))
    
    if len(physical_roots) == 0:
        print('Oh dear, no physical roots!')
        quit()
    elif len(physical_roots) == 1:
        print('One physical root found, using this...')
        zero = physical_roots[0]
    else:
        print('Multiple physical roots found. Please choose one...')
        x_min = min(physical_roots) - 0.5*(max(physical_roots) -
                                           min(physical_roots))
        x_max = max(physical_roots) + 0.2*(max(physical_roots) -
                                           min(physical_roots))
        x_dummy = np.linspace(x_min, x_max, 10000)
    
        fig, ax = plt.subplots()
        
        plt.plot(x_dummy, poly8(x_dummy, a, b, c), 'k.', ms=1)
        plt.axhline(y=0, color='r', linestyle='--')
        
        plt.title('Please select a region close to a zero...')
        plt.xlabel('x')
        plt.ylabel('F(x)')
        
        toggleSelector.RS = RectangleSelector(ax, 
                                              onSelect,
                                              drawtype='box',
                                              interactive=True)
        
        plt.connect('key_press_event', toggleSelector)
        plt.show()
        
        zero = newton(poly8, zero, poly8prime, (a, b, c))
    
    print('Zero: {}'.format(str(zero))) # geocentric radius r_2
    
    # step 9 - obtain slant ranges
    num_1 = (6*(d_31*(tau_1 / tau_3) + d_21*(tau / tau_3))*zero**3 +
                           MU*d_31*(tau**2 - tau_1**2)*(tau_1 / tau_3))
    den_1 = 6*zero**3 + MU*(tau**2 - tau_3**2)
    num_3 = (6*(d_13*(tau_3 / tau_1) - d_23*(tau / tau_1))*zero**3 +
                           MU*d_13*(tau**2 - tau_3**2)*(tau_3 / tau_1))
    den_3 = 6*zero**3 + MU*(tau**2 - tau_1**2)
    
    rho_1 = (1 / d_0)*(num_1 / den_1 - d_11)
    rho_2 = A + MU*B / zero**3
    rho_3 = (1 / d_0)*(num_3 / den_3 - d_33)
    
    # step 10 - orbiting body geocentric position vectors
    r_1 = R_1 + rho_1*rho_hat_1
    r_2 = R_2 + rho_2*rho_hat_2
    r_3 = R_3 + rho_3*rho_hat_3
    
    # step 11 - Lagrange coefficients
    f_1 = 1 - (1 / 2)*(MU / zero**3)*tau_1**2
    f_3 = 1 - (1 / 2)*(MU / zero**3)*tau_3**2
    g_1 = tau_1 - (1 / 6)*(MU / zero**3)*tau_1**3
    g_3 = tau_3 - (1 / 6)*(MU / zero**3)*tau_3**3
    
    # step 12 - velocity vector for obs2
    v_2 = (1 / (f_1*g_3 - f_3*g_1))*(-f_3*r_1 + f_1*r_3)
    
    # step 13 - orbital elements
    orb_elements = orbElementsAlgorithm(r_2, v_2)
    
    # Improve the state vector 
    ## Using example from Curtis textbook Chapter 5.10
    r_2 = np.array([5659.1, 6533.8, 3270.1])
    v_2 = np.array([-3.880, 5.1156, -2.2387])
    tau_1 = -118.10
    tau_3 = 119.47
    rho_1 = 3639.1
    rho_2 = 3864.8
    rho_3 = 4172.8
    R_1 = np.array([3489.8, 3430.2, 4078.5])
    R_2 = np.array([3460.1, 3460.1, 4078.5])
    R_3 = np.array([3429.9, 3490.1, 4078.5])
    rho_hat_1 = np.array([0.71643, 0.68074, -0.15270])
    rho_hat_2 = np.array([0.56897, 0.79531, -0.20917])
    rho_hat_3 = np.array([0.41841, 0.87007, -0.26059])
    d_0 = -0.0015198
    d_11 = 782.15
    d_12 = 1646.5
    d_13 = 887.10
    d_21 = 784.72
    d_22 = 1651.5
    d_23 = 889.60
    d_31 = 787.31
    d_32 = 1656.6
    d_33 = 892.13
    
    params = {}
    params.update({'rho_1':rho_1, 'rho_2':rho_2, 'rho_3':rho_3,
                   'R_1':R_1, 'R_2':R_2, 'R_3':R_3,
                   'rho_hat_1':rho_hat_1, 
                   'rho_hat_2':rho_hat_2,
                   'rho_hat_3':rho_hat_3,
                   'd_0':d_0, 
                   'd_11':d_11, 'd_12':d_12, 'd_13':d_13,
                   'd_21':d_21, 'd_22':d_22, 'd_23':d_23,
                   'd_31':d_31, 'd_32':d_32, 'd_33':d_33})
    _ = improveOrbit(r_2, v_2, tau_1, tau_3, params)
    
    return orb_elements

if __name__ == "__main__":
    
    args = argParse()
    
    orb_elements = gaussAlgorithm(args)
    
    
