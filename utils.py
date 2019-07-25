"""
Useful functions for the Orbits module
"""

import math

def fractional_yday(epoch):
    """
    Convert a datetime object to fractional day of the year
    
    Parameters
    ----------
    epoch : datetime object
        Epoch to convert to year day
    
    Returns
    -------
    yday : float
        Fractional year day corresponding to input epoch
    """
    frac = epoch.hour / 24. + epoch.minute / 60. + epoch.second / 60.
    
    return epoch.timetuple().tm_yday + frac

def n_digits(integer):
	"""
	Determine the number of digits in a given integer
	
	Parameters
	----------
	integer : int
	    Input integer
	
	Returns
	-------
	n_digits : int
	    Number of digits in input integer
	"""
	return len(str(abs(int(integer))))

def tle_standard_form(number):
	"""
	Convert a number to a scientific notation string consistent 
	with the NORAD Two Line Element data format
	
	Parameters
	----------
	number : float
	    Number to convert to TLE-friendly scientific notation
	
	Returns
	-------
	sci_not : str
	    Correctly formatted string giving TLE-friendly scientific 
	    notation for input number
	"""
	sci_not = '{:.4e}'.format(number)
	e_pos = sci_not.find('e')
	
	mantissa = float(sci_not[:e_pos])
	if mantissa >= 0:
		mantissa_sign = '+'
	else:
		mantissa_sign = '-'
	mantissa_str = '{:.5f}'.format(mantissa / 10.).lstrip('+-0.')
	
	exponent = int(sci_not[e_pos+1:]) + 1
	if exponent >= 0:
		exponent_sign = '+'
	else:
		exponent_sign = '-'
	exponent_str = str(exponent).lstrip('+-')
	
	return mantissa_sign + mantissa_str + exponent_sign + exponent_str
	
