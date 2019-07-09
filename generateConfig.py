"""
Script for generating tlemcee config file
"""

from params import Geo
import argparse as ap
import json

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
    
    parser.add_argument('--verbose',
                        help='print result?',
                        action='store_true')
    
    return parser.parse_args()

if __name__ == "__main__":
	
	args = argParse()
	
	
