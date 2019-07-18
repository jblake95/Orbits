"""
Script for generating tlemcee config file
"""

import argparse as ap
import json
import pprint

## TODO: Generalise for more orbital types
TYPES = ['geo']

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
    
    parser.add_argument('out_dir',
                        help='output directory for config json file',
                        type=str)
    
    parser.add_argument('orb_type',
                        help='desired orbit type\n'
                             'Circular geosyncronous - "geo"',
                        type=str)
    
    parser.add_argument('n_steps',
                        help='Number of steps',
                        type=int)
    
    parser.add_argument('n_walkers',
                        help='Number of walkers',
                        type=int)
    
    parser.add_argument('--verbose',
                        help='print result?',
                        action='store_true')
    
    return parser.parse_args()

def generateConfig(args):
	"""
	Generate config dictionary
	"""
	config = {}
	
	config.update({'ephem_path':args.ephem_path})
	config.update({'out_dir':args.out_dir})
	config.update({'n_steps':args.n_steps})
	config.update({'n_walkers':args.n_walkers})
	
	if args.orb_type in TYPES:
		fixed = {}
		uniform = {}
		if args.orb_type.lower() == 'geo':
			from params import Geo
			g = Geo()
			# fixed params
			ecc = {}
			ecc.update({'value':g.ecc_val})
			fixed.update({'ecc':ecc})
			
			# uniform prior params
			argp = {}
			argp.update({'init':g.argp_init})
			argp.update({'llim':g.argp_llim})
			argp.update({'ulim':g.argp_ulim})
			argp.update({'wght':g.argp_wght})
			uniform.update({'argp':argp})
			
			raan = {}
			raan.update({'init':g.raan_init})
			raan.update({'llim':g.raan_llim})
			raan.update({'ulim':g.raan_ulim})
			raan.update({'wght':g.raan_wght})
			uniform.update({'raan':raan})
			
			mmdot = {}
			mmdot.update({'init':g.mmdot_init})
			mmdot.update({'llim':g.mmdot_llim})
			mmdot.update({'ulim':g.mmdot_ulim})
			mmdot.update({'wght':g.mmdot_wght})
			uniform.update({'mmdot':mmdot})
			
			mm = {}
			mm.update({'init':g.mm_init})
			mm.update({'llim':g.mm_llim})
			mm.update({'ulim':g.mm_ulim})
			mm.update({'wght':g.mm_wght})
			uniform.update({'mm':mm})
			
			incl = {}
			incl.update({'init':g.incl_init})
			incl.update({'llim':g.incl_llim})
			incl.update({'ulim':g.incl_ulim})
			incl.update({'wght':g.incl_wght})
			uniform.update({'incl':incl})
			
			anom = {}
			anom.update({'init':g.anom_init})
			anom.update({'llim':g.anom_llim})
			anom.update({'ulim':g.anom_ulim})
			anom.update({'wght':g.anom_wght})
			uniform.update({'anom':anom})	
	else:
		print('Please provide a valid orbital type:\n'
		      'Circular geosynchronous - "geo"')
		quit()
	
	config.update({'fixed':fixed})
	config.update({'uniform':uniform})
	
	if args.verbose:
		print('Config:\n'
		      '------')
		pp = pprint.PrettyPrinter(indent=4)
		pp.pprint(config)
	
	# save config to file
	with open('{}/config.json'.format(args.out_dir), 'w') as f:
		json.dump(config, f)
	
	return config

if __name__ == "__main__":
	
	args = argParse()
	
	config = generateConfig(args)
