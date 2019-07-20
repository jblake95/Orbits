"""
Script for generating tlemcee config file
"""

import argparse as ap
import json
import pprint

## TODO: Generalise for more orbital types
TYPES = ['geo']
POS_CHECK = ['y', 'yes']
NEG_CHECK = ['n', 'no']

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
    
    parser.add_argument('--n_steps',
                        help='number of steps, default=1000',
                        type=int)
    
    parser.add_argument('--walker_scaling',
                        help='scaling factor for n_walkers, default=1',
                        type=int)
    
    parser.add_argument('--thinning_factor',
                        help='thinning factor for MCMC, default=1',
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
	
	defaults = 0
	if args.n_steps:
		config.update({'n_steps':args.n_steps})
	else:
		# resort to default=1000
		config.update({'n_steps':1000})
		print('No n_steps provided. Using default=1000')
		defaults += 1
	if args.walker_scaling:
		config.update({'walker_scaling':args.walker_scaling})
	else:
		# resort to default=1
		config.update({'walker_scaling':1})
		print('No walker_scaling provided. Using default=1')
		defaults += 1
	if args.thinning_factor:
		config.update({'thinning_factor':args.thinning_factor})
	else:
		# resort to default=1
		config.update({'thinning_factor':1})
		print('No thinning_factor provided. Using default=1')
		defaults += 1
	if defaults > 0:
		# check if defaults are desired
		while True:
			check = input('Accept the above defaults? (y|n) ')
			if check in POS_CHECK:
				print('Keeping defaults...')
				break
			elif check in NEG_CHECK:
				print('Quitting...')
				quit()
			else:
				print('Please specify "y" to keep or "n" to quit...')
	
	if args.orb_type in TYPES:
		if args.orb_type.lower() == 'geo':
			from params import Geo
			g = Geo()
			config.update(g.__dict__)	
	else:
		print('Please provide a valid orbital type:\n'
		      'Circular geosynchronous - "geo"')
		quit()
	
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
