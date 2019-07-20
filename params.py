"""
Config params for different orbit types
"""

class Geo:
	"""
	Circular geosynchronous orbit
	"""
	def __init__(self):
		"""
		Initiate Geo object
		"""
		# fixed params
		self.fixed = {}
		self.fixed.update({'ecc':{'value':0.}})
		
		# uniform prior params
		self.uniform = {}
		self.uniform.update({'argp':{'init':180.,
		                             'llim':0.,
		                             'ulim':360.,
		                             'wght':0.5}})
		self.uniform.update({'raan':{'init':180.,
		                             'llim':0.,
		                             'ulim':360.,
		                             'wght':0.5}})
		self.uniform.update({'mmdot':{'init':0.,
		                              'llim':-0.000004,
		                              'ulim':0.000002,
		                              'wght':0.00000001}})
		self.uniform.update({'mm':{'init':1.,
		                           'llim':0.99,
		                           'ulim':1.01,
		                           'wght':0.0001}})
		self.uniform.update({'incl':{'init':15.,
		                             'llim':0.,
		                             'ulim':60.,
		                             'wght':0.5}})
		self.uniform.update({'anom':{'init':180.,
		                             'llim':0.,
		                             'ulim':360.,
		                             'wght':0.5}})
