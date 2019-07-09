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
		self.ecc_val = 0.
		
		# uniform prior params
		self.argp_init = 180.
		self.argp_llim = 0.
		self.argp_ulim = 360.
		self.argp_wght = 0.5
		
		self.raan_init = 180.
		self.raan_llim = 0.
		self.raan_ulim = 360.
		self.raan_wght = 0.5
		
		self.mmdot_init = 0.
		self.mmdot_llim = -0.000004
		self.mmdot_ulim = 0.000002
		self.mmdot_wght = 0.00000001
		
		self.mm_init = 1.
		self.mm_llim = 0.99
		self.mm_ulim = 1.01
		self.mm_wght = 0.0001
		
		self.incl_init = 15.
		self.incl_llim = 0.
		self.incl_ulim = 60.
		self.incl_wght = 0.5
		
		self.anom_init = 180.
		self.anom_llim = 0.
		self.anom_ulim = 360.
		self.anom_wght = 0.5
