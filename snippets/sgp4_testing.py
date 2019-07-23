"""
Script for testing the effects of TLE params on propagation
"""

from skyfield.sgp4lib import EarthSatellite
from skyfield.api import (
    load, 
    Topos,
    utc,
    )
from datetime import datetime
from astropy.coordinates import (
    Latitude,
    Longitude,
    )
from astropy import units as u

TS = load.timescale() # save repeated use in iterative loops

# location of RASA, La Palma
SITE_LATITUDE = 28.7603135
SITE_LONGITUDE = -17.8796168
SITE_ELEVATION = 2387
TOPOS_LOCATION = Topos(SITE_LATITUDE, 
                       SITE_LONGITUDE, 
                       elevation_m=SITE_ELEVATION)

class TLE:
    """
    Two Line Element
    """
    def __init__(self, line1, line2, name=None):
        """
        Initiates the TLE
        
        Parameters
        ----------
        line1, line2 : str
            First and second lines of the TLE
        name : str, optional
            Name of the object to which the TLE is attributed
        """
        self.line1 = line1
        self.line2 = line2
        if name is not None:
            self.name = name[2:]
        else:
            self.name = 'UNKNOWN'
        
        # ephemeris info
        self.obs = TOPOS_LOCATION
        self.obj = EarthSatellite(line1, line2, name)
        self.ts = TS
    
    def radec(self, epoch):
        """
        Determine radec coords for a given epoch
        
        Parameters
        ----------
        epoch : datetime object
            datetime object [utc] with replaced tzinfo as utc
        
        Returns
        -------
        ra, dec : Longitude and Latitude objects
            Right ascension and declination of object at given epoch
        """
        ra, dec, _ = (self.obj-self.obs).at(self.ts.utc(epoch)).radec()
        
        return Longitude(ra.hours, u.hourangle), Latitude(dec.degrees, u.deg)

if __name__ == "__main__":
	
	## GEO object - SBS 3
	# Investigating effect of t0
	LINE1 = '1 13651U 82110B   19202.77479931 -.00000287 +00000-0 +00000-0 0  9997'
	LINE2 = '2 13651 014.9046 354.9232 0009216 043.0022 143.8103 00.99935313021330'
	tle = TLE(LINE1, LINE2, name='SBS 3')
	
	epoch = datetime.strptime('2019-07-24T12:00:00',
	                          '%Y-%m-%dT%H:%M:%S').replace(tzinfo=utc)
	
	ra, dec = tle.radec(epoch)
	print(ra)
	print(dec)
	
	LINE1 = '1 13651U 82110B   19203.27511181 -.00000287 +00000-0 +00000-0 0  9997'
	LINE2 = '2 13651 014.9046 354.9232 0009216 043.0022 323.8103 00.99935313021330'
	tle = TLE(LINE1, LINE2, name='SBS 3')
	
	ra, dec = tle.radec(epoch)
	print(ra)
	print(dec)
	
	## LEO object - ISS
	LINE1 = '1 25544U 98067A   19203.55917916  .00000561  00000-0  17341-4 0  9990'
	LINE2 = '2 25544  51.6424 185.7766 0006739 167.1766 299.1181 15.50994811180747'
	tle = TLE(LINE1, LINE2, name='ISS')
	
	epoch = datetime.strptime('2019-07-24T12:00:00',
	                          '%Y-%m-%dT%H:%M:%S').replace(tzinfo=utc)
	
	ra, dec = tle.radec(epoch)
	print(ra)
	print(dec)
	
	LINE1 = '1 25544U 98067A   19203.59140138  .00000561  00000-0  17341-4 0  9990'
	LINE2 = '2 25544  51.6424 185.7766 0006739 167.1766 119.1181 15.50994811180747'
	tle = TLE(LINE1, LINE2, name='ISS')
	
	ra, dec = tle.radec(epoch)
	print(ra)
	print(dec)
