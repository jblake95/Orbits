"""
Module for dealing with two-line element sets
"""

import getpass as gp
from spacetrack import SpaceTrackClient
from skyfield.sgp4lib import EarthSatellite
from skyfield.api import (
    load, 
    Topos,
    )
from datetime import datetime
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import (
    SkyCoord, 
    Longitude, 
    Latitude,
    EarthLocation,
    AltAz,
    )

TS = load.timescale() # save repeated use in iterative loops
LE_FORMAT = '3le'     # catalog storage is set up to use 3le format

# location of RASA, La Palma
SITE_LATITUDE = 28.7603135
SITE_LONGITUDE = -17.8796168
SITE_ELEVATION = 2387
SITE_LOCATION = EarthLocation(lat=SITE_LATITUDE*u.deg,
                              lon=SITE_LONGITUDE*u.deg,
                              height=SITE_ELEVATION*u.m)

# topocentric location of RASA, La Palma
TOPOS_LOCATION = Topos(SITE_LATITUDE, 
                       SITE_LONGITUDE, 
                       elevation_m=SITE_ELEVATION)

# checks for orbital type requests
GEO_CHECK = ['g', 'geo']
LEO_CHECK = ['l', 'leo']
MEO_CHECK = ['m', 'meo']
HEO_CHECK = ['h', 'heo']
ALL_CHECK = ['a', 'all']

class ST:
    """
    Space-Track Interface
    """
    def __init__(self):
        """
        Connects user to Space-Track account
        """
        un, pw = self.requestAccess()
        self.username = un
        self.password = pw
        self.client = SpaceTrackClient(identity=un, password=pw)
    
    def requestAccess(self):
        """
        Obtain user access details - requests password from user
        """
        st_un = 'J.Blake@warwick.ac.uk'
        st_pw = gp.getpass('Space-Track password: ')
        
        return st_un, st_pw
    
    def getLatestTLE(self, norad_id):
        """
        Obtain latest TLE for a NORAD object
        
        Parameters
        ----------
        norad_id : int
            NORAD catalogue ID for the object of interest
        
        Returns
        -------
        tle : TLE object
            Latest TLE object for the object of interest
        """
        elset = self.client.tle_latest(norad_cat_id=norad_id,
                                       iter_lines=True,
                                       ordinal=1,
                                       format=LE_FORMAT)
        tle = [line for line in elset]
        
        return TLE(tle[1], tle[2], name=tle[0])
    
    def getSatCat(self, norad_id):
        """
        Obtain the satcat entry for the object of interest
        
        Parameters
        ----------
        norad_id : int | list
            NORAD catalogue ID(s) for the object(s) of interest
        
        Returns
        -------
        satcat : list
            List of directories of information for the object(s) 
            of interest
        """
        return self.client.satcat(norad_cat_id=norad_id)
    
    def getLatestCatalog(self, orb):
        """
        Obtain latest TLE catalog for an orbital type
        
        Parameters
        ----------
        orb : Orbit object
            Orbit object encoding the desired orbital type
        
        Returns
        -------
        catalog : dict
            Directory of latest TLE objects, look up by NORAD ID
        """
        print('Pulling TLEs from Space-Track...')
        if orb.type in GEO_CHECK + LEO_CHECK:
            result = self.client.tle(iter_lines=True,
                                     eccentricity=orb.e_lim,
                                     mean_motion=orb.mm_lim,
                                     epoch='>now-30',
                                     limit=200000,
                                     format=LE_FORMAT)
        elif orb.type in MEO_CHECK:
            result = self.client.tle(iter_lines=True,
                                     eccentricity=orb.e_lim,
                                     period=orb.p_lim,
                                     epoch='>now-30',
                                     limit=200000,
                                     format=LE_FORMAT)
        elif orb.type in HEO_CHECK:
            result = self.client.tle(iter_lines=True,
                                     eccentricity=orb.e_lim,
                                     epoch='>now-30',
                                     limit=200000,
                                     format=LE_FORMAT)
        elif orb.type in ALL_CHECK:
            result = self.client.tle(iter_lines=True,
                                     epoch='>now-30',
                                     limit=200000,
                                     format=LE_FORMAT)
        else:
            print('Incorrect format! Please supply a valid' 
                  'orbit type... \n'
                  'GEO - "g" \n'
                  'LEO - "l" \n'
                  'MEO - "m" \n'
                  'HEO - "h" \n'
                  'ALL - "a" \n')
            quit()
        result = [line for line in result]
        
        # organize into user-friendly format
        i = 0
        catalog = {}
        while i < len(result):
            tle = TLE(result[i+1], 
                      result[i+2], 
                      name=result[i])
            catalog.update({tle.norad_id:tle})
            i += 3
        
        return catalog

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
        
        # book-keeping
        self.norad_id = int(self.line1[2:7])
        self.yday = float(self.line1[20:32])
        
        # orbital properties
        self.inclination = float(self.line2[8:16])
        self.eccentricity = float(self.line2[26:33])
        self.raan = float(self.line2[17:25])
        self.argperigree = float(self.line2[34:42])
        self.mean_anomaly = float(self.line2[43:51])
        self.mean_motion = float(self.line2[52:63])
    
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
    
    def hourangle(self, epoch):
        """
        Determine hourangle for a given epoch
        
        Parameters
        ----------
        epoch : datetime object
            datetime object [utc] with replaced tzinfo as utc
        
        Returns
        -------
        hourangle : Longitude object
            Hourangle of object at given epoch
        """
        ra, _ = self.radec(epoch)
        lst = Time(datetime.utcnow(), 
                   scale='utc', 
                   location=SITE_LOCATION).sidereal_time('apparent')
        
        return (lst - ra).wrap_at(12*u.hourangle)
    
    def altaz(self, epoch):
        """
        Determine altaz coords for a given epoch
        
        Parameters
        ----------
        epoch : datetime object
            datetime object [utc] with replaced tzinfo as utc
        
        Returns
        -------
        alt, az : Angle objects
            Altitude and azimuth of object at given epoch
        """
        ra, dec = self.radec(epoch)
        
        aa = AltAz(location=SITE_LOCATION, obstime=epoch)
        altaz = SkyCoord(ra, dec).transform_to(aa)
        
        return altaz.alt, altaz.az

class SatCat:
    """
    Satellite catalog entry
    """
    def __init__(self, entry):
        """
        Initiate satcat entry
        
        Parameters
        ----------
        entry : dict
            Directory of satcat information for the object of interest
        """
        # launch info
        self.launchyear = entry['LAUNCH_YEAR']
        self.launchdate = entry['LAUNCH']
        self.launchsite = entry['SITE']
        
        # object info
        self.norad_id = entry['NORAD_CAT_ID']
        self.name = entry['OBJECT_NAME']
        self.objtype = entry['OBJECT_TYPE']
        self.country = entry['COUNTRY']
        self.size = entry['RCS_SIZE']
        
        # orbit info
        self.period = entry['PERIOD']
        self.apogee = entry['APOGEE']
        self.perigee = entry['PERIGEE']
        self.inclination = entry['INCLINATION']

class Orbit:
    """
    Convenience class for orbit-specific searches
    """
    def __init__(self, orb_type):
        """
        Initiate Orbit object using SpaceTrack definitions
        
        Parameters
        ----------
        orb_type : str
            Desired type of orbit
            'g' - GEO
            'l' - LEO
            'm' - MEO
            'h' - HEO
            'a' - ALL
        """
        self.type = orb_type.lower()
        if self.type in GEO_CHECK:
            self.e_lim = '<0.01'
            self.mm_lim = '0.99--1.01'
        elif self.type in LEO_CHECK:
            self.e_lim = '<0.25'
            self.mm_lim = '>11.25'
        elif self.type in MEO_CHECK:
            self.e_lim = '<0.25'
            self.p_lim = '600--800'
        elif self.type in HEO_CHECK:
            self.e_lim = '>0.25'
        elif self.type in ALL_CHECK:
            print('Full catalogue specified; no limits placed.')
        else:
            print('Incorrect format! Please provide a valid' 
                  'orbit type... \n'
                  'GEO - "g" \n'
                  'LEO - "l" \n'
                  'MEO - "m" \n'
                  'HEO - "h" \n'
                  'ALL - "a" \n')
            quit()
