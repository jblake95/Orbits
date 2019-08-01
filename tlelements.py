"""
Convenience classes for manipulation of TLE parameters

A note on default values:
-------------------------
Each parameter class has an in-built default
Norad_id     - 0 (next field set to Unclassified 'U' by default)
Launch       - empty (Intl designator requires all three components)
Epoch        - 01/01/2000 (WARNING: TLEs quickly deteriorate)
Mmdot        - 0 (not used by SGP4/SDP4, only SGP)
Mmdot2       - 0 (not used by SGP4/SDP4, only SGP)
Drag         - 0 (ballistic coefficient adjusted for atm. density)
Housekeeping - 0 (element set number and revolution number)
Checksum     - 0 (modulo-10 checksum summing all numbers with '-'=1)
Inclination  - 0
RAAN         - 0
Eccentricity - 0
ArgPerigee   - 0
MeanAnomaly  - 0
Mm           - 0
"""

from utils import (
    fractional_yday,
    n_digits,
    tle_standard_form,
    )
from datetime import datetime
from astropy.coordinates import Angle
from astropy import units as u

class Norad_id:
    """
    NORAD identification number
    """
    def __init__(self, norad_id=0):
        """
        Initiate Norad_id object
        
        Parameters
        ----------
        norad_id : int, optional
            NORAD identification number
            default = 0
        """
        self._validate(norad_id)
        self.entry = self._formatTLE()
    
    def __str__(self):
        """
        Printable string for Norad_id object
        """
        return 'Norad_id <{}>'.format(self.entry)
    
    def _validate(self, norad_id):
        """
        Validate the input
        """
        assert isinstance(norad_id, int), "NORAD ID must be int"
        assert norad_id >= 0, "Norad ID cannot be negative"
        assert n_digits(norad_id) < 6, "NORAD ID exceeds limit"
        
        self.value = norad_id
    
    def _formatTLE(self):
        """
        TLE-friendly format
        """
        return '{:5.0f}'.format(self.value)

class Launch:
    """
    Launch information for international designator
    """
    def __init__(self, launch_yr=None, launch_no=None, launch_id=None):
        """
        Initiate Launch object - international designator can be left
        blank but all three inputs are required if either is given
        
        Parameters
        ----------
        launch_yr : datetime object
            Year of launch in datetime format
        launch_no : int
            Launch number for that year
        launch_id : str
            Identifier for piece of launch (e.g. 'ABC')
        """
        if (launch_yr is not None and
            launch_no is not None and
            launch_id is not None):
                    self._validate(launch_yr, launch_no, launch_id)
                    [self.yr_entry, 
                     self.no_entry, 
                     self.id_entry] = self._formatTLE()
        else:
            [self.yr_entry, 
             self.no_entry, 
             self.id_entry] = self._formatTLE(empty=True)
    
    def __str__(self):
        """
        Printable string for Launch object
        """
        return 'International designator <{}{}{}>'.format(self.yr_entry,
                                                          self.no_entry,
                                                          self.id_entry)
    
    def _validate(self, launch_yr, launch_no, launch_id):
        """
        Validate the input
        """
        assert isinstance(launch_yr, datetime), "Launch date must be datetime"
        assert isinstance(launch_no, int), "Launch no must be int"
        assert n_digits(launch_no) < 4, "Launch no exceeds limit"
        assert isinstance(launch_id, str), "Launch ID must be str"
        assert launch_id.isalpha(), "Invalid launch ID; ABC not 123"
        
        self.yr_value = launch_yr.year
        self.no_value = launch_no
        self.id_value = launch_id
        return None
    
    def _formatTLE(self, empty=False):
        """
        TLE-friendly format
        
        Parameters
        ----------
        empty : bool, optional
            Toggle to select default in case of insufficient info
            default = False
        """
        if not empty:
            return ['{:02.0f}'.format(self.yr_value % 100),
                    '{:03.0f}'.format(self.no_value),
                    '{:3}'.format(self.id_value)]
        else:
            return['{:2}'.format(''),
                   '{:3}'.format(''),
                   '{:3}'.format('')]

class Epoch:
    """
    Reference epoch
    """
    def __init__(self, epoch=None):
        """
        Initiate Epoch object
        
        Parameters
        ----------
        epoch : datetime object
            Reference epoch in datetime format
        """
        if epoch is not None:
            self._validate(epoch)
            self.yr_entry, self.yday_entry = self._formatTLE()
        else:
            self.yr_entry, self.yday_entry = self._formatTLE(empty=True)
    
    def __str__(self):
        """
        Printable string for Epoch object
        """
        return 'Reference epoch <{}{}>'.format(self.yr_entry,
                                               self.yday_entry)
    
    def _validate(self, epoch):
        """
        Validate input
        """
        assert isinstance(epoch, datetime), "Epoch must be datetime"
        
        self.yr_value = epoch.year
        self.yday_value = fractional_yday(epoch)
    
    def _formatTLE(self, empty=False):
        """
        TLE-friendly format
        
        Parameters
        ----------
        empty : bool, optional
            Toggle to select default in case of insufficient info
            default = False (01/01/2000)
        """
        if not empty:
            return ['{:02.0f}'.format(self.yr_value % 100),
                    '{:012.8f}'.format(self.yday_value)]
        else:
            print('WARNING: TLE propagation requires reference epoch!')
            return['{:02}'.format(0),
                   '{:012.8f}'.format(1)]

class Mmdot:
    """
    First time derivative of mean motion
    """
    def __init__(self, mmdot=0):
        """
        Initiate Mmdot object
        
        Parameters
        ----------
        mmdot : float, optional
            First time derivative of mean motion
            default = 0
        """
        self._validate(mmdot)
        self.entry = self._formatTLE()
    
    def __str__(self):
        """
        Printable string for Mmdot object
        """
        return 'Mmdot <{}>'.format(self.entry)
    
    def _validate(self, mmdot):
        """
        Validate input
        """
        assert isinstance(mmdot, (float,int)), "Mmdot must be float|int"
        assert abs(mmdot) < 1, "Invalid mmdot, exceeds limit"
        
        self.value = mmdot
    
    def _formatTLE(self):
        """
        TLE-friendly format
        """
        if self.value >= 0 :
            return '+' + '{:9.8f}'.format(self.value).lstrip('0')
        else:
            return '-' + '{:9.8f}'.format(self.value).lstrip('-0')

class Mmdot2:
    """
    Second time derivative of mean motion
    """
    def __init__(self, mmdot2=0):
        """
        Initiate Mmdot2 object
        
        Parameters
        ----------
        mmdot2 : float, optional
            Second time derivative of mean motion
            default=0
        """
        self._validate(mmdot2)
        self.entry = self._formatTLE()
    
    def __str__(self):
        """
        Printable string for Mmdot2 object
        """
        return 'Mmdot2 <{}>'.format(self.entry)
    
    def _validate(self, mmdot2):
        """
        Validate input
        """
        assert isinstance(mmdot2, (float,int)), "Mmdot2 must be float|int"
        
        self.value = mmdot2
    
    def _formatTLE(self):
        """
        TLE-friendly format
        """
        return tle_standard_form(self.value)

class Drag:
    """
    B-star drag term
    """
    def __init__(self, drag=0):
        """
        Initiate Drag object
        
        Parameters
        ----------
        drag : float, optional
            B-star drag term
            default = 0
        """
        self._validate(drag)
        self.entry = self._formatTLE()
    
    def __str__(self):
        """
        Printable string for Drag object
        """
        return 'Drag <{}>'.format(self.entry)
    
    def _validate(self, drag):
        """
        Validate input
        """
        assert isinstance(drag, (float,int)), "Drag must be float|int"
        
        self.value = drag
    
    def _formatTLE(self):
        """
        TLE-friendly format
        """
        return tle_standard_form(self.value)

class Housekeeping:
    """
    Elements for general housekeeping
    """
    def __init__(self, set_no=0, rev_no=0):
        """
        Initiate Housekeeping object
        
        Parameters
        ----------
        set_no : int, optional
            Element set number
            default = 0
        rev_no : int, optional
            Revolution number
            default = 0
        """
        self._validate(set_no, rev_no)
        self.set_entry, self.rev_entry = self._formatTLE()
    
    def __str__(self):
        """
        Printable string for Housekeeping object
        """
        return 'Set No <{}> Rev No <{}>'.format(self.set_entry,
                                                self.rev_entry)
    
    def _validate(self, set_no, rev_no):
        """
        Validate input
        """
        assert isinstance(set_no, int), "Set no must be int"
        assert isinstance(rev_no, int), "Rev no must be int"
        assert n_digits(set_no) < 5, "Set no exceeds limit"
        assert n_digits(rev_no) < 6, "Rev no exceeds limit"
        
        self.set_value = set_no
        self.rev_value = rev_no
    
    def _formatTLE(self):
        """
        TLE-friendly format
        """
        return ['{:4}'.format(self.set_value),
                '{:5}'.format(self.rev_value)]

class Checksum:
    """
    Checksum
    """
    def __init__(self, line1=None, line2=None):
        """
        Initiate Checksum object
        
        Parameters
        ----------
        line1 : list, optional
            Line 1 of TLE in list format
            default = None
        line2 : list, optional
            Line 2 of TLE in list format,
            default = None
        """
        if line1 is not None and line2 is not None:
            self.line1_value = self._calculate(line1)
            self.line2_value = self._calculate(line2)
        else:
            self.line1_value = 0
            self.line2_value = 0
        self.line1_entry, self.line2_entry = self._formatTLE()
    
    def __str__(self):
        """
        Printable string for Checksum object
        """
        return 'Checksums <{},{}>'.format(self.line1_value,
                                          self.line2_value)
    
    def _calculate(self, line):
        """
        Calculate checksum
        """
        checksum = 0
        for character in line:
            if character.isnumeric():
                checksum += int(character)
            elif character.isalpha():
                continue
            elif character.isspace():
                continue
            else:
                if character == '-':
                    checksum += 1
                else:
                    continue
        return checksum % 10
    
    def _formatTLE(self):
        """
        TLE-friendly format
        """
        return ['{:1}'.format(self.line1_value),
                '{:1}'.format(self.line2_value)]

class Inclination:
    """
    Orbital inclination
    """
    def __init__(self, incl=0):
        """
        Initiate Inclination object
        
        Parameters
        ----------
        incl : float, optional
            Orbital inclination [deg]
            default = 0
        """
        self._validate(incl)
        self.entry = self._formatTLE()
    
    def __str__(self):
        """
        Printable string for Inclination object
        """
        return 'Inclination <{}>'.format(self.entry)
    
    def _validate(self, incl):
        """
        Validate input
        """
        assert isinstance(incl, (float,int)), "Inclination must be float|int"
        
        self.value = Angle(incl, u.deg)
    
    def _formatTLE(self):
        """
        TLE-friendly format
        """
        return '{:8.4f}'.format(self.value.deg)

class RAAN:
    """
    Right ascension of the ascending node
    """
    def __init__(self, raan=0):
        """
        Initiate RAAN object
        
        Parameters
        ----------
        raan : float, optional
            Right ascension of the ascending node [deg]
            default = 0
        """
        self._validate(raan)
        self.entry = self._formatTLE()
    
    def __str__(self):
        """
        Printable string for RAAN object
        """
        return 'RAAN <{}>'.format(self.entry)
    
    def _validate(self, raan):
        """
        Validate input
        """
        assert isinstance(raan, (float,int)), "Raan must be float|int"
        
        self.value = Angle(raan, u.deg)
    
    def _formatTLE(self):
        """
        TLE-friendly format
        """
        return '{:8.4f}'.format(self.value.deg)

class Eccentricity:
    """
    Orbital eccentricity
    """
    def __init__(self, ecc = 0):
        """
        Initiate Eccentricity object
        
        Parameters
        ----------
        ecc : float, optional
            Orbital eccentricity
            default = 0
        """
        self._validate(ecc)
        self.entry = self._formatTLE()
    
    def __str__(self):
        """
        Printable string for Eccentricity object
        """
        return 'Eccentricity <{}>'.format(self.entry)
    
    def _validate(self, ecc):
        """
        Validate input
        """
        assert isinstance(ecc, (float,int)), "Eccentricity must be float|int"
        
        self.value = ecc
    
    def _formatTLE(self):
        """
        TLE-friendly format
        """
        return '{:.7f}'.format(self.value).split('.')[1]

class ArgPerigee:
    """
    Argument of perigee
    """
    def __init__(self, argp=0):
        """
        Initiate ArgPerigee object
        
        Parameters
        ----------
        argp : float, optional
            Argument of perigee [deg]
            default = 0
        """
        self._validate(argp)
        self.entry = self._formatTLE()
    
    def __str__(self):
        """
        Printable string for ArgPerigee object
        """
        return 'ArgPerigee <{}>'.format(self.entry)
    
    def _validate(self, argp):
        """
        Validate input
        """
        assert isinstance(argp, (float,int)), "Argp must be float|int"
        
        self.value = argp
    
    def _formatTLE(self):
        """
        TLE-friendly format
        """
        return '{:8.4f}'.format(self.value)

class MeanAnomaly:
    """
    Mean Anomaly
    """
    def __init__(self, anom=0):
        """
        Initiate MeanAnomaly object
        
        Parameters
        ----------
        anom : float, optional
            Mean anomaly [deg]
            default = 0
        """
        self._validate(anom)
        self.entry = self._formatTLE()
    
    def __str__(self):
        """
        Printable string for MeanAnomaly object
        """
        return 'MeanAnomaly <{}>'.format(self.entry)
    
    def _validate(self, anom):
        """
        Validate input
        """
        assert isinstance(anom, (float,int)), "Anom must be float|int"
        
        self.value = anom
    
    def _formatTLE(self):
        """
        TLE-friendly format
        """
        return '{:8.4f}'.format(self.value)

class Mm:
    """
    Mean motion
    """
    def __init__(self, mm=0):
        """
        Initiate Mm object
        
        Parameters
        ----------
        mm : float, optional 
            Mean motion [revs per day]
            default = 0
        """
        self._validate(mm)
        self.entry = self._formatTLE()
    
    def __str__(self):
        """
        Printable string for Mm object
        """
        return 'Mm <{}>'.format(self.entry)
    
    def _validate(self, mm):
        """
        Validate input 
        """
        assert isinstance(mm, (float,int)), "Mm must be float|int"
        
        self.value = mm
    
    def _formatTLE(self):
        """
        TLE-friendly format
        """
        return '{:11.8f}'.format(self.value)
