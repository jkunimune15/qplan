#
# entity.py -- various entities used by queue system
#
#  Eric Jeschke (eric@naoj.org)
#
from datetime import tzinfo, timedelta, datetime
import string
import math
import pytz

# 3rd party imports
import ephem
import numpy

import astroplan as apn
import astropy.units as u
# ZOPE imports
try:
    import persistent

    class PersistentEntity(persistent.Persistent):

        def __init__(self):
            super(PersistentEntity, self).__init__()

except ImportError:
    # No ZOPE, so define as an ordinary base object
    class PersistentEntity(object):

        def __init__(self):
            super(PersistentEntity, self).__init__()
            self._p_changed = False

from ginga.misc import Bunch

# Subaru defines a dark night as one that is 2-3 days before or
# after a new moon (0%).  Since a half moon (50%) occurs just 7 days
# prior to a new moon, we can roughly calculate a dark night as
# being < 25% illumination
dark_moon_lim = 0.25


class Program(PersistentEntity):
    """
    Program
    Defines a program that has been accepted for observation.
    """
    def __init__(self, proposal, pi='', observers='', rank=1.0,
                 propid=None, grade=None, partner=None, hours=None,
                 category=None, instruments=[], description=None,
                 skip=False):
        super(Program, self).__init__()

        self.proposal = proposal
        if propid == None:
            # TODO: is there an algorithm to turn proposals
            # into propids?
            propid = proposal
        self.propid = propid
        self.pi = pi
        self.observers = observers
        self.rank = rank
        self.grade = grade
        self.partner = partner
        self.category = category.lower()
        self.instruments = map(string.upper, instruments)
        self.total_time = hours * 3600.0
        # TODO: eventually this will contain all the relevant info
        # pertaining to a proposal
        self.skip = skip

    def __repr__(self):
        return self.proposal

    __str__ = __repr__


class BaseTarget(object):
    pass

class StaticTarget(BaseTarget):
    def __init__(self, name=None, ra=None, dec=None, equinox=2000.0,
                 comment=''):
        super(StaticTarget, self).__init__()
        self.name = name
        self.ra = ra
        self.dec = dec
        self.equinox = equinox
        self.comment = comment

        if self.ra is not None:
            self._recalc_body()

    def _recalc_body(self):
        self.xeph_line = "%s,f|A,%s,%s,0.0,%s" % (
            self.name[:20], self.ra, self.dec, self.equinox)
        self.body = ephem.readdb(self.xeph_line)


    def import_record(self, rec):
        code = rec.code.strip()
        self.name = rec.name
        self.ra = rec.ra
        self.dec = rec.dec

        # transform equinox, e.g. "J2000" -> 2000
        eq = rec.eq
        if isinstance(eq, str):
            eq = eq.upper()
            if eq[0] in ('B', 'J'):
                eq = eq[1:]
                eq = float(eq)
        eq = int(eq)
        self.equinox = eq
        self.comment = rec.comment.strip()

        self._recalc_body()
        return code

    def calc(self, observer, time_start):
        return CalculationResult(self, observer, time_start)

    # for pickling

    def __getstate__(self):
        d = self.__dict__.copy()
        # ephem objects can't be pickled
        d['body'] = None
        return d

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.body = ephem.readdb(self.xeph_line)


class TelescopeConfiguration(object):
    """Contains information about the telescope that must be true during an OB"""

    def __init__(self, focus=None, dome=None, comment=''):
        super(TelescopeConfiguration, self).__init__()
        self.focus = focus
        if dome is None:
            dome = 'open'
        else:
            dome = dome.lower()
        self.dome = dome
        self.comment = comment

    def import_record(self, rec):
        code = rec.code.strip()
        self.focus = rec.focus.upper()
        self.dome = rec.dome.lower()
        self.comment = rec.comment.strip()
        return code
    
    def get_constraints(self):	# returns a list of Constraints representing this cfg
        return []
    
    def get_configuration(self):	# return a dictionary representing this cfg
        return {'focus':self.focus, 'dome':self.dome}


class InstrumentConfiguration(object):

    def __init__(self):
        super(InstrumentConfiguration, self).__init__()

        self.insname = None
        self.mode = None
        self.comment = ''

class SPCAMConfiguration(InstrumentConfiguration):	#TODO: find out if I need all these configurations

    def __init__(self, filter=None, guiding=False, num_exp=1, exp_time=10,
                 mode='IMAGE', offset_ra=0, offset_dec=0, pa=90,
                 dith1=60, dith2=None):
        super(SPCAMConfiguration, self).__init__()

        self.insname = 'SPCAM'
        if filter is not None:
            filter = filter.lower()
        self.filter = filter
        self.dither = dither
        self.guiding = guiding
        self.num_exp = num_exp
        self.exp_time = exp_time
        self.mode = mode
        self.offset_ra = offset_ra
        self.offset_dec = offset_dec
        self.pa = pa
        self.dith1 = dith1
        if dith2 == None:
            # TODO: defaults for this depends on mode
            dith2 = 0
        self.dith2 = dith2

    def calc_filter_change_time(self):
        # TODO: this needs to become more accurate
        filter_change_time_sec = 10.0 * 60.0
        return filter_change_time_sec

class HSCConfiguration(InstrumentConfiguration):
    """Contains information about the ??? that must be true during an OB"""
			#TODO: Figure out what HSC stands for
    def __init__(self, filter=None, guiding=False, num_exp=1, exp_time=10,
                 mode='IMAGE', dither=1, offset_ra=0, offset_dec=0, pa=90,
                 dith1=60, dith2=None, skip=0, stop=None, comment=''):
        super(HSCConfiguration, self).__init__()

        self.insname = 'HSC'
        self.mode = mode
        if filter is not None:
            filter = filter.lower()
        self.filter = filter
        self.dither = dither
        self.guiding = guiding
        self.num_exp = int(num_exp)
        self.exp_time = float(exp_time)
        self.offset_ra = offset_ra
        self.offset_dec = offset_dec
        self.pa = pa
        self.dith1 = dith1
        if dith2 is None:
            # TODO: defaults for this depends on mode
            dith2 = 0
        self.dith2 = dith2
        self.skip = skip
        if stop is None:
            stop = num_exp
        self.stop = stop
        self.comment = comment

    def calc_filter_change_time(self):
        # TODO: this needs to become more accurate
        filter_change_time_sec = 35.0 * 60.0
        return filter_change_time_sec

    def import_record(self, rec):
        code = rec.code.strip()
        self.insname = 'HSC'
        self.filter = rec.filter.lower()
        self.mode = rec.mode
        self.dither = rec.dither
        self.guiding = rec.guiding in ('y', 'Y', 'yes', 'YES')
        self.num_exp = int(rec.num_exp)
        self.exp_time = float(rec.exp_time)
        self.pa = float(rec.pa)
        self.offset_ra = float(rec.offset_ra)
        self.offset_dec = float(rec.offset_dec)
        self.dith1 = float(rec.dith1)
        self.dith2 = float(rec.dith2)
        self.skip = int(rec.skip)
        self.stop = int(rec.stop)
        self.comment = rec.comment.strip()
        return code
    
    def get_constraints(self):	# return a list of Constraints representing this cfg
        return []
    
    def get_configuration(self):	# return a dictionary representing this cfg
        return {'mode':self.mode, 'filter':self.filter, 'dither':self.dither,
                'guiding':self.guiding, 'num_exp':self.num_exp, 'exp_time':self.exp_time,
                'offset_ra':self.offset_ra, 'offset_dec':self.offset_dec, 'pa':self.pa,
                'dith1':self.dith1, 'dith2':self.dith2, 'skip':self.skip, 'stop':self.stop}

class FOCASConfiguration(InstrumentConfiguration):

    def __init__(self, filter=None, guiding=False, num_exp=1, exp_time=10,
                 mode='IMAGE', binning='1x1', offset_ra=0, offset_dec=0,
                 pa=0, dither_ra=5, dither_dec=5, dither_theta=0.0):
        super(FOCASConfiguration, self).__init__()

        self.insname = 'FOCAS'
        self.mode = mode
        if filter is not None:
            filter = filter.lower()
        self.filter = filter
        self.guiding = guiding
        self.num_exp = int(num_exp)
        self.exp_time = float(exp_time)
        self.pa = float(pa)
        self.binning = binning
        self.offset_ra = float(offset_ra)
        self.offset_dec = float(offset_dec)
        self.dither_ra = float(dither_ra)
        self.dither_dec = float(dither_dec)
        self.dither_theta = float(dither_theta)

    def calc_filter_change_time(self):
        # TODO: this needs to become more accurate
        filter_change_time_sec = 30.0
        return filter_change_time_sec

    def import_record(self, rec):
        code = rec.code.strip()
        self.insname = 'FOCAS'
        self.mode = rec.mode
        self.filter = rec.filter.lower()
        self.guiding = rec.guiding in ('y', 'Y', 'yes', 'YES')
        self.num_exp = int(rec.num_exp)
        self.exp_time = float(rec.exp_time)
        self.pa = float(rec.pa)
        self.offset_ra = float(rec.offset_ra)
        self.offset_dec = float(rec.offset_dec)
        self.dither_ra = float(rec.dither_ra)
        self.dither_dec = float(rec.dither_dec)
        self.dither_theta = float(rec.dither_theta)
        self.binning = rec.binning
        return code

class EnvironmentConfiguration(object):
    """Contains information about the environment that must be true during an OB"""

    def __init__(self, seeing=None, airmass=None, moon='any',
                 transparency=None, moon_sep=None, comment=''):
        super(EnvironmentConfiguration, self).__init__()
        self.seeing = seeing			# ???
        self.airmass = airmass			# maximum allowable air mass
        self.transparency = transparency	# minimum allowable air transparency
        self.moon_sep = moon_sep		# minimum allowable distance from moon
        if (moon == None) or (len(moon) == 0):
            moon = 'any'
        self.moon = moon.lower()		# desired moon phase
        self.comment = comment

    def import_record(self, rec):
        code = rec.code.strip()

        seeing = rec.seeing.strip()
        if len(seeing) != 0:
            self.seeing = float(seeing)
        else:
            self.seeing = None

        airmass = rec.airmass.strip()
        if len(airmass) != 0:
            self.airmass = float(airmass)
        else:
            self.airmass = None

        self.moon = rec.moon
        self.moon_sep = float(rec.moon_sep)*u.degree
        self.transparency = float(rec.transparency)	#TODO: Account for transparency
        self.comment = rec.comment.strip()
        return code
    
    def get_constraints(self):	# return a list of Constraints representing this cfg
        output = []
        output.append(apn.AirmassConstraint(self.airmass))
        #TODO: TransparencyConstraint
        output.append(apn.MoonSeparationConstraint(self.moon_sep))
        if self.moon == 'dark':
            output.append(apn.MoonIlluminationConstraint(dark_moon_lim))
        return output
    
    def get_configuration(self):	# return a dictionary representing this cfg
        return {}

#END
