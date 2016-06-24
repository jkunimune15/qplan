#
# constraints.py -- custom extensions to the astroplan.Constraint class
#
#  Justin Kunimune (justinku@naoj.org)
#
import numpy as np

# 3rd party imports
import astropy.units as u
from astroplan import Constraint
from astroplan.moon import moon_illumination


class TransparencyConstraint(Constraint):
    """
    Constraint the transparency of the sky
    """
    def __init__(self, min=None, max=None, sdlr=None):
        self.sdlr = sdlr	# this scheudler is used to calculate transparency
        self.min = min
        self.max = max

    def set_sdlr(self, sdlr):
        self.sdlr = sdlr
    
    def compute_constraint(self, times, observer, targets):
        transparency = np.array([self.sdlr.calc_transparency(times)]*len(targets))
        if self.min is None and self.max is not None:
            mask = self.max >= transparency
        elif self.max is None and self.min is not None:
            mask = self.min <= transparency
        elif self.min is not None and self.max is not None:
            mask = ((self.min <= transparency) &
                    (illumination <= self.max))
        else:
            raise ValueError("No max and/or min specified in "
                             "TransparencyConstraint.")
        return mask


class MoonIlluminationConstraint(Constraint):
    """
    Constrain the fractional illumination of the Earth's moon.
    Almsot identical to astrolpan.MoonIlluminationConstraint, but with one bug fixed for.
    """
    def __init__(self, min=None, max=None):
        """
        Parameters
        ----------
        min : float or `None` (optional)
            Minimum acceptable fractional illumination (inclusive). `None`
            indicates no limit.
        max : float or `None` (optional)
            Maximum acceptable fractional illumination (inclusive). `None`
            indicates no limit.
        """
        self.min = min
        self.max = max

    def compute_constraint(self, times, observer, targets):
        illumination = np.array([moon_illumination(times, observer.location)]*len(targets))
	#                       ^                                           ^
	# I added these brackets and the *len(targets) to force illumination to have the correct shape
        if self.min is None and self.max is not None:
            mask = self.max >= illumination
        elif self.max is None and self.min is not None:
            mask = self.min <= illumination
        elif self.min is not None and self.max is not None:
            mask = ((self.min <= illumination) &
                    (illumination <= self.max))
        else:
            raise ValueError("No max and/or min specified in "
                             "MoonSeparationConstraint.")
        return mask

#END
