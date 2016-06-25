#! /usr/bin/env python
#
# qsim.py -- Various functions used for qplan
#
#  Justin Kunimune (justinku@naoj.org)
#
import numpy as np

# 3rd party imports
from ginga.misc import Bunch

import astropy.units as u
import astropy.coordinates
from astropy.coordinates.angle_utilities import angular_separation

import astroplan
from astroplan.constraints import _get_altaz

# the speed at which the telescope slews
slew_rate = 0.5*u.degree/u.second

# the time it takes to change a filter
filter_change_time = 1*u.hour

# define weights (see score() method)
weights = Bunch.Bunch(w_rank=0.3, w_delay=0.2,
                      w_slew=0.2, w_priority=0.1,
                      w_filterchange = 0.3)


def transition(ob1, ob2, start_time, site):
    """
    returns a new TransitionBlock that will fit between the two OBs
    and acount for slewing and filter changes
    """
    components = {}

    aaz = _get_altaz(astropy.time.Time([start_time]), site,
                     [ob1.target, ob2.target])['altaz']
    sep = aaz[0].separation(aaz[1])[0]
    components['slew_time'] = sep / slew_rate

    if ob1.configuration['filter'] != ob2.configuration['filter']:
        components['filter_change'] = filter_change_time
    
    output = astroplan.TransitionBlock(components=components, start_time=start_time)
    output.components = components
    return output


def score(tb, ob):
    """
    calculates a score for an OB thar represents how well it would fit here
    Remember: Lower numbers are better
    """
    score = 0
    score += ob.priority*weights.w_priority	# the marked priority
    score += ob.program.rank*weights.w_rank	# the rank of the program
    if tb is not None:
	# the time to slew
        score += tb.components.get('slew_time', 0*u.second).value*weights.w_slew
	# the time to change filter
        score += tb.components.get('filter_change',
                                   0*u.second).value*weights.w_filterchange
        
    return score

    
def const_score(tb, ob, times, site):
    """
    calculates a score for an OB thar represents how well it would fit here,
    infinity if constraints make it unobservable
    Remember: Lower numbers are better
    """
    moon = astroplan.get_moon(times, site.location, site.pressure)

    altaz = site.altaz(times, ob.target)
    if np.any(altaz.alt < ob.constraints[1].min) or np.any(altaz.alt > ob.constraints[1].max):
        return float('inf')	# altitude constraint
    if np.any(site.moon_illumination(times) > ob.constraints[2].max):
        return float('inf')	# moon illumination constraint
    moon_sep = astropy.coordinates.Angle([angular_separation(moon.spherical.lon,
                                                             moon.spherical.lat,
                                                             ob.target.coord.spherical.lon,
                                                             ob.target.coord.spherical.lat)])
    if np.any(moon_sep < ob.constraints[0].min):
        return float('inf')	# moon separation constraint
    score = 0.
    score += ob.priority*weights.w_priority	# the marked priority
    score += ob.program.rank*weights.w_rank	# the rank of the program
    if tb is not None:
	# the time to slew
        score += tb.components.get('slew_time', 0*u.second).value*weights.w_slew	
	# the time to change filter
        score += tb.components.get('filter_change',
                                   0*u.second).value*weights.w_filterchange
        
    return score


def eval_schedule(schedule):
    """
    counts the number of filter changes and the wasted seconds
    and returns them in a Bunch
    """

    current_filter = None
    num_filter_exchanges = 0
    time_waste_sec = 0.0

    for ob in schedule:
        if type(ob).__name__ == "TransitionBlock":
            delta = (ob.end_time-ob.start_time).to(u.second).value
            time_waste_sec += delta

        elif ((ob.settings['filter'] != None) and
            (ob.settings['filter'] != current_filter)):
            num_filter_exchanges += 1
            current_filter = ob.settings['filter']

    res = Bunch.Bunch(num_filter_exchanges=num_filter_exchanges,
                      time_waste_sec=time_waste_sec)
    return res

#END
