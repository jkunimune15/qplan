#! /usr/bin/env python
#
# qsim.py -- Observing Queue Planner
#
#  Eric Jeschke (eric@naoj.org)
#
from datetime import timedelta
import time

# Gen2 imports
from ginga.misc import Bunch

# local imports
import misc
#import constraints
import entity


# maximum rank for a program
max_rank = 10.0

# minimum slot size in sec
minimum_slot_size = 60.0
#minimum_slot_size = 10.0

# telescope parked position
#parked_az_deg = -90.0
parked_az_deg = 270.0
parked_alt_deg = 90.0

# Subaru defines a dark night as one that is 2-3 days before or
# after a new moon (0%).  Since a half moon (50%)occurs just 7 days
# prior to a new moon, we can roughly calculate a dark night as
# being < 25% illumination
dark_night_moon_pct_limit = 0.25


def filterchange_ob(ob, total_time):
    """
    copies ob and adds a filter change comment
    """
    new_ob = entity.OB(program=ob.program, target=ob.target,
                       telcfg=ob.telcfg,
                       inscfg=ob.inscfg, envcfg=ob.envcfg,
                       total_time=total_time, derived=True,
                       comment="Filter change for %s" % (ob))
    return new_ob


def longslew_ob(prev_ob, ob, total_time):
    """
    Copies ob with a Long slew comment and a time of total_time. Also removes the filter of the inscfg if prev_ob is None
    """
    inscfg = ob.inscfg
    new_ob = entity.OB(program=ob.program, target=ob.target,
                       telcfg=ob.telcfg,
                       inscfg=inscfg, envcfg=ob.envcfg,
                       total_time=total_time, derived=True,
                       comment="Long slew for %s" % (ob))
    return new_ob


def calibration_ob(ob, sdss_target, total_time):
    """
    Copies ob with a target of sdss_target, a time of total_time, and a comment about SDSS.
    """
    new_ob = entity.OB(program=ob.program, target=sdss_target,
                       telcfg=ob.telcfg,
                       inscfg=ob.inscfg, envcfg=ob.envcfg,
                       total_time=total_time, derived=True,
                       comment="SDSS calibration for %s" % (ob))
    return new_ob


def delay_ob(ob, total_time):
    """
    Copies ob with a total time of total_time
    """
    new_ob = entity.OB(program=ob.program, target=ob.target,
                       telcfg=ob.telcfg,
                       inscfg=ob.inscfg, envcfg=ob.envcfg,
                       total_time=total_time, derived=True,
                       comment="Delay for %s visibility" % (ob))
    return new_ob


def setup_ob(ob, total_time):
    """
    Copies ob with a time of total_time and a comment about setting up
    """
    d = dict(obid=str(ob), obname=ob.name,
             comment=ob.comment,    # root OB's comment
             proposal=ob.program.proposal)
    # make this derived OB's comment include root OB comment
    comment = "%(proposal)s %(obname)s: %(comment)s" % d
    new_ob = entity.OB(program=ob.program, target=ob.target,
                       telcfg=ob.telcfg,
                       inscfg=ob.inscfg, envcfg=ob.envcfg,
                       total_time=total_time, derived=True,
                       comment="Setup OB: %s" % (comment))
    return new_ob


def obs_to_slots(logger, slots, site, obs, check_moon=False, check_env=False):
    """
    returns a dictionary where each slot has a list of obs that fit into it
    """
    obmap = {}
    for slot in slots:
        key = str(slot)
        obmap[key] = []
        if slot.size() < minimum_slot_size:
            continue
        for ob in obs:
            # this OB OK for this slot at this site?
            res = check_slot(site, None, slot, ob,
                             check_moon=check_moon, check_env=check_env)
            if res.obs_ok:
                obmap[key].append(ob)
            else:
                logger.debug("OB %s no good for slot because: %s" % (
                    ob, res.reason))

    return obmap

def calc_slew_time(cur_alt_deg, cur_az_deg, to_alt_deg, to_az_deg):
    """
    calculates the slew time between two sets of coordinates
    """
    delta_alt, delta_az = to_alt_deg - cur_alt_deg, to_az_deg - cur_az_deg

    slew_sec = misc.calc_slew_time(delta_az, delta_alt)
    return slew_sec


def check_schedule_invariant(site, schedule, ob):
    """
    returns a Bunch answering the question, Can this ob fit into this schedule at this site?
    """
    res = Bunch.Bunch(ob=ob, obs_ok=False, reason="No good reason!")

    # check if instrument will be installed
    if not (ob.inscfg.insname in schedule.data.instruments):
        res.setvals(obs_ok=False, reason="Instrument '%s' not installed" % (
            ob.inscfg.insname))
        return res

    # check if filter will be installed
    if not (ob.inscfg.filter in schedule.data.filters):
        res.setvals(obs_ok=False, reason="Filter '%s' not installed [%s]" % (
            ob.inscfg.filter, schedule.data.filters))
        return res

    # check if this schedule can take this category
    if not ob.program.category in schedule.data.categories:
        res.setvals(obs_ok=False,
                    reason="Slot cannot take category '%s'" % (
            ob.program.category))
        return res

    res.setvals(obs_ok=True)
    return res


def check_night_visibility(site, schedule, ob):
    """
    returns a bunch representing the answer to,
    Will this ob be visible during this schedule at this site, and when?
    """
    res = Bunch.Bunch(ob=ob, obs_ok=False, reason="No good reason!")

    min_el, max_el = ob.telcfg.get_el_minmax()

    # is this target visible during this night, and when?
    (obs_ok, t_start, t_stop) = site.observable(ob.target,
                                                schedule.start_time,
                                                schedule.stop_time,
                                                min_el, max_el, ob.total_time,
                                                airmass=ob.envcfg.airmass,
                                                moon_sep=ob.envcfg.moon_sep)
    if not obs_ok:
        res.setvals(reason="Time or visibility of target")

    res.setvals(obs_ok=obs_ok, start_time=t_start, stop_time=t_stop)
    return res


def check_moon_cond(site, start_time, stop_time, ob, res):
    """Check whether the moon is at acceptable darkness for this OB
    and an acceptable distance from the target.
    """
    c1 = ob.target.calc(site, start_time)
    c2 = ob.target.calc(site, stop_time)

    # is this a dark night? check moon illumination
    is_dark_night = c1.moon_pct <= dark_night_moon_pct_limit

    desired_moon_sep = ob.envcfg.moon_sep

    # if the moon is down for entire exposure, override illumination
    # and consider this a dark night
    horizon_deg = 0.0   # change as necessary
    if (c1.moon_alt < horizon_deg) and (c2.moon_alt < horizon_deg):
        #print("moon down, dark night")
        is_dark_night = True

    # if observer specified a moon phase, check it now
    if ob.envcfg.moon == 'dark':
        ## print "moon pct=%f moon alt=%f moon_sep=%f" % (
        ##     c1.moon_pct, c1.moon_alt, c1.moon_sep)
        if not is_dark_night:
            return False

    # override the observer's desired separation if it is a dark night
    limit_sep = min(30.0, desired_moon_sep)
    if (desired_moon_sep is not None) and is_dark_night:
        desired_moon_sep = min(desired_moon_sep, limit_sep)

    # if observer specified a moon separation from target, check it now
    # TODO: do we need to check this at the end of the exposure as well?
    # If so, then we may need to do it in observable() method
    if desired_moon_sep is not None:
        if ((c1.moon_sep < desired_moon_sep) or
            (c2.moon_sep < desired_moon_sep)):
            return False

    # moon looks good!
    return True


def check_slot(site, prev_slot, slot, ob, check_moon=True, check_env=True):
    """
    returns a Bunch saying whether ob will fit into slot
    """
    res = Bunch.Bunch(ob=ob, obs_ok=False, reason="No good reason!")

    # Check whether OB will fit in this slot
    delta = (slot.stop_time - slot.start_time).total_seconds()
    if ob.total_time > delta:
        res.setvals(obs_ok=False,
                    reason="Slot duration (%d) too short for OB (%d)" % (
            delta, ob.total_time))
        return res

    ## # check if instrument will be installed
    ## if not (ob.inscfg.insname in slot.data.instruments):
    ##     res.setvals(obs_ok=False, reason="Instrument '%s' not installed" % (
    ##         ob.inscfg.insname))
    ##     return res

    ## # check if filter will be installed
    ## if not (ob.inscfg.filter in slot.data.filters):
    ##     res.setvals(obs_ok=False, reason="Filter '%s' not installed [%s]" % (
    ##         ob.inscfg.filter, slot.data.filters))
    ##     return res

    ## # check if this slot can take this category
    ## if not ob.program.category in slot.data.categories:
    ##     res.setvals(obs_ok=False,
    ##                 reason="Slot cannot take category '%s'" % (
    ##         ob.program.category))
    ##     return res

    filterchange = False
    cur_filter = None
    filterchange_sec = 0.0

    # get immediately previous ob and filter
    if (prev_slot == None) or (prev_slot.ob == None):
        prev_ob = None
        cur_filter = slot.data.cur_filter

    else:
        prev_ob = prev_slot.ob
        cur_filter = prev_ob.inscfg.filter

    # calculate cost of filter exchange
    if cur_filter != ob.inscfg.filter:
        # filter exchange necessary
        filterchange = True
        filterchange_sec = ob.inscfg.calc_filter_change_time()
    #print "filter change time for new ob is %f sec" % (filterchange_sec)

    # for adding up total preparation time for new OB
    prep_sec = filterchange_sec

    # check dome status
    if slot.data.dome != ob.telcfg.dome:
        res.setvals(obs_ok=False, reason="Dome status OB(%s) != slot(%s)" % (
            ob.telcfg.dome, slot.data.dome))
        return res

    start_time = slot.start_time + timedelta(0, prep_sec)

    if slot.data.dome == 'closed':
        # <-- dome closed

        stop_time = start_time + timedelta(0, ob.total_time)

        # Check whether OB will fit in this slot
        if slot.stop_time < stop_time:
            res.setvals(obs_ok=False, reason="Not enough time in slot")
            return res

        res.setvals(obs_ok=True, prev_ob=prev_ob,
                    prep_sec=prep_sec, slew_sec=0.0, slew2_sec=0.0,
                    filterchange=filterchange,
                    filterchange_sec=filterchange_sec,
                    calibration_sec=0.0,
                    start_time=start_time, stop_time=stop_time,
                    delay_sec=0.0)
        return res

    # <-- dome open, need to check visibility and other criteria

    if check_env:
        # check seeing on the slot is acceptable to this ob
        if (slot.data.seeing > ob.envcfg.seeing):
            res.setvals(obs_ok=False,
                        reason="Seeing (%f > %f) not acceptable" % (
                slot.data.seeing, ob.envcfg.seeing))
            return res

        # check sky condition on the slot is acceptable to this ob
        if ob.envcfg.transparency is not None:
            if slot.data.transparency < ob.envcfg.transparency:
                res.setvals(obs_ok=False,
                            reason="Transparency (%f < %f) not acceptable" % (
                    slot.data.transparency, ob.envcfg.transparency))
                return res

    # Calculate cost of slew to this target
    # Assume that we want to do the SDSS calibration target first
    target = ob.target.sdss_calib
    if target is None:
        # No SDSS calibration target specified--go with OB main target
        target = ob.target

    if prev_ob is None:
        # no previous target--calculate cost from ...
        if slot.data.cur_az is not None:
            # ... current telescope position
            cur_alt_deg, cur_az_deg = slot.data.cur_el, slot.data.cur_az
        else:
            # ... parked position
            cur_alt_deg, cur_az_deg = parked_alt_deg, parked_az_deg

    else:
        # assume telescope is at previous target
        c0 = prev_ob.target.calc(site, start_time)
        cur_alt_deg, cur_az_deg = c0.alt_deg, c0.az_deg

    c1 = target.calc(site, start_time)

    slew_sec = calc_slew_time(cur_alt_deg, cur_az_deg, c1.alt_deg, c1.az_deg)
    #print("first slew time for new ob is %f sec" % (slew_sec))

    prep_sec += slew_sec
    # adjust on-target start time
    start_time += timedelta(0, slew_sec)

    # Is there an SDSS calibration target?  If so, then calculate in
    # calibration exposure and slew to main OB target
    calibration_sec = 0.0
    slew2_sec = 0.0
    if ob.target.sdss_calib is not None:

        calibration_sec = 30.0

        prep_sec += calibration_sec
        # adjust on-target start time
        start_time += timedelta(0, calibration_sec)

        # add slew time from calibration target to main target
        c2 = ob.target.calc(site, start_time)
        slew2_sec = calc_slew_time(c1.alt_deg, c1.az_deg, c2.alt_deg, c2.az_deg)
        #print("slew time from calib tgt to ob target is %f sec" % (slew2_sec))

        prep_sec += slew2_sec
        # adjust on-target start time
        start_time += timedelta(0, slew2_sec)

    # Check whether OB will fit in this slot
    ## delta = (slot.stop_time - start_time).total_seconds()
    ## if ob.total_time > delta:
    ##     return False

    min_el, max_el = ob.telcfg.get_el_minmax()

    # find the time that this object begins to be visible
    # TODO: figure out the best place to split the slot
    (obs_ok, t_start, t_stop) = site.observable(ob.target,
                                                start_time, slot.stop_time,
                                                min_el, max_el, ob.total_time,
                                                airmass=ob.envcfg.airmass,
                                                moon_sep=ob.envcfg.moon_sep)

    if not obs_ok:
        res.setvals(obs_ok=False,
                    reason="Time or visibility of target")
        return res

    # calculate delay until we could actually start observing the object
    # in this slot
    delay_sec = (t_start - start_time).total_seconds()

    stop_time = t_start + timedelta(0, ob.total_time)

    # check moon constraints between start and stop time
    if check_moon:
        obs_ok = check_moon_cond(site, start_time, stop_time, ob, res)
    else:
        obs_ok = True

    res.setvals(obs_ok=obs_ok, prev_ob=prev_ob,
                prep_sec=prep_sec, slew_sec=slew_sec,
                slew2_sec=slew2_sec,
                filterchange=filterchange,
                filterchange_sec=filterchange_sec,
                calibration_sec=calibration_sec,
                start_time=t_start, stop_time=stop_time,
                delay_sec=delay_sec)
    return res


def eval_schedule(schedule):
    """
    Calculates the number of filter exchanges and the time wasted in schedule
    """
    current_filter = None
    num_filter_exchanges = 0
    time_waste_sec = 0.0

    for slot in schedule.slots:
        ob = slot.ob
        # TODO: fix up a more solid check for delays
        if (ob == None) or ob.comment.startswith('Delay'):
            delta = (slot.stop_time - slot.start_time).total_seconds()
            time_waste_sec += delta
            continue

        if ((ob.inscfg.filter is not None) and
            (ob.inscfg.filter != current_filter)):
            num_filter_exchanges += 1
            current_filter = ob.inscfg.filter

    res = Bunch.Bunch(num_filter_exchanges=num_filter_exchanges,
                      time_waste_sec=time_waste_sec)
    return res


# END
