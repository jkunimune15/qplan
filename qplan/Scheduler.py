#! /usr/bin/env python
#
# Scheduler.py -- Observing Queue Scheduler
#
#  Eric Jeschke (eric@naoj.org)
#
import os
import time
from datetime import timedelta
import pytz
import numpy as np
import StringIO

# 3rd party imports
from ginga.misc import Callback, Bunch
import astroplan
import astropy.units as u
import astropy.time as aptime

# maximum rank for a program
max_rank = 10.0

# time (sec) beyond which we breakout a slew into it's own OB
slew_breakout_limit = 30.0


class Scheduler(Callback.Callbacks):

    def __init__(self, logger):
        Callback.Callbacks.__init__(self)

        self.logger = logger

        self.site = astroplan.Observer.at_site("Subaru", timezone="US/Hawaii")

        # these are the main data structures used to schedule
        self.oblist = []
        self.schedule_recs = None
        self.programs = {}

        # FOR SCALING PURPOSES ONLY, define limits
        # (see cmp_res() )
        self.max_slew = 20.0*u.minute		# max slew
        self.max_rank = 10.0*u.minute		# max rank
        self.max_delay = 10.0*u.hour		# max wait for visibility
        self.min_delay = 30.0*u.minute		# min gap to try to place an OB
        self.max_filterchange = 35*u.minute  	# max filter exchange time
        self.alt_limits = (15.0*u.degree, 89.0*u.degree)

	# other telescope-related constants
        self.slew_rate = 0.5*u.degree/u.second	# rate of slew
        self.inst_reconfig_times = calc_reconfig_time()		# a mapping of variables to a mapping of pairs of states to times

        # define weights (see cmp_res() method)	TODO: Use these
        self.weights = Bunch.Bunch(w_rank=0.3, w_delay=0.2,
                                   w_slew=0.2, w_priority=0.1,
                                   w_filterchange = 0.3)

        # For callbacks
        for name in ('schedule-cleared', 'schedule-added', 'schedule-completed',):
            self.enable_callback(name)

    def set_weights(self, weights):
        self.weights = weights

    def set_programs_info(self, info):
        self.programs = {}
        for key, rec in info.items():
            if not rec.skip:
                self.programs[key] = rec

    def set_oblist_info(self, info):	#import some data for our oblist
        self.oblist = info
        #for ob in self.oblist:
        #    for c in ob.constraints:
        #        try:
        #            c.set_sdlr(self)	# assisgn yourself to the constraints while you're at it
        #        except AttributeError:
        #            continue

    def set_schedule_info(self, info):
        # Set our schedule_recs attribute to the supplied data
        # structure.
        self.schedule_recs = info


    def schedule_all(self):
	""" call the main scheduling algorithm and deduce stats about the new schedule
	"""
	# Get ready to time yourself
        t_t1 = time.time()
        
        # figure out the start and stop times
        start_str = self.schedule_recs[0].date+" "+self.schedule_recs[0].starttime
        end_str = self.schedule_recs[-1].date+" "+self.schedule_recs[-1].stoptime
        start_time = aptime.Time(start_str)
        end_time = aptime.Time(end_str) + 24*u.hour
        num_days = (end_time.to_datetime()-start_time.to_datetime()).days + 1
	
	# do some final assignments
        self.logger.info("scheduling %d OBs (from %d programs) for %d nights" % (
            len(self.oblist), len(self.programs), num_days))
        self.transitioner = astroplan.Transitioner(self.slew_rate, self.inst_reconfig_times)
        self.constraints = [astroplan.AtNightConstraint(), astroplan.AltitudeConstraint(*self.alt_limits)]

	# call the algorithm
        blk_lsts = self(start_time, end_time)
        self.schedule = blk_lsts[0]
        self.unscheduled_obs = blk_lsts[1]
        self.logger.info("schedule generated in %.2f sec" (t_elapsed))

	#TODO: delete this
        for ob in self.schedule:
            try:
                print ob.start_time,"%-12.12s"%ob.target.name,ob.end_time
            except Exception:
                print ob.start_time,"%-12.12s"%"Transition","I can't even"

        # build a lookup table of programs -> OBs
        props = {}
        total_program_time = 0
        for key in self.programs:
            total_time = self.programs[key].total_time
            props[key] = Bunch.Bunch(pgm=self.programs[key], obs=[],
                                     obcount=0, sched_time=0.0,
                                     total_time=total_time)
            total_program_time += total_time

        # count OBs in each program
        total_ob_time = 0
        for ob in self.oblist:
            pgmname = str(ob.program)
            props[pgmname].obs.append(ob)
            props[pgmname].obcount += 1
            # New policy is not to charge any overhead to the client,
            # including readout time
            obtime_no_overhead = ob.settings['exp_time'] * ob.settings['num_exp']
            total_ob_time += obtime_no_overhead

        # Note oversubscribed time
        self.logger.info("total program time=%d  total ob time=%d" % (
            total_program_time, total_ob_time))
        diff = total_ob_time - total_program_time
        if diff > 0:
            hrs = float(diff) / 3600.0
            self.logger.info("oversubscribed by %.2f hours" % (hrs))
        elif diff < 0:
            hrs = float(-diff) / 3600.0
            self.logger.info("undersubscribed by %.2f hours" % (hrs))

        # check time
        t_elapsed = time.time() - t_t1
        self.logger.info("%.2f sec to schedule all" % (t_elapsed))

        # print a summary
        out_f = StringIO.StringIO()
        num_obs = len(self.oblist)
        if num_obs > 0:
            pct = float(num_obs - len(unscheduled_obs)) / float(num_obs)
        else:
            pct = 0.
        out_f.write("%5.2f %% of OBs scheduled\n" % (pct*100.0))

	# check how many requested programs we completed
        completed, uncompleted = [], []
        for key in self.programs:
            bnch = props[key]
            for i in range(len(bnch.obs), 0, -1):
                ob = bnch.obs[i-1]
                if ob in unscheduled_obs:
                    if not bnch in uncompleted:
                        uncompleted.append(bnch)
                else:
                    bnch.obs.remove(ob)
            if not bnch in uncompleted:
                completed.append(bnch)

        completed = sorted(completed,
                           key=lambda bnch: max_rank - bnch.pgm.rank)
        uncompleted = sorted(uncompleted,
                             key=lambda bnch: max_rank - bnch.pgm.rank)
        self.make_callback('schedule-completed',
                           completed, uncompleted, self.schedule)
	
	# print out the completed programs
        out_f.write("Completed programs\n")
        for bnch in completed:
            out_f.write("%-12.12s   %5.2f  %d/%d  100%%\n" % (
                str(bnch.pgm), bnch.pgm.rank,
                bnch.obcount, bnch.obcount))
        out_f.write("\n")

	# calculate the amount of wasted time
        total_waste = 0*u.minute
	total_used = 0*u.minute
        for block in self.schedule:
            if type(block).__name__ == 'TransitionBlock':
                total_waste += (block.end_time-block.start_time).to(u.minute)
            else:
                total_used += (block.end_time-block.start_time).to(u.minute)

	# print out the uncompleted ones
        out_f.write("Uncompleted programs\n")
        for bnch in uncompleted:
            pct = float(bnch.obcount-len(bnch.obs)) / float(bnch.obcount) * 100.0
            uncompleted_s = "["+", ".join(map(lambda ob: ob.target.name, props[str(bnch.pgm)].obs))+"]"

            out_f.write("%-12.12s   %5.2f  %d/%d  %5.2f%%  %-142.142s\n" % (
                str(bnch.pgm), bnch.pgm.rank,
                bnch.obcount-len(bnch.obs), bnch.obcount, pct,
                uncompleted_s))
        out_f.write("\n")
        out_f.write("Total time: avail={} sched={} unsched={} min\n".format(
            (end_time-start_time).to(u.minute).value, total_used.value, total_waste.value))
        self.summary_report = out_f.getvalue()
        out_f.close()
        self.logger.info(self.summary_report)


    def select_schedule(self, schedule):
        self.selected_schedule = schedule
        self.make_callback('schedule-selected', schedule)


    def calc_transparency(self, times):	# finds the sky transparency at specified times
        """
        finds the sky transparency at specified times and place
        times is an astropy.Time object, place is an astroplan.site object, and return val is a list of floats
        """
        output = np.ones(times.size)*0.9	# initialize output to an array of zeros
        times.format = 'iso'
        for i, time in enumerate(times):
            for rec in self.schedule_recs:	# assume transparency changes on a daily basis
                if time.value[0:10] == rec.date[0:10]:
                    rec.data.transparency
        return output


    def __call__(self, start_time, end_time):
        """
        Compose a schedule made out of the OBs in self.oblist.
        Parameters
            start_time: an astropy.Time object indicating the beginning of this schedule
            end_time:   an astropy.Time object indicating the end of this schedule
        returns
            a list of two lists, the unused OBs, and the organized OBs and TBs.
	NOTE:
	    self.oblist will not be altered, nor will any attributes of self.
	    The OBs in the returned lists will not be copies, though, and the
	    OBs that make it into the final schedule will be altered
        """
        # start by combining universal constraints with OB-specific constraints
        for ob in self.oblist:
            if ob.constraints is None:
                ob._all_constraints = self.constraints
            else:
                ob._all_constraints = self.constraints + ob.constraints
            ob._time_scale = u.Quantity([0*u.minute, ob.duration/2, ob.duration])
        
        # define the main variables and start scheduling
        final_blocks = []
        unschedulable = []
        remaining_blocks = list(self.oblist)
        current_time = start_time
        while len(remaining_blocks) > 0 and current_time < end_time:
	    # score each potential ob based on how well it would fit here
            block_transitions = []
            block_scores = []
            for ob in reversed(remaining_blocks):
                # first calculate transition
                if len(final_blocks) > 0 and type(final_blocks[-1]) == astroplan.ObservingBlock:
                    tb = self.transitioner(final_blocks[-1], ob, current_time, self.site)
                    transition_time = tb.duration
                else:
                    tb = None
                    transition_time = 0*u.minute
                block_transitions.append(tb)
                
		# now verify that it is observable during this time
                times = current_time + transition_time + ob._time_scale
                if times[-1] <= end_time:
                    observable = astroplan.is_always_observable(ob._all_constraints, self.site,
                                                                [ob.target], times)
                    if observable[0]:
                        block_scores.append(ob.priority)
                    else:
                        block_scores.append(0)
		# if it would run over the end of the schedule, then assume it is unschedulable
                else:
                    unschedulable.append(ob)
                    block_transitions.pop()
                    remaining_blocks.remove(ob)
                    continue

	    # now that all that's been calculated, pick the best block
            best_block_idx = np.argmax(block_scores)
        
	    # if the best block is unobservable, then we obviously need a delay
            if block_scores[best_block_idx] == 0.:
                self.logger.info("nothing observable at %s" % (current_time))
                final_blocks.append(astroplan.TransitionBlock({'nothing_observable': self.min_delay}, current_time))
                current_time += self.min_delay
	    # otherwise, go ahead and add it to the schedule; don't forget the TransitionBlock!
            else:
                tb = block_transitions.pop(best_block_idx)
                ob = remaining_blocks.pop(best_block_idx)
                if tb is not None:
                    final_blocks.append(tb)
                    current_time += tb.duration
                ob.start_time = current_time
                current_time += ob.duration
                ob.end_time = current_time
                ob.constraints_value = block_scores[best_block_idx]
                final_blocks.append(ob)
        
        # return the scheduled blocks and the unscheduled blocks
        return [final_blocks, remaining_blocks+unschedulable]   


def eval_schedule(schedule):	# counts the number of filter changes and the wasted seconds

    current_filter = None
    num_filter_exchanges = 0
    time_waste_sec = 0.0

    for ob in schedule:
        if type(ob).__name__ == "TransitionBlock" and ob.components['AtNightConstraint']:
            delta = (ob.end_time-ob.start_time).to(u.second).value
            time_waste_sec += delta

        elif ((ob.settings['filter'] != None) and
            (ob.settings['filter'] != current_filter)):
            num_filter_exchanges += 1
            current_filter = ob.settings['filter']

    res = Bunch.Bunch(num_filter_exchanges=num_filter_exchanges,
                      time_waste_sec=time_waste_sec)
    return res


def calc_reconfig_time():	# builds a nested dictionary of reconfiguration times
    output = {'filter':{}}
    list_of_filters = ['g','r','i','z','Y','SH','nb656','nb515',None]
    for filter1 in list_of_filters:
        for filter2 in list_of_filters:
            if filter1 == filter2:
                output['filter'][(filter1, filter2)] = 0*u.hour
            else:
                output['filter'][(filter1, filter2)] = 1*u.hour
    return output

# END
