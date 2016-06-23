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
import numpy
import StringIO

# 3rd party imports
from ginga.misc import Callback, Bunch
import astropy.units as u
import astropy.time as aptime
import astroplan

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

        # define weights (see cmp_res() method)
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

    def set_oblist_info(self, info):
        self.oblist = info

    def set_schedule_info(self, info):
        # Set our schedule_recs attribute to the supplied data
        # structure.
        self.schedule_recs = info


    def schedule_all(self):
	""" The main method of Scheduler: take self.oblist and sort it into usable schedule,
	which is saved in self.schedule
	"""
	# Get ready to time yourself
        t_t1 = time.time()
        
        # figure out the start and stop times
        start_str = self.schedule_recs[0].date+"T"+self.schedule_recs[0].starttime
        stop_str = self.schedule_recs[-1].date+"T"+self.schedule_recs[-1].stoptime
        start_time = aptime.Time(start_str, format='isot')
        stop_time = aptime.Time(stop_str, format='isot') + 24*u.hour
        num_days = (stop_time.to_datetime()-start_time.to_datetime()).days + 1
	
	# call the astroplan scheduling algorithm
        self.logger.info("preparing to schedule {}".format(map(lambda r: r.date, self.schedule_recs)))
        transitioner = astroplan.scheduling.Transitioner(self.slew_rate, self.inst_reconfig_times)
        constraints = [astroplan.AtNightConstraint(), astroplan.AltitudeConstraint(*self.alt_limits)]
        astroSdlr = astroplan.PriorityScheduler(start_time, stop_time, constraints, self.site,
                                                transitioner, self.min_delay, self.max_slew)
        self.schedule = astroSdlr(self.oblist)	# TODO: make this sequential

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

        self.logger.info("scheduling %d OBs (from %d programs) for %d nights" % (
            len(self.oblist), len(self.programs), num_days))

        # check time
        t_elapsed = time.time() - t_t1
        self.logger.info("%.2f sec to schedule all" % (t_elapsed))

        # print a summary
        out_f = StringIO.StringIO()
        num_obs = len(self.oblist)
        pct = 0.0

        # count all the obs that didn't make it into the final schedule
        unscheduled_obs = []
        for proposed_ob in self.oblist:
            scheduled = False
            for chosen_ob in self.schedule:
                if chosen_ob.target.name == proposed_ob.target.name:
                    scheduled = True
                    break
            if not scheduled:
                unscheduled_obs.append(proposed_ob)    
        if num_obs > 0:
            pct = float(num_obs - len(unscheduled_obs)) / float(num_obs)
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
            (start_time-stop_time).to(u.minute).value, total_used.value, total_waste.value))
        self.summary_report = out_f.getvalue()
        out_f.close()
        self.logger.info(self.summary_report)


    def select_schedule(self, schedule):
        self.selected_schedule = schedule
        self.make_callback('schedule-selected', schedule)


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
