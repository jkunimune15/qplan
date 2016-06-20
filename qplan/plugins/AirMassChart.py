#
# AirMassChart.py -- AirMass chart plugin
#
# Eric Jeschke (eric@naoj.org)
#
from __future__ import print_function
from datetime import timedelta
#import pytz

from ginga.gw import Widgets, Plot
from ginga.misc import Bunch

import PlBase
from plots.airmass import AirMassPlot


class AirMassChart(PlBase.Plugin):

    def __init__(self, model, view, controller, logger):
        super(AirMassChart, self).__init__(model, view, controller, logger)

        self.schedules = {}
        self.initialized = False

        # Set preferred timezone for plot
        #self.tz = pytz.utc

        sdlr = model.get_scheduler()
        self.tz = sdlr.timezone
        sdlr.add_callback('schedule-cleared', self.clear_schedule_cb)
        sdlr.add_callback('schedule-added', self.new_schedule_cb)

        model.add_callback('schedule-selected', self.show_schedule_cb)

    def build_gui(self, container):

        self.plot = AirMassPlot(700, 500, logger=self.logger)

        self.canvas = Plot.PlotWidget(self.plot, width=700, height=500)

        container.set_margins(2, 2, 2, 2)
        container.set_spacing(4)

        container.add_widget(self.canvas, stretch=1)

    def show_schedule_cb(self, qmodel, schedule):
        try:
            info = self.schedules[schedule]

            if not self.initialized:
                self.plot.setup()
                self.initialized = True

            if info.num_tgts == 0:
                self.logger.debug("no targets for plotting airmass")
                self.view.gui_do(self.plot.clear)
            else:
                self.logger.debug("plotting airmass")
                self.view.gui_do(self.plot.clear)
                site = info.site
                target_data = info.target_data

                # Plot a subset of the targets
                idx = int((self.controller.idx_tgt_plots / 100.0) * len(target_data))
                num_tgts = self.controller.num_tgt_plots
                target_data = target_data[idx:idx+num_tgts]

                self.view.gui_do(self.plot.plot_altitude, site, target_data,
                                 self.tz)
            self.view.gui_do(self.plot.draw)
        ## except KeyError:
        ##     pass
        except Exception as e:
            self.logger.error("Error plotting airmass: %s" % (str(e)))

        return True

    def add_schedule(self, schedule):
        self.logger.debug("adding schedule %s" % (schedule))

        start_time = schedule.start_time
        sdlr = self.model.get_scheduler()
        t = start_time.astimezone(sdlr.timezone)
        # if schedule starts after midnight, change start date to the
        # day before, this is due to the way the Observer module charts
        # airmass
        if 0 <= t.hour < 12:
            t -= timedelta(0, 3600*12)
        ndate = t.strftime("%Y/%m/%d")

        targets = []
        site = sdlr.site

        for slot in schedule.slots:
            ob = slot.ob
            if (ob is not None) and (not ob.derived):
                # not an OB generated to serve another OB
                # TODO: make sure targets are unique in pointing
                targets.append(ob.target)

        # make airmass plot
        num_tgts = len(targets)
        target_data = []
        lengths = []
        if num_tgts > 0:
            for tgt in targets:
                info_list = get_target_info(site,tgt,(schedule.start_time+schedule.stop_time)/2.0)
                target_data.append(Bunch.Bunch(history=info_list, target=tgt))
                lengths.append(len(info_list))

        # clip all arrays to same length
        min_len = 0
        if len(lengths) > 0:
            min_len = min(*lengths)
        for il in target_data:
            il.history = il.history[:min_len]

        self.schedules[schedule] = Bunch.Bunch(site=site, num_tgts=num_tgts,
                                               target_data=target_data)

    def new_schedule_cb(self, qscheduler, schedule):
        self.add_schedule(schedule)

    def clear_schedule_cb(self, qscheduler):
        #self.view.gui_do(self.plot.clear)
        self.logger.info("cleared plot")


def get_target_info(site, target, date,
                    time_interval=5):
    """Compute various values for a target from sunrise to sunset"""

    def _round_time(dtime):
        # Sets time to nice rounded value
        y, m ,d, hh, mm, ss = dtime.tuple()
        mm = mm - (mm % 5)
        return ephem.Date(datetime(y, m , d, hh, mm, 5, 0))

    def _set_data_range(sunset, sunrise, tint):
        # Returns numpy array of dates ranging from 15 minutes before sunset
        # to 15 minutes after sunrise
        ss = _round_time(ephem.Date(sunset - 15*ephem.minute))
        sr = _round_time(ephem.Date(sunrise + 15*ephem.minute))
        return numpy.arange(ss, sr, tint)

    # default for start time is sunset on the current date
    time_start = site.sun_set_time(date, u'previous')
    # default for stop time is sunrise on the current date
    time_stop = site.sun_rise_time(date, u'next')

    t_range = _set_data_range(time_start, time_stop,
                                  time_interval*ephem.minute)
    #print('computing airmass history...')
    history = []

    # TODO: this should probably return a generator
    for ut in t_range:
        # ugh
        tup = ephem.Date(ut).tuple()
        args = tup[:-1] + (int(tup[-1]),)
        ut_with_tz = datetime(*args,
                                  tzinfo=self.tz_utc)
        info = target.calc(self, ut_with_tz)
        history.append(info)
    #print(('computed airmass history', self.history))
    return history


#END
