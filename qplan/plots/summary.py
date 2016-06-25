#
# summary.py -- creates summary plots
#
# Russell Kackley (rkackley@naoj.org)
#
import numpy as np
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
import matplotlib.patches as mpatches

from ginga.util import plots
from ginga.misc import Bunch

import astropy.units as u

from qsim import eval_schedule

class BaseSumPlot(plots.Plot):
    def __init__(self, width, height, logger=None):
        # create matplotlib figure
        super(BaseSumPlot, self).__init__(width=width, height=height,
                                          logger=logger)
        self.logger = logger
        self.barWidth = 0.5
        self.legendFont = FontProperties('sans-serif', 'normal', 'normal', 'normal', 'normal', 'small')
        self.bottomMargin = 0.15
        self.heightSpace = 0.6 # space between plots
        self.fig.subplots_adjust(bottom=self.bottomMargin, hspace=self.heightSpace)
        self.grades = ('A', 'B', 'C', 'F')
        self.grade_colors = {
            'A': 'green',
            'B': 'violet',
            'C': 'coral',
            'F': 'aqua',
            }
        self.activity_colors = {
            'Long slew':     'blue',
            'Filter change': 'cyan',
            'Delay':         'orchid',
            'Science':       'green',
            'Unscheduled':   'darkred'
            }

        self.ob_types = ('Long slew', 'Filter change', 'Delay', 'Science', 'Unscheduled')

    def clear(self):
        self.fig.clf()

class NightSumPlot(BaseSumPlot):
    # Make a bar chart to show which types of OB's will be executed
    # during the night

    def plot(self, full_sched):
        # Create a plot that shows the activities (e.g., slews, filter
        # changes, etc.) for all the nights.
        plt = self.fig.add_subplot(111)
        plt.set_title('Nightly Activity')
        plt.set_xlabel('Minutes from start of night')
        
        schedules = each_night(full_sched)	# split full_sched into nights

        # Iterate through all the dates in the schedules list. Note
        # that we iterate in reverse order so that the oldest dates
        # appear at the top of the plot, which is the same order as is
        # shown in the "Schedule" section on the left-hand side of
        # qplan/qexec.
        date_list = []
        for i, schedule in enumerate(list(reversed(schedules))):
            if len(schedule) <= 0:
                continue
            date_list.append(start_of(schedule).strftime('%Y-%m-%d'))
            y = [i]
            previous_slot_right = np.array([0.0])
            for ob in schedule:
                dt = ob.end_time.to_datetime() - ob.start_time.to_datetime()
                dt_minutes = dt.total_seconds() / 60.0
                width = np.array([dt_minutes])
                if type(ob).__name__ == 'ObservingBlock':
                    ob_type = 'Science'
                elif 'slew' in ob.component.keys():
                    ob_type = 'Long slew'
                elif 'filter' in ob.component.keys():
                    ob_type = 'Filter change'
                else:
                    ob_type = 'Delay'

                bar = plt.barh(y, width, self.barWidth, left=previous_slot_right, color=self.activity_colors[ob_type])
                previous_slot_right += width

        # Add the y-axis titles, which are the dates from the
        # schedules list.
        y = np.arange(len(schedules))
        plt.set_yticks(y+self.barWidth/2.)
        plt.set_yticklabels(date_list)

        # Reduce the plot area by a little bit so that we have room
        # for the legend outside the plot.
        box = plt.get_position()
        plt.set_position([box.x0, box.y0, box.width * 0.9, box.height])

        # Create some matplotlib "Patches" so that we can use them in
        # the legend
        legend_patches = []
        legend_titles = []
        for ob_type in self.ob_types:
            legend_patches.append(mpatches.Patch(color=self.activity_colors[ob_type]))
            legend_titles.append(ob_type)

        # Add a legend to the plot. We put the legend outside the plot
        # area so that we don't obscure any of the bars.
        plt.legend(legend_patches, legend_titles, prop=self.legendFont, loc='center left', bbox_to_anchor=(1, 0.5), handlelength=1)

        self.draw()

class ProposalSumPlot(BaseSumPlot):
    # Make a bar chart to show the completed OB percentage for each
    # proposal.
    def plot(self, completed, uncompleted):
        plt = self.fig.add_subplot(111)
        plt.set_title('Proposal Completion Percentage')
        plt.set_ylabel('OB Percent Complete')
        ind = np.arange(len(completed)+len(uncompleted))
        propID_comp_percent = {}
        grades_dict = {}
        for grade in self.grades:
            grades_dict[grade] = []
        for i, proposal in enumerate(completed+uncompleted):
            #x = [i]
            total_ob_count = float(proposal.obcount)
            uncompleted_count = len(proposal.obs)
            completed_count = float(total_ob_count - uncompleted_count)
            propID = str(proposal.pgm)
            completed_pct = 0.0
            if total_ob_count > 0:
                completed_pct = completed_count / total_ob_count
            propID_comp_percent[propID] = completed_pct * 100.0
            if propID not in grades_dict[proposal.pgm.grade]:
                grades_dict[proposal.pgm.grade].append(propID)

        # For the bar chart, we want all proposals grouped into their
        # "grade" category and then, within that category, sort the
        # proposals by their proposal ID.
        propid_list = []
        comp_percent = []
        colors = []
        for grade in self.grades:
            for propID in sorted(grades_dict[grade]):
                propid_list.append(propID)
                comp_percent.append(propID_comp_percent[propID])
                colors.append(self.grade_colors[grade])

        plt.bar(ind, comp_percent, self.barWidth, color=colors)
        plt.set_xticks(ind+self.barWidth/2.)
        plt.set_xticklabels(propid_list, rotation=45, ha='right')

        # Create some matplotlib "Patches" so that we can use them in
        # the legend
        legend_patches = []
        legend_titles = []
        for grade in self.grades:
            legend_patches.append(mpatches.Patch(color=self.grade_colors[grade]))
            legend_titles.append(grade)

        plt.legend(legend_patches, legend_titles, prop=self.legendFont, title='Grades', loc='center left', bbox_to_anchor=(1, 0.5), handlelength=1)

        self.draw()

class ScheduleSumPlot(BaseSumPlot):
    # Makes a bar chart to show scheduled/unscheduled minutes for each
    # night
    def plot(self, full_sched):
        schedules = each_night(full_sched)
        plt = self.fig.add_subplot(111)
        plt.set_title('Nightly Schedules')
        plt.set_ylabel('Time (min)')
        ind = np.arange(len(schedules))
        date_list = []
        sched_minutes = []
        unsched_minutes = []
        for schedule in schedules:
            if len(schedule) <= 0:
                continue
            date_list.append(start_of(schedule).strftime('%Y-%m-%d'))
            time_avail = length_of(schedule)
            time_avail_minutes = time_avail.total_seconds() / 60.0
            time_waste_minutes = eval_schedule(schedule).time_waste_sec / 60.0
            sched_minutes.append(time_avail_minutes - time_waste_minutes)
            unsched_minutes.append(time_waste_minutes)
        self.logger.debug('ind %s' % ind)
        self.logger.debug('date_list %s' % date_list)
        self.logger.debug('sched_minutes %s' % sched_minutes)
        self.logger.debug('unsched_minutes %s' % unsched_minutes)
        try:
            sched_bar = plt.bar(ind, sched_minutes, self.barWidth, color='g')
            unsched_bar = plt.bar(ind, unsched_minutes, self.barWidth, color='darkred', bottom=sched_minutes)
            plt.set_xticks(ind+self.barWidth/2.)
            plt.set_xticklabels(date_list, rotation=45, ha='right')
            plt.legend((unsched_bar, sched_bar), ('Delay+Unscheduled', 'Scheduled'), prop=self.legendFont)
 
            self.draw()
        except Exception:
            self.logger.error('could not produce bar chart due to lack of data')

class SemesterSumPlot(BaseSumPlot):
    # Makes a pie chart to show percentage of available time allocated
    # to each proposal and also the unscheduled time.
    def plot(self, full_sched):
        plt = self.fig.add_subplot(111)
        total_time_avail = 0.
        total_time_waste = 0.
        propID_alloc_minutes = {}
        grades_dict = {}
        for grade in self.grades:
            grades_dict[grade] = []
        time_avail = length_of(full_sched)
        total_time_avail = time_avail.total_seconds() / 60.0
        total_time_waste = eval_schedule(full_sched).time_waste_sec / 60.0

        for ob in full_sched:
            if type(ob).__name__ != "TransitionBlock":
                propID = str(ob.program)
                if propID not in grades_dict[ob.program.grade]:
                    grades_dict[ob.program.grade].append(propID)
                if propID in propID_alloc_minutes:
                    propID_alloc_minutes[propID] += (ob.end_time-ob.start_time).to(u.minute).value
                else:
                    propID_alloc_minutes[propID] = (ob.end_time-ob.start_time).to(u.minute).value

        total_time_sched = total_time_avail - total_time_waste
        self.logger.debug('propID_alloc_minutes %s' % propID_alloc_minutes)
        self.logger.debug('total_time_sched %f' % total_time_sched)
        self.logger.debug('total_time_waste %f' % total_time_waste)

        # For the pie chart, we want all proposals grouped into their
        # "grade" category and then, within that category, sort the
        # proposals by their proposal ID.
        labels = []
        sizes = []
        colors = []
        for grade in self.grades:
            for propID in sorted(grades_dict[grade]):
                labels.append(propID)
                sizes.append(propID_alloc_minutes[propID])
                colors.append(self.grade_colors[grade])

        labels.append('Unscheduled')
        sizes.append(total_time_waste)
        colors.append('darkred')
        plt.pie(sizes, labels=labels, colors=colors, autopct='%1.0f%%', shadow=True)
        if '-' in labels[0]:
            # TODO: title if nothing can be scheduled
            semester, ident = labels[0].split('-')
            plt.set_title('Total for Semester %s = %5.0f Hours' % (semester, total_time_avail / 60.0))

        # Create some matplotlib "Patches" so that we can use them in
        # the legend
        legend_patches = []
        legend_titles = []
        for grade in self.grades:
            legend_patches.append(mpatches.Patch(color=self.grade_colors[grade]))
            legend_titles.append(grade)

        plt.legend(legend_patches, legend_titles, prop=self.legendFont, title='Grades', loc='center left', bbox_to_anchor=(1, 0.5), handlelength=1)

        self.draw()


def each_night(schedule):	# breaks schedule up into several schedules, separated by daytime
    schedules = []		# schedules: a list of lists of obs, each inner list being a night
    night_start = 0
    for i, block in enumerate(schedule):		# parse through the schedule
        if type(block).__name__ == "TransitionBlock" and 'day' in block.components.keys():	# if there is a (day) transitionblock at i
            if night_start < i:					# split the list,
                schedules.append(schedule[night_start:i])	# add the piece to schedules
            night_start = i+1					# and move on
    if night_start < len(schedule):
        schedules.append(schedule[night_start:i])
    return schedules


def start_of(schedule):	#calculates the start time of a list of OBs; returns a datetime.time
    if len(schedule) <= 0:
        return None
    return schedule[0].start_time.to_datetime()


def end_of(schedule):	#calculates the end time of a list of OBs; returns a datetime.time
    if len(schedule) <= 0:
        return None
    return schedule[-1].end_time.to_datetime()


def length_of(schedule):#calculates the duration of a list of OBs; returns a datetime.timedelta
    return end_of(schedule) - start_of(schedule)

