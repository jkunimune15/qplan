#
# ControlPanel.py -- Controller plugin for operating scheduler
#
# Eric Jeschke (eric@naoj.org)
#
import os.path

from ginga.misc import Bunch
from ginga.gw import Widgets

import PlBase
import filetypes

have_qdb = False
try:
    from qplan import q_db, q_query
    from Gen2.db.db_config import qdb_addr
    have_qdb = True

except ImportError:
    pass


class ControlPanel(PlBase.Plugin):

    def __init__(self, model, view, controller, logger):
        super(ControlPanel, self).__init__(model, view, controller, logger)

        self.input_dir = "inputs"

        self.weights_qf = None
        self.schedule_qf = None
        self.programs_qf = None
        self.ob_qf_dict = None
        self.tgtcfg_qf_dict = None
        self.envcfg_qf_dict = None
        self.inscfg_qf_dict = None
        self.telcfg_qf_dict = None

        self.qdb = None
        self.qa = None
        self.qq = None
        if have_qdb:
            self.connect_qdb()

    def connect_qdb(self):
        # Set up Queue database access
        self.qdb = q_db.QueueDatabase(self.logger, qdb_addr)
        self.qa = q_db.QueueAdapter(self.qdb)
        self.qq = q_query.QueueQuery(self.qa)

    def build_gui(self, container):
        vbox = Widgets.VBox()
        vbox.set_border_width(4)
        vbox.set_spacing(2)
        vbox.cfg_expand(8, 8)

        sw = Widgets.ScrollArea()
        sw.set_widget(vbox)

        fr = Widgets.Frame("Files")

        captions = (('Inputs:', 'label', 'Input dir', 'entry'),
                    ('Load Info', 'button'),
                    ('Build Schedule', 'button', 'Use QDB', 'checkbutton'))
        w, b = Widgets.build_info(captions, orientation='vertical')
        self.w = b

        b.input_dir.set_length(128)
        b.input_dir.set_text(self.controller.input_dir)
        b.load_info.add_callback('activated', self.initialize_model_cb)
        b.build_schedule.add_callback('activated', self.build_schedule_cb)

        fr.set_widget(w)
        vbox.add_widget(fr, stretch=0)

        spacer = Widgets.Label('')
        vbox.add_widget(spacer, stretch=1)

        hbox = Widgets.HBox()

        adj = Widgets.Slider(orientation='horizontal', track=True)
        adj.set_limits(0, 100, incr_value=1)
        idx = self.controller.idx_tgt_plots
        adj.set_value(idx)
        #adj.resize(200, -1)
        adj.set_tooltip("Choose subset of targets plotted")
        #self.w.plotgrp = adj
        adj.add_callback('value-changed', self.set_plot_pct_cb)
        hbox.add_widget(adj, stretch=1)

        sb = Widgets.SpinBox()
        sb.set_limits(1, 100)
        num = self.controller.num_tgt_plots
        sb.set_value(num)
        sb.set_tooltip("Adjust size of subset of targets plotted")
        sb.add_callback('value-changed', self.set_plot_limit_cb)
        hbox.add_widget(sb, stretch=0)

        vbox.add_widget(hbox, stretch=0)

        ## btns = Widgets.HBox()
        ## btns.set_spacing(3)

        ## btn = Widgets.Button("Close")
        ## #btn.add_callback('activated', lambda w: self.close())
        ## btns.add_widget(btn, stretch=0)
        ## btns.add_widget(Widgets.Label(''), stretch=1)
        ## vbox.add_widget(btns, stretch=0)

        self.sw = sw
        container.add_widget(sw, stretch=1)

    def initialize_model_cb(self, widget):
        self.input_dir = self.w.input_dir.get_text().strip()
        self.input_fmt = self.controller.input_fmt

        try:
            # read weights
            self.weights_qf = filetypes.WeightsFile(self.input_dir, self.logger, file_ext=self.input_fmt)
            # Load "Weights" Tab
            if 'weightstab' not in self.view.plugins:
                self.view.load_plugin('weightstab', 'WeightsTab', 'WeightsTab', 'report', 'Weights')
            self.model.set_weights_qf(self.weights_qf)

            # read schedule
            self.schedule_qf = filetypes.ScheduleFile(self.input_dir, self.logger, file_ext=self.input_fmt)
            # Load "Schedule" Tab
            if 'scheduletab' not in self.view.plugins:
                self.view.load_plugin('scheduletab', 'ScheduleTab', 'ScheduleTab', 'report', 'Schedule')
            self.model.set_schedule_qf(self.schedule_qf)

            # read proposals
            self.programs_qf = filetypes.ProgramsFile(self.input_dir, self.logger, file_ext=self.input_fmt)
            if 'programstab' not in self.view.plugins:
                self.view.load_plugin('programstab', 'ProgramsTab', 'ProgramsTab', 'report', 'Programs')
            self.model.set_programs_qf(self.programs_qf)

            # read observing blocks
            self.ob_qf_dict = {}
            self.tgtcfg_qf_dict = {}
            self.envcfg_qf_dict = {}
            self.inscfg_qf_dict = {}
            self.telcfg_qf_dict = {}
            self.oblist_info = []

            propnames = list(self.programs_qf.programs_info.keys())
            propnames.sort()

            for propname in propnames:
                pf = filetypes.ProgramFile(self.input_dir, self.logger, propname, self.programs_qf.programs_info, file_ext=self.input_fmt)

                # Set telcfg
                telcfg_qf = pf.cfg['telcfg']
                self.telcfg_qf_dict[propname] = telcfg_qf
                self.model.set_telcfg_qf_dict(self.telcfg_qf_dict)

                # Set inscfg
                inscfg_qf = pf.cfg['inscfg']
                self.inscfg_qf_dict[propname] = inscfg_qf
                self.model.set_inscfg_qf_dict(self.inscfg_qf_dict)

                # Set envcfg
                envcfg_qf = pf.cfg['envcfg']
                self.envcfg_qf_dict[propname] = envcfg_qf
                self.model.set_envcfg_qf_dict(self.envcfg_qf_dict)

                # Set targets
                tgtcfg_qf = pf.cfg['targets']
                self.tgtcfg_qf_dict[propname] = tgtcfg_qf
                self.model.set_tgtcfg_qf_dict(self.tgtcfg_qf_dict)

                # Finally, set OBs
                self.ob_qf_dict[propname] = pf.cfg['ob']
                #self.oblist_info.extend(self.oblist[propname].obs_info)
                self.model.set_ob_qf_dict(self.ob_qf_dict)

        except Exception as e:
            self.logger.error("Error initializing: %s" % (str(e)))


    def update_scheduler(self, use_db=False):
        sdlr = self.model.get_scheduler()
        try:
            sdlr.set_weights(self.weights_qf.weights)
            sdlr.set_schedule_info(self.schedule_qf.schedule_info)
            pgms = self.programs_qf.programs_info
            sdlr.set_programs_info(pgms)

            # TODO: this maybe should be done in the Model
            ob_keys = set([])
            propnames = list(self.programs_qf.programs_info.keys())
            okprops = []
            ob_dict = {}
            for propname in propnames:
                if pgms[propname].skip:
                    continue
                okprops.append(propname)

                # get all OB keys for this program
                for ob in self.ob_qf_dict[propname].obs_info:
                    key = (propname, ob.target.name)
                    ob_keys.add(key)
                    ob_dict[key] = ob

            self.logger.info("%s OBs after excluding skipped programs." % (
                len(ob_keys)))

            # TODO: remove keys for OBs that are already executed
            if use_db:
                do_not_execute = set(self.qq.get_do_not_execute_ob_keys())
                ob_keys -= do_not_execute
                self.logger.info("%s OBs after removing executed OBs." % (
                    len(ob_keys)))

            elif self.model.completed_keys is not None:
                do_not_execute = set(self.model.completed_keys)
                ob_keys -= do_not_execute
                self.logger.info("%s OBs after removing executed OBs." % (
                    len(ob_keys)))

            # for a deterministic result
            ob_keys = list(ob_keys)
            ob_keys.sort()

            # Now turn the keys back into actual OB list
            oblist_info = []
            for key in ob_keys:
                oblist_info.append(ob_dict[key])

            self.oblist_info = oblist_info

            # TODO: only needed if we ADD or REMOVE programs

            sdlr.set_oblist_info(self.oblist_info)

        except Exception as e:
            self.logger.error("Error storing into scheduler: %s" % (str(e)))

        self.logger.info("Scheduler initialized")

    def build_schedule_cb(self, widget):
        # update the model with any changes from GUI
        use_db = self.w.use_qdb.get_state()
        self.update_scheduler(use_db=use_db)

        sdlr = self.model.get_scheduler()
        self.view.nongui_do(sdlr.schedule_all)

    def set_plot_pct_cb(self, w, val):
        #print(('pct', val))
        self.controller.idx_tgt_plots = val
        self.model.select_schedule(self.model.selected_schedule)
        return True

    def set_plot_limit_cb(self, w, val):
        #print(('limit', val))
        self.controller.num_tgt_plots = val
        self.model.select_schedule(self.model.selected_schedule)
        return True

#END
