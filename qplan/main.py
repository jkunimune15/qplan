#
# main.py -- Queue planner/scheduler main program logic
#
# eric@naoj.org
#
# stdlib imports
import sys, os
import threading
import logging

# Subaru python stdlib imports
from ginga.misc import ModuleManager, Settings, Task, Bunch
from ginga.misc import log
import ginga.toolkit as ginga_toolkit

# Local application imports
from .Control import Controller
from .Model import QueueModel
from qplan import version

moduleHome = os.path.split(sys.modules['qplan.version'].__file__)[0]
sys.path.insert(0, moduleHome)
pluginHome = os.path.join(moduleHome, 'plugins')
sys.path.insert(0, pluginHome)

default_layout = ['seq', {},
                   ['vbox', dict(name='top', width=1440, height=900),
                    dict(row=['hbox', dict(name='menu')],
                         stretch=0),
                    dict(row=['hpanel', {},
                     ['ws', dict(name='left', width=100, show_tabs=False),
                      # (tabname, layout), ...
                      ],
                     ['vpanel', dict(width=700),
                      ['hpanel', dict(height=400),
                       ['vbox', dict(name='main', width=700),
                        dict(row=['ws', dict(name='report', group=1)], stretch=1)],
                       ['ws', dict(name='right', width=350, group=2),
                        # (tabname, layout), ...
                        ],
                       ],
                      ['hpanel', {},
                       ['ws', dict(name='sub1', width=700, height=520,
                                   group=1)],
                       ['ws', dict(name='sub2', width=500, group=1)],
                       ],
                      ],
                     ], stretch=1),
                    dict(row=['hbox', dict(name='status')], stretch=0),
                    ]]


plugins = [
    # pluginName, moduleName, className, workspaceName, tabName
    ('slewchart', 'SlewChart', 'SlewChart', 'sub2', 'Slew Chart'),
    ('airmasschart', 'AirMassChart', 'AirMassChart', 'sub1', 'AirMass Chart'),
    ('schedule', 'Schedule', 'Schedule', 'left', 'Schedule'),
    ('report', 'Report', 'Report', 'report', 'Report'),
    ('logger', 'Logger', 'Logger', 'report', 'Log'),
    ('cp', 'ControlPanel', 'ControlPanel', 'right', 'Control Panel'),
    ('night_activity', 'SumChart', 'NightSumChart', 'sub1', 'Night Activity Chart'),
    ('night_sched', 'SumChart', 'SchedSumChart', 'sub1', 'Schedules Chart'),
    ('proposals', 'SumChart', 'ProposalSumChart', 'sub1', 'Proposals Chart'),
    ('semester', 'SumChart', 'SemesterSumChart', 'sub1', 'Semester Chart'),
    ]


class QueuePlanner(object):
    """
    This class exists solely to be able to customize the queue planner
    startup/application.
    """
    def __init__(self, layout=default_layout):
        self.plugins = []
        self.layout = layout

    def add_plugins(self, plugins):
        self.plugins.extend(plugins)

    def add_default_options(self, optprs):
        """
        Adds the default reference viewer startup options to an
        OptionParser instance `optprs`.
        """
        optprs.add_option("-c", "--completed", dest="completed", default=None,
                          metavar="FILE",
                          help="Specify FILE of completed OB keys")
        optprs.add_option("--date-start", dest="date_start", default=None,
                          help="Define the start of the schedule ('YYYY-MM-DD HH:MM')")
        optprs.add_option("--date-stop", dest="date_stop", default=None,
                          help="Define the end of the schedule ('YYYY-MM-DD HH:MM')")
        optprs.add_option("--debug", dest="debug", default=False, action="store_true",
                          help="Enter the pdb debugger on main()")
        optprs.add_option("--display", dest="display", metavar="HOST:N",
                          help="Use X display on HOST:N")
        optprs.add_option("-g", "--geometry", dest="geometry",
                          metavar="GEOM", default=None,
                          help="X geometry for initial size and placement")
        optprs.add_option("-i", "--input", dest="input_dir", default=".",
                          metavar="DIRECTORY",
                          help="Read input files from DIRECTORY")
        optprs.add_option("-f", "--format", dest="input_fmt", default=None,
                          metavar="FILE_FORMAT",
                          help="Specify input file format (csv, xls, or xlsx)")
        ## optprs.add_option("--modules", dest="modules", metavar="NAMES",
        ##                   help="Specify additional modules to load")
        ## optprs.add_option("--monitor", dest="monitor", metavar="NAME",
        ##                   default='monitor',
        ##                   help="Synchronize from monitor named NAME")
        ## optprs.add_option("--monchannels", dest="monchannels",
        ##                   default='status', metavar="NAMES",
        ##                   help="Specify monitor channels to subscribe to")
        ## optprs.add_option("--monport", dest="monport", type="int",
        ##                   help="Register monitor using PORT", metavar="PORT")
        optprs.add_option("--numthreads", dest="numthreads", type="int",
                          default=30,
                          help="Start NUM threads in thread pool", metavar="NUM")
        optprs.add_option("-o", "--output", dest="output_dir", default="output",
                          metavar="DIRECTORY",
                          help="Write output files to DIRECTORY")
        ## optprs.add_option("--plugins", dest="plugins", metavar="NAMES",
        ##                   help="Specify additional plugins to load")
        ## optprs.add_option("--port", dest="port", type="int", default=None,
        ##                   help="Register using PORT", metavar="PORT")
        optprs.add_option("--profile", dest="profile", action="store_true",
                          default=False,
                          help="Run the profiler on main()")
        optprs.add_option("-t", "--toolkit", dest="toolkit", metavar="NAME",
                          default='qt4',
                          help="Prefer GUI toolkit (gtk|qt)")
        log.addlogopts(optprs)


    def main(self, options, args):
        # Create top level logger.
        svcname = 'qplan'
        logger = log.get_logger(name=svcname, options=options)

        ev_quit = threading.Event()

        thread_pool = Task.ThreadPool(logger=logger, ev_quit=ev_quit,
                                     numthreads=options.numthreads)

        ginga_toolkit.use(options.toolkit)

        from qplan.View import Viewer
        # must import AFTER Viewer
        from ginga.Control import GuiLogHandler

        class QueuePlanner(Controller, Viewer):

            def __init__(self, logger, thread_pool, module_manager, preferences,
                         ev_quit, model):

                Viewer.__init__(self, logger, ev_quit)
                Controller.__init__(self, logger, thread_pool, module_manager,
                                    preferences, ev_quit, model)

        # Get settings folder
        ## if os.environ.has_key('CONFHOME'):
        ##     basedir = os.path.join(os.environ['CONFHOME'], svcname)
        ## else:
        basedir = os.path.join(os.environ['HOME'], '.' + svcname)
        if not os.path.exists(basedir):
            os.mkdir(basedir)
        prefs = Settings.Preferences(basefolder=basedir, logger=logger)

        mm = ModuleManager.ModuleManager(logger)

        ## # Add any custom modules
        ## if options.modules:
        ##     modules = options.modules.split(',')
        ##     for mdlname in modules:
        ##         #self.mm.loadModule(name, pfx=pluginconfpfx)
        ##         self.mm.loadModule(name)

        model = QueueModel(logger=logger)

        if options.completed is not None:
            # specify a list of completed OB keys
            with open(options.completed, 'r') as in_f:
                buf = in_f.read()
            import ast
            model.completed_keys = ast.literal_eval(buf)

        # Start up the control/display engine
        qplanner = QueuePlanner(logger, thread_pool, mm, prefs, ev_quit, model)
        qplanner.set_input_dir(options.input_dir)
        qplanner.set_input_fmt(options.input_fmt)

        # Build desired layout
        qplanner.build_toplevel(default_layout)
        for w in qplanner.ds.toplevels:
            w.show()

        # load plugins
        for pluginName, moduleName, className, wsName, tabName in self.plugins:
            qplanner.load_plugin(pluginName, moduleName, className,
                                 wsName, tabName)

        guiHdlr = GuiLogHandler(qplanner)
        #guiHdlr.setLevel(options.loglevel)
        guiHdlr.setLevel(logging.INFO)
        fmt = logging.Formatter(log.LOG_FORMAT)
        guiHdlr.setFormatter(fmt)
        logger.addHandler(guiHdlr)

        qplanner.update_pending()

        # Did user specify geometry
        if options.geometry:
            qplanner.set_geometry(options.geometry)

        # Raise window
        w = qplanner.w.root
        w.show()

        server_started = False

        # Create threadpool and start it
        try:
            # Startup monitor threadpool
            thread_pool.startall(wait=True)

            try:
                # if there is a network component, start it
                if hasattr(qplanner, 'start'):
                    task = Task.FuncTask2(qplanner.start)
                    thread_pool.addTask(task)

                # Main loop to handle GUI events
                qplanner.mainloop(timeout=0.001)

            except KeyboardInterrupt:
                logger.error("Received keyboard interrupt!")

        finally:
            logger.info("Shutting down...")
            thread_pool.stopall(wait=True)

        sys.exit(0)


def planner(sys_argv):

    viewer = QueuePlanner(layout=default_layout)
    viewer.add_plugins(plugins)

    # Parse command line options with optparse module
    from optparse import OptionParser

    usage = "usage: %prog [options] cmd [args]"
    optprs = OptionParser(usage=usage,
                          version=('%%prog %s' % version.version))
    viewer.add_default_options(optprs)

    (options, args) = optprs.parse_args(sys_argv[1:])

    if options.display:
        os.environ['DISPLAY'] = options.display

    # Are we debugging this?
    if options.debug:
        import pdb

        pdb.run('viewer.main(options, args)')

    # Are we profiling this?
    elif options.profile:
        import profile

        print(("%s profile:" % sys_argv[0]))
        profile.run('viewer.main(options, args)')

    else:
        viewer.main(options, args)


# END
