import logging
import sys
import warnings
import numpy as np

global logger
logger = logging.getLogger('simnibs')
sh = logging.StreamHandler()
formatter = logging.Formatter('[ %(name)s ] %(levelname)s: %(message)s')
sh.setFormatter(formatter)
sh.setLevel(logging.INFO)
logger.addHandler(sh)
logger.setLevel(logging.DEBUG)
logging.addLevelName(25, 'SUMMARY')
logging.addLevelName(26, 'SUMMARY')

def log_warnings(message, category, filename, lineno, file=None, line=None):
    logger.warn(warnings.formatwarning(message, category, filename, lineno))


# This is causing errors in pytest
#warnings.showwarning = log_warnings


def register_excepthook(logger):
    def log_excep(exc_type, exc_value, exc_traceback):
        if issubclass(exc_type, KeyboardInterrupt):
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return
        logger.debug(
            "Traceback",
            exc_info=(exc_type, exc_value, exc_traceback)
        )
        logger.critical(
            "Uncaught exception",
            exc_info=(exc_type, exc_value, None)
        )
    sys.excepthook = log_excep


def unregister_excepthook():
    sys.excepthook = sys.__excepthook__


def format_time(running_time):
    """Format time in seconds as hours:minutes:seconds.
    
    PARAMETERS
    ----------
    running_time : float
        Time in seconds.
    
    RETURNS
    ----------
    running_time : str
        The time formatted as hours:minutes:seconds.
    """
    hrs = np.uint16(np.floor(running_time/(60.**2)))
    mts = np.uint16(np.floor(running_time/60.-hrs*60))
    sec = np.uint16(np.round(running_time-hrs*60.**2-mts*60.))

    return "{:02d}:{:02d}:{:02d}".format(hrs,mts,sec)