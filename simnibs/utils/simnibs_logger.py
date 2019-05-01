import logging
import sys
import warnings
import traceback
from multiprocessing import Lock

global logger
logger = logging.getLogger('simnibs')
sh = logging.StreamHandler()
formatter = logging.Formatter('[ %(name)s ]%(levelname)s: %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)
logger.setLevel(logging.INFO)
logging.addLevelName(25, 'SUMMARY')


def log_warnings(message, category, filename, lineno, file=None, line=None):
    logger.warn(warnings.formatwarning(message, category, filename, lineno))


warnings.showwarning = log_warnings


# Redirect the exceptions to the logger
def register_handler(orig_excepthook=sys.excepthook):
    def log_excep(type, value, tback):
        """Log uncaught exceptions. When an exception occurs, sys.exc_info()
        returns a tuple of three variables (exception class, exception value,
        traceback). Setting
            sys.excepthook = log_excep
        will replace the standard way of handling exceptions but that of log_excep.
        log_excep takes the sys.exc_info() as input and prints the exception to 
        "logger" at level error.
        """
        #logger.debug("Traceback:", exc_info=(type, value, tback))
        logger.critical("Unhandled exception:", exc_info=(type, value, tback))
        #orig_excepthook(*exc_info)
    sys.excepthook = log_excep

register_handler()



class PicklebleFileHandler(logging.FileHandler):
    def __getstate__(self):
        d = dict(self.__dict__)
        d['stream'] = None
        return d

if __name__ == '__main__':
    warnings.warn('aaaa')
