import logging
import sys
import warnings

global logger
logger = logging.getLogger('simnibs')
sh = logging.StreamHandler()
formatter = logging.Formatter('[ %(name)s ]%(levelname)s: %(message)s')
sh.setFormatter(formatter)
sh.setLevel(logging.INFO)
logger.addHandler(sh)
logger.setLevel(logging.DEBUG)
logging.addLevelName(25, 'SUMMARY')


def log_warnings(message, category, filename, lineno, file=None, line=None):
    logger.warn(warnings.formatwarning(message, category, filename, lineno))


warnings.showwarning = log_warnings


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
