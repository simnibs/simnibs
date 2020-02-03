# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 16:27:07 2020

@author: axthi
"""

import configparser
from simnibs.utils.simnibs_logger import logger


def read_ini(config_file=None):
    """ read .ini-file and return its content as dictionary
    
    NOTE: This function is potentially UNSAFE, as it allows to
          enter adversarial code via the .ini-file. It should be
          updated as some point.
    """
    logger.info("Reading settings from "+config_file)
    cfg = configparser.ConfigParser()
    fn=cfg.read(config_file)
    if not len(fn):
        logger.error("error when reading "+config_file)
        raise IOError('error when reading config file')
    
    # convert settings from strings to bool, int, float, string, ...
    # Note: the entry bla = 1 will be converted to int, bla = 1.0 to float!
    config={}
    for section in cfg.sections():
        config[section]={}
        for key in cfg[section]:
            try:
                config[section][key]=eval(cfg[section][key])
            except:
                logger.error("error when formatting "+section+" "+key+": "+cfg[section][key])
                raise IOError('error when parsing config file')
    del cfg
    return config

