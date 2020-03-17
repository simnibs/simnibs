# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 16:27:07 2020

@author: axthi
"""

import configparser
import json
from simnibs.utils.simnibs_logger import logger


def read_ini(config_file):
    """ read .ini-file and return its content as dictionary
        left-hand side should be in JSON style
        int = 1
        float = 1.2
        bool = true/false (lower case)
        str = "string" (double qotes)
        list = [1, 2, 3]
        dict = {"a": 1, "b": 2} (key must be string)

    Parameters
    -------------
    config_file: str
        Name of configuration file

    Returns
    ---------
    config: dict
        Nested dictionary with parsed values for each section
    """
    logger.info("Reading settings from "+config_file)
    cfg = configparser.ConfigParser()
    fn = cfg.read(config_file)
    if not len(fn):
        logger.error("error when reading " + config_file)
        raise IOError('error when reading config file')

    # convert settings from strings to bool, int, float, string, ...
    # Note: the entry bla = 1 will be converted to int, bla = 1.0 to float!
    config = {}
    for section in cfg.sections():
        config[section] = {}
        for key in cfg[section]:
            try:
                config[section][key] = json.loads(cfg[section][key])
            except json.decoder.JSONDecodeError:
                logger.error(
                    f"error when formatting {section} {key}="
                    f"{cfg[section][key]}"
                )
                raise IOError('error when parsing config file')

    return config
