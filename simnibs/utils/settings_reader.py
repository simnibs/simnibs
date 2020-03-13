# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 16:27:07 2020

@author: axthi
"""

import configparser
from simnibs.utils.simnibs_logger import logger


def convert_str(s):
    """ try to convert string to integer, float or boolean;
        if this doesn't work, return the string
    """    
    try:
        val = int(s)
    except:
        try:
            val = float(s)
        except:
            if s.lower()== 'true': val=True
            elif s.lower()== 'false': val=False
            else: val = s
    return val


def split_str(s):
    """ splits a string delimited by brackets at commas
        and returns a list of sub-strings
    """    
    s=s[1:-1]
    in_quotes=False
    in_dquotes=False
    in_brackets=0
    idx_comma=-1
    idx_substr=[]
    for i in range(len(s)):
        if s[i] == '\'': in_quotes = not in_quotes
        if s[i] == '\"': in_dquotes = not in_dquotes
        if s[i] == '[' and not (in_quotes or in_dquotes):
            in_brackets += 1
        if s[i] == ']' and not (in_quotes or in_dquotes):
            in_brackets -= 1 
        if s[i] == ',' and not (in_quotes or in_dquotes or in_brackets != 0):
            idx_substr.append([min(idx_comma+1,len(s)), i])
            idx_comma=i
            
    if idx_comma<len(s):
        idx_substr.append([min(idx_comma+1,len(s)), len(s)])
    
    if in_quotes or in_dquotes or in_brackets != 0:
        logger.error("error when parsing "+s)
        raise IOError('error when parsing config file')
    
    s_split=[]
    for i in range(len(idx_substr)):
        s_split.append(s[idx_substr[i][0]:idx_substr[i][1]])
    return s_split


def parse_str(s):
    """ convert string into list by locating pairs of [ and ]
        and splitting the content based on the commas;
        convert content from string to integer or float if possible
    """ 
    s = s.strip()
    if not len(s): return s
    if (s[0] == '\'') and (s[-1] == '\''):
        s = s[1:-1]
        if not len(s): return s
    if (s[0] == '\"') and (s[-1] == '\"'): 
        s = s[1:-1]
        if not len(s): return s
    if (s[0] == '[') and (s[-1] == ']'):
        # split s and parse subparts separately
        s=split_str(s)
        for i in range(len(s)):
            s[i] = parse_str(s[i])
    else:
        s = convert_str(s)
    return s


def read_ini(config_file=None):
    """ read .ini-file and return its content as dictionary
        supported variable types are bool, int, float and string,
        as well lists [var1,var2,var3]
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
                config[section][key]=parse_str(cfg[section][key])
            except:
                logger.error("error when formatting "+section+" "+key+": "+cfg[section][key])
                raise IOError('error when parsing config file')
    del cfg
    return config

