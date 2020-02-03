# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 17:12:24 2020

@author: axthi
"""


import logging
import nibabel as nib
import os
import shutil
import time
from simnibs import utils
from simnibs import SIMNIBSDIR
from simnibs.utils.simnibs_logger import logger



def view(subject_dir):
    print('charm viewer not yet implemented, sorry...')



def run(subject_dir=None, T1=None, T2=None,
        registerT2=False, initatlas=False, segment=False, mesh=False,
        skipregisterT2=False, usesettings=None, options_str=None):
    """charm pipeline
    
    PARAMETERS
    ----------
    subject_dir : str, mandatory
        output directory
    T1 : str
        filename of T1 image
    T2 : str
        filename of T2 image
    
    --> parameters to control the workflow:
    registerT2 : bool
        run T2-to-T1 registration (default = False)
    initatlas : bool
        run affine registration of atlas to input images (default = False)
    segment : bool
        run volume and surface segmentation (default = False)
    mesh : bool
        run tetrahedral meshing (default = False)
        
    --> further parameters: 
    skipregisterT2 : bool
        skip T2-to-T1 registration, just copy T2-weighted image (default = False)
    usesettings : str
        filename of alternative settings-file (default = None)
    options_str : str
        string of command line options to add to logging (default = None)
        
    RETURNS
    ----------
        None
    """  
    

    # ------------------------START UP-----------------------------------------
    start = time.time()
        
    # make subject_dir if not existent
    if not os.path.exists(subject_dir):
            os.mkdir(subject_dir)
    
    # start logging ...
    logfile = os.path.join(subject_dir, "charm_log.html")
    with open(logfile, 'a') as f:
        f.write('<HTML><HEAD><TITLE>charm report</TITLE></HEAD><BODY><pre>')
        f.close()
    fh = logging.FileHandler(logfile, mode='a')
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    fh.setFormatter(formatter)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    utils.simnibs_logger.register_excepthook(logger)
        
    if options_str is not None:
        logger.debug('options: '+options_str)
    logger.info('charm run started: '+time.asctime())
     
    # get subID
    subID=os.path.basename(subject_dir)
    idx=[i for i in range(len(subID)) if subID[i] == '_']
    if not len(idx) or idx[0]+1==len(subID):
        raise RuntimeError("ERROR: subID could not be determined from subject directory:  "+subject_dir)
    subID=subID[idx[0]+1:]
    
    # copy T1 (as nii.gz) if supplied
    if T1 is not None:
        if not os.path.exists(T1):
            raise FileNotFoundError(T1)
        if len(T1)>7 and T1[-7:].lower()=='.nii.gz':
            shutil.copyfile(T1,os.path.join(subject_dir, 'T1.nii.gz'))
        else:
            nib.save(nib.load(T1),os.path.join(subject_dir, 'T1.nii.gz'))
    
    if skipregisterT2 and T2 is not None:
        # skip T2-to-T1 registration, just copy T2 image (as nii.gz)
        registerT2=False
        if not os.path.exists(T2):
            raise FileNotFoundError(T2)
        if len(T2)>7 and T2[-7:].lower()=='.nii.gz':
            shutil.copyfile(T2,os.path.join(subject_dir, 'T2_reg.nii.gz'))
        else:
            nib.save(nib.load(T2),os.path.join(subject_dir, 'T2_reg.nii.gz'))
    
    # read settings and copy settings file
    if usesettings is None:
        fn_settings=os.path.join(SIMNIBSDIR,'charm.ini')
    else:
        fn_settings=usesettings
    settings=utils.settings_reader.read_ini(fn_settings)
    shutil.copyfile(fn_settings,os.path.join(subject_dir, 'settings.ini'))    


    # -------------------------PIPELINE STEPS----------------------------------
    if registerT2:
    # register T2 to T1
        logger.info('starting registerT2')
        
        # get local settings
        local_settings=settings['registerT2']
        
        # check input files
        T1=os.path.join(subject_dir, 'T1.nii.gz')
        if not os.path.exists(T1):
            raise FileNotFoundError(T1)       
        if not os.path.exists(T2):
            raise FileNotFoundError(T2)
            
        # do your stuff
            
        # write QA results


    if initatlas:
    # initial affine registration of atlas to input images, including break neck
        logger.info('starting initatlas')
        
        
    if segment:
    # nonlinear registration of atlas, reconstruct cortical surfaces, 
    # register to fsaverage and MNI space, create label image
        logger.info('starting segment')
        
        
    if mesh:
    # create mesh from label image
        logger.info('starting mesh')

    
    
    # -------------------------TIDY UP-----------------------------------------
        
    # log stopping time and total duration ...
    logger.info('charm run finished: '+time.asctime())
    logger.info('Total running time: '+utils.simnibs_logger.format_time(time.time()-start))
    
    # stop logging ...
    while logger.hasHandlers():
        logger.removeHandler(logger.handlers[0])
    utils.simnibs_logger.unregister_excepthook()
    logging.shutdown()
    with open(logfile, 'a') as f:
        f.write('</pre></BODY></HTML>')
        f.close()
    
    