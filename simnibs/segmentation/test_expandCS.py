import configparser
from glob import glob
import logging
import nibabel as nib
import os
import shutil
import timeit
import copy

from simnibs.utils.simnibs_logger import logger
import simnibs.segmentation.brain_surface as sf
from simnibs import mesh_io

subject_dir=""
config_file="surface_reco.ini" # content would be added to charm.ini somewhere in the simnibs install folder
ref_fs="C:/Users/axthi/simnibs/simnibs_examples/ernie/m2m_ernie/ref_FS.nii.gz"
# ref_fs is needed to set correct header in freesurfer meshes in expandCS for debugging only;
# it can likely be dropped unless its also needed for the gifti headers


# =============================================================================
# start logging
# =============================================================================

logger.setLevel(logging.DEBUG) 

# log to console and file
#logger.addHandler(sf.make_handler("stream", logging.INFO))
#    f_log.write('<HTML><HEAD><TITLE>headreco report</TITLE></HEAD><BODY><pre>')

# log command line call
# sf.log(" ".join(argv), level=logging.INFO)



# =============================================================================
# read config file
# =============================================================================
logger.info("Reading settings from "+config_file)
shutil.copy(config_file, os.path.join(subject_dir, "charm_settings.ini"))
cfg_raw = configparser.ConfigParser()
try:
    cfg_raw.read(config_file)
except:
    logger.error("error when reading "+config_file)
    raise IOError('error whenreading config file')
    
# convert settings from strings to bool, int, float, string, ...
# Note: the entry bla = 1 will be converted to int, bla = 1.0 to float!
config={}
for section in cfg_raw.sections():
    config[section]={}
    for key in cfg_raw[section]:
        try:
            config[section][key]=eval(cfg_raw[section][key])
        except:
            logger.error("error when formatting "+section+" "+key+": "+cfg_raw[section][key])
            raise IOError('error when parsing config file')
del cfg_raw           
        


# =============================================================================
# call of expandCS --> to be integrated into cat_surf_createCS_recoSurf.m
# =============================================================================

# get relevant settings for surface creation
opt=config['surface_reco']
opt['surffolder'] = os.path.join(subject_dir,"surf")
                                                            
actualsurf = 'lh' # 'lh' or 'rh'
Ppial = os.path.join(
    opt['surffolder'],actualsurf+".pial.gii") # has to be replaced by a .gii

# load surface and thickness data
central = glob(os.path.join(opt['surffolder'], actualsurf+".central.gii"))[0]
thickness = glob(os.path.join(opt['surffolder'], actualsurf+".thickness"))[0]

central = mesh_io.read_gifti_surface(central)
cv = central.nodes[:]
cf = central.elm[:, :3] - 1
th = mesh_io.read_curv(thickness)

# expand and save as .gii
tic = timeit.default_timer()
logger.info(actualsurf+": Expanding ...")
cv_pial = sf.expandCS(cv, cf, th/2, actualsurf=actualsurf, ref_fs=ref_fs) # ref_fs only for debugging
toc = timeit.default_timer()
logger.info(actualsurf+": Expanding done after "+str(round(toc-tic))+" seconds")
pial = copy.deepcopy(central)
pial.nodes.node_coord = cv_pial
mesh_io.write_gifti_surface(pial, Ppial) # replace by saving gifti


# =============================================================================
# stop logging
# =============================================================================

# end HTML file
#with open(logfile, 'a') as f_log:
#    f_log.write('</pre></BODY></HTML>')

logging.shutdown()
