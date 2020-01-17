import configparser
from glob import glob
import logging
import nibabel as nib
import os
import shutil
import surface_functions as sf
import timeit

subject_dir=""
config_file="surface_reco.ini" # content would be added to charm.ini somewhere in the simnibs install folder
ref_fs="C:/Users/axthi/simnibs/simnibs_examples/ernie/m2m_ernie/ref_FS.nii.gz"
# ref_fs is needed to set correct header in freesurfer meshes in expandCS for debugging only;
# it can likely be dropped unless its also needed for the gifti headers


# =============================================================================
# start logging
# =============================================================================

# logging
logging.captureWarnings(True)
logger = logging.getLogger("py.warnings")
# use logging.INFO for normal use; logging.DEBUG for development:
logger.setLevel(logging.DEBUG) 

# log to console and file
logfile = os.path.join(subject_dir, "charm_log.html")
logger.addHandler(sf.make_handler("stream", logging.INFO))
logger.addHandler(sf.make_handler("file", logging.DEBUG, logfile ))
with open(logfile, 'a') as f_log:
    f_log.write('<HTML><HEAD><TITLE>headreco report</TITLE></HEAD><BODY><pre>')

# log command line call
# sf.log(" ".join(argv), level=logging.INFO)



# =============================================================================
# read config file
# =============================================================================
sf.log("Reading settings from "+config_file, level=logging.INFO)
shutil.copy(config_file, os.path.join(subject_dir, "charm_settings.ini"))
cfg_raw = configparser.ConfigParser()
try:
    cfg_raw.read(config_file)
except:
    sf.log("error when reading "+config_file, level=logging.ERROR)
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
            sf.log("error when formatting "+section+" "+key+": "+cfg_raw[section][key], level=logging.ERROR)
            raise IOError('error when parsing config file')
del cfg_raw           
        


# =============================================================================
# call of expandCS --> to be integrated into cat_surf_createCS_recoSurf.m
# =============================================================================

# get relevant settings for surface creation
opt=config['surface_reco']
opt['surffolder'] = os.path.join(subject_dir,"surf")
                                                            
actualsurf='lh' # 'lh' or 'rh'
Ppial=os.path.join(opt['surffolder'],actualsurf+".pial.stl") # has to be replaced by a .gii

# load surface and thickness data
central = glob(os.path.join(opt['surffolder'], actualsurf+".central.gii"))[0]
thickness = glob(os.path.join(opt['surffolder'], actualsurf+".thickness"))[0]

central = nib.load(central)
cv = central.darrays[0].data
cf = central.darrays[1].data
th = sf.read_curv(thickness)

# expand and save as .gii
tic = timeit.default_timer()
sf.log(actualsurf+": Expanding ...", level=logging.INFO)
cv_pial = sf.expandCS(cv,cf,th/2, actualsurf=actualsurf, ref_fs=ref_fs) # ref_fs only for debugging
toc = timeit.default_timer()
sf.log(actualsurf+": Expanding done after "+str(round(toc-tic))+" seconds", level=logging.INFO)
sf.mesh_save(cv_pial, cf, Ppial) # replace by saving gifti


# =============================================================================
# stop logging
# =============================================================================

# end HTML file
with open(logfile, 'a') as f_log:
    f_log.write('</pre></BODY></HTML>')

# remove handlers, otherwise I get duplicate logging lines when calling the script
# repeatedaly in the IDLE
while logger.hasHandlers():
    logger.handlers[0].close()
    logger.removeHandler(logger.handlers[0])
logging.shutdown()


