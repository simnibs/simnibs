import sys
import os
import shutil
import subprocess



def _get_default_dir():
    if sys.platform == 'win32':
        return os.path.join(os.environ['LOCALAPPDATA'], 'SimNIBS')
    elif sys.platform == 'linux':
       return os.path.join(os.environ['HOME'], 'SimNIBS')
    elif sys.platform == 'darwin':
       return os.path.join(os.environ['HOME'], 'Applications', 'SimNIBS')
    else:
        raise OSError('OS not supported')


def install(prefix):
    # Copy the simnibs_env folder 
    # pip-install simnibs from the new simnibs_env folder
    # run postinstall

def distribute():
    import conda_pack
    # Create a new environment
    subprocess.run(f'conda env create -y --name simnibs_env_tmp --clone {os.environ['CONDA_PREFIX']}', shell=True)
    subprocess.run('conda activate simnibs_env_tmp && pip -y uninstall simnibs', shell=True)
    # Remove existing SimNIBS install before??
    conda_pack.pack(prefix=os.environ['simnibs_env_tmp'], dest_prefix='simnibs_env')
    subprocess.run(f'conda env remove -y --name simnibs_env_tmp', shell=True)
