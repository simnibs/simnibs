import os
import sys
import shutil

def add_symlinks_or_cmd(external_progs,script_dir):
     ''' add symbolic links or .cmd '''
     for s in external_progs:
        if not os.path.exists(s):
            raise IOError('Could not find '+s)
        s = os.path.abspath(s)
        bash_name = os.path.join(script_dir, os.path.basename(s))
        if sys.platform == 'win32':
            bash_name=os.path.splitext(bash_name)[0] + '.cmd'
            print('making cmd link '+bash_name+' --> '+s)
            with open(bash_name, 'w') as f:
                f.write("@echo off\n")
                f.write(f'"{s}" %*')
        else:
            if os.path.lexists(bash_name):
                os.remove(bash_name)
            print('making sym link '+bash_name+' --> '+s)
            os.symlink(s, bash_name)

########################################################################################################
# external stuff for which symlinks or .cmd should be added to the scripts folder
########################################################################################################
external_progs = ['gmsh','meshfix']

bin_dir = os.path.join('simnibs', 'external', 'bin')
ending=''
if sys.platform == 'darwin':
    bin_dir = os.path.join(bin_dir, 'osx')
elif sys.platform == 'linux':
    bin_dir = os.path.join(bin_dir, 'linux')
elif sys.platform == 'win32':
    bin_dir = os.path.join(bin_dir, 'win')
    ending='.exe'
else:
    raise OSError('OS not supported!')
for i in range(len(external_progs)):
    external_progs[i] = os.path.join(bin_dir, external_progs[i]+ending)

if not sys.platform == 'win32':
    external_progs.append(os.path.join('simnibs','external','dwi2cond'))


script_dir = shutil.which('simnibs')
if script_dir is None:
    raise IOError('could not locate folder with console-scripts')
else:
    script_dir = os.path.dirname(script_dir)
    add_symlinks_or_cmd(external_progs,script_dir)
