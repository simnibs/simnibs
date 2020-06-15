''' Has to be run from the conda root environment!!
'''
import stat
import glob
import sys
import argparse
import os
import subprocess
import shutil

from jinja2 import Template


def build():
    import conda_pack
    version = open("../simnibs/_version.py").readlines()[-1].split()[-1].strip("\"'")
    pack_dir = os.path.abspath('pack')
    # Create a new environment
    if os.path.dirname(__file__):
        os.chdir(os.path.dirname(__file__))
    if not os.path.isdir(pack_dir):
        os.mkdir(pack_dir)
    if sys.platform == 'linux':
        os_name = 'linux'
    elif sys.platform == 'darwin':
        os_name = 'macOS'
    elif sys.platform == 'win32':
        os_name = 'win'
    else:
        raise OSError('OS not supported!')
    # Create temporary environment
    env = os.path.join(
        os.path.dirname(__file__), '..', f'environment_{os_name}.yml'
    )
    subprocess.run(
        f'conda env create -n simnibs_env_tmp -f {env}',
        check=True,
        shell=True
    )
    # Pack
    conda_pack.pack(
        name='simnibs_env_tmp',
        dest_prefix='simnibs_env',
        output=os.path.join(pack_dir, 'simnibs_env.tar.gz'),
        compress_level=0,
        force=True
    )
    shutil.unpack_archive(
        os.path.join(pack_dir, 'simnibs_env.tar.gz'),
        os.path.join(pack_dir, 'simnibs_env'),
    )
    os.remove(os.path.join(pack_dir, 'simnibs_env.tar.gz'))
    # Remove temporary env
    subprocess.run(
        'conda env remove -y --name simnibs_env_tmp',
        check=True,
        shell=True
    )

    # Copy wheel
    wheels = glob.glob(f'../dist/simnibs-{version}*.whl')
    if len(wheels) == 0:
        raise FileNotFoundError(f'Did not find any wheels for version {version}')
    for f in wheels:
        shutil.copy(f, pack_dir)

    # Create bash or bat file for installation
    if sys.platform == 'win32':
        print('TODO')
        '''
        with open(os.path.join('pack/install.cmd'), 'w') as f:
            f.write(f'SET INSTALL_DIR=%USERPROFILE%\SimNIBS\n')
            f.write('mkdir "%INSTALL_DIR%\simnibs_env"\n')
            f.write('powershell.exe -nologo -noprofile -command "& '
                    '{ Add-Type -A \'System.IO.Compression.FileSystem\'; '
                    '[IO.Compression.ZipFile]::ExtractToDirectory(\'%~dp0simnibs_env.zip\', \'%INSTALL_DIR%\simnibs_env\'); }"\n'
            )
            f.write('call "%INSTALL_DIR%\simnibs_env\\Scripts\\activate"\n')
            f.write('python -m pip install simnibs --no-cache-dir --no-index --upgrade --find-links=./\n')
            f.write('postinstall_simnibs -d "%INSTALL_DIR%" --copy-matlab --setup-links --no-extra-coils')
        '''
    else:
        if sys.platform == 'darwin':
            print('TODO')
            '''
            fn_script = os.path.join('pack', 'install')
            install_dir = "$HOME/Applications/SimNIBS"
            '''
        else:
            fn_script = os.path.join(pack_dir, 'install')
            install_dir = "$HOME/SimNIBS"
            with open('install', 'r') as f:
                install_script = Template(f.read()).render(
                    version='.'.join(version.split('.')[:2])
                )
            with open(fn_script, 'w') as f:
                f.write(install_script)

            os.chmod(fn_script,
                     os.stat(fn_script).st_mode |
                     stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    # zip the whole thing
    shutil.make_archive(f'simnibs-{version}-{os_name}', 'zip', pack_dir)

if __name__ == '__main__':
    build()
