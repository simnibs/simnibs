''' Has to be run from the conda root environment!!
'''
import stat
import glob
import sys
import argparse
import os
import shutil
import subprocess


def build():
    import conda_pack
    version = open("simnibs/_version.py").readlines()[-1].split()[-1].strip("\"'")
    # Create a new environment
    if os.path.dirname(__file__):
        os.chdir(os.path.dirname(__file__))
    if not os.path.isdir('pack'):
        os.mkdir('pack')
    env = os.path.join(os.path.dirname(__file__), 'environment_{}.yml')
    if sys.platform == 'linux':
        env = env.format('linux')
        out_pack = os.path.join('pack', 'simnibs_env.tar.gz')
        os_name = 'linux'
    elif sys.platform == 'darwin':
        env = env.format('macOS')
        out_pack = os.path.join('pack', 'simnibs_env.tar.gz')
        os_name = 'macOS'
    elif sys.platform == 'win32':
        env = env.format('win')
        out_pack = os.path.join('pack', 'simnibs_env.zip')
        os_name = 'win'
    else:
        raise OSError('OS not supported!')
    # Create temporary environment
    subprocess.run(
        f'conda env create -n simnibs_env_tmp -f {env}',
        check=True,
        shell=True
    )
    # Pack
    if os.path.isfile(out_pack):
        os.remove(out_pack)
    conda_pack.pack(
        name='simnibs_env_tmp',
        dest_prefix='simnibs_env',
        output=out_pack,
    )
    # Remove temporary env
    subprocess.run(
        'conda env remove -y --name simnibs_env_tmp',
        check=True,
        shell=True
    )

    # Copy wheel
    wheels = glob.glob(f'dist/simnibs-{version}*.whl')
    if len(wheels) == 0:
        raise FileNotFoundError(f'Did not find any wheels for version {version}')
    for f in wheels:
        shutil.copy(f, 'pack')

    # Create bash or bat file for installation
    if sys.platform == 'win32':
        with open(os.path.join('pack/install.cmd'), 'w') as f:
            f.write(f'SET INSTALL_DIR=%LOCALAPPDATA%\SimNIBS\n')
            f.write('mkdir "%INSTALL_DIR%\simnibs_env"\n')
            f.write('powershell.exe -nologo -noprofile -command "& '
                    '{ Add-Type -A \'System.IO.Compression.FileSystem\'; '
                    '[IO.Compression.ZipFile]::ExtractToDirectory(\'%~dp0simnibs_env.zip\', \'%INSTALL_DIR%\simnibs_env\'); }"\n'
            )
            f.write('call "%INSTALL_DIR%\simnibs_env\\Scripts\\activate"\n')
            f.write('python -m pip install simnibs --no-cache-dir --no-index --upgrade --find-links=./\n')
            f.write('postinstall_simnibs -d "%INSTALL_DIR%" --copy-matlab --setup-links --no-extra-coils')

    else:
        if sys.platform == 'darwin':
            fn_script = os.path.join('pack', 'install')
            install_dir = "$HOME/Applications/SimNIBS"
        else:
            fn_script = os.path.join('pack', 'install')
            install_dir = "$HOME/SimNIBS"
        with open(fn_script, 'w') as f:
            # Some problem with Qt
            f.write('#! /bin/bash -e\n')
            f.write(f'INSTALL_DIR={install_dir}\n')
            f.write('CUR_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)"\n')
            f.write('mkdir -p "$INSTALL_DIR/simnibs_env"\n')
            f.write('tar -zxf "$CUR_DIR/simnibs_env.tar.gz" -C "$INSTALL_DIR/simnibs_env"\n')
            f.write('source "$INSTALL_DIR/simnibs_env/bin/activate"\n')
            f.write('python -m pip install simnibs --no-cache-dir --no-index --upgrade --find-links="$CUR_DIR"\n')
            f.write('python -m pip install pyqt5 --force-reinstall --no-cache-dir\n') # I need to re-install pyqt, this requires internet
            f.write('postinstall_simnibs -d "$INSTALL_DIR" --copy-matlab --setup-links --no-extra-coils')

        os.chmod(fn_script,
                 os.stat(fn_script).st_mode |
                 stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    # zip the whole thing
    shutil.make_archive(f'simnibs-{version}-{os_name}', 'zip', 'pack')


def _get_default_dir():
    if sys.platform == 'win32':
        return os.path.join(os.environ['LOCALAPPDATA'], 'SimNIBS')
    elif sys.platform == 'linux':
       return os.path.join(os.environ['HOME'], 'SimNIBS')
    elif sys.platform == 'darwin':
       return os.path.join(os.environ['HOME'], 'Applications', 'SimNIBS')
    else:
        raise OSError('OS not supported')


if __name__ == '__main__':
    build()
