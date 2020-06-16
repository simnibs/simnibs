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
        output=os.path.join(pack_dir, 'simnibs_env.zip'),
        compress_level=0,
        force=True
    )
    shutil.unpack_archive(
        os.path.join(pack_dir, 'simnibs_env.zip'),
        os.path.join(pack_dir, 'simnibs_env'),
    )
    os.remove(os.path.join(pack_dir, 'simnibs_env.zip'))
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

    # Copy documentation
    shutil.copytree('../docs/build/html', os.path.join(pack_dir, 'documentation'))

    # Create bash or bat file for installation
    if sys.platform == 'win32':
        shutil.copy('../simnibs/resources/gui_icon.ico', os.path.join(pack_dir, 'gui_icon.ico'))
        fn_script = os.path.join(pack_dir, 'installer.nsi')
        with open('installer.nsi', 'r') as f:
            install_script = Template(f.read()).render(
                version='.'.join(version.split('.')[:2]),
                full_version=version
            )
        with open(fn_script, 'w') as f:
            f.write(install_script)
        print('Creating ')
        subprocess.run(
            fr'"%programfiles(x86)%\NSIS\makensis.exe" {fn_script}',
            check=True,
            shell=True
        )
    if sys.platform == 'darwin':
        print('TODO')
    elif sys.platform=='linux':
        fn_script = os.path.join(pack_dir, 'install')
        with open('install', 'r') as f:
            install_script = Template(f.read()).render(
                version='.'.join(version.split('.')[:2])
            )
        with open(fn_script, 'w') as f:
            f.write(install_script)
        os.chmod(
            fn_script,
            os.stat(fn_script).st_mode |
            stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
        )
        # zip the whole thing
        shutil.make_archive(f'simnibs-{version}-{os_name}', 'zip', pack_dir)

if __name__ == '__main__':
    build()
