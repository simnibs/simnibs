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
    env_prefix = os.path.join(pack_dir, 'simnibs_env_tmp')
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
    # Install requirements
    subprocess.run(
        f'conda env create -p {env_prefix} -f {env} --force',
        check=True,
        shell=True
    )
    # Install SimNIBS
    wheels = glob.glob(f'../dist/simnibs-{version}*.whl')
    if len(wheels) == 0:
        raise FileNotFoundError(f'Did not find any wheels for version {version}')
    if sys.platform == 'win32':
        env_pip = os.path.join(env_prefix, 'Scripts', 'pip.exe')
    else:
        env_pip = os.path.join(env_prefix, 'bin', 'pip') 
    subprocess.run(
        f'{env_pip} install simnibs --no-cache-dir --no-index --upgrade --find-links=../dist',
        check=True,
        shell=True
    )
    # Pack
    # I use .tar because MacOS erases the execute permission in .zip
    conda_pack.pack(
        prefix=env_prefix,
        dest_prefix='simnibs_env',
        output=os.path.join(pack_dir, 'simnibs_env.tar'),
        compress_level=0,
        force=True
    )
    shutil.unpack_archive(
        os.path.join(pack_dir, 'simnibs_env.tar'),
        os.path.join(pack_dir, 'simnibs_env'),
    )
    os.remove(os.path.join(pack_dir, 'simnibs_env.tar'))
    # Remove temporary env
    subprocess.run(
        f'conda env remove -y -p {env_prefix}',
        check=True,
        shell=True
    )

    # Copy documentation
    shutil.copytree('../docs/build/html', os.path.join(pack_dir, 'documentation'))

    # Copy postinstall script
    shutil.copy('../simnibs/cli/postinstall_simnibs.py', pack_dir) 

    # Create OS-specific installer
    if sys.platform == 'win32':
        #Use the installer.nsi template to create an NSIS installer
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
        subprocess.run([
            'pkgbuild',
            '--root', pack_dir,
            '--identifier', f'org.SimNIBS.{version}',
            '--version', version,
            '--scripts', 'macOS_scripts/', 
            '--install-location', 
            '/Applications/SimNIBS-'+ '.'.join(version.split('.')[:2]),
            'simnibs_installer_macos.pkg'
            ],
            check=True,
        )
    elif sys.platform=='linux':
        fn_script = os.path.join(pack_dir, 'install')
        with open('install', 'r') as f:
            install_script = Template(f.read()).render(
                version='.'.join(version.split('.')[:2]),
                full_version=version
            )
        with open(fn_script, 'w') as f:
            f.write(install_script)
        os.chmod(
            fn_script,
            os.stat(fn_script).st_mode |
            stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
        )
        # zip the whole thing
        shutil.make_archive(f'simnibs_installer_linux', 'zip', pack_dir)

if __name__ == '__main__':
    build()
