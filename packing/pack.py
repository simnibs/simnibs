''' Has to be run from the conda root environment!!
'''
import stat
import glob
import sys
import argparse
import os
import subprocess
import shutil
import tempfile

from jinja2 import Template


def build():
    import conda_pack
    version = open("../simnibs/_version.py").readlines()[-1].split()[-1].strip("\"'")
    pack_dir = os.path.abspath('simnibs_installer')
    env_prefix = os.path.join(pack_dir, 'simnibs_env_tmp')
    # Create a new environment
    if os.path.dirname(__file__):
        os.chdir(os.path.dirname(__file__))
    if os.path.isdir(pack_dir):
        shutil.rmtree(pack_dir)
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

    # Copy the fix_entrypoints script and the postinstall script
    shutil.copy('fix_entrypoints.py', os.path.join(pack_dir, 'simnibs_env'))

    # Create OS-specific installer
    if sys.platform == 'win32':
        # Move the sitecustomize.py file to the site-packages directory
        # This should allow for using the python interpreter without activating the environment
        shutil.copy('../simnibs/utils/sitecustomize.py', os.path.join(pack_dir, 'simnibs_env', 'Lib', 'site-packages'))
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
        print('Creating NSIS installer')
        subprocess.run(
            fr'"%programfiles(x86)%\NSIS\makensis.exe" {fn_script}',
            check=True,
            shell=True
        )
    if sys.platform == 'darwin':
        with tempfile.TemporaryDirectory() as tmpdir:
            for fn in glob.glob('macOS_installer/*'):
                fn_out = os.path.join(tmpdir, os.path.basename(fn))
                with open(fn, 'r') as f:
                    template = Template(f.read()).render(
                        version='.'.join(version.split('.')[:2]),
                        full_version=version
                    )
                with open(fn_out, 'w') as f:
                    f.write(template)
                os.chmod(fn_out, os.stat(fn).st_mode)
            print('Running pkgbuild')
            subprocess.run([
                'pkgbuild',
                '--root', pack_dir,
                '--identifier', f'org.SimNIBS.{version}',
                '--version', version,
                '--scripts', tmpdir, 
                '--install-location', 
                '/Applications/SimNIBS-'+ '.'.join(version.split('.')[:2]),
                os.path.join(tmpdir, 'simnibs_installer_macos.pkg')
                ],
                check=True,
            )
            print('Running productbuid')
            subprocess.run([
                'productbuild',
                '--distribution', os.path.join(tmpdir, 'Distribution'),
                '--package-path', tmpdir,
                '--resources', tmpdir,
                'simnibs_installer_macos.pkg'
                ],
                check=True
            )
    elif sys.platform=='linux':
        # Write the install script
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
        print('compressing')
        shutil.make_archive(
            'simnibs_installer_linux',
            'gztar',
            # I use root_dir and base_dir so that it decompresses into a folder called
            # simnibs_installer
            root_dir='.',
            base_dir=os.path.relpath(pack_dir)
        )

if __name__ == '__main__':
    build()
