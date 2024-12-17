import stat
import glob
import sys
import argparse
import os
import subprocess
import shutil
import tempfile
import argparse
import re


from jinja2 import Template
import conda_pack
from setuptools_scm import get_version


def build(simnibs_dist_dir, developer_id=None):
    simnibs_root_dir = os.path.normpath(os.path.join(
        os.path.abspath(os.path.dirname(__file__)),
        '..'
    ))
    version = get_version(git_describe_command="git describe --tags --abbrev=0")

    pack_dir = os.path.abspath('simnibs_installer')
    env_prefix = os.path.join(pack_dir, 'simnibs_env_tmp')
    simnibs_dist_dir = os.path.abspath(simnibs_dist_dir)

    wheels = glob.glob(
        os.path.join(simnibs_dist_dir, f'simnibs-{version}*.whl')
    )
    if len(wheels) == 0:
        raise FileNotFoundError(f'Did not find any wheels for version {version}')

    # Create temporary environment
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
    env = os.path.join(
        simnibs_root_dir, f'environment_{os_name}.yml'
    )
    subprocess.run(
        f'conda env create -p {env_prefix} -f {env} --yes',
        check=True,
        shell=True
    )
    if sys.platform == 'win32':
        env_pip = os.path.join(env_prefix, 'Scripts', 'pip.exe')
    else:
        env_pip = os.path.join(env_prefix, 'bin', 'pip')

    # Install SimNIBS
    subprocess.run(
        f'{env_pip} install simnibs=={version} --no-deps --no-index --find-links={simnibs_dist_dir}',
        check=True,
        shell=True
    )

    # Pack
    # I use .tar because MacOS erases the execute permission in .zip
    conda_pack.pack(
        prefix=env_prefix,
        output=os.path.join(pack_dir, 'simnibs_env.tar'),
        compress_level=0,
        force=True,
        ignore_missing_files=True
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
    shutil.copytree(
        os.path.join(simnibs_root_dir, 'docs', 'build', 'html'),
        os.path.join(pack_dir, 'documentation')
    )

    # Copy the fix_entrypoints script and the postinstall script
    shutil.copy(
        os.path.join(simnibs_root_dir, 'packing', 'fix_entrypoints.py'),
        os.path.join(pack_dir, 'simnibs_env')
    )

    # Create OS-specific installer
    if sys.platform == 'win32':
        # Move the sitecustomize.py file to the site-packages directory
        # This should allow for using the python interpreter without activating the environment
        shutil.copy(
            os.path.join(simnibs_root_dir, 'simnibs', '_internal_resources', 'sitecustomize.py'),
            os.path.join(pack_dir, 'simnibs_env', 'Lib', 'site-packages')
        )
        #Use the installer.nsi template to create an NSIS installer
        shutil.copy(
            os.path.join(simnibs_root_dir, 'simnibs', '_internal_resources',
                         'icons', 'simnibs', 'gui_icon.ico'),
            os.path.join(pack_dir, 'gui_icon.ico')
        )
        fn_script = os.path.join(pack_dir, 'installer.nsi')
        with open(os.path.join(simnibs_root_dir, 'packing', 'installer.nsi'), 'r') as f:
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
        shutil.move(
            os.path.join(pack_dir, 'simnibs_installer_windows.exe'),
            'simnibs_installer_windows.exe'
        )
    if sys.platform == 'darwin':
        with tempfile.TemporaryDirectory() as tmpdir:
            for fn in glob.glob(os.path.join(simnibs_root_dir, 'packing', 'macOS_installer', '*')):
                fn_out = os.path.join(tmpdir, os.path.basename(fn))
                with open(fn, 'r') as f:
                    template = Template(f.read()).render(
                        version='.'.join(version.split('.')[:2]),
                        full_version=version
                    )
                with open(fn_out, 'w') as f:
                    f.write(template)
                os.chmod(fn_out, os.stat(fn).st_mode)

            # Workaroud for Notarization
            # Instead of signing all binaries, I zip the enironment with a password
            # The postinstall script will unzip it in the user's computer
            orig_folder = os.path.abspath(os.curdir)
            os.chdir(pack_dir)
            subprocess.run([
                'zip', '-q', '-P', 'password', '-r',
                'simnibs_env.zip',
                'simnibs_env'
            ])
            os.chdir(orig_folder)
            shutil.rmtree(os.path.join(pack_dir, 'simnibs_env'))

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
            print('Running productbuild')
            if developer_id is not None:
                sign = ['--sign', developer_id]
            else:
                sign = []
            subprocess.run([
                'productbuild',
                '--distribution', os.path.join(tmpdir, 'Distribution'),
                '--package-path', tmpdir,
                '--resources', tmpdir,
                'simnibs_installer_macos.pkg'
                ] + sign,
                check=True
            )
    elif sys.platform=='linux':
        # Write the install script
        fn_script = os.path.join(pack_dir, 'install')
        with open(os.path.join(simnibs_root_dir, 'packing', 'install'), 'r') as f:
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
    parser = argparse.ArgumentParser(
        prog="simnibs-pack",
        description="Prepare a SimNIBS intaller"
    )
    parser.add_argument("dist_dir", help="Directory with the SimNIBS wheels to be packed")
    parser.add_argument("--developer-id", default=None, help="Developer ID for signing in MacOS, DOES NOT SUPPORT NOTARIZATION (optional)")
    args = parser.parse_args(sys.argv[1:])
    build(args.dist_dir, args.developer_id)
