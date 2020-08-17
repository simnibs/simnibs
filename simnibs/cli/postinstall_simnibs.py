'''
    Optional post-install procedures in SimNIBS
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2020 Guilherme B Saturnino

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

'''
import sys
import os
import shutil
import argparse
import glob
import stat
import re
import subprocess
import tempfile
import time
import functools
import zipfile

import requests

from simnibs import SIMNIBSDIR
from simnibs import __version__
from simnibs import file_finder
try:
    from PyQt5 import QtCore, QtWidgets, QtGui
    GUI = True
except ImportError:
    GUI = False

if sys.platform == 'win32':
    import winreg

MINOR_VERSION = '.'.join(__version__.split('.')[:2])

def create_scripts(dest_dir):
    ''' Create scripts to call SimNIBS
    We need to write sh/cmd scripts due to the way python handles DLLs
    Additionaly, creates a 'simnibs' command in the matlab folder
    '''
    # On windows, copy the sitecustomize script
    # in order to be able to use the python interpreter without activating the environment
    if sys.platform == 'win32':
        simnibs_sitecustomize = os.path.join(SIMNIBSDIR, '_internal_resources', 'sitecustomize.py')
        env_sitecustomize = os.path.join(os.path.dirname(sys.executable), 'Lib', 'site-packages', 'sitecustomize.py')
        write_sitecustomize =True
        with open(simnibs_sitecustomize, 'r') as f:
            simnibs_sitecustomize_contents = f.read()
        # Check if there is a sitecustomize file alread present and if it is identical to the SimNIBS one
        if os.path.isfile(env_sitecustomize):
            with open(env_sitecustomize, 'r') as f:
                env_sitecustomize_contents = f.read()
            # If it's alteady there, will not append the PATH
            write_sitecustomize = not(simnibs_sitecustomize_contents in env_sitecustomize_contents)
        if write_sitecustomize:
            with open(env_sitecustomize, 'a') as f:
                f.write('\n')
                f.write(simnibs_sitecustomize_contents)
    scripts = glob.glob(os.path.join(SIMNIBSDIR, 'cli', '[!_]*.py'))
    if not os.path.isdir(dest_dir):
        os.makedirs(dest_dir)
    # Create scripts
    for s in scripts:
        basename = os.path.splitext(os.path.basename(s))[0]
        if basename == 'run_simnibs':
            basename = 'simnibs'
        if basename == 'simnibs_gui':
            gui = True
        else:
            gui = False
        # Normal things
        bash_name = os.path.join(dest_dir, basename)

        # Special treatment to meshfix and gmsh
        if basename in ['meshfix', 'gmsh']:
            if sys.platform == 'win32':
                with open(bash_name + '.cmd', 'w') as f:
                    f.write("@echo off\n")
                    f.write(f'"{file_finder.path2bin(basename)}" %*')
            else:
                if os.path.lexists(bash_name):
                    os.remove(bash_name)
                os.symlink(file_finder.path2bin(basename), bash_name)
        # Other stuff
        else:
            if sys.platform == 'win32':
                _write_windows_cmd(s, bash_name, gui)
            else:
                _write_unix_sh(s, bash_name)
            
    # simnibs_python interpreter
    if sys.platform == 'win32':
        _write_windows_cmd(None, os.path.join(dest_dir, 'simnibs_python'))

    else:
        if os.path.lexists(os.path.join(dest_dir, 'simnibs_python')):
            os.remove(os.path.join(dest_dir, 'simnibs_python'))
        os.symlink(
            sys.executable,
            os.path.join(dest_dir, 'simnibs_python'))


def _write_unix_sh(python_cli, bash_cli, commands='"$@"'):
    ''' Writes a bash script to evoke a python program '''
    print(f'Writing {bash_cli}')
    with open(bash_cli, 'w') as f:
        f.write('#! /bin/bash -e\n')
        f.write(f'"{sys.executable}" -E -u "{python_cli}" {commands}')
    os.chmod(bash_cli,
             os.stat(bash_cli).st_mode |
             stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

def _write_windows_cmd(python_cli, bash_cli, gui=False, commands='%*'):
    bash_cli = bash_cli + '.cmd'
    executables_dir = os.path.dirname(sys.executable)
    if gui:
        python_interpreter = f'start "Loading SimNIBS" "{os.path.join(executables_dir, "pythonw.exe")}"'
    else:
        python_interpreter = f'"{os.path.join(executables_dir, "python.exe")}"'
    with open(bash_cli, 'w') as f:
        f.write("@echo off\n")
        if python_cli is None:
            f.write(f'"{python_interpreter}" %*')
        else:
            f.write(f'{python_interpreter} -E -u "{python_cli}"  {commands}')


def _get_activate_bin():
    activate_bin = os.path.abspath(os.path.join(
        os.path.dirname(sys.executable),
        '..', '..', 'Scripts', 'activate'))
    if os.path.isfile(activate_bin):
        return activate_bin
    activate_bin = os.path.abspath(os.path.join(
        os.path.dirname(sys.executable), 'Scripts', 'activate.bat'
    ))
    if os.path.isfile(activate_bin):
        return activate_bin
    else:
        raise OSError("Can't run postinstall script (not a conda environment?)")


def _get_conda_env():
    try:
        return os.environ['CONDA_DEFAULT_ENV']
    except KeyError:
        return ''

def setup_gmsh_options(force=False, silent=False):
    ''' Copies the gmsh_options file to the appropriate place '''
    if sys.platform in ['linux', 'darwin']:
        target = os.path.expanduser('~/.gmsh-options')
    else:
        target = os.path.join(os.getenv('APPDATA'), 'gmsh-options')
    copy = True
    if os.path.isfile(target):
        if force:
            copy = True
        else:
            copy = _get_input(
                'Found a gmsh configuration file, do you whish to overwite it?',
                silent)

    if copy:
        gmsh_options_path = os.path.join(
            SIMNIBSDIR, '_internal_resources', 'gmsh-options_simnibsdefault')
        _copy_and_log(gmsh_options_path, target)


def _get_input(message, silent):
    '''Simple function to get user input via command line or GUI '''
    if silent:
        answer = input(
            f'{message} [Y/n]')
        if answer in ['n', 'N', 'no', 'No']:
            return False
        elif answer in ['', 'y', 'Y', 'yes', 'Yes']:
            return True
        else:
            raise ValueError(f'Unrecognized answer: {answer}')
    else:
        answer = QtWidgets.QMessageBox.question(
            None,'SimNIBS Postinstall', message,
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
            QtWidgets.QMessageBox.Yes)
        return answer == QtWidgets.QMessageBox.Yes



def _copy_and_log(src, dst):
    print(f'Copying {src} -> {dst}')
    shutil.copy(src, dst)


def _get_bashrc():
    ''' Get the bash startup file '''
    if sys.platform == 'linux':
        bashrc = os.path.expanduser('~/.bashrc')
        backup_file = os.path.expanduser('~/.bashrc_simnibs_bk')
    else:
        bashrc = os.path.expanduser('~/.bash_profile')
        backup_file = os.path.expanduser('~/.bash_profile_simnibs_bk')
    return bashrc, backup_file

def _get_win_simnibs_env_vars():
    ''' 'Creates a disctionary with environment names and values '''
    simnibs_env_vars = {}
    with winreg.OpenKey(winreg.HKEY_CURRENT_USER, 'Environment') as reg:
        _, num_values, _ = winreg.QueryInfoKey(reg)
        for i in range(num_values):
            var_name, value, _ = winreg.EnumValue(reg, i)
            if 'SIMNIBS' in var_name:
                simnibs_env_vars[var_name] = value
    return simnibs_env_vars

def path_setup(scripts_dir, force=False, silent=False, shell_type='bash'):
    ''' Modifies the bash startup path and postpends SimNIBS to the PATH '''
    scripts_dir = os.path.abspath(scripts_dir)
    if sys.platform in ['linux', 'darwin']:
        if shell_type == 'bash':
            bashrc, _ = _get_bashrc()
        elif shell_type =='zsh':
            bashrc = os.path.expanduser('~/.zprofile')
        else:
            raise OSError('Invalid shell type')
        if os.path.exists(bashrc):
            has_simnibs = (
                re.search('simnibs', open(bashrc, 'r').read(), re.IGNORECASE)
                is not None)
        else:
            has_simnibs = False

    if sys.platform == 'win32':
        has_simnibs = False

    if has_simnibs:
        if force:
             overwrite=True
        else:
            overwrite = _get_input(
                'Found another SimNIBS install, overwite it from the PATH?',
                silent)

        if not overwrite:
            print(f'Not Adding the current SimNIBS install to the {shell_type} PATH')
            return False

        path_cleanup(scripts_dir, shell_type=shell_type)

    print(f'Postpending {scripts_dir} to the {shell_type} PATH')
    if sys.platform in ['linux', 'darwin']:
        with open(bashrc, 'a') as f:
            f.write('\n')
            f.write('## Added by SimNIBS\n')
            f.write(f'SIMNIBS_BIN="{scripts_dir}"\n')
            f.write('export PATH=${PATH}:${SIMNIBS_BIN}')

    else:
        from simnibs.utils import _system_path
        _system_path.add_to_system_path(scripts_dir, allusers=False)
        _system_path.broadcast_environment_settings_change()

    return True

def path_cleanup(scripts_dir, shell_type='bash'):
    ''' Removes SIMNIBS from PATH '''
    if sys.platform in ['linux', 'darwin']:
        if shell_type == 'bash':
            bashrc, backup_file = _get_bashrc()
        if shell_type == 'zsh':
            bashrc = os.path.expanduser('~/.zprofile')
            backup_file = os.path.expanduser('~/.zprofile_simnibs_bk')
        if not os.path.isfile(bashrc):
            print('Could not find bashrc file')
            return

        print(f'Removing SimNIBS install from {shell_type} PATH')
        print(f'Backing up the bashrc file at {backup_file}')
        _copy_and_log(bashrc, backup_file)
        with open(backup_file, 'r') as fin:
            with open(bashrc, 'w') as fout:
                for line in fin:
                    if not re.search('simnibs', line, re.IGNORECASE):
                        fout.write(line)
    else:
        simnibs_env_vars = _get_win_simnibs_env_vars()
        for key, value in simnibs_env_vars.items():
            # Remove environments variables with SimNIBS in their names.
            # These are leftovers from previous (3.0, 3.1) installs
            with winreg.OpenKey(winreg.HKEY_CURRENT_USER, 'Environment', access=winreg.KEY_WRITE) as reg:
                winreg.DeleteValue(reg, key)
                
        from simnibs.utils import _system_path
        _system_path.remove_from_system_path(scripts_dir, allusers=False)
        _system_path.broadcast_environment_settings_change()


def matlab_prepare():
    with open(os.path.join(SIMNIBSDIR, 'matlab', 'SIMNIBSDIR.m'), 'w') as f:
        f.write("function path=SIMNIBSDIR\n")
        f.write("% Function writen by SimNIBS postinstaller\n")
        f.write(f"path='{os.path.join(SIMNIBSDIR)}';\n")
        f.write("end\n")

    with open(os.path.join(SIMNIBSDIR, 'matlab', 'SIMNIBSPYTHON.m'), 'w') as f:
        python_call = f'"{sys.executable}" -E -u '
        f.write("function python_call=SIMNIBSPYTHON\n")
        f.write("% Function writen by SimNIBS postinstaller\n")
        f.write(f"python_call='{python_call}';\n")
        f.write("end\n")

def links_setup(install_dir):
    if sys.platform == 'win32':
        _create_shortcut(
            os.path.join(install_dir, 'examples'),
            os.path.join(SIMNIBSDIR, 'examples')
        )
        _create_shortcut(
            os.path.join(install_dir, 'simnibs'),
            SIMNIBSDIR
        )
        _create_shortcut(
            os.path.join(install_dir, 'resources'),
            os.path.join(SIMNIBSDIR, 'resources')
        )
        _create_shortcut(
            os.path.join(install_dir, 'matlab_tools'),
            os.path.join(SIMNIBSDIR, 'matlab_tools')
        )
    else:
        def _new_symlink(link_name, target):
            if os.path.islink(link_name):
                os.remove(link_name)
            os.symlink(target, link_name)

        _new_symlink(
            os.path.join(install_dir, 'examples'),
            os.path.join(SIMNIBSDIR, 'examples')
        )
        _new_symlink(
            os.path.join(install_dir, 'simnibs'),
            SIMNIBSDIR
        )
        _new_symlink(
            os.path.join(install_dir, 'resources'),
            os.path.join(SIMNIBSDIR, 'resources')
        )
        _new_symlink(
            os.path.join(install_dir, 'matlab_tools'),
            os.path.join(SIMNIBSDIR, 'matlab_tools')
        )

def setup_shortcut_icons(scripts_dir, force=False, silent=False):
    ''' Creates shortcut icons for the gui_scripts '''
    if sys.platform == 'darwin':
        _create_apps(os.path.abspath(os.path.join(scripts_dir, '..')))
        return
    elif sys.platform == 'win32':
        shortcut_folder = os.path.join(
            os.environ['APPDATA'],
            "Microsoft", "Windows", "Start Menu",
            "Programs", f"SimNIBS {MINOR_VERSION}"
        )
        gmsh_icon = None
        simnibs_icon = os.path.join(SIMNIBSDIR, '_internal_resources', 'icons', 'simnibs',  'gui_icon.ico')

    elif sys.platform == 'linux':
        shortcut_folder = os.path.expanduser(f'~/.local/share/applications/SimNIBS-{MINOR_VERSION}')
        gmsh_icon = os.path.join(SIMNIBSDIR, '_internal_resources', 'icons', 'gmsh', 'logo.png')
        simnibs_icon = os.path.join(SIMNIBSDIR, '_internal_resources', 'icons', 'simnibs', 'gui_icon.png')


    if os.path.isdir(shortcut_folder):
        if force:
            overwrite = True
        else:
            overwrite = _get_input(
                'Found other SimNIBS menu icons, overwrite them?',
                silent)
        if not overwrite:
            print('Not adding shortucts to the current SimNIBS install')
            return
        
        else:
            shortcut_icons_clenup()

    os.makedirs(shortcut_folder, exist_ok=True)

    _create_shortcut(
        os.path.join(shortcut_folder, 'Gmsh'),
        file_finder.path2bin('gmsh'),
        'Gmsh is a free 3D finite element mesh generator with a built-in CAD engine and'
        ' post-processor',
        mime_type='model/x.stl-binary',
        icon=gmsh_icon
    )
    fn_gui_script = os.path.join(scripts_dir, 'simnibs_gui')
    if sys.platform == 'win32':
        fn_gui_script += '.cmd'
    _create_shortcut(
        os.path.join(shortcut_folder, 'SimNIBS GUI'),
        fn_gui_script,
        'SimNIBS is software for simulating electric fields caused by NIBS',
        icon=simnibs_icon
    )
    if sys.platform == 'win32':
        _create_shortcut(
            os.path.join(shortcut_folder, 'SimNIBS Directory'),
            os.path.abspath(os.path.join(scripts_dir, '..')),
        )
        _create_shortcut(
            os.path.join(shortcut_folder, 'SimNIBS Documentation'),
            os.path.abspath(os.path.join(scripts_dir, '..', 'documentation', 'index.html')),
        )
        _create_shortcut(
             os.path.join(shortcut_folder, 'SimNIBS Prompt'),
             r'%windir%\System32\cmd.exe',
             arguments=f'/K ""{_get_activate_bin()}"" {_get_conda_env()}')
    if sys.platform == 'linux':
        try:
            subprocess.run(
                ['update-desktop-database',
                 os.path.expanduser('~/.local/share/applications')])
        except:
            print('Could not update desktop database')


def _create_shortcut(shortcut_name, target_path, comment='', icon=None, mime_type=None, arguments=None):
    if sys.platform == 'win32':
        with tempfile.NamedTemporaryFile('w', delete=False, suffix='.ps1') as f:
            f.write(f'$objShell = New-Object -ComObject ("WScript.Shell")\n')
            f.write(f'$objShortCut = $objShell.CreateShortcut("{shortcut_name}.lnk")\n')
            f.write(f'$objShortCut.TargetPath="{target_path}"\n')
            f.write(f'$objShortCut.WorkingDirectory="%HOMEPATH%"\n')
            if icon:
                f.write(f'$objShortCut.IconLocation="{icon}"\n')
            if arguments:
                f.write(f'$objShortCut.Arguments="{arguments}"\n')
            f.write(f'$objShortCut.Save()')
            temp_fn = f.name
        subprocess.run(f'powershell.exe -noprofile -executionpolicy bypass -file "{temp_fn}"', shell=True).check_returncode()
        os.remove(temp_fn)
    elif sys.platform == 'linux':
        with open(shortcut_name + '.desktop', 'w') as f:
            f.write('[Desktop Entry]\n')
            f.write(f'Name={os.path.basename(shortcut_name)}\n')
            f.write('Comment={comment}\n')
            f.write(f'Exec={target_path} %f\n')
            if icon:
                f.write(f'Icon={icon}\n')
            if mime_type:
                f.write(f'MimeType={mime_type}\n')
            f.write('Terminal=false\n')
            f.write('Encoding=UTF-8\n')
            f.write('Type=Application\n')


def _create_apps(install_dir):
    """ Creates an apps for MacOS
    """

    plist = dict(
        CFBundleDisplayName="SimNIBS GUI",
        CFBundleName="SimNIBS GUI",
        CFBundleIdentifier="org.simnibs",
        CFBundleShortVersionString=__version__,
        CFBundleGetInfoString=f'SimNIBS GUI {__version__}',
        CFBundleIconFile="gui_icon.icns",
        CFBundleExecutable="simnibs_gui",
        CFBundleInfoDictionaryVersion='6.0'
    ) 
    _create_app(
        os.path.join(install_dir, 'SimNIBS GUI.app'),
        os.path.join(SIMNIBSDIR, 'cli', 'simnibs_gui.py'),
        os.path.join(SIMNIBSDIR, '_internal_resources', 'icons', 'simnibs', 'gui_icon.icns'),
        plist)

    # Gmsh app setup
    target_dir = os.path.join(install_dir, 'Gmsh.app')
    if os.path.isdir(target_dir):
        shutil.rmtree(target_dir)
    print('Installing Gmsh')
    contents_dir = os.path.join(target_dir, 'Contents')
    resouces_dir = os.path.join(contents_dir, 'Resources')
    macos_dir = os.path.join(contents_dir, 'MacOS')
    os.mkdir(target_dir)
    os.mkdir(contents_dir)
    os.mkdir(resouces_dir)
    os.mkdir(macos_dir)

    _copy_and_log(
        os.path.join(SIMNIBSDIR, '_internal_resources', 'icons', 'gmsh', 'Info.plist'),
        os.path.join(contents_dir))

    for icns in glob.glob(os.path.join(SIMNIBSDIR, '_internal_resources', 'icons', 'gmsh', '*.icns')):
        _copy_and_log(icns, resouces_dir)

    _copy_and_log(
        file_finder.path2bin('gmsh'),
        os.path.join(macos_dir))

def _create_app(app_name, executable, icon, plist):
    '''
    https://stackoverflow.com/questions/1596945/building-osx-app-bundle
    https://developer.apple.com/library/archive/documentation/CoreFoundation/Conceptual/CFBundles/BundleTypes/BundleTypes.html#//apple_ref/doc/uid/10000123i-CH101-SW16
    Only works if app_name finishes in ".app"
    '''
    import plistlib 
    # SimNIBS app setup
    if os.path.isdir(app_name):
        shutil.rmtree(app_name)
    os.makedirs(app_name)
    contents_dir = os.path.join(app_name, 'Contents')
    resouces_dir = os.path.join(contents_dir, 'Resources')
    macos_dir = os.path.join(contents_dir, 'MacOS')
    if not os.path.isdir(contents_dir):
        os.mkdir(contents_dir)
    if not os.path.isdir(resouces_dir):
        os.mkdir(resouces_dir)
    if not os.path.isdir(macos_dir):
        os.mkdir(macos_dir)
    _copy_and_log(
        icon,
        os.path.join(resouces_dir))
    # Write a simnibs_gui and a gmsh executable
    _write_unix_sh(
        executable,
        os.path.join(macos_dir, 'simnibs_gui'))
    with open(os.path.join(contents_dir, 'Info.plist'), 'wb') as fp:
        plistlib.dump(plist, fp)

def shortcut_icons_clenup():
    if sys.platform == 'win32':
        shortcut_folder=os.path.join(
            os.environ['APPDATA'],
            "Microsoft", "Windows", "Start Menu",
            "Programs", f"SimNIBS {MINOR_VERSION}"
        )
    elif sys.platform == 'linux':
        shortcut_folder = os.path.expanduser('~/.local/share/applications/SimNIBS')
    else:
        return
    if os.path.isdir(shortcut_folder):
        shutil.rmtree(shortcut_folder)
        while os.path.exists(shortcut_folder):
            pass

def setup_file_association(force=False, silent=False):
    # Linux and OSX file associations are done together with desktop items
    if sys.platform != 'win32':
        return
    gmsh_bin = file_finder.path2bin('gmsh')
    extensions = ['.msh', '.geo', '.stl']
    associate = dict.fromkeys(extensions)
    for ext in extensions:
        if _is_associated(ext):
            if force:
                associate[ext] = True
            else:
                associate[ext] = _get_input(
                    f'Found other association for "{ext}" files, overwrite it?',
                    silent)
        else:
            associate[ext] = True
    # If all rejected, return
    if not any(associate.values()):
        return

    if sys.platform == 'win32':
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, r'Software\Classes', access=winreg.KEY_WRITE) as reg:
            winreg.CreateKey(reg, rf'SimNIBS.Gmsh.v{MINOR_VERSION}\shell\open\command')
            winreg.SetValue(reg, rf'SimNIBS.Gmsh.v{MINOR_VERSION}\shell\open\command', winreg.REG_SZ, f'"{gmsh_bin}" "%1"')
            for ext in extensions:
                try:
                    value = winreg.QueryValue(reg, ext)
                except FileNotFoundError:
                    register = True
                else:
                    if value:
                        register = _get_input(
                            f'Found other association for "{ext}" files, overwrite it?',
                            silent
                        )
                    else:
                        register = True
                if register:
                    winreg.CreateKey(reg, ext)
                    winreg.SetValue(reg, ext, winreg.REG_SZ, fr'SimNIBS.Gmsh.v{MINOR_VERSION}')

def _is_associated(ext):
    if sys.platform == 'win32':
        # Also needs to be run as shell
        ret = subprocess.run('assoc '+ ext, shell=True, capture_output=True, text=True,
                             errors='replace')
        if ret.returncode == 0:
            assoc = re.findall(f'{ext}=(.*)\n', ret.stdout)
            if len(assoc) == 0:
                return False
            return assoc[0]
        else:
            return False

def file_associations_cleanup():
    extensions = ['.msh', '.geo', '.stl']
    # Linux file associations are done together with desktop items
    # MacOS file associations are set using the .app files
    if sys.platform in ['linux', 'darwin']:
        return
    
    if sys.platform == 'win32':
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER, r'Software\Classes', access=winreg.KEY_WRITE) as reg:
            # Remove SimNIBS Gmsh call from the registry
            try:
                winreg.QueryValue(reg, rf'SimNIBS.Gmsh.v{MINOR_VERSION}\shell\open\command')
            except FileNotFoundError:
                pass
            else:
                # Delete recursivelly
                paths = rf'SimNIBS.Gmsh.v{MINOR_VERSION}\shell\open\command'.split('\\')
                for i in reversed(range(len(paths))):
                    try:
                        winreg.DeleteKey(reg, '\\'.join(paths[:i+1]))
                    except OSError:
                        break
            # Remove the extensions from the registry
            for ext in extensions:
                try:
                    entry = winreg.QueryValue(reg, ext)
                except FileNotFoundError:
                    pass
                else:
                    if entry == fr'SimNIBS.Gmsh.v{MINOR_VERSION}':
                        winreg.SetValue(reg, ext, winreg.REG_SZ, '')

def uninstaller_setup(install_dir, force, silent):
    uninstaller = os.path.join(install_dir, 'uninstall_simnibs')
    simnibs_env_dir = os.path.join(install_dir, 'simnibs_env')
    if sys.platform == 'win32':
        _write_windows_cmd(
            os.path.join(SIMNIBSDIR, 'cli', 'postinstall_simnibs.py'),
            uninstaller, commands=f'-u %* -d "{install_dir}"',
            gui=False)
    else:
        _write_unix_sh(
            os.path.join(SIMNIBSDIR, 'cli', 'postinstall_simnibs.py'),
            uninstaller, commands=f'-u "$@" -d "{install_dir}"')
        with open(uninstaller, 'a') as f:
            f.write(
                f'&& rm -rf {simnibs_env_dir} '
                f'; rm "{uninstaller}" '
                f'; rm -rf "{install_dir}"')

def activator_setup(install_dir):
    activator = os.path.join(install_dir, 'activate_simnibs')
    if sys.platform == 'win32':
        _write_windows_cmd(
            os.path.join(SIMNIBSDIR, 'cli', 'postinstall_simnibs.py'),
            activator, gui=True, commands=f'-d "{install_dir}"')

        _create_shortcut(
            os.path.join(install_dir, 'Activate SimNIBS'),
            activator,
            icon=os.path.join(SIMNIBSDIR, '_internal_resources', 'icons', 'simnibs', 'gui_icon.ico')
        )
    else:
        _write_unix_sh(
            os.path.join(SIMNIBSDIR, 'cli', 'postinstall_simnibs.py'),
            activator, commands=f'-d "{install_dir}"')

# from https://stackoverflow.com/questions/16694907/download-large-file-in-python-with-requests
def download_file(url, local_filename, timeout=None):
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raw.read = functools.partial(r.raw.read, decode_content=True)
        with open(local_filename, 'wb') as f:
            shutil.copyfileobj(r.raw, f)

def download_extra_coils(timeout=None):
    version = 'master'
    url = f'https://github.com/simnibs/simnibs-coils/archive/{version}.zip'
    with tempfile.NamedTemporaryFile('wb', delete=False) as tmpf:
        tmpname = tmpf.name
    download_file(url, tmpf.name, timeout)
    with zipfile.ZipFile(tmpname) as z:
        z.extractall(file_finder.coil_models)
    os.remove(tmpname)
    src = os.path.join(file_finder.coil_models, f'simnibs-coils-{version}')
    for f in glob.glob(os.path.join(src, '*')):
        d = os.path.join(file_finder.coil_models, os.path.basename(f))
        if os.path.isdir(d):
            shutil.rmtree(d)
        if os.path.isfile(d):
            os.remove(d)
        shutil.move(f, d)
    shutil.rmtree(
        os.path.join(SIMNIBSDIR, 'ccd-files', f'simnibs-coils-{version}'))

def run_tests(args):
    ''' run tests on pytest '''
    import pytest
    exitcode = pytest.main(args)
    return exitcode

if GUI:
    class PostInstallGUI(QtWidgets.QDialog):
        def __init__(self,
                     install_dir,
                     copy_gmsh_options=True,
                     add_to_path=True,
                     extra_coils=True,
                     add_shortcut_icons=True,
                     associate_files=True):
            super().__init__()
            self.install_dir = os.path.abspath(os.path.expanduser(install_dir))
            self.copy_gmsh_options = copy_gmsh_options
            self.add_to_path = add_to_path
            self.extra_coils = extra_coils
            self.add_shortcut_icons = add_shortcut_icons
            self.associate_files = associate_files

            install_box = QtWidgets.QGroupBox('Installing SimNIBS to:')
            layout = QtWidgets.QHBoxLayout()
            self.install_line_edit = QtWidgets.QLineEdit(self.install_dir)
            self.install_line_edit.setEnabled(False)
            layout.addWidget(self.install_line_edit)
            install_box.setLayout(layout)

            ## Set check boxes
            layout = QtWidgets.QVBoxLayout()
            options_box = QtWidgets.QGroupBox()
            # gmsh options
            gmsh_options_cb = QtWidgets.QCheckBox('Copy gmsh_options file')
            if self.copy_gmsh_options: gmsh_options_cb.toggle()
            gmsh_options_cb.toggled.connect(self.set_copy_gmsh_options)
            layout.addWidget(gmsh_options_cb, 0, QtCore.Qt.Alignment(1))

            # path options
            path_cb = QtWidgets.QCheckBox('Add SimNIBS to the system PATH')
            if self.add_to_path: path_cb.toggle()
            path_cb.toggled.connect(self.set_add_to_path)
            layout.addWidget(path_cb, 0, QtCore.Qt.Alignment(1))

            # extra coils options
            coils_cb = QtWidgets.QCheckBox('Download additional coil files')
            if self.extra_coils: coils_cb.toggle()
            coils_cb.toggled.connect(self.set_extra_coils)
            layout.addWidget(coils_cb, 0, QtCore.Qt.Alignment(1))


            # shortcut options
            shortcut_cb = QtWidgets.QCheckBox('Add SimNIBS shortcut icons')
            if self.add_shortcut_icons: shortcut_cb.toggle()
            shortcut_cb.toggled.connect(self.set_add_shortcut_icons)
            layout.addWidget(shortcut_cb, 0, QtCore.Qt.Alignment(1))

            # file association options
            # Only separate from the shortcut stuff in Windows
            if sys.platform == 'win32':
                associate_cb = QtWidgets.QCheckBox('Associate .msh, .stl and .geo files')
                if self.associate_files: associate_cb.toggle()
                associate_cb.toggled.connect(self.set_associate_files)
                layout.addWidget(associate_cb, 0, QtCore.Qt.Alignment(1))

            options_box.setLayout(layout)


            # Button Boxes
            button_box = QtWidgets.QDialogButtonBox(
                QtWidgets.QDialogButtonBox.Ok|QtWidgets.QDialogButtonBox.Cancel)
            button_box.accepted.connect(self.accept)
            button_box.rejected.connect(self.reject)

            # Join Everything
            mainLayout = QtWidgets.QVBoxLayout()
            mainLayout.addWidget(install_box)
            mainLayout.addWidget(options_box)
            mainLayout.addWidget(button_box)

            self.setLayout(mainLayout)

            self.setWindowTitle(f'SimNIBS {__version__} Post-Install Options')
            gui_icon = os.path.join(SIMNIBSDIR, '_internal_resources', 'icons', 'simnibs', 'gui_icon.ico')
            self.setWindowIcon(QtGui.QIcon(gui_icon))

        def set_copy_gmsh_options(self, new_value):
            self.copy_gmsh_options = new_value
            if not new_value:
                QtWidgets.QMessageBox.warning(
                    self, 'SimNIBS',
                    'Visualization of head models and simulation results will be affected')

        def set_add_to_path(self, new_value):
            self.add_to_path = new_value
            if not new_value:
                QtWidgets.QMessageBox.warning(
                    self, 'SimNIBS',
                    'SimNIBS functions will not be callable from the terminal')

        def set_add_shortcut_icons(self, new_value):
            self.add_shortcut_icons = new_value
            if not new_value:
                QtWidgets.QMessageBox.warning(
                    self, 'SimNIBS',
                    'SimNIBS icons will not be found in the Start Menu')

        def set_associate_files(self, new_value):
            self.associate_files = new_value
            if not new_value:
                QtWidgets.QMessageBox.warning(
                    self, 'SimNIBS',
                    'Models and simulation results will not be associated with Gmsh')

        def set_extra_coils(self, new_value):
            self.extra_coils = new_value
            if not new_value:
                QtWidgets.QMessageBox.warning(
                    self, 'SimNIBS',
                    'Less coil files will be availiable')


    def start_gui(simnibsdir,
                  setup_links,
                  copy_gmsh_options,
                  add_to_path,
                  extra_coils,
                  add_shortcut_icons,
                  associate_files):
        app = QtWidgets.QApplication(sys.argv)
        ex = PostInstallGUI(
            simnibsdir,
            copy_gmsh_options,
            add_to_path,
            extra_coils,
            add_shortcut_icons,
            associate_files
        )
        ex.show()
        app.exec_()
        if ex.result():
            install(ex.install_dir,
                    False, False,
                    ex.copy_gmsh_options,
                    ex.add_to_path,
                    ex.extra_coils,
                    ex.add_shortcut_icons,
                    ex.associate_files,
                    setup_links)
        else:
            raise Exception('Post-installation cancelled by user')

    class UnintallerGUI(QtWidgets.QDialog):
        def __init__(self):
            super().__init__()
            button_box = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok|QtWidgets.QDialogButtonBox.Cancel)
            button_box.accepted.connect(self.accept)
            button_box.rejected.connect(self.reject)
            mainLayout = QtWidgets.QVBoxLayout()
            mainLayout.addWidget(
                QtWidgets.QLabel(
                    f'SimNIBS version {__version__} will be uninstalled. '
                    'Are you sure?'))
            mainLayout.addWidget(button_box)
            self.setLayout(mainLayout)
            self.setWindowTitle('SimNIBS Uninstaller')
            gui_icon = os.path.join(SIMNIBSDIR, '_internal_resources', 'icons', 'simnibs', 'gui_icon.ico')
            self.setWindowIcon(QtGui.QIcon(gui_icon))

    def start_uninstall_gui(install_dir):
        app = QtWidgets.QApplication(sys.argv)
        ex = UnintallerGUI()
        ex.show()
        app.exec_()
        if ex.result():
            uninstall(install_dir)
        else:
            raise Exception('uninstall cancelled by user')

def install(install_dir,
            force,
            silent,
            copy_gmsh_options=True,
            add_to_path=True,
            extra_coils=True,
            add_shortcut_icons=False,
            associate_files=False,
            setup_links=False):
    install_dir = os.path.abspath(os.path.expanduser(install_dir))
    os.makedirs(install_dir, exist_ok=True)
    scripts_dir = os.path.join(install_dir, 'bin')
    matlab_prepare()
    if extra_coils:
        print('Downloading extra coils, this might take some time', flush=True)
        download_extra_coils(timeout=10*60)
    if copy_gmsh_options:
        print('Copying Gmsh Options')
        setup_gmsh_options(force, silent)
    if add_shortcut_icons:
        print('Adding Shortcut Icons', flush=True)
        setup_shortcut_icons(scripts_dir, force, silent)
    if associate_files:
        print('Associating Files', flush=True)
        setup_file_association(force, silent)
    if setup_links:
        links_setup(install_dir)
    activator_setup(install_dir)
    uninstaller_setup(install_dir, force, silent)
    test_call = [
        os.path.join(SIMNIBSDIR, 'tests', 'simulation', 'test_fem.py'),
        '-k', 'TestSolve', '-q', '-W', 'ignore'
    ]
    run_tests(test_call)
    if sys.platform == 'win32':
        pythonw = os.path.join(os.path.dirname(sys.executable), 'pythonw')
        subprocess.run([pythonw, '-m', 'pytest'] + test_call)
    shutil.rmtree(os.path.join(install_dir, '.pytest_cache'), True)
    create_scripts(scripts_dir)
    if add_to_path:
        try:
            added_to_path = path_setup(scripts_dir, force, silent)
            if sys.platform == 'darwin' and added_to_path:
                path_setup(scripts_dir, force=True, silent=True, shell_type='zsh')
        except:
            print('Could not add SimNIBS to the system PATH')

def uninstall(install_dir):
    path_cleanup(os.path.join(install_dir, 'bin'))
    if sys.platform == 'darwin':
        path_cleanup(os.path.join(install_dir, 'bin'), shell_type='zsh')
    shortcut_icons_clenup()
    file_associations_cleanup()
    shutil.rmtree(os.path.join(install_dir, 'documentation'), True)
    shutil.rmtree(os.path.join(install_dir, 'bin'), True)
    def try_remove(f):
        try:
            os.remove(f)
        except OSError:
            pass

    if sys.platform == 'win32':
        try_remove(os.path.join(install_dir, 'activate_simnibs.cmd'))
        try_remove(os.path.join(install_dir, 'examples.lnk'))
        try_remove(os.path.join(install_dir, 'simnibs.lnk'))
        try_remove(os.path.join(install_dir, 'matlab_tools.lnk'))
        try_remove(os.path.join(install_dir, 'Activate SimNIBS.lnk'))
        try_remove(os.path.join(install_dir, 'Uninstall SimNIBS.lnk'))
        conda_uninstaller = os.path.join(
             install_dir, 'miniconda3',
             'Uninstall-Miniconda3.exe')
        if os.path.isfile(conda_uninstaller):
            subprocess.run(f'"{conda_uninstaller}" /S', shell=True)

    else:
        try_remove(os.path.join(install_dir, 'activate_simnibs'))
        try_remove(os.path.join(install_dir, 'simnibs'))
        try_remove(os.path.join(install_dir, 'examples'))
        try_remove(os.path.join(install_dir, 'matlab_tools'))
        if sys.platform == 'darwin':
            shutil.rmtree(os.path.join(install_dir, 'SimNIBS GUI.app'), True)
            shutil.rmtree(os.path.join(install_dir, 'Gmsh.app'), True)

def main():
    parser = argparse.ArgumentParser(prog="postinstall_simnibs",
                                     description="Optional post-installation procedures "
                                     "for SimNIBS ")
    parser.add_argument('-d', "--target_dir", required=True,
                        help="SimNIBS install directory")
    parser.add_argument('-f', "--force", required=False, action='store_true',
                        help="Perform all install steps")
    parser.add_argument('-s', "--silent", required=False, action='store_true',
                        help="Silent mode, will install without the GUI")
    parser.add_argument('--setup-links', action='store_true',
                        help='Setups links or shortcuts (on windows) to the simnibs '
                        'and example folders')
    parser.add_argument('--no-copy-gmsh-options', dest='copy_gmsh_options',
                        action='store_false', help='Do not copy gmsh options')
    parser.add_argument('--no-add-to-path', dest='add_to_path',
                        action='store_false', help='Do not add SimNBIS to system PATH')
    parser.add_argument('--no-extra-coils', dest='extra_coils',
                        action='store_false', help='Do not download extra coil files')
    parser.add_argument('--no-shortcut-icons', dest='add_shortcut_icons',
                        action='store_false', help='Do not create start menu shortcuts')
    parser.add_argument('--no-associate-files', dest='associate_files',
                        action='store_false', help='Do not create file associations')
    parser.add_argument('-u', "--uninstall", required=False, action='store_true',
                        help="Ignores all other arguments and uninstall SimNIBS")
    parser.add_argument('--version', action='version', version=__version__)
    args = parser.parse_args(sys.argv[1:])
    install_dir = os.path.abspath(os.path.expanduser(args.target_dir))
    if args.uninstall:
        if args.silent:
            uninstall(install_dir)
        else:
            start_uninstall_gui(install_dir)
        return
    if not args.silent:
        if not GUI:
            raise ImportError(
                'Trying to run post-install script without PyQt istalled, '
                'please use the silent mode (postinstall_simnibs --help for '
                'more information')
        start_gui(
            install_dir,
            args.setup_links,
            args.copy_gmsh_options,
            args.add_to_path,
            args.extra_coils,
            args.add_shortcut_icons,
            args.associate_files
        )

    else:
        install(
            install_dir,
            args.force,
            args.silent,
            args.copy_gmsh_options,
            args.add_to_path,
            args.extra_coils,
            args.add_shortcut_icons,
            args.associate_files,
            args.setup_links
        )

if __name__ == '__main__':
    main()

