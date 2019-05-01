'''
    Optional post-install procedures in SimNIBS
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2019  Guilherme B Saturnino

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
from simnibs import SIMNIBSDIR
from simnibs import __version__
from simnibs import file_finder
try:
    from PyQt5 import QtCore, QtWidgets, QtGui
    GUI = True
except ImportError:
    GUI = False

def copy_scripts(dest_dir):
    ''' Create scripts to call SimNIBS
    We need to write sh/cmd scripts due to the way python handles DLLs
    Additionaly, creates a 'simnibs' command in the matlab folder
    '''
    scripts = glob.glob(os.path.join(SIMNIBSDIR, 'cli', '[!_]*.py'))
    script_folder = os.path.basename(sys.executable)
    if not os.path.isdir(dest_dir):
        os.makedirs(dest_dir)
    # Create scripts
    for s in scripts:
        basename = os.path.splitext(os.path.basename(s))[0]
        if basename == 'run_simnibs':
            basename = 'simnibs'
        if basename in ['simnibs_gui', 'gmsh']:
            gui = True
        else:
            gui = False
        # Norma things
        bash_name = os.path.join(dest_dir, basename)
        if sys.platform == 'win32':
            _write_windows_cmd(s, bash_name, gui)
        else:
            _write_unix_sh(s, bash_name)
            
    # Create an additionl simnibs script in the matlab folder
    if sys.platform == 'win32':
        _write_windows_cmd(os.path.join(SIMNIBSDIR, 'cli', 'run_simnibs.py'),
                           os.path.join(SIMNIBSDIR, 'matlab', 'simnibs'))
    else:
        _write_unix_sh(os.path.join(SIMNIBSDIR, 'cli', 'run_simnibs.py'),
                       os.path.join(SIMNIBSDIR, 'matlab', 'simnibs'))

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
        f.write(f'{sys.executable} -E -u {python_cli} {commands}')
    os.chmod(bash_cli,
             os.stat(bash_cli).st_mode |
             stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

def _write_windows_cmd(python_cli, bash_cli, gui=False, commands='%*'):
    bash_cli = bash_cli + '.cmd'
    # I need to activate the environment first
    activate_bin = os.path.abspath(os.path.join(
        os.path.dirname(sys.executable),
        '..', '..', 'Scripts', 'activate'))
    if not os.path.isfile(activate_bin):
        raise OSError("Can't run postinstall script (not a conda environment)")
    if gui:
        python_interpreter = 'start pythonw'
    else:
        python_interpreter = 'python'
    with open(bash_cli, 'w') as f:
        f.write("@echo off\n")
        f.write(f'call "{activate_bin}" simnibs_env\n')
        if python_cli is None:
            f.write(f'{python_interpreter} %*')
        else:
            f.write(f'{python_interpreter} -E -u "{python_cli}"  {commands}')

def setup_gmsh_options(force=False, silent=False):
    ''' Copies the gmsh_options file to the appropriate place '''
    if sys.platform in ['linux', 'darwin']:
        target = os.path.expanduser('~/.gmsh-options')
    else:
        target = os.path.join(os.getenv('APPDATA'), 'gmsh-options')

    silent = silent and GUI
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
            SIMNIBSDIR, 'resources', 'gmsh-options_simnibsdefault')
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
    res = subprocess.run(
        'reg query HKEY_CURRENT_USER\Environment',
        capture_output=True, shell=True)
    res.check_returncode()
    simnibs_env_vars = re.findall('\s+(SIMNIBS\w+)\s+REG_SZ\s+(\S+)', res.stdout.decode())
    simnibs_env_vars = {s[0]: s[1] for s in simnibs_env_vars}
    return simnibs_env_vars

def _get_win_path():
    res = subprocess.run(
        'reg query HKEY_CURRENT_USER\Environment '
        '/f Path', capture_output=True, shell=True)
    try:
        res.check_returncode()
    except subprocess.CalledProcessError:
        return ''
    path = re.findall('REG\S+\s+(.*)\r', res.stdout.decode())
    if len(path) == 0:
        raise OSError('Could not derermine the system PATH')
    return path[0]

def path_setup(scripts_dir, force=False, silent=False):
    ''' Modifies the bash startup path and postpends SimNIBS to the PATH '''
    scripts_dir = os.path.abspath(scripts_dir)
    silent = silent and GUI
    if sys.platform in ['linux', 'darwin']:
        bashrc, backup_file = _get_bashrc()
        if os.path.exists(bashrc):
            has_simnibs = (
                re.search('simnibs', open(bashrc, 'r').read(), re.IGNORECASE)
                is not None)
        else:
            has_simnibs = False

    if sys.platform == 'win32':
        simnibs_env_vars = _get_win_simnibs_env_vars()
        has_simnibs = len(simnibs_env_vars) != 0

    if has_simnibs:
        if force:
             overwrite=True
        else:
            overwrite = _get_input(
                'Found another SimNIBS install, overwite it from the PATH?',
                silent)

        if not overwrite:
            print('Not Adding the current SimNIBS install to the PATH')
            return

        path_cleanup()

    print(f'Postpending {scripts_dir} to the PATH')
    if sys.platform in ['linux', 'darwin']:
        with open(bashrc, 'a') as f:
            f.write('\n')
            f.write('## Added by SimNIBS\n')
            f.write(f'SIMNIBS_BIN="{scripts_dir}"\n')
            f.write('export PATH=${PATH}:${SIMNIBS_BIN}')

    else:
        res = subprocess.run(
            'reg add HKEY_CURRENT_USER\Environment '
            f'/v SIMNIBS_BIN /d "{scripts_dir}" /f',
            shell=True)
        path = _get_win_path()
        path = scripts_dir + ';' + path
        subprocess.run(
            'reg add HKEY_CURRENT_USER\Environment '
            f'/v Path /t REG_EXPAND_SZ /d "{path}" /f',
            shell=True).check_returncode()
        # Signal that the variables changed
        #subprocess.run('setx simnibs_throwaway "trash"', shell=True)
        #subprocess.run(
        #    'reg delete HKEY_CURRENT_USER\Environment',
        #    f'/v simnibs_throwaway /f', shell=True)

def path_cleanup():
    ''' Removes SIMNIBS from PATH '''
    if sys.platform in ['linux', 'darwin']:
        bashrc, backup_file = _get_bashrc()

        if not os.path.isfile(bashrc):
            print('Could not find bashrc file')
            return

        print('Removing SimNIBS install from PATH')
        print(f'Backing up the bashrc file at {backup_file}')
        _copy_and_log(bashrc, backup_file)
        with open(backup_file, 'r') as fin:
            with open(bashrc, 'w') as fout:
                for i, line in enumerate(fin):
                    if not re.search('simnibs', line, re.IGNORECASE):
                        fout.write(line)
    else:
        # Requires admin privileges
        #registry_bk = os.path.abspath('registry_bk.hiv')
        #print(f'Backing up registry at {registry_bk}')
        #res = subprocess.run(
        #    ['reg', 'save', 'HKEY_CURRENT_USER\Environment',
        #     registry_bk, '/y'])
        #res.check_returncode()
        simnibs_env_vars = _get_win_simnibs_env_vars()
        path = _get_win_path()
        path = path.split(';')
        path = [p for p in path if len(p) > 0]
        for key, value in simnibs_env_vars.items():
            # If the directory is in the PATH variable, remove it
            path = [
                os.path.normpath(p) for p in path if not (
                os.path.normpath(value) in os.path.normpath(p))]
            # Remove environment variable
            res = subprocess.run(
                'reg delete HKEY_CURRENT_USER\Environment '
                f'/v  "{key}" /f', shell=True)
            res.check_returncode()
        path = ';'.join(path) + ';'
        # write out the PATH with SimNIBS removed
        res = subprocess.run(
            'reg add HKEY_CURRENT_USER\Environment '
            f'/v Path /t REG_EXPAND_SZ /d "{path}" /f',
            shell=True)
        res.check_returncode()

def matlab_setup(install_dir):
    destdir =  os.path.abspath(os.path.join(install_dir, 'matlab'))
    if not os.path.isdir(destdir):
        os.mkdir(destdir)
    for m in glob.glob(os.path.join(SIMNIBSDIR, 'matlab', '*')):
        _copy_and_log(m, destdir)

def setup_shortcut_icons(scripts_dir, force=False, silent=False):
    ''' Creates shortcut icons for the gui_scripts '''
    if sys.platform == 'win32':
        shortcut_folder = os.path.join(
            os.environ['APPDATA'],
            "Microsoft\Windows\Start Menu\Programs\SimNIBS"
        )
    elif sys.platform == 'linux':
        shortcut_folder = os.path.expanduser('~/.local/share/applications/SimNIBS')

    if os.path.isdir(shortcut_folder):
        overwrite = _get_input(
            'Found other SimNIBS menu icons, overwrite them?',
            silent)

        if not overwrite:
            print('Not adding shortucts to the current SimNIBS install')
            return
        
        else:
            shortcut_icons_clenup()

    os.mkdir(shortcut_folder)

    _create_shortcut(
        os.path.join(shortcut_folder, 'Gmsh'),
        file_finder.path2bin('gmsh'),
        mime_type='model/x.stl-binary'
    )

    _create_shortcut(
        os.path.join(shortcut_folder, 'SimNIBS GUI'),
        os.path.join(scripts_dir, 'simnibs_gui'),
        icon=os.path.join(SIMNIBSDIR, 'resources', 'gui_icon.bmp')
    )
    if sys.platform == 'linux':
        subprocess.run(
            ['update-desktop-database',
             os.path.expanduser('~/.local/share/applications')])

def _create_shortcut(shortcut_name, target_path, icon=None, mime_type=None):
    if sys.platform == 'win32':
        with tempfile.NamedTemporaryFile('w', delete=False, suffix='.ps1') as f:
            f.write(f'$objShell = New-Object -ComObject ("WScript.Shell")\n')
            f.write(f'$objShortCut = $objShell.CreateShortcut("{shortcut_name}.lnk")\n')
            f.write(f'$objShortCut.TargetPath="{target_path}"\n')
            f.write(f'$objShortCut.WorkingDirectory="%HOMEPATH%"\n')
            if icon:
                f.write(f'$objShortCut.IconLocation="{icon}"\n')
            f.write(f'$objShortCut.Save()')
            temp_fn = f.name
        subprocess.run(f'powershell.exe {temp_fn}', shell=True).check_returncode()
        os.remove(temp_fn)
    elif sys.platform == 'linux':
        with open(shortcut_name + '.desktop', 'w') as f:
            f.write('[Desktop Entry]\n')
            f.write(f'Name={os.path.basename(shortcut_name)}\n')
            f.write('Comment=SimNIBS is software '
                    'for simulating electric fields caused by NIBS\n')
            f.write(f'Exec={target_path} %f\n')
            if icon:
                f.write(f'Icon={icon}\n')
            if mime_type:
                f.write(f'MimeType={mime_type}\n')
            f.write('Terminal=false\n')
            f.write('Type=Application\n')

    # TODO: Mac version

def shortcut_icons_clenup():
    if sys.platform == 'win32':
        shortcut_folder=os.path.join(
            os.environ['APPDATA'],
            "Microsoft\Windows\Start Menu\Programs\SimNIBS"
        )
        shutil.rmtree(shortcut_folder)
        # wait
        while os.path.exists(shortcut_folder):
            pass
    elif sys.platform == 'linux':
        shortcut_folder = os.path.expanduser('~/.local/share/applications/SimNIBS')
        shutil.rmtree(shortcut_folder)

def setup_file_association(force=False, silent=False):
    # Linux file associations are done together with desktop items
    if sys.platform == 'linux':
        return
    gmsh_bin = file_finder.path2bin('gmsh')
    extensions = ['.msh', '.geo', '.stl']
    associate = dict.fromkeys(extensions)
    for ext in extensions:
        if _is_associated(ext):
            associate[ext] = _get_input(
                f'Found other association for "{ext}" files, overwrite it?',
                silent)
        else:
            associate[ext] = True
    # If all rejected, return
    if not any(associate.values()):
        return

    if sys.platform == 'win32':
        # We need to run with admin privileges
        # So a write a .cmd script and run it with administrative privileges
        with tempfile.NamedTemporaryFile('w', delete=False, suffix='.cmd') as f:
            [f.write(f'assoc {ext}=gmsh.simnibs\n') for ext, val in associate.items() if val]
            f.write(f'ftype gmsh.simnibs="{gmsh_bin}" "%%1"')
            temp_fn = f.name
        # I need to run as shell for some reason
        ret = subprocess.run(
            'powershell.exe -Command '
            f'"Start-Process -Wait -WindowStyle Hidden -Verb RunAs -FilePath {temp_fn}"',
            shell=True)
        try:
            ret.check_returncode()
        except subprocess.CalledProcessError:
            print('Could not associate files')
        finally:
            os.remove(temp_fn)

def _is_associated(ext):
    if sys.platform == 'win32':
        # Also needs to be run as shell
        ret = subprocess.run('assoc '+ ext, shell=True, capture_output=True)
        if ret.returncode == 0:
            assoc = re.findall(f'{ext}=(.*)\r', ret.stdout.decode())
            if len(assoc) == 0:
                return False
            return assoc[0]
        else:
            return False

def file_associations_cleanup():
    # Linux file associations are done together with desktop items
    if sys.platform == 'linux':
        return
    extensions = ['.msh', '.geo', '.stl']
    associate = dict.fromkeys(extensions)
    for ext in extensions:
        ass = _is_associated(ext)
        associate[ext] = (ass and 'simnibs' in ass)
    if not any(associate.values()):
        return
    
    if sys.platform == 'win32':
        with tempfile.NamedTemporaryFile('w', delete=False, suffix='.cmd') as f:
            [f.write(f'assoc {ext}=\n') for ext, val in associate.items() if val]
            f.write(f'ftype gmsh.simnibs=')
            temp_fn = f.name

        # I need to run as shell for some reason
        ret = subprocess.run(
        'powershell.exe -Command '
            f'"Start-Process -Wait -WindowStyle Hidden -Verb RunAs -FilePath {temp_fn}"',
            shell=True)
        try:
            ret.check_returncode()
        except subprocess.CalledProcessError:
            print('Could not cleanup file associations')
        finally:
            os.remove(temp_fn)

def uninstaller_setup(install_dir):
    uninstaller = os.path.join(install_dir, 'uninstall_simnibs')
    if sys.platform == 'win32':
        _write_windows_cmd(
            os.path.join(SIMNIBSDIR, 'cli', 'simnibs_postinstall.py'),
            uninstaller, commands='-u')

        conda_uninstaller = os.path.join(
            install_dir,'miniconda3',
            'Uninstall-Miniconda3.exe')


        with open(uninstaller + '.cmd', 'a') as f:
            f.write('\n')
            f.write(f'"{conda_uninstaller}"" /S\n')
            f.write(f'rd /s /q "{install_dir}"\n')

        _create_shortcut(
            os.path.join(install_dir, 'Uninstall SimNIBS'),
            uninstaller,
            os.path.join(SIMNIBSDIR, 'resources', 'gui_icon.bmp')
        )
    else:
        _write_unix_sh(
            os.path.join(SIMNIBSDIR, 'cli', 'simnibs_postinstall.py'),
            uninstaller, commands='-u')
        with open(uninstaller + '.cmd', 'a') as f:
            f.write('\n')
            f.write(f'rm -rf {install_dir}')


def activator_setup(install_dir):
    activator = os.path.join(install_dir, 'activate_simnibs')
    if sys.platform == 'win32':
        _write_windows_cmd(
            os.path.join(SIMNIBSDIR, 'cli', 'simnibs_postinstall.py'),
            activator, gui=True, commands=f'-d {install_dir}')

        _create_shortcut(
            os.path.join(install_dir, 'Activate SimNIBS'),
            activator,
            os.path.join(SIMNIBSDIR, 'resources', 'gui_icon.bmp')
        )
    else:
        _write_unix_sh(
            os.path.join(SIMNIBSDIR, 'cli', 'simnibs_postinstall.py'),
            activator, commands=f'-d {install_dir}')


if GUI:
    class PostInstallGUI(QtWidgets.QDialog):
        def __init__(self,
                     install_dir,
                     copy_gmsh_options=True,
                     add_to_path=True,
                     add_shortcut_icons=True,
                     associate_files=True):
            super().__init__()
            self.install_dir = os.path.abspath(os.path.expanduser(install_dir))
            self.copy_gmsh_options = copy_gmsh_options
            self.add_to_path = add_to_path
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

            self.setWindowTitle('Post-Install Options')
            gui_icon = os.path.join(SIMNIBSDIR,'resources', 'gui_icon.ico')
            self.setWindowIcon(QtGui.QIcon(gui_icon))

        def set_copy_gmsh_options(self, new_value):
            self.copy_gmsh_options = new_value
            if not new_value:
                QtWidgets.QMessageBox.warning(
                    self, 'SimNIBS',
                    'Visualization of models and simulations will be affected')

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
                    'Models and Simulation results will not automatically open in Gmsh')


    def start_gui(simnibsdir):
        app = QtWidgets.QApplication(sys.argv)
        ex = PostInstallGUI(simnibsdir)
        ex.show()
        app.exec_()
        if ex.result():
            install(ex.install_dir,
                    False, False,
                    ex.copy_gmsh_options,
                    ex.add_to_path,
                    ex.add_shortcut_icons,
                    ex.associate_files)
        else:
            raise Exception('Post-installation cancelled by user')


def install(install_dir, force, silent,
            copy_gmsh_options=True,
            add_to_path=True,
            add_shortcut_icons=False,
            associate_files=False):
    install_dir = os.path.abspath(os.path.expanduser(install_dir))
    scripts_dir = os.path.join(install_dir, 'bin')
    copy_scripts(scripts_dir)
    if copy_gmsh_options:
        setup_gmsh_options(force, silent)
    if add_to_path:
        path_setup(scripts_dir, force, silent)
    if add_shortcut_icons:
        setup_shortcut_icons(scripts_dir, force, silent)
    if associate_files:
        setup_file_association(force, silent)
    matlab_setup(install_dir)
    activator_setup(install_dir)
    uninstaller_setup(install_dir)

def uninstall():
    path_cleanup()
    shortcut_icons_clenup()
    file_associations_cleanup()


def _get_default_dir():
    if sys.platform == 'win32':
        return os.path.join(os.environ['LOCALAPPDATA'], 'SimNIBS')
    elif sys.platform == 'linux':
        return os.path.join(os.environ['HOME'], 'SimNIBS')

def main():
    parser = argparse.ArgumentParser(prog="simnibs_postinstall",
                                     description="Optional post-installation procedures "
                                     "for SimNIBS ")
    parser.add_argument('-d', "--target_dir", required=False,
                        help="SimNIBS install directory",
                        default=_get_default_dir())
    parser.add_argument('-f', "--force", required=False, action='store_true',
                        help="Perform all install steps")
    parser.add_argument('-s', "--silent", required=False, action='store_true',
                        help="Silent mode, will install without the GUI")
    parser.add_argument('-u', "--uninstall", required=False, action='store_true',
                        help="Ignores all other arguments and uninstall SimNIBS")
    args = parser.parse_args(sys.argv[1:])
    if args.uninstall:
        uninstall()
        return
    install_dir = os.path.abspath(os.path.expanduser(args.target_dir))
    if not args.silent:
        if not GUI:
            raise ImportError(
                'Trying to run post-install script without PyQt istalled, '
                'please use the silent mode (simnibs_postinstall --help for '
                'more information')
        start_gui(install_dir)

    else:
        install(install_dir, args.force, args.silent)

if __name__ == '__main__':
    main()

