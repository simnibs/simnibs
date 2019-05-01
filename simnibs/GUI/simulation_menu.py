#!/usr/bin/python2.7

'''
    Menu with simulation options for SimNIBS
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.
    
    Copyright (C) 2018  Guilherme B Saturnino

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''


from PyQt5 import QtCore, QtGui, QtWidgets
import os
import sys
from .. import SIMNIBSDIR



##Pop-up menu for selecting options for the simulation
class SimulationOptionsDialog (QtWidgets.QDialog):
    def __init__(self, parent, session):
        super(SimulationOptionsDialog, self).__init__(parent)

        self.session = session


        self.fields_box = self.selectFields()

        self.options_box = self.selectOtherOptions()
        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok|QtWidgets.QDialogButtonBox.Cancel)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)


        mainLayout = QtWidgets.QVBoxLayout()

        mainLayout.addWidget(self.fields_box)
        mainLayout.addWidget(self.options_box)
        mainLayout.addWidget(self.button_box)

        self.setLayout(mainLayout)  

        self.setWindowTitle('Simulation Options')

    def selectFields(self):
        fields_box = QtWidgets.QGroupBox("Fields:")
        layout = QtWidgets.QHBoxLayout()

        v_check = QtWidgets.QCheckBox()
        v_check.setText('v')
        if 'v' in self.session.fields: v_check.toggle()
        v_check.toggled.connect(self.changeFields)
        layout.addWidget(v_check)

        E_check = QtWidgets.QCheckBox()
        E_check.setText('vector E')
        if 'E' in self.session.fields: E_check.toggle()
        E_check.toggled.connect(self.changeFields)
        layout.addWidget(E_check)

        e_check = QtWidgets.QCheckBox()
        e_check.setText('norm E')
        if 'e' in self.session.fields: e_check.toggle()
        e_check.toggled.connect(self.changeFields)
        layout.addWidget(e_check)

        J_check = QtWidgets.QCheckBox()
        J_check.setText('vector J')
        if 'J' in self.session.fields: J_check.toggle()
        J_check.toggled.connect(self.changeFields)
        layout.addWidget(J_check)

        j_check = QtWidgets.QCheckBox()
        j_check.setText('norm J')
        if 'j' in self.session.fields: j_check.toggle()
        j_check.toggled.connect(self.changeFields)
        layout.addWidget(j_check)

        s_check = QtWidgets.QCheckBox()
        s_check.setText('Conductivities')
        if 's' in self.session.fields: s_check.toggle()
        s_check.toggled.connect(self.changeFields)
        layout.addWidget(s_check)

        A_check = QtWidgets.QCheckBox()
        A_check.setText('dA/dt (TMS only)')
        if 'D' in self.session.fields: A_check.toggle()
        A_check.toggled.connect(self.changeFields)
        layout.addWidget(A_check)


        fields_box.setLayout(layout)

        fields_box.v_check = v_check
        fields_box.E_check = E_check
        fields_box.e_check = e_check
        fields_box.J_check = J_check
        fields_box.j_check = j_check
        fields_box.s_check = s_check
        fields_box.A_check = A_check

        return fields_box


    def changeFields(self):
        self.session.fields  = ''
        if self.fields_box.v_check.isChecked():
            self.session.fields +='v'
        if self.fields_box.e_check.isChecked():
            self.session.fields  +='e'
        if self.fields_box.E_check.isChecked():
            self.session.fields  +='E'
        if self.fields_box.j_check.isChecked():
            self.session.fields  +='j'
        if self.fields_box.J_check.isChecked():
            self.session.fields  +='J'
        if self.fields_box.s_check.isChecked():
            self.session.fields  +='s'
        if self.fields_box.A_check.isChecked():
            self.session.fields  +='D'

    def selectOtherOptions(self):
        options_box = QtWidgets.QGroupBox("Additional Options:")
        layout = QtWidgets.QHBoxLayout()

        open_gmsh_cb = QtWidgets.QCheckBox()
        open_gmsh_cb.setText('Open in Gmsh')
        if self.session.open_in_gmsh: open_gmsh_cb.toggle()
        open_gmsh_cb.toggled.connect(self.changeOptions)
        layout.addWidget(open_gmsh_cb,0 , QtCore.Qt.Alignment(1))

        map_to_surf_cb = QtWidgets.QCheckBox()
        map_to_surf_cb.setText('Interpolate to cortical surface')
        if self.session.map_to_surf: map_to_surf_cb.toggle()
        map_to_surf_cb.toggled.connect(self.changeOptions)
        layout.addWidget(map_to_surf_cb,0, QtCore.Qt.Alignment(1))

        map_to_fsavg_cb = QtWidgets.QCheckBox()
        map_to_fsavg_cb.setText('Transform to fsaverage space')
        if self.session.map_to_fsavg: map_to_fsavg_cb.toggle()
        map_to_fsavg_cb.toggled.connect(self.changeOptions)
        layout.addWidget(map_to_fsavg_cb,0, QtCore.Qt.Alignment(1))

        map_to_vol_cb = QtWidgets.QCheckBox()
        map_to_vol_cb.setText('Interpolate to a nifti volume')
        if self.session.map_to_vol: map_to_vol_cb.toggle()
        map_to_vol_cb.toggled.connect(self.changeOptions)
        layout.addWidget(map_to_vol_cb,0, QtCore.Qt.Alignment(1))

        map_to_MNI_cb = QtWidgets.QCheckBox()
        map_to_MNI_cb.setText('Transform to MNI space')
        if self.session.map_to_MNI: map_to_MNI_cb.toggle()
        map_to_MNI_cb.toggled.connect(self.changeOptions)
        layout.addWidget(map_to_MNI_cb,0, QtCore.Qt.Alignment(1))



        options_box.setLayout(layout)

        options_box.open_gmsh_cb = open_gmsh_cb
        options_box.map_to_MNI_cb = map_to_MNI_cb
        options_box.map_to_vol_cb = map_to_vol_cb
        options_box.map_to_fsavg_cb = map_to_fsavg_cb
        options_box.map_to_surf_cb = map_to_surf_cb

        return options_box

    def changeOptions(self):
        self.session.open_in_gmsh = self.options_box.open_gmsh_cb.isChecked()
        self.session.map_to_surf = self.options_box.map_to_surf_cb.isChecked()
        self.session.map_to_fsavg = self.options_box.map_to_fsavg_cb.isChecked()
        self.session.map_to_MNI = self.options_box.map_to_MNI_cb.isChecked()
        self.session.map_to_vol = self.options_box.map_to_vol_cb.isChecked()


    def getOptions(self):
        result = self.exec_()
        return (self.session, result == QtWidgets.QDialog.Accepted)

        
 
#Dialog for selecting tensor files
class tensorFilesDialog(QtWidgets.QDialog):
    def __init__(self, parent, fname_tensor):
        super(tensorFilesDialog, self).__init__(parent)

        self.fname = fname_tensor

        groupBox = QtWidgets.QGroupBox('Tensor File')
        layout = QtWidgets.QHBoxLayout()

        groupBox.lineEdit = QtWidgets.QLineEdit()
        if self.fname is not None:
            groupBox.lineEdit.setText(self.fname)
        layout.addWidget(groupBox.lineEdit)

        groupBox.selectFile = QtWidgets.QPushButton('&Browse')
        groupBox.selectFile.clicked.connect(lambda: self.selectFile())
        layout.addWidget(groupBox.selectFile)

        groupBox.setLayout(layout)
        self.group_box = groupBox

        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok|QtWidgets.QDialogButtonBox.Cancel)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)

        mainLayout = QtWidgets.QGridLayout()

        mainLayout.addWidget(self.group_box)
        mainLayout.addWidget(self.button_box)

        self.setLayout(mainLayout)  

        self.setWindowTitle('Tensor file names')
        self.resize(400,200)

    def selectFile(self):
        if self.fname is not None and os.path.isfile(self.fname):
            directory = self.fname
        else:
            directory = QtCore.QDir.currentPath()
        dialog = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select tensor conductivity file', directory,  'Tensor files (*.nii *.nii.gz)')
        if dialog[0] != 0:
            self.fname = str(dialog[0])
            self.group_box.lineEdit.setText(str(dialog[0]))

    def getFileNames(self):
        result = self.exec_()
        return (self.fname, result == QtWidgets.QDialog.Accepted)


class EEGFileDialog(QtWidgets.QDialog):
    def __init__(self, parent, eeg_cap):
        super(EEGFileDialog, self).__init__(parent)

        self.fname = eeg_cap

        groupBox = QtWidgets.QGroupBox('EEG Cap File')
        layout = QtWidgets.QHBoxLayout()

        groupBox.lineEdit = QtWidgets.QLineEdit()
        if self.fname is not None:
            groupBox.lineEdit.setText(self.fname)
        layout.addWidget(groupBox.lineEdit)

        groupBox.selectFile = QtWidgets.QPushButton('&Browse')
        groupBox.selectFile.clicked.connect(self.selectFile)
        layout.addWidget(groupBox.selectFile)

        groupBox.setLayout(layout)
        self.group_box = groupBox

        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok|QtWidgets.QDialogButtonBox.Cancel)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)

        mainLayout = QtWidgets.QGridLayout()

        mainLayout.addWidget(self.group_box)
        mainLayout.addWidget(self.button_box)

        self.setLayout(mainLayout)  

        self.setWindowTitle('Tensor file names')
        self.resize(400,200)

    def selectFile(self):
        if self.fname is not None and os.path.isfile(self.fname):
            eeg_cap_dir = os.path.dirname(self.fname)
        else:
            eeg_cap_dir = QtCore.QDir.currentPath()
        dialog = QtWidgets.QFileDialog(self)
        dialog.setWindowTitle('Open EEG Position file')
        dialog.setNameFilter('(*.csv)')
        dialog.setDirectory(eeg_cap_dir)
        dialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        filename = None
        if dialog.exec_() == QtWidgets.QDialog.Accepted:
            filename = dialog.selectedFiles()
        if filename:
            self.fname = str(filename[0])
            self.group_box.lineEdit.setText(self.fname)


    def getFileNames(self):
        result = self.exec_()
        return (self.fname, result == QtWidgets.QDialog.Accepted)


@QtCore.pyqtSlot(int)
class SimProgressScreen (QtWidgets.QMainWindow):
    terminate_signal = QtCore.pyqtSignal()

    def __init__(self):
        super(SimProgressScreen, self).__init__()

        self.text = ''
        self.simFinished = False

        self.textBox = QtWidgets.QTextEdit()
        self.textBox.setReadOnly(True)
        self.textBox.setAcceptRichText(True)


        self.terminate_btn = QtWidgets.QPushButton('Terminate')
        self.terminate_btn.clicked.connect(self.close)

        mainLayout = QtWidgets.QGridLayout()

        mainLayout.addWidget(self.textBox)
        mainLayout.addWidget(self.terminate_btn)

        self.central_widget = QtWidgets.QWidget()
        self.central_widget.setLayout(mainLayout)


        self.setCentralWidget(self.central_widget)


        self.resize(800,500)

        self.setWindowTitle('Simulation Progress')
        try:
            gui_icon = os.path.join(SIMNIBSDIR,'resources', 'gui_icon.gif')
            self.setWindowIcon(QtGui.QIcon(gui_icon))
        except:
            pass
            

    def appendText(self, text):
        self.textBox.append(text)
        QtWidgets.QApplication.processEvents()


    def showSimProgress(self):
        self.show()

    def terminate(self):
        if self.simFinished:
            self.terminate_signal.emit()

        msgBox = QtWidgets.QMessageBox(
            QtWidgets.QMessageBox.Warning, 'Warning',
            "Are you sure?", QtWidgets.QMessageBox.NoButton, self)
        msgBox.addButton("Terminate", QtWidgets.QMessageBox.AcceptRole)
        msgBox.addButton("Continue", QtWidgets.QMessageBox.RejectRole)
        if msgBox.exec_() == QtWidgets.QMessageBox.AcceptRole:
            self.terminate_signal.emit()
            return True
        else:
            return False


    def setSimFinished(self):
        self.simFinished = True
        self.terminate_btn.setText('Close')



    def closeEvent(self, event):
        if self.simFinished:
            event.accept()
        else:
            if self.terminate():
                event.accept()
            else:
                event.ignore()
