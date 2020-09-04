# -*- coding: utf-8 -*-
'''
    Main GUI for SimNIBS

    Copyright (C) 2018 Guilherme B Saturnino

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
from PyQt5 import QtCore, QtGui, QtWidgets, QtPrintSupport
import copy
import sys
import os
import time
import logging
import numpy as np

from .. import SIMNIBSDIR
from ..simulation import sim_struct
from ..simulation.cond import standard_cond
from ..simulation import coil_numpy
from ..simulation.run_simnibs import run_simnibs
from ..msh import transformations
from .. import msh
from . import electrodeGUI
from . import head_model_OGL
from . import simulation_menu
from .. import __version__
from ..utils.simnibs_logger import logger
from ..utils.file_finder import SubjectFiles, path2bin
from ..utils.matlab_read import read_mat


class TDCS_GUI(QtWidgets.QMainWindow):

    def __init__(self):

        super(TDCS_GUI, self).__init__()

        self.selectFileLayout()
        #self.out_folder_box, self.out_folder_lineEdit = self.selectOutFolderLayout()
        self.table_rows = []
        self.session = sim_struct.SESSION()
        self.saveFn = ''



        self.simThreads = []
        self.progressScreens = []

        self.button_box = self.setButtonBox()
        self.createMenus()

        self.headModelWidget = head_model_OGL.HEADMODEL_UI()

        #self.electrodeTable()
        poslits_tabs = self.poslistTabs()

        left_widgets = QtWidgets.QGroupBox()
        left_widget_layout = QtWidgets.QGridLayout()

        left_widget_layout.addWidget(self.select_file,0,0,1,4)
        #left_widget_layout.addWidget(self.out_folder_box,1,0,1,4)
        #left_widget_layout.addWidget(self.electrodes)
        left_widget_layout.addWidget(poslits_tabs,2,0,4,4)
        left_widget_layout.addWidget(self.button_box)

        left_widgets.setLayout(left_widget_layout)

        #left_widgets.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.MinimumExpanding)
        #left_widgets.setFixedWidth(500)
        left_widgets.setMaximumWidth(500)
        left_widgets.setMinimumWidth(300)


        mainLayout = QtWidgets.QHBoxLayout()
        mainLayout.addWidget(left_widgets)
        mainLayout.addWidget(self.headModelWidget)

        self.central_widget = QtWidgets.QWidget()
        self.central_widget.setLayout(mainLayout)


        self.setCentralWidget(self.central_widget)

        self.setWindowTitle(f'SimNIBS {__version__}')

        try:
            gui_icon = os.path.join(SIMNIBSDIR,'resources', 'gui_icon.gif')
            self.setWindowIcon(QtGui.QIcon(gui_icon))
        except:
            pass

        self.resize(1200,1000)

    def sizeHint(self):
        return QtCore.QSize(1200,1000)

    def minumumSizeHint(self):
        return QtCore.QSize(500, 500)

    def selectFileLayout(self):

        self.select_file = QtWidgets.QGroupBox('')
        layout = QtWidgets.QGridLayout()

        tag_m2m = QtWidgets.QLabel('<b>m2m Folder:<\\b>')
        layout.addWidget(tag_m2m, 1, 0, 1, 3)

        self.m2m_folder_lineEdit = QtWidgets.QLineEdit()
        layout.addWidget(self.m2m_folder_lineEdit, 2, 0, 1, 3)
 
        file_browse_m2m = QtWidgets.QPushButton('Browse')
        file_browse_m2m.clicked.connect(self.m2mFolderDialog)
        layout.addWidget(file_browse_m2m,2,3,1,1)

        tag = QtWidgets.QLabel('<b>Head Mesh:<\\b>')
        layout.addWidget(tag,3,0,1,3)

        self.file_name = QtWidgets.QLineEdit()
        layout.addWidget(self.file_name, 4, 0, 1, 3)
 
        file_browse = QtWidgets.QPushButton('Browse')
        file_browse.clicked.connect(self.fileDialog)
        layout.addWidget(file_browse, 4, 3, 1, 1)

        tag_Out = QtWidgets.QLabel('<b>Output Folder:<\\b>')
        layout.addWidget(tag_Out,5,0,1,3)

        self.out_folder_lineEdit = QtWidgets.QLineEdit()
        layout.addWidget(self.out_folder_lineEdit,6,0,1,3)
 
        file_browse_out = QtWidgets.QPushButton('Browse')
        file_browse_out.clicked.connect(self.outFolderDialog)
        layout.addWidget(file_browse_out,6,3,1,1)



        self.select_file.setLayout(layout)

    def fileDialog(self):
        dialog = QtWidgets.QFileDialog(self)
        dialog.setWindowTitle('Open Mesh File')
        dialog.setNameFilter('gmsh files (*.msh)')
        dialog.setDirectory(QtCore.QDir.currentPath())
        dialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        filename = None
        if dialog.exec_() == QtWidgets.QDialog.Accepted:
            filename = dialog.selectedFiles()
        if filename:
            fn = str(filename[0])
            self.file_name.setText(fn)
            self.loadHeadModel(fn)
            self.lookForTensors(fn)

    def m2mFolderDialog(self):
        folder = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Directory"))
        if folder == '':
            return None
        else:
            self.m2m_folder_lineEdit.setText(folder)
            self.session.subpath = folder
            sub_files = SubjectFiles(subpath=folder)
            fn_mesh = sub_files.fnamehead
            if fn_mesh and not self.session.fnamehead:
                self.loadHeadModel(fn_mesh)



    #Defines the dialog for selecting output folder
    def outFolderDialog(self):
        folder = str( QtWidgets.QFileDialog.getExistingDirectory(self, "Select Directory"))
        if folder == '':
            return None
        else:
            self.out_folder_lineEdit.setText(folder)
            self.session.pathfem = folder
    #Loads the head model in a separate thread
    def loadHeadModel(self, fn):
        fn = os.path.abspath(fn)
        self.file_name.setText(fn)
        self.lookForTensors(fn)
        self.session.fnamehead = fn
        sub_files = SubjectFiles(fn, self.session.subpath)

        if os.path.isdir(sub_files.subpath):
            self.m2m_folder_lineEdit.setText(sub_files.subpath)
            self.session.subpath = sub_files.subpath

        if not self.session.eeg_cap:
            if os.path.isfile(sub_files.eeg_cap_1010):
                self.session.eeg_cap = sub_files.eeg_cap_1010

        if not self.session.pathfem:
            pathfem = os.path.join(sub_files.basedir, 'simnibs_simulation')
            self.out_folder_lineEdit.setText(pathfem)
            self.session.pathfem = pathfem

        self.loadThread = loadMeshThread(fn, self.headModelWidget.glHeadModel)
        self.loadThread.start()
        self.loadThread.finished.connect(self.drawModel)
        self.loadThread.finished.connect(self.changeGlStimulators)


    #Updates the OpenGL head model. The  function needs to be called here.
    def drawModel(self):
        QtWidgets.QApplication.processEvents()
        self.headModelWidget.glHeadModel.drawSkinAndGm()
        self.headModelWidget.glHeadModel.setEEG(self.session.eeg_cap)
        #self.headModelWidget.glHeadModel.selectRenderSurface('Scalp')

    #Defines the tabs
    def poslistTabs(self):
        poslistsTabs = QtWidgets.QGroupBox("Position Lists:")
        layout = QtWidgets.QGridLayout()
        self.poslistTabWidget = QtWidgets.QTabWidget()
        #elc_table = ElcTable(0,3,self.headModelWidget.glHeadModel,self)
        #self.poslistTabWidget.addTab(elc_table.generateGroupBox(), "tDCS")
        layout.addWidget(self.poslistTabWidget,0,0,3,0)

        add_tDCS_button = QtWidgets.QPushButton("Add tDCS Poslist")
        add_tDCS_button.clicked.connect(self.addTdcsPoslistTab)
        layout.addWidget(add_tDCS_button,3,0)

        add_TMS_button = QtWidgets.QPushButton("Add TMS Poslist")
        add_TMS_button.clicked.connect(self.addTmsPoslistTab)
        layout.addWidget(add_TMS_button,3,1)

        copy_button = QtWidgets.QPushButton("Copy Poslist")
        copy_button.clicked.connect(self.copyPoslist)
        layout.addWidget(copy_button,3,2)

        self.poslistTabWidget.setTabsClosable(True)
        self.poslistTabWidget.tabCloseRequested.connect(self.removePoslistTab)

        self.poslistTabWidget.currentChanged.connect(self.changeGlStimulators)
        #self.pos

        poslistsTabs.setLayout(layout)

        return poslistsTabs

    #Adds a tdcs-type tab to the poslist tabs
    def addTdcsPoslistTab(self, tdcslist=None):
        if not tdcslist:
            tdcslist = sim_struct.TDCSLIST()
            self.session.poslists.append(tdcslist)
        elc_table = ElcTable(self.headModelWidget.glHeadModel,
                             tdcslist, self,
                             eeg_cap=self.session.eeg_cap)
        self.poslistTabWidget.addTab(elc_table, "tDCS")
        # This command is buggy: it changes the first tab
        self.poslistTabWidget.setCurrentWidget(elc_table)

    #adds a tms-type tab to the poslist tabs
    def addTmsPoslistTab(self, tmslist=None):
        if not tmslist:
            tmslist = sim_struct.TMSLIST()
            self.session.poslists.append(tmslist)
        try:
            while not self.loadThread.isFinished():
                time.sleep(1)
        except:
            pass
        coil_table = CoilTable(self.headModelWidget.glHeadModel,
                               tmslist, self, eeg_cap=self.session.eeg_cap)
        self.poslistTabWidget.addTab(coil_table, "TMS")
        self.poslistTabWidget.setCurrentWidget(coil_table)


    def copyPoslist(self):
        ret = QtWidgets.QInputDialog.getInt(self, 'Select Poslist Number', 'poslist', 1,
                                            1, self.poslistTabWidget.count())
        if ret[1] == False:
            return None
        else:
            tab = self.poslistTabWidget.widget(ret[0]-1)
            if tab.type == 'tDCS':
                tdcslist = copy.deepcopy(tab.returnElAndCond())
                self.addTdcsPoslistTab(tdcslist)
                self.session.poslists.append(tdcslist)
            elif tab.type == 'TMS':
                tmslist = copy.deepcopy(tab.returnCoilAndConds())
                self.addTmsPoslistTab(tmslist)
                self.session.poslists.append(tmslist)
            self.changeGlStimulators()


    def removePoslistTab(self, index):
        '''
        ret = QtWidgets.QMessageBox.warning(self, "Warning",
                'Are you sure you want to delete poslist ' + str(index+1) + '?')
        if ret == False:
            return
        '''
        widget = self.poslistTabWidget.widget(index)
        widget.deleteLater()
        self.poslistTabWidget.removeTab(index)
        del self.session.poslists[index]

    #Changes the stimulators (coils and electrodes) being shown in the OpenGL window
    def changeGlStimulators(self):
        self.headModelWidget.glHeadModel.cleardAdtFields()
        widget = self.poslistTabWidget.currentWidget()
        if widget is not None:
            widget.updateStimulatorModels()

    #Creates the buttons "Run", "Save" and "Close"
    def setButtonBox(self):
        runButton = QtWidgets.QPushButton("Run")
        runButton.setDefault(True)
        runButton.clicked.connect(self.runSimnibs)

        button_box = QtWidgets.QDialogButtonBox()
        button_box.addButton(runButton, QtWidgets.QDialogButtonBox.ActionRole)

        return button_box

    #Creates the menu bar
    def createMenus(self):
        menu = self.menuBar()


        openAction = QtWidgets.QAction('Open', self)
        openAction.triggered.connect(self.openSimnibsFile)

        #Opens a simulation in GMSH
        #openSimAction = QtWidgets.QAction('Open Simulation', self)
        #openSimAction.triggered.connect(self.openSimulation)

        saveAction = QtWidgets.QAction('Save', self)
        saveAction.triggered.connect(self.saveSimFile)

        exitAction = QtWidgets.QAction('Exit', self)
        exitAction.triggered.connect(self.checkAndClose)

        fileMenu = menu.addMenu('File')
        fileMenu.addAction(openAction)
        #fileMenu.addAction(openSimAction)
        fileMenu.addAction(saveAction)
        fileMenu.addSeparator()
        fileMenu.addAction(exitAction)

        sim_options = QtWidgets.QAction('Simulation Options', self)
        sim_options.triggered.connect(self.setSimOptions)

        tensor_fns = QtWidgets.QAction('Select Tensor Files', self)
        tensor_fns.triggered.connect(self.selectTensorFiles)

        eeg_fn = QtWidgets.QAction('Select EEG Cap', self)
        eeg_fn.triggered.connect(self.selectEEGCap)

        select_colors = QtWidgets.QAction('Select Model Colors', self)
        select_colors.triggered.connect(self.selectColors)

        editMenu = menu.addMenu('Edit')
        editMenu.addAction(sim_options)
        editMenu.addAction(tensor_fns)
        editMenu.addAction(eeg_fn)
        editMenu.addAction(select_colors)

        licenseAction = QtWidgets.QAction('License', self)
        licenseAction.triggered.connect(self.licencePopup)

        aboutMenu = menu.addMenu('About')
        aboutMenu.addAction(licenseAction)


    def licencePopup(self):

        QtWidgets.QMessageBox.information(self, "Licensing", 
    '<p><center>     <b>SimNIBS</b>    </center></p>' +\
    '<p><center>      version ' + __version__ +'      </center></p>'+\
    '<p><center>Simulation of eletromagnetic fields generated by tDCS and TMS</center></p>' +\
    '<p><center>Copyright (C) 2019 SimNIBS Developers </center></p>'+\
    '<p>This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by'+\
    ' the Free Software Foundation, either version 3 of the License, or any later version.</p>'+\
    '<p>This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'+\
    ' See the GNU General Public License for more details.</p>'+\
    '<p>You should have received a copy of the GNU General Public License along with this program.  If not, see http://www.gnu.org/licenses/ </p>')
    

    #Calls the simulation option GUI
    def setSimOptions(self):
        SimGUI = simulation_menu.SimulationOptionsDialog(
            self, copy.deepcopy(self.session))
        options, ok = SimGUI.getOptions()
        if ok:
            self.session = options


    def selectTensorFiles(self):
        tensorDialog = simulation_menu.tensorFilesDialog(self, self.session.fname_tensor)
        fname_tensor, ok = tensorDialog.getFileNames()
        if ok:
            self.session.fname_tensor = fname_tensor

    def selectEEGCap(self):
        eeg_cap_dialog = simulation_menu.EEGFileDialog(self, self.session.eeg_cap)
        eeg_cap, ok = eeg_cap_dialog.getFileNames()
        if ok:
            self.session.eeg_cap = eeg_cap
            self.drawModel()

    def lookForTensors(self, fn):
        if os.path.isfile(str(self.session.fname_tensor)):
            return
        sub_files = SubjectFiles(fn)
        if os.path.isfile(sub_files.tensor_file):
            self.session.fname_tensor = sub_files.tensor_file


    def selectColors(self):
        items = ("Scalp", "Gray Matter")
        item, ok = QtWidgets.QInputDialog.getItem(self, "Color",
                "Select Surface:", items, 0, False)
        if not ok:
            return
        else:
            color = QtWidgets.QColorDialog.getColor()
            if color.isValid():
                tup = color.getRgb()
                rgb = [float(v)/255 for v in tup[0:3]]
                if  item == "Scalp":
                    self.headModelWidget.glHeadModel.changeSurfaceColors('Skin', rgb)
                if item == "Gray Matter":
                    self.headModelWidget.glHeadModel.changeSurfaceColors('GM', rgb)
            else:
                return

    #Open a SimNIBS file and fills out the fields
    def openSimnibsFile(self, fn=None):
        if fn:
            if not os.path.isfile(fn):
                print("ERROR: could not open file", fn)
                return None

            file_full_path = os.path.abspath(fn)
            path = os.path.split(file_full_path)[0]
            fn_no_extension, file_extension = os.path.splitext(file_full_path)

        else:
            dialog = QtWidgets.QFileDialog(self)
            dialog.setWindowTitle('Open SimNIBS File')
            dialog.setNameFilters(['MATLAB SimNIBS Configuration files (*.mat)'])
            dialog.setDirectory(QtCore.QDir.currentPath())
            dialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)
            if dialog.exec_() == QtWidgets.QDialog.Accepted:
                file_full_path = str(dialog.selectedFiles()[0])
                fn_no_extension, file_extension = os.path.splitext(file_full_path)

            else:
                return None

        self.saveFn = file_full_path
        
        S = read_mat(file_full_path)
        self.session = S

        count = self.poslistTabWidget.count()

        for idx in range(count):
            self.removePoslistTab(0)

        self.file_name.setText(S.fnamehead)
        if S.pathfem:
            self.out_folder_lineEdit.setText(S.pathfem)
        if S.subpath:
            self.m2m_folder_lineEdit.setText(S.subpath)
        if self.session.fnamehead:
            if os.path.isfile(self.session.fnamehead):
                self.loadHeadModel(S.fnamehead)
            else:
                QtWidgets.QMessageBox.critical(self, "Warning",
                    'Could not open file\n'+str(self.file_name.text()))

        for poslist in S.poslists:
            if poslist.type == 'TDCSLIST':
                self.addTdcsPoslistTab(poslist)
            elif poslist.type == 'TMSLIST':
                self.addTmsPoslistTab(poslist)


    #Opens a .msh file in gmsh
    def openSimulation(self):
        dialog = QtWidgets.QFileDialog(self)
        dialog.setWindowTitle('Open GMSH File')
        dialog.setNameFilter('GMSH files (*.msh)')
        dialog.setDirectory(QtCore.QDir.currentPath())
        dialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        if dialog.exec_() == QtWidgets.QDialog.Accepted:
            file_full_path = str(dialog.selectedFiles()[0])

        else:
            return None

        self.thread = openGmshThread(file_full_path)
        self.thread.start()



    #Generates a sim_struct.session() structure, for saving and running
    def generateNNAVSession(self):
        self.session.fnamehead = str(self.file_name.text())
        self.session.subpath = str(self.m2m_folder_lineEdit.text())
        self.session.pathfem = str(self.out_folder_lineEdit.text())

        self.session.poslists = []
        tab_count = self.poslistTabWidget.count()
        for index in range(tab_count):
            widget = self.poslistTabWidget.widget(index)
            if widget.type == 'tDCS':
                self.session.poslists.append(widget.returnElAndCond())
            if widget.type == 'TMS':
                self.session.poslists.append(widget.returnCoilAndConds())
        return self.session

    def saveSimFile(self, fn=None):
        S = self.generateNNAVSession()
        if not S:
            return None, None

        try:
            full_path = os.path.abspath(fn)
        except TypeError:
            if os.path.isfile(self.saveFn):
                directory = self.saveFn
            else:
                directory = QtCore.QDir.currentPath()

            dialog = QtWidgets.QFileDialog.getSaveFileName(
                self,'Save SimNIBS configuration file', directory,
                'MATLAB file (*.mat)')

            if dialog[0] != '':
                fn, extension = dialog
            else:
                return None, None

            if ' ' in fn:
                QtWidgets.QMessageBox.critical(
                        self, 'warning',  'Invalid file name:\n' +
                        'There are white space in the path:' +
                        '\n{0}'.format(fn))
                return None, None

            path, fn = os.path.split(fn)
            name, ext = os.path.splitext(fn)
            if ext == '':
                fn = name+'.mat'
            full_path = os.path.join(path, fn)
            if not os.path.exists(path):
                os.makedirs(path)
            self.saveFn = full_path

        sim_struct.save_matlab_sim_struct(S, full_path)

        return full_path, S

    #Checks if the user really wants to exit without saving
    def checkAndClose(self):
        msgBox = QtWidgets.QMessageBox()
        msgBox.setWindowTitle('SimNIBS')
        msgBox.setText("Save changes before exiting?")
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Save |  QtWidgets.QMessageBox.Close)
        msgBox.setDefaultButton(QtWidgets.QMessageBox.Save)
        ret = msgBox.exec_()
        if ret == QtWidgets.QMessageBox.Save:
            self.saveSimFile()
            self.close()
        elif ret == QtWidgets.QMessageBox.Close:
            self.close()


    #Runs Simnibs in a new thread
    def runSimnibs(self):
        S = self.generateNNAVSession()
        if self.warnings(S):
            return None

        self.simThreads.append(runSimThread(S))

        self.progressScreens.append(simulation_menu.SimProgressScreen())
        self.progressScreens[-1].show()

        self.simThreads[-1].start()
        self.simThreads[-1].output_signal.connect(self.progressScreens[-1].appendText)
        self.simThreads[-1].finished.connect(self.progressScreens[-1].setSimFinished)
        self.progressScreens[-1].terminate_signal.connect(self.simThreads[-1].set_stop)
        #self.progressScreens[-1].terminate_signal.connect(self.simThreads[-1].terminate)

    #Looks for problems before saving the poslist
    def warnings(self, S):

        problems = False
        if not os.path.isfile(S.fnamehead):
            QtWidgets.QMessageBox.critical(self, "Warning",
                'Invalid head mesh file')
            problems = True

        full_path = S.fnamehead
        path,fn = os.path.split(full_path)
        name, extension = os.path.splitext(fn)
        if extension != '.msh':
            QtWidgets.QMessageBox.critical(self, "Warning",
                'Please select a .msh file')
            problems = True

        if len(S.poslists) == 0:
            QtWidgets.QMessageBox.critical(self, "Warning",
                'Please define at least one TMS or tDCS Poslist')
            problems = True


        for ii, poslist in enumerate(S.poslists):
            if (poslist.anisotropy_type in ['vn', 'dir', 'mc'] and not
                    os.path.isfile(self.session.fname_tensor)):
                    QtWidgets.QMessageBox.critical(self, "Warning",
                            "<b><font color = red>Error in poslist "+str(ii+1) + "</b></font>"+\
                            '<p>Could not find tensor file corresponding to anisotropy type</p>' + \
                            "<p> For setting an anisotropy file go to Edit-> Select Tensor File</p>")

                    problems = True


            if poslist.type == 'TDCSLIST':
                if not np.isclose(np.sum(poslist.currents), 0, atol=1e-3):
                    QtWidgets.QMessageBox.critical(self, "Warning",
                            "Error in poslist "+str(ii+1) + "\nThe sum of the currents going through each electrode must be 0")
                    problems = True

                if np.any(np.isclose(poslist.currents, 0, atol=1e-5)):
                    QtWidgets.QMessageBox.critical(self, "Warning",
                            "Error in poslist "+ str(ii+1)+ "\nOne of the electrodes has no currents assigned")
                    problems = True


                for i in range (len(poslist.electrode)):
                    if not poslist.electrode[i].dimensions:
                        QtWidgets.QMessageBox.critical(self, "Warning",
                            "Error in poslist "+ str(ii+1)+ "\nDouble click on the third column cell in order to define an electrode")
                        problems = True

                    if len(poslist.electrode[i].centre) == 0:
                        QtWidgets.QMessageBox.critical(self, "Warning",
                            "Error in poslist "+ str(ii+1)+ "\nDouble click on the second column cell in order to define the electrode's position")
                        problems = True

                if len(poslist.electrode) < 2:
                    QtWidgets.QMessageBox.critical(self, "Warning",
                        "Error in poslist "+ str(ii+1)+ "\nAt least 2 electrodes must be defined")
                    problems = True


            if poslist.type == 'TMSLIST':
                if len(poslist.pos) == 0:
                    QtWidgets.QMessageBox.critical(self, "Warning",
                        "Error in poslist "+ str(ii+1)+ "\nDefine at least one coil position by double clicking the position coloumn")
                    problems = True

                if not os.path.isfile(poslist.fnamecoil):
                    QtWidgets.QMessageBox.critical(self, "Warning",
                        "Error in poslist "+ str(ii+1)+ "\nInvalid coil definition file!")
                    problems = True



        return problems





#tDCS Poslist tab
class ElcTable(QtWidgets.QWidget):
    def __init__(self, glHeadModel, tdcslist, parent, eeg_cap=None):
        super(ElcTable, self).__init__(parent)

        self.type = 'tDCS'
        self.glHeadModel = glHeadModel
        self.table_rows = []
        self.shapesOn = False

        self.colors =  [QtGui.QColor.fromCmykF(0., 0., 1., 0.),
                        QtGui.QColor.fromCmykF(0.72, 0.52, 0., 0.),
                        QtGui.QColor.fromCmykF(0., 1., 1., 0.),
                        QtGui.QColor.fromCmykF(1., 0., 1., 0.),
                        QtGui.QColor.fromCmykF(0., 0., 0., 1.0),
                        QtGui.QColor.fromCmykF(0., 0., 0., 0.34),
                        QtGui.QColor.fromCmykF(0., 1.0, 0., 0.5),
                        QtGui.QColor.fromCmykF(0., 1., 1., 0.45),
                        QtGui.QColor.fromCmykF(0., 0.09, 0.27, 0.04)]

        self.table = self.createTable()
        self.addNewRow()

        self.add_elt_btn = QtWidgets.QPushButton('Add Electrode')
        self.add_elt_btn.clicked.connect(self.addNewRow)

        self.rem_elt_btn = QtWidgets.QPushButton('Remove Electrode')
        self.rem_elt_btn.clicked.connect(self.deleteRow)

        self.preview_btn = QtWidgets.QPushButton('Preview Shapes')
        self.preview_btn.clicked.connect(self.showShapes)

        self.hide_btn = QtWidgets.QPushButton('Hide Shapes')
        self.hide_btn.clicked.connect(self.hideShapes)

        self.conduct_btn = QtWidgets.QPushButton("Set Conductivities")
        self.conduct_btn.clicked.connect(self.setConductivities)



        self.setRightClickMenu()

        layout = QtWidgets.QGridLayout()

        layout.addWidget(self.table,0,0,3,0)
        layout.addWidget(self.add_elt_btn,3,0)
        layout.addWidget(self.rem_elt_btn,3,1)
        layout.addWidget(self.preview_btn,3,2)
        layout.addWidget(self.hide_btn, 4,2)
        layout.addWidget(self.conduct_btn,4,1)  

        self.setLayout(layout)

        self.loadStruct(tdcslist, eeg_cap)

    #Defines the table itself
    def createTable(self):
        table =  QtWidgets.QTableWidget(0,4)
        table.setHorizontalHeaderLabels(("Current", "Position", "Shape", "Name"))
        #table.setColumnWidth(3,1)
        table.setShowGrid(False)
        #table.horizontalHeader().setResizeMode(2, QtWidgets.QHeaderView.Stretch)
        table.verticalHeader().hide()
        table.cellDoubleClicked.connect(self.tableLeftClick)
        return table

    def tableLeftClick(self, row, column):
        if column == 1:
            self.definePosition(row)
        if column == 2:
            self.defineElectrode(row)

    def addNewRow(self):
        row_nbr = self.table.rowCount()

        self.table_rows.append(ELC_TABLE_ROW())
        self.table.insertRow(row_nbr)

        self.table.setCellWidget(row_nbr, 0, self.table_rows[row_nbr].current_box)
        self.table.setItem(row_nbr, 1, self.table_rows[row_nbr].position_item)
        self.table.setItem(row_nbr, 2, self.table_rows[row_nbr].shape_size_item)
        self.table.setItem(row_nbr, 3, self.table_rows[row_nbr].name_item)
        self.table_rows[row_nbr].name_item.setBackground(self.colors[row_nbr%len(self.colors)])


    def deleteRow(self):
        currentRow = self.table.currentRow()
        if currentRow == -1:
            currentRow = self.table.rowCount()-1
        self.table.removeRow(currentRow)
        self.table_rows.pop(currentRow)
        self.updateStimulatorModels()
        for ii,row in enumerate(self.table_rows):
            row.name_item.setBackground(self.colors[ii%len(self.colors)])


    #Gets the electrode positions, uses the Position_GUI struct
    def definePosition(self, row):
         pos_gui = Position_GUI(self.table_rows[row].electrode.centre,
                                self.table_rows[row].electrode.pos_ydir,
                                str(self.table_rows[row].name_item.text()),
                                self.glHeadModel)
         pos_gui.show()
         centre, pos_ydir, name, ok = pos_gui.GetPositions()

         if ok and pos_ydir == [0.0,0.0,0.0]:
            QtWidgets.QMessageBox.critical(self, "Warning",
                'You must select a direction by checking the box and clicking the head model')
            self.definePosition(row)
            #centre, pos_ydir,ok = pos_gui.GetPositions()

         elif ok:
             self.table_rows[row].centre = centre
             self.table_rows[row].name = name
             self.table_rows[row].pos_ydir = pos_ydir
             self.table_rows[row].electrode.centre = centre
             self.table_rows[row].electrode.pos_ydir = pos_ydir
             self.writeOnPosColumn(row, centre)

             scalp_surf = self.glHeadModel.getSurface('Scalp')
             self.table_rows[row].transf_matrix = scalp_surf.calculateMatSimnibs(centre, pos_ydir)

             self.table_rows[row].name_item.setText(name)
             self.updateStimulatorModels()


    #Writes the position in the appropriate column
    def writeOnPosColumn(self, row, pos):
        if pos is not None and len(pos) == 3:
            pos_text = "{0:0.1f}, {1:0.1f}, {2:0.1f}".format(*pos)
            self.table_rows[row].position_item.setText(pos_text)

    #Uses the Ui_Electrode sctuct to retrieve information about the electrode
    def defineElectrode(self, row):
        electrode_gui = electrodeGUI.Ui_Electrode(self.table_rows[row].electrode)
        electrode, ok = electrode_gui.return_el_struct()
        electrode.channelnr = row+1

        if ok:
            self.table_rows[row].electrode = electrode
            self.writeOnElColumn(row, electrode)


    def writeOnElColumn(self, row, electrode):
        shape = ''
        if self.table_rows[row].electrode.shape == 'rect':
            shape = 'Rectangular'

        elif self.table_rows[row].electrode.shape == 'ellipse':
            shape = 'Elliptical'

        else:
            return
        string = str(self.table_rows[row].electrode.dimensions[0]/10)\
                + 'x' + str(self.table_rows[row].electrode.dimensions[1]/10) + \
                ' ' + shape
        self.table_rows[row].shape_size_item.setText(string)

    def add_conductivity_to_list(self):
        for i in range(len(self.table_rows)):
            self.tdcslist.cond[99+i].name = 'Electrode_%d_rubber'%(i+1)
            if not self.tdcslist.cond[99+i].value:
                self.tdcslist.cond[99+i].value = self.tdcslist.cond[99].value

            self.tdcslist.cond[499+i].name = 'Electrode_%d_saline_or_gel'%(i+1)
            if not self.tdcslist.cond[499+i].value:
                self.tdcslist.cond[499+i].value = self.tdcslist.cond[499].value

    #Calls the conduvtivitiesGui
    def setConductivities(self):
        self.add_conductivity_to_list()
        CondGui = ConductivitiesGui(copy.deepcopy(self.tdcslist))
        tdcslist, ok = CondGui.getCond()
        if ok:
            self.tdcslist = tdcslist


    #Loads form poslist
    def loadStruct(self, tdcslist, eeg_cap):
        self.tdcslist = tdcslist
        nbr_electrodes = len(tdcslist.electrode)
        self.deleteRow()
        for i in range(nbr_electrodes):
            self.addNewRow()
            tdcslist.electrode[i].substitute_positions_from_cap(eeg_cap)
            if tdcslist.electrode[i].pos_ydir is None or len(
                tdcslist.electrode[i].pos_ydir) == 0:
                while self.glHeadModel.getSurface('Scalp') is 'Loading':
                    time.sleep(1)
                tdcslist.electrode[i].pos_ydir = _get_posy(
                    tdcslist.electrode[i].centre,
                    self.glHeadModel.getSurface('Scalp'))
            self.table_rows[i].centre = tdcslist.electrode[i].centre
            self.table_rows[i].pos_ydir = tdcslist.electrode[i].pos_ydir
            self.table_rows[i].electrode = tdcslist.electrode[i]
            self.table_rows[i].current_box.setValue(tdcslist.currents[i]*1000)
            self.table_rows[i].name_item.setText(tdcslist.electrode[i].name)

            self.writeOnPosColumn(i, tdcslist.electrode[i].centre)
            self.writeOnElColumn(i, tdcslist.electrode[i])
        self.updateStimulatorModels()

    #Returns currents and electrode
    def returnElAndCond(self):
        nbr_electrodes = self.table.rowCount()
        self.tdcslist.currents = [0]*nbr_electrodes
        self.tdcslist.electrode = [0]*nbr_electrodes
        for i in range(nbr_electrodes):
            self.tdcslist.currents[i] = self.table_rows[i].current_box.value()/1000
            self.tdcslist.electrode[i] = self.table_rows[i].electrode
            self.tdcslist.electrode[i].name = str(self.table_rows[i].name_item.text())
        return self.tdcslist

    #Copying and pasting electrodes
    def setRightClickMenu (self):
        self.menu = QtWidgets.QMenu(self)
        self.copyAction = QtWidgets.QAction('Copy', self)
        self.copyAction.triggered.connect(self.copyElectrode)

        self.pasteAction = QtWidgets.QAction('Paste', self)
        self.pasteAction.triggered.connect(self.pasteElectrode)
        self.pasteAction.setEnabled(False)

        self.deleteAction = QtWidgets.QAction('Delete', self)
        self.deleteAction.triggered.connect(self.deleteElectrode)

        self.menu.addAction(self.copyAction)
        self.menu.addAction(self.pasteAction)
        self.menu.addAction(self.deleteAction)


    def contextMenuEvent(self, event):
        if self.table.currentColumn() == 2:
            self.menu.popup(QtGui.QCursor.pos())

    def copyElectrode(self):
        self.electrode_copy = copy.deepcopy(self.table_rows[self.table.currentRow()].electrode)
        self.pasteAction.setEnabled(True)

    def pasteElectrode(self):
        self.electrode_copy.channelnr = self.table.currentRow() + 1
        self.table_rows[self.table.currentRow()].electrode = copy.deepcopy(self.electrode_copy)
        self.table_rows[self.table.currentRow()].electrode.centre = \
            self.table_rows[self.table.currentRow()].centre
        self.table_rows[self.table.currentRow()].electrode.pos_ydir = \
            self.table_rows[self.table.currentRow()].pos_ydir
        self.writeOnElColumn(self.table.currentRow(), self.electrode_copy)
        self.updateStimulatorModels()

    def deleteElectrode(self):
        self.table_rows[self.currentRow()].electrode = sim_struct.ELECTRODE()
        self.table_rows[self.currentRow()].shape_size_item.setText('')

    #update the electrode objects in openGL
    def updateStimulatorModels(self):
        objects = []
        i = 0
        for row in self.table_rows:
            if row.transf_matrix != []:
                ogl_object = self.glHeadModel.drawPointAndDirs(row.transf_matrix, self.colors[i%len(self.colors)])
                objects.append(ogl_object)
            else:
                try:
                    scalp_surf = self.glHeadModel.getSurface('Scalp')
                    row.transf_matrix = scalp_surf.calculateMatSimnibs(row.electrode.centre, row.electrode.pos_ydir)
                    ogl_object = self.glHeadModel.drawPointAndDirs(row.transf_matrix, self.colors[i%len(self.colors)])
                    objects.append(ogl_object)
                except:
                    pass

            if self.shapesOn:
                center = row.electrode.centre
                if row.electrode.dimensions_sponge not in [None, []]:
                    dimensions = row.electrode.dimensions_sponge
                else:
                    dimensions = row.electrode.dimensions
                if center is not None and len(center) == 3 and dimensions is not None:
                    objects.append(self.glHeadModel.drawElectrode(row.electrode, self.colors[i%len(self.colors)]))
            i += 1

        self.glHeadModel.stimulatorList(objects)
        self.glHeadModel.update()

    def showShapes(self):
        self.shapesOn = True
        self.updateStimulatorModels()

    def hideShapes(self):
        self.shapesOn = False
        self.updateStimulatorModels()




#tDCS poslist tab
class ELC_TABLE_ROW:
    def __init__(self):
        self.electrode = sim_struct.ELECTRODE()
        #self.electrode.definition = 'plane'
        self.centre = None
        self.pos_ydir = None
        self.current_box = QtWidgets.QDoubleSpinBox()
        self.current_box.setDecimals(3)
        self.current_box.setSuffix("mA")
        self.current_box.setMinimum(-10)
        self.current_box.setMaximum(10)
        self.current_box.setSingleStep(0.25)
        self.current = self.current_box.value()

        self.position_item = QtWidgets.QTableWidgetItem('')
        self.position_item.setFlags(self.position_item.flags() ^ QtCore.Qt.ItemIsEditable)

        self.shape_size_item = QtWidgets.QTableWidgetItem('')
        self.shape_size_item.setFlags(self.shape_size_item.flags() ^ QtCore.Qt.ItemIsEditable)

        #self.color_item = QtWidgets.QTableWidgetItem('')

        self.name_item = QtWidgets.QTableWidgetItem('')

        self.transf_matrix = []


#Defies a row in the table on the TMS poslist tab
class COIL_TABLE_ROW:
    def __init__(self):
        #self.position_struc = sim_struct.POSITION()
        self.p1 = []
        self.p2 = []
        #self.ogl_object = None
        #self.electrode.definition = 'plane'
        self.didt_box = QtWidgets.QDoubleSpinBox()
        self.didt_box.setSuffix("x1e6 A/s")
        self.didt_box.setMinimum(0)
        self.didt_box.setMaximum(200)
        self.didt_box.setSingleStep(0.5)
        self.didt_box.setValue(1)

        self.dist_box = QtWidgets.QDoubleSpinBox()
        self.dist_box.setSuffix("mm")
        self.dist_box.setMinimum(0.)
        self.dist_box.setMaximum(200)
        self.dist_box.setSingleStep(1)
        self.dist_box.setValue(4)

        self.position_item = QtWidgets.QTableWidgetItem('')
        self.position_item.setFlags(self.position_item.flags() ^ QtCore.Qt.ItemIsEditable)

        self.name_item = QtWidgets.QTableWidgetItem('')


    def calc_matsimnibs(self, surface):
        return surface.calculateMatSimnibs(self.p1, self.p2, skin_distance=self.dist_box.value())

    def calc_distance(self, surface):
        return surface.calculateDistance(self.p1)


#TMS Poslist Tab
class CoilTable (QtWidgets.QWidget):
    def __init__(self, glHeadModel, tmslist, parent, eeg_cap=None):
        super(CoilTable, self).__init__(parent)

        self.type = 'TMS'
        self.glHeadModel = glHeadModel
        self.table_rows = []

        self.colors =  [QtGui.QColor.fromCmykF(0., 0., 1., 0.),
                        QtGui.QColor.fromCmykF(0.72, 0.52, 0., 0.),
                        QtGui.QColor.fromCmykF(0., 1., 1., 0.),
                        QtGui.QColor.fromCmykF(1., 0., 1., 0.),
                        QtGui.QColor.fromCmykF(0., 0., 0., 1.0),
                        QtGui.QColor.fromCmykF(0., 0., 0., 0.34),
                        QtGui.QColor.fromCmykF(0., 0., 0., 0.),
                        QtGui.QColor.fromCmykF(0., 1.0, 0., 0.5),
                        QtGui.QColor.fromCmykF(0., 1., 1., 0.45),
                        QtGui.QColor.fromCmykF(0., 0.09, 0.27, 0.04)]


        self.coil_box, self.coil_line_edit = self.coilBox()

        self.table = self.createTable()
        self.addNewRow()

        self.add_pos_btn = QtWidgets.QPushButton('Add Position')
        self.add_pos_btn.clicked.connect(self.addNewRow)

        self.rem_pos_btn = QtWidgets.QPushButton('Remove Position')
        self.rem_pos_btn.clicked.connect(self.deleteRow)

        self.conduct_btn = QtWidgets.QPushButton("Set Conductivities")
        self.conduct_btn.clicked.connect(self.setConductivities)


        self.show_dAdt_btn = QtWidgets.QPushButton("Show dA/dt field")
        self.show_dAdt_btn.clicked.connect(self.showdAdt)

        layout = QtWidgets.QGridLayout()

        layout.addWidget(self.coil_box,0,0,1,3)
        layout.addWidget(self.table,1,0,3,3)
        layout.addWidget(self.add_pos_btn,4,0)
        layout.addWidget(self.rem_pos_btn,4,1)
        layout.addWidget(self.conduct_btn,5,1)    
        layout.addWidget(self.show_dAdt_btn, 4,2)


        self.setLayout(layout)

        self.loadStruct(tmslist, eeg_cap)

    #Box for defining the electrode file
    def coilBox(self):
        box = QtWidgets.QGroupBox("Coil Definition File:")
        layout = QtWidgets.QGridLayout()

        coil_line_edit = QtWidgets.QLineEdit()
        coil_line_edit.textChanged.connect(self.set_coil_fn)
        layout.addWidget(coil_line_edit,0,0,1,3)

        file_browse_btn = QtWidgets.QPushButton('Browse')
        file_browse_btn.clicked.connect(self.coilDialog)
        layout.addWidget(file_browse_btn,0,3,1,1)

        box.setLayout(layout)

        return box, coil_line_edit



    def coilDialog(self):
        #get folder with ccd files
        try:
            ccd_folder = os.path.join(SIMNIBSDIR, 'ccd-files')
        except:
            ccd_folder = './'


        dialog = QtWidgets.QFileDialog(self)
        dialog.setWindowTitle('Open Coil Definition File')
        dialog.setNameFilter('Coil Definition files (*.ccd *.nii *.gz)')
        dialog.setDirectory(ccd_folder)
        dialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        if dialog.exec_() == QtWidgets.QDialog.Accepted:
            fn = str(dialog.selectedFiles()[0])
        else:
            return None

        self.tmslist.fnamecoil = fn
        self.coil_line_edit.setText(fn)

    def set_coil_fn(self):
        self.tmslist.fnamecoil = str(self.coil_line_edit.text())

    def createTable(self):
        table =  QtWidgets.QTableWidget(0,4)
        table.setHorizontalHeaderLabels(("dI/dt", "Skin Distance", "  Position  ", " Name "))
        table.setShowGrid(False)
        table.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
        table.verticalHeader().hide()
        table.cellDoubleClicked.connect(self.tableLeftClick)
        return table


    def tableLeftClick(self, row, column):
        if column == 2:
            self.definePosition(row)


    def addNewRow(self):
        row_nbr = self.table.rowCount()

        self.table_rows.append(COIL_TABLE_ROW())
        self.table.insertRow(row_nbr)

        self.table.setCellWidget(row_nbr, 0, self.table_rows[row_nbr].didt_box)
        self.table.setCellWidget(row_nbr, 1, self.table_rows[row_nbr].dist_box)
        self.table.setItem(row_nbr, 2, self.table_rows[row_nbr].position_item)
        self.table.setItem(row_nbr,3, self.table_rows[row_nbr].name_item)
        self.table_rows[row_nbr].name_item.setBackground(self.colors[row_nbr%len(self.colors)])
        self.table_rows[row_nbr].dist_box.valueChanged.connect(self.updateStimulatorModels)
        #self.table.setItem(row_nbr, 2, self.table_rows[row_nbr].coil)


    def deleteRow(self):
        currentRow = self.table.currentRow()
        if currentRow == -1:
            currentRow = self.table.rowCount()-1
        self.table.removeRow(currentRow)
        self.table_rows.pop(currentRow)
        self.updateStimulatorModels()
        for ii,row in enumerate(self.table_rows):
            row.name_item.setBackground(self.colors[ii%len(self.colors)])


    def definePosition(self, row):
         pos_gui = Position_GUI(self.table_rows[row].p1, self.table_rows[row].p2,
                                str(self.table_rows[row].name_item.text()),
                                self.glHeadModel)
         pos_gui.show()
         p1, p2, name, ok = pos_gui.GetPositions()

         #to make sure the reference is given
         while ok and p2 == [0.0,0.0,0.0]:
            QtWidgets.QMessageBox.critical(self, "Warning",
                'You must select a direction by checking the box and clicking the head model')
            p1, p2, name, ok = pos_gui.GetPositions()

         if ok:
             self.table_rows[row].p1 = p1
             self.table_rows[row].p2 = p2
             scalp_surf = self.glHeadModel.getSurface('Scalp')
             self.writeOnPosColumn(row, p1)
             self.table_rows[row].name_item.setText(name)
             self.updateStimulatorModels()

    def writeOnPosColumn(self, row, pos):
        if pos is not None:
            pos_text = "{0:0.1f}, {1:0.1f}, {2:0.1f}".format(*pos)
            self.table_rows[row].position_item.setText(pos_text)


    #Calls the conduvtivitiesGui
    def setConductivities(self):
        CondGui = ConductivitiesGui(copy.deepcopy(self.tmslist))
        tmslist, ok = CondGui.getCond()
        if ok:
            self.tmslist = tmslist


    #Loads poslist
    def loadStruct(self, poslist, eeg_cap):
        self.tmslist = poslist
        self.coil_line_edit.setText(self.tmslist.fnamecoil)
        nbr_positions = len(poslist.pos)
        self.deleteRow()
        for i in range(nbr_positions):
            self.addNewRow()
            try:
                matsimnibs = np.array(poslist.pos[i].matsimnibs)
            except:
                matsimnibs = None
            if matsimnibs is not None and matsimnibs.shape == (4, 4):
                self.table_rows[i].p1 = matsimnibs[:3, 3]
                self.table_rows[i].p2 = matsimnibs[:3, 1] + matsimnibs[:3, 3]
                surf = self.glHeadModel.getSurface('Scalp')
                if surf is not None:
                    self.table_rows[i].dist_box.setValue(
                        self.table_rows[i].calc_distance(surf))
                self.writeOnPosColumn(i, self.table_rows[i].p1)
            elif poslist.pos[i].centre is not None and len(poslist.pos[i].centre) > 0:
                poslist.pos[i].substitute_positions_from_cap(eeg_cap)
                self.table_rows[i].p1 = poslist.pos[i].centre
                self.table_rows[i].p2 = poslist.pos[i].pos_ydir
                self.table_rows[i].dist_box.setValue(
                    poslist.pos[i].distance)
                self.writeOnPosColumn(i, self.table_rows[i].p1)

            if poslist.pos[i].name is not None:
                self.table_rows[i].name_item.setText(poslist.pos[i].name)
            if poslist.pos[i].didt is not None:
                self.table_rows[i].didt_box.setValue(poslist.pos[i].didt*1e-6)



    def returnCoilAndConds(self):
        self.tmslist.pos = []
        for row in self.table_rows:
            self.tmslist.pos.append(sim_struct.POSITION())
            if len(row.p1) == 3:
                self.tmslist.pos[-1].matsimnibs = row.calc_matsimnibs(self.glHeadModel.getSurface('Scalp'))
            self.tmslist.pos[-1].didt = row.didt_box.value()*1e6
            self.tmslist.pos[-1].name = str(row.name_item.text())

        return self.tmslist


    def updateStimulatorModels(self):
        c_list = []
        i = 0
        for row in self.table_rows:
            if len(row.p1) != 0:
                ogl_object = self.glHeadModel.drawPointAndDirs(row.calc_matsimnibs(self.glHeadModel.getSurface('Scalp')),
                                                            self.colors[i%len(self.colors)])
                c_list.append(ogl_object)
                i += 1
        self.glHeadModel.stimulatorList(c_list)


    def showdAdt(self):
        if str(self.coil_line_edit.text()) == '' or str(self.coil_line_edit.text()).endswith('.ccd'):
            QtWidgets.QMessageBox.critical(self, "Warning",
                'For previewing dA/dt fields, please select a nifti coil file')
        row = self.table_rows[self.table.currentRow()]
        if row.p1 is None or len(row.p1) != 3:
            return
        self.glHeadModel.get_dAdtField(row.calc_matsimnibs(self.glHeadModel.getSurface('Scalp')),
                                       str(self.coil_line_edit.text()))
        #self.glHeadModel.drawModel()
        self.glHeadModel.update()


#GUI where the user clicks and chooses the position
class Position_GUI(QtWidgets.QDialog):
    def __init__(self, centre, reference, name, glHeadModel):
        super(Position_GUI, self).__init__()
        self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.centre = centre
        self.reference = reference
        self.glHeadModel = glHeadModel
        self.ogl_object = None
        self.name = name
        glHeadModel.windowClicked.connect(self.getClickedPosition)

        mainLayout = QtWidgets.QGridLayout()

        self.set_position_layout()
        self.set_direction_layout()
        #self.getClickedPosition()

        mainLayout.addWidget(self.position_box,0,0,1,2)
        mainLayout.addWidget(self.direction_box,1,0,1,2)

        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok|QtWidgets.QDialogButtonBox.Cancel)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)

        mainLayout.addWidget(self.button_box,2,1,1,3)


        self.project_button = QtWidgets.QPushButton('Project on scalp') 
        self.project_button.clicked.connect(self.project_points)
        mainLayout.addWidget(self.project_button,2,0,1,1)

        self.setLayout(mainLayout)
        self.setWindowTitle("Position")
        self.el_type = []
        self.el_coords = []
        self.el_name = []
        self.load_cap()
        #self.resize(200,300)

    def load_cap(self):
        if self.glHeadModel.eeg_coordinates is None:
            self.el_coords = []
            self.eeg_pos_selector.clear()
            self.eeg_pos_selector.setCurrentIndex(0)
            self.eeg_pos_selector.setEnabled(False)
        else:
            coords = np.array(self.glHeadModel.eeg_coordinates)
            eeg_names = np.array(self.glHeadModel.eeg_names)
            order = np.argsort(eeg_names)
            self.glHeadModel.eeg_names = eeg_names[order].tolist()
            self.glHeadModel.eeg_coordinates = coords[order].tolist()
            self.el_coords = [None] + list(self.glHeadModel.eeg_coordinates)
            self.el_name = [None] + self.glHeadModel.eeg_names
            if self.name in self.el_name:
                curr_idx = self.el_name.index(self.name)
            else:
                curr_idx = 0
            self.eeg_pos_selector.clear()
            self.eeg_pos_selector.addItems(self.el_name)
            self.eeg_pos_selector.setCurrentIndex(curr_idx)
            self.eeg_pos_selector.setEnabled(True)


    def set_position_layout(self):
        self.position_box = QtWidgets.QGroupBox("Position")
        load_pos = self.centre is not None and len(self.centre) == 3

        Layout = QtWidgets.QGridLayout()
        self.eeg_pos_selector = QtWidgets.QComboBox()
        self.eeg_pos_selector.activated.connect(self.select_eeg_pos)
        Layout.addWidget(self.eeg_pos_selector, 1,0)
        self.eeg_pos_selector.setEnabled(False)

        self.label_x = QtWidgets.QLabel("X:")
        Layout.addWidget(self.label_x,0,1)
        self.pos_x = QtWidgets.QDoubleSpinBox()
        self.pos_x.setRange(-1000,1000)
        if load_pos: self.pos_x.setValue(self.centre[0])
        #self.pos_x.valueChanged.connect(self.update_center)
        Layout.addWidget(self.pos_x,1,1)

        self.label_y = QtWidgets.QLabel("Y:")
        Layout.addWidget(self.label_y,0,2)
        self.pos_y = QtWidgets.QDoubleSpinBox()
        self.pos_y.setRange(-1000,1000)
        if load_pos: self.pos_y.setValue(self.centre[1])
       # self.pos_y.valueChanged.connect(self.update_center)
        Layout.addWidget(self.pos_y,1,2)

        self.label_z = QtWidgets.QLabel("Z:")
        Layout.addWidget(self.label_z,0,3)
        self.pos_z = QtWidgets.QDoubleSpinBox()
        self.pos_z.setRange(-1000,1000)
        if load_pos: self.pos_z.setValue(self.centre[2])
        #self.pos_z.valueChanged.connect(self.update_center)
        Layout.addWidget(self.pos_z,1,3)

        self.position_box.setLayout(Layout)

    def select_eeg_pos(self):
        i = self.eeg_pos_selector.currentIndex()
        c = self.el_coords[i]
        if c is not None:
            self.pos_x.setValue(c[0])
            self.pos_y.setValue(c[1])
            self.pos_z.setValue(c[2])
            self.name = self.el_name[i]
            if not self.check.isChecked():
                self.check.toggle()
        self.update_center()

    def set_direction_layout(self):
        box_text = "Direction Reference"
        self.direction_box = QtWidgets.QGroupBox(box_text)

        load_ref = self.reference is not None and len(self.reference) == 3

        Layout = QtWidgets.QGridLayout()

        self.check = QtWidgets.QCheckBox('Define Reference Coordinates')

        Layout.addWidget(self.check, 0,0,1,3)

        self.label_x2 = QtWidgets.QLabel("X:")
        Layout.addWidget(self.label_x2,1,1)
        self.ref_x = QtWidgets.QDoubleSpinBox()
        self.ref_x.setRange(-1000, 1000)
        if load_ref: self.ref_x.setValue(self.reference[0])
        #self.ref_x.valueChanged.connect(self.update_stimulator)
        Layout.addWidget(self.ref_x,2,1)

        self.label_y2 = QtWidgets.QLabel("Y:")
        Layout.addWidget(self.label_y2,1,2)
        self.ref_y = QtWidgets.QDoubleSpinBox()
        self.ref_y.setRange(-1000, 1000)
        if load_ref: self.ref_y.setValue(self.reference[1])
        #self.ref_y.valueChanged.connect(self.update_stimulator)
        Layout.addWidget(self.ref_y,2,2)

        self.label_z2 = QtWidgets.QLabel("Z:")
        Layout.addWidget(self.label_z2,1,3)
        self.ref_z = QtWidgets.QDoubleSpinBox()
        self.ref_z.setRange(-1000, 1000)
        if load_ref: self.ref_z.setValue(self.reference[2])
        #self.ref_z.valueChanged.connect(self.update_stimulator)
        Layout.addWidget(self.ref_z,2,3)

        self.direction_box.setLayout(Layout)

    #Returns the positions
    def getPositions(self):
        self.centre = [0,0,0]
        self.centre[0] = self.pos_x.value()
        self.centre[1] = self.pos_y.value()
        self.centre[2] = self.pos_z.value()
        #if self.mode == 'TMS' or (self.mode == 'TDCS' and self.Use_ref.isChecked()): 
        self.reference = [0,0,0]
        self.reference[0] = self.ref_x.value()
        self.reference[1] = self.ref_y.value()
        self.reference[2] = self.ref_z.value()

        return self.centre, self.reference 

    #Gets the positions from clicking in the head model
    def getClickedPosition(self):
        #Get ref => clicling the model will change the numbers in the reference box
        get_ref = False
        if self.check.isChecked():
            get_ref = True

        Position = self.glHeadModel.getIntersectPoint()
        if Position is not None and not get_ref:
            self.pos_x.setValue(Position[0])
            self.pos_y.setValue(Position[1])
            self.pos_z.setValue(Position[2])
            #if self.mode == 'TMS':
            self.name = ''
            self.check.toggle()
            self.update_center()
        if Position is not None and get_ref:
            self.ref_x.setValue(Position[0])
            self.ref_y.setValue(Position[1])
            self.ref_z.setValue(Position[2])
            #if self.mode == 'TMS':
            self.check.toggle()
            self.update_stimulator()


    def project_points(self):
        reference = np.array((self.ref_x.value(),self.ref_y.value(),self.ref_z.value()))
        center = np.array((self.pos_x.value(),self.pos_y.value(),self.pos_z.value()))
        projected_center,_ = self.glHeadModel.getSurface('Scalp').projectPoint(center)
        projected_reference = projected_center + (reference-center)
        self.ref_x.setValue(projected_reference[0])
        self.ref_y.setValue(projected_reference[1])
        self.ref_z.setValue(projected_reference[2])

        self.pos_x.setValue(projected_center[0])
        self.pos_y.setValue(projected_center[1])
        self.pos_z.setValue(projected_center[2])
        self.update_stimulator()

    #Executes the GUI
    def update_center(self):
        centre, _ = self.getPositions()
        pos_y = _get_posy(centre, self.glHeadModel.getSurface('Scalp'))
        self.ref_x.setValue(pos_y[0])
        self.ref_y.setValue(pos_y[1])
        self.ref_z.setValue(pos_y[2])

        self.update_stimulator()

    def update_stimulator(self):
        centre, reference = self.getPositions()
        scalp_surf = self.glHeadModel.getSurface('Scalp')
        transf_matrix = scalp_surf.calculateMatSimnibs(centre, reference)
        ogl_object = self.glHeadModel.drawPointAndDirs(
            transf_matrix, QtGui.QColor.fromCmykF(0., 0., 0., 0.74))
        self.glHeadModel.tmpObjectList([ogl_object])
        self.glHeadModel.update()

    def GetPositions(self):
        result = self.exec_()
        self.glHeadModel.clearTmpObjects()
        centre, reference = self.getPositions()
        return (centre, reference, self.name, result == QtWidgets.QDialog.Accepted)


def _get_posy(centre, surf):
    _, normal = surf.projectPoint(centre)
    pos_y = [centre[0], centre[1], centre[2]]
    if np.abs(normal[1]) < .8:
        pos_y[1] += 10
    else:
        pos_y[0] += 10
    return pos_y

#GUI for definin conductivities
class ConductivitiesGui (QtWidgets.QDialog):
    def __init__(self, simulist):
        super(ConductivitiesGui, self).__init__()

        self.simulist = simulist

        mainLayout = QtWidgets.QVBoxLayout()

        self.setCondTable()
        self.tensorBox = self.setTensorThings()
        self.changeCondType()

        mainLayout.addWidget(self.condTable)
        mainLayout.addWidget(self.tensorBox)

        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.button_box.accepted.connect(self.checkAndAccept)
        #self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)
        self.resetButton = QtWidgets.QPushButton("Reset")
        self.resetButton.clicked.connect(self.reset)
        self.button_box.addButton(self.resetButton, QtWidgets.QDialogButtonBox.ResetRole)

        mainLayout.addWidget(self.button_box)

        self.setLayout(mainLayout)    
        self.setWindowTitle("Conductivities")
        self.resize(450,400)

    def setCondTable(self):
        self.condTable = QtWidgets.QTableWidget(0,2)
        self.condTable.setHorizontalHeaderLabels(("Name", "Value")) 
        self.condTable.setShowGrid(False)
        self.condTable.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        self.condTable.verticalHeader().hide()


        self.setRowValues()

    #Sets up each row
    def setRowValues(self):
        self.conds_dict = {}
        self.cond_values = []
        self.cond_names = []
        i = 0
        j = 0
        for i,cnd in enumerate(self.simulist.cond):
            if cnd.name != '' and cnd.value is not None:
                self.condTable.insertRow(j)
                self.conds_dict[j] = i
                self.condTable.setCellWidget(j, 0, QtWidgets.QLabel(cnd.name))
                ##Setting up spin box for conduictivity
                self.cond_values.append(QtWidgets.QDoubleSpinBox())
                self.cond_values[-1].setSuffix(" S/m")
                self.cond_values[-1].setSingleStep(0.01)
                self.cond_values[-1].setDecimals(3)
                self.cond_values[-1].setValue(cnd.value)
                self.cond_values[-1].setMinimum(1e-5)
                self.condTable.setCellWidget(j, 1, self.cond_values[-1])
                j += 1



    def setTensorThings(self):
        tensorBox = QtWidgets.QGroupBox('Brain Anisotropy')
        layout = QtWidgets.QVBoxLayout()

        cond_type_CBox = QtWidgets.QComboBox()
        cond_type_CBox.addItems(['scalar','volume normalized','direct mapping','mean conductivity'])
        cond_type_CBox.currentIndexChanged.connect(self.changeCondType)
        if self.simulist.anisotropy_type == 'vn':
            cond_type_CBox.setCurrentIndex(1)
        elif self.simulist.anisotropy_type == 'dir':
            cond_type_CBox.setCurrentIndex(2)
        elif self.simulist.anisotropy_type == 'mc':
            cond_type_CBox.setCurrentIndex(3)
        #layout.addWidget(cond_type_CBox, 0, QtCore.Qt.Alignment(1))
        layout.addWidget(cond_type_CBox, 0)

        layout.addWidget(QtWidgets.QLabel('Maximum ratio between eigenvalues:'), 1)

        aniso_maxratio_sbox = QtWidgets.QDoubleSpinBox()
        aniso_maxratio_sbox.setSingleStep(0.1)
        aniso_maxratio_sbox.setDecimals(2)
        aniso_maxratio_sbox.setMinimum(1)
        aniso_maxratio_sbox.setValue(self.simulist.aniso_maxratio)
        aniso_maxratio_sbox.valueChanged.connect(self.changeCondType)

        layout.addWidget(aniso_maxratio_sbox, 0)

        layout.addWidget(QtWidgets.QLabel('Maximum eigenvalue:'), 1)

        aniso_maxcond_sbox = QtWidgets.QDoubleSpinBox()
        aniso_maxcond_sbox.setSingleStep(0.1)
        aniso_maxcond_sbox.setDecimals(2)
        aniso_maxcond_sbox.setMinimum(1e-2)
        aniso_maxcond_sbox.setValue(self.simulist.aniso_maxcond)
        aniso_maxcond_sbox.valueChanged.connect(self.changeCondType)

        layout.addWidget(aniso_maxcond_sbox, 0)
 
        tensorBox.setLayout(layout)
        tensorBox.cond_type_CBox = cond_type_CBox
        tensorBox.aniso_maxratio_sbox = aniso_maxratio_sbox
        tensorBox.aniso_maxcond_sbox = aniso_maxcond_sbox

        return tensorBox


    def changeCondType(self):
        if self.tensorBox.cond_type_CBox.currentIndex() == 0:
            self.simulist.anisotropy_type = 'scalar'
            self.tensorBox.aniso_maxratio_sbox.setEnabled(False)
            self.tensorBox.aniso_maxcond_sbox.setEnabled(False)

        elif self.tensorBox.cond_type_CBox.currentIndex() == 1:
            self.simulist.anisotropy_type = 'vn'
            self.tensorBox.aniso_maxratio_sbox.setEnabled(True)
            self.tensorBox.aniso_maxcond_sbox.setEnabled(True)

        elif self.tensorBox.cond_type_CBox.currentIndex() == 2:
            self.simulist.anisotropy_type = 'dir'
            self.tensorBox.aniso_maxratio_sbox.setEnabled(True)
            self.tensorBox.aniso_maxcond_sbox.setEnabled(True)

        elif self.tensorBox.cond_type_CBox.currentIndex() == 3:
            self.simulist.anisotropy_type = 'mc'
            self.tensorBox.aniso_maxratio_sbox.setEnabled(False)
            self.tensorBox.aniso_maxcond_sbox.setEnabled(True)

        self.simulist.aniso_maxratio = self.tensorBox.aniso_maxratio_sbox.value()
        self.simulist.aniso_maxcond = self.tensorBox.aniso_maxcond_sbox.value()

    #Resets everything to the standard conductivities
    def reset(self):
        standard = standard_cond()
        for j in range(self.condTable.rowCount()):
            if self.conds_dict[j] > 99 and self.conds_dict[j] < 499:
                self.cond_values[j].setValue(standard[99].value)
            elif self.conds_dict[j] > 499:
                self.cond_values[j].setValue(standard[499].value)
            else:
                self.cond_values[j].setValue(standard[self.conds_dict[j]].value)

        self.simulist.aniso_maxratio = sim_struct.SimuList().aniso_maxratio
        self.simulist.aniso_maxcond = sim_struct.SimuList().aniso_maxcond
        self.simulist.anisotropy_type = sim_struct.SimuList().anisotropy_type
        self.tensorBox.cond_type_CBox.setCurrentIndex(0)
        self.tensorBox.aniso_maxratio_sbox.setValue(self.simulist.aniso_maxratio)
        self.tensorBox.aniso_maxcond_sbox.setValue(self.simulist.aniso_maxcond)


    #Returns the conductivities
    def checkAndAccept(self):
        for j in range(self.condTable.rowCount()):
            self.simulist.cond[self.conds_dict[j]].value = self.cond_values[j].value()

        self.accept()


    def getCond(self):
        result = self.exec_()
        #cond = self.returnCond()
        return (self.simulist, result == QtWidgets.QDialog.Accepted)


#Thread for loading meshes
class loadMeshThread(QtCore.QThread):

    def __init__(self, mesh_fn, glHeadModel):
        QtCore.QThread.__init__(self)
        self.fn = mesh_fn
        self.glHeadModel = glHeadModel

    def run(self):
        self.glHeadModel.loadMesh(self.fn)
        self.exit(0)
        self.exec_()

#Thread for runing simulation
@QtCore.pyqtSlot(str)
class runSimThread(QtCore.QThread):
    output_signal = QtCore.pyqtSignal(str)

    def __init__(self, session):
        QtCore.QThread.__init__(self)
        self.session = session
        self._stop = False

    def run(self):
        class WriteToBoxHandler(logging.StreamHandler):
            def __init__(self, out_signal, stop_signal):
                super().__init__()
                self.out_signal = out_signal
                self.stop_signal = stop_signal

            def emit(self, record):
                msg = self.format(record)
                if self.stop_signal():
                    raise Exception('Execution Stopped')
                self.out_signal.emit(msg)

        w2b_handler = WriteToBoxHandler(
            self.output_signal, self.get_stop)
        w2b_handler.setFormatter(
            logging.Formatter('%(levelname)s: %(message)s'))
        w2b_handler.setLevel(logging.INFO)
        logger.addHandler(w2b_handler)
        run_simnibs(self.session)
        logger.removeHandler(w2b_handler)
        self.exit(0)

    def set_stop(self):
        self._stop = True
        self.exit(1)

    def get_stop(self):
        return self._stop

class openGmshThread(QtCore.QThread):
    def __init__(self, fn):
        QtCore.QThread.__init__(self)
        self.fn = fn

    def run(self):
        gmsh = path2bin('gmsh')
        gmsh_return = subprocess.run([gmsh, self.fn])
        self.exit(gmsh_return.returncode)
        self.exec_()

def except_hook(cls, exception, traceback):
    QtWidgets.QMessageBox.critical(None, 'SimNIBS GUI Error', str(exception))
    sys.__excepthook__(cls, exception, traceback)
    #sys.exit(app.exec_())

def start_gui(args):
    app = QtWidgets.QApplication(args)
    #app.setAttribute(QtCore.Qt.AA_UseDesktopOpenGL)
    #app.setAttribute(QtCore.Qt.AA_UseSoftwareOpenGL)
    #app.setAttribute(QtCore.Qt.AA_UseOpenGLES)
    sys.excepthook = except_hook
    ex = TDCS_GUI()
    ex.show()
    if len(args) > 1:
        if args[1].endswith(".mat"):
            ex.openSimnibsFile(args[1])
        elif args[1].endswith(".msh"):
            ex.loadHeadModel(args[1])
        else:
            raise IOError('simnibs_gui can only load .mat and .msh files')
    sys.exit(app.exec_())


        
