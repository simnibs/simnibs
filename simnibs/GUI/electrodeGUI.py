#!/usr/bin/python2.7
# -*- coding: utf-8 -*-\
'''

    Window for defining electrode proprieties
    Part of the main SimNIBS GUI 
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.
    
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


from PyQt5 import QtCore, QtGui, QtOpenGL, QtWidgets
from OpenGL import GL, GLU
import math
import copy
import sys
import numpy

from ..simulation import sim_struct


class Ui_Electrode(QtWidgets.QDialog):

    def __init__(self, electrode_struct):
        super(Ui_Electrode, self).__init__()

        self.electrode_struct = copy.deepcopy(electrode_struct)
        mainLayout = QtWidgets.QGridLayout()

        self.glElectrode = GLElectrode()
        self.set_shape_size_layout()
        self.set_connector_layout()
        self.set_type_tick_layout()
        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok|QtWidgets.QDialogButtonBox.Cancel)
        self.button_box.accepted.connect(self.checkAndAccept)
        #self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)

        mainLayout.addWidget(self.shape_size, 0,4)
        mainLayout.addWidget(self.type_thick, 1,4)
        mainLayout.addWidget(self.connector_box,2,4)
        mainLayout.addWidget(self.glElectrode,0,0,4,4)
        mainLayout.addWidget(self.button_box, 4,4)


        self.update_electrode()
        self.glElectrode.electrode_object, self.glElectrode.plug_object = self.glElectrode.drawElectrode(self.electrode_struct)
        self.glElectrode.electrode_struct = self.electrode_struct
        self.glElectrode.updateGL()

        self.setLayout(mainLayout)
        self.resize(900,450)

        self.setWindowTitle('Electrode')


    #Shape and size box
    def set_shape_size_layout(self):

        self.shape_size = QtWidgets.QGroupBox("Shape and Size")
        layout = QtWidgets.QGridLayout()

        self.shape_selector = QtWidgets.QComboBox()
        self.shape_selector.addItem("Rectangular")
        self.shape_selector.addItem("Elliptical")
        if self.electrode_struct.shape == 'ellipse':
            self.shape_selector.setCurrentIndex(1)
        self.shape_selector.activated.connect(self.update_electrode)

        layout.addWidget(self.shape_selector,0,0)

        self.dim_x = QtWidgets.QDoubleSpinBox()
        if self.electrode_struct.dimensions is not None:
             self.dim_x.setValue(self.electrode_struct.dimensions[0]/10.0)
        else:
            self.dim_x.setValue(5.0)
        self.dim_x.setSuffix(" cm")
        self.dim_x.setSingleStep(0.5)
        self.dim_x.setMinimum(0.5)

        layout.addWidget(self.dim_x, 0, 1)
        self.dim_x.valueChanged.connect(self.update_electrode)


        self.dim_y = QtWidgets.QDoubleSpinBox()
        if self.electrode_struct.dimensions is not None:
            self.dim_y.setValue(self.electrode_struct.dimensions[1]/10)
        else:
            self.dim_y.setValue(5.0)
        self.dim_y.setSuffix(" cm")
        self.dim_y.setSingleStep(0.5)
        self.dim_y.setMinimum(0.5)

        layout.addWidget(self.dim_y, 0, 2)
        self.dim_y.valueChanged.connect(self.update_electrode)


        self.shape_size.setLayout(layout)

    #Type and thickness box
    def set_type_tick_layout(self):

        self.type_thick = QtWidgets.QGroupBox("Type and Thickness")
        layout = QtWidgets.QGridLayout()

        self.type_selector = QtWidgets.QComboBox()
        self.type_selector.addItem("Simple")
        self.type_selector.addItem("Electrode+Gel")
        self.type_selector.addItem("Electrode+Sponge")
        if len(self.electrode_struct.thickness) == 2:
            self.type_selector.setCurrentIndex(1)
        if len(self.electrode_struct.thickness) == 3:
            self.type_selector.setCurrentIndex(2)
        self.type_selector.activated.connect(self.update_labels)
        self.type_selector.activated.connect(self.update_electrode)
        layout.addWidget(self.type_selector,0,0,1,2)


        self.label_thick1= QtWidgets.QLabel("Electrode Thickness:")
        layout.addWidget(self.label_thick1, 1, 0,1,1)
        self.thick1 = QtWidgets.QDoubleSpinBox()
        if len(self.electrode_struct.thickness) == 1:
            self.thick1.setValue(self.electrode_struct.thickness[0])
        if len(self.electrode_struct.thickness) == 2 or len(self.electrode_struct.thickness) == 3:
            self.thick1.setValue(self.electrode_struct.thickness[1])
        self.thick1.valueChanged.connect(self.update_electrode)
        self.thick1.setSuffix(" mm")
        self.thick1.setSingleStep(0.5)
        layout.addWidget(self.thick1,1,1,1,1)

        self.label_thick2= QtWidgets.QLabel("Gel/Sponge Thickness:")
        layout.addWidget(self.label_thick2, 2, 0,1,1)
        self.thick2 = QtWidgets.QDoubleSpinBox()
        if len(self.electrode_struct.thickness) == 2:
            self.thick2.setValue(self.electrode_struct.thickness[0])
        if len(self.electrode_struct.thickness) == 3:
            self.thick2.setValue(2*self.electrode_struct.thickness[0])
        if len(self.electrode_struct.thickness) == 0:
            self.thick1.setValue(5.0)
            self.thick2.setValue(5.0)
        self.thick2.valueChanged.connect(self.update_electrode)
        self.thick2.setSuffix(" mm")
        self.thick2.setSingleStep(0.5)
        layout.addWidget(self.thick2,2,1,1,1)

        self.label_sporatio= QtWidgets.QLabel("Electrode/Sponge size:")
        layout.addWidget(self.label_sporatio, 3, 0,1,1)

        self.spo_x = QtWidgets.QDoubleSpinBox()
        self.spo_x.setMinimum(1.0 * self.dim_x.value())
        self.spo_x.setSingleStep(0.5)
        layout.addWidget(self.spo_x,3,1,1,1)


        self.spo_y = QtWidgets.QDoubleSpinBox()
        self.spo_y.setMinimum(1.0 * self.dim_y.value())
        self.spo_y.setSingleStep(0.5)
        layout.addWidget(self.spo_y,3,2,1,1)

        if len(self.electrode_struct.thickness) == 3 and\
           self.electrode_struct.dimensions_sponge is not None:
            self.spo_x.setValue(self.electrode_struct.dimensions_sponge[0] / 10.)
            self.spo_y.setValue(self.electrode_struct.dimensions_sponge[1] / 10.)

        self.spo_y.valueChanged.connect(self.update_electrode)
        self.spo_x.valueChanged.connect(self.update_electrode)
        self.update_labels()
        self.type_thick.setLayout(layout)

    #Box for defining connector proprieties
    def set_connector_layout(self):

        self.connector_box = QtWidgets.QGroupBox("Connector information")

        layout1 = QtWidgets.QGridLayout()

        has_connector = False
        if len(self.electrode_struct.plug) == 1:
            has_connector = True

        self.No_Connector = QtWidgets.QRadioButton("Whole Surface")
        if not(has_connector):
            self.No_Connector.toggle()
        self.No_Connector.toggled.connect(self.enable_disable_connector)
        self.No_Connector.toggled.connect(self.update_electrode)
        layout1.addWidget(self.No_Connector, 0, 0,1,1)


        self.Yes_Connector = QtWidgets.QRadioButton("Define Connector")
        if (has_connector):
            self.Yes_Connector.toggle()
        self.Yes_Connector.toggled.connect(self.enable_disable_connector)
        self.Yes_Connector.toggled.connect(self.update_electrode)
        layout1.addWidget(self.Yes_Connector, 0, 1,1,1)

        self.label_con_pos = QtWidgets.QLabel("Position in relation to electrode center:")
        layout1.addWidget(self.label_con_pos, 1, 0,1,3)

        self.con_position_x = QtWidgets.QDoubleSpinBox()
        self.con_position_x.setMinimum(-20.0)
        if has_connector:
            self.con_position_x.setValue(self.electrode_struct.plug[0].centre[0]/10.0)
        self.con_position_x.valueChanged.connect(self.update_electrode)
        self.con_position_x.setSingleStep(0.5)
        layout1.addWidget(self.con_position_x,2,0)

        self.con_position_y = QtWidgets.QDoubleSpinBox()
        self.con_position_y.setMinimum(-20.0)
        if has_connector:
            self.con_position_y.setValue(self.electrode_struct.plug[0].centre[1]/10.0)
        self.con_position_y.valueChanged.connect(self.update_electrode)
        self.con_position_y.setSingleStep(0.5)
        layout1.addWidget(self.con_position_y,2,1)

        self.label_con_shape = QtWidgets.QLabel("Shape:")
        layout1.addWidget(self.label_con_shape, 3, 0)

        self.con_shape_selector = QtWidgets.QComboBox()
        self.con_shape_selector.addItem("Rectangular")
        self.con_shape_selector.addItem("Elliptical")
        if has_connector:
            if self.electrode_struct.plug[0].shape == 'ellipse':
                self.con_shape_selector.setCurrentIndex(1)
        self.con_shape_selector.activated.connect(self.update_electrode)
        layout1.addWidget(self.con_shape_selector,4,0)

        self.con_dim_x = QtWidgets.QDoubleSpinBox()
        if has_connector:
            self.con_dim_x.setValue(self.electrode_struct.plug[0].dimensions[0]/10.0)
        else:
            self.con_dim_x.setValue(2.0)
        self.con_dim_x.valueChanged.connect(self.update_electrode)
        self.con_dim_x.setSuffix(" cm")
        self.con_dim_x.setSingleStep(0.5)
        self.con_dim_x.setMinimum(0.5)

        layout1.addWidget(self.con_dim_x, 4, 1)

        self.con_dim_y = QtWidgets.QDoubleSpinBox()
        if has_connector:
            self.con_dim_y.setValue(self.electrode_struct.plug[0].dimensions[1]/10.0)
        else:
            self.con_dim_y.setValue(1.0)
        self.con_dim_y.valueChanged.connect(self.update_electrode)
        self.con_dim_y.setSuffix(" cm")
        self.con_dim_y.setSingleStep(0.5)
        self.con_dim_y.setMinimum(0.5)

        layout1.addWidget(self.con_dim_y, 4, 2)


        self.enable_disable_connector()

        self.connector_box.setLayout(layout1)


    #Updates electrode structure and OpenGL model
    def update_electrode(self) :

        if self.shape_selector.currentIndex() == 0:
            self.electrode_struct.definition = 'plane'
            self.electrode_struct.shape = 'rect'
        if self.shape_selector.currentIndex() == 1:
            self.electrode_struct.definition = 'plane'
            self.electrode_struct.shape = 'ellipse'


        #write dimensions
        self.electrode_struct.dimensions = [self.dim_x.value()*10,
                                            self.dim_y.value()*10]

        #write thickness, for different electrode types
            #electrode in sponge
        if self.type_selector.currentIndex() == 2:
            self.electrode_struct.thickness = [0.0,0.0,0.0]
            self.electrode_struct.thickness[0] = self.thick2.value()/2.0
            self.electrode_struct.thickness[1] = self.thick1.value()
            self.electrode_struct.thickness[2] = self.thick2.value()/2.0
            self.electrode_struct.dimensions_sponge = [self.spo_x.value() * 10,
                                                       self.spo_y.value() * 10]
            self.spo_x.setMinimum(1.0 * self.dim_x.value())
            self.spo_y.setMinimum(1.0 * self.dim_y.value())
            self.dim_x.setMaximum(1./1.0 * self.spo_x.value())
            self.dim_y.setMaximum(1./1.0 * self.spo_y.value())


            #electrode+gel
        if self.type_selector.currentIndex() == 1:
            self.electrode_struct.thickness = [0.0,0.0]
            self.electrode_struct.thickness[0] = self.thick2.value()
            self.electrode_struct.thickness[1] = self.thick1.value()
            self.electrode_struct.dimensions_sponge = None
            #simple
        if self.type_selector.currentIndex() == 0:
            self.electrode_struct.thickness = [0.0]
            self.electrode_struct.thickness[0] = self.thick1.value()
            self.electrode_struct.dimensions_sponge = None

        ##connector
        if self.Yes_Connector.isChecked():
            if len(self.electrode_struct.plug) == 0 :
                self.electrode_struct.plug.append(sim_struct.ELECTRODE())

            if self.con_shape_selector.currentIndex() == 0:
                self.electrode_struct.plug[0].definition = 'plane'
                self.electrode_struct.plug[0].shape = 'rect'
            if self.con_shape_selector.currentIndex() == 1:
                self.electrode_struct.plug[0].definition = 'plane'
                self.electrode_struct.plug[0].shape = 'ellipse'

            self.electrode_struct.plug[0].dimensions = [self.con_dim_x.value()*10,
                                                        self.con_dim_y.value()*10]
            self.electrode_struct.plug[0].centre = [0.0,0.0]
            self.electrode_struct.plug[0].centre[0]= self.con_position_x.value()*10
            self.electrode_struct.plug[0].centre[1]= self.con_position_y.value()*10


        if self.No_Connector.isChecked():
            self.electrode_struct.plug = []



        self.glElectrode.electrode_object, self.glElectrode.plug_object= self.glElectrode.drawElectrode(self.electrode_struct)
        self.glElectrode.electrode_struct = self.electrode_struct
        self.glElectrode.updateGL()
        QtWidgets.QApplication.processEvents()


    def update_labels(self):
        if self.type_selector.currentIndex() == 0:
            self.thick2.setEnabled(False)
            self.label_thick2.setEnabled(False)
            self.spo_x.setEnabled(False)
            self.spo_y.setEnabled(False)
            self.label_sporatio.setEnabled(False)

            if self.electrode_struct.thickness is not None and\
                len(self.electrode_struct.thickness) == 1:
                self.thick1.setValue(self.electrode_struct.thickness[0])
            else:
                self.thick1.setValue(5.00)
            self.thick2.setValue(0.00)
            self.spo_x.setValue(0)
            self.spo_y.setValue(0)


        if self.type_selector.currentIndex() == 1:
            self.thick2.setEnabled(True)
            self.label_thick2.setEnabled(True)
            self.spo_x.setEnabled(False)
            self.spo_y.setEnabled(False)
            self.label_sporatio.setEnabled(False)

            if self.electrode_struct.thickness is not None and\
                len(self.electrode_struct.thickness) == 2:
                self.thick2.setValue(self.electrode_struct.thickness[0])
                self.thick1.setValue(self.electrode_struct.thickness[1])
            else:
                self.thick1.setValue(1.00)
                self.thick2.setValue(5.00)

            self.spo_x.setValue(0.0)
            self.spo_y.setValue(0.0)


        if self.type_selector.currentIndex() == 2:
            self.thick2.setEnabled(True)
            self.label_thick2.setEnabled(True)
            self.spo_x.setEnabled(True)
            self.spo_y.setEnabled(True)
            self.label_sporatio.setEnabled(True)

            if self.electrode_struct.thickness is not None and\
                len(self.electrode_struct.thickness) == 3:
                self.thick2.setValue(self.electrode_struct.thickness[0] * 2)
                self.thick1.setValue(self.electrode_struct.thickness[1])
            else:
                self.thick1.setValue(1.00)
                self.thick2.setValue(8.00)

            self.spo_x.setMinimum(1.0 * self.dim_x.value())
            self.spo_y.setMinimum(1.0 * self.dim_y.value())
            self.dim_x.setMaximum(1./1.0 * self.spo_x.value())
            self.dim_y.setMaximum(1./1.0 * self.spo_y.value())



    #Disable/enable things envolving the connector
    def enable_disable_connector(self):
        has_con = self.Yes_Connector.isChecked()

        self.con_position_x.setEnabled(has_con)
        self.con_position_y.setEnabled(has_con)
        self.con_dim_x.setEnabled(has_con)
        self.con_dim_y.setEnabled(has_con)
        self.con_shape_selector.setEnabled(has_con)
        self.label_con_pos.setEnabled(has_con)
        self.label_con_shape.setEnabled(has_con)


    #checks if connector is completelly inside the electrode
    def checkAndAccept(self):
        '''
        if len(self.electrode_struct.plug) == 1:
            if abs(self.electrode_struct.plug[0].centre[0])+self.electrode_struct.plug[0].dimX/2.0 > self.electrode_struct.dimX/2.0:
                reply = QtWidgets.QMessageBox.critical(self, "Warning",
                "Plug outside electrode in X dimension")
                return None
            elif abs(self.electrode_struct.plug[0].centre[1])+self.electrode_struct.plug[0].dimY/2.0 > self.electrode_struct.dimY/2.0:
                reply = QtWidgets.QMessageBox.critical(self, "Warning",
                "Plug outside electrode in Y dimension")
                return None
        '''
        self.accept()




    #Returns electrode structure
    def return_el_struct(self):
        dialog = Ui_Electrode(self.electrode_struct)
        result = dialog.exec_()
        struct = dialog.electrode_struct
        return (struct, result == QtWidgets.QDialog.Accepted)






#####################################################################
############Electrode OpenGL Model###################################
#####################################################################
class GLElectrode(QtOpenGL.QGLWidget):
    def __init__(self, parent=None):
        QtOpenGL.QGLWidget.__init__(self, parent)

        self.electrode_object = 0
        self.plug_object = 0
        self.center_dir = 0

        self.xRot = 3160
        self.yRot = 5536
        self.zRot = 0
        self.zoom = 1.0

        self.electrode_placeholder = sim_struct.ELECTRODE()
        self.electrode_placeholder.dimensions = [0.0, 0.0]
        self.electrode_placeholder.thickness = [0.0]

        self.lastPos = QtCore.QPoint()

    def minimumSizeHint(self):
        return QtCore.QSize(50, 50)

    def sizeHint(self):
        return QtCore.QSize(450, 450)

    def setXRotation(self, angle):
        angle = self.normalizeAngle(angle)
        if angle != self.xRot:
            self.xRot = angle
            #self.emit(QtCore.SIGNAL("xRotationChanged(int)"), angle)
            self.updateGL()

    def setYRotation(self, angle):
        angle = self.normalizeAngle(angle)
        if angle != self.yRot:
            self.yRot = angle
            #self.emit(QtCore.SIGNAL("yRotationChanged(int)"), angle)
            self.updateGL()

    def setZRotation(self, angle):
        angle = self.normalizeAngle(angle)
        if angle != self.zRot:
            self.zRot = angle
            #self.emit(QtCore.SIGNAL("zRotationChanged(int)"), angle)
            self.updateGL()

    def initializeGL(self):
        self.qglClearColor(QtGui.QColor.fromCmykF(0.1, .0, 0.2, 1.0))
        self.electrode_object, self.plug_object = \
            self.drawElectrode(self.electrode_struct)
        self.center_dir = self.drawPointAndDirs()
        GL.glShadeModel(GL.GL_FLAT)
        GL.glEnable(GL.GL_DEPTH_TEST)
        
    def paintGL(self):
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glLoadIdentity()
        GL.glTranslated(0.0, 0.0, -10.0)
        GL.glRotated(self.xRot / 16.0, 1.0, 0.0, 0.0)
        GL.glRotated(self.yRot / 16.0, 0.0, 1.0, 0.0)
        GL.glRotated(self.zRot / 16.0, 0.0, 0.0, 1.0)
        GL.glScalef(self.zoom,self.zoom,self.zoom)

        GL.glCallList(self.electrode_object)
        GL.glCallList(self.center_dir)
        if self.plug_object != 0:
            GL.glCallList(self.plug_object)

    def resizeGL(self, width, height):
        side = min(width, height)
        GL.glViewport((width - side) // 2, (height - side) // 2, side, side)
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        #GL.glOrtho(-0.5, +0.5, +0.5, -0.5, 4.0, 15.0)
        GL.glFrustum(-0.5, +0.5, +0.5, -0.5, 4.0, 15.0)
        GL.glMatrixMode(GL.GL_MODELVIEW)

    def mousePressEvent(self, event):
        self.lastPos = QtCore.QPoint(event.pos())

    def mouseMoveEvent(self, event):
        dx = event.x() - self.lastPos.x()
        dy = event.y() - self.lastPos.y()

        if event.buttons() & QtCore.Qt.LeftButton:
            self.setXRotation(self.xRot - 8 * dy)
            self.setYRotation(self.yRot - 8 * dx)
        elif event.buttons() & QtCore.Qt.RightButton:
            self.setXRotation(self.xRot - 8 * dy)
            self.setZRotation(self.zRot - 8 * dx)

        self.lastPos = QtCore.QPoint(event.pos())


    def wheelEvent(self, event):
        delta = event.angleDelta().y()
        zoom = self.zoom+delta/1200.0
        if zoom > 0.1 and zoom < 10:
            self.zoom = zoom
            self.updateGL()


    def drawElectrode(self, electrode_struct):
        yellow = QtGui.QColor.fromCmykF(0.00, 0.0, 1.0, 0.0)
        grey = QtGui.QColor.fromCmykF(0.00, 0.0, 0.0, 0.7)
        blue = QtGui.QColor.fromCmykF(0.83, 0.54, 0.0, 0.0)
        red = QtGui.QColor.fromCmykF(0.00, 0.99, 1.0, 0.0)
        scale = 0.02

        try:
            size_x = electrode_struct.dimensions[0] * scale
        except:
            size_x = 0
        try:
            size_y = electrode_struct.dimensions[1] * scale
        except: 
            size_y = 0
            
        thick1 = 0
        thick2 = 0
        thick_sponge = 0
        z_center = 0
        model = ''
        if len(electrode_struct.thickness) == 1:
            thick1 = electrode_struct.thickness[0]*scale
            model = 'simple'
        if len(electrode_struct.thickness) == 2:
            thick1 = electrode_struct.thickness[1]*scale
            thick2 = electrode_struct.thickness[0]*scale
            model = '2layer'
        if len(electrode_struct.thickness) == 3:
            thick1 = electrode_struct.thickness[1]*scale
            thick2 = electrode_struct.thickness[0]*scale + \
                electrode_struct.thickness[1]*scale + electrode_struct.thickness[2]*scale
            model = 'sponge'

        electrode_list = GL.glGenLists(1)
        GL.glNewList(electrode_list, GL.GL_COMPILE)

        if electrode_struct.shape == 'rect':
            if model == 'simple' :
                self.cuboid(0,0,0,size_x,size_y,thick1,blue, False)
            if model == '2layer' and thick2 > 0:
                self.cuboid(0,0,0,size_x,size_y,thick1,grey, False)
                self.cuboid(0,0,(thick1+thick2)/2,size_x,size_y,thick2,blue, False)
            if model == 'sponge':
                self.cuboid(0,0,0,size_x,size_y,thick1,grey, False)
                sponge_size_x = electrode_struct.dimensions_sponge[0] * scale
                sponge_size_y = electrode_struct.dimensions_sponge[1] * scale
                self.cuboid(0,0,0,sponge_size_x,sponge_size_y,thick2,yellow, True)

        if electrode_struct.shape == 'ellipse':
            if model == 'simple' :
                self.cylinder(0,0,0,size_x,size_y,thick1,blue, False)
            if model == '2layer' and thick2 > 0:
                self.cylinder(0,0,0,size_x,size_y,thick1,grey, False)
                self.cylinder(0,0,(thick1+thick2)/2,size_x,size_y,thick2,blue, False)
            if model == 'sponge':
                self.cylinder(0,0,0,size_x,size_y,thick1,grey, False)
                sponge_size_x = electrode_struct.dimensions_sponge[0] * scale
                sponge_size_y = electrode_struct.dimensions_sponge[1] * scale
                self.cylinder(0,0,0,sponge_size_x,sponge_size_y,thick2,yellow, True)

        GL.glEndList()

        plug_list = GL.glGenLists(1)
        GL.glNewList(plug_list, GL.GL_COMPILE)

        ##Model Plug
        if len(electrode_struct.plug) > 0 :
            con_size_x = electrode_struct.plug[0].dimensions[0]*scale
            con_size_y = electrode_struct.plug[0].dimensions[1]*scale
            con_center_x = electrode_struct.plug[0].centre[0]*scale
            con_center_y = electrode_struct.plug[0].centre[1]*scale
            dz = 0.001
            if electrode_struct.plug[0].shape == 'rect':
                self.plane_xy(con_center_x, con_center_y, -(thick1/2+dz), con_size_x, con_size_y, red)
            if electrode_struct.plug[0].shape == 'ellipse':
                self.ellipse_xy(con_center_x, con_center_y, -(thick1/2+dz), con_size_x, con_size_y, red)


        GL.glEndList()

        return electrode_list, plug_list



    def cuboid(self, center_x, center_y,center_z, size_x, size_y, size_z, color, wireframe):

        if wireframe :
            GL.glPolygonMode( GL.GL_FRONT_AND_BACK, GL.GL_LINE ) 

        GL.glBegin(GL.GL_QUADS)
        self.qglColor(color)

        px = center_x+size_x/2
        mx =center_x-size_x/2
        py = center_y+size_y/2
        my =center_y-size_y/2
        pz = center_z+size_z/2
        mz =center_z-size_z/2


        GL.glVertex3d(px,py, pz)
        GL.glVertex3d(px,my, pz)
        GL.glVertex3d(mx,my, pz)
        GL.glVertex3d(mx,py, pz)

        GL.glVertex3d(px,py, mz)
        GL.glVertex3d(mx,py, mz)
        GL.glVertex3d(mx,my, mz)
        GL.glVertex3d(px,my, mz)

        self.qglColor(color.darker(250))

        GL.glVertex3d(px,py, pz)
        GL.glVertex3d(mx,py, pz)
        GL.glVertex3d(mx,py, mz)
        GL.glVertex3d(px,py, mz)

        GL.glVertex3d(px,my, pz)
        GL.glVertex3d(mx,my, pz)
        GL.glVertex3d(mx,my, mz)
        GL.glVertex3d(px,my, mz)
        
        GL.glVertex3d(px,py, pz)
        GL.glVertex3d(px,my, pz)
        GL.glVertex3d(px,my, mz)
        GL.glVertex3d(px,py, mz)

        GL.glVertex3d(mx,py, pz)
        GL.glVertex3d(mx,my, pz)
        GL.glVertex3d(mx,my, mz)
        GL.glVertex3d(mx,py, mz)


        GL.glEnd()
        if wireframe:
            GL.glPolygonMode( GL.GL_FRONT_AND_BACK, GL.GL_FILL );

    def cylinder(self, center_x,center_y, center_z, x_axis, y_axis, size_z, color, wireframe) :
        Pi = 3.14159265358979323846
        if wireframe:
            nbr_sides = 15
        else:
            nbr_sides = 100

        pz = center_z+size_z/2
        mz = center_z-size_z/2

        #top
        if wireframe :
            GL.glPolygonMode( GL.GL_FRONT_AND_BACK, GL.GL_LINE ) 

        GL.glBegin(GL.GL_POLYGON)
        self.qglColor(color)
        for i in range(nbr_sides):
             angle1 = (i * 2 * Pi) / nbr_sides
             x1 = center_x + x_axis/2 * math.sin(angle1)
             y1 = center_y + y_axis/2 * math.cos(angle1)
             GL.glVertex3d(x1, y1, pz)
 
        GL.glEnd()
        #bottom
        GL.glBegin(GL.GL_POLYGON)
        self.qglColor(color)
        for i in range(nbr_sides):
             angle1 = (i * 2 * Pi) / nbr_sides
             x1 = center_x + x_axis/2 * math.sin(angle1)
             y1 = center_y + y_axis/2 * math.cos(angle1)
             GL.glVertex3d(x1, y1, mz)
        GL.glEnd()

        #Sides
        GL.glBegin(GL.GL_QUADS)
        self.qglColor(color.darker(250))
        for i in range(nbr_sides):
             angle1 = (i * 2 * Pi) / nbr_sides
             angle2 = ((i+1) * 2 * Pi) / nbr_sides
             x1 = center_x + x_axis/2 * math.sin(angle1)
             y1 = center_y + y_axis/2 * math.cos(angle1)
             x2 = center_x + x_axis/2 * math.sin(angle2)
             y2 = center_y + y_axis/2 * math.cos(angle2)
             GL.glVertex3d(x1, y1, pz)
             GL.glVertex3d(x1, y1, mz)
             GL.glVertex3d(x2, y2, mz)
             GL.glVertex3d(x2, y2, pz)
        GL.glEnd()

        #if wireframe: lines in top and bottom
        if wireframe:
            GL.glBegin(GL.GL_LINES)
            self.qglColor(color)
            for i in range(nbr_sides/2):
                angle1 = (i * 2 * Pi) / nbr_sides
                angle2 = (i * 2 * Pi) / nbr_sides + Pi
                x1 = center_x + x_axis/2 * math.sin(angle1)
                y1 = center_y + y_axis/2 * math.cos(angle1)
                x2 = center_x + x_axis/2 * math.sin(angle2)
                y2 = center_y + y_axis/2 * math.cos(angle2)
                GL.glVertex3d(x1, y1, pz)
                GL.glVertex3d(x2, y2, pz)
                GL.glVertex3d(x1, y1, mz)
                GL.glVertex3d(x2, y2, mz)
            GL.glEnd()

        if wireframe:
            GL.glPolygonMode( GL.GL_FRONT_AND_BACK, GL.GL_FILL );


    def plane_xy (self, center_x,center_y, center_z, size_x, size_y,color) :
        GL.glBegin(GL.GL_QUADS)
        self.qglColor(color)

        px = center_x+size_x/2
        mx =center_x-size_x/2
        py = center_y+size_y/2
        my =center_y-size_y/2

        GL.glVertex3d(px,py, center_z)
        GL.glVertex3d(mx,py, center_z)
        GL.glVertex3d(mx,my, center_z)
        GL.glVertex3d(px,my, center_z)

        GL.glEnd()

    def ellipse_xy (self, center_x, center_y, center_z, x_axis, y_axis, color):
        Pi = 3.14159265358979323846
        nbr_sides = 100

        GL.glBegin(GL.GL_POLYGON)
        self.qglColor(color)
        for i in range(nbr_sides):
             angle1 = (i * 2 * Pi) / nbr_sides
             x1 = center_x + x_axis/2 * math.sin(angle1)
             y1 = center_y + y_axis/2 * math.cos(angle1)
             GL.glVertex3d(x1, y1, center_z)
        GL.glEnd()

    def normalizeAngle(self, angle):
        while angle < 0:
            angle += 360 * 16
        while angle > 360 * 16:
            angle -= 360 * 16
        return angle

    def drawPointAndDirs(self):

        xAxis = numpy.array([1,0,0], 'float')
        yAxis = numpy.array([0,1,0], 'float')
        zAxis = numpy.array([0,0,1], 'float')
        center =numpy.array([0,0,0], 'float')



        qobj = GLU.gluNewQuadric()
        sphere_radius = 0.1
        cyl_height = 1
        cyl_radius = 0.01
        z_dir = numpy.array([0.,0.,1.], dtype = 'float64')

        genList = GL.glGenLists(1)
        GL.glNewList(genList, GL.GL_COMPILE)

        #Cylinder along z, we want it to be allong xAxis
        #Solution: use the appropriate rotation matrix
        #We have to choose the right axis, perpendicular to both, and rotate by the appropriate angle, the angle between the 2 vectors

        self.qglColor(QtGui.QColor.fromCmykF(0., 1., 1., 0.))

        rotation_dir = numpy.cross(z_dir,xAxis)
        angle = math.acos(z_dir.dot(xAxis))*180/3.14159

        GL.glPushMatrix()
        GL.glTranslatef(center[0],center[1],center[2])
        GL.glRotatef(angle, rotation_dir[0],rotation_dir[1],rotation_dir[2])
        GLU.gluCylinder(qobj,cyl_radius,cyl_radius,cyl_height,20,20)
        GL.glPopMatrix()

        self.qglColor(QtGui.QColor.fromCmykF(1., 0., 1., 0.))


        rotation_dir = numpy.cross(z_dir,yAxis)
        angle = math.acos(z_dir.dot(yAxis))*180/3.14159

        GL.glPushMatrix()
        GL.glTranslatef(center[0],center[1],center[2])
        GL.glRotatef(angle, rotation_dir[0],rotation_dir[1],rotation_dir[2])
        GLU.gluCylinder(qobj,cyl_radius,cyl_radius,cyl_height,20,20)
        GL.glPopMatrix()


        GL.glEndList()

        return genList
if __name__ == '__main__':
     app = QtWidgets.QApplication(sys.argv)
     electrode_struct = sim_struct.ELECTRODE()
     ex = Ui_Electrode(electrode_struct)
     ex.show()
     sys.exit(app.exec_())

