#!/usr/bin/python
# -*- coding: utf-8 -*-\
'''
    Head Model in OpenGL for the SimNIBS GUI
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2018 Guilherme B Saturnnino

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

import math
import sys
import time
import numpy
from PyQt5 import QtCore, QtGui, QtWidgets

try:
    from OpenGL import GLUT
except NotImplementedError as e:
    print(f"{e}, setting PYOPENGL_PLATFORM=osmesa")
    import os
    os.environ['PYOPENGL_PLATFORM'] = 'osmesa'
    # clean up opengl imports and try again
    for module in tuple(filter(lambda m: m.startswith("OpenGL"), sys.modules)):
        del sys.modules[module]
    from OpenGL import GLUT

from OpenGL import GL, GLU
import OpenGL

from simnibs.simulation.tms_coil.tms_coil import TmsCoil
from ..mesh_tools import surface, mesh_io
from ..utils.csv_reader import read_csv_positions

global YELLOW
global BLUE
global RED
global GREEN
global BLACK
global GRAY
global WHITE
global PURPLE
global DARK_RED
global WHEAT

YELLOW = QtGui.QColor.fromCmykF(0., 0., 1., 0.)
BLUE = QtGui.QColor.fromCmykF(0.72, 0.52, 0., 0.)
RED = QtGui.QColor.fromCmykF(0., 1., 1., 0.)
GREEN = QtGui.QColor.fromCmykF(1., 0., 1., 0.)
BLACK = QtGui.QColor.fromCmykF(0., 0., 0., 1.0)
GRAY = QtGui.QColor.fromCmykF(0., 0., 0., 0.34)
WHITE = QtGui.QColor.fromCmykF(0., 0., 0., 0.)
PURPLE = QtGui.QColor.fromCmykF(0., 1.0, 0., 0.5)
DARK_RED = QtGui.QColor.fromCmykF(0., 1., 1., 0.45)
WHEAT = QtGui.QColor.fromCmykF(0., 0.09, 0.27, 0.04)


@QtCore.pyqtSlot(int)
@QtCore.pyqtSlot(int)
class GLHeadModel(QtWidgets.QOpenGLWidget):
    windowClicked = QtCore.pyqtSignal(int)
    loadStage = QtCore.pyqtSignal(int)

    def __init__(self, parent=None):
        super(GLHeadModel, self).__init__(parent)

        self.mesh_fn = ''
        #Objects
        self.model = 0
        self.indicator = 0
        self.stimulator_objects = []
        self.tmp_objects = []
        self.eegReferences = 0
        self.eegPositions = 0
        self.eeg_coordinates = None
        self.eeg_names = None
        self.eeg_cap = None
        #self.axis = self.drawAxis()


        #Scene parameters
        self.near = [0,0,0]
        self.far = [0,0,0]

        self.xTran = 0.0
        self.yTran = -20.0
        self.zTran = -400

        self.xRot = 90*16
        #self.xRot = 0
        self.yRot = 0
        self.zRot = 0

        self.zoom = 1.0


        #self.skin_triangles = []
        self.skin_surf = None
        self.skin_model = None
        self.skin_model_field = None
        #self.gm_triangles = []
        self.gm_surf = None
        self.gm_model = None
        self.gm_model_field = None
        self.skin_color = [135./255, 204./255, 238./255]
        self.gm_color = [160./255, 160./255, 160./255]



        self.currenSurface = ''


        self.figure_centre = None
        self.mesh_name = None

        self.hasField = False

        #holds the point in surface which was clicked
        self.intersect_point = None
        self.intersect_normal = []



        self.lastPos = QtCore.QPoint()

        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)

    #initializes the visualization from the file name
    #call processEvents in order to maintain responsiveness
    def loadMesh(self, mesh_fn):
        print("reading ", mesh_fn)
        self.skin_surf = 'Loading'
        self.gm_surf = 'Loading'
        self.mesh_fn = mesh_fn
        self.loadStage.emit(0)
        QtWidgets.QApplication.processEvents()
        mesh_struct = mesh_io.read_msh(mesh_fn)

        self.loadStage.emit(1)
        QtWidgets.QApplication.processEvents()
        self.skin_surf = surface.Surface(mesh_struct, [5,1005])

        self.loadStage.emit(2)
        QtWidgets.QApplication.processEvents()
        self.gm_surf = surface.Surface(mesh_struct, [2,1002])
        self.loadStage.emit(3)

        QtWidgets.QApplication.processEvents()

        self.skin_surf.mesh_name = mesh_fn

    def drawSkinAndGm(self):
        self.loadStage.emit(3)
        self.skin_model = self.drawModel('Skin')
        #QtWidgets.QApplication.processEvents()
        self.gm_model = self.drawModel('GM')
        self.loadStage.emit(4)
        #QtWidgets.QApplication.processEvents()
        self.selectRenderSurface('Skin')





    def getSurface(self, surf):
        if surf == 'Scalp':
            return self.skin_surf
        elif surf == 'GM':
            return self.gm_surf
        else:
            print('Invalid Surface Name!!')


    def minimumSizeHint(self):
        return QtCore.QSize(300, 300)

    def sizeHint(self):
        #return QtCore.QSize(450, 450)
        return QtCore.QSize(700, 700)

    def setXRotation(self, angle):
        angle = self.normalizeAngle(angle)
        if angle != self.xRot:
            self.xRot = angle
            self.update()

    def setYRotation(self, angle):
        angle = self.normalizeAngle(angle)
        if angle != self.yRot:
            self.yRot = angle
            self.update()

    def setZRotation(self, angle):
        angle = self.normalizeAngle(angle)
        if angle != self.zRot:
            self.zRot = angle
            self.update()

    def setXTranslation(self, distance):
        if distance != 0.0:
            self.xTran = distance
            self.update()

    def setYTranslation(self, distance):
        if distance != 0.0:
            self.yTran = distance
            self.update()

    def getOpenglInfo(self):
        info = """
            Vendor: {0}
            Renderer: {1}
            OpenGL Version: {2}
            Shader Version: {3}
            """.format(
            GL.glGetString(GL.GL_VENDOR).decode(),
            GL.glGetString(GL.GL_RENDERER).decode(),
            GL.glGetString(GL.GL_VERSION).decode(),
            GL.glGetString(GL.GL_SHADING_LANGUAGE_VERSION).decode()
        )
        return info

    def initializeGL(self):
        self.setClearColor(WHITE)
        GL.glShadeModel(GL.GL_SMOOTH)
        light_ambient =  [0.0, 0.0, 0.0, 1.0]
        GL.glLightfv(GL.GL_LIGHT0, GL.GL_AMBIENT, light_ambient)
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_LIGHTING)
        GL.glEnable(GL.GL_LIGHT0)
        GL.glEnable(GL.GL_COLOR_MATERIAL)
        GL.glEnable(GL.GL_NORMALIZE)
        GL.glColorMaterial(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE)

    def setClearColor(self, c):
        GL.glClearColor(c.redF(), c.greenF(), c.blueF(), c.alphaF())

    def paintGL(self):
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        width = float(self.width())
        height = float(self.height())
        frustrumx = 60*width/self.sizeHint().width()
        frustrumy = 60*height/self.sizeHint().height()
        GL.glFrustum(-frustrumx, frustrumx, frustrumy, -frustrumy, 200, 500)
        GL.glMatrixMode(GL.GL_MODELVIEW)

        self.drawAxis()

        GL.glLoadIdentity()
        GL.glTranslatef(self.xTran, self.yTran, self.zTran)
        GL.glRotatef(self.xRot / 16.0, 1.0, 0.0, 0.0)
        GL.glRotatef(self.yRot / 16.0, 0.0, 1.0, 0.0)
        GL.glRotatef(self.zRot / 16.0, 0.0, 0.0, 1.0)
        GL.glScalef(-self.zoom,self.zoom,self.zoom)


        if self.model != 0:
            GL.glCallList(self.model)

        if self.indicator:
            GL.glCallList(self.indicator)

        for tmp in self.tmp_objects:
            try:
                GL.glCallList(tmp)
            except:
                pass

        try:
            GL.glCallList(self.eegReferences)
        except:
            pass

        try:
            GL.glCallList(self.eegPositions)
        except:
            pass

        for stimulator in self.stimulator_objects:
            try:
                GL.glCallList(stimulator)
            except:
                pass

        self.model_matrix = GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX)
        self.projection_matrix = GL.glGetDoublev(GL.GL_PROJECTION_MATRIX)

    def resizeGL(self, width, height):
        GL.glViewport(0, 0, width, height)
        self.view = [0, 0, width, height]


    def getIntersectPoint(self):
        return self.intersect_point

    #double clickng gets the position
    def mouseDoubleClickEvent(self, event):
        if isinstance(self.skin_surf,  surface.Surface):
            self.lastPos = QtCore.QPoint(event.pos())
            x = float(self.lastPos.x())
            y = float(self.view[3] - self.lastPos.y())
            Near = GLU.gluUnProject(
                x, y, 0.,
                self.model_matrix, self.projection_matrix, self.view)
            Far = GLU.gluUnProject(
                x, y, 1.,
                self.model_matrix, self.projection_matrix, self.view)
            self.intersect_point, self.intersect_normal = self.skin_surf.interceptRay(Near, Far)
            if self.intersect_point is not None:
                self.indicator = self.drawIndicator(self.intersect_point, self.intersect_normal)
                self.update()
            self.windowClicked.emit(1)


    def wheelEvent(self, event):
        delta = event.angleDelta().y()
        zoom = self.zoom+delta/1200.0
        if zoom > 0.1 and zoom < 10:
            self.zoom = zoom
            self.update()

    #gets a new mouse position, for a smoother rotation/translation
    def mousePressEvent(self,event):
        self.lastPos = QtCore.QPoint(event.pos())

    #Rotates / translates model
    def mouseMoveEvent(self, event):
        dx = event.x() - self.lastPos.x()
        dy = event.y() - self.lastPos.y()
        size = self.size()
        size_hint = self.sizeHint()
        dx = dx*size_hint.width()/float(size.width())
        dy = dy*size_hint.height()/float(size.height())

        if event.buttons() & QtCore.Qt.LeftButton:
            self.setXRotation(self.xRot - 8 * dy)
            self.setZRotation(self.zRot - 8 * dx)
        elif event.buttons() & QtCore.Qt.RightButton:
            self.setXTranslation(self.xTran + dx*0.5)
            self.setYTranslation(self.yTran + dy*0.5)

        self.lastPos = QtCore.QPoint(event.pos())



    def normalizeAngle(self, angle):
        while angle < 0:
            angle += 360 * 16
        while angle > 360 * 16:
            angle -= 360 * 16
        return angle



    #Selects if GM or Skin will be rendered
    def selectRenderSurface(self, surf_name):
        if surf_name == 'Skin':
            self.currenSurface = 'Skin'
            if not self.hasField:
                self.model = self.skin_model
            else:
                self.model = self.skin_model_field
        if surf_name == 'GM':
            self.currenSurface = 'GM'
            if not self.hasField:
                self.model = self.gm_model
            if self.hasField:
                self.model = self.gm_model_field

        self.update()


    ##Renders the surface
    #if it has fields, also creates a heatMap
    def drawModel(self, surf_name, field= None):

        if surf_name == 'Skin':
            rendered_surf = self.skin_surf
            color = self.skin_color
        elif surf_name == 'GM':
            rendered_surf = self.gm_surf
            color = self.gm_color


        else:
            print('Invalid argument at drawModel')
            sys.exit(1)

        if not isinstance(rendered_surf, surface.Surface):
            return 0


        if len(rendered_surf.nodes) == 0:
            print("ERROR: Could not find surface", surf_name, "in mesh file\n")
            return 0



        #Buffers the node positions and normals
        nodes_pos = numpy.array(rendered_surf.nodes, dtype='float32')
        node_normals = numpy.array([normal for normal in rendered_surf.nodes_normals], dtype = 'float32')
        nr_nodes = len(nodes_pos)

        #Creates heat map
        if self.hasField:
            blue = numpy.array([0,0,1], dtype='float32')
            red = numpy.array([1,0,0], dtype='float32')
            green = numpy.array([0,1,0], dtype='float32')

            #Kernel: lowest value: [0,0,255], middle value: [0,255,0], highest value:[0,0,255]
            #The rest are interpolated
            #http://stackoverflow.com/questions/20792445/calculate-rgb-value-for-a-range-of-values-to-create-heat-map
            max_v = numpy.max(field)
            min_v = numpy.min(field)
            #norm_factor = max_v-min_v
            norm_field = (field-min_v)/(max_v-min_v) #normalized field
            b = numpy.zeros(len(field))
            g = numpy.zeros(len(field))
            r = numpy.zeros(len(field))
            '''
            b = 1-2*(field-min_v)/norm_factor
            r = -1+2*(field-min_v)/norm_factor
            b = numpy.maximum(b, numpy.zeros([1,nr_nodes]))
            r = numpy.maximum(r, numpy.zeros([1,nr_nodes]))
            g = 1-b-r
            '''
            octile = norm_field<0.125
            b[octile] = 4*norm_field[octile]+.5

            octile = (norm_field>0.125) & (norm_field<=3*0.125)
            g[octile] = 4*norm_field[octile]-.5
            b[octile] = 1

            octile = (norm_field>3*0.125) & (norm_field<=5*0.125)
            b[octile] = -4*norm_field[octile]+2.5
            g[octile] = 1
            r[octile] = 4*norm_field[octile]-1.5


            octile = (norm_field>5*0.125) & (norm_field<=7*0.125)
            g[octile] = -4*norm_field[octile]+3.5
            r[octile] = 1

            octile = (norm_field>7*0.125)
            r[octile] = -4*norm_field[octile]+4.5

            node_colors = (r*red[:,numpy.newaxis] + b*blue[:,numpy.newaxis] + g*green[:,numpy.newaxis]).T
        else:
            node_colors = numpy.array([color, ]*nr_nodes,  dtype = 'float32')

        GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
        GL.glEnableClientState(GL.GL_NORMAL_ARRAY)
        GL.glEnableClientState(GL.GL_COLOR_ARRAY)

        GL.glVertexPointerf(nodes_pos)
        GL.glNormalPointerf(node_normals)
        GL.glColorPointerf(node_colors)


        genList = GL.glGenLists(1)
        GL.glNewList(genList, GL.GL_COMPILE)
        #self.qglColor(BLUE)

        #Call the triangles
        for vertices in rendered_surf.tr_nodes:
            GL.glDrawElements(GL.GL_TRIANGLES, 3, GL.GL_UNSIGNED_INT, vertices)



        if not self.hasField:
            GL.glEndList()
            GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
            GL.glDisableClientState(GL.GL_NORMAL_ARRAY)
            GL.glDisableClientState(GL.GL_COLOR_ARRAY)
            return genList

        else:
            #additionally, draws the heat map
            self.drawHeatMapScale(field)
            GL.glEndList()
            GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
            GL.glDisableClientState(GL.GL_NORMAL_ARRAY)
            GL.glDisableClientState(GL.GL_COLOR_ARRAY)

            return genList

    def qglColor(self, c):
        GL.glColor3f(c.redF(), c.greenF(), c.blueF())


    #Draws an "indicator" showing the clicked position
    def drawIndicator(self, intersect_point, normal):
        genList = GL.glGenLists(1)
        GL.glNewList(genList, GL.GL_COMPILE)
        u, v = self.skin_surf.getTangentCoordinates(normal)


        self.qglColor(RED)
        self.ellipse(intersect_point,normal, u,v,[6,6]);


        GL.glEndList()


        return genList

    #draws a simple model of an electrode
    def drawElectrode(self, elec, color = GREEN):

        if len(elec.centre)< 3:
            print('invalid center!')
            return None
        if elec.dimensions is None:
            print('invalid dimensions!')
            return None

        center = numpy.array(elec.centre, 'float')
        normal = self.skin_surf.projectPoint(center, smooth=True)[1]
        center = center+2*normal

        if normal is None:
            print('Invalid Point!')
            return None

        if len(elec.pos_ydir) != 3:
            x_axis,y_axis = self.skin_surf.getTangentCoordinates(normal)

        else:
            y_dir = numpy.array(elec.pos_ydir, 'float')
            y_axis = y_dir-center
            y_axis = y_axis-y_axis.dot(normal)*normal
            y_axis = y_axis/numpy.linalg.norm(y_axis)

            x_axis = numpy.cross(y_axis, -normal)



        genList = GL.glGenLists(1)
        GL.glNewList(genList, GL.GL_COMPILE)


        self.qglColor(color)

        if elec.shape == 'rect':
            self.rectangle(center,normal, x_axis, y_axis, elec.dimensions)


        elif elec.shape == 'ellipse':
            self.ellipse(center,normal, x_axis, y_axis, elec.dimensions)

        else:
            GL.glEndList()
            print("Unrecognized electrode shape!")
            return None

        #draws the plug
        self.qglColor(RED)
        if  len(elec.plug) != 0:
            plug_center = center+elec.plug[0].centre[1]*y_axis+elec.plug[0].centre[0]*x_axis+normal
            if elec.plug[0].shape == 'rect':
                self.rectangle(plug_center,normal, x_axis, y_axis,
                               elec.plug[0].dimensions)
            elif elec.plug[0].shape == 'ellipse':
                 self.ellipse(plug_center,normal, x_axis, y_axis,
                              elec.plug[0].dimensions)
            else:
                GL.glEndList()
                print("Unrecognized electrode shape!")
                return None

        GL.glEndList()




        return genList

    def ellipse(self, center, normal, u, v, dimensions):
        dim_x = dimensions[0]
        dim_y = dimensions[1]
        GL.glBegin(GL.GL_POLYGON)
        Pi = 3.14159265358979323846
        nbr_sides = 20
        GL.glNormal3f(normal[0],normal[1],normal[2])
        for i in range(nbr_sides):
            angle = (i * 2 * Pi) / nbr_sides
            vect = center+normal+( dim_x/2)*u*math.sin(angle)+ \
                +( dim_y/2)*v*math.cos(angle)
            GL.glVertex3d(vect[0], vect[1], vect[2])

        GL.glEnd()

    def rectangle(self, center, normal, x_axis, y_axis, dimensions):
        GL.glBegin(GL.GL_QUADS)
        GL.glNormal3f(normal[0],normal[1],normal[2])
        x_dim = dimensions[0]/2.0
        y_dim = dimensions[1]/2.0
        v1 = center+normal+x_dim*x_axis+y_dim*y_axis
        v2 = center+normal+x_dim*x_axis-y_dim*y_axis
        v3 = center+normal-x_dim*x_axis-y_dim*y_axis
        v4 = center+normal-x_dim*x_axis+y_dim*y_axis
        GL.glVertex3d(v1[0], v1[1], v1[2])
        GL.glVertex3d(v2[0], v2[1], v2[2])
        GL.glVertex3d(v3[0], v3[1], v3[2])
        GL.glVertex3d(v4[0], v4[1], v4[2])
        GL.glEnd()

    def setEEG(self, cap_fn):
        if cap_fn is None:
            self.eegReferences = 0
            self.eegPositions = 0
            self.eeg_coordinates = None
            self.eeg_names = None
            self.eeg_cap = None
            self.update()
            return
        try:
            type_, coordinates, _, name, _, _ = read_csv_positions(cap_fn)
        except:
            raise IOError('Could not read EEG position file: ' + self.eeg_cap)
        self.eeg_cap = cap_fn
        self.eeg_coordinates = coordinates
        self.eeg_names = name
        self.drawEEGPositions(list(coordinates), name)

    #Rest of the positions
    def drawEEGPositions(self, points, names):
        genList = GL.glGenLists(1)
        GL.glNewList(genList, GL.GL_COMPILE)
        for point, name in zip(points, names):
            normal = self.skin_surf.projectPoint(point, smooth=False)[1]
            if normal is None or len(normal) == 0 or point is None or len(point) == 0:
                print('Invalid Point!')
            else:
                u, v = self.skin_surf.getTangentCoordinates(normal)
                pos_el = point + 2*normal
                GL.glEnable(GL.GL_LIGHTING)
                self.qglColor(GREEN)
                self.ellipse(pos_el, normal, u, v, [6, 6])
                if bool(GLUT.glutBitmapString):
                    pos_text = point + 5*normal
                    GL.glDisable(GL.GL_LIGHTING)
                    GL.glColor3f(0.0, 0.0, 0.0)
                    GL.glRasterPos3f(*pos_text)
                    GLUT.glutBitmapString(
                        GLUT.GLUT_BITMAP_HELVETICA_12,
                        name.encode())
                    GL.glEnable(GL.GL_LIGHTING)

        GL.glEndList()
        self.eegPositions = genList
        self.update()

    def clear_eeg_positions(self):
        self.eegPositions = 0

    #Draws the points and directions from a 4x4 matrice, like matsimnibs
    def drawPointAndDirs(self, matrix, color=GREEN, fn_coil=''):
        try:

            xAxis = numpy.array([matrix[0][0], matrix[1][0], matrix[2][0]], 'float')
            yAxis = numpy.array([matrix[0][1], matrix[1][1], matrix[2][1]], 'float')
            zAxis = numpy.array([matrix[0][2], matrix[1][2], matrix[2][2]], 'float')
            center =numpy.array([matrix[0][3], matrix[1][3], matrix[2][3]], 'float')
        except:
            print('invalid matrix!')
            return None

        if fn_coil != '':
            coil = TmsCoil.from_file(fn_coil)
            coil_mesh = coil.get_mesh(numpy.array(matrix), include_optimization_points=False, include_coil_elements=False)
            if coil_mesh.elm.nr == 0:
                coil_genList = GL.glGenLists(1)
                GL.glNewList(coil_genList, GL.GL_COMPILE)
                GL.glEndList()
            else:
                QtWidgets.QApplication.processEvents()
                coil_mesh.elm.tag1[:] = 1
                coil_mesh.elm.tag2[:] = 1
                rendered_surf = surface.Surface(coil_mesh, [1])
                coil_color = [color.redF(), color.greenF(), color.blueF(), 0.4]
                nodes_pos = numpy.array(rendered_surf.nodes, dtype='float32')
                node_normals = numpy.array([normal for normal in rendered_surf.nodes_normals], dtype = 'float32')
                nr_nodes = len(nodes_pos)
                node_colors = numpy.array([coil_color, ]*nr_nodes,  dtype = 'float32')

                GL.glEnableClientState(GL.GL_NORMAL_ARRAY)
                GL.glEnableClientState(GL.GL_COLOR_ARRAY)
                GL.glEnableClientState(GL.GL_VERTEX_ARRAY)

                GL.glNormalPointerf(node_normals)
                GL.glColorPointer(4, GL.GL_FLOAT, 0, node_colors)
                GL.glVertexPointerf(nodes_pos)

                coil_genList = GL.glGenLists(1)
                GL.glNewList(coil_genList, GL.GL_COMPILE)
                GL.glPushAttrib (GL.GL_ALL_ATTRIB_BITS)
                GL.glEnable(GL.GL_CULL_FACE)
                GL.glCullFace(GL.GL_BACK)
                GL.glEnable(GL.GL_BLEND)
                GL.glBlendFunc( GL.GL_SRC_ALPHA,  GL.GL_ONE_MINUS_SRC_ALPHA)

                #for vertices in rendered_surf.tr_nodes:
                GL.glDrawElements(GL.GL_TRIANGLES, len(rendered_surf.tr_nodes) * 3, GL.GL_UNSIGNED_INT, numpy.ravel(rendered_surf.tr_nodes))

                GL.glPopAttrib()
                GL.glEndList()
                GL.glDisableClientState(GL.GL_NORMAL_ARRAY)
                GL.glDisableClientState(GL.GL_COLOR_ARRAY)
                GL.glDisableClientState(GL.GL_VERTEX_ARRAY)


        qobj = GLU.gluNewQuadric()
        sphere_radius = 3
        cyl_height = 20
        cyl_radius = 0.7
        z_dir = numpy.array([0.,0.,1.], dtype = 'float64')

        genList = GL.glGenLists(1)
        GL.glNewList(genList, GL.GL_COMPILE)

        self.qglColor(color)
        GL.glPushMatrix()
        GL.glTranslatef(center[0],center[1],center[2])
        GLU.gluSphere(qobj, sphere_radius, 20,20)
        GL.glPopMatrix()

        #Cylinder along z, we want it to be allong xAxis
        #Solution: use the appropriate rotation matrix
        #We have to choose the right axis, perpendicular to both, and rotate by the appropriate angle, the angle between the 2 vectors

        self.qglColor(RED)

        rotation_dir = numpy.cross(z_dir,xAxis)
        angle = math.acos(z_dir.dot(xAxis))*180/3.14159

        GL.glPushMatrix()
        GL.glTranslatef(center[0],center[1],center[2])
        GL.glRotatef(angle, rotation_dir[0],rotation_dir[1],rotation_dir[2])
        GLU.gluCylinder(qobj,cyl_radius,cyl_radius,cyl_height,20,20)
        GL.glPopMatrix()

        self.qglColor(GREEN)


        rotation_dir = numpy.cross(z_dir,yAxis)
        angle = math.acos(z_dir.dot(yAxis))*180/3.14159

        GL.glPushMatrix()
        GL.glTranslatef(center[0],center[1],center[2])
        GL.glRotatef(angle, rotation_dir[0],rotation_dir[1],rotation_dir[2])
        GLU.gluCylinder(qobj,cyl_radius,cyl_radius,cyl_height,20,20)
        GL.glPopMatrix()


        self.qglColor(BLUE)

        rotation_dir = numpy.cross(z_dir,zAxis)
        angle = math.acos(z_dir.dot(zAxis))*180/3.14159
        z_axis_center = center-zAxis*10

        GL.glPushMatrix()
        GL.glTranslatef(z_axis_center[0],z_axis_center[1],z_axis_center[2])
        GL.glRotatef(angle, rotation_dir[0],rotation_dir[1],rotation_dir[2])
        GLU.gluCylinder(qobj,cyl_radius,cyl_radius,cyl_height*2,20,20)
        GL.glPopMatrix()


        GL.glEndList()

        if fn_coil == '':
            return genList
        else:
            return genList, coil_genList



    def addElectrodeToList(self, electrode_struct, number):
        while len(self.electrodes_objects) < number+1:
            self.electrodes_objects.append(None)
        color = QtGui.QColor.fromCmykF(0.00, 0.99, 0.1, 0.0)
        center = electrode_struct.centre
        dimensions = electrode_struct.dimensions
        shape = electrode_struct.shape
        self.electrodes_objects[number] = self.drawElectrode(center, dimensions, shape, color)


    #List of stimulator objects (Electrodes or Coils)
    def stimulatorList(self, objects):
        self.stimulator_objects = objects
        #for obj in objects:
        #    self.coil_objects.append(obj)
        self.update()

    def tmpObjectList(self, objects):
        self.tmp_objects = objects
        self.update()

    def clearTmpObjects(self):
        self.tmp_objects = []
        self.update()


    #Gets the dAdt for the point in gray matter surface in the position defined by matsimnibs using the data from fn_coil
    def get_dAdtField(self, matsimnibs, fn_coil):
        #Try getting the values for plotting dA/dt
        if fn_coil.endswith('.nii') or fn_coil.endswith('.nii.gz'):
            M = numpy.array(matsimnibs)
            if self.skin_surf is not None:
                field_skin = numpy.linalg.norm(TmsCoil.from_file(fn_coil).get_a_field(numpy.array(self.skin_surf.nodes), M), axis=1)*1e7
            if self.gm_surf is not None:
                field_gm = numpy.linalg.norm(TmsCoil.from_file(fn_coil).get_a_field(numpy.array(self.gm_surf.nodes), M), axis=1)*1e7

            self.hasField = True
            if self.skin_surf is not None:
                self.skin_model_field = self.drawModel('Skin', field_skin)
            QtWidgets.QApplication.processEvents()
            if self.gm_surf is not None:
                self.gm_model_field = self.drawModel('GM', field_gm)
            QtWidgets.QApplication.processEvents()
            self.selectRenderSurface(self.currenSurface)
            self.update()

        else:
            return None



    def cleardAdtFields(self):
        self.hasField = False
        self.selectRenderSurface(self.currenSurface)


    #Static X,Y and Z Axis
    #X = Red
    #Y = Green
    #Z = Blue
    def drawAxis(self):
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glPushMatrix()
        GL.glLoadIdentity()
        width = float(self.width())
        height = float(self.height())
        frustrumx = 60*width/self.sizeHint().width()
        frustrumy = 60*height/self.sizeHint().height()
        GL.glFrustum(-frustrumx, frustrumx, frustrumy, -frustrumy, 200, 500)
        GL.glMatrixMode(GL.GL_MODELVIEW)

        GL.glPushMatrix()
        GL.glLoadIdentity()

        qobj = GLU.gluNewQuadric()
        GL.glTranslatef(1.7*frustrumx, 1.7*frustrumy, -400)
        GL.glRotatef(self.xRot / 16.0, 1.0, 0.0, 0.0)
        GL.glRotatef(self.yRot / 16.0, 0.0, 1.0, 0.0)
        GL.glRotatef(self.zRot / 16.0, 0.0, 0.0, 1.0)
        GL.glDisable(GL.GL_LIGHTING)
        GL.glBegin(GL.GL_LINES)

        self.qglColor(RED)
        GL.glVertex3f(-20.0,0,0) #X is inverted
        GL.glVertex3f(0.,0.,0.)
        self.qglColor(GREEN)
        GL.glVertex3f(0,20,0)
        GL.glVertex3f(0.,0.,0.)
        self.qglColor(BLUE)
        GL.glVertex3f(0,0,20)
        GL.glVertex3f(0.,0.,0.)

        GL.glEnd()
        GL.glEnable(GL.GL_LIGHTING)
        GL.glPopMatrix()

        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glPopMatrix()

        GL.glMatrixMode(GL.GL_MODELVIEW)




    #Scale for HeatMap
    def drawHeatMapScale(self, field):
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glPushMatrix()
        GL.glLoadIdentity()
        GL.glOrtho(-60,60,-60,60,0,200)
       # GL.glPopMatrix()

        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glPushMatrix()
        GL.glLoadIdentity()

        GL.glDisable(GL.GL_LIGHTING)

        if bool(GLUT.glutBitmapString):
            GL.glColor3f(0.0, 0.0, 0.0)
            GL.glRasterPos2i(26,-58)
            GLUT.glutBitmapString(
                GLUT.GLUT_BITMAP_HELVETICA_12,
                f"{max(field) : .2f}".encode())

            GL.glRasterPos2i(-32,-58)
            GLUT.glutBitmapString(
                GLUT.GLUT_BITMAP_HELVETICA_12,
                f"{min(field) : .2f}".encode())

            GL.glRasterPos2i(-5,-58)
            GLUT.glutBitmapString(
                GLUT.GLUT_BITMAP_HELVETICA_12,
                "dA/dt (V/m)".encode())


        top = -50
        bot = -55

        GL.glBegin(GL.GL_QUAD_STRIP)
        GL.glColor3f(0,0,0.5)
        GL.glVertex2f(-30,top)
        GL.glVertex2f(-30,bot)

        GL.glColor3f(0,0,1)
        GL.glVertex2f(-22.5,top)
        GL.glVertex2f(-22.5,bot)

        GL.glColor3f(0,1,1)
        GL.glVertex2f(-7.5,top)
        GL.glVertex2f(-7.5,bot)

        GL.glColor3f(0.5,1,0.5)
        GL.glVertex2f(0,top)
        GL.glVertex2f(0,bot)

        GL.glColor3f(1,1,0)
        GL.glVertex2f(7.5,top)
        GL.glVertex2f(7.5,bot)

        GL.glColor3f(1,0,0)
        GL.glVertex2f(22,top)
        GL.glVertex2f(22,bot)

        GL.glColor3f(0.5,0,0)
        GL.glVertex2f(30,top)
        GL.glVertex2f(30,bot)

        GL.glEnd()

        GL.glEnable(GL.GL_LIGHTING)
        GL.glPopMatrix()


        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glPopMatrix()

        GL.glMatrixMode(GL.GL_MODELVIEW)



    def changeSurfaceColors(self, surf_name, color):
        if surf_name == 'Skin':
            self.skin_color = color
            self.skin_model = self.drawModel('Skin')
            QtWidgets.QApplication.processEvents()
        elif surf_name == 'GM':
            self.gm_color = color
            QtWidgets.QApplication.processEvents()
            self.gm_model = self.drawModel('GM')
        else:
            print("ERROR in head_model_OGL: Invalid surface name:", surf_name)

        self.selectRenderSurface(self.currenSurface)
        self.update()


class HEADMODEL_UI(QtWidgets.QWidget):

    def __init__(self):
        super(HEADMODEL_UI, self).__init__()
        try:
            GLUT.glutInit()
        except OpenGL.error.NullFunctionError:
            pass
        self.referencePoints = [None]*5
        mainLayout = QtWidgets.QVBoxLayout()

        #Radio buttons for surfaces
        self.radio_buttons_box = QtWidgets.QGroupBox('')
        radio_buttons_layout = QtWidgets.QHBoxLayout()
        self.view_scalp = QtWidgets.QRadioButton("Scalp")
        self.view_scalp.toggle()
        self.view_scalp.toggled.connect(self.changeSurface)
        radio_buttons_layout.addWidget(self.view_scalp)

        self.view_gm = QtWidgets.QRadioButton("Gray Matter")
        radio_buttons_layout.addWidget(self.view_gm)

        self.radio_buttons_box.setLayout(radio_buttons_layout)



        self.glHeadModel = GLHeadModel()
        self.glHeadModel.windowClicked.connect(self.writePosition)

        self.text_label = QtWidgets.QLabel("Current Position:")
        self.pos_label = QtWidgets.QLabel("")

        self.setProgressBar()
        #self.eeg_grop_box = self.eegPosStuf()



        mainLayout.addWidget(self.radio_buttons_box)
        mainLayout.addWidget(self.glHeadModel)
        mainLayout.addWidget(self.text_label)
        mainLayout.addWidget(self.pos_label)
        mainLayout.addWidget(self.progressBar)
        #mainLayout.addWidget(self.eeg_grop_box)


        self.setLayout(mainLayout)

    def changeSurface(self):
        if self.view_gm.isChecked():
            if self.glHeadModel.gm_surf is None:
                QtWidgets.QMessageBox.critical(
                    self, 'warning',  'No GM surface found in model')
                return
            else:
                self.glHeadModel.selectRenderSurface('GM')
        else:
            if self.glHeadModel.skin_model is None:
                QtWidgets.QMessageBox.critical(
                    self, 'warning',  'No skin surface found in model')
                return
            else:
                self.glHeadModel.selectRenderSurface('Skin')


        self.glHeadModel.update()

    def writePosition(self):
        if self.glHeadModel.intersect_point is not None:
            string = self.Coords2String(self.glHeadModel.intersect_point)
            self.pos_label.setText(string)

    def Coords2String(self, Coords):
        if Coords is not None:
            strings = ['%.3f' % x for x in Coords ]
            string = ", ".join(strings)
            return string
        else:
            return None

    def setProgressBar(self):

        self.progressBar = QtWidgets.QProgressBar(self)
        self.progressBar.setValue(0)
        self.progressBar.setVisible(False)

        self.glHeadModel.loadStage.connect(self.updateProgressBar)


    #Updates the progress bar using the signals (stage) from the glHeadModel
    def updateProgressBar(self, stage):
        self.progressBar.setValue(25*stage)
        if stage == 0:
            self.progressBar.setVisible(True)
            self.progressBar.setFormat("Reading Mesh File")
        if stage == 1:
            self.progressBar.setFormat("Generating Scalp Surface")
        if stage == 2:
            self.progressBar.setFormat("Generating GM Surface")
        if stage == 3:
            self.progressBar.setFormat("Drawing Surfaces")
        if stage == 4:
            self.progressBar.setVisible(False)


if __name__ == '__main__':
     app = QtWidgets.QApplication(sys.argv)
     ex = HEADMODEL_UI()
     ex.show()
     #ex.glHeadModel.loadMesh('sphere_mm.msh')
     #ex.glHeadModel.loadMesh('almi5-binary.msh')
     #ex.glHeadModel.selectRenderSurface('Scalp')
     ex.glHeadModel.drawAxis()
     ex.glHeadModel.update()
     sys.exit(app.exec_())
