# -*- coding: utf-8 -*-\
'''
    Calculates the standard 10/10 electrode position in a head model, given the references
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.
    
    Copyright (C) 2015 - 2018  Guilherme B Saturnino, Axel Thielscher

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
import numpy

from simnibs import __version__
import csv  # csv is in the standard library


class EEG:
    ''' Class defing EEG electrode positions

    Parameters
    --------------
    surface: simnibs.msh.surface.Surface
        Surface structure with skin triangles
    pos: dict
        Dictionary with electrode name and positions
    references:
        List of fiducial names
    '''
    def __init__(self, surface):
        self.pos = {}
        self.surface = surface
        self.references = ['Nz', 'Iz', 'LPA', 'RPA']

    def print2csv(self, fn):
        with open(fn, 'w') as f:
            writer = csv.writer(f)
            for n, p in self.pos.items():
                if n in self.references:
                    type_ = 'Fiducial'
                elif n == 'Cz':
                    type_ = 'ReferenceElectrode'
                else:
                    type_ = 'Electrode'
                writer.writerow([type_] + list(p) + [n])

    def print2geo(self, fn):
        with open(fn, 'w') as f:
            f.write('View"' + 'EEG' + '"{\n')
            for p in self.pos.values():
                f.write("SP(" + ", ".join([str(i) for i in p]) + "){0};\n")
            f.write("};\n")
            f.write("myView = PostProcessing.NbViews-1;\n")
            f.write("View[myView].PointType=1; // spheres\n")
            f.write("View[myView].PointSize=" + str(6) + ";")
            
            normals = self.get_normals()
            distance = 2.5
            f.write('View"' + 'EEG_names' + '"{\n')
            for t, p in self.pos.items():
                text_position = p + distance * normals[t]
                f.write('T3(' + ', '.join([str(i) for i in text_position]) + ', TextAttributes("FontSize", "24")){"' + t + '"};\n')
            f.write("};\n")
            f.write("myView = PostProcessing.NbViews-1;\n")
            f.write("View[myView].PointType=1; // spheres\n")
            f.write("View[myView].PointSize=" + str(6) + ";")

    def get_normals(self):
        normals = {}
        for name in self.pos.keys():
            _, normals[name] = self.surface.projectPoint(self.pos[name])
        return normals


class Arch:
    def __init__(self, R, C, begin, end, matrix):
        self.R = R
        # original_space
        # begin and end of circle. the center is also the origin of the transform space
        self.C = C
        self.begin = begin
        self.end = end
        self.transf_matrix = matrix
        # polar coordinates (just angle)
        self.begin_angle = None
        self.end_angle = None
        self.findBeginEndAngles()

    def findPolarAngle(self, point):
        point_transf = point - self.C
        point_transf = point_transf.dot(self.transf_matrix)
        x = point_transf[0]
        y = point_transf[1]
        # print x,y
        # if abs(y)<1e-3: y = 1e-3
        angle = math.atan2(y, x)
        # 3rd quadrant
        if angle < (-3.1415 / 2.0): angle += 2 * 3.1415

        return angle

    # finds the coordinate that's a given fraction of the arch
    def findFraction(self, fraction):
        angle = self.begin_angle + fraction * (self.end_angle - self.begin_angle)
        coor_circle = numpy.array([self.R * math.cos(angle), self.R * math.sin(angle), 0.], 'float')
        inv_transf = numpy.transpose(self.transf_matrix)
        coor_transformed = coor_circle.dot(inv_transf)
        coor_transformed = coor_transformed + self.C
        return coor_transformed

    def findBeginEndAngles(self):
        self.begin_angle = self.findPolarAngle(self.begin)
        self.end_angle = self.findPolarAngle(self.end)


# Find positions, according toi
def FindPositions(Nazion, Inion, LPA, RPA, surface, NE_cap=False):
    ''' Find the EEG positions based in the fiducials

    Parameters
    -----------
    Nazion: 3x1 ndarray
        Position of Nazion fiducial
    Inion: 3x1 ndarray
        Position of inion fiducial
    LPA: 3x1 ndarray
        Position of left preauricular point
    RPA: 3x1 ndarray
        Position of right preauricular point
    surface: simnibs.msh.surface.Surface
        Surface structure containing the skin triangles
    NE_cap: bool (optional)
        Wether to use the NeuroElectrics cap
    
    Returns
    -----------
    eeg: EEG
        EEG structure with the electrode positions
    '''
    # eeg = EEG_TenTen()
    eeg = EEG(surface)
    eeg.pos['Nz'] = numpy.array(Nazion, 'float')
    eeg.pos['Iz'] = numpy.array(Inion, 'float')
    eeg.pos['LPA'] = numpy.array(LPA, 'float')
    eeg.pos['RPA'] = numpy.array(RPA, 'float')

    findCz(eeg)
    if NE_cap:
        sagittalPositions_NE(eeg)
        coronalPositions_NE(eeg)
        otherPositions_NE(eeg)

    else:
        sagittalPositions(eeg)
        coronalPositions2(eeg)
        otherPositions(eeg)


    return eeg


def avgSideSize(surface):
    avg_side = 0.0
    nr_triangles = 0.0
    for tr in surface.triangles:
        side = tr.V1 - tr.V0
        avg_side += numpy.linalg.norm(side)

        nr_triangles += 1

    # print avg_side
    avg_side /= nr_triangles

    return avg_side


# Traces a circle that goes through 3 given points
def traceArch2(Begin_Point, End_Point, Mid_Point):
    # define plane with the 3 points
    u = Begin_Point - End_Point
    v = Mid_Point - End_Point
    if numpy.linalg.norm(u) < 1e-7 or numpy.linalg.norm(v) < 1e-7:
        raise ValueError('Invalid Points!')
    u /= numpy.linalg.norm(u)
    v /= numpy.linalg.norm(v)
    # orthogonalize
    # if the vectors are paralell
    if u.dot(v) > 0.9999:
        v = v - u.dot(v) * u
        v /= numpy.linalg.norm(v)
        raise ValueError('Orthogonalization Error!!\n')

    v = v - u.dot(v) * u
    v /= numpy.linalg.norm(v)
    normal = numpy.cross(u, v)


    # transform the 3 points in the new coordinase system:
    transformation = numpy.transpose(numpy.array([u, v, normal]))
    inv_transfo = numpy.array([u, v, normal])

    P1 = numpy.array([0., 0., 0.], 'float')
    P2 = Begin_Point - End_Point
    P2 = P2.dot(transformation)
    P3 = Mid_Point - End_Point
    P3 = P3.dot(transformation)

    # gets the parameters of the circle that goes through the 3 points
    # http://mathworld.wolfram.com/Circle.html

    D = 3 * [0.]
    D[0] = P1.dot(P1)
    D[1] = P2.dot(P2)
    D[2] = P3.dot(P3)

    A = numpy.array([[P1[0], P1[1], 1], [P2[0], P2[1], 1], [P3[0], P3[1], 1]])
    B = numpy.array([[D[0], P1[1], 1], [D[1], P2[1], 1], [D[2], P3[1], 1]])
    E = numpy.array([[D[0], P1[0], 1], [D[1], P2[0], 1], [D[2], P3[0], 1]])
    F = numpy.array([[D[0], P1[0], P1[1]], [D[1], P2[0], P2[1]], [D[2], P3[0], P3[1]]])

    a = numpy.linalg.det(A)
    b = -numpy.linalg.det(B)
    e = numpy.linalg.det(E)
    f = -numpy.linalg.det(F)

    if abs(a) < 1e-7:
        print("Fit failed!")
        return None

    R = math.sqrt((b * b + e * e) / (4 * a * a) - f / a)
    Center = numpy.array([-b / (2 * a), -e / (2 * a), 0.])
    Center = Center.dot(inv_transfo) + End_Point

    arch = Arch(R, Center, Begin_Point, End_Point, transformation)

    return arch


# Finds Cz iteractively
# Cz is defined as the mid-point of both curves Nz-Iz and LPA-RPA
def findCz(eeg):
    # first estimate: projection of the central point in the head surface, perpendicular to LPA-RPA and Cz-Iz lines
    central_point = (eeg.pos['LPA'] + eeg.pos['RPA'] + eeg.pos['Nz'] + eeg.pos['Iz']) * 0.25
    direction = numpy.cross(eeg.pos['LPA'] - eeg.pos['RPA'], eeg.pos['Iz'] - eeg.pos['Nz'])
    direction /= numpy.linalg.norm(direction)
    # Iz-Nz is just a size estimate
    out_point = numpy.linalg.norm(eeg.pos['Iz'] - eeg.pos['Nz']) * direction + central_point
    eeg.pos['Cz'], tmp = eeg.surface.projectPoint(out_point, central_point)

    iterations = 0
    while True:
        Cz_Old = eeg.pos['Cz']

        ArchNz_Iz = traceArch2(eeg.pos['Nz'], eeg.pos['Iz'], eeg.pos['Cz'])

        eeg.pos['Cz'] = ArchNz_Iz.findFraction(0.5)

        ArchLPA_RPA = traceArch2(eeg.pos['LPA'], eeg.pos['RPA'], eeg.pos['Cz'])
        eeg.pos['Cz'] = ArchLPA_RPA.findFraction(0.5)
        iterations += 1
        Cz_distance = numpy.linalg.norm(eeg.pos['Cz'] - Cz_Old)

        if Cz_distance < 0.1 or iterations == 20:
            break


def sagittalPositions(eeg):
    positions = ['Fpz', 'AFz', 'Fz', 'FCz', 'Cz', 'CPz', 'Pz', 'POz', 'Oz']
    equalySpaced(eeg, eeg.pos['Nz'], eeg.pos['Iz'], eeg.pos['Cz'], positions)


def coronalPositions2(eeg):
    LPA_RPA = ['T7', 'C5', 'C3', 'C1', 'Cz', 'C2', 'C4', 'C6', 'T8']
    equalySpaced(eeg, eeg.pos['LPA'], eeg.pos['RPA'], eeg.pos['Cz'], LPA_RPA)

    coronal_0_left = ['N1', 'AF9', 'F9', 'FT9', 'T9', 'TP9', 'P9', 'PO9', 'I1']
    equalySpaced(eeg, eeg.pos['Nz'], eeg.pos['Iz'], eeg.pos['LPA'], coronal_0_left)

    coronal_0_right = ['N2', 'AF10', 'F10', 'FT10', 'T10', 'TP10', 'P10', 'PO10', 'I2']
    equalySpaced(eeg, eeg.pos['Nz'], eeg.pos['Iz'], eeg.pos['RPA'], coronal_0_right)

    coronal_10_left = ['Fp1', 'AF7', 'F7', 'FT7', 'T7', 'TP7', 'P7', 'PO7', 'O1']
    equalySpaced(eeg, eeg.pos['Fpz'], eeg.pos['Oz'], eeg.pos['T7'], coronal_10_left)

    coronal_10_right = ['Fp2', 'AF8', 'F8', 'FT8', 'T8', 'TP8', 'P8', 'PO8', 'O2']
    equalySpaced(eeg, eeg.pos['Fpz'], eeg.pos['Oz'], eeg.pos['T8'], coronal_10_right)


def otherPositions(eeg):
    AF = ['AF5', 'AF3', 'AF1', 'AFz', 'AF2', 'AF4', 'AF6']
    equalySpaced(eeg, eeg.pos['AF7'], eeg.pos['AF8'], eeg.pos['AFz'], AF)

    F = ['F5', 'F3', 'F1', 'Fz', 'F2', 'F4', 'F6']
    equalySpaced(eeg, eeg.pos['F7'], eeg.pos['F8'], eeg.pos['Fz'], F)

    FC = ['FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6']
    equalySpaced(eeg, eeg.pos['FT7'], eeg.pos['FT8'], eeg.pos['FCz'], FC)

    CP = ['CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6']
    equalySpaced(eeg, eeg.pos['TP7'], eeg.pos['TP8'], eeg.pos['CPz'], CP)

    P = ['P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6']
    equalySpaced(eeg, eeg.pos['P7'], eeg.pos['P8'], eeg.pos['Pz'], P)

    PO = ['PO5', 'PO3', 'PO1', 'POz', 'PO2', 'PO4', 'PO6']
    equalySpaced(eeg, eeg.pos['PO7'], eeg.pos['PO8'], eeg.pos['POz'], PO)


# Gets equaly spaced points
def equalySpaced(eeg, begin, end, middle, position_names):
    interval = 1.0 / (len(position_names) + 1)
    curve = traceArch2(begin, end, middle)

    # print middle, '\n',curveA.R, curveA.C, '\n', curve.R, curve.C,'\n\n'
    fraction = interval
    for name in position_names:
        in_arch = curve.findFraction(fraction)
        fraction += interval
        eeg.pos[name], _ = eeg.surface.projectPoint(in_arch, curve.C)
        if eeg.pos[name] is None:
            raise Exception("Error finding EEG positions: references invalid?")

def sagittalPositions_NE(eeg):
    positions = ['Fpz', 'AFz', 'Fz', 'FCz', 'Cz', 'CPz', 'Pz', 'POz', 'Oz']
    spacing = [.113, .099, .105, .091, .098, .103, .091, .098, .101]
    unequalySpaced(eeg, eeg.pos['Nz'], eeg.pos['Iz'], eeg.pos['Cz'], positions, spacing)


def coronalPositions_NE(eeg):
    LPA_RPA = ['T7', 'C5', 'C3', 'C1', 'Cz', 'C2', 'C4', 'C6', 'T8']
    equalySpaced(eeg, eeg.pos['LPA'], eeg.pos['RPA'], eeg.pos['Cz'], LPA_RPA)

    # 10% reference
    coronal_10_left = ['Fp1', 'AF7', 'F7', 'FT7', 'TP7', 'P7', 'PO7', 'O1']
    spacing = [.12, .118, .113, .097, .093 + .082, .086, .094, .101]
    unequalySpaced(eeg, eeg.pos['Fpz'], eeg.pos['Oz'], eeg.pos['T7'],
                   coronal_10_left, spacing)

    coronal_10_right = ['Fp2', 'AF8', 'F8', 'FT8', 'TP8', 'P8', 'PO8', 'O2']
    spacing = [.12, .122, .123, .090, .087 + .081, .085, .096, .103]
    unequalySpaced(eeg, eeg.pos['Fpz'], eeg.pos['Oz'], eeg.pos['T8'],
                   coronal_10_right, spacing)

    # 0% reference
    unequalySpaced(eeg, eeg.pos['Nz'], eeg.pos['Iz'], eeg.pos['LPA'], ['P9'], [.7])
    unequalySpaced(eeg, eeg.pos['Nz'], eeg.pos['Iz'], eeg.pos['RPA'], ['P10'], [.7])



def otherPositions_NE(eeg):
    AF = ['AF3', 'AF4']
    spacing = [.19, .31+.30]
    unequalySpaced(eeg, eeg.pos['AF7'], eeg.pos['AF8'], eeg.pos['AFz'],
                   AF, spacing)

    PO = ['PO3', 'PO4']
    spacing = [.30, .20+.20]
    unequalySpaced(eeg, eeg.pos['PO7'], eeg.pos['PO8'], eeg.pos['POz'],
                 PO, spacing)

    F = ['F5', 'F3', 'F1', 'Fz', 'F2', 'F4', 'F6']
    equalySpaced(eeg, eeg.pos['F7'], eeg.pos['F8'], eeg.pos['Fz'], F)

    FC = ['FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6']
    equalySpaced(eeg, eeg.pos['FT7'], eeg.pos['FT8'], eeg.pos['FCz'], FC)

    CP = ['CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6']
    equalySpaced(eeg, eeg.pos['TP7'], eeg.pos['TP8'], eeg.pos['CPz'], CP)

    P = ['P7', 'P5', 'P3', 'P1', 'P2', 'P4', 'P6', 'P8']
    spacing = [.125, .112, .0853, .105, .082+.083, .0952, .0952, .1]
    unequalySpaced(eeg, eeg.pos['P9'], eeg.pos['P10'], eeg.pos['Pz'],
                   P, spacing)


def unequalySpaced(eeg, begin, end, middle, position_names, spacing):
    assert len(position_names) == len(spacing)
    curve = traceArch2(begin, end, middle)
    fraction = 0.
    for sp, name in zip(spacing, position_names):
        fraction += sp
        in_arch = curve.findFraction(fraction)
        eeg.pos[name], _= eeg.surface.projectPoint(in_arch, curve.C)
        if eeg.pos[name] is None:
            raise Exception("Error finding EEG positions: references invalid?")


