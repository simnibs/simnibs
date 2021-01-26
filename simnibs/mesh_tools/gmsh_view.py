import os
import tempfile
import subprocess
import threading

import numpy as np
from ..utils.file_finder import path2bin


class Visualization:
    ''' Defines a visualization for a 3D mesh

    Parameters
    ------------
    mesh: simnibs.msh.Msh
        Mesh to be visualized

    Attributes
    -----------
    General: simnibs.gmsh_view.General
        General visualization settings

    Mesh: simnibs.gmsh_view.Mesh
        Mesh visualization settings

    View: simnibs.gmsh_view.View of list
        Settings for views (fields), either general or one per field in the mesh.

    mesh: simnibs.msh.Msh
        Mesh to be visualized

    merge: list
        Files to be merged
    '''
    def __init__(self, mesh):
        if isinstance(mesh, str):
            if not mesh.endswith('.msh'):
                raise ValueError('mesh file name should end with .msh')
        self.General = General()
        self.Mesh = Mesh()
        self.View = View()
        self.visibility = None
        self.mesh = mesh
        self.merge = []

    def __str__(self):
        string = str(self.General) + '\n'
        string += str(self.Mesh) + '\n'
        try:
            string += '\n'.join([str(v) for v in self.View])
        except:
            string += str(self.View)
        return string

    def _apply(self):
        # Not implemented yet
        self.General._apply()
        self.Mesh._apply()
        try:
            self.View._apply()
        except AttributeError:
            [v._apply() for v in self.View]

    def show(self, new_thread=False):
        ''' Shows the mesh in gmsh

        Parameters
        -----------
        new_thread: bool (optional)
            Wether to show the mesh in a new thread. Default: False
        '''
        if self.mesh is None:
            raise ValueError('Mesh not set!')

        with tempfile.NamedTemporaryFile(suffix='.msh') as f:
            mesh_fn = f.name
        self.mesh.write(mesh_fn)
        geo_fn = mesh_fn + '.opt'
        self._write(geo_fn)
        command = [path2bin('gmsh'), mesh_fn]
        if new_thread:
            t = threading.Thread(
                target=_run, args=(command, [mesh_fn, geo_fn]))
            t.start()
        else:
            _run(command, [mesh_fn, geo_fn])

    def write_opt(self, fn_mesh):
        ''' Writes a .opt file

        Parameters
        -----------
        fn_mesh: str
            Name of mesh file ('msh'). Will add the '.opt' suffix
        '''
        fn_out = fn_mesh + ".opt"
        self._write(fn_out)

    def add_view(self, **kwargs):
        ''' adds a View entry to the current visualization object

        Parameters
        -----------
        kwargs:
            Argumets to be passed to the View object
        '''
        try:
            len(self.View)
        except TypeError:
            self.View = []
        new_view = View(len(self.View), **kwargs)
        self.View.append(new_view)
        return new_view

    def add_merge(self, fn):
        ''' adds a file to be merged to the current Visualization

        Parameters
        -----------
        fn: str
            file name
        '''
        self.merge.append(fn)


    def _write(self, fn, mode='w'):
        with open(fn, mode) as f:
            f.write('// Visualization File Created by SimNIBS\n')
            dirname = os.path.dirname(fn)
            [f.write(f'Merge "{os.path.relpath(m, dirname)}";\n') for m in self.merge]
            f.write(str(self.General))
            f.write(str(self.Mesh))
            try:
                [f.write(str(v)) for v in self.View]
            except TypeError:
                f.write(str(self.View))
            f.write(self._visibility_str())

    def _visibility_str(self):
        if self.visibility is None:
            return ''
        vis_str = ','.join([str(v) for v in self.visibility])
        s = 'Hide "*";\n'
        s += 'Show {\n'
        s += 'Volume{' + vis_str + '};\n'
        s += 'Surface{' + vis_str + '};\n'
        s += '}\n'
        return s

    def screenshot_views(self, fn_out, sleep=1):
        try:
            self.View[0]
            view = self.View
        except IndexError:
            view = [self.View]
        assert len(view) == len(fn_out), \
                'Please define one name per view'
        with tempfile.NamedTemporaryFile(suffix='.msh') as f:
            mesh_fn = f.name
        self.mesh.write(mesh_fn)
        with tempfile.NamedTemporaryFile(suffix='.geo', delete=False, mode='w') as f:
            f.write('// Visualization File Created by SimNIBS\n')
            f.write('Merge "{0}";\n'.format(mesh_fn))
            f.write(str(self.General))
            f.write(str(self.Mesh))
            for fn, v in zip(fn_out, view):
                if v.indx is None:
                    raise ValueError('Please assign an index to all views')
                fn = os.path.abspath(fn)
                f.write(str(v))
                f.write('View[{0}].Visible = 1;\n'.format(v.indx))
                f.write('Draw;\n')
                f.write('Print "{0}";\n'.format(fn))
                f.write('Sleep {0};\n'.format(sleep))
                f.write('View[{0}].Visible = 0;\n'.format(v.indx))
            f.write('Exit;\n')
            geo_fn = f.name

        _run([path2bin('gmsh'), geo_fn], [geo_fn])



class General(object):
    ''' General Gmsh visualization options.
    For more information see http://gmsh.info/doc/texinfo/gmsh.html

    Attributes
    -----------
    FieldWidth: int
        Width (in pixels) of the field window
    MaxX: float
        Maximum model coordinate along the X-axis (read-only)
    MaxY: float
        Maximum model coordinate along the Y-axis (read-only)
    MaxZ: float
        Maximum model coordinate along the Z-axis (read-only)
    MinX: float
        Maximum model coordinate along the X-axis (read-only)
    MinY: float
        Maximum model coordinate along the Y-axis (read-only)
    MinZ: float
        Maximum model coordinate along the Z-axis (read-only)
    RotationX: float
        First Euler angle (used if Trackball=0)
    RotationY: float
        Second Euler angle (used if Trackball=0)
    RotationZ: float
        Third Euler angle (used if Trackball=0)
    ScaleX: float
        X-axis scale factor
    ScaleY: float
        Y-axis scale factor
    ScaleZ: float
        Z-axis scale factor
    Trackball: int
        Use trackball rotation mode
    TrackballQuaternion0: float
        First trackball quaternion component (used if General.Trackball=1)
    TrackballQuaternion1: float
        Second trackball quaternion component (used if General.Trackball=1)
    TrackballQuaternion2: float
        Third trackball quaternion component (used if General.Trackball=1)
    TrackballQuaternion3: float
        Fourth trackball quaternion component (used if General.Trackball=1)
    TranslationX: float
        X-axis translation (in model units)
    TranslationY: float
        Y-axis translation (in model units)
    TranslationZ: float
        Z-axis translation (in model units)
    VisibilityPositionX: int
        Horizontal position (in pixels) of the upper left corner of the visibility window
    VisibilityPositionY: int
        Vertical position (in pixels) of the upper left corner of the visibility window
    VectorType: int
        Default vector display type (for normals, etc.)
    SmallAxes: int
        Display the small axes

    References
    ------------
     `Gmsh documentation <http://gmsh.info/doc/texinfo/gmsh.html>`_
    '''
    def __init__(self, **kwargs):
        self.FieldWidth = 449
        self.MaxX = 15.0
        self.MaxY = 150.0
        self.MaxZ = 150.0
        self.MinX = -150.0
        self.MinY = -150.0
        self.MinZ = -150.0
        self.RotationX = 291.7866150042338
        self.RotationY = 359.3666166267545
        self.RotationZ = 153.1795186816626
        self.ScaleX = 1.0
        self.ScaleY = 1.0
        self.ScaleZ = 1.0
        self.Trackball = 1
        self.TrackballQuaternion0 = 0.1344966115037936
        self.TrackballQuaternion1 = -0.5443772043350141
        self.TrackballQuaternion2 = -0.8061256037170713
        self.TrackballQuaternion3 = 0.1890122533757437
        self.TranslationX = 3.844376848011358
        self.TranslationY = 3.699782618579281
        self.VisibilityPositionX = 10
        self.VisibilityPositionY = 443
        self.VectorType = 1
        self.SmallAxes = 1
        self.__dict__.update(kwargs)

    def __str__(self):
        string = ''
        for k, v in self.__dict__.items():
            string += 'General.{0} = {1};\n'.format(k, v)
        return string

    def apply(self):
        from . import gmsh
        for k, v in self.__dict__.items():
            name = 'General.{0}'.format(k)
            gmsh.option.setNumber(name, v)

class Mesh(object):
    ''' Mesh Gmsh visualization options.
    For more information see http://gmsh.info/doc/texinfo/gmsh.html

    Attributes
    -----------
    AngleSmoothNormals: float
        Threshold angle below which normals are not smoothed
    SmoothNormals: int
        Smooth the mesh normals?
    SurfaceEdges: int
        Display edges of surface mesh?
    SurfaceFaces: int
        Display faces of surface mesh?
    VolumeEdges: int
        Display edges of volume mesh?
    VolumeFaces: int
        Display faces of volume mesh?
    Color: simnibs.gmsh_view.Color
        Color options

    References
    ------------
     `Gmsh documentation <http://gmsh.info/doc/texinfo/gmsh.html>`_

    '''
    def __init__(self, **kwargs):
        self.AngleSmoothNormals = 180  # Threshold angle below which normals are not smoothed
        self.SmoothNormals = 1  # Smooth the mesh normals?
        self.SurfaceEdges = 0  # Display edges of surface mesh?
        self.SurfaceFaces = 1  # Display faces of surface mesh?
        self.VolumeEdges = 0  # Display edges of volume mesh?
        self.VolumeFaces = 0  # Display faces of volume mesh?
        self.Color = Color()
        self.__dict__.update(kwargs)

    def __str__(self):
        string = ''
        for k, v in self.__dict__.items():
            if k == 'Color':
                c_st = str(v)
                st = c_st.replace('Color', 'Mesh.Color')
            else:
                st = 'Mesh.{0} = {1};\n'.format(k, v) 
            # st = st.replace('[', '{').replace(']', '}') 
            string += st

        return string 

    def apply(self):
        from . import gmsh
        for k, v in self.__dict__.items():
            name = 'Mesh.{0}'.format(k)
            gmsh.option.setNumber(name, v)


class Color(object):
    ''' Gmsh visualization Color options.
    For more information see http://gmsh.info/doc/texinfo/gmsh.html

    Attributes
    -----------
    One: list of ints
        Color 1 in color carousel
    Two: list of ints
        Color 2 in color carousel
    Three: list of ints
        Color 3 in color carousel
    Four: list of ints
        Color 4 in color carousel
    Five: list of ints
        Color 5 in color carousel


    References
    ------------
     `Gmsh documentation <http://gmsh.info/doc/texinfo/gmsh.html>`_

    '''

    def __init__(self, **kwargs):
        self.One = [230, 230, 230]  # Color 1 in color carousel
        self.Two = [129, 129, 129]  # Color 2 in color carousel
        self.Three = [104, 163, 255]  # Color 3 in color carousel
        self.Four = [255, 239, 179]  # Color 4 in color carousel
        self.Five = [255, 166, 133]  # Color 5 in color carousel
        self.__dict__.update(kwargs)

    def __str__(self):
        string = ''
        for k, v in self.__dict__.items():
            st = 'Color.{0} = {1};\n'.format(k, v)
            st = st.replace('[', '{').replace(']', '}')
            string += st
        return string

    def apply(self):
        from . import gmsh
        # NOT WORKING
        for k, v in self.__dict__.items():
            name = 'Mesh.Color.{0}'.format(k)
            string = str(v).replace('[', '{').replace(']', '}')
            gmsh.option.setString(name, string)


class View(object):
    ''' Gmsh visualization View options.
    For more information see http://gmsh.info/doc/texinfo/gmsh.html

    Parameters
    -----------
    indx: int (optional)
        Index of the view in the mesh

    Attributes
    -----------
    CenterGlyphs: int
        Center glyphs (arrows, numbers, etc.)? (0: left, 1: centered, 2: right)
    GlyphLocation: int
        Glyph (arrow, number, etc.) location (1: center of gravity, 2: node)
    VectorType: int
        Vector display type (1: segment, 2: arrow, 3: pyramid, 4: 3D arrow, 5:
        displacement, 6: comet)
    Visible: int
        Is the view visible?
    CustomMax: float
        User-defined maximum value to be displayed
    CustomMin: float
        User-defined minimum value to be displayed
    SaturateValues: int
        Saturate the view values to custom min and max (1: true, 0: false)
    ShowScale: int
        Show value scale?
    ColormapNumber: int
        Default colormap number (0: black, 1: vis5d, 2: jet, 3: lucie, 4: rainbow, 5:
        emc2000, 6: incadescent, 7: hot, 8: pink, 9: grayscale, 10: french, 11: hsv,
        12: spectrum, 13: bone, 14: spring, 15: summer, 16: autumm, 17: winter, 18:
        cool, 19: copper, 20: magma, 21: inferno, 22: plasma, 23: viridis)
    ColorTable: 255x3 list of ints
        Color table used to draw the view
    indx: int
        Index of the view in the mesh
    References
    ------------
     `Gmsh documentation <http://gmsh.info/doc/texinfo/gmsh.html>`_

    '''
    def __init__(self, indx=None, **kwargs):
        self.CenterGlyphs = 1
        self.GlyphLocation = 1
        self.VectorType = 1
        self.Visible = 0
        self.CustomMax = 0
        self.CustomMin = 0
        self.SaturateValues = 0
        self.RangeType = 1
        self.ShowScale = 1
        self.indx = indx
        self.ColormapNumber = 2
        self.ColorTable = None
        self.__dict__.update(kwargs)

    def _color_table_string(self):
        st = '{'
        st += ', '.join(
            ['{' + ', '.join([str(i) for i in c]) + '}'
             for c in self.ColorTable])
        st += '}'
        return st

    def __str__(self):
        string = ''
        exclude = ['indx', 'ColorTable']
        if self.indx is None:
            add = 'View.'
        else:
            add = 'View[%d].' % self.indx

        if self.ColorTable is not None:
            st = add + 'ColorTable = {0};\n'.format(
                        self._color_table_string())
            string += st
            exclude += ['ColormapNumber']

        for k, v in self.__dict__.items():
            if k not in exclude:
                st = add + '{0} = {1};\n'.format(k, v)
                string += st
        return string

    def apply(self):
        from . import gmsh
        for k, v in self.__dict__.items():
            if self.indx is None:
                name = 'View.{0}'.format(k)
            else:
                name = 'View[{0}].{1}'.format(self.indx, k)

            if k == 'ColorTable':
                if v is not None:
                    string = self._color_table_string()
                    gmsh.option.setString(name, string)

            elif k == 'indx':
                pass

            else:
                gmsh.option.setNumber(name, v)


def _run(command, remove=[]):
    try:
        subprocess.call(command)
    finally:
        [os.remove(r) for r in remove]


def _coolwarm_cm():
    # from http://www.kennethmoreland.com/color-maps/ 
    cm = np.array([
        [59,76,192],
        [60,78,194],
        [61,80,195],
        [62,81,197],
        [63,83,198],
        [64,85,200],
        [66,87,201],
        [67,88,203],
        [68,90,204],
        [69,92,206],
        [70,93,207],
        [71,95,209],
        [73,97,210],
        [74,99,211],
        [75,100,213],
        [76,102,214],
        [77,104,215],
        [79,105,217],
        [80,107,218],
        [81,109,219],
        [82,110,221],
        [84,112,222],
        [85,114,223],
        [86,115,224],
        [87,117,225],
        [89,119,226],
        [90,120,228],
        [91,122,229],
        [93,123,230],
        [94,125,231],
        [95,127,232],
        [96,128,233],
        [98,130,234],
        [99,131,235],
        [100,133,236],
        [102,135,237],
        [103,136,238],
        [104,138,239],
        [106,139,239],
        [107,141,240],
        [108,142,241],
        [110,144,242],
        [111,145,243],
        [112,147,243],
        [114,148,244],
        [115,150,245],
        [116,151,246],
        [118,153,246],
        [119,154,247],
        [120,156,247],
        [122,157,248],
        [123,158,249],
        [124,160,249],
        [126,161,250],
        [127,163,250],
        [129,164,251],
        [130,165,251],
        [131,167,252],
        [133,168,252],
        [134,169,252],
        [135,171,253],
        [137,172,253],
        [138,173,253],
        [140,174,254],
        [141,176,254],
        [142,177,254],
        [144,178,254],
        [145,179,254],
        [147,181,255],
        [148,182,255],
        [149,183,255],
        [151,184,255],
        [152,185,255],
        [153,186,255],
        [155,187,255],
        [156,188,255],
        [158,190,255],
        [159,191,255],
        [160,192,255],
        [162,193,255],
        [163,194,255],
        [164,195,254],
        [166,196,254],
        [167,197,254],
        [168,198,254],
        [170,199,253],
        [171,199,253],
        [172,200,253],
        [174,201,253],
        [175,202,252],
        [176,203,252],
        [178,204,251],
        [179,205,251],
        [180,205,251],
        [182,206,250],
        [183,207,250],
        [184,208,249],
        [185,208,248],
        [187,209,248],
        [188,210,247],
        [189,210,247],
        [190,211,246],
        [192,212,245],
        [193,212,245],
        [194,213,244],
        [195,213,243],
        [197,214,243],
        [198,214,242],
        [199,215,241],
        [200,215,240],
        [201,216,239],
        [203,216,238],
        [204,217,238],
        [205,217,237],
        [206,217,236],
        [207,218,235],
        [208,218,234],
        [209,219,233],
        [210,219,232],
        [211,219,231],
        [213,219,230],
        [214,220,229],
        [215,220,228],
        [216,220,227],
        [217,220,225],
        [218,220,224],
        [219,220,223],
        [220,221,222],
        [221,221,221],
        [222,220,219],
        [223,220,218],
        [224,219,216],
        [225,219,215],
        [226,218,214],
        [227,218,212],
        [228,217,211],
        [229,216,209],
        [230,216,208],
        [231,215,206],
        [232,215,205],
        [232,214,203],
        [233,213,202],
        [234,212,200],
        [235,212,199],
        [236,211,197],
        [236,210,196],
        [237,209,194],
        [238,209,193],
        [238,208,191],
        [239,207,190],
        [240,206,188],
        [240,205,187],
        [241,204,185],
        [241,203,184],
        [242,202,182],
        [242,201,181],
        [243,200,179],
        [243,199,178],
        [244,198,176],
        [244,197,174],
        [245,196,173],
        [245,195,171],
        [245,194,170],
        [245,193,168],
        [246,192,167],
        [246,191,165],
        [246,190,163],
        [246,188,162],
        [247,187,160],
        [247,186,159],
        [247,185,157],
        [247,184,156],
        [247,182,154],
        [247,181,152],
        [247,180,151],
        [247,178,149],
        [247,177,148],
        [247,176,146],
        [247,174,145],
        [247,173,143],
        [247,172,141],
        [247,170,140],
        [247,169,138],
        [247,167,137],
        [247,166,135],
        [246,164,134],
        [246,163,132],
        [246,161,131],
        [246,160,129],
        [245,158,127],
        [245,157,126],
        [245,155,124],
        [244,154,123],
        [244,152,121],
        [244,151,120],
        [243,149,118],
        [243,147,117],
        [242,146,115],
        [242,144,114],
        [241,142,112],
        [241,141,111],
        [240,139,109],
        [240,137,108],
        [239,136,106],
        [238,134,105],
        [238,132,103],
        [237,130,102],
        [236,129,100],
        [236,127,99],
        [235,125,97],
        [234,123,96],
        [233,121,95],
        [233,120,93],
        [232,118,92],
        [231,116,90],
        [230,114,89],
        [229,112,88],
        [228,110,86],
        [227,108,85],
        [227,106,83],
        [226,104,82],
        [225,102,81],
        [224,100,79],
        [223,98,78],
        [222,96,77],
        [221,94,75],
        [220,92,74],
        [218,90,73],
        [217,88,71],
        [216,86,70],
        [215,84,69],
        [214,82,67],
        [213,80,66],
        [212,78,65],
        [210,75,64],
        [209,73,62],
        [208,71,61],
        [207,69,60],
        [205,66,59],
        [204,64,57],
        [203,62,56],
        [202,59,55],
        [200,57,54],
        [199,54,53],
        [198,51,52],
        [196,49,50],
        [195,46,49],
        [193,43,48],
        [192,40,47],
        [190,37,46],
        [189,34,45],
        [188,30,44],
        [186,26,43],
        [185,22,41],
        [183,17,40],
        [181,11,39],
        [180,4,38]
    ], dtype=np.int)
    alpha = np.abs(np.arange(len(cm)) - len(cm)/2)**1.5
    alpha = (alpha * 255/alpha.max()).astype(int)
    cm = np.append(cm, alpha[:, None], axis=1)
    return cm


def _gray_red_lightblue_blue_cm():
   # simple colormap to visualize coils    
    cm = np.zeros((255,4),dtype = np.int)
    cm[0:64,0:3] = 132 # first quarter is gray 
    cm[64:128,0] = 255 # second quarter is red
    cm[128:192,2] = 255 # third quarter is lightblue
    cm[128:192,0:2] = 127
    cm[192:,2] = 255 # forth quarter is blue
    cm[0:128,3] = 60 # the coil is transparent
    cm[128:,3] = 114 # the logo is a bit less transparent
    return cm
