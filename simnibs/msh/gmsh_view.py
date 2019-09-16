import os
import tempfile
import subprocess
import threading
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
        self.View.append(View(len(self.View), **kwargs))

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
        self.CenterGlyphs = 0
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
            st = add + 'ColorTable = {1};\n'.format(
                        self._color_table_string())
            string += st

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
