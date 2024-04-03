from setuptools import setup, Extension
import os
import sys
import shutil
from setuptools.command.build_ext import build_ext
from distutils.dep_util import newer_group
import numpy as np


''' C extensions

CGAL Compilation
-----------------

CGAL >= 5 is a header-only library, so we download it right before compiling.

Compilation requires:
GCC >= 6.3 or Apple Clang == 10.0.1 or MSVC >= 14.0
    conda install gcc_linux-64 gxx_linux-64 gfortran_linux-64
Boost >= 1.57

Boost can be instaled with
    Ubuntu: sudo apt install libboost-all-dev
    MacOS: brew install boost
    Windows: conda install boost
    Boost is also header-only, so we only need it during compile time

For more info, refer to https://doc.cgal.org/latest/Manual/thirdparty.html

'''

cgal_mesh_macros = [
    ('CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX', None),
    ('CGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS', None),
    ('CGAL_CONCURRENT_MESH_3', None),
    ('CGAL_EIGEN3_ENABLED', None),
    ('CGAL_USE_ZLIB', 1),
    ('CGAL_LINKED_WITH_TBB', None)
]

is_conda = 'CONDA_PREFIX' in os.environ
# No conda, no setup
if not is_conda:
    raise Exception("Cannot run setup without conda")

# Information for eigen library
eigen_headers = os.path.abspath(os.path.join(os.environ['CONDA_PREFIX'], 'include','eigen3'))

#### Setup compilation arguments

if sys.platform == 'win32':
    petsc_libs = ['libpetsc', 'msmpi']
    petsc_include = [
        np.get_include(),
        'simnibs/external/include/win/petsc',
        'simnibs/external/include/win/hypre',
        'simnibs/external/include/win/mpi'
    ]
    petsc_dirs = ['simnibs/external/lib/win']
    petsc_runtime = None
    petsc_extra_link_args = None

    cgal_libs = ['mpfr', 'gmp', 'zlib', 'tbb', 'tbbmalloc']
    cgal_dirs = [os.path.join(os.environ['CONDA_PREFIX'], 'Library', 'lib')]
    cgal_include = [os.path.join(os.environ['CONDA_PREFIX'], 'Library', 'include')]
    cgal_include += [
        np.get_include(),
        eigen_headers,
    ]
    # Find boost headers if installed with conda
    cgal_runtime = None
    # Got those arguments from compiling a CGAL program following the instructions in the website
    cgal_compile_args = [
        '/Zi', '/WX-', '/diagnostics:classic', '/Ob0', '/Oy',
        '/D WIN32', '/D _WINDOWS', '/D _SCL_SECURE_NO_DEPRECATE',
        '/D _SCL_SECURE_NO_WARNINGS', '/D BOOST_ALL_DYN_LINK=1',
        '/D _MBCS'
    ]
    cgal_link_args = None

    cat_compile_args = None

elif sys.platform == 'linux':
    petsc_libs = ['petsc']
    petsc_include = [
        np.get_include(),
        'simnibs/external/include/linux/petsc'
    ]
    petsc_dirs = ['simnibs/external/lib/linux']
    petsc_runtime = ['$ORIGIN/../external/lib/linux']
    petsc_extra_link_args = None

    cgal_libs = ['mpfr', 'gmp', 'z', 'tbb', 'tbbmalloc', 'pthread']
    # To find the boost headers if installed with conda
    cgal_include = []
    if is_conda:
        cgal_include += [os.path.join(os.environ['CONDA_PREFIX'], 'include')]

    cgal_include += [
        np.get_include(),
        eigen_headers,
    ]

    cgal_dirs = [os.path.join(os.environ['CONDA_PREFIX'], 'Library', 'lib')]
    cgal_runtime = ['$ORIGIN/../../external/lib/linux']
    # Add -Os -flto for much smaller binaries
    cgal_compile_args = [
        '-Os', '-flto',
        '-frounding-math',
        '-std=gnu++14',
    ]
    cgal_mesh_macros += [('NOMINMAX', None)]
    cgal_link_args = None

    cat_compile_args = [
      '-std=gnu99',
    ]

elif sys.platform == 'darwin':
    petsc_libs = ['petsc']
    petsc_include = [
        np.get_include(),
        'simnibs/external/include/osx/petsc'
    ]
    petsc_dirs = ['simnibs/external/lib/osx']
    petsc_runtime = None
    # add RPATH as the _runtime argument does not work in MacOS, likely bug in setuptools
    petsc_extra_link_args = ['-Wl,-rpath,@loader_path/../external/lib/osx']

    cgal_libs = ['mpfr', 'gmp', 'z', 'tbb', 'tbbmalloc']
    cgal_dirs = [os.path.join(os.environ['CONDA_PREFIX'], 'Library', 'lib')]
    cgal_include = [os.path.join(os.environ['CONDA_PREFIX'], 'include')]
    cgal_include += [
        np.get_include(),
        eigen_headers,
    ]

    cgal_runtime = None
    cgal_compile_args = [
        '-std=gnu++14',
        '-stdlib=libc++',
    ]
    cgal_mesh_macros += [('NOMINMAX', None)]
    cgal_link_args = [
        '-stdlib=libc++',
        '-Wl,-rpath,@loader_path/../../external/lib/osx'
    ]

    cat_compile_args = None

else:
    raise OSError('OS not supported!')

cython_msh = Extension(
    'simnibs.mesh_tools.cython_msh',
    ["simnibs/mesh_tools/cython_msh.pyx"],
    include_dirs=[np.get_include()]
)
marching_cubes_lewiner_cy = Extension(
    'simnibs.segmentation._marching_cubes_lewiner_cy',
    ["simnibs/segmentation/_marching_cubes_lewiner_cy.pyx"],
    include_dirs=[np.get_include()]
)
cat_c_utils = Extension(
    'simnibs.segmentation._cat_c_utils',
    ["simnibs/segmentation/_cat_c_utils.pyx", "simnibs/segmentation/cat_c_utils/genus0.c"],
    include_dirs=[np.get_include(), 'simnibs/segmentation/cat_c_utils'],
    extra_compile_args=cat_compile_args
)
thickness = Extension(
    'simnibs.segmentation._thickness',
    ["simnibs/segmentation/_thickness.pyx"],
    include_dirs=[np.get_include()]
)
petsc_solver = Extension(
    'simnibs.simulation.petsc_solver',
    sources=["simnibs/simulation/petsc_solver.pyx"],
    depends=["simnibs/simulation/_solver.c"],
    include_dirs=petsc_include,
    library_dirs=petsc_dirs,
    libraries=petsc_libs,
    runtime_library_dirs=petsc_runtime,
    extra_link_args=petsc_extra_link_args
)
# I separated the CGAL functions into several files for two reasons
# 1. Reduce memory consumption during compilation in Linux
# 2. Fix some compilation problems in Windows
create_mesh_surf = Extension(
    'simnibs.mesh_tools.cgal.create_mesh_surf',
    sources=["simnibs/mesh_tools/cgal/create_mesh_surf.pyx"],
    depends=["simnibs/mesh_tools/cgal/_mesh_surfaces.cpp"],
    language='c++',
    include_dirs=cgal_include,
    libraries=cgal_libs,
    library_dirs=cgal_dirs,
    runtime_library_dirs=cgal_runtime,
    extra_compile_args=cgal_compile_args,
    extra_link_args=cgal_link_args,
    define_macros=cgal_mesh_macros
)
create_mesh_vol = Extension(
    'simnibs.mesh_tools.cgal.create_mesh_vol',
    sources=["simnibs/mesh_tools/cgal/create_mesh_vol.pyx"],
    depends=["simnibs/mesh_tools/cgal/_mesh_volumes.cpp"],
    language='c++',
    include_dirs=cgal_include,
    libraries=cgal_libs,
    library_dirs=cgal_dirs,
    runtime_library_dirs=cgal_runtime,
    extra_compile_args=cgal_compile_args,
    extra_link_args=cgal_link_args,
    define_macros=cgal_mesh_macros
)
cgal_misc = Extension(
    'simnibs.mesh_tools.cgal.cgal_misc',
    sources=["simnibs/mesh_tools/cgal/cgal_misc.pyx"],
    depends=["simnibs/mesh_tools/cgal/_cgal_intersect.cpp"],
    language='c++',
    include_dirs=cgal_include,
    libraries=cgal_libs,
    library_dirs=cgal_dirs,
    runtime_library_dirs=cgal_runtime,
    extra_compile_args=cgal_compile_args,
    extra_link_args=cgal_link_args,
)
cgal_pmp = Extension(
    "simnibs.mesh_tools.cgal.polygon_mesh_processing",
    sources = ["simnibs/mesh_tools/cgal/polygon_mesh_processing.pyx"],
    depends = ["simnibs/mesh_tools/cgal/polygon_mesh_processing_src.cpp"],
    language='c++',
    include_dirs=cgal_include,
    libraries=cgal_libs,
    library_dirs=cgal_dirs,
    runtime_library_dirs=cgal_runtime,
    extra_compile_args=cgal_compile_args,
    extra_link_args=cgal_link_args,
)

extensions = [
    cython_msh,
    marching_cubes_lewiner_cy,
    cat_c_utils,
    thickness,
    petsc_solver,
    create_mesh_surf,
    create_mesh_vol,
    cgal_misc,
    cgal_pmp,
]


class build_ext_(build_ext):
    '''
        Build the extension, download some dependencies and remove stuff from other OS
    '''
    def run(self):
        from Cython.Build import cythonize
        ## Cythonize
        self.extension = cythonize(self.extensions)
        ## Download requirements
        changed_meshing = (
            newer_group(
                create_mesh_surf.sources + create_mesh_surf.depends,
                self.get_ext_fullpath(create_mesh_surf.name),
                'newer'
            ) or
            newer_group(
                create_mesh_vol.sources + create_mesh_vol.depends,
                self.get_ext_fullpath(create_mesh_vol.name),
                'newer'
            ) or
            newer_group(
                cgal_misc.sources + cgal_misc.depends,
                self.get_ext_fullpath(cgal_misc.name),
                'newer'
            )
        )

        # Compile
        build_ext.run(self)

        # Remove unescessary binary files
        linux_folders = [
            os.path.join(self.build_lib, 'simnibs', 'external', 'bin', 'linux'),
            os.path.join(self.build_lib, 'simnibs', 'external', 'include', 'linux'),
            os.path.join(self.build_lib, 'simnibs', 'external', 'lib', 'linux'),
        ]
        osx_folders = [
            os.path.join(self.build_lib, 'simnibs', 'external', 'bin', 'osx'),
            os.path.join(self.build_lib, 'simnibs', 'external', 'include', 'osx'),
            os.path.join(self.build_lib, 'simnibs', 'external', 'lib', 'osx'),
        ]
        win_folders = [
            os.path.join(self.build_lib, 'simnibs', 'external', 'bin', 'win'),
            os.path.join(self.build_lib, 'simnibs', 'external', 'include', 'win'),
            os.path.join(self.build_lib, 'simnibs', 'external', 'lib', 'win'),
        ]
        if sys.platform == 'linux':
            [shutil.rmtree(f, True) for f in osx_folders]
            [shutil.rmtree(f, True) for f in win_folders]

        if sys.platform == 'darwin':
            [shutil.rmtree(f, True) for f in linux_folders]
            [shutil.rmtree(f, True) for f in win_folders]

        if sys.platform == 'win32':
            [shutil.rmtree(f, True) for f in linux_folders]
            [shutil.rmtree(f, True) for f in osx_folders]


setup(ext_modules=extensions,
      cmdclass={
          'build_ext': build_ext_
          },
      )
