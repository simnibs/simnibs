from setuptools import setup, find_packages, Extension
import os
import sys
import glob
import shutil
import urllib
import tempfile
import zipfile
import tarfile
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

# Information for CGAL download
CGAL_version = '5.0.2'
CGAL_headers = os.path.abspath(f'CGAL-{CGAL_version}/include')
CGAL_url = (
    f'https://github.com/CGAL/cgal/releases/download/'
    f'releases/CGAL-{CGAL_version}/'
    f'CGAL-{CGAL_version}-library.zip'
)
cgal_mesh_macros = [
    ('CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX', None),
    ('CGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS', None),
    ('CGAL_CONCURRENT_MESH_3', None),
    ('CGAL_EIGEN3_ENABLED', None),
    ('CGAL_USE_ZLIB', 1),
    ('CGAL_LINKED_WITH_TBB', None)
]

# Information for eigen library
# I don't download it because gitlab does not allow it
eigen_version = '3.3.7'
eigen_headers = os.path.abspath(f'simnibs/external/include/eigen-{eigen_version}')

# Information for Intel TBB download
tbb_version = '2020.1'
tbb_path = os.path.abspath('tbb')
tbb_headers = os.path.join(tbb_path, 'tbb', 'include')
if sys.platform == 'win32':
    tbb_url = (
        f'https://github.com/intel/tbb/releases/download/'
        f'v{tbb_version}/tbb-{tbb_version}-win.zip'
    )
    tbb_libs = [
        os.path.join(tbb_path, 'tbb', 'bin', 'intel64', 'vc14', 'tbb.dll'),
        os.path.join(tbb_path, 'tbb', 'lib', 'intel64', 'vc14', 'tbb.lib'),
        os.path.join(tbb_path, 'tbb', 'bin', 'intel64', 'vc14', 'tbbmalloc.dll'),
        os.path.join(tbb_path, 'tbb', 'lib', 'intel64', 'vc14', 'tbbmalloc.lib'),
    ]
elif sys.platform == 'linux':
    tbb_url = (
        f'https://github.com/intel/tbb/releases/download/'
        f'v{tbb_version}/tbb-{tbb_version}-lin.tgz'
    )
    tbb_libs = [
        os.path.join(tbb_path, 'tbb', 'lib', 'intel64', 'gcc4.8', 'libtbb.so'),
        os.path.join(tbb_path, 'tbb', 'lib', 'intel64', 'gcc4.8', 'libtbb.so.2'),
        os.path.join(tbb_path, 'tbb', 'lib', 'intel64', 'gcc4.8', 'libtbbmalloc.so'),
        os.path.join(tbb_path, 'tbb', 'lib', 'intel64', 'gcc4.8', 'libtbbmalloc.so.2'),
    ]
elif sys.platform == 'darwin':
    tbb_url = (
        f'https://github.com/intel/tbb/releases/download/'
        f'v{tbb_version}/tbb-{tbb_version}-mac.tgz'
    )
    tbb_libs = [
        os.path.join(tbb_path, 'tbb', 'lib', 'libtbb.dylib'),
        os.path.join(tbb_path, 'tbb', 'lib', 'libtbbmalloc.dylib'),
    ]
else:
    raise OSError('OS not supported!')


#### Setup compilation arguments
is_conda = 'CONDA_PREFIX' in os.environ

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

    cgal_libs = ['libmpfr-4', 'libgmp-10', 'zlib', 'tbb', 'tbbmalloc']
    cgal_include = [
        np.get_include(),
        CGAL_headers,
        eigen_headers,
        tbb_headers,
        'simnibs/external/include/win/mpfr',
        'simnibs/external/include/win/gmp'
    ]
    # Find boost headers if installed with conda
    if is_conda:
        cgal_include += [os.path.join(os.environ['CONDA_PREFIX'], 'Library', 'include')]
    cgal_dirs = ['simnibs/external/lib/win']
    cgal_runtime = None
    # Got those arguments from compiling a CGAL program following the instructions in the website
    cgal_compile_args = [
        '/Zi', '/WX-', '/diagnostics:classic', '/Ob0', '/Oy',
        '/D WIN32', '/D _WINDOWS', '/D _SCL_SECURE_NO_DEPRECATE',
        '/D _SCL_SECURE_NO_WARNINGS', '/D BOOST_ALL_DYN_LINK=1',
        '/D _MBCS'
    ]
    cgal_link_args = None

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
    cgal_include = [
        np.get_include(),
        CGAL_headers,
        eigen_headers,
        tbb_headers,
        'simnibs/external/include/linux/mpfr',
        'simnibs/external/include/linux/gmp'
    ]
    # To find the boost headers if installed with conda
    if is_conda:
        cgal_include += [os.path.join(os.environ['CONDA_PREFIX'], 'include')]
    cgal_dirs = ['simnibs/external/lib/linux']
    cgal_runtime = ['$ORIGIN/../../external/lib/linux']
    # Add -Os -flto for much smaller binaries
    cgal_compile_args = [
        '-Os', '-flto',
        '-frounding-math',
        '-std=gnu++14',
    ]
    cgal_mesh_macros += [('NOMINMAX', None)]
    cgal_link_args = None

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
    cgal_include = [
        np.get_include(),
        CGAL_headers,
        eigen_headers,
        tbb_headers,
        'simnibs/external/include/osx/mpfr',
        'simnibs/external/include/osx/gmp'
    ]
    if is_conda:
        cgal_include += [os.path.join(os.environ['CONDA_PREFIX'], 'include')]
    cgal_dirs = ['simnibs/external/lib/osx']
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
    include_dirs=[np.get_include(), 'simnibs/segmentation/cat_c_utils']
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


extensions = [
    cython_msh,
    marching_cubes_lewiner_cy,
    cat_c_utils,
    thickness,
    petsc_solver,
    create_mesh_surf,
    create_mesh_vol,
    cgal_misc
]
####################################################
# Add all scripts in the cli folder and a bit more #
####################################################

# IMPORTANT: For the postinstall script to also work
# ALL scripts should be in the simnibs/cli folder and have
# a if __name__ == '__main__' clause

script_names = [os.path.splitext(os.path.basename(s))[0]
                for s in glob.glob('simnibs/cli/*.py')]

console_scripts = []

for s in script_names:
    if s not in ['__init__', 'simnibs_gui']:
        console_scripts.append(f'{s}=simnibs.cli.{s}:main')

console_scripts.append(f'simnibs=simnibs.cli.run_simnibs:main')

gui_scripts = [
    'simnibs_gui=simnibs.cli.simnibs_gui:main',
]


def download_and_extract(url, path='.'):
    ''' Downloads and extracts a zip or tar-gz folder '''
    print('Downloading:', url)
    with urllib.request.urlopen(url) as response:
        with tempfile.NamedTemporaryFile('wb', delete=False) as tmpf:
            shutil.copyfileobj(response, tmpf)
            tmpname = tmpf.name

    if url.endswith('.zip'):
        with zipfile.ZipFile(tmpname) as z:
            z.extractall(path)

    elif url.endswith('.tgz') or url.endswith('.tar.gz'):
        with tarfile.open(tmpname, 'r:gz') as z:
            z.extractall(path)
    else:
        raise IOError('Could not extract file, unrecognized extension')

    os.remove(tmpname)


def install_lib(url, path, libs):
    ''' Downloads a compiled library from the internet and move to "lib" folder '''
    download_and_extract(url, path)
    if sys.platform == 'darwin':
        folder_name = 'osx'
    elif sys.platform == 'linux':
        folder_name = 'linux'
    elif sys.platform == 'win32':
        folder_name = 'win'
    for l in libs:
        shutil.copy(
            l, f'simnibs/external/lib/{folder_name}',
            follow_symlinks=False
        )

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
        if self.force or changed_meshing:
            download_and_extract(CGAL_url)
            install_lib(tbb_url, tbb_path, tbb_libs)
        # Compile
        build_ext.run(self)
        # cleanup downloads
        if self.force or changed_meshing:
            shutil.rmtree(f'CGAL-{CGAL_version}', ignore_errors=True)
            shutil.rmtree(tbb_path, ignore_errors=True)
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

setup(name='simnibs',
      version=open("simnibs/_version.py").readlines()[-1].split()[-1].strip("\"'"),
      description='www.simnibs.org',
      author='SimNIBS developers',
      author_email='support@simnibs.org',
      packages=find_packages(),
      license='GPL3',
      ext_modules=extensions,
      include_package_data=True,
      cmdclass={
          'build_ext': build_ext_,
        },
      entry_points={
          'console_scripts': console_scripts,
          'gui_scripts': gui_scripts
      },
      install_requires=[
          'numpy>=1.16',
          'scipy>=1.2',
          'h5py>=2.9',
          'nibabel>=2.3',
          'packaging',
          'requests',
          'charm-gems'
      ],
      extras_require={
          'GUI': ['pyqt5', 'pyopengl']
      },
      setup_requires=[
          'numpy>=1.16',
          'cython'
      ],
      tests_require=['pytest', 'mock'],
      zip_safe=False)
