from setuptools import setup, find_packages, Extension
import os
import sys
import glob
import shutil
import urllib
import tempfile
import zipfile
from setuptools.command.build_ext import build_ext
from setuptools.command.develop import develop
from distutils.core import setup, Command
import numpy as np

''' C extensions

PETSc Linking
-------------

If you want to setup SimNIBS to use your own locally compiled PETSC, you can call
setup.py with the SIMNIBS_PETSC_INCLUDE and SIMNIBS_PETSC_DIR environment variables, for
example

export SIMNIBS_PETSC_DIR="/path/to/petsc/simnibs_petsc_arch/lib"
export SIMNIBS_PETSC_INCLUDE="/path/to/petsc/simnibs_direct_petsc_arch/include"
python setup.py develop

CGAL Compilation
-----------------

CGAL >= 5 is a header-only library, so we download it right before compiling.

Compilation requires:
GCC >= 6.3 or Apple Clang == 10.0.1 or MSVC >= 14.0

Boost >= 1.57

Boost can be instaled with
    Ubuntu: sudo apt install libboost-all-dev
    MacOS: brew install boost
    Windows: conda install boost

    Boost is also header-only, so we only need it during compile time

For more info, refer to https://doc.cgal.org/latest/Manual/thirdparty.html

'''

CGAL_version = '5.0'  # I tried 5.0.1 but tests fail!
CGAL_headers = os.path.abspath(f'CGAL-{CGAL_version}/include')

eigen_version = '3.3.7'
eigen_headers = os.path.abspath(f'eigen-{eigen_version}')

is_conda = 'CONDA_PREFIX' in os.environ

if sys.platform == 'win32':
    petsc_libs = ['libpetsc', 'msmpi']
    if 'SIMNIBS_PETSC_INCLUDE' in os.environ:
        petsc_include = [np.get_include()] + os.environ['SIMNIBS_PETSC_INCLUDE'].split(';')
    else:
        petsc_include = [np.get_include(),
                         'simnibs/external/include/win/petsc',
                         'simnibs/external/include/win/hypre',
                         'simnibs/external/include/win/mpi']
    if 'SIMNIBS_PETSC_DIR' in os.environ:
        petsc_dirs = [os.environ['SIMNIBS_PETSC_DIR']]
    else:
        petsc_dirs = ['simnibs/external/lib/win']
    petsc_runtime = None

    cgal_libs = ['libmpfr-4', 'libgmp-10', 'zlib']
    cgal_include = [
        np.get_include(),
        CGAL_headers,
        eigen_headers,
        'simnibs/external/include/win/mpfr',
        'simnibs/external/include/win/gmp'
    ]
    if is_conda:
        cgal_include += [os.path.join(os.environ['CONDA_PREFIX'], 'Library', 'include')]
    cgal_dirs = ['simnibs/external/lib/win']
    cgal_runtime = None
    cgal_compile_args = None

elif sys.platform == 'linux':
    petsc_libs = ['petsc']
    if 'SIMNIBS_PETSC_INCLUDE' in os.environ:
        petsc_include = [np.get_include()] + os.environ['SIMNIBS_PETSC_INCLUDE'].split(':')
    else:
        petsc_include = [np.get_include(), 'simnibs/external/include/linux/petsc']
    if 'SIMNIBS_PETSC_DIR' in os.environ:
        petsc_dirs = [os.environ['SIMNIBS_PETSC_DIR']]
        petsc_runtime = petsc_dirs
    else:
        petsc_dirs = ['simnibs/external/lib/linux']
        petsc_runtime = ['$ORIGIN/../external/lib/linux']

    cgal_libs = ['mpfr', 'gmp', 'z']
    cgal_include = [
        np.get_include(),
        CGAL_headers,
        eigen_headers,
        'simnibs/external/include/linux/mpfr',
        'simnibs/external/include/linux/gmp'
    ]
    if is_conda:
        cgal_include += [os.path.join(os.environ['CONDA_PREFIX'], 'include')]
    cgal_dirs = ['simnibs/external/lib/linux']
    cgal_runtime = ['$ORIGIN/../external/lib/linux']
    # TODO: setup compile args in other platforms
    cgal_compile_args = [
        '-Os', '-flto',
        '-DCGAL_CONCURRENT_MESH_3',
        '-DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS',
        '-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX',
        '-DNOMINMAX',
        '-DCGAL_EIGEN3_ENABLED',
        '-DCGAL_LINKED_WITH_TBB',
        '-frounding-math',
        '-std=gnu++14',

    ]
    # TODO: Add libtbb
    cgal_link_args = [
        '/usr/lib/x86_64-linux-gnu/libtbb.so',
        '/usr/lib/x86_64-linux-gnu/libtbbmalloc.so',
        '-lpthread'
    ]

elif sys.platform == 'darwin':
    petsc_libs = ['petsc']
    if 'SIMNIBS_PETSC_INCLUDE' in os.environ:
        petsc_include = [np.get_include()] + os.environ['SIMNIBS_PETSC_INCLUDE'].split(':')
    else:
        petsc_include = [np.get_include(), 'simnibs/external/include/osx/petsc']
    if 'SIMNIBS_PETSC_DIR' in os.environ:
        petsc_dirs = [os.environ['SIMNIBS_PETSC_DIR']]
    else:
        petsc_dirs = ['simnibs/external/lib/osx']
    petsc_runtime = None

    cgal_libs = ['mpfr', 'gmp', 'z']
    cgal_include = [
        np.get_include(),
        CGAL_headers,
        eigen_headers,
        'simnibs/external/include/osx/mpfr',
        'simnibs/external/include/osx/gmp'
    ]
    if is_conda:
        cgal_include += [os.path.join(os.environ['CONDA_PREFIX'], 'include')]
    cgal_dirs = ['simnibs/external/lib/osx']
    cgal_runtime = None
    cgal_compile_args = ['-march=native', '-stdlib=libc++', '-std=c++14']


else:
    raise OSError('OS not supported!')

extension = [
    Extension('simnibs.mesh_tools.cython_msh',
              ["simnibs/mesh_tools/cython_msh.pyx"],
              include_dirs=[np.get_include()]),
    Extension('simnibs.mesh_tools._marching_cubes_lewiner_cy',
              ["simnibs/mesh_tools/_marching_cubes_lewiner_cy.pyx"],
              include_dirs=[np.get_include()]),
    Extension('simnibs.segmentation._cs_utils',
              ["simnibs/segmentation/_cs_utils.pyx",
               "simnibs/external/cat12/genus0.c"],
              include_dirs=[np.get_include(),
                            'simnibs/external/cat12']),
    Extension('simnibs._compiled.petsc_solver',
              ["simnibs/_compiled/petsc_solver.pyx"],
              include_dirs=petsc_include,
              library_dirs=petsc_dirs,
              libraries=petsc_libs,
              runtime_library_dirs=petsc_runtime),
    Extension('simnibs._compiled.create_mesh',
              sources=["simnibs/_compiled/create_mesh.pyx"],
              language='c++',
              include_dirs=cgal_include,
              libraries=cgal_libs,
              library_dirs=cgal_dirs,
              runtime_library_dirs=cgal_runtime,
              extra_compile_args=cgal_compile_args,
              extra_link_args=cgal_link_args),
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


def download_and_extract(url):
    with urllib.request.urlopen(url) as response:
        with tempfile.NamedTemporaryFile('wb', delete=False) as tmpf:
            shutil.copyfileobj(response, tmpf)
            tmpname = tmpf.name

    with zipfile.ZipFile(tmpname) as z:
        z.extractall()

    os.remove(tmpname)

def move_libraries(build_folder, operation=shutil.move):
    if sys.platform == 'darwin':
        folder_name = 'osx'
    elif sys.platform == 'linux':
        folder_name = 'linux'
    elif sys.platform == 'win32':
        folder_name = 'win'

    lib_folder = os.path.join(build_folder, 'simnibs', 'external', 'lib', folder_name)
    compliled_folder = os.path.join(build_folder, 'simnibs', '_compiled')
    for fn in glob.glob(os.path.join(lib_folder, '*')):
        operation(
            fn,
            os.path.join(compliled_folder, os.path.basename(fn))
        )

class build_ext_(build_ext):
    def run(self):
        from Cython.Build import cythonize
        self.extension = cythonize(self.extensions)
        download_and_extract(
            f'https://github.com/CGAL/cgal/releases/download/'
            f'releases/CGAL-{CGAL_version}/'
            f'CGAL-{CGAL_version}-library.zip'
        )

        download_and_extract(
            f'https://gitlab.com/libeigen/eigen/-/'
            f'archive/{eigen_version}/eigen-{eigen_version}.zip'
        )
        build_ext.run(self)
        # Remove unescessary binary files
        shutil.rmtree(f'CGAL-{CGAL_version}', ignore_errors=True)
        shutil.rmtree(f'eigen-{eigen_version}', ignore_errors=True)
        linux_folders = [
            os.path.join(self.build_lib, 'simnibs', 'extenal', 'bin', 'linux'),
            os.path.join(self.build_lib, 'simnibs', 'extenal', 'include', 'linux'),
            os.path.join(self.build_lib, 'simnibs', 'extenal', 'lib', 'linux'),
            os.path.join(self.build_lib, 'simnibs', 'resources',
                         'spm12', 'toolbox', 'cat12', 'CAT.glnx86'),
        ]
        osx_folders = [
            os.path.join(self.build_lib, 'simnibs', 'extenal', 'bin', 'osx'),
            os.path.join(self.build_lib, 'simnibs', 'extenal', 'include', 'osx'),
            os.path.join(self.build_lib, 'simnibs', 'extenal', 'lib', 'osx'),
            os.path.join(self.build_lib, 'simnibs', 'resources',
                         'spm12', 'toolbox', 'cat12', 'CAT.maci64'),
        ]
        win_folders = [
            os.path.join(self.build_lib, 'simnibs', 'extenal', 'bin', 'win'),
            os.path.join(self.build_lib, 'simnibs', 'extenal', 'include', 'win'),
            os.path.join(self.build_lib, 'simnibs', 'extenal', 'lib', 'win'),
            os.path.join(self.build_lib, 'simnibs', 'resources',
                         'spm12', 'toolbox', 'cat12', 'CAT.w32'),
        ]
        if sys.platform == 'linux':
            [shutil.rmtree(f, True) for f in osx_folders]
            [shutil.rmtree(f, True) for f in win_folders]

        if sys.platform == 'darwin':
            [shutil.rmtree(f, True) for f in linux_folders]
            [shutil.rmtree(f, True) for f in win_folders]
            move_libraries(self.build_lib)

        if sys.platform == 'win32':
            [shutil.rmtree(f, True) for f in linux_folders]
            [shutil.rmtree(f, True) for f in osx_folders]
            move_libraries(self.build_lib)


class develop_(develop):
    def run(self):
        develop.run(self)
        if sys.platform == 'win32':
            move_libraries('.', shutil.copy)
        if sys.platform == 'darwin':
            move_libraries('.', shutil.copy)

setup(name='simnibs',
      version=open("simnibs/_version.py").readlines()[-1].split()[-1].strip("\"'"),
      description='www.simnibs.org',
      author='SimNIBS developers',
      author_email='support@simnibs.org',
      packages=find_packages(),
      license='GPL3',
      ext_modules=extension,
      include_package_data=True,
      cmdclass={
          'build_ext': build_ext_,
          'develop': develop_
        },
      entry_points={
          'console_scripts': console_scripts,
          'gui_scripts': gui_scripts
      },
      install_requires=[
          'numpy>=1.16',
          'scipy>=1.2',
          'h5py>=2.9',
          'nibabel>=2.3'
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
