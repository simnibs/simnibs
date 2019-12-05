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

CGAL_version = '5.0'
CGAL_headers = os.path.abspath(f'CGAL-{CGAL_version}/include')

if sys.platform == 'win32':
    petsc_libs = ['libpetsc', 'msmpi']
    if 'SIMNIBS_PETSC_INCLUDE' in os.environ:
        petsc_include = [np.get_include()] + os.environ['SIMNIBS_PETSC_INCLUDE'].split(';')
    else:
        petsc_include = [np.get_include(),
                         'simnibs/include/win/petsc',
                         'simnibs/include/win/hypre',
                         'simnibs/include/win/mpi']
    if 'SIMNIBS_PETSC_DIR' in os.environ:
        petsc_dirs = [os.environ['SIMNIBS_PETSC_DIR']]
    else:
        petsc_dirs = ['simnibs/lib/win']
    petsc_runtime = None

    cgal_libs = ['libmpfr-4', 'libgmp-10', 'zlib']
    cgal_include = [
        np.get_include(),
        CGAL_headers,
        'simnibs/include/win/mpfr',
        'simnibs/include/win/gmp',
        # Assuming conda install boost
        os.path.join(os.environ['CONDA_PREFIX'], 'Library', 'include')
    ]
    cgal_dirs = ['simnibs/lib/win']
    cgal_runtime = None
    cgal_compile_args = None

elif sys.platform == 'linux':
    petsc_libs = ['petsc']
    if 'SIMNIBS_PETSC_INCLUDE' in os.environ:
        petsc_include = [np.get_include()] + os.environ['SIMNIBS_PETSC_INCLUDE'].split(':')
    else:
        petsc_include = [np.get_include(), 'simnibs/include/linux/petsc']
    if 'SIMNIBS_PETSC_DIR' in os.environ:
        petsc_dirs = [os.environ['SIMNIBS_PETSC_DIR']]
        petsc_runtime = petsc_dirs
    else:
        petsc_dirs = ['simnibs/lib/linux']
        petsc_runtime = ['$ORIGIN/../lib/linux']

    cgal_libs = ['mpfr', 'gmp', 'z']
    cgal_include = [
        np.get_include(),
        CGAL_headers,
        'simnibs/include/linux/mpfr',
        'simnibs/include/linux/gmp',
        # Assuming conda install boost
        os.path.join(os.environ['CONDA_PREFIX'], 'include')
    ]
    cgal_dirs = ['simnibs/lib/linux']
    cgal_runtime = ['$ORIGIN/../lib/linux']
    cgal_compile_args = ['-Os', '-flto']

elif sys.platform == 'darwin':
    petsc_libs = ['petsc']
    if 'SIMNIBS_PETSC_INCLUDE' in os.environ:
        petsc_include = [np.get_include()] + os.environ['SIMNIBS_PETSC_INCLUDE'].split(':')
    else:
        petsc_include = [np.get_include(), 'simnibs/include/osx/petsc']
    if 'SIMNIBS_PETSC_DIR' in os.environ:
        petsc_dirs = [os.environ['SIMNIBS_PETSC_DIR']]
    else:
        petsc_dirs = ['simnibs/lib/osx']
    petsc_runtime = None

    cgal_libs = ['mpfr', 'gmp', 'z']
    cgal_include = [
        np.get_include(),
        CGAL_headers,
        'simnibs/include/osx/mpfr',
        'simnibs/include/osx/gmp',
    ]
    cgal_dirs = ['simnibs/lib/osx']
    cgal_runtime = None
    cgal_compile_args = None


else:
    raise OSError('OS not supported!')

extension = [
    Extension('simnibs.msh.cython_msh',
              ["simnibs/msh/cython_msh.pyx"],
              include_dirs=[np.get_include()]),
    Extension('simnibs.msh._marching_cubes_lewiner_cy',
              ["simnibs/msh/_marching_cubes_lewiner_cy.pyx"],
              include_dirs=[np.get_include()]),
    Extension('simnibs.segmentation._cs_utils',
              ["simnibs/segmentation/_cs_utils.pyx",
               "simnibs/segmentation/genus0.c"],
              include_dirs=[np.get_include()]),
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
              extra_link_args=cgal_compile_args),
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


def download_cgal():
    CGAL_url = (
        'https://github.com/CGAL/cgal/releases/download/'
        'releases/CGAL-{0}/CGAL-{0}-library.zip'.format(CGAL_version)
    )
    with urllib.request.urlopen(CGAL_url) as response:
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

    lib_folder = os.path.join(build_folder, 'simnibs', 'lib', folder_name)
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
        download_cgal()
        build_ext.run(self)
        # Remove unescessary binary files
        shutil.rmtree(f'CGAL-{CGAL_version}', ignore_errors=True)
        linux_folders = [
            os.path.join(self.build_lib, 'simnibs', 'bin', 'linux'),
            os.path.join(self.build_lib, 'simnibs', 'include', 'linux'),
            os.path.join(self.build_lib, 'simnibs', 'lib', 'linux'),
            os.path.join(self.build_lib, 'simnibs', 'resources',
                         'spm12', 'toolbox', 'cat12', 'CAT.glnx86'),
        ]
        osx_folders = [
            os.path.join(self.build_lib, 'simnibs', 'bin', 'osx'),
            os.path.join(self.build_lib, 'simnibs', 'include', 'osx'),
            os.path.join(self.build_lib, 'simnibs', 'lib', 'osx'),
            os.path.join(self.build_lib, 'simnibs', 'resources',
                         'spm12', 'toolbox', 'cat12', 'CAT.maci64'),
        ]
        win_folders = [
            os.path.join(self.build_lib, 'simnibs', 'bin', 'win'),
            os.path.join(self.build_lib, 'simnibs', 'include', 'win'),
            os.path.join(self.build_lib, 'simnibs', 'lib', 'win'),
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
