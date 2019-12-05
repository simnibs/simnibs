from setuptools import setup, find_packages, Extension
import os
import sys
import glob
import shutil
from setuptools.command.build_ext import build_ext
from setuptools.command.develop import develop
from distutils.core import setup, Command
import numpy as np

########################
# Compile C extensions #
########################
'''
If you want to setup SimNIBS to use your own locally compiled PETSC, you can call
setup.py with the SIMNIBS_PETSC_INCLUDE and SIMNIBS_PETSC_DIR environment variables, for
example

SIMNIBS_PETSC_DIR=/path/to/petsc/simnibs_petsc_arch/lib \
SIMNIBS_PETSC_INCLUDE="/path/to/petsc/simnibs_direct_petsc_arch/include:/path/to/petsc/include" \
python setup.py develop
'''
if sys.platform == 'win32':
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
    petsc_libs = ['libpetsc', 'msmpi']
    petsc_runtime = None

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

else:
    raise OSError('OS not supported!')

extension = [
    Extension('simnibs.msh.cython_msh',
              ["simnibs/msh/cython_msh.pyx"],
              include_dirs=[np.get_include()]),
    Extension('simnibs.msh._marching_cubes_lewiner_cy',
              ["simnibs/msh/_marching_cubes_lewiner_cy.pyx"],
              include_dirs=[np.get_include()]),
    Extension('simnibs._compiled.petsc_solver',
              ["simnibs/_compiled/petsc_solver.pyx"],
              include_dirs=petsc_include,
              library_dirs=petsc_dirs,
              libraries=petsc_libs,
              runtime_library_dirs=petsc_runtime)
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
        build_ext.run(self)
        # Remove unescessary binary files
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
