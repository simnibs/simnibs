from setuptools import setup, find_packages, Extension
import os
import sys
import glob
import shutil
from setuptools.command.build_ext import build_ext
from setuptools.command.develop import develop
import numpy as np

########################
# Compile C extensions #
########################
try:
    from Cython.Build import cythonize
    ext = '.pyx'
except ImportError:
    ext = '.c'

if sys.platform == 'win32':
    petsc_include = [np.get_include(),
                     'simnibs/include/win/petsc',
                     'simnibs/include/win/hypre',
                     'simnibs/include/win/mpi']
    petsc_libs = ['libpetsc', 'msmpi']
    petsc_dirs = ['simnibs/lib/win']
    petsc_runtime = None

elif sys.platform == 'linux':
    petsc_include = [np.get_include(), 'simnibs/include/linux/petsc']
    petsc_libs = ['petsc']
    petsc_dirs = ['simnibs/lib/linux']
    petsc_runtime = ['$ORIGIN/../lib/linux']

elif sys.platform == 'darwin':
    petsc_include = [np.get_include(), 'simnibs/include/osx/petsc']
    petsc_libs = ['petsc']
    petsc_dirs = ['simnibs/lib/osx']
    petsc_runtime = None

else:
    raise OSError('OS not supported!')

extension = [
    Extension('simnibs.cython_code.cython_msh',
              ["simnibs/cython_code/cython_msh" + ext],
              include_dirs=[np.get_include()]),
    Extension('simnibs.cython_code._marching_cubes_lewiner_cy',
              ["simnibs/cython_code/_marching_cubes_lewiner_cy" + ext],
              include_dirs=[np.get_include()]),
    Extension('simnibs.cython_code.petsc_solver',
              ["simnibs/cython_code/petsc_solver" + ext],
              include_dirs=petsc_include,
              library_dirs=petsc_dirs,
              libraries=petsc_libs,
              runtime_library_dirs=petsc_runtime,
    )]

if ext == '.pyx':
    extension = cythonize(extension)


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
    if s not in ['__init__', 'gmsh', 'simnibs_gui']:
        console_scripts.append(f'{s}=simnibs.cli.{s}:main')

console_scripts.append(f'simnibs=simnibs.cli.run_simnibs:main')

gui_scripts = [
    'simnibs_gui=simnibs.cli.simnibs_gui:main',
    'gmsh=simnibs.cli.gmsh:main',
]


class build_ext_(build_ext):
    def run(self):
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
            for f in glob.glob(os.path.join(self.build_lib, 'simnibs', 'lib', 'osx', '*')):
                shutil.move(
                            f,
                            os.path.join(self.build_lib, 'simnibs', 'cython_code', os.path.basename(f)))
        
        if sys.platform == 'win32':
            [shutil.rmtree(f, True) for f in linux_folders]
            [shutil.rmtree(f, True) for f in osx_folders]
            for f in glob.glob(os.path.join(self.build_lib, 'simnibs', 'lib', 'win', '*')):
                shutil.move(
                    f,
                    os.path.join(self.build_lib, 'simnibs', 'cython_code', os.path.basename(f)))


class develop_(develop):
    def run(self):
        develop.run(self)
        if sys.platform == 'win32':
            for f in glob.glob(os.path.join('simnibs', 'lib', 'win', '*')):
                shutil.copy(
                    f,
                    os.path.join('simnibs', 'cython_code', os.path.basename(f)))
        if sys.platform == 'darwin':
            for f in glob.glob(os.path.join('simnibs', 'lib', 'osx', '*')):
                shutil.copy(
                    f,
                    os.path.join('simnibs', 'cython_code', os.path.basename(f)))


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
          'nibabel>=2.3', 'pandas', 'nose'
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
