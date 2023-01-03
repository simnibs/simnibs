''' Automated testing of the examples
'''
import sys
import os
import stat
import subprocess
import shutil
import glob
import tempfile
import requests
import functools
import zipfile

import pytest
from simnibs.utils.file_finder import path2bin
import simnibs

EXAMPLES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

@pytest.fixture(scope="module")
def example_dataset():
    url = (
        f'https://github.com/simnibs/example-dataset/releases/'
        'download/v4.0-lowres/ernie_lowres_V2.zip'
    )
    fn_folder = os.path.abspath(os.path.join(os.path.dirname(__file__),'test_data'))
    tmpname = tempfile.mktemp(".zip")
    # Download the dataset
    with requests.get(url, stream=True) as r:
        r.raw.read = functools.partial(r.raw.read, decode_content=True)
        with open(tmpname, 'wb') as f:
            shutil.copyfileobj(r.raw, f)
    # Unzip the dataset
    with zipfile.ZipFile(tmpname) as z:
        z.extractall(fn_folder, )
    os.remove(tmpname)
    yield fn_folder
    try:
        shutil.rmtree(fn_folder)
        pass
    except:
        print('Could not remove example dataset folder')


@pytest.fixture
def replace_gmsh():
    fn_gmsh = path2bin('gmsh')
    fn_gmsh_tmp = path2bin('gmsh_tmp')
    # move
    if sys.platform == 'win32': fn_gmsh += '.exe'
    shutil.move(fn_gmsh, fn_gmsh_tmp)
    # replace
    if sys.platform == 'win32':
        # replace gmsh with an .exe that does not to anything
        shutil.copy(r'C:\Windows\System32\rundll32.exe', fn_gmsh)
    else:
        with open(fn_gmsh, 'w') as f:
            f.write('#! /bin/bash -e\n')
            f.write(f'"echo" "$@"')
        os.chmod(
            fn_gmsh,
            os.stat(fn_gmsh).st_mode |
            stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
        )
    yield
    shutil.move(fn_gmsh_tmp, fn_gmsh)

def octave_call(script):
    if sys.platform == 'win32':
        octave_exe = glob.glob(os.path.join(
            r'C:\Octave\Octave-*\mingw64\bin\octave-cli.exe'
        ))
        if len(octave_exe) == 0:
            raise OSError('Did not find Octave executable')
        else:
            octave_exe = octave_exe[0]
    else:
        octave_exe = 'octave'

    matlab_dir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), '..', '..', 'matlab_tools')
    )
    octave_dir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), 'octave_compatibility')
    )
    cmd = f"{octave_exe} -W --eval \""
    cmd += "addpath('{0}');".format(matlab_dir)
    cmd += "addpath('{0}');".format(octave_dir)
    cmd += "try,run('{}');catch ME,rethrow(ME);end,exit;\"".format(script)
    return cmd

@pytest.fixture
def replace_show_surface():
    fn_orig = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            '..', '..', 'matlab_tools',
            'mesh_show_surface.m'
        )
    )
    fn_tmp = fn_orig + '.bk'
    shutil.move(fn_orig, fn_tmp)
    # replace
    with open(fn_orig, 'w') as f:
        f.write('function varargout = mesh_show_surface(m, varargin)\n')
        f.write('end')
    yield
    shutil.move(fn_tmp, fn_orig)

class TestPythonErnie:
    def run_script(self, script_folder, script_name, clean=None):
        if clean is not None and os.path.exists(clean):
            shutil.rmtree(clean)
        fn = os.path.join(EXAMPLES_DIR, script_folder, script_name)
        print(fn)
        return subprocess.run([sys.executable, fn])

    def test_transform_coordinates(self, example_dataset):
        os.chdir(example_dataset)
        ret = self.run_script('analysis', 'transform_coordinates.py')
        assert ret.returncode == 0

    def test_tDCS(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('simulations', 'tDCS.py', 'tdcs_simu')
        assert ret.returncode == 0

    def test_roi_analysis_surf(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        simnibs.utils.transformations.middle_gm_interpolation(
            'tdcs_simu/ernie_TDCS_1_scalar.msh',
            'm2m_ernie',
            'tdcs_simu/subject_overlays'
        )
        ret = self.run_script('analysis', 'roi_analysis_surf.py')
        assert ret.returncode == 0

    def test_roi_analysis_mni(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('analysis', 'roi_analysis_mni.py')
        assert ret.returncode == 0

    def test_TMS(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('simulations', 'TMS.py', 'tms_simu')
        assert ret.returncode == 0

    def test_tDCS_ring(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('simulations', 'tDCS_ring.py', 'tdcs_ring')
        assert ret.returncode == 0

    def test_TMS_MNI(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('simulations', 'TMS_MNI.py', 'tms_hand')
        assert ret.returncode == 0

    def test_leadfield(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'leadfield.py', 'leadfield')
        assert ret.returncode == 0

    def test_tdcs_optimize(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'teslf_optimize.py')
        assert ret.returncode == 0

    def test_tdcs_optimize_multi_target(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'teslf_optimize_multi_target.py')
        assert ret.returncode == 0

    def test_tdcs_optimize_avoid(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'teslf_optimize_avoid.py')
        assert ret.returncode == 0

    def test_tdcs_optimize_mni(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'teslf_optimize_mni.py')
        assert ret.returncode == 0

    def test_tdcs_optimize_strength(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'teslf_optimize_strength.py')
        assert ret.returncode == 0

    def test_tdcs_optimize_distributed(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        shutil.copy(
            os.path.join(
                os.path.dirname(simnibs.SIMNIBSDIR), 'docs',
                'data', 'ID03_MOTOR_ICA.nii.gz'),
            example_dataset
        )
        ret = self.run_script('optimization', 'teslf_optimize_distributed.py')
        assert ret.returncode == 0

    def test_tms_optimize_ADM(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'tms_optimization_adm.py', 'tms_optimization_adm')
        assert ret.returncode == 0

    def test_tes_electrode_setup_init(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'tes_electrode_setup_init.py', clean=None)
        assert ret.returncode == 0


class TestMatlabErnie:
    def run_script(self, script_folder, script_name, clean=None):
        if clean is not None and os.path.exists(clean):
            shutil.rmtree(clean)
        shutil.copy(
            os.path.join(EXAMPLES_DIR, script_folder, script_name),
            script_name
        )
        return subprocess.run(octave_call(script_name), shell=True)

    def test_transform_coordinates(self, example_dataset):
        os.chdir(example_dataset)
        ret = self.run_script('analysis', 'transform_coordinates.m')
        assert ret.returncode == 0

    def test_tDCS(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('simulations', 'tDCS.m', 'tdcs_simu')
        assert ret.returncode == 0

    def test_roi_analysis_surf(self, example_dataset, replace_show_surface):
        os.chdir(example_dataset)
        simnibs.transformations.middle_gm_interpolation(
            'tdcs_simu/ernie_TDCS_1_scalar.msh',
            'm2m_ernie',
            'tdcs_simu/subject_overlays'
        )
        ret = self.run_script('analysis', 'roi_analysis_surf.m')
        assert ret.returncode == 0

    def test_roi_analysis_mni(self, example_dataset, replace_show_surface):
        os.chdir(example_dataset)
        ret = self.run_script('analysis', 'roi_analysis_mni.m')
        assert ret.returncode == 0

    def test_TMS(self, example_dataset, replace_show_surface):
        os.chdir(example_dataset)
        ret = self.run_script('simulations', 'TMS.m', 'tms_simu')
        assert ret.returncode == 0

    def test_tDCS_ring(self, example_dataset, replace_show_surface):
        os.chdir(example_dataset)
        ret = self.run_script('simulations', 'tDCS_ring.m', 'tdcs_ring')
        assert ret.returncode == 0

    def test_TMS_MNI(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('simulations', 'TMS_MNI.m', 'tms_hand')
        assert ret.returncode == 0

    def test_leadfield(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'leadfield.m', 'leadfield')
        assert ret.returncode == 0

    def test_tdcs_optimize(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'teslf_optimize.m')
        assert ret.returncode == 0

    def test_tdcs_optimize_multi_target(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'teslf_optimize_multi_target.m')
        assert ret.returncode == 0

    def test_tdcs_optimize_avoid(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'teslf_optimize_avoid.m')
        assert ret.returncode == 0

    def test_tdcs_optimize_mni(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'teslf_optimize_mni.m')
        assert ret.returncode == 0

    def test_tdcs_optimize_strength(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'teslf_optimize_strength.m')
        assert ret.returncode == 0

    def test_tdcs_optimize_distributed(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        shutil.copy(
            os.path.join(
                os.path.dirname(simnibs.SIMNIBSDIR), 'docs',
                'data', 'ID03_MOTOR_ICA.nii.gz'),
            example_dataset
        )
        ret = self.run_script('optimization', 'teslf_optimize_distributed.m')
        assert ret.returncode == 0

    def test_tms_optimize_ADM(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'tms_optimization_adm.m', 'tms_optimization_adm')
        assert ret.returncode == 0
