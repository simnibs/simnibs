import sys
import os
import stat
import subprocess
import zipfile
import urllib.request
import tempfile
import shutil
import pytest

from simnibs.utils.file_finder import path2bin
import simnibs

EXAMPLES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'examples'))

@pytest.fixture(scope="module")
def example_dataset():
    url = (
        f'https://github.com/simnibs/example-dataset/releases/'
        'download/v3.0-lowres/ernie_lowres.zip'
    )
    fn_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), 'ernie_lowres'))
    tmpname = tempfile.mktemp(".zip")
    urllib.request.urlretrieve(url, tmpname)
    with zipfile.ZipFile(tmpname) as z:
        z.extractall(os.path.dirname(__file__), )
    os.remove(tmpname)
    yield fn_folder
    shutil.rmtree(fn_folder)


@pytest.fixture
def replace_gmsh():
    fn_gmsh = path2bin('gmsh')
    fn_gmsh_tmp = path2bin('gmsh_tmp')
    # move
    shutil.move(fn_gmsh, fn_gmsh_tmp)
    # replace
    if sys.platform == 'win32':
        fn_script = fn_gmsh[:4] + '.cmd'
        with open(fn_script, 'w') as f:
            f.write('echo "GMSH"')
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
    cmd = "octave -W --eval \""
    matlab_dir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), '..', '..', 'matlab_tools')
    )
    octave_dir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), 'octave_compatibility')
    )
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
        return subprocess.run([sys.executable, fn])

    def test_transform_coordinates(self, example_dataset):
        os.chdir(example_dataset)
        ret = self.run_script('analysis', 'transform_coordinates.py')
        assert ret.returncode == 0

    def test_tDCS(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('simulations', 'tDCS.py', 'tdcs')
        assert ret.returncode == 0

    def test_roi_analysis_surf(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        simnibs.transformations.middle_gm_interpolation(
            'tdcs/ernie_TDCS_1_scalar.msh',
            'm2m_ernie',
            'tdcs/subject_overlays'
        )
        ret = self.run_script('analysis', 'roi_analysis_surf.py')
        assert ret.returncode == 0

    def test_roi_analysis_mni(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('analysis', 'roi_analysis_mni.py')
        assert ret.returncode == 0

    def test_TMS(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('simulations', 'TMS.py', 'tms')
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
        ret = self.run_script('optimization', 'tdcs_optimize.py')
        assert ret.returncode == 0

    def test_tdcs_optimize_multi_target(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'tdcs_optimize_multi_target.py')
        assert ret.returncode == 0

    def test_tdcs_optimize_avoid(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'tdcs_optimize_avoid.py')
        assert ret.returncode == 0

    def test_tdcs_optimize_mni(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'tdcs_optimize_mni.py')
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
        ret = self.run_script('simulations', 'tDCS.m', 'tdcs')
        assert ret.returncode == 0

    def test_roi_analysis_surf(self, example_dataset, replace_show_surface):
        os.chdir(example_dataset)
        simnibs.transformations.middle_gm_interpolation(
            'tdcs/ernie_TDCS_1_scalar.msh',
            'm2m_ernie',
            'tdcs/subject_overlays'
        )
        ret = self.run_script('analysis', 'roi_analysis_surf.m')
        assert ret.returncode == 0

    def test_roi_analysis_mni(self, example_dataset, replace_show_surface):
        os.chdir(example_dataset)
        ret = self.run_script('analysis', 'roi_analysis_mni.m')
        assert ret.returncode == 0

    def test_TMS(self, example_dataset, replace_show_surface):
        os.chdir(example_dataset)
        ret = self.run_script('simulations', 'TMS.m', 'tms')
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
        ret = self.run_script('optimization', 'tdcs_optimize.m')
        assert ret.returncode == 0

    def test_tdcs_optimize_multi_target(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'tdcs_optimize_multi_target.m')
        assert ret.returncode == 0

    def test_tdcs_optimize_avoid(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'tdcs_optimize_avoid.m')
        assert ret.returncode == 0

    def test_tdcs_optimize_mni(self, example_dataset, replace_gmsh):
        os.chdir(example_dataset)
        ret = self.run_script('optimization', 'tdcs_optimize_mni.m')
        assert ret.returncode == 0
