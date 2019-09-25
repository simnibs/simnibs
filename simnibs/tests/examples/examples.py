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

@pytest.fixture
def example_dataset(scope="module"):
    url = (f'https://github.com/simnibs/example-dataset/'
            'releases/download/v3.0/simnibs_examples.zip')
    fn_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), 'simnibs_examples'))
    tmpname = tempfile.mktemp(".zip")
    urllib.request.urlretrieve(url, tmpname)
    with zipfile.ZipFile(tmpname) as z:
        z.extractall(os.path.dirname(__file__))
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
    cmd = "octave -W --eval "
    matlab_dir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), '..', '..', 'matlab')
    )
    cmd += "\"addpath('{0}');".format(matlab_dir)
    cmd += "try,run('{}');catch ME,rethrow(ME);end,exit;\"".format(script)
    return cmd

@pytest.fixture
def replace_show_surface():
    fn_orig = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            '..', '..', 'matlab',
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
    def run_script(self, script_name, clean=None):
        if clean is not None and os.path.exists(clean):
            shutil.rmtree(clean)
        fn = os.path.join(EXAMPLES_DIR, script_name)
        return subprocess.run([sys.executable, fn])

    def test_transform_coordinates(self, example_dataset):
        os.chdir(os.path.join(example_dataset, 'ernie'))
        ret = self.run_script('transform_coordinates.py')
        assert ret.returncode == 0

    def test_tDCS(self, example_dataset, replace_gmsh):
        os.chdir(os.path.join(example_dataset, 'ernie'))
        ret = self.run_script('tDCS.py', 'tdcs')
        assert ret.returncode == 0

    def test_roi_analysis_atlas(self, example_dataset, replace_gmsh):
        os.chdir(os.path.join(example_dataset, 'ernie'))
        simnibs.msh.transformations.middle_gm_interpolation(
            'tdcs/ernie_TDCS_1_scalar.msh',
            'm2m_ernie',
            'tdcs/subject_overlays'
        )
        ret = self.run_script('roi_analysis_atlas.py')
        assert ret.returncode == 0

    def test_roi_analysis_mni(self, example_dataset, replace_gmsh):
        os.chdir(os.path.join(example_dataset, 'ernie'))
        ret = self.run_script('roi_analysis_mni.py')
        assert ret.returncode == 0

    def test_TMS(self, example_dataset, replace_gmsh):
        os.chdir(os.path.join(example_dataset, 'ernie'))
        ret = self.run_script('TMS.py', 'tms')
        assert ret.returncode == 0

    def test_tDCS_ring(self, example_dataset, replace_gmsh):
        os.chdir(os.path.join(example_dataset, 'ernie'))
        ret = self.run_script('tDCS_ring.py', 'tdcs_ring')
        assert ret.returncode == 0


class TestMatlabErnie:
    def run_script(self, script_name, clean=None):
        if clean is not None and os.path.exists(clean):
            shutil.rmtree(clean)
        shutil.copy(
            os.path.join(EXAMPLES_DIR, script_name),
            script_name
        )
        return subprocess.run(octave_call(script_name), shell=True)

    def test_transform_coordinates(self, example_dataset):
        os.chdir(os.path.join(example_dataset, 'ernie'))
        ret = self.run_script('transform_coordinates.m')
        assert ret.returncode == 0

    def test_tDCS(self, example_dataset, replace_show_surface):
        os.chdir(os.path.join(example_dataset, 'ernie'))
        ret = self.run_script('tDCS.m', 'tdcs')
        assert ret.returncode == 0

    def test_roi_analysis_atlas(self, example_dataset, replace_show_surface):
        os.chdir(os.path.join(example_dataset, 'ernie'))
        simnibs.msh.transformations.middle_gm_interpolation(
            'tdcs/ernie_TDCS_1_scalar.msh',
            'm2m_ernie',
            'tdcs/subject_overlays'
        )
        ret = self.run_script('roi_analysis_atlas.m')
        assert ret.returncode == 0

    def test_roi_analysis_mni(self, example_dataset, replace_show_surface):
        os.chdir(os.path.join(example_dataset, 'ernie'))
        ret = self.run_script('roi_analysis_mni.m')
        assert ret.returncode == 0

    def test_TMS(self, example_dataset, replace_show_surface):
        os.chdir(os.path.join(example_dataset, 'ernie'))
        ret = self.run_script('TMS.m', 'tms')
        assert ret.returncode == 0

    def test_tDCS_ring(self, example_dataset, replace_show_surface):
        os.chdir(os.path.join(example_dataset, 'ernie'))
        ret = self.run_script('tDCS_ring.m', 'tdcs_ring')
        assert ret.returncode == 0

    def test_tDCS_advanced(self, example_dataset, replace_show_surface):
        os.chdir(os.path.join(example_dataset, 'ernie'))
        ret = self.run_script('tDCS_advanced.m', 'tdcs_simulation')
        assert ret.returncode == 0
