import os
import numpy as np
import pytest
import tempfile
import pathlib

from .. import file_finder


class TestTemplates:
    def test_find(self):
        templates = file_finder.templates
        for k, fn in templates.__dict__.items():
            is_file = os.path.isfile(fn)
            is_dir = os.path.isdir(fn)
            assert is_file or is_dir

@pytest.mark.parametrize('atlas_name', ['a2009s', 'DK40', 'HCP_MMP1'])
@pytest.mark.parametrize('hemi', ['lh', 'rh', 'both'])
def test_get_atlas(atlas_name, hemi):
    atlas = file_finder.get_atlas(atlas_name, hemi)

    if atlas_name == 'a2009':
        if hemi == 'both':
            assert len(atlas) == 77*2
        else:
            assert len(atlas) == 77
    elif atlas_name == 'DK40':
        if hemi == 'both':
            assert len(atlas) == 36*2
        else:
            assert len(atlas) == 36
    elif atlas_name == 'HCP_MMP1':
        if hemi == 'both':
            assert len(atlas) == 181*2
        else:
            assert len(atlas) == 181

    for name, mask in atlas.items():
        if hemi in ['lh', 'rh']:
            assert len(mask) == 163842
        else:
            assert len(mask) == 2*163842
            if name.startswith('lh'):
                assert not np.any(mask[163842:])
            if name.startswith('rh'):
                assert not np.any(mask[:163842])

@pytest.mark.parametrize('region', ['lh', 'rh', 'lc', 'rc'])
@pytest.mark.parametrize('surf_type', ['central', 'sphere', 'inflated'])
def test_get_reference_surf(region, surf_type):
    if (region in ['lc', 'rc']) and surf_type == 'inflated':
        with pytest.raises(FileNotFoundError):
            file_finder.get_reference_surf(region, surf_type)
    else:
        file_finder.get_reference_surf(region, surf_type)


class TestSubjectFiles:
    def test_define_fnamehead(self):
        s = file_finder.SubjectFiles(
            os.path.join('path', 'to', 'sub.msh'))
        assert s.fnamehead == os.path.join('path', 'to', 'sub.msh')
        assert s.subid == 'sub'
        assert s.basedir == os.path.join('path', 'to')
        assert s.subpath == os.path.join('path', 'to', 'm2m_sub')

    def test_define_subpath(self):
        s = file_finder.SubjectFiles(subpath=os.path.join('path', 'to', 'm2m_sub'))
        assert s.subpath == os.path.join('path', 'to', 'm2m_sub')
        assert s.subid == 'sub'
        assert s.basedir == os.path.join('path', 'to')
        assert s.fnamehead == os.path.join('path', 'to', 'sub.msh')

    def test_define_fnamehead_subpath(self):
        s = file_finder.SubjectFiles(
            os.path.join('some', 'random', 'file.msh'),
            subpath=os.path.join('path', 'to', 'm2m_sub'))

        assert s.fnamehead == os.path.join('some', 'random', 'file.msh')
        assert s.subid == 'sub'
        assert s.basedir == os.path.join('path', 'to')
        assert s.subpath == os.path.join('path', 'to', 'm2m_sub')

    def test_tensor_file(self):
        s = file_finder.SubjectFiles(os.path.join('path', 'to', 'sub.msh'))
        assert s.tensor_file == os.path.join('path', 'to', 'd2c_sub',
                                             'dti_results_T1space',
                                             'DTI_conf_tensor.nii.gz')

    def test_cap_files(self):
        s = file_finder.SubjectFiles(os.path.join('path', 'to', 'sub.msh'))
        assert s.eeg_cap_folder == os.path.join('path', 'to', 'm2m_sub',
                                                'eeg_positions')
        assert s.eeg_cap_1010 == os.path.join('path', 'to', 'm2m_sub',
                                              'eeg_positions',
                                              'EEG10-10_UI_Jurak_2007.csv')

        assert s.get_eeg_cap('test.csv') == os.path.join(
            'path', 'to', 'm2m_sub', 'eeg_positions', 'test.csv')
    
    def test_surfaces(self):
        with tempfile.TemporaryDirectory(prefix='m2m_') as tmpdir:
            surf_dir = file_finder.SubjectFiles(subpath=tmpdir).surf_dir
            os.makedirs(surf_dir)
            for region in ['lh', 'rh']:
                for surf_type in ['central', 'sphere', 'sphere.reg']:
                    pathlib.Path(
                        os.path.join(surf_dir, f'{region}.{surf_type}.SOMETHING.gii')
                    ).touch()
            pathlib.Path(os.path.join(surf_dir, 'lc.sphere.reg.SOMETHING.gii')).touch()
            pathlib.Path(os.path.join(surf_dir, 'rc.sphere.SOMETHING.gii')).touch()
            pathlib.Path(os.path.join(surf_dir, 'rh.thickness.SOMETHING')).touch()
            s = file_finder.SubjectFiles(subpath=tmpdir)

            assert s.get_surface('lh', 'central') == \
                    os.path.join(surf_dir, 'lh.central.SOMETHING.gii')
            assert s.get_surface('rh', 'central') == \
                    os.path.join(surf_dir, 'rh.central.SOMETHING.gii')
            assert s.get_surface('lc', 'sphere_reg') == \
                    os.path.join(surf_dir, 'lc.sphere.reg.SOMETHING.gii')
            assert s.get_surface('rc', 'sphere') == \
                    os.path.join(surf_dir, 'rc.sphere.SOMETHING.gii')
            with pytest.raises(FileNotFoundError):
                s.get_surface('lc', 'sphere')
            assert s.regions == ['lh', 'rh']
