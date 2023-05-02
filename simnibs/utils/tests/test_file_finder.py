import os
import numpy as np
import pytest
import tempfile
from pathlib import Path

from .. import file_finder


class TestTemplates:
    def test_find(self):
        templates = file_finder.templates
        for k, fn in templates.__dict__.items():
            if k not in ("fsaverage_resolutions", ):
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

@pytest.mark.parametrize('surf_type', file_finder.fs_surfaces + ["spheree"])
@pytest.mark.parametrize('region', file_finder.HEMISPHERES + ["lhh"])
@pytest.mark.parametrize('resolution', file_finder.fs_resolutions + [20])
def test_get_reference_surf(surf_type, region, resolution):
    if resolution not in file_finder.fs_resolutions:
        with pytest.raises(AssertionError):
            file_finder.get_reference_surf(surf_type, region, resolution)
    elif (
        region not in file_finder.HEMISPHERES
        or surf_type not in file_finder.fs_surfaces
    ):
        with pytest.raises(FileNotFoundError):
            file_finder.get_reference_surf(surf_type, region, resolution)
    else:
        file_finder.get_reference_surf(surf_type, region, resolution)


class TestSubjectFiles:
    def test_define_fnamehead(self):
        s = file_finder.SubjectFiles(
            os.path.join('path', 'm2m_sub', 'sub.msh'))
        assert s.fnamehead == os.path.join('path', 'm2m_sub', 'sub.msh')
        assert s.subid == 'sub'
        assert s.subpath == os.path.join('path', 'm2m_sub')

    def test_define_subpath(self):
        s = file_finder.SubjectFiles(subpath=os.path.join('path', 'to', 'm2m_sub'))
        assert s.subpath == os.path.join('path', 'to', 'm2m_sub')
        assert s.subid == 'sub'
        assert s.fnamehead == os.path.join('path', 'to', 'm2m_sub', 'sub.msh')

    def test_define_fnamehead_subpath(self):
        s = file_finder.SubjectFiles(
            os.path.join('some', 'random', 'file.msh'),
            subpath=os.path.join('path', 'to', 'm2m_sub'))

        assert s.fnamehead == os.path.join('some', 'random', 'file.msh')
        assert s.subid == 'sub'
        assert s.subpath == os.path.join('path', 'to', 'm2m_sub')

    def test_tensor_file(self):
        s = file_finder.SubjectFiles(os.path.join('path', 'm2m_sub', 'sub.msh'))
        assert s.tensor_file == os.path.join('path', 'm2m_sub',
                                             'DTI_coregT1_tensor.nii.gz')

    def test_cap_files(self):
        s = file_finder.SubjectFiles(os.path.join('path', 'm2m_sub', 'sub.msh'))
        assert s.eeg_cap_folder == os.path.join('path', 'm2m_sub', 'eeg_positions')
        assert s.eeg_cap_1010 == os.path.join('path', 'm2m_sub', 'eeg_positions',
                                              'EEG10-10_UI_Jurak_2007.csv')

        assert s.get_eeg_cap('test.csv') == os.path.join(
            'path', 'm2m_sub', 'eeg_positions', 'test.csv')

    def test_surfaces(self):
        with tempfile.TemporaryDirectory(prefix='m2m_') as tmpdir:
            m2m = file_finder.SubjectFiles(subpath=tmpdir)
            surface_folder = Path(m2m.surface_folder)

            for s in ("central", "pial", "sphere", "sphere.reg"):
                for h in ("lh", "rh"):
                    assert m2m.surfaces[s][h] == surface_folder / f"{h}.{s}.gii"

            # surface_folder = Path(file_finder.SubjectFiles(subpath=tmpdir).surface_folder)
            # surface_folder.mkdir(parents=True)
            # regions = ['lc', 'lh', 'rh']
            # surf_types = ('central', 'sphere_reg')
            # surf_to_name = dict(central='central', sphere_reg='sphere.reg')
            # subsamplings = (None, 10000)
            # fnames = {}
            # for region in regions:
            #     for surf_type in surf_types:
            #         for subsampling in subsamplings:
            #             fileparts = [region, surf_to_name[surf_type]]
            #             if subsampling is not None:
            #                 fileparts.append(str(subsampling))
            #             fileparts.append('gii')
            #             filename = '.'.join(fileparts)
            #             fn = surface_folder / filename
            #             fn.touch()
            #             fnames[(region, surf_type, subsampling)] = str(fn)

            # # Needs to be reinitialized
            # s = file_finder.SubjectFiles(subpath=tmpdir)
            # assert s.regions == regions
            # for region in regions:
            #     for surf_type in surf_types:
            #         for subsampling in subsamplings:
            #             assert s.get_surface(region, surf_type, subsampling) \
            #                    == fnames[region, surf_type, subsampling]

            # with pytest.raises(FileNotFoundError):
            #     s.get_surface('rc', 'central')
            #     s.get_surface('rc', 'sphere_reg')
            #     s.get_surface('lh', 'central', 20000)


    @pytest.mark.parametrize('surface', ["central", "layer0"])
    @pytest.mark.parametrize('hemi', ["lh", "rh"])
    @pytest.mark.parametrize('subsampling', [None, 12345])
    def test_get_surface(self, surface, hemi, subsampling):
        print(f"surface: {surface}")
        print(f"hemi: {hemi}")
        print(f"subsampling: {subsampling}")
        m2m = file_finder.SubjectFiles(subpath="/path/to/m2m_subid")
        surface_dir = Path(m2m.surface_folder)

        p = m2m.get_surface(hemi, surface, subsampling)
        if subsampling is None:
            assert p == surface_dir / f"{hemi}.{surface}.gii"
        else:
            assert p == surface_dir / str(subsampling) / f"{hemi}.{surface}.gii"

    @pytest.mark.parametrize('morph_data', ["thickness", "curv"])
    @pytest.mark.parametrize('hemi', ["lh", "rh"])
    @pytest.mark.parametrize('subsampling', [None, 12345])
    def test_get_morph_data(self, hemi, morph_data, subsampling):
        m2m = file_finder.SubjectFiles(subpath="/path/to/m2m_subid")
        surface_dir = Path(m2m.surface_folder)

        p = m2m.get_morph_data(hemi, morph_data, subsampling)
        if subsampling is None:
            assert p == surface_dir / f"{hemi}.{morph_data}"
        else:
            assert p == surface_dir / str(subsampling) / f"{hemi}.{morph_data}"
