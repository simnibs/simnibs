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
            surface_folder = Path(file_finder.SubjectFiles(subpath=tmpdir).surface_folder)
            surface_folder.mkdir(parents=True)
            regions = ['lc', 'lh', 'rh']
            surf_types = ('central', 'sphere_reg')
            surf_to_name = dict(central='central', sphere_reg='sphere.reg')
            subsamplings = (None, 10000)
            fnames = {}
            for region in regions:
                for surf_type in surf_types:
                    for subsampling in subsamplings:
                        fileparts = [region, surf_to_name[surf_type]]
                        if subsampling is not None:
                            fileparts.append(str(subsampling))
                        fileparts.append('gii')
                        filename = '.'.join(fileparts)
                        fn = surface_folder / filename
                        fn.touch()
                        fnames[(region, surf_type, subsampling)] = str(fn)

            # Needs to be reinitialized
            s = file_finder.SubjectFiles(subpath=tmpdir)
            assert s.regions == regions
            for region in regions:
                for surf_type in surf_types:
                    for subsampling in subsamplings:
                        assert s.get_surface(region, surf_type, subsampling) \
                               == fnames[region, surf_type, subsampling]

            with pytest.raises(FileNotFoundError):
                s.get_surface('rc', 'central')
                s.get_surface('rc', 'sphere_reg')
                s.get_surface('lh', 'central', 20000)

    def test_v2v2_get_surface_and_morph_files(self):
        """Check that the """
        m2m_dir = "/path/to/m2m_subid"
        sf = file_finder.SubjectFiles(subpath=m2m_dir)

        # surface
        surf = "mysurf"
        subsamp = 12345

        res = sf.v2v2_get_surface_file(surf)
        assert res["lh"] == Path(f"{m2m_dir}/surfaces/lh.{surf}.gii")
        assert res["rh"] == Path(f"{m2m_dir}/surfaces/rh.{surf}.gii")

        res = sf.v2v2_get_surface_file(surf, subsamp)
        assert res["lh"] == Path(f"{m2m_dir}/surfaces/{subsamp}/lh.{surf}.gii")
        assert res["rh"] == Path(f"{m2m_dir}/surfaces/{subsamp}/rh.{surf}.gii")

        res = sf.v2v2_get_surface_file(surf, subsamp, "lh")
        assert res["lh"] == Path(f"{m2m_dir}/surfaces/{subsamp}/lh.{surf}.gii")
        with pytest.raises(KeyError):
            res["rh"]

        # morph data
        data = "mydata"
        subsamp = 10000

        res = sf.v2v2_get_morph_data_file(data)
        assert res["lh"] == Path(f"{m2m_dir}/surfaces/lh.{data}")
        assert res["rh"] == Path(f"{m2m_dir}/surfaces/rh.{data}")

        res = sf.v2v2_get_morph_data_file(data, subsamp)
        assert res["lh"] == Path(f"{m2m_dir}/surfaces/{subsamp}/lh.{data}")
        assert res["rh"] == Path(f"{m2m_dir}/surfaces/{subsamp}/rh.{data}")

        res = sf.v2v2_get_morph_data_file(data, subsamp, "rh")
        with pytest.raises(KeyError):
            res["lh"]
        assert res["rh"] == Path(f"{m2m_dir}/surfaces/{subsamp}/rh.{data}")

    def test_v2v2_read_write_files(self):
        """Write a few files, read them again, and check equality."""
        rng = np.random.default_rng()

        # surface
        surf_hemi = dict(
            points=rng.uniform(0, 10, size=(10,3)),
            tris=rng.integers(0, 10, size=(16,3))
        )
        surf_dict = dict(lh=surf_hemi, rh=surf_hemi)
        surf_name = "mysurf"

        # data
        d = rng.uniform(0,10,size=10)
        data_dict = dict(lh=d, rh=d)
        data_name = "mydata"

        with tempfile.TemporaryDirectory(prefix='m2m_') as tmpdir:

            sf = file_finder.SubjectFiles(subpath=tmpdir)

            sf.v2v2_write_surface(surf_dict, surf_name)
            surf_read = sf.v2v2_read_surface(surf_name)
            for hemi in surf_dict:
                for field in surf_dict[hemi]:
                    np.testing.assert_allclose(surf_dict[hemi][field], surf_read[hemi][field])

            sf.v2v2_write_morph_data(data_dict, data_name)
            data_read = sf.v2v2_read_morph_data(data_name)
            for hemi in surf_dict:
                np.testing.assert_allclose(data_dict[hemi], data_read[hemi])
