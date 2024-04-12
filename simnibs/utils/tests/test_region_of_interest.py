import os
import numpy as np
import nibabel as nib
import pytest
from simnibs.mesh_tools import mesh_io
from simnibs.mesh_tools.mesh_io import Msh
from simnibs.utils import file_finder
from simnibs.utils.file_finder import SubjectFiles
import simnibs.utils.region_of_interest as region_of_interest
from simnibs.utils.region_of_interest import RegionOfInterest


class TestLoadVolume:
    def test_load_from_msh_instance(self, sphere3_msh: Msh):
        roi = RegionOfInterest()
        roi.load_mesh(mesh=sphere3_msh)
        np.testing.assert_allclose(
            roi.get_roi_mesh().nodes.node_coord, sphere3_msh.nodes.node_coord
        )
        np.testing.assert_allclose(
            roi.get_roi_mesh().elm.node_number_list, sphere3_msh.elm.node_number_list
        )
        np.testing.assert_allclose(
            roi.get_nodes(), sphere3_msh.elements_baricenters().value
        )
        assert np.all(roi._mask)
        assert len(roi._mask) == sphere3_msh.elm.nr

    def test_load_from_msh_path(self, sphere3_msh: Msh):
        roi = RegionOfInterest()
        roi.load_mesh(mesh=str(sphere3_msh.fn))
        np.testing.assert_allclose(
            roi.get_roi_mesh().nodes.node_coord, sphere3_msh.nodes.node_coord
        )
        np.testing.assert_allclose(
            roi.get_roi_mesh().elm.node_number_list, sphere3_msh.elm.node_number_list
        )
        np.testing.assert_allclose(
            roi.get_nodes(), sphere3_msh.elements_baricenters().value
        )
        assert np.all(roi._mask)
        assert len(roi._mask) == sphere3_msh.elm.nr

    @pytest.mark.slow
    def test_load_from_m2m_folder(self, example_dataset):
        m2m_path = os.path.join(example_dataset, "m2m_ernie")

        roi = RegionOfInterest()
        roi.load_mesh(subpath=m2m_path)
        ernie = mesh_io.read_msh(os.path.join(m2m_path, "ernie.msh"))
        np.testing.assert_allclose(
            roi.get_roi_mesh().nodes.node_coord, ernie.nodes.node_coord
        )
        np.testing.assert_allclose(
            roi.get_roi_mesh().elm.node_number_list, ernie.elm.node_number_list
        )
        np.testing.assert_allclose(roi.get_nodes(), ernie.elements_baricenters().value)
        assert np.all(roi._mask)
        assert len(roi._mask) == ernie.elm.nr

    def test_load_none(self):
        with pytest.raises(Exception) as e_info:
            roi = RegionOfInterest()
            roi.load_mesh()


class TestLoadSurface:
    @pytest.mark.slow
    def test_load_central_gm_surface(self, example_dataset):
        m2m_path = os.path.join(example_dataset, "m2m_ernie")
        subject_files = SubjectFiles(subpath=m2m_path)

        lh_central = mesh_io.read_gifti_surface(
            subject_files.get_surface("lh", "central")
        )
        rh_central = mesh_io.read_gifti_surface(
            subject_files.get_surface("rh", "central")
        )
        central = lh_central.join_mesh(rh_central)

        roi = RegionOfInterest()
        roi.load_surfaces(surface_type="central", subpath=m2m_path)
        np.testing.assert_allclose(
            roi.get_roi_mesh().nodes.node_coord, central.nodes.node_coord
        )
        np.testing.assert_allclose(
            roi.get_roi_mesh().elm.node_number_list, central.elm.node_number_list
        )
        np.testing.assert_allclose(roi.get_nodes(), central.nodes.node_coord)
        assert np.all(roi._mask)
        assert len(roi._mask) == central.nodes.nr
        assert roi._surface_divide == lh_central.nodes.nr

    def test_load_custom_surface(self, sphere3_msh: Msh, tmp_path):
        surface = sphere3_msh.crop_mesh(1004)
        surface_path = os.path.join(tmp_path, "surface.msh")
        surface.write(surface_path)

        roi = RegionOfInterest()
        roi.load_surfaces(surface_type="custom", surface_path=surface_path)
        np.testing.assert_allclose(
            roi.get_roi_mesh().nodes.node_coord, surface.nodes.node_coord
        )
        np.testing.assert_allclose(
            roi.get_roi_mesh().elm.node_number_list, surface.elm.node_number_list
        )
        np.testing.assert_allclose(roi.get_nodes(), surface.nodes.node_coord)
        assert np.all(roi._mask)
        assert len(roi._mask) == surface.nodes.nr


class TestLoadSurfaceFromFile:
    def test_load_msh(self, sphere3_msh: Msh, tmp_path):
        surface = sphere3_msh.crop_mesh(1004)
        surface_path = os.path.join(tmp_path, "surface.msh")
        surface.write(surface_path)
        loaded_surface = region_of_interest.load_surface_from_file(surface_path)

        np.testing.assert_allclose(
            loaded_surface.nodes.node_coord, surface.nodes.node_coord
        )
        np.testing.assert_allclose(
            loaded_surface.elm.node_number_list, surface.elm.node_number_list
        )

    def test_load_msh_with_tetrahedrons(self, sphere3_msh: Msh, tmp_path):
        surface = sphere3_msh
        surface_path = os.path.join(tmp_path, "surface.msh")
        surface.write(surface_path)
        with pytest.raises(Exception) as e_info:
            region_of_interest.load_surface_from_file(surface_path)

    def test_load_gifti(self, sphere3_msh: Msh, tmp_path):
        surface = sphere3_msh.crop_mesh(1004)
        surface_path = os.path.join(tmp_path, "surface.gii")
        mesh_io.write_gifti_surface(surface, surface_path)

        loaded_surface = region_of_interest.load_surface_from_file(surface_path)

        np.testing.assert_allclose(
            loaded_surface.nodes.node_coord, surface.nodes.node_coord
        )
        np.testing.assert_allclose(
            loaded_surface.elm.node_number_list, surface.elm.node_number_list
        )

    def test_load_fs_surface(self, sphere3_msh: Msh, tmp_path):
        surface = sphere3_msh.crop_mesh(1004)
        surface_path = os.path.join(tmp_path, "surface")
        mesh_io.write_freesurfer_surface(surface, surface_path)

        loaded_surface = region_of_interest.load_surface_from_file(surface_path)

        np.testing.assert_allclose(
            loaded_surface.nodes.node_coord, surface.nodes.node_coord
        )
        np.testing.assert_allclose(
            loaded_surface.elm.node_number_list, surface.elm.node_number_list
        )


class TestLoadSurfaceMaskFromFile:
    def test_load_label(self, sphere3_msh: Msh, tmp_path):
        label_file_path = os.path.join(tmp_path, "test.label")
        vertex_indexes = []
        for i in range(sphere3_msh.nodes.nr):
            if i % 3:
                vertex_indexes.append(i)
        with open(label_file_path, "w") as file:
            file.write(f"{len(vertex_indexes)}{os.linesep}")
            for vertex in vertex_indexes:
                file.write(
                    f"{vertex} {1.345 * vertex} {0.23 * vertex} {12.3213 * vertex} {0.000}{os.linesep}"
                )

        index_mask = region_of_interest.load_surface_mask_from_file(label_file_path)
        assert np.all(index_mask % 3)

    def test_load_annot(self, sphere3_msh: Msh, tmp_path):
        annot_file_path = os.path.join(tmp_path, "test.annot")
        vertex_labels = []
        names = []
        for i in range(sphere3_msh.nodes.nr):
            vertex_labels.append(np.rint(np.sqrt(i)))
        for i in np.unique(vertex_labels):
            names.append(f"{i}-name")

        nib.freesurfer.io.write_annot(
            annot_file_path,
            np.array(vertex_labels, dtype=np.int_),
            np.full((len(np.unique(vertex_labels)), 4), 122),
            names,
        )

        index_mask = region_of_interest.load_surface_mask_from_file(annot_file_path, 3)
        assert np.all(np.rint(np.sqrt(index_mask)) == 3)

    def test_load_curv(self, sphere3_msh: Msh, tmp_path):
        curv_file_path = os.path.join(tmp_path, "test")
        vertex_labels = []
        for i in range(sphere3_msh.nodes.nr):
            vertex_labels.append(np.rint(np.sqrt(i)))

        nib.freesurfer.io.write_morph_data(curv_file_path, np.array(vertex_labels))

        index_mask = region_of_interest.load_surface_mask_from_file(curv_file_path, 3)
        assert np.all(np.rint(np.sqrt(index_mask)) == 3)


class TestMaskImageToIndexMask:
    def test_transform_node_mask(self, sphere3_msh: Msh):
        mask_array = np.zeros((500, 500, 500))
        mask_array[:, :, :250] = 1

        mask_affine = np.eye(4)
        mask_affine[:3, 3] = [-250, -250, -250]
        mask_img = nib.Nifti1Image(mask_array, mask_affine)

        index_mask = region_of_interest.mask_image_to_index_mask("node", mask_img, sphere3_msh, mask_value=1)

        bool_mask = np.ones((sphere3_msh.nodes.nr), dtype=np.bool_)
        bool_mask[index_mask] = False
        assert len(index_mask) > 0
        assert np.any(bool_mask)

        assert np.all(sphere3_msh.nodes.node_coord[bool_mask][:,2] > -0.5)
        assert np.all(sphere3_msh.nodes.node_coord[index_mask][:,2] < 0.5)

    def test_transform_element_mask(self, sphere3_msh: Msh):
        mask_array = np.zeros((500, 500, 500))
        mask_array[:, :, :250] = 1

        mask_affine = np.eye(4)
        mask_affine[:3, 3] = [-250, -250, -250]
        mask_img = nib.Nifti1Image(mask_array, mask_affine)

        index_mask = region_of_interest.mask_image_to_index_mask("elm_center", mask_img, sphere3_msh, mask_value=1)

        bool_mask = np.ones((sphere3_msh.elm.nr), dtype=np.bool_)
        bool_mask[index_mask] = False
        assert len(index_mask) > 0
        assert np.any(bool_mask)
        assert np.all(sphere3_msh.elements_baricenters().value[bool_mask][:,2] > -0.5)
        assert np.all(sphere3_msh.elements_baricenters().value[index_mask][:,2] < 0.5)
        
class TestMniMaskToSub:
    @pytest.mark.slow
    @pytest.mark.skip(reason="Example dataset is missing label_prep/T1_upsampled.nii.gz")
    def test_transform_final_tissues(self, example_dataset):
        m2m_path = os.path.join(example_dataset, "m2m_ernie")
        subject_files = SubjectFiles(subpath=m2m_path)

        mask_img = nib.load(subject_files.final_labels_MNI)
        mni_to_sub_final_tissues = region_of_interest.mni_mask_to_sub(mask_img, subject_files)

        sub_final_tissues = nib.load(subject_files.final_labels)

        np.testing.assert_allclose(sub_final_tissues.affine, mni_to_sub_final_tissues.affine)
        np.testing.assert_allclose(sub_final_tissues.get_fdata(), mni_to_sub_final_tissues.get_fdata())

class TestCombineMask:
    def test_union_distinct(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[:20] = True

        index_mask = np.arange(30, 50)

        result_bool_mask = region_of_interest.combine_mask(bool_mask, index_mask, "union")

        assert np.all(result_bool_mask[:20])
        assert np.all(result_bool_mask[30:])
        assert not np.any(result_bool_mask[20:30])

    def test_union_overlap(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[10:30] = True

        index_mask = np.arange(20, 40)

        result_bool_mask = region_of_interest.combine_mask(bool_mask, index_mask, "union")

        assert np.all(result_bool_mask[10:40])
        assert not np.any(result_bool_mask[:10])
        assert not np.any(result_bool_mask[40:])

    def test_union_equal(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[20:30] = True

        index_mask = np.arange(20, 30)

        result_bool_mask = region_of_interest.combine_mask(bool_mask, index_mask, "union")

        assert np.all(result_bool_mask[20:30])
        assert not np.any(result_bool_mask[:20])
        assert not np.any(result_bool_mask[30:])

    def test_union_inside(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[10:40] = True

        index_mask = np.arange(20, 30)

        result_bool_mask = region_of_interest.combine_mask(bool_mask, index_mask, "union")

        assert np.all(result_bool_mask[10:40])
        assert not np.any(result_bool_mask[:10])
        assert not np.any(result_bool_mask[40:])

    def test_intersection_distinct(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[:20] = True

        index_mask = np.arange(30, 50)

        result_bool_mask = region_of_interest.combine_mask(bool_mask, index_mask, "intersection")

        assert not np.any(result_bool_mask)

    def test_intersection_overlap(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[10:30] = True

        index_mask = np.arange(20, 40)

        result_bool_mask = region_of_interest.combine_mask(bool_mask, index_mask, "intersection")

        assert np.all(result_bool_mask[20:30])
        assert not np.any(result_bool_mask[:20])
        assert not np.any(result_bool_mask[30:])

    def test_intersection_equal(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[20:30] = True

        index_mask = np.arange(20, 30)

        result_bool_mask = region_of_interest.combine_mask(bool_mask, index_mask, "intersection")

        assert np.all(result_bool_mask[20:30])
        assert not np.any(result_bool_mask[:20])
        assert not np.any(result_bool_mask[30:])

    def test_intersection_inside(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[10:40] = True

        index_mask = np.arange(20, 30)

        result_bool_mask = region_of_interest.combine_mask(bool_mask, index_mask, "intersection")

        assert np.all(result_bool_mask[20:30])
        assert not np.any(result_bool_mask[:20])
        assert not np.any(result_bool_mask[30:])

    def test_difference_distinct(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[:20] = True

        index_mask = np.arange(30, 50)

        result_bool_mask = region_of_interest.combine_mask(bool_mask, index_mask, "difference")

        assert np.all(result_bool_mask[:20])
        assert not np.any(result_bool_mask[20:])

    def test_difference_overlap(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[10:30] = True

        index_mask = np.arange(20, 40)

        result_bool_mask = region_of_interest.combine_mask(bool_mask, index_mask, "difference")

        assert np.all(result_bool_mask[10:20])
        assert not np.any(result_bool_mask[:10])
        assert not np.any(result_bool_mask[20:])

    def test_difference_equal(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[20:30] = True

        index_mask = np.arange(20, 30)

        result_bool_mask = region_of_interest.combine_mask(bool_mask, index_mask, "difference")

        assert not np.any(result_bool_mask)

    def test_difference_inside(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[10:40] = True

        index_mask = np.arange(20, 30)

        result_bool_mask = region_of_interest.combine_mask(bool_mask, index_mask, "difference")

        assert np.all(result_bool_mask[10:20])
        assert np.all(result_bool_mask[30:40])
        assert not np.any(result_bool_mask[:10])
        assert not np.any(result_bool_mask[20:30])
        assert not np.any(result_bool_mask[40:])

class TestFsAvrMaskToSub:
    @pytest.mark.slow        
    def test_full_mask(self, example_dataset):
        m2m_path = os.path.join(example_dataset, "m2m_ernie")
        subject_files = SubjectFiles(subpath=m2m_path)

        # 163842 nodes in fs avr per hemi
        index_mask = np.arange(163842, dtype=np.int_)

        result_index_mask = region_of_interest.fs_avr_mask_to_sub(index_mask, "lh", subject_files)

        central_lh = mesh_io.read_gifti_surface(subject_files.get_surface('lh', 'central'))
        np.testing.assert_equal(result_index_mask, np.arange(central_lh.nodes.nr))

    @pytest.mark.slow        
    def test_empty_mask(self, example_dataset):
        m2m_path = os.path.join(example_dataset, "m2m_ernie")
        subject_files = SubjectFiles(subpath=m2m_path)

        # 163842 nodes in fs avr per hemi
        index_mask = np.array([], dtype=np.int_)

        result_index_mask = region_of_interest.fs_avr_mask_to_sub(index_mask, "lh", subject_files)

        assert len(result_index_mask) == 0

    @pytest.mark.slow        
    def test_half_mask(self, example_dataset):
        m2m_path = os.path.join(example_dataset, "m2m_ernie")
        subject_files = SubjectFiles(subpath=m2m_path)

        sphere: Msh = mesh_io.read_gifti_surface(
            file_finder.get_reference_surf('lh', "sphere")
        )

        # 163842 nodes in fs avr per hemi
        index_mask = np.arange(163842, dtype=np.int_)[sphere.nodes.node_coord[:, 2] > 0]

        result_index_mask = region_of_interest.fs_avr_mask_to_sub(index_mask, "lh", subject_files)

        sphere_reg_lh = mesh_io.read_gifti_surface(subject_files.get_surface('lh', 'sphere_reg'))
        assert np.all(np.isin(result_index_mask, np.arange(sphere_reg_lh.nodes.nr)[sphere_reg_lh.nodes.node_coord[:, 2] > 0]))
        assert not np.any(np.isin(result_index_mask, np.arange(sphere_reg_lh.nodes.nr)[sphere_reg_lh.nodes.node_coord[:, 2] < 0]))

