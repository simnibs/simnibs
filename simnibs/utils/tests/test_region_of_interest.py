import os
import numpy as np
import nibabel as nib
import pytest
import scipy
from simnibs.mesh_tools import mesh_io
from simnibs.mesh_tools.mesh_io import Msh
from simnibs.utils import file_finder
from simnibs.utils.file_finder import SubjectFiles
from simnibs.utils.matlab_read import dict_from_matlab
from simnibs.utils.mesh_element_properties import ElementTags
import simnibs.utils.region_of_interest as region_of_interest
from simnibs.utils.region_of_interest import RegionOfInterest


class TestLoadVolume:
    def test_load_from_msh_instance(self, sphere3_msh: Msh):
        roi = RegionOfInterest()
        roi.load_mesh(mesh=sphere3_msh)
        assert not np.any(roi._mask)
        np.testing.assert_allclose(
            roi._mesh.nodes.node_coord, sphere3_msh.nodes.node_coord
        )
        np.testing.assert_allclose(
            roi._mesh.elm.node_number_list, sphere3_msh.elm.node_number_list
        )
        assert len(roi.get_nodes()) == 0

    def test_load_from_msh_path(self, sphere3_msh: Msh):
        roi = RegionOfInterest()
        roi.load_mesh(mesh=str(sphere3_msh.fn))
        assert not np.any(roi._mask)
        np.testing.assert_allclose(
            roi._mesh.nodes.node_coord, sphere3_msh.nodes.node_coord
        )
        np.testing.assert_allclose(
            roi._mesh.elm.node_number_list, sphere3_msh.elm.node_number_list
        )
        assert len(roi.get_nodes()) == 0

    @pytest.mark.slow
    def test_load_from_m2m_folder(self, example_dataset):
        m2m_path = os.path.join(example_dataset, "m2m_ernie")

        roi = RegionOfInterest()
        roi.load_mesh(subpath=m2m_path)
        ernie = mesh_io.read_msh(os.path.join(m2m_path, "ernie.msh"))
        assert not np.any(roi._mask)
        np.testing.assert_allclose(roi._mesh.nodes.node_coord, ernie.nodes.node_coord)
        np.testing.assert_allclose(
            roi._mesh.elm.node_number_list, ernie.elm.node_number_list
        )
        assert len(roi.get_nodes()) == 0

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
        assert not np.any(roi._mask)
        np.testing.assert_allclose(roi._mesh.nodes.node_coord, central.nodes.node_coord)
        np.testing.assert_allclose(
            roi._mesh.elm.node_number_list, central.elm.node_number_list
        )
        assert roi._surface_divide == lh_central.nodes.nr
        assert len(roi.get_nodes()) == 0

    def test_load_custom_surface(self, sphere3_msh: Msh, tmp_path):
        surface = sphere3_msh.crop_mesh(1004)
        surface_path = os.path.join(tmp_path, "surface.msh")
        surface.write(surface_path)

        roi = RegionOfInterest()
        roi.load_surfaces(surface_type="custom", surface_path=surface_path)
        assert not np.any(roi._mask)
        np.testing.assert_allclose(roi._mesh.nodes.node_coord, surface.nodes.node_coord)
        np.testing.assert_allclose(
            roi._mesh.elm.node_number_list, surface.elm.node_number_list
        )
        assert len(roi.get_nodes()) == 0


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
        assert len(index_mask) > 0
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
            np.array([[i + 1, i + 1, i + 1, i + 1] for i in np.unique(vertex_labels)]),
            names,
        )

        index_mask = region_of_interest.load_surface_mask_from_file(annot_file_path, 3)
        assert len(index_mask) > 0
        assert np.all(np.rint(np.sqrt(index_mask)) == 3)

    def test_load_curv(self, sphere3_msh: Msh, tmp_path):
        curv_file_path = os.path.join(tmp_path, "test")
        vertex_labels = []
        for i in range(sphere3_msh.nodes.nr):
            vertex_labels.append(np.rint(np.sqrt(i)))

        nib.freesurfer.io.write_morph_data(curv_file_path, np.array(vertex_labels))

        index_mask = region_of_interest.load_surface_mask_from_file(curv_file_path, 3)
        assert len(index_mask) > 0
        assert np.all(np.rint(np.sqrt(index_mask)) == 3)


class TestMaskImageToIndexMask:
    def test_transform_node_mask(self, sphere3_msh: Msh):
        mask_array = np.zeros((500, 500, 500))
        mask_array[:, :, :250] = 1

        mask_affine = np.eye(4)
        mask_affine[:3, 3] = [-250, -250, -250]
        mask_img = nib.Nifti1Image(mask_array, mask_affine)

        index_mask = region_of_interest.mask_image_to_index_mask(
            "node", mask_img, sphere3_msh, mask_value=1
        )

        bool_mask = np.ones((sphere3_msh.nodes.nr), dtype=np.bool_)
        bool_mask[index_mask] = False
        assert len(index_mask) > 0
        assert np.any(bool_mask)

        assert np.all(sphere3_msh.nodes.node_coord[bool_mask][:, 2] > -0.5)
        assert np.all(sphere3_msh.nodes.node_coord[index_mask][:, 2] < 0.5)

    def test_transform_element_mask(self, sphere3_msh: Msh):
        mask_array = np.zeros((500, 500, 500))
        mask_array[:, :, :250] = 1

        mask_affine = np.eye(4)
        mask_affine[:3, 3] = [-250, -250, -250]
        mask_img = nib.Nifti1Image(mask_array, mask_affine)

        index_mask = region_of_interest.mask_image_to_index_mask(
            "elm_center", mask_img, sphere3_msh, mask_value=1
        )

        bool_mask = np.ones((sphere3_msh.elm.nr), dtype=np.bool_)
        bool_mask[index_mask] = False
        assert len(index_mask) > 0
        assert np.any(bool_mask)
        assert np.all(sphere3_msh.elements_baricenters().value[bool_mask][:, 2] > -0.5)
        assert np.all(sphere3_msh.elements_baricenters().value[index_mask][:, 2] < 0.5)


class TestMniMaskToSub:
    @pytest.mark.slow
    @pytest.mark.skip(reason="Final labels MNI has artifacts")
    def test_transform_final_tissues(self, example_dataset):
        m2m_path = os.path.join(example_dataset, "m2m_ernie")
        subject_files = SubjectFiles(subpath=m2m_path)

        mask_img = nib.load(subject_files.final_labels_MNI)
        mni_to_sub_final_tissues = region_of_interest.mni_mask_to_sub(
            mask_img, subject_files
        )

        sub_final_tissues = nib.load(subject_files.final_labels)

        np.testing.assert_allclose(
            sub_final_tissues.affine, mni_to_sub_final_tissues.affine
        )

        np.testing.assert_allclose(
            np.squeeze(sub_final_tissues.get_fdata()),
            mni_to_sub_final_tissues.get_fdata(),
        )


class TestCombineMask:
    def test_union_distinct(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[:20] = True

        index_mask = np.arange(30, 50)

        result_bool_mask = region_of_interest.combine_mask(
            bool_mask, index_mask, "union"
        )

        assert np.all(result_bool_mask[:20])
        assert np.all(result_bool_mask[30:])
        assert not np.any(result_bool_mask[20:30])

    def test_union_overlap(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[10:30] = True

        index_mask = np.arange(20, 40)

        result_bool_mask = region_of_interest.combine_mask(
            bool_mask, index_mask, "union"
        )

        assert np.all(result_bool_mask[10:40])
        assert not np.any(result_bool_mask[:10])
        assert not np.any(result_bool_mask[40:])

    def test_union_equal(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[20:30] = True

        index_mask = np.arange(20, 30)

        result_bool_mask = region_of_interest.combine_mask(
            bool_mask, index_mask, "union"
        )

        assert np.all(result_bool_mask[20:30])
        assert not np.any(result_bool_mask[:20])
        assert not np.any(result_bool_mask[30:])

    def test_union_inside(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[10:40] = True

        index_mask = np.arange(20, 30)

        result_bool_mask = region_of_interest.combine_mask(
            bool_mask, index_mask, "union"
        )

        assert np.all(result_bool_mask[10:40])
        assert not np.any(result_bool_mask[:10])
        assert not np.any(result_bool_mask[40:])

    def test_intersection_distinct(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[:20] = True

        index_mask = np.arange(30, 50)

        result_bool_mask = region_of_interest.combine_mask(
            bool_mask, index_mask, "intersection"
        )

        assert not np.any(result_bool_mask)

    def test_intersection_overlap(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[10:30] = True

        index_mask = np.arange(20, 40)

        result_bool_mask = region_of_interest.combine_mask(
            bool_mask, index_mask, "intersection"
        )

        assert np.all(result_bool_mask[20:30])
        assert not np.any(result_bool_mask[:20])
        assert not np.any(result_bool_mask[30:])

    def test_intersection_equal(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[20:30] = True

        index_mask = np.arange(20, 30)

        result_bool_mask = region_of_interest.combine_mask(
            bool_mask, index_mask, "intersection"
        )

        assert np.all(result_bool_mask[20:30])
        assert not np.any(result_bool_mask[:20])
        assert not np.any(result_bool_mask[30:])

    def test_intersection_inside(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[10:40] = True

        index_mask = np.arange(20, 30)

        result_bool_mask = region_of_interest.combine_mask(
            bool_mask, index_mask, "intersection"
        )

        assert np.all(result_bool_mask[20:30])
        assert not np.any(result_bool_mask[:20])
        assert not np.any(result_bool_mask[30:])

    def test_difference_distinct(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[:20] = True

        index_mask = np.arange(30, 50)

        result_bool_mask = region_of_interest.combine_mask(
            bool_mask, index_mask, "difference"
        )

        assert np.all(result_bool_mask[:20])
        assert not np.any(result_bool_mask[20:])

    def test_difference_overlap(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[10:30] = True

        index_mask = np.arange(20, 40)

        result_bool_mask = region_of_interest.combine_mask(
            bool_mask, index_mask, "difference"
        )

        assert np.all(result_bool_mask[10:20])
        assert not np.any(result_bool_mask[:10])
        assert not np.any(result_bool_mask[20:])

    def test_difference_equal(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[20:30] = True

        index_mask = np.arange(20, 30)

        result_bool_mask = region_of_interest.combine_mask(
            bool_mask, index_mask, "difference"
        )

        assert not np.any(result_bool_mask)

    def test_difference_inside(self):
        bool_mask = np.zeros((50), dtype=np.bool_)
        bool_mask[10:40] = True

        index_mask = np.arange(20, 30)

        result_bool_mask = region_of_interest.combine_mask(
            bool_mask, index_mask, "difference"
        )

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

        result_index_mask = region_of_interest.fs_avr_mask_to_sub(
            index_mask, "lh", subject_files
        )

        central_lh = mesh_io.read_gifti_surface(
            subject_files.get_surface("lh", "central")
        )
        np.testing.assert_equal(result_index_mask, np.arange(central_lh.nodes.nr))

    @pytest.mark.slow
    def test_empty_mask(self, example_dataset):
        m2m_path = os.path.join(example_dataset, "m2m_ernie")
        subject_files = SubjectFiles(subpath=m2m_path)

        # 163842 nodes in fs avr per hemi
        index_mask = np.array([], dtype=np.int_)

        result_index_mask = region_of_interest.fs_avr_mask_to_sub(
            index_mask, "lh", subject_files
        )

        assert len(result_index_mask) == 0

    @pytest.mark.slow
    def test_half_mask(self, example_dataset):
        m2m_path = os.path.join(example_dataset, "m2m_ernie")
        subject_files = SubjectFiles(subpath=m2m_path)

        sphere: Msh = mesh_io.read_gifti_surface(
            file_finder.get_reference_surf("lh", "sphere")
        )

        # 163842 nodes in fs avr per hemi
        index_mask = np.arange(163842, dtype=np.int_)[sphere.nodes.node_coord[:, 2] > 0]

        result_index_mask = region_of_interest.fs_avr_mask_to_sub(
            index_mask, "lh", subject_files
        )

        sphere_reg_lh = mesh_io.read_gifti_surface(
            subject_files.get_surface("lh", "sphere_reg")
        )
        assert np.all(
            np.isin(
                result_index_mask,
                np.arange(sphere_reg_lh.nodes.nr)[
                    sphere_reg_lh.nodes.node_coord[:, 2] > 0
                ],
            )
        )
        assert not np.any(
            np.isin(
                result_index_mask,
                np.arange(sphere_reg_lh.nodes.nr)[
                    sphere_reg_lh.nodes.node_coord[:, 2] < 0
                ],
            )
        )


class TestApplySurfaceMask:
    @pytest.mark.slow
    def test_apply_single_minimal_surface_mask(self, example_dataset, tmp_path):
        m2m_path = os.path.join(example_dataset, "m2m_ernie")
        subject_files = SubjectFiles(subpath=m2m_path)

        central_lh = mesh_io.read_gifti_surface(
            subject_files.get_surface("lh", "central")
        )
        curv_file_path = os.path.join(tmp_path, "test")
        vertex_labels = np.zeros((central_lh.nodes.nr))
        vertex_labels[: central_lh.nodes.nr // 2] = 1
        nib.freesurfer.io.write_morph_data(curv_file_path, np.array(vertex_labels))

        roi = RegionOfInterest()
        roi.load_surfaces("central", m2m_path)
        roi.apply_surface_mask(
            "central", "subject_lh", curv_file_path, subpath=m2m_path
        )

        assert len(roi.get_nodes()) == central_lh.nodes.nr // 2
        assert np.count_nonzero(roi._mask) == central_lh.nodes.nr // 2
        assert np.all(roi._mask[: central_lh.nodes.nr // 2])

    @pytest.mark.slow
    def test_apply_multiple_full_surface_mask(self, example_dataset, tmp_path):
        m2m_path = os.path.join(example_dataset, "m2m_ernie")
        subject_files = SubjectFiles(subpath=m2m_path)

        ref_sphere_surface_lh: Msh = mesh_io.read_gifti_surface(
            file_finder.get_reference_surf("lh", "sphere")
        )

        ref_sphere_surface_rh: Msh = mesh_io.read_gifti_surface(
            file_finder.get_reference_surf("rh", "sphere")
        )

        central_lh = mesh_io.read_gifti_surface(
            subject_files.get_surface("lh", "central")
        )
        sub_lh_path = os.path.join(tmp_path, "sub_lh")
        vertex_labels = np.ones((central_lh.nodes.nr))
        vertex_labels[0] = 0
        vertex_labels[central_lh.nodes.nr // 2] = 0
        vertex_labels[-1] = 0
        nib.freesurfer.io.write_morph_data(sub_lh_path, np.array(vertex_labels))

        fs_lh_path = os.path.join(tmp_path, "fs_lh.annot")
        names = []
        # 163842 nodes in fs avr per hemi
        vertex_labels = np.zeros((163842))
        labels = [0, 1, 2, 3]
        vertex_labels[ref_sphere_surface_lh.nodes.node_coord[:, 2] > 0] = 3
        for i in labels:
            names.append(f"{i}-name")
        nib.freesurfer.io.write_annot(
            fs_lh_path,
            np.array(vertex_labels, dtype=np.int_),
            np.array([[1, 1, 1, 1], [2, 2, 2, 2], [3, 3, 3, 3], [4, 4, 4, 4]]),
            names,
        )

        central_rh = mesh_io.read_gifti_surface(
            subject_files.get_surface("rh", "central")
        )
        sub_rh_path = os.path.join(tmp_path, "sub_rh")
        vertex_labels = np.zeros((central_rh.nodes.nr))
        vertex_labels[0] = 2
        vertex_labels[central_rh.nodes.nr // 2] = 2
        vertex_labels[-1] = 2
        nib.freesurfer.io.write_morph_data(sub_rh_path, np.array(vertex_labels))

        fs_rh_path = os.path.join(tmp_path, "fs_rh.annot")
        names = []
        # 163842 nodes in fs avr per hemi
        vertex_labels = np.zeros((163842))
        vertex_labels[ref_sphere_surface_rh.nodes.node_coord[:, 2] < 0] = 4
        labels = [0, 1, 2, 3, 4]

        for i in labels:
            names.append(f"{i}-name")
        nib.freesurfer.io.write_annot(
            fs_rh_path,
            np.array(vertex_labels, dtype=np.int_),
            np.array(
                [[1, 1, 1, 1], [2, 2, 2, 2], [3, 3, 3, 3], [4, 4, 4, 4], [5, 5, 5, 5]]
            ),
            names,
        )

        sub_mask_path = os.path.join(tmp_path, "sub.nii.gz")
        mask_array = np.zeros((500, 500, 500))
        mask_array[:, 250:, :] = 5

        mask_affine = np.eye(4)
        mask_affine[:3, 3] = [-250, -250, -250]
        mask_img = nib.Nifti1Image(mask_array, mask_affine)
        nib.save(mask_img, sub_mask_path)

        roi = RegionOfInterest()
        roi.load_surfaces("central", m2m_path)
        roi.apply_surface_mask(
            "central",
            ["subject_lh", "subject_rh", "fs_avg_lh", "fs_avg_rh", "subject"],
            [
                sub_lh_path,
                sub_rh_path,
                fs_lh_path,
                fs_rh_path,
                sub_mask_path,
            ],
            [
                1,
                2,
                3,
                4,
                5,
            ],
            ["union", "union", "intersection", "union", "difference"],
            subpath=m2m_path,
        )

        joined = central_lh.join_mesh(central_rh)

        assert np.all(
            roi._mask[
                (joined.nodes.node_coord[:, 1] < -1)
                & (joined.nodes.node_coord[:, 2] > 39)
                & (joined.nodes.node_coord[:, 0] < 0)
            ]
        )
        assert (
            np.count_nonzero(
                roi._mask[
                    (joined.nodes.node_coord[:, 1] < -1)
                    & (joined.nodes.node_coord[:, 2] > 39)
                    & (joined.nodes.node_coord[:, 0] < 0)
                ]
            )
            > 0
        )

        assert np.all(
            roi._mask[
                (joined.nodes.node_coord[:, 1] < -1)
                & (joined.nodes.node_coord[:, 2] < 20)
                & (joined.nodes.node_coord[:, 0] > 3)
            ]
        )

        assert (
            np.count_nonzero(
                roi._mask[
                    (joined.nodes.node_coord[:, 1] < -1)
                    & (joined.nodes.node_coord[:, 2] < 20)
                    & (joined.nodes.node_coord[:, 0] > 3)
                ]
            )
            > 0
        )

        assert not roi._mask[0]
        assert roi._mask[-1]
        assert not roi._mask[central_lh.nodes.nr // 2]
        assert roi._mask[central_lh.nodes.nr + (central_rh.nodes.nr // 2)]

        mask = np.copy(roi._mask)

        mask[0] = False
        mask[-1] = False
        mask[central_lh.nodes.nr // 2] = False
        mask[central_lh.nodes.nr + (central_rh.nodes.nr // 2)] = False

        assert not np.any(
            mask[
                (joined.nodes.node_coord[:, 1] > 0)
                & (joined.nodes.node_coord[:, 2] < 14)
                & (joined.nodes.node_coord[:, 0] < 0)
            ]
        )

        assert (
            np.count_nonzero(
                mask[
                    (joined.nodes.node_coord[:, 1] > 0)
                    & (joined.nodes.node_coord[:, 2] < 14)
                    & (joined.nodes.node_coord[:, 0] < 0)
                ]
                == False
            )
            > 0
        )

        assert not np.any(
            mask[
                (joined.nodes.node_coord[:, 1] > 0)
                & (joined.nodes.node_coord[:, 2] > 46)
                & (joined.nodes.node_coord[:, 0] > 3)
            ]
        )

        assert (
            np.count_nonzero(
                mask[
                    (joined.nodes.node_coord[:, 1] > 0)
                    & (joined.nodes.node_coord[:, 2] > 46)
                    & (joined.nodes.node_coord[:, 0] > 3)
                ]
                == False
            )
            > 0
        )

    def test_apply_surface_mask_mismatch(self):
        roi = RegionOfInterest()
        with pytest.raises(Exception) as e_info:
            roi.apply_surface_mask(
                "central",
                mask_space=["subject_lh", "subject_lh"],
                mask_path="path/to/mask",
                mask_value=1,
                mask_operator=["union"],
                subpath="path/to/sub",
            )
        assert "same amount of elements" in str(e_info)

        with pytest.raises(Exception) as e_info:
            roi.apply_surface_mask(
                "central",
                mask_space=["subject_lh"],
                mask_path=["path/to/mask", "path/to/mask"],
                mask_value=1,
                mask_operator=["union"],
                subpath="path/to/sub",
            )
        assert "same amount of elements" in str(e_info)

        with pytest.raises(Exception) as e_info:
            roi.apply_surface_mask(
                "central",
                mask_space=["subject_lh"],
                mask_path=["path/to/mask"],
                mask_value=[1, 2],
                mask_operator=["union"],
                subpath="path/to/sub",
            )
        assert "same amount of elements" in str(e_info)

        with pytest.raises(Exception) as e_info:
            roi.apply_surface_mask(
                "central",
                mask_space=["subject_lh"],
                mask_path=["path/to/mask"],
                mask_value=1,
                mask_operator=["union", "union"],
                subpath="path/to/sub",
            )
        assert "same amount of elements" in str(e_info)


class TestApplyTissueMask:
    def test_apply_tissue_mask_union(self, sphere3_msh: Msh):
        roi = RegionOfInterest()
        roi.load_mesh(sphere3_msh)
        roi.apply_tissue_mask(ElementTags.BONE, "union")

        cropped_sphere = sphere3_msh.crop_mesh(ElementTags.BONE)

        np.testing.assert_allclose(
            roi.get_roi_mesh().nodes.node_coord, cropped_sphere.nodes.node_coord
        )
        np.testing.assert_allclose(
            roi.get_roi_mesh().elm.node_number_list, cropped_sphere.elm.node_number_list
        )

    def test_apply_tissue_mask_intersection(self, sphere3_msh: Msh):
        roi = RegionOfInterest()
        roi.load_mesh(sphere3_msh)
        roi.apply_tissue_mask([ElementTags.SCALP, ElementTags.BONE], "union")
        roi.apply_tissue_mask(ElementTags.BONE, "intersection")

        cropped_sphere = sphere3_msh.crop_mesh(ElementTags.BONE)

        np.testing.assert_allclose(
            roi.get_roi_mesh().nodes.node_coord, cropped_sphere.nodes.node_coord
        )
        np.testing.assert_allclose(
            roi.get_roi_mesh().elm.node_number_list, cropped_sphere.elm.node_number_list
        )

    def test_apply_tissue_mask_difference(self, sphere3_msh: Msh):
        roi = RegionOfInterest()
        roi.load_mesh(sphere3_msh)
        roi.apply_tissue_mask([ElementTags.SCALP, ElementTags.BONE], "union")
        roi.apply_tissue_mask(ElementTags.SCALP, "difference")

        cropped_sphere = sphere3_msh.crop_mesh(ElementTags.BONE)

        np.testing.assert_allclose(
            roi.get_roi_mesh().nodes.node_coord, cropped_sphere.nodes.node_coord
        )
        np.testing.assert_allclose(
            roi.get_roi_mesh().elm.node_number_list, cropped_sphere.elm.node_number_list
        )


class TestApplyVolumeMask:
    def test_apply_simple_volume_mask(self, sphere3_msh, tmp_path):
        roi = RegionOfInterest()
        roi.load_mesh(sphere3_msh)

        sub_mask_path = os.path.join(tmp_path, "sub.nii.gz")
        mask_array = np.zeros((500, 500, 500))
        mask_array[:, :, :250] = 1
        mask_affine = np.eye(4)
        mask_affine[:3, 3] = [-250, -250, -250]
        mask_img = nib.Nifti1Image(mask_array, mask_affine)
        nib.save(mask_img, sub_mask_path)

        roi.apply_volume_mask(mask_space="subject", mask_path=sub_mask_path)

        assert np.all(roi.get_nodes()[:, 2] < 0)
        assert len(roi.get_nodes()) > 0
        assert np.all(roi.get_roi_mesh().elements_baricenters().value[:, 2] < 0)
        assert len(roi.get_roi_mesh().elements_baricenters().value) > 0

    @pytest.mark.slow
    def test_apply_full_volume_mask(self, example_dataset, tmp_path):
        m2m_path = os.path.join(example_dataset, "m2m_ernie")

        roi = RegionOfInterest()
        roi.load_mesh(subpath=m2m_path)

        sub_mask_path = os.path.join(tmp_path, "sub.nii.gz")
        mask_array = np.ones((500, 500, 500), dtype=np.int32)
        mask_array[:, :, :250] = 42
        mask_affine = np.eye(4)
        mask_affine[:3, 3] = [-250, -250, -250]
        mask_img = nib.Nifti1Image(mask_array, mask_affine)
        nib.save(mask_img, sub_mask_path)

        mni_mask_path = os.path.join(tmp_path, "mni.nii.gz")
        mask_array = np.full((500, 500, 500), 76, dtype=np.int32)
        mask_array[:250, :, :] = 31
        mask_affine = np.eye(4)
        mask_affine[:3, 3] = [-250, -250, -250]
        mask_img = nib.Nifti1Image(mask_array, mask_affine)
        nib.save(mask_img, mni_mask_path)

        sub2_mask_path = os.path.join(tmp_path, "sub2.nii.gz")
        mask_array = np.zeros((500, 500, 500), dtype=np.int32)
        mask_array[:, :250, :] = 86
        mask_affine = np.eye(4)
        mask_affine[:3, 3] = [-250, -250, -250]
        mask_img = nib.Nifti1Image(mask_array, mask_affine)
        nib.save(mask_img, sub2_mask_path)

        roi.apply_volume_mask(
            node_type="elm_center",
            mask_space=["subject", "mni", "subject"],
            mask_path=[sub_mask_path, mni_mask_path, sub2_mask_path],
            mask_value=[42, 31, 86],
            mask_operator=["union", "intersection", "difference"],
            subpath=m2m_path,
        )

        assert np.all(roi.get_nodes()[:, 0] < 2)
        assert np.all(roi.get_nodes()[:, 1] > -1)
        assert np.all(roi.get_nodes()[:, 2] < 0)
        assert len(roi.get_nodes()) > 0

        assert np.all(roi.get_roi_mesh().elements_baricenters().value[:, 0] < 2)
        assert np.all(roi.get_roi_mesh().elements_baricenters().value[:, 1] > -1)
        assert np.all(roi.get_roi_mesh().elements_baricenters().value[:, 2] < 0)
        assert len(roi.get_roi_mesh().elements_baricenters().value) > 0

    def test_apply_surface_mask_mismatch(self):
        roi = RegionOfInterest()
        with pytest.raises(Exception) as e_info:
            roi.apply_volume_mask(
                mask_space=["subject", "subject"],
                mask_path="path/to/mask",
                mask_value=1,
                mask_operator=["union"],
            )
        assert "same amount of elements" in str(e_info)

        with pytest.raises(Exception) as e_info:
            roi.apply_volume_mask(
                mask_space=["subject"],
                mask_path=["path/to/mask", "path/to/mask"],
                mask_value=1,
                mask_operator=["union"],
            )
        assert "same amount of elements" in str(e_info)

        with pytest.raises(Exception) as e_info:
            roi.apply_volume_mask(
                mask_space=["subject"],
                mask_path=["path/to/mask"],
                mask_value=[1, 2],
                mask_operator=["union"],
            )
        assert "same amount of elements" in str(e_info)

        with pytest.raises(Exception) as e_info:
            roi.apply_volume_mask(
                mask_space=["subject"],
                mask_path=["path/to/mask"],
                mask_value=1,
                mask_operator=["union", "union"],
            )
        assert "same amount of elements" in str(e_info)


class TestApplyVolumeMaskFromSurfaceRoi:
    def test_apply_volume_mask_from_surface_roi(self, sphere3_msh: Msh, tmp_path):
        surface_roi = RegionOfInterest()
        surface_path = os.path.join(tmp_path, "surf.msh")
        sphere3_msh.crop_mesh(tags=[1003]).write(surface_path)
        surface_roi.load_surfaces("custom", surface_path=surface_path)
        surface_roi.invert()

        roi = RegionOfInterest()
        roi.load_mesh(sphere3_msh)
        roi.apply_volume_mask_from_surface_roi(surface_roi, 5, "union")

        # 1003 has distance 85 from center
        roi_distance = np.sqrt(np.sum(roi.get_nodes() ** 2, axis=1))
        assert np.all(roi_distance < 90)
        assert np.all(roi_distance > 80)


class TestApplySphereMask:
    def test_apply_simple_sphere_mask(self, sphere3_msh: Msh):
        roi = RegionOfInterest()
        roi.load_mesh(sphere3_msh)
        roi.apply_sphere_mask(
            node_type="elm_center",
            roi_sphere_center=[0, 0, 0],
            roi_sphere_radius=85,
            roi_sphere_center_space="subject",
        )

        roi_mesh = roi.get_roi_mesh()
        assert np.all(roi_mesh.elm.tag1[roi_mesh.elm.elm_type == 4] == 3)

    def test_apply_full_sphere_mask(self, sphere3_msh: Msh):
        roi = RegionOfInterest()
        roi.load_mesh(sphere3_msh)
        roi.apply_sphere_mask(
            node_type="elm_center",
            roi_sphere_center=[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
            roi_sphere_radius=[95, 90, 85],
            roi_sphere_center_space=["subject", "subject", "subject"],
            roi_sphere_operator=["union", "intersection", "difference"],
        )

        roi_mesh = roi.get_roi_mesh()
        assert np.all(roi_mesh.elm.tag1[roi_mesh.elm.elm_type == 4] == 4)


class TestToFromDict:
    def test_surface_roi_write_read_mat(self, tmp_path):
        roi = RegionOfInterest()

        roi.method = "surface"

        roi.surface_type = "central"
        roi.subpath = "path/to/m2m"

        roi.mask_space = "subject_lh"
        roi.mask_path = "path_to_file"
        roi.mask_value = 2
        roi.mask_operator = "intersection"

        roi.roi_sphere_center = [[1, 2, 3], [4, 5, 6]]
        roi.roi_sphere_radius = [3, 45]
        roi.roi_sphere_center_space = ["subject", "mni"]
        roi.roi_sphere_operator = ["union", "intersection"]

        mat_path = os.path.join(tmp_path, "test.mat")

        scipy.io.savemat(mat_path, roi.to_dict())
        roi_loaded = RegionOfInterest(dict_from_matlab(scipy.io.loadmat(mat_path)))
        assert roi.__dict__ == roi_loaded.__dict__

        scipy.io.savemat(mat_path, {'instance': roi.to_dict()})
        roi_loaded = RegionOfInterest(dict_from_matlab(scipy.io.loadmat(mat_path)['instance']))
        assert roi.__dict__ == roi_loaded.__dict__

    def test_custom_roi_write_read_mat(self, tmp_path):
        roi = RegionOfInterest()

        roi.method = "custom"

        roi.nodes = [[1.1, 2.2, 3.3], [4.4, 5.5, 6.6]]

        mat_path = os.path.join(tmp_path, "test.mat")

        scipy.io.savemat(mat_path, roi.to_dict())
        roi_loaded = RegionOfInterest(dict_from_matlab(scipy.io.loadmat(mat_path)))
        assert roi.__dict__ == roi_loaded.__dict__

        scipy.io.savemat(mat_path, {'instance': roi.to_dict()})
        roi_loaded = RegionOfInterest(dict_from_matlab(scipy.io.loadmat(mat_path)['instance']))
        assert roi.__dict__ == roi_loaded.__dict__

    def test_volume_roi_write_read_mat(self, tmp_path):
        roi = RegionOfInterest()

        roi.method = "volume"
        roi.mesh = "path/to/mesh"
        roi.tissues = [1, 2]

        roi.mask_space = ["mni", "subject"]
        roi.mask_path = ["path_to_file", "path_to_file_2"]
        roi.mask_value = [2, 33]
        roi.mask_operator = ["intersection", "difference"]

        roi.roi_sphere_center = [1, 2, 3]
        roi.roi_sphere_radius = 3
        roi.roi_sphere_center_space = "subject"
        roi.roi_sphere_operator = "union"

        mat_path = os.path.join(tmp_path, "test.mat")

        scipy.io.savemat(mat_path, roi.to_dict())
        roi_loaded = RegionOfInterest(dict_from_matlab(scipy.io.loadmat(mat_path)))
        assert roi.__dict__ == roi_loaded.__dict__

        scipy.io.savemat(mat_path, {'instance': roi.to_dict()})
        roi_loaded = RegionOfInterest(dict_from_matlab(scipy.io.loadmat(mat_path)['instance']))
        assert roi.__dict__ == roi_loaded.__dict__

    def test_surface_to_volume_roi_write_read_mat(self, tmp_path):
        roi = RegionOfInterest()

        roi.method = "surface"

        roi.surface_type = "central"
        roi.subpath = "path/to/m2m"

        roi.mask_space = ["subject_lh", "subject_rh"]
        roi.mask_path = ["path_to_file", "path_to_file_2"]
        roi.mask_value = [2, 4]
        roi.mask_operator = ["intersection", "union"]

        roi.roi_sphere_center = [1, 2, 3]
        roi.roi_sphere_radius = 3
        roi.roi_sphere_center_space = "mni"
        roi.roi_sphere_operator = "intersection"

        roi.tissues = 1
        roi.surface_inclusion_radius = 1.5

        mat_path = os.path.join(tmp_path, "test.mat")

        scipy.io.savemat(mat_path, roi.to_dict())
        roi_loaded = RegionOfInterest(dict_from_matlab(scipy.io.loadmat(mat_path)))
        assert roi.__dict__ == roi_loaded.__dict__

        scipy.io.savemat(mat_path, {'instance': roi.to_dict()})
        roi_loaded = RegionOfInterest(dict_from_matlab(scipy.io.loadmat(mat_path)['instance']))
        assert roi.__dict__ == roi_loaded.__dict__

    def test_mesh_and_mask_node_roi_write_read_mat(self, tmp_path):
        roi = RegionOfInterest()

        roi.method = "mesh+mask"
        roi.mesh = "path/to/mesh"
        roi.node_mask = [True, False, True, True, False]

        mat_path = os.path.join(tmp_path, "test.mat")

        scipy.io.savemat(mat_path, roi.to_dict())
        roi_loaded = RegionOfInterest(dict_from_matlab(scipy.io.loadmat(mat_path)))
        assert roi.__dict__ == roi_loaded.__dict__

        scipy.io.savemat(mat_path, {'instance': roi.to_dict()})
        roi_loaded = RegionOfInterest(dict_from_matlab(scipy.io.loadmat(mat_path)['instance']))
        assert roi.__dict__ == roi_loaded.__dict__

    def test_mesh_and_mask_elm_roi_write_read_mat(self, tmp_path):
        roi = RegionOfInterest()

        roi.method = "mesh+mask"
        roi.subpath = "path/to/m2m"
        roi.elm_mask = [True, False, True, True, False, False, False, True]

        mat_path = os.path.join(tmp_path, "test.mat")

        scipy.io.savemat(mat_path, roi.to_dict())
        roi_loaded = RegionOfInterest(dict_from_matlab(scipy.io.loadmat(mat_path)))
        assert roi.__dict__ == roi_loaded.__dict__

        scipy.io.savemat(mat_path, {'instance': roi.to_dict()})
        roi_loaded = RegionOfInterest(dict_from_matlab(scipy.io.loadmat(mat_path)['instance']))
        assert roi.__dict__ == roi_loaded.__dict__
