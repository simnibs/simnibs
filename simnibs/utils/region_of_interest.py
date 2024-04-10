from copy import deepcopy
import os
import numpy as np
import nibabel as nib
import scipy

from simnibs.mesh_tools import mesh_io
from simnibs.mesh_tools.mesh_io import Elements, Msh, Nodes

import simnibs.utils.file_finder as file_finder
from simnibs.utils.mesh_element_properties import ElementTags
from .file_finder import SubjectFiles
from .transformations import (
    mni2subject_coords,
)
from simnibs.utils import transformations


class RegionOfInterest:
    method: str  # ("manual", "custom", "surface", "volume", "volume_from_surface", "mesh+mask")

    subpath: str | None
    mesh: str | Msh | None

    mask_space: (
        str | list[str] | None
    )  # method = "surface" : ("subject", "subject_lh", "fs_avg_lh", "subject_rh", "fs_avg_rh", "mni") | method = "volume" : ("subject", "mni")
    mask_path: (
        str | list[str] | None
    )  # method = "surface" : (label, annot, curv, nifti) | method = "volume" : (nifti)
    mask_value: int | list[int] | None  # default 1
    mask_operator: str | list[str] | None  # default "union" ("union", "intersection", "difference")

    roi_sphere_center: list[float] | list[list[float]] | None
    roi_sphere_radius: float | list[float] | None
    roi_sphere_center_space: str | list[str] | None  # ("subject", "mni")
    roi_sphere_operator: str | list[str] | None # default "union" ("union", "intersection", "difference")

    # method = "custom"
    nodes: list[list[float]] | None

    # method = "surface"
    surface_type: str | None  # ("central", "custom")
    surface_path: str | None  # surface_type = "custom" : (msh, freesurfer, gifti)

    # method = "volume"
    tissues: int | list[int] | ElementTags | list[ElementTags] | None  # default 2

    # method = "volume_from_surface"
    surface_inclusion_radius: float | None

    # method = "mesh+mask"
    node_mask: list[bool] | None  # exclusive -> nodes
    elm_mask: list[bool] | None  # exclusive -> elm_center

    def __init__(self):
        self.method = "manual"

        self.subpath = None
        self.mesh = None

        self.mask_space = None
        self.mask_path = None
        self.mask_value = None
        self.mask_operator = None

        self.roi_sphere_center = None
        self.roi_sphere_radius = None
        self.roi_sphere_center_space = None
        self.roi_sphere_operator = None

        self.nodes = None

        self.surface_type = None
        self.surface_path = None

        self.tissues = None

        self.surface_inclusion_radius = None

        self.node_mask = None
        self.elm_mask = None

        self._prepared = False

    def _prepare(self):

        match self.method:
            case "custom":
                self._mesh = Msh(Nodes(np.array(self.nodes)), Elements())
                self._mask = np.ones((self._mesh.nodes.nr), dtype=np.bool_)
                self._mask_type = "node"

            case "surface":
                self.load_surfaces(self.surface_type, self.subpath, self.surface_path)
                self.apply_surface_mask(
                    self.surface_type,
                    self.mask_space,
                    self.mask_path,
                    self.mask_value,
                    self.mask_operator,
                    self.subpath,
                )
                self.apply_sphere_mask(
                    "node",
                    self.roi_sphere_center,
                    self.roi_sphere_radius,
                    self.roi_sphere_center_space,
                    self.roi_sphere_operator,
                    self.subpath,
                )
            case "volume":
                self.load_mesh(self.mesh, self.subpath)
                self.apply_tissue_mask(self.tissues)
                self.apply_volume_mask(
                    "elm_center",
                    self.mask_space,
                    self.mask_path,
                    self.mask_value,
                    self.mask_operator,
                    self.subpath,
                )
                self.apply_sphere_mask(
                    "elm_center",
                    self.roi_sphere_center,
                    self.roi_sphere_radius,
                    self.roi_sphere_center_space,
                    self.roi_sphere_operator,
                    self.subpath,
                )

            case "volume_from_surface":
                surface_roi = RegionOfInterest()
                surface_roi.load_surfaces(
                    self.surface_type, self.subpath, self.surface_path
                )
                surface_roi.apply_surface_mask(
                    self.surface_type,
                    self.mask_space,
                    self.mask_path,
                    self.mask_value,
                    self.mask_operator,
                    self.subpath,
                )
                surface_roi.apply_sphere_mask(
                    "node",
                    self.roi_sphere_center,
                    self.roi_sphere_radius,
                    self.roi_sphere_center_space,
                    self.roi_sphere_operator,
                    self.subpath,
                )

                self.load_mesh(self.mesh, self.subpath)
                self.apply_tissue_mask(self.tissues)
                self.apply_volume_mask_from_surface_roi(
                    surface_roi, self.surface_inclusion_radius
                )
            case "mesh+mask":
                self.load_mesh(self.mesh, self.subpath)
                if self.node_mask is not None and self.elm_mask is None:
                    self._mask = self.node_mask
                    self._mask_type = "node"
                elif self.elm_mask is not None and self.node_mask is None:
                    self._mask = self.elm_mask
                    self._mask_type = "elm_center"
                else:
                    raise ValueError(
                        f'For method = "mesh+mask" either node_mask or elm_mask need to be set'
                    )
            case _:
                raise ValueError(
                    f'surface_type needs to be one of ["custom", "surface", "volume", "volume_from_surface", "mesh+mask"] (was {self.method})'
                )
        self._prepared = True

    def to_mat(self):
        return {'a': 'b'}

    def get_nodes(self):
        if self.method != "manual" and not self._prepared:
            self._prepare()

        match self._mask_type:
            case "node":
                return self._mesh.nodes.node_coord[self._mask]
            case "elm_center":
                return self._mesh.elements_baricenters().value[self._mask]
            case _:
                raise ValueError(f"No mesh or surface was loaded")

    def get_roi_mesh(self):
        if self.method != "manual" and not self._prepared:
            self._prepare()

        if self.method == "custom":
            return self._mesh

        match self._mask_type:
            case "node":
                node_indexes = self._mesh.nodes.node_number[self._mask]
                if len(node_indexes) == 0:
                    return Msh()
                return self._mesh.crop_mesh(nodes=node_indexes)
            case "elm_center":
                elm_indexes = self._mesh.elm.elm_number[self._mask]
                if len(elm_indexes) == 0:
                    return Msh()
                return self._mesh.crop_mesh(elements=elm_indexes)
            case _:
                raise ValueError(f"No mesh or surface was loaded")

    def write_visualization(self, folder_path: str, base_file_name: str):
        if self.method != "manual" and not self._prepared:
            self._prepare()

        geo_file_path = os.path.join(folder_path, f"{base_file_name}.geo")
        if os.path.isfile(geo_file_path):
            os.remove(geo_file_path)

        roi_colormap = np.zeros((255, 4), dtype=int)
        roi_colormap[:-1, :] = [100, 100, 255, 255]
        roi_colormap[-1, :] = [255, 0, 0, 255]

        if self._mesh.elm.nr > 0:
            match self._mask_type:
                case "node":
                    v = self._mesh.view()
                    roi_node_data = np.zeros(self._mesh.nodes.nr)
                    roi_node_data[self._mask] = 1
                    self._mesh.add_node_field(roi_node_data, "ROI")
                    v.add_view(
                        ColorTable=roi_colormap,
                        Visible=1,
                        ShowScale=0,
                        CustomMin=0,
                        CustomMax=1,
                        RangeType=2,
                    )
                case "elm_center":
                    v = self._mesh.view(visible_tags=np.unique(self._mesh.elm.tag1[self._mask]))
                    elm_node_data = np.zeros(self._mesh.elm.nr)
                    elm_node_data[self._mask] = 1
                    self._mesh.add_element_field(elm_node_data, "ROI")

                    v.add_view(
                        ColorTable=roi_colormap,
                        Visible=1,
                        ShowScale=0,
                        CustomMin=0.9,
                        CustomMax=1,
                        RangeType=2,
                    )
        else:
            v = self._mesh.view()

        v.Mesh.SurfaceFaces = 0
        v.Mesh.VolumeFaces = 0

        nodes = self.get_nodes()
        mesh_io.write_geo_spheres(
            nodes,
            geo_file_path,
            np.full((len(nodes)), 1),
            name=f"roi-nodes",
            mode="ba",
        )
        node_view = v.add_view(
            ShowScale=0, PointType=0, PointSize=3.0, ColormapNumber=2
        )
        if self._mesh.elm.nr == 0:
            node_view.Visible = 1

        if ElementTags.SCALP_TH_SURFACE in self._mesh.elm.tag1:
            skin_mesh = self._mesh.crop_mesh(tags=[ElementTags.SCALP_TH_SURFACE])
            mesh_io.write_geo_triangles(
                skin_mesh.elm.node_number_list - 1,
                skin_mesh.nodes.node_coord,
                geo_file_path,
                name="scalp",
                mode="ba",
            )

            v.add_view(
                ColormapNumber=8, ColormapAlpha=0.3, Visible=1, ShowScale=0
            )  # scalp

        v.add_merge(geo_file_path)
        v.mesh.write(os.path.join(folder_path, f"{base_file_name}.msh"))
        v.write_opt(os.path.join(folder_path, f"{base_file_name}.msh"))

        match self._mask_type:
            case "node":
                del self._mesh.nodedata[-1]
            case "elm_center":
                del self._mesh.elmdata[-1]

    def load_surfaces(
        self, surface_type: str | None, subpath: str | None, surface_path: str | None
    ):
        surfaces: list[Msh] = []
        match surface_type:
            case "central":
                if subpath is None:
                    raise ValueError(
                        f'If surface_type = "central", subpath needs to be set (was {subpath})'
                    )
                subject_files = _init_subject_files(subpath)
                surfaces.append(
                    mesh_io.read_gifti_surface(
                        subject_files.get_surface("lh", "central")
                    )
                )
                surfaces.append(
                    mesh_io.read_gifti_surface(
                        subject_files.get_surface("rh", "central")
                    )
                )

            case "custom":
                if surface_path is None:
                    raise ValueError(
                        f'If surface_type = "custom", surface_path needs to be set (was {surface_path})'
                    )
                surfaces.append(load_surface_from_file(surface_path))
            case _:
                raise ValueError(
                    f'surface_type needs to be one of ["central", "custom"] (was {surface_type})'
                )
        surface = surfaces[0]
        if len(surfaces) == 2:
            self._surface_divide = surface.nodes.nr
            surface = surface.join_mesh(surfaces[1])
        self._mesh = surface
        self._mask_type = "node"
        self._mask = np.ones((self._mesh.nodes.nr), dtype=np.bool_)

    def load_mesh(self, mesh: str | Msh | None, subpath: str | None):
        if mesh is not None:
            if isinstance(mesh, Msh):
                self._mesh = mesh
            else:
                self._mesh = mesh_io.read_msh(mesh)

        if subpath is not None:
            self._mesh = mesh_io.read_msh(_init_subject_files(subpath).fnamehead)

        if self._mesh is None:
            raise ValueError(f"mesh or subpath needs to be set (was {mesh}, {subpath})")

        self._mask_type = "elm_center"
        self._mask = np.ones((self._mesh.elm.nr), dtype=np.bool_)

    def apply_surface_mask(
        self,
        surface_type: str,
        mask_space: str | list[str] | None,
        mask_path: str | list[str] | None,
        mask_value: int | list[int] | None,
        mask_operator: str | list[str] | None,
        subpath: str | None,
    ):
        if mask_path is not None or mask_space is not None:
            if mask_path is None:
                raise ValueError(f"mask_path needs to be set (was {mask_space})")

            if not isinstance(mask_path, list):
                mask_path = [mask_path]

            if mask_space is None:
                raise ValueError(
                    f'elements of mask_space needs to be one of ["subject_lh", "fs_avg_lh", "subject_rh", "fs_avg_rh", "subject", "mni"] (was {mask_space})'
                )

            if not isinstance(mask_space, list):
                mask_space = [mask_space]

            if mask_value is None:
                mask_value = 1

            if isinstance(mask_value, int):
                mask_value = [mask_value] * len(mask_path)

            if mask_operator is None:
                mask_operator = "intersection"

            if isinstance(mask_operator, str):
                mask_operator = [mask_operator] * len(mask_path)

            if not (
                len(mask_path)
                == len(mask_space)
                == len(mask_value)
                == len(mask_operator)
            ):
                raise ValueError(
                    "mask_path, mask_space, mask_value and mask_operator need the same amount of elements"
                )

            for mask_path, mask_space, mask_value, mask_operator in zip(
                mask_path, mask_space, mask_value, mask_operator
            ):
                if os.path.splitext(mask_path)[1] in ["nii", "nii.gz"]:
                    self.apply_volume_mask(
                        "node",
                        mask_space,
                        mask_path,
                        mask_value,
                        mask_operator,
                        subpath,
                    )
                else:
                    index_mask = load_surface_mask_from_file(mask_path, mask_value)
                    if surface_type == "custom":
                        match mask_space:
                            case "subject":
                                self._mask = combine_mask(
                                    self._mask, index_mask, mask_operator
                                )
                            case _:
                                raise ValueError(
                                    f'mask_space needs to be one of ["subject", "mni"] for a custom surface (was {mask_space})'
                                )
                    else:
                        match mask_space:
                            case "subject_lh":
                                self._mask[: self._surface_divide] = combine_mask(
                                    self._mask[: self._surface_divide],
                                    index_mask,
                                    mask_operator,
                                )
                            case "subject_rh":
                                self._mask[self._surface_divide :] = combine_mask(
                                    self._mask[self._surface_divide :],
                                    index_mask,
                                    mask_operator,
                                )
                            case "fs_avg_lh":
                                index_mask = fs_avr_mask_to_sub(
                                    index_mask, "lh", _init_subject_files(subpath)
                                )
                                self._mask[: self._surface_divide] = combine_mask(
                                    self._mask[: self._surface_divide],
                                    index_mask,
                                    mask_operator,
                                )
                            case "fs_avg_rh":
                                index_mask = fs_avr_mask_to_sub(
                                    index_mask, "rh", _init_subject_files(subpath)
                                )
                                self._mask[self._surface_divide :] = combine_mask(
                                    self._mask[self._surface_divide :],
                                    index_mask,
                                    mask_operator,
                                )
                            case _:
                                raise ValueError(
                                    f'elements of mask_space need to be one of ["subject_lh", "fs_avg_lh", "subject_rh", "fs_avg_rh", "mni"] (was {mask_space})'
                                )

    def apply_volume_mask(
        self,
        node_type: str,
        mask_space: str | list[str] | None,
        mask_path: str | list[str] | None,
        mask_value: int | list[int] | None,
        mask_operator: str | list[str] | None,
        subpath: str | None,
    ):
        if mask_path is not None or mask_space is not None:
            if mask_path is None:
                raise ValueError(f"mask_path needs to be set (was {mask_space})")

            if not isinstance(mask_path, list):
                mask_path = [mask_path]

            if mask_space is None:
                raise ValueError(
                    f'elements of mask_space needs to be one of ["subject", "mni"] (was {mask_space})'
                )

            if not isinstance(mask_space, list):
                mask_space = [mask_space]

            if mask_value is None:
                mask_value = 1

            if isinstance(mask_value, int):
                mask_value = [mask_value] * len(mask_path)

            if mask_operator is None:
                mask_operator = "intersection"

            if isinstance(mask_operator, str):
                mask_operator = [mask_operator] * len(mask_path)

            if not (
                len(mask_path)
                == len(mask_space)
                == len(mask_value)
                == len(mask_operator)
            ):
                raise ValueError(
                    "mask_path, mask_space, mask_value and mask_operator need the same amount of elements"
                )

            for mask_path, mask_space, mask_value, mask_operator in zip(
                mask_path, mask_space, mask_value, mask_operator
            ):
                mask_img = nib.load(mask_path)
                match mask_space:
                    case "subject":
                        index_mask = mask_image_to_index_mask(
                            node_type, mask_img, self._mesh, mask_value
                        )
                        self._mask = combine_mask(self._mask, index_mask, mask_operator)
                    case "mni":
                        mask_img = mni_mask_to_sub(
                            mask_img, _init_subject_files(subpath)
                        )
                        index_mask = mask_image_to_index_mask(
                            node_type, mask_img, self._mesh, mask_value
                        )
                        self._mask = combine_mask(self._mask, index_mask, mask_operator)
                    case _:
                        raise ValueError(
                            f'mask_space needs to be one of ["subject", "mni"] for a NIfTI file (was {mask_space})'
                        )

    def apply_volume_mask_from_surface_roi(
        self,
        surface_roi: "RegionOfInterest",
        surface_inclusion_radius: float | None,
        surface_roi_operator: str | None = None,
    ):
        if surface_inclusion_radius is None:
            raise ValueError(
                f"surface_inclusion_radius needs to be set (was {surface_inclusion_radius})"
            )

        if surface_roi_operator is None:
            surface_roi_operator = "intersection"

        kd_tree = scipy.spatial.cKDTree(self._mesh.elements_baricenters().value)
        index_mask = kd_tree.query_ball_point(
            surface_roi.get_nodes(), surface_inclusion_radius
        )
        index_mask = np.unique(np.concatenate(index_mask))
        self._mask = combine_mask(self._mask, index_mask, surface_roi_operator)

    def apply_sphere_mask(
        self,
        node_type: str,
        roi_sphere_center: list[float] | list[list[float]] | None,
        roi_sphere_radius: float | list[float] | None,
        roi_sphere_center_space: str | list[str] | None,
        roi_sphere_operator: str | list[str] | None,
        subpath: str | None,
    ):
        if (
            roi_sphere_center is not None
            or roi_sphere_radius
            or roi_sphere_center_space is not None
        ):
            if roi_sphere_center is None:
                raise ValueError(
                    f"roi_sphere_center needs to be set (was {roi_sphere_center})"
                )

            if not isinstance(roi_sphere_center[0], list):
                roi_sphere_center = [roi_sphere_center]

            if roi_sphere_radius is None:
                raise ValueError(
                    f"roi_sphere_radius needs to be set (was {roi_sphere_radius})"
                )

            if not isinstance(roi_sphere_radius, list):
                roi_sphere_radius = [roi_sphere_radius]

            if roi_sphere_center_space is None:
                raise ValueError(
                    f"roi_sphere_center_space needs to be set (was {roi_sphere_center_space})"
                )

            if not isinstance(roi_sphere_center_space, list):
                roi_sphere_center_space = [roi_sphere_center_space]

            if roi_sphere_operator is None:
                roi_sphere_operator = "intersection"

            if isinstance(roi_sphere_operator, str):
                roi_sphere_operator = [roi_sphere_operator] * len(roi_sphere_center)

            if not (
                len(roi_sphere_center)
                == len(roi_sphere_radius)
                == len(roi_sphere_center_space)
                == len(roi_sphere_operator)
            ):
                raise ValueError(
                    "roi_sphere_center, roi_sphere_radius, roi_sphere_center_space and roi_sphere_operator need the same amount of elements"
                )

            match node_type:
                case "node":
                    kd_tree = scipy.spatial.cKDTree(self._mesh.nodes.node_coord)
                case "elm_center":
                    kd_tree = scipy.spatial.cKDTree(
                        self._mesh.elements_baricenters().value
                    )
                case _:
                    raise ValueError(
                        f'node_type needs to be one of ["node", "elm_center"] (was {node_type})'
                    )

            for (
                roi_sphere_center,
                roi_sphere_radius,
                roi_sphere_center_space,
                roi_sphere_operator,
            ) in zip(
                roi_sphere_center,
                roi_sphere_radius,
                roi_sphere_center_space,
                roi_sphere_operator,
            ):
                match roi_sphere_center_space:
                    case "subject":
                        center = roi_sphere_center
                    case "mni":
                        center = mni2subject_coords(
                            roi_sphere_center, _init_subject_files(subpath).subpath
                        )
                    case _:
                        raise ValueError(
                            f'elements of roi_sphere_center_space need to be one of ["subject", "mni"] (was {roi_sphere_center_space})'
                        )
                index_mask = kd_tree.query_ball_point(center, roi_sphere_radius)
                self._mask = combine_mask(self._mask, index_mask, roi_sphere_operator)

    def apply_tissue_mask(
        self,
        tissues: int | list[int] | ElementTags | list[ElementTags] | None,
        tissue_mask_operator: str | None = None,
    ):
        if tissues is None:
            tissues = 2

        if not isinstance(tissues, list):
            tissues = [tissues]

        if tissue_mask_operator is None:
            tissue_mask_operator = "intersection"

        tissues = np.array(tissues)

        self._mask = combine_mask(
            self._mask,
            self._mesh.elm.elm_number[np.in1d(self._mesh.elm.tag1, tissues)] - 1,
            tissue_mask_operator,
        )


def _init_subject_files(subpath: str | None):
    if subpath is None:
        raise ValueError(f"subpath needs to be set (was {subpath})")
    subpath = os.path.abspath(os.path.expanduser(subpath))
    if not os.path.isdir(subpath):
        raise IOError(f"Cannot locate subjects m2m folder: {subpath}")

    return SubjectFiles(subpath=subpath)


def load_surface_from_file(surface_path: str) -> Msh:
    if not os.path.isfile(surface_path):
        raise ValueError("surface_path needs to be a file")

    _, file_extension = os.path.splitext(surface_path)
    surface = None
    match file_extension:
        case "msh":
            surface: Msh = mesh_io.read_msh(surface_path)
            if np.any(surface.elm.elm_type != 2):
                raise ValueError(".msh file contains non triangle elements")
        case "gii":
            surface = mesh_io.read_gifti_surface(surface_path)
        case _:
            try:
                surface = mesh_io.read_freesurfer_surface(surface_path, True)
            except:
                pass
    if surface is None:
        raise ValueError(
            "surface_path needs to be in one of the following file formats: [msh, Gifti, FreeSurfer surface]"
        )
    return surface


def load_surface_mask_from_file(mask_path: str, mask_value: int):
    if not os.path.isfile(mask_path):
        raise ValueError("surface_path needs to be a file")

    _, file_extension = os.path.splitext(mask_path)
    index_mask = None
    match file_extension:
        case "label":
            index_mask = nib.freesurfer.io.read_label(mask_path)
        case "annot":
            labels, _, _ = nib.freesurfer.io.read_annot(mask_path)
            index_mask = np.where(labels == mask_value)[0]
        case _:
            try:
                data = nib.freesurfer.io.read_morph_data(mask_path)
                index_mask = np.where(np.rint(data) == mask_value)[0]
            except:
                pass

    if index_mask is None:
        ValueError(
            "mask_path needs to be in one of the following file formats: [FreeSurfer label, FreeSurfer annot, FreeSurfer curv, NIfTI]"
        )

    return index_mask


def fs_avr_mask_to_sub(index_mask, hemi, subject_files: SubjectFiles):
    sphere_surface: Msh = mesh_io.read_gifti_surface(
        file_finder.get_reference_surf(hemi, "sphere")
    )
    registration_surface: Msh = mesh_io.read_gifti_surface(
        subject_files.get_surface(hemi, "sphere_reg")
    )

    morph = transformations.SurfaceMorph(
        sphere_surface, registration_surface, method="nearest"
    )

    index_mask = morph.transform(index_mask) > 0.0001

    return index_mask


def combine_mask(bool_mask, index_mask, operator):
    match operator:
        case "union":
            bool_mask[index_mask] = True
        case "intersection":
            bool_index_mask = np.zeros_like(bool_mask, dtype=np.bool_)
            bool_index_mask[index_mask] = True
            bool_mask = bool_mask & bool_index_mask
        case "difference":
            bool_index_mask = np.ones_like(bool_mask, dtype=np.bool_)
            bool_index_mask[index_mask] = False
            bool_mask = bool_mask & bool_index_mask
        case _:
            raise ValueError('elements of mask_space needs to be one of ["union", "intersection", "difference"]')

    return bool_mask


def mni_mask_to_sub(mask, subject_files: SubjectFiles):
    target_image = nib.load(subject_files.T1_upsampled)
    target_dim = list(target_im.get_fdata().shape)

    image_deformation = nib.load(subject_files.conf2mni_nonl)

    transformed_mask = transformations.volumetric_nonlinear(
        (mask.get_fdata(), mask.affine),
        (image_deformation.get_fdata(), image_deformation.affine),
        target_space_affine=target_image.affine,
        target_dimensions=target_dim,
        intorder=0,
    )
    transformed_mask = np.squeeze(transformed_mask)

    return nib.Nifti1Image(transformed_mask, target_image.affine)


def mask_image_to_index_mask(node_type: str, mask_img, surface, mask_value):
    match node_type:
        case "node":
            data = mesh_io.NodeData.from_data_grid(
                surface, mask_img.get_fdata(), mask_img.affine, "", order=0
            ).value
        case "elm_center":
            data = mesh_io.ElementData.from_data_grid(
                surface, mask_img.get_fdata(), mask_img.affine, "", order=0
            ).value
        case _:
            raise ValueError(
                f'node_type needs to be one of ["node", "elm_center"] (was {node_type})'
            )

    return np.where(np.rint(data) == mask_value)[0]
