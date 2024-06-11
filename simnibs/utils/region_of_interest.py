import os
from typing import get_type_hints
import numpy as np
import numpy.typing as npt
import nibabel as nib
import scipy

from simnibs.mesh_tools import mesh_io
from simnibs.mesh_tools.mesh_io import Elements, Msh, Nodes

import simnibs.utils.file_finder as file_finder
from simnibs.utils.matlab_read import matlab_field_to_list, remove_None, try_to_read_matlab_field
from simnibs.utils.mesh_element_properties import ElementTags
from .file_finder import SubjectFiles
from .transformations import (
    mni2subject_coords,
)
from simnibs.utils import transformations


class RegionOfInterest:
    """A class describing a region of Interest in a volume mesh or a surface mesh.
    
    Parameters
    ------------------------
    matlab_struct: (optional) scipy.io.loadmat()
        matlab structure
    """

    method: str
    """ The method to create the ROI {"manual", "custom", "surface", "volume", "volume_from_surface", "mesh+mask"} """

    subpath: str | None
    """ Path to the m2m folder of the subject (example: "path/to/m2m")"""
    mesh: str | Msh | None
    """ Path to a mesh or a mesh instance (example: "path/to/msh" | mesh_instance)"""

    mask_space: str | list[str] | None
    """ The space the mask is defined in, method = "surface" : {"subject", "subject_lh", "fs_avg_lh", "subject_rh", "fs_avg_rh", "mni"} | method = "volume" : {"subject", "mni"} """
    mask_path: str | list[str] | None
    """The path to the mask, method = "surface" : (label, annot, curv, nifti) | method = "volume" : (nifti) (example: "path/to/file")"""
    mask_value: int | list[int] | None
    """ The values that are considered as the mask in the mask files, default 1 (example: 1 | [1, 2])"""
    mask_operator: str | list[str] | None
    """The operator to combine the mask with the ROI {"union", "intersection", "difference"}, default "union" """

    roi_sphere_center: list[float] | list[list[float]] | None
    """ Sphere center coordinates for spherical masks in mm (example: [0,0,0] | [[1,2,3], [4,5,6]]) """
    roi_sphere_radius: float | list[float] | None
    """ The radius of the spherical masks in mm (example: 5 | [3, 45])"""
    roi_sphere_center_space: str | list[str] | None
    """ The space the center coordinates of the spheres are defined in {"subject", "mni"} """

    roi_sphere_operator: str | list[str] | None
    """The operator to combine the mask with the ROI {"union", "intersection", "difference"}, default "union" """

    nodes: list[list[float]] | None
    """ Only for method = "custom" -> a custom list of node coordinates (example: [[1.1, 2.2, 3.3], [4.4, 5.5, 6.6]])"""

    surface_type: str | None
    """ Only for method = "surface" -> Weather the surface is the subject specific central gm surface or a custom surface {"central", "custom"} """
    surface_path: str | None
    """ Only for method = "surface" -> Only for surface_type = "custom" -> The path to a custom surface (msh, freesurfer, gifti) """

    tissues: int | list[int] | ElementTags | list[ElementTags] | None
    """ Only for method = "volume" -> a number of volume tissue tags, default 2 (example: ElementTags.GM | [ElementTags.WM, ElementTags.GM]) """

    surface_inclusion_radius: float | None
    """ Only for method = "volume_from_surface" -> The radius from the surface nodes at which the volume elements should be included in mm (example: 5)"""

    node_mask: list[bool] | None  
    """ Only for method = "mesh+mask" -> a boolean node mask (exclusive with elm_mask) (example: [True, ..., False])"""
    elm_mask: list[bool] | None  
    """ Only for method = "mesh+mask" -> a boolean node mask (exclusive with node_mask) (example: [True, ..., False])"""

    def __init__(self, matlab_struct=None):
        self.method = "manual"
        self.subpath = None
        self.mesh = None
        self.fname_visu = None

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

        if matlab_struct:
           self.read_mat_struct(matlab_struct)

    def _prepare(self):
        """Prepares the Region of Interest based on the scripting parameters of the class.

        Raises
        ------
        ValueError
            If method = "mesh+mask" and node_mask and elm_mask are None
        ValueError
            method not in {"custom", "surface", "volume", "volume_from_surface", "mesh+mask"}
        """

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
                    self._mask = np.array(self.node_mask)
                    self._mask_type = "node"
                elif self.elm_mask is not None and self.node_mask is None:
                    self._mask = np.array(self.elm_mask)
                    self._mask_type = "elm_center"
                else:
                    raise ValueError(
                        'For method = "mesh+mask" either node_mask or elm_mask need to be set'
                    )
            case _:
                raise ValueError(
                    f'method needs to be one of ["custom", "surface", "volume", "volume_from_surface", "mesh+mask"] (was {self.method})'
                )
        self._prepared = True

    def to_mat(self):
        """ Makes a dictionary for saving a matlab structure with scipy.io.savemat()

        Returns
        --------------------
        dict
            Dictionaty for usage with scipy.io.savemat
        """
        # Generate dict from instance variables (excluding variables starting with _ or __)
        mat = {
            key:remove_None(value) for key, value in self.__dict__.items() 
            if not key.startswith('__')
            and not key.startswith('_')
            and not callable(value) 
            and not callable(getattr(value, "__get__", None))
        }

        # Add class name as type (type is protected in python so it cannot be a instance variable)
        mat['type'] = 'RegionOfInterest'

        return mat
    
    def read_mat_struct(self, mat):
        """ Reads parameters from matlab structure

        Parameters
        ----------
        mat: scipy.io.loadmat
            Loaded matlab structure
        """
        self.method = try_to_read_matlab_field(mat, 'method', str)
        self.subpath = try_to_read_matlab_field(mat, 'subpath', str)
        self.mesh = try_to_read_matlab_field(mat, 'mesh', str)
        self.fname_visu = try_to_read_matlab_field(mat, 'fname_visu', str)
        
        self.mask_space = matlab_field_to_list(mat, 'mask_space', 1) 
        self.mask_path = matlab_field_to_list(mat, 'mask_path', 1) 
        self.mask_value = matlab_field_to_list(mat, 'mask_value', 1)
        self.mask_operator = matlab_field_to_list(mat, 'mask_operator', 1)

        self.roi_sphere_center = matlab_field_to_list(mat, 'roi_sphere_center', 2)
        self.roi_sphere_radius = matlab_field_to_list(mat, 'roi_sphere_radius', 1)
        self.roi_sphere_center_space = matlab_field_to_list(mat, 'roi_sphere_center_space', 1)
        self.roi_sphere_operator = matlab_field_to_list(mat, 'roi_sphere_operator', 1)

        self.nodes = matlab_field_to_list(mat, 'nodes', 2)

        self.surface_type = try_to_read_matlab_field(mat, 'surface_type', str)
        self.surface_path = try_to_read_matlab_field(mat, 'surface_path', str)

        self.tissues = matlab_field_to_list(mat, 'tissues', 1)

        self.surface_inclusion_radius = try_to_read_matlab_field(mat, 'surface_inclusion_radius', float)

        self.node_mask = matlab_field_to_list(mat, 'node_mask', 1)
        self.elm_mask = matlab_field_to_list(mat, 'elm_mask', 1)

    def get_nodes(self, node_type: str | None = None) -> npt.NDArray[np.float_]:
        """Returns the nodes which are part of the Region of Interest.
        Element center coordinates in the case of node_type="elm_center", node coordinates in the case of node_type="node"

        Parameters
        ----------
        node_type : str | None, optional {"node", "elm_center"}
            Weather to return element center coordinates or node coordinates, by default elm_center if volume was loaded, node if surface was loaded

        Returns
        -------
        npt.NDArray[np.float_]
            The node coordinates which are part of the Region of Interest

        Raises
        ------
        ValueError
            If no mesh or surface was loaded
        ValueError
            node_type not in {"node", "elm_center"}
        """
        if self.method != "manual" and not self._prepared:
            self._prepare()

        if node_type is None:
            match self._mask_type:
                case "node":
                    return self._mesh.nodes.node_coord[self._mask]
                case "elm_center":
                    return self._mesh.elements_baricenters().value[self._mask]
                case _:
                    raise ValueError("No mesh or surface was loaded")
        else:
            roi_mesh = self.get_roi_mesh()
            match node_type:
                case "node":
                    return roi_mesh.nodes.node_coord
                case "elm_center":
                    return roi_mesh.elements_baricenters().value
                case _:
                    raise ValueError(
                        f'node_type needs to be one of ["node", "elm_center"] (was {node_type})'
                    )

    def get_roi_mesh(self) -> Msh:
        """Returns the Region of Interest as a mesh

        Returns
        -------
        Msh
            The Region of Interest as a mesh

        Raises
        ------
        ValueError
            If no mesh or surface was loaded
        """
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
                raise ValueError("No mesh or surface was loaded")

    def write_visualization(self, folder_path: str, base_file_name: str):
        """Writes a visualization of the Region of Interest to a folder

        Parameters
        ----------
        folder_path : str
            Folder to write the visualization to
        base_file_name : str
            The base file name of the visualization files
        """
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
                    v = self._mesh.view(
                        visible_tags=np.unique(self._mesh.elm.tag1[self._mask])
                    )
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
            name="roi-nodes",
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
                if len(self._mesh.nodedata):
                    del self._mesh.nodedata[-1]
            case "elm_center":
                del self._mesh.elmdata[-1]

    def load_surfaces(
        self,
        surface_type: str | None = None,
        subpath: str | None = None,
        surface_path: str | None = None,
    ):
        """Initializes the Region of Interest with a surface either with a surface from a m2m folder or a custom surface.
        The surface will then be subject to masking to define a Region of Interest in the surface.

        Parameters
        ----------
        surface_type : str | None, optional {"central", "custom"}
            Weather to load a subject specific central surface or a custom surface, by default None
        subpath : str | None, optional
            The m2m folder to load surfaces from, by default None
        surface_path : str | None, optional
            The path to a custom surface, by default None

        Raises
        ------
        ValueError
            If surface_type == "central" and subpath is not set
        ValueError
            If surface_type == "custom" and surface_path is not set
        ValueError
            If surface_type not in {"central", "custom"}
        """
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
        #self._mask = np.zeros((self._mesh.nodes.nr), dtype=np.bool_)
        self._mask = np.ones((self._mesh.nodes.nr), dtype=np.bool_)

    def load_mesh(self, mesh: str | Msh | None = None, subpath: str | None = None):
        """Initializes the Region of Interest with a mesh either with a .msh file or a m2m folder.
        The mesh will then be subject to masking to define a Region of Interest in the mesh.

        Parameters
        ----------
        mesh : str | Msh | None, optional
            The mesh to initilize the ROI, can be a path to an .msh file or a Msh instance, by default None
        subpath : str | None, optional
            m2m folder to load the .msh file from, by default None

        Raises
        ------
        ValueError
            If mesh and subpath are None
        """
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
        #self._mask = np.zeros((self._mesh.elm.nr), dtype=np.bool_)
        self._mask = np.ones((self._mesh.elm.nr), dtype=np.bool_)

    def apply_surface_mask(
        self,
        surface_type: str,
        mask_space: str | list[str] | None = None,
        mask_path: str | list[str] | None = None,
        mask_value: int | list[int] | None = None,
        mask_operator: str | list[str] | None = None,
        subpath: str | None = None,
    ):
        """Applies a surface mask based.

        Parameters
        ----------
        surface_type : str {"central", "custom"}
            The type of the surface
        mask_space : str | list[str] | None, optional {"subject", "subject_lh", "fs_avg_lh", "subject_rh", "fs_avg_rh", "mni"}
            The space the mask is defined in, by default None
        mask_path : str | list[str] | None, optional
            The path to the mask file, by default None
        mask_value : int | list[int] | None, optional {.label, .annot, curv, nifti}
            The value that in the mask file that describes the mask, by default None
        mask_operator : str | list[str] | None, optional {"union", "intersection", "difference"}
            The operator to be used, by default "union"
            union: current_mask or surface_mask
            intersection: current_mask and surface_mask
            difference: current_mask and not surface_mask
        subpath : str | None, optional
            The path to the m2m folder of the subject, used for transformation from MNI / fs average to subject space, by default None

        Raises
        ------
        ValueError
            mask_path is not set
        ValueError
            mask_space is not in {"subject", "subject_lh", "fs_avg_lh", "subject_rh", "fs_avg_rh", "mni"}
        ValueError
            mask_path, mask_space, mask_value and mask_operator have different counts
        ValueError
            If surface_type == "custom" but mask_space not in {"subject", "mni"}
        """
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
                #mask_operator = "union"
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
                if os.path.splitext(mask_path)[1] in [".nii", ".gz"]:
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
        node_type: str | None = None,
        mask_space: str | list[str] | None = None,
        mask_path: str | list[str] | None = None,
        mask_value: int | list[int] | None = None,
        mask_operator: str | list[str] | None = None,
        subpath: str | None = None,
    ):
        """Applies a volume mask based on Nifti images.

        Parameters
        ----------
        node_type : str | None, optional {"node", "elm_center"}
            Weather to use the nodes or the element center, by default elm_center
        mask_space : str | list[str] | None, optional {"subject", "mni"}
            The space the Nifti mask images are defined in, by default None
        mask_path : str | list[str] | None, optional
            The path to the Nifti mask images, by default None
        mask_value : int | list[int] | None, optional
            The value in the Nifti image that describes the mask, by default 1
        mask_operator : str | list[str] | None, optional {"union", "intersection", "difference"}
            The operator to be used, by default "union"
            union: current_mask or volume_mask
            intersection: current_mask and volume_mask
            difference: current_mask and not volume_mask
        subpath : str | None, optional
            The path to the m2m folder of the subject, used for transformation from MNI to subject space, by default None

        Raises
        ------
        ValueError
            If mask_path is not set
        ValueError
            If mask_space is not in {"subject", "mni"}
        ValueError
            If mask_path, mask_space, mask_value and mask_operator have different counts
        """
        if node_type is None:
            node_type = "elm_center"
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
                #mask_operator = "union"
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
        """Applies a mask based on a Region of Interest that is defined by a surface.
        The volume will be masked based on the distance to the surface ROI.

        Parameters
        ----------
        surface_roi : RegionOfInterest
            The surface ROI
        surface_inclusion_radius : float | None
            The radius from each node of the surface roi in which the volume is included.
        surface_roi_operator : str | None, optional {"union", "intersection", "difference"}
            The operator to be used, by default "union"
            union: current_mask or surface_roi_mask
            intersection: current_mask and surface_roi_mask
            difference: current_mask and not surface_roi_mask

        Raises
        ------
        ValueError
            If surface_inclusion_radius is not set
        """
        if surface_inclusion_radius is None:
            raise ValueError(
                f"surface_inclusion_radius needs to be set (was {surface_inclusion_radius})"
            )

        if surface_roi_operator is None:
            #surface_roi_operator = "union"
            surface_roi_operator = "intersection"

        kd_tree = scipy.spatial.cKDTree(self._mesh.elements_baricenters().value)
        index_mask = kd_tree.query_ball_point(
            surface_roi.get_nodes(), surface_inclusion_radius
        )
        index_mask = np.unique(np.concatenate(index_mask)).astype(np.int_)
        self._mask = combine_mask(self._mask, index_mask, surface_roi_operator)

    def apply_sphere_mask(
        self,
        node_type: str | None = None,
        roi_sphere_center: list[float] | list[list[float]] | None = None,
        roi_sphere_radius: float | list[float] | None = None,
        roi_sphere_center_space: str | list[str] | None = None,
        roi_sphere_operator: str | list[str] | None = None,
        subpath: str | None = None,
    ):
        """Applies a mask based on one or multiple spheres. Masked are applied in the input order.

        Parameters
        ----------
        node_type : str | None, optional {"node", "elm_center"}
            Weather to use the nodes or the element center, by default elm_center
        roi_sphere_center : list[float] | list[list[float]] | None, optional
            Center coordinates of the spheres that are used to create masks, by default None
        roi_sphere_radius : float | list[float] | None, optional
            The radii of the spheres used to create masks, by default None
        roi_sphere_center_space : str | list[str] | None, optional {"subject", "mni"}
            The space the sphere centers are defined in, by default None
        roi_sphere_operator : str | list[str] | None, optional {"union", "intersection", "difference"}
            The operator to be used, by default "union"
            union: current_mask or sphere_mask
            intersection: current_mask and sphere_mask
            difference: current_mask and not sphere_mask
        subpath : str | None, optional
            The path to the m2m folder of the subject, used for transformation from MNI to subject space, by default None

        Raises
        ------
        ValueError
            If roi_sphere_center is not set
        ValueError
            If roi_sphere_radius is not set
        ValueError
            If roi_sphere_center_space is not set
        ValueError
            If roi_sphere_center, roi_sphere_radius, roi_sphere_center_space and roi_sphere_operator have different counts
        ValueError
            If node_type is not in {"node", "elm_center"}
        ValueError
            If roi_sphere_center_space is not in {"subject", "mni"}
        """
        if node_type is None:
            node_type = "elm_center"
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
                #roi_sphere_operator = "union"
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
        tissues: int | list[int] | ElementTags | list[ElementTags] | None = None,
        tissue_mask_operator: str | None = None,
    ):
        """Applies a mask based on the tissue tags in the mesh using the tissue_mask_operator.

        Parameters
        ----------
        tissues : int | list[int] | ElementTags | list[ElementTags] | None
            The tissues that are supposed to be masked, by default GM
        tissue_mask_operator : str | None, optional {"union", "intersection", "difference"}
            The operator to be used, by default "union"
            union: current_mask or tissue_mask
            intersection: current_mask and tissue_mask
            difference: current_mask and not tissue_mask
        """
        if tissues is None:
            tissues = ElementTags.GM

        if not isinstance(tissues, list):
            tissues = [tissues]

        if tissue_mask_operator is None:
            #tissue_mask_operator = "union"
            tissue_mask_operator = "intersection"

        tissues = np.array(tissues)

        self._mask = combine_mask(
            self._mask,
            self._mesh.elm.elm_number[np.in1d(self._mesh.elm.tag1, tissues)] - 1,
            tissue_mask_operator,
        )

    def invert(self):
        """Inverts the current mask"""
        self._mask = ~self._mask
        
    def run(self, cpus=1):
        """writes visualization for visual control"""
        if self.fname_visu is None:
            raise ValueError(
                    'fname_visu (mesh filename for ROI visualization) has to be set'
            )
        folder_path, base_file_name = os.path.split(os.path.abspath(os.path.expanduser(self.fname_visu)))
        base_file_name = os.path.splitext(base_file_name)[0]
        self.write_visualization(folder_path, base_file_name)


def _init_subject_files(subpath: str | None) -> SubjectFiles:
    """Initializes subject files from a path to a m2m folder

    Parameters
    ----------
    subpath : str | None
        The path to a m2m folder

    Returns
    -------
    SubjectFiles
        The initialized subject files

    Raises
    ------
    ValueError
        If subpath is None
    IOError
        If subpath is not a folder
    """
    if subpath is None:
        raise ValueError(f"subpath needs to be set (was {subpath})")
    subpath = os.path.abspath(os.path.expanduser(subpath))
    if not os.path.isdir(subpath):
        raise IOError(f"Cannot locate subjects m2m folder: {subpath}")

    return SubjectFiles(subpath=subpath)


def load_surface_from_file(surface_path: str) -> Msh:
    """Loads a surface from a .msh (Gmsh), .gii (GIfTI surface) or FreeSurfer surface file

    Parameters
    ----------
    surface_path : str
        The path to the surface file

    Returns
    -------
    Msh
        The loaded surface

    Raises
    ------
    ValueError
        If surface_path is not a file
    ValueError
        If a .msh file contains elements other then triangles
    ValueError
        If the file type is not one of .msh (Gmsh), .gii (GIfTI surface) or FreeSurfer surface
    """
    if not os.path.isfile(surface_path):
        raise ValueError("surface_path needs to be a file")

    _, file_extension = os.path.splitext(surface_path)
    surface = None
    match file_extension:
        case ".msh":
            surface: Msh = mesh_io.read_msh(surface_path)
            if np.any(surface.elm.elm_type != 2):
                raise ValueError(".msh file contains non triangle elements")
        case ".gii":
            surface = mesh_io.read_gifti_surface(surface_path)
        case _:
            try:
                surface = mesh_io.read_freesurfer_surface(surface_path, True)
            except:
                pass
    if surface is None:
        raise ValueError(
            f"surface_path needs to be in one of the following file formats: [.msh (Gmsh), .gii (GIfTI surface), FreeSurfer surface] (was {file_extension})"
        )
    return surface


def load_surface_mask_from_file(
    mask_path: str, mask_value: int = 1
) -> npt.NDArray[np.int_]:
    """Loads a surface mask from a .label, .annot or freesurfer curv file

    Parameters
    ----------
    mask_path : str
        The path to the mask file
    mask_value : int, optional
        The value in the mask file that describes the mask, by default 1
        Ignored if .label file

    Returns
    -------
    npt.NDArray[np.int_]
        The index mask loaded from the file

    Raises
    ------
    ValueError
        If surface_path is not a file
    ValueError
        If the file type is not one of .label, .annot or freesurfer curv
    """
    if not os.path.isfile(mask_path):
        raise ValueError("surface_path needs to be a file")

    _, file_extension = os.path.splitext(mask_path)
    index_mask = None
    match file_extension:
        case ".label":
            index_mask = nib.freesurfer.io.read_label(mask_path)
        case ".annot":
            labels, _, _ = nib.freesurfer.io.read_annot(mask_path)
            index_mask = np.where(labels == mask_value)[0]
        case _:
            try:
                data = nib.freesurfer.io.read_morph_data(mask_path)
                index_mask = np.where(np.rint(data) == mask_value)[0]
            except:
                pass

    if index_mask is None:
        raise ValueError(
            f"mask_path needs to be in one of the following file formats: [.label (FreeSurfer label), .annot (FreeSurfer annot), FreeSurfer curv] (was {file_extension})"
        )

    return index_mask.astype(np.int_)


def fs_avr_mask_to_sub(
    index_mask: npt.NDArray[np.int_], hemi: str, subject_files: SubjectFiles
) -> npt.NDArray[np.int_]:
    """Turns a mask, defined on the fs avr space into subject space

    Parameters
    ----------
    index_mask : npt.NDArray[np.int_]
        The index mask in fs avr space
    hemi : str {"lh", "rh"}
        The hemisphere used
    subject_files : SubjectFiles
        The subject files of the subject which describes the subject space to transform into

    Returns
    -------
    npt.NDArray[np.int_]
        The index mask in subject space
    """
    sphere_surface: Msh = mesh_io.read_gifti_surface(
        file_finder.get_reference_surf(hemi, "sphere")
    )
    registration_surface: Msh = mesh_io.read_gifti_surface(
        subject_files.get_surface(hemi, "sphere_reg")
    )

    morph = transformations.SurfaceMorph(
        sphere_surface, registration_surface, method="nearest"
    )

    bool_mask = np.zeros((sphere_surface.nodes.nr), dtype=np.int_)
    bool_mask[index_mask] = 1

    index_mask = np.where(morph.transform(bool_mask) > 0.0001)[0]

    return index_mask


def combine_mask(
    bool_mask: npt.NDArray[np.bool_], index_mask: npt.NDArray[np.int_], operator: str
) -> npt.NDArray[np.bool_]:
    """Combines a boolean mask and an index mask based on an operator
    bool_mask *(operator)* index_mask

    Parameters
    ----------
    bool_mask : npt.NDArray[np.bool_]
        The boolean mask that the index mask is applied on top of
    index_mask : npt.NDArray[np.int_]
        The index mask that is combined with the boolean mask
    operator : str {"union", "intersection", "difference"}
        The operator to be used.
        union: bool_mask or index_mask
        intersection: bool_mask and index_mask
        difference: bool_mask and not index_mask

    Returns
    -------
    npt.NDArray[np.bool_]
        The resulting boolean mask

    Raises
    ------
    ValueError
        If mask_space is not in {"union", "intersection", "difference"}
    """
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
            raise ValueError(
                'mask_space needs to be one of ["union", "intersection", "difference"]'
            )

    return bool_mask


def mni_mask_to_sub(
    mask_img: nib.Nifti1Image, subject_files: SubjectFiles
) -> nib.Nifti1Image:
    """Turns a Nifti image that is in MNI space into subject space

    Parameters
    ----------
    mask_img : nib.Nifti1Image
        The input Nifti image
    subject_files : SubjectFiles
        The subject files of the subject which describes the subject space to transform into

    Returns
    -------
    nib.Nifti1Image
        The Nifti image in subject space

    Raises
    ------
    ValueError
        If the subject folder does not contain a T1 image as a reference
    """
    if os.path.isfile(subject_files.T1_upsampled):
        target_image = nib.load(subject_files.T1_upsampled)
    elif os.path.isfile(subject_files.reference_volume):
        target_image = nib.load(subject_files.reference_volume)
    else:
        raise ValueError(
            "Subject folder does not contain a T1 image for the MNI transformation"
        )

    target_dim = list(target_image.get_fdata().shape)

    image_deformation = nib.load(subject_files.conf2mni_nonl)

    transformed_mask = transformations.volumetric_nonlinear(
        (mask_img.get_fdata(), mask_img.affine),
        (image_deformation.get_fdata(), image_deformation.affine),
        target_space_affine=target_image.affine,
        target_dimensions=target_dim,
        intorder=0,
    )
    transformed_mask = np.squeeze(transformed_mask)

    return nib.Nifti1Image(transformed_mask, target_image.affine)


def mask_image_to_index_mask(
    node_type: str, mask_img: nib.Nifti1Image, mesh: Msh, mask_value: int
) -> npt.NDArray[np.int_]:
    """Turns a Nifti image mask (mask_img) into a list of mesh element or node indexes.

    Parameters
    ----------
    node_type : str {"node", "elm_center"}
        Weather to return the node indexes or the element indexes
    mask_img : nib.Nifti1Image
        The Nifti image mask
    mesh : Msh
        The mesh that is supposed to be masked
    mask_value : int
        The value in the mask_img that describes the mask

    Returns
    -------
    npt.NDArray[np.int_]
        A list of mesh or node indexes that are within the image mask

    Raises
    ------
    ValueError
        If node_type not in {"node", "elm_center"}
    """
    match node_type:
        case "node":
            data = mesh_io.NodeData.from_data_grid(
                mesh, mask_img.get_fdata(), mask_img.affine, "", order=0
            ).value
        case "elm_center":
            data = mesh_io.ElementData.from_data_grid(
                mesh, mask_img.get_fdata(), mask_img.affine, "", order=0
            ).value
        case _:
            raise ValueError(
                f'node_type needs to be one of ["node", "elm_center"] (was {node_type})'
            )

    return np.where(np.rint(data) == mask_value)[0]
