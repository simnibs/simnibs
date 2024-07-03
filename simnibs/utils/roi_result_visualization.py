from copy import deepcopy
import os
import re

import numpy as np
from simnibs.mesh_tools import gmsh_view
from simnibs.mesh_tools.mesh_io import Msh
from simnibs.utils.mesh_element_properties import ElementTags
from simnibs.utils.region_of_interest import RegionOfInterest


class RoiResultVisualization:

    def __init__(
        self,
        rois: list[RegionOfInterest],
        sim_result_filenames: list[str],
        folder_path: str,
        base_file_name: str,
        roi_names: list[str] | None = None,
        result_prefixes: list[str] | None = None,
        base_head_mesh: Msh | None = None,
    ) -> None:
        self.head_mesh: Msh | None = None
        self.surface_mesh: Msh | None = None

        self.sim_result_filenames = sim_result_filenames
        self.folder_path = folder_path
        self.base_file_name = base_file_name
        self.geo_file_name = os.path.join(
            self.folder_path, f"{self.base_file_name}.geo"
        )

        if roi_names is None:
            self.roi_names = [f"roi_{i}" for i in np.arange(len(rois))]
        else:
            self.roi_names = roi_names

        if result_prefixes is None:
            self.result_prefixes = [f"{i}__" for i in np.arange(len(rois))]
        else:
            self.result_prefixes = result_prefixes
            prepared_result_prefixes = []
            for prefix in self.result_prefixes:
                if prefix == "":
                    prepared_result_prefixes.append(prefix)
                else:
                    prepared_result_prefixes.append(f"{prefix}__")
            self.result_prefixes = prepared_result_prefixes

        if base_head_mesh is not None:
            self.head_mesh = base_head_mesh

        self._head_mesh_roi_names = []
        self._surface_mesh_roi_names = []
        self._head_mesh_roi_counts = {'ed': 0, 'nd': 0}
        self._surface_mesh_roi_counts = {'ed': 0, 'nd': 0}
        surface_names_to_surface: dict[str, Msh] = {}
        for i, roi in enumerate(rois):
            if 4 in roi._mesh.elm.elm_type:
                self._head_mesh_roi_names.append(self.roi_names[i])
                if self.head_mesh is None:
                    self.head_mesh = roi._mesh
                match roi._mask_type:
                    case "node":
                        self.head_mesh.add_node_field(roi._mask, self.roi_names[i])
                        self._head_mesh_roi_counts['nd'] += 1
                    case "elm_center":
                        self.head_mesh.add_element_field(roi._mask, self.roi_names[i])
                        self.head_mesh.field[self.roi_names[i]].assign_triangle_values()
                        self._head_mesh_roi_counts['ed'] += 1
                    case _:
                        raise NotImplementedError()
            else:
                self._surface_mesh_roi_names.append(self.roi_names[i])
                match roi.surface_type:
                    case "central":
                        surface_name = "central"
                        if surface_name not in surface_names_to_surface:
                            surface_names_to_surface[surface_name] = roi._mesh
                            surface_names_to_surface[surface_name].elm.tag1[
                                roi._surface_divide_elements :
                            ] = ElementTags.LH_CENTRAL_GM
                            surface_names_to_surface[surface_name].elm.tag2[
                                roi._surface_divide_elements :
                            ] = ElementTags.LH_CENTRAL_GM
                            surface_names_to_surface[surface_name].elm.tag1[
                                : roi._surface_divide_elements
                            ] = ElementTags.RH_CENTRAL_GM
                            surface_names_to_surface[surface_name].elm.tag2[
                                : roi._surface_divide_elements
                            ] = ElementTags.RH_CENTRAL_GM
                    case "custom":
                        surface_name = f"surface_{i}"
                        if roi.surface_path is not None:
                            surface_name = roi.surface_path
                        if surface_name not in surface_names_to_surface:
                            surface_names_to_surface[surface_name] = roi._mesh
                            surface_names_to_surface[surface_name].elm.tag1[:] = (
                                ElementTags.RH_CENTRAL_SURFACE_END + i
                            )
                            surface_names_to_surface[surface_name].elm.tag2[:] = (
                                ElementTags.RH_CENTRAL_SURFACE_END + i
                            )
                    case _:
                        raise NotImplementedError()

                match roi._mask_type:
                    case "node":
                        surface_names_to_surface[surface_name].add_node_field(
                            roi._mask, self.roi_names[i]
                        )
                        self._surface_mesh_roi_counts['nd'] += 1
                    case "elm_center":
                        surface_names_to_surface[surface_name].add_element_field(
                            roi._mask, self.roi_names[i]
                        )
                        self._surface_mesh_roi_counts['ed'] += 1
                    case _:
                        raise NotImplementedError()

        if len(surface_names_to_surface) > 0:
            self.surface_mesh = Msh()
            for surface_name in surface_names_to_surface:
                self.surface_mesh = self.surface_mesh.join_mesh(
                    surface_names_to_surface[surface_name]
                )

        self.elm_field_count_per_file = []
        self.node_field_count_per_file = []
        for i, filename in enumerate(self.sim_result_filenames):
            result_mesh = Msh(fn=filename)
            if self.head_mesh is not None:
                result_mesh = result_mesh.crop_mesh(
                    tags=np.unique(self.head_mesh.elm.tag1)
                )
            self.elm_field_count_per_file.append(len(result_mesh.elmdata))
            self.node_field_count_per_file.append(len(result_mesh.nodedata))

            for elm_data in result_mesh.elmdata:
                if self.head_mesh is not None:
                    self.head_mesh.add_element_field(
                        elm_data.value[: self.head_mesh.elm.nr],
                        f"{self.result_prefixes[i]}{elm_data.field_name}",
                    )

                if self.surface_mesh is not None:
                    self.surface_mesh.add_element_field(
                        elm_data.interpolate_scattered(
                            self.surface_mesh.elements_baricenters().value
                        ),
                        f"{self.result_prefixes[i]}{elm_data.field_name}",
                    )
                
            for node_data in result_mesh.nodedata:
                if self.head_mesh is not None:
                    self.head_mesh.add_node_field(
                        node_data.value[: self.head_mesh.nodes.nr],
                        f"{self.result_prefixes[i]}{node_data.field_name}",
                    )

                if self.surface_mesh is not None:
                    self.surface_mesh.add_node_field(
                        node_data.interpolate_scattered(
                            self.surface_mesh.nodes.node_coord
                        ),
                        f"{self.result_prefixes[i]}{node_data.field_name}",
                    )

    def has_surface_mesh(self) -> bool:
        return self.surface_mesh is not None

    def get_surface_mesh(self) -> Msh:
        if self.has_surface_mesh():
            return self.surface_mesh
        else:
            raise AttributeError()

    def has_head_mesh(self) -> bool:
        return self.head_mesh is not None

    def get_head_mesh(self) -> Msh:
        if self.has_head_mesh():
            return self.head_mesh
        else:
            raise AttributeError()
        
    def remove_field_from_surface_mesh(self, field_name: str):
        if self.surface_mesh is not None:
            RoiResultVisualization._remove_field(self.surface_mesh, self.surface_mesh_opt, field_name)

    def remove_field_from_head_mesh(self, field_name: str):
        if self.head_mesh is not None:
            RoiResultVisualization._remove_field(self.head_mesh, self.head_mesh_opt, field_name)
    
    @staticmethod
    def _remove_field(mesh:Msh, gmsh_options: gmsh_view.Visualization, field_name):
        view_index = 0
        found = False
        for node_data_field in mesh.nodedata:
            if field_name == node_data_field.field_name:
                found = True
                mesh.nodedata.remove(node_data_field)
                break
            view_index += 1

        if not found:
            for elm_data_field in mesh.elmdata:
                if field_name == elm_data_field.field_name:
                    found = True
                    mesh.elmdata.remove(elm_data_field)
                    break
                view_index += 1

        if found:
            del gmsh_options.View[view_index]

            for i, current_gmsh_view in enumerate(gmsh_options.View):
                current_gmsh_view.indx = i

    def create_visualization(self):
        gmsh_options: list[gmsh_view.Visualization] = []
        for sim_result_filename in self.sim_result_filenames:
            gmsh_options.append(
                RoiResultVisualization._read_opt(f"{sim_result_filename}.opt")
            )

        self.head_mesh_data_name_to_gmsh_view = {}
        self.surface_mesh_data_name_to_gmsh_view = {}
        meshes_and_opt = []
        if self.head_mesh is not None:
            self.head_mesh_opt = deepcopy(gmsh_options[0])
            self.head_mesh_opt.View = []
            self.head_mesh_opt.merge = []
            self.head_mesh_opt.add_merge(self.geo_file_name)
            unique_tags = np.unique(self.head_mesh.elm.tag1)
            for tag in list(self.head_mesh_opt.PhysicalNames.PhysicalSurfaces.keys()):
                if tag not in unique_tags:
                    del self.head_mesh_opt.PhysicalNames.PhysicalSurfaces[tag]

            for tag in list(self.head_mesh_opt.PhysicalNames.PhysicalVolumes.keys()):
                if tag not in unique_tags:
                    del self.head_mesh_opt.PhysicalNames.PhysicalVolumes[tag]

            meshes_and_opt.append(
                (
                    self.head_mesh,
                    self.head_mesh_opt,
                    self.head_mesh_data_name_to_gmsh_view,
                    self._head_mesh_roi_counts
                )
            )

        if self.surface_mesh is not None:
            self.surface_mesh_opt = deepcopy(gmsh_options[0])
            self.surface_mesh_opt.View = []
            self.surface_mesh_opt.merge = []
            self.surface_mesh_opt.add_merge(self.geo_file_name)
            self.surface_mesh_opt.PhysicalNames = gmsh_view.PhysicalNames(None, None)
            self.surface_mesh_opt.visibility = np.unique(self.surface_mesh.elm.tag1)
            meshes_and_opt.append(
                (
                    self.surface_mesh,
                    self.surface_mesh_opt,
                    self.surface_mesh_data_name_to_gmsh_view,
                    self._surface_mesh_roi_counts
                )
            )

        for mesh, mesh_opt, mesh_data_name_to_gmsh_view, roi_counts in meshes_and_opt:
            
            # 1. node data of rois
            data_count_sum = 0
            current_view_index = 0
            for i in range(roi_counts['nd']):
                mesh_opt.add_view()
                current_view_index += 1
                mesh_data_name_to_gmsh_view[mesh.nodedata[i].field_name] = (
                    mesh_opt.View[-1]
                )
                data_count_sum += 1

            # 2. node data of all input files
            offset = 0
            for i, result_filename in enumerate(self.sim_result_filenames):
                node_data_count = self.node_field_count_per_file[i]
                opt = gmsh_options[i]
                for j in range(
                    roi_counts['nd'],
                    node_data_count + roi_counts['nd'],
                ):
                    view: gmsh_view.View = deepcopy(
                        opt.View[j - roi_counts['nd']]
                    )
                    view.indx = current_view_index
                    current_view_index += 1
                    mesh_opt.View.append(view)
                    mesh_data_name_to_gmsh_view[mesh.nodedata[j+offset].field_name] = view     
                offset += node_data_count
                data_count_sum += node_data_count

            # 3. node data added inbetween
            for i in range(len(mesh.nodedata) - data_count_sum):
                mesh_opt.add_view()
                current_view_index += 1
                mesh_data_name_to_gmsh_view[
                    mesh.nodedata[data_count_sum + i].field_name
                ] = mesh_opt.View[-1]

            # 4. elm data of rois
            data_count_sum = 0
            for i in range(roi_counts['ed']):
                mesh_opt.add_view()
                current_view_index += 1
                mesh_data_name_to_gmsh_view[mesh.elmdata[i].field_name] = mesh_opt.View[
                    -1
                ]
                data_count_sum += 1

            # 5. elm data of all input files
            offset = 0
            for i, result_filename in enumerate(self.sim_result_filenames):
                node_data_count = self.node_field_count_per_file[i]
                elm_data_count = self.elm_field_count_per_file[i]
                opt = gmsh_options[i]
                for j in range(
                    node_data_count + roi_counts['ed'],
                    node_data_count + roi_counts['ed'] + elm_data_count,
                ):
                    view: gmsh_view.View = deepcopy(
                        opt.View[j - roi_counts['ed']]
                    )
                    view.indx = current_view_index
                    current_view_index += 1
                    mesh_opt.View.append(view)
                    mesh_data_name_to_gmsh_view[mesh.elmdata[j+offset].field_name] = view        
                offset += elm_data_count
                data_count_sum += node_data_count + elm_data_count

            # 6. elm data added inbetween
            for i in range(len(mesh.elmdata) - data_count_sum):
                mesh_opt.add_view()
                current_view_index += 1
                mesh_data_name_to_gmsh_view[
                    mesh.elmdata[data_count_sum + i].field_name
                ] = mesh_opt.View[-1]

        # add roi view settings
        if self.head_mesh is not None:
            for roi_name in self._head_mesh_roi_names:
                roi_view = self.head_mesh_data_name_to_gmsh_view[roi_name]
                roi_view.Visible = 0
                roi_view.ShowScale = 0
                roi_view.CustomMin = 0.9
                roi_view.CustomMax = 1
                roi_view.RangeType = 2

        if self.surface_mesh is not None:
            for roi_name in self._surface_mesh_roi_names:
                roi_view = self.surface_mesh_data_name_to_gmsh_view[roi_name]
                roi_view.Visible = 0
                roi_view.ShowScale = 0
                roi_view.CustomMin = 0
                roi_view.CustomMax = 1
                roi_view.RangeType = 2

        geo_filenames = []
        for sim_result_filename in self.sim_result_filenames:
            match os.path.basename(sim_result_filename).split("_")[1]:
                case "TMS":
                    geo_filenames.append(
                        f"{'_'.join(sim_result_filename.split('_')[:-1])}_coil_pos.geo"
                    )
                case "TDCS":
                    geo_filenames.append(
                        f"{'_'.join(sim_result_filename.split('_')[:-1])}_el_currents.geo"
                    )
                case _:
                    raise NotImplementedError

        dont_duplicate_views = ["scalp"]
        added_views = []
        geo_content = []
        added_geo_views_indexes = []
        post_geo_content = []
        post_added_geo_views_indexes = []
        for i, geo_filename in enumerate(geo_filenames):
            geo_views = RoiResultVisualization._read_geo(geo_filename)
            for j, geo_view_name in enumerate(geo_views):
                if geo_view_name in dont_duplicate_views:
                    if geo_view_name not in added_views:
                        added_views.append(geo_view_name)
                        post_geo_content.append(geo_views[geo_view_name])
                        post_added_geo_views_indexes.append([i,j])
                else:
                    new_geo_view_name = f"{self.result_prefixes[i]}{geo_view_name}"
                    added_views.append(new_geo_view_name)
                    geo_content.append(geo_views[geo_view_name].replace(
                        geo_view_name, new_geo_view_name
                    ))
                    added_geo_views_indexes.append([i,j])
        
        self.geo_content = ''.join(geo_content) + ''.join(post_geo_content)
        added_geo_views_indexes.extend(post_added_geo_views_indexes)
        
        for idx_fn, idx_v in added_geo_views_indexes:
            view = gmsh_options[idx_fn].View[
                self.node_field_count_per_file[idx_fn]
                + self.elm_field_count_per_file[idx_fn]
                + idx_v
            ]
            if self.head_mesh is not None:
                view2: gmsh_view.View = deepcopy(view)
                view2.indx = self.head_mesh_opt.View[-1].indx + 1
                self.head_mesh_opt.View.append(view2)
                
            if self.surface_mesh is not None:
                view2: gmsh_view.View = deepcopy(view)
                view2.indx = self.surface_mesh_opt.View[-1].indx + 1
                self.surface_mesh_opt.View.append(view2)

    def write_visualization(self):
        """Writes a visualization of the results and region of interests to a folder

        Parameters
        ----------
        folder_path : str
            Folder to write the visualization to
        base_file_name : str
            The base file name of the visualization files
        """
        file_names = []
        if self.head_mesh is not None:
            file_names.append(
                os.path.join(self.folder_path, f"{self.base_file_name}_head_mesh.msh")
            )
            self.head_mesh.write(file_names[-1])

        if self.surface_mesh is not None:
            file_names.append(
                os.path.join(
                    self.folder_path, f"{self.base_file_name}_surface_mesh.msh"
                )
            )
            self.surface_mesh.write(file_names[-1])

        self.write_gmsh_options()

        with open(
            os.path.join(self.folder_path, f"{self.base_file_name}.geo"), "w"
        ) as file:
            file.write(self.geo_content)

        return file_names

    def write_gmsh_options(self):
        if self.head_mesh is not None:
            self.head_mesh_opt.write_opt(
                os.path.join(self.folder_path, f"{self.base_file_name}_head_mesh.msh")
            )

        if self.surface_mesh is not None:
            self.surface_mesh_opt.write_opt(
                os.path.join(
                    self.folder_path, f"{self.base_file_name}_surface_mesh.msh"
                )
            )

    @staticmethod
    def _read_geo(geo_filename: str) -> dict[str, str]:
        with open(geo_filename, "r") as file:
            data = file.readlines()
        views = {}
        current_view_name: str = None
        for line in data:
            if current_view_name is None:
                current_view_name = re.findall(r'View"([^"]+)"\{', line)[0]
                views[current_view_name] = []

            views[current_view_name].append(line)

            if line.startswith("};"):
                views[current_view_name] = ''.join(views[current_view_name])
                current_view_name = None

        return views

    @staticmethod
    def _read_opt(opt_filename: str) -> gmsh_view.Visualization:
        with open(opt_filename, "r") as file:
            data = file.readlines()

        opt = gmsh_view.Visualization(None)
        opt.View = []
        opt.visibility = [ElementTags.GM_TH_SURFACE.value]
        for line in data:
            if line.startswith("Physical"):
                pattern = (
                    r'Physical (Volume|Surface) \("([^"]+)",(\d+)\) = \{ (\d+) \};'
                )
                match = re.search(pattern, line)
                if match is not None:
                    if "Volume" in line:
                        opt.PhysicalNames.PhysicalVolumes[int(match.group(3))] = (
                            match.group(2)
                        )
                    elif "Surface" in line:
                        opt.PhysicalNames.PhysicalSurfaces[int(match.group(3))] = (
                            match.group(2)
                        )
                else:
                    raise IOError()

            elif line.startswith("General"):
                pattern = r"General\.([A-Za-z0-9_\.]+) = ([^;]+);"
                match = re.search(pattern, line)
                if match is not None:
                    setattr(opt.General, match.group(1), match.group(2))
                else:
                    raise IOError()

            elif line.startswith("Mesh"):
                pattern = r"Mesh\.([A-Za-z0-9_\.]+) = ([^;]+);"
                match = re.search(pattern, line)
                if match is not None:
                    if line.startswith('Mesh.Color'):
                        data_str = match.group(2).strip("{}")
                        data_list = list(map(int, data_str.split(",")))
                        setattr(opt.Mesh.Color, match.group(1).removeprefix('Color.'), data_list)
                    else:
                        setattr(opt.Mesh, match.group(1), match.group(2))
                else:
                    raise IOError()
            elif line.startswith("View"):
                pattern = r"View\[(\d+)\]\.([A-Za-z0-9_\.]+) = ([^;]+);"
                match = re.search(pattern, line)
                if match is not None:
                    view_index = int(match.group(1))
                    name = match.group(2)
                    value = match.group(3)

                    if len(opt.View) <= view_index:
                        opt.add_view()

                    if name == "ColorTable":
                        data_str = value.removeprefix("{")
                        data_str = data_str.removesuffix("}")
                        matches = re.findall(r"\{([^}]+)\}", data_str)
                        data_list = [
                            list(map(int, match.split(","))) for match in matches
                        ]
                        value = np.array(data_list)
                    setattr(opt.View[view_index], name, value)
                else:
                    raise IOError()

        return opt
