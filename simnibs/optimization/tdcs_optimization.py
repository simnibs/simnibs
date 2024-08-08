import copy
import csv
import re
import os
import h5py
import functools
import logging

import nibabel
import numpy as np
import numpy.typing as npt
import scipy.linalg
import scipy.optimize
import scipy.spatial

from ..mesh_tools import mesh_io, gmsh_view
from ..utils import transformations
from ..utils.simnibs_logger import logger
from ..utils.matlab_read import try_to_read_matlab_field, remove_None
from ..utils.mesh_element_properties import ElementTags


class TDCSoptimize:
    """Defines a tdcs optimization problem

    Parameters
    --------------
    leadfield_hdf: str (optional)
        Name of file with leadfield
    max_total_current: float (optional)
        Maximum current across all electrodes (in Amperes). Default: 2e-3
    max_individual_current: float (optional)
        Maximum current for any single electrode (in Amperes). Default: 1e-3
    max_active_electrodes: int (optional)
        Maximum number of active electrodes. Default: no maximum
    name: str (optional)
        Name of optimization problem. Default: optimization
    target: list of TDCStarget objects (optional)
        Targets for the optimization. Default: no target
    avoid: list of TDCSavoid objects
        list of TDCSavoid objects defining regions to avoid


    Attributes
    --------------
    leadfield_hdf: str
        Name of file with leadfield
    max_total_current: float (optional)
        Maximum current across all electrodes (in Amperes). Default: 2e-3
    max_individual_current: float
        Maximum current for any single electrode (in Amperes). Default: 1e-3
    max_active_electrodes: int
        Maximum number of active electrodes. Default: no maximum

    ledfield_path: str
        Path to the leadfield in the hdf5 file. Default: '/mesh_leadfield/leadfields/tdcs_leadfield'
    mesh_path: str
        Path to the mesh in the hdf5 file. Default: '/mesh_leadfield/'

    The two above are used to define:

    mesh: simnibs.msh.mesh_io.Msh
        Mesh with problem geometry

    leadfield: np.ndarray
        Leadfield matrix (N_elec -1 x M x 3) where M is either the number of nodes or the
        number of elements in the mesh. We assume that there is a reference electrode

    Alternatively, you can set the three attributes above and not leadfield_path,
    mesh_path and leadfield_hdf

    lf_type: None, 'node' or 'element'
        Type of leadfield.

    name: str
        Name for the optimization problem. Defaults tp 'optimization'

    target: list of TDCStarget objects
        list of TDCStarget objects defining the targets of the optimization

    avoid: list of TDCSavoid objects (optional)
        list of TDCSavoid objects defining regions to avoid

    open_in_gmsh: bool (optional)
        Whether to open the result in Gmsh after the calculations. Default: False

    Warning
    -----------
    Changing leadfield_hdf, leadfield_path and mesh_path after constructing the class
    can cause unexpected behaviour
    """

    def __init__(
        self,
        leadfield_hdf=None,
        max_total_current=2e-3,
        max_individual_current=1e-3,
        max_active_electrodes=None,
        name="optimization/tdcs",
        target=None,
        avoid=None,
        open_in_gmsh=True,
    ):
        self.leadfield_hdf = leadfield_hdf
        self.max_total_current = max_total_current
        self.max_individual_current = max_individual_current
        self.max_active_electrodes = max_active_electrodes
        self.leadfield_path = "/mesh_leadfield/leadfields/tdcs_leadfield"
        self.mesh_path = "/mesh_leadfield/"
        self.open_in_gmsh = open_in_gmsh
        self._mesh = None
        self._leadfield = None
        self._field_name = None
        self._field_units = None
        self.name = name
        self.target = target if target is not None else []
        self.avoid = avoid if avoid is not None else []

        self._assert_valid_currents()

    @property
    def lf_type(self):
        if self.mesh is None or self.leadfield is None:
            return None
        if self.leadfield.shape[1] == self.mesh.nodes.nr:
            return "node"
        elif self.leadfield.shape[1] == self.mesh.elm.nr:
            return "element"
        else:
            raise ValueError(
                "Could not find if the leadfield is node- or " "element-based"
            )

    @property
    def leadfield(self):
        """Reads the leadfield from the HDF5 file"""
        if self._leadfield is None and self.leadfield_hdf is not None:
            with h5py.File(self.leadfield_hdf, "r") as f:
                self.leadfield = f[self.leadfield_path][:]

        return self._leadfield

    @leadfield.setter
    def leadfield(self, leadfield):
        if leadfield is not None:
            assert leadfield.ndim == 3, "leadfield should be 3 dimensional"
            assert (
                leadfield.shape[2] == 3
            ), "Size of last dimension of leadfield should be 3"
        self._leadfield = leadfield

    @property
    def mesh(self):
        if self._mesh is None and self.leadfield_hdf is not None:
            self.mesh = mesh_io.Msh.read_hdf5(self.leadfield_hdf, self.mesh_path)

        return self._mesh

    @mesh.setter
    def mesh(self, mesh):
        elm_type = np.unique(mesh.elm.elm_type)
        if len(elm_type) > 1:
            raise ValueError("Mesh has both tetrahedra and triangles")
        else:
            self._mesh = mesh

    @property
    def field_name(self):
        if self.leadfield_hdf is not None and self._field_name is None:
            try:
                with h5py.File(self.leadfield_hdf, "r") as f:
                    self.field_name = f[self.leadfield_path].attrs["field"]
            except:
                return "Field"

        if self._field_name is None:
            return "Field"
        else:
            return self._field_name

    @field_name.setter
    def field_name(self, field_name):
        self._field_name = field_name

    @property
    def field_units(self):
        if self.leadfield_hdf is not None and self._field_units is None:
            try:
                with h5py.File(self.leadfield_hdf, "r") as f:
                    self.field_units = f[self.leadfield_path].attrs["units"]
            except:
                return "Au"

        if self._field_units is None:
            return "Au"
        else:
            return self._field_units

    @field_units.setter
    def field_units(self, field_units):
        self._field_units = field_units

    def to_mat(self):
        """Makes a dictionary for saving a matlab structure with scipy.io.savemat()

        Returns
        --------------------
        dict
            Dictionaty for usage with scipy.io.savemat
        """
        mat = {}
        mat["type"] = "TDCSoptimize"
        mat["leadfield_hdf"] = remove_None(self.leadfield_hdf)
        mat["max_total_current"] = remove_None(self.max_total_current)
        mat["max_individual_current"] = remove_None(self.max_individual_current)
        mat["max_active_electrodes"] = remove_None(self.max_active_electrodes)
        mat["open_in_gmsh"] = remove_None(self.open_in_gmsh)
        mat["name"] = remove_None(self.name)
        mat["target"] = _save_TDCStarget_mat(self.target)
        mat["avoid"] = _save_TDCStarget_mat(self.avoid)
        return mat

    @classmethod
    def read_mat_struct(cls, mat):
        """Reads a .mat structure

        Parameters
        -----------
        mat: dict
            Dictionary from scipy.io.loadmat

        Returns
        ----------
        p: TDCSoptimize
            TDCSoptimize structure
        """
        t = cls()
        leadfield_hdf = try_to_read_matlab_field(
            mat, "leadfield_hdf", str, t.leadfield_hdf
        )
        max_total_current = try_to_read_matlab_field(
            mat, "max_total_current", float, t.max_total_current
        )
        max_individual_current = try_to_read_matlab_field(
            mat, "max_individual_current", float, t.max_individual_current
        )
        max_active_electrodes = try_to_read_matlab_field(
            mat, "max_active_electrodes", int, t.max_active_electrodes
        )
        open_in_gmsh = try_to_read_matlab_field(
            mat, "open_in_gmsh", bool, t.open_in_gmsh
        )
        name = try_to_read_matlab_field(mat, "name", str, t.name)
        target = []
        if len(mat["target"]) > 0:
            for t in mat["target"][0]:
                target_struct = TDCStarget.read_mat_struct(t)
                if target_struct is not None:
                    target.append(target_struct)
        if len(target) == 0:
            target = None

        avoid = []
        if len(mat["avoid"]) > 0:
            avoid = []
            for t in mat["avoid"][0]:
                avoid_struct = TDCSavoid.read_mat_struct(t)
                if avoid_struct is not None:
                    avoid.append(avoid_struct)
        if len(avoid) == 0:
            avoid = None

        return cls(
            leadfield_hdf,
            max_total_current,
            max_individual_current,
            max_active_electrodes,
            name,
            target,
            avoid,
            open_in_gmsh,
        )

    def get_weights(self):
        """Calculates the volumes or areas of the mesh associated with the leadfield"""
        assert self.mesh is not None, "Mesh not defined"
        if self.lf_type == "node":
            weights = self.mesh.nodes_volumes_or_areas().value
        elif self.lf_type == "element":
            weights = self.mesh.elements_volumes_and_areas().value
        else:
            raise ValueError("Cant calculate weights: mesh or leadfield not set")

        weights *= self._get_avoid_field()
        return weights

    def _get_avoid_field(self):
        fields = []
        for a in self.avoid:
            a.mesh = self.mesh
            a.lf_type = self.lf_type
            fields.append(a.avoid_field())

        if len(fields) > 0:
            total_field = np.ones_like(fields[0])
            for f in fields:
                total_field *= f
            return total_field

        else:
            return 1.0

    def add_target(self, target=None):
        """Adds a target to the current tDCS optimization

        Parameters:
        ------------
        target: TDCStarget (optional)
            TDCStarget structure to be added. Default: empty TDCStarget

        Returns:
        -----------
        target: TDCStarget
            TDCStarget added to the structure
        """
        if target is None:
            target = TDCStarget(mesh=self.mesh, lf_type=self.lf_type)
        self.target.append(target)
        return target

    def add_avoid(self, avoid=None):
        """Adds an avoid structure to the current tDCS optimization

        Parameters:
        ------------
        target: TDCStarget (optional)
            TDCStarget structure to be added. Default: empty TDCStarget

        Returns:
        -----------
        target: TDCStarget
            TDCStarget added to the structure
        """
        if avoid is None:
            avoid = TDCSavoid(mesh=self.mesh, lf_type=self.lf_type)
        self.avoid.append(avoid)
        return avoid

    def _assign_mesh_lf_type_to_target(self):
        for t in self.target:
            if t.mesh is None:
                t.mesh = self.mesh
            if t.lf_type is None:
                t.lf_type = self.lf_type
        for a in self.avoid:
            if a.mesh is None:
                a.mesh = self.mesh
            if a.lf_type is None:
                a.lf_type = self.lf_type

    def _assert_valid_currents(
            self,
            max_total_current: float | None = None,
            max_individual_current: float | None = None,
        ):
        if max_total_current is None:
            max_total_current = self.max_total_current
        if max_individual_current is None:
            max_individual_current = self.max_individual_current

        assert max_total_current > 0, f"`max_total_current` must be positive (got {max_total_current})"
        assert max_individual_current > 0, f"`max_individual_current` must be positive (got {max_individual_current})"
        assert max_total_current >= max_individual_current

    def optimize(self, fn_out_mesh=None, fn_out_csv=None):
        """Runs the optimization problem

        Parameters
        -------------
        fn_out_mesh: str
            If set, will write out the electric field and currents to the mesh

        fn_out_mesh: str
            If set, will write out the currents and electrode names to a CSV file


        Returns
        ------------
        currents: N_elec x 1 ndarray
            Optimized currents. The first value is the current in the reference electrode
        """
        assert len(self.target) > 0, "No target defined"
        assert self.leadfield is not None, "Leadfield not defined"
        assert self.mesh is not None, "Mesh not defined"
        if self.max_active_electrodes is not None:
            assert (
                self.max_active_electrodes > 1
            ), "The maximum number of active electrodes should be at least 2"

        self._assert_valid_currents()
        max_total_current = self.max_total_current
        max_individual_current = self.max_individual_current

        self._assign_mesh_lf_type_to_target()
        weights = self.get_weights()
        norm_constrained = [t.directions is None for t in self.target]

        # Angle-constrained optimization
        if any([t.max_angle is not None for t in self.target]):
            if len(self.target) > 1:
                raise ValueError("Can't apply angle constraints with multiple target")
            t = self.target[0]
            max_angle = t.max_angle
            indices, directions = t.get_indexes_and_directions()

            assert max_angle >= 0, "`max_angle` must be >= 0"
            if self.max_active_electrodes is None:
                opt_problem = TESLinearAngleConstrained(
                    indices,
                    directions,
                    t.intensity,
                    max_angle,
                    self.leadfield,
                    max_total_current,
                    max_individual_current,
                    weights=weights,
                    target_weights=t.get_weights(),
                )
            else:
                opt_problem = TESLinearAngleElecConstrained(
                    self.max_active_electrodes,
                    indices,
                    directions,
                    t.intensity,
                    max_angle,
                    self.leadfield,
                    max_total_current,
                    max_individual_current,
                    weights,
                    target_weights=t.get_weights(),
                )

        # Norm-constrained optimization
        elif any(norm_constrained):
            if not all(norm_constrained):
                raise ValueError("Can't mix norm and linear constrained optimization")
            if self.max_active_electrodes is None:
                opt_problem = TESNormConstrained(
                    self.leadfield, max_total_current, max_individual_current, weights
                )
            else:
                opt_problem = TESNormElecConstrained(
                    self.max_active_electrodes,
                    self.leadfield,
                    max_total_current,
                    max_individual_current,
                    weights,
                )
            for t in self.target:
                if t.intensity < 0:
                    raise ValueError("Intensity must be > 0")
                opt_problem.add_norm_constraint(
                    t.get_indexes_and_directions()[0], t.intensity, t.get_weights()
                )

        # Simple QP-style optimization
        else:
            if self.max_active_electrodes is None:
                opt_problem = TESLinearConstrained(
                    self.leadfield, max_total_current, max_individual_current, weights
                )

            else:
                opt_problem = TESLinearElecConstrained(
                    self.max_active_electrodes,
                    self.leadfield,
                    max_total_current,
                    max_individual_current,
                    weights,
                )

            for t in self.target:
                opt_problem.add_linear_constraint(
                    *t.get_indexes_and_directions(), t.intensity, t.get_weights()
                )

        currents = opt_problem.solve()

        logger.log(25, "\n" + self.summary(currents))

        if fn_out_mesh is not None:
            fn_out_mesh = os.path.abspath(fn_out_mesh)
            m = self.field_mesh(currents)
            m.write(fn_out_mesh)
            v = m.view()
            ## Configure view
            v.Mesh.SurfaceFaces = 0
            v.View[0].Visible = 1
            # Change vector type for target field
            offset = 2
            if self.lf_type == "node":
                offset = 3
            for i, t in enumerate(self.target):
                v.View[offset + i].VectorType = 4
                v.View[offset + i].ArrowSizeMax = 60
                v.View[offset + i].Visible = 1
            # Electrode geo file
            el_geo_fn = os.path.splitext(fn_out_mesh)[0] + "_el_currents.geo"
            self.electrode_geo(el_geo_fn, currents)
            v.add_merge(el_geo_fn)
            max_c = np.max(np.abs(currents))
            v.add_view(
                Visible=1,
                RangeType=2,
                ColorTable=gmsh_view._coolwarm_cm(),
                CustomMax=max_c,
                CustomMin=-max_c,
                PointSize=10,
            )
            v.write_opt(fn_out_mesh)
            if self.open_in_gmsh:
                mesh_io.open_in_gmsh(fn_out_mesh, True)

        if fn_out_csv is not None:
            self.write_currents_csv(currents, fn_out_csv)

        return currents

    def field(self, currents):
        """Outputs the electric fields caused by the current combination

        Parameters
        -----------
        currents: N_elec x 1 ndarray
            Currents going through each electrode, in A. Usually from the optimize
            method. The sum should be approximately zero

        Returns
        ----------
        E: simnibs.mesh.NodeData or simnibs.mesh.ElementData
            NodeData or ElementData with the field caused by the currents
        """

        assert np.isclose(np.sum(currents), 0, atol=1e-5), "Currents should sum to zero"
        E = np.einsum("ijk,i->jk", self.leadfield, currents[1:])

        if self.lf_type == "node":
            E = mesh_io.NodeData(E, self.field_name, mesh=self.mesh)

        if self.lf_type == "element":
            E = mesh_io.ElementData(E, self.field_name, mesh=self.mesh)

        return E

    def electrode_geo(
        self, fn_out, currents=None, mesh_elec=None, elec_tags=None, elec_positions=None
    ):
        """Creates a mesh with the electrodes and their currents

        Parameters
        ------------
        currents: N_elec x 1 ndarray (optional)
            Electric current values per electrode. Default: do not print currents
        mesh_elec: simnibs.mesh.Msh (optional)
            Mesh with the electrodes. Default: look for a mesh called mesh_electrodes in
            self.leadfield_hdf
        elec_tags: N_elec x 1 ndarray of ints (optional)
            Tags of the electrodes corresponding to each leadfield column. The first is
            the reference electrode. Default: load at the attribute electrode_tags in the
            leadfield dataset
        elec_positions: N_elec x 3 ndarray of floats (optional)
            Positions of the electrodes in the head. If mesh_elec is not defined, will
            create small sphres at those positions instead.
            Default: load at the attribute electrode_pos in the leadfield dataset

        """
        # First try to set the electrode visualizations using meshed electrodes
        if mesh_elec is None:
            if self.leadfield_hdf is not None:
                try:
                    mesh_elec = mesh_io.Msh.read_hdf5(
                        self.leadfield_hdf, "mesh_electrodes"
                    )
                except KeyError:
                    pass
            else:
                raise ValueError("Please define a mesh with the electrodes")

        if elec_tags is None and mesh_elec is not None:
            if self.leadfield_hdf is not None:
                with h5py.File(self.leadfield_hdf, "r") as f:
                    elec_tags = f[self.leadfield_path].attrs["electrode_tags"]
            else:
                raise ValueError("Please define the electrode tags")

        # If not, use point electrodes
        if mesh_elec is None and elec_positions is None:
            if self.leadfield_hdf is not None:
                with h5py.File(self.leadfield_hdf, "r") as f:
                    elec_positions = f[self.leadfield_path].attrs["electrode_pos"]
            else:
                raise ValueError("Please define the electrode positions")

        if mesh_elec is not None:
            elec_pos = self._electrode_geo_triangles(
                fn_out, currents, mesh_elec, elec_tags
            )
            # elec_pos is used for writing electrode names
        elif elec_positions is not None:
            self._electrode_geo_points(fn_out, currents, elec_positions)
            elec_pos = elec_positions
        else:
            raise ValueError("Neither mesh_elec nor elec_positions defined")
        if self.leadfield_hdf is not None:
            with h5py.File(self.leadfield_hdf, "r") as f:
                try:
                    elec_names = f[self.leadfield_path].attrs["electrode_names"]
                    elec_names = [
                        n.decode() if isinstance(n, bytes) else n for n in elec_names
                    ]
                except KeyError:
                    elec_names = None

            if elec_names is not None:
                mesh_io.write_geo_text(
                    elec_pos, elec_names, fn_out, name="electrode_names", mode="ba"
                )

    def _electrode_geo_triangles(self, fn_out, currents, mesh_elec, elec_tags):
        if currents is None:
            currents = np.ones(len(elec_tags))

        assert len(elec_tags) == len(currents), "Define one current per electrode"

        triangles = []
        values = []
        elec_pos = []
        bar = mesh_elec.elements_baricenters()
        norms = mesh_elec.triangle_normals()
        for t, c in zip(elec_tags, currents):
            triangles.append(mesh_elec.elm[mesh_elec.elm.tag1 == t, :3])
            values.append(c * np.ones(len(triangles[-1])))
            avg_norm = np.average(norms[mesh_elec.elm.tag1 == t], axis=0)
            pos = np.average(bar[mesh_elec.elm.tag1 == t], axis=0)
            pos += avg_norm * 4
            elec_pos.append(pos)

        triangles = np.concatenate(triangles, axis=0)
        values = np.concatenate(values, axis=0)
        elec_pos = np.vstack(elec_pos)
        mesh_io.write_geo_triangles(
            triangles - 1,
            mesh_elec.nodes.node_coord,
            fn_out,
            values,
            "electrode_currents",
        )

        return elec_pos

    def _electrode_geo_points(self, fn_out, currents, elec_positions):
        if currents is None:
            currents = np.ones(len(elec_positions))

        assert len(elec_positions) == len(currents), "Define one current per electrode"
        mesh_io.write_geo_spheres(
            elec_positions, fn_out, currents, "electrode_currents"
        )

    def field_mesh(self, currents):
        """Creates showing the targets and the field
        Parameters
        -------------
        currents: N_elec x 1 ndarray
            Currents going through each electrode, in A. Usually from the optimize
            method. The sum should be approximately zero

        Returns
        ---------
        results: simnibs.msh.mesh_io.Msh
            Mesh file
        """
        target_fields = [
            t.as_field("target_{0}".format(i + 1)) for i, t in enumerate(self.target)
        ]
        weight_fields = [
            t.as_field("avoid_{0}".format(i + 1)) for i, t in enumerate(self.avoid)
        ]
        e_field = self.field(currents)
        e_magn_field = e_field.norm()
        if self.lf_type == "node":
            normals = -self.mesh.nodes_normals()[:]
            e_normal_field = np.sum(e_field[:] * normals, axis=1)
            e_normal_field = mesh_io.NodeData(
                e_normal_field, "normal" + e_field.field_name, mesh=self.mesh
            )
        m = copy.deepcopy(self.mesh)
        if self.lf_type == "node":
            m.nodedata = (
                [e_magn_field, e_field, e_normal_field] + target_fields + weight_fields
            )
        elif self.lf_type == "element":
            m.elmdata = [e_magn_field, e_field] + target_fields + weight_fields
        return m

    def write_currents_csv(self, currents, fn_csv, electrode_names=None):
        """Writes the currents and the corresponding electrode names to a CSV file

        Parameters
        ------------
        currents: N_elec x 1 ndarray
            Array with electrode currents
        fn_csv: str
            Name of CSV file to write
        electrode_names: list of strings (optional)
            Name of electrodes. Default: will read from the electrode_names attribute in
            the leadfield dataset
        """
        if electrode_names is None:
            if self.leadfield_hdf is not None:
                with h5py.File(self.leadfield_hdf, "r") as f:
                    electrode_names = f[self.leadfield_path].attrs["electrode_names"]
                    electrode_names = [
                        n.decode() if isinstance(n, bytes) else n
                        for n in electrode_names
                    ]
            else:
                raise ValueError("Please define the electrode names")

        assert len(electrode_names) == len(currents)
        with open(fn_csv, "w", newline="") as f:
            writer = csv.writer(f)
            for n, c in zip(electrode_names, currents):
                writer.writerow([n, c])

    def run(self, cpus=1):
        """Interface to use with the run_simnibs function

        Parameters
        ---------------
        cpus: int (optional)
            Does not do anything, it is just here for the common interface with the
            simulation's run function
        """
        if not self.name:
            if self.leadfield_hdf is not None:
                try:
                    name = re.search(r"(.+)_leadfield_", self.leadfield_hdf).group(1)
                except AttributeError:
                    name = "optimization"
            else:
                name = "optimization"
        else:
            name = self.name
        out_folder = os.path.dirname(name)
        os.makedirs(out_folder, exist_ok=True)

        # Set-up logger
        fh = logging.FileHandler(name + ".log", mode="w")
        formatter = logging.Formatter(
            "[ %(name)s - %(asctime)s - %(process)d ]%(levelname)s: %(message)s"
        )
        fh.setFormatter(formatter)
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)

        fn_summary = name + "_summary.txt"
        fh_s = logging.FileHandler(fn_summary, mode="w")
        fh_s.setFormatter(logging.Formatter("%(message)s"))
        fh_s.setLevel(25)
        logger.addHandler(fh_s)

        fn_out_mesh = name + ".msh"
        fn_out_csv = name + ".csv"
        logger.info("Optimizing")
        logger.log(25, str(self))
        self.optimize(fn_out_mesh, fn_out_csv)
        logger.log(
            25,
            "\n=====================================\n"
            "SimNIBS finished running optimization\n"
            "Mesh file: {0}\n"
            "CSV file: {1}\n"
            "Summary file: {2}\n"
            "=====================================".format(
                fn_out_mesh, fn_out_csv, fn_summary
            ),
        )

        logger.removeHandler(fh)
        logger.removeHandler(fh_s)

        return fn_out_mesh

    def __str__(self):
        s = "Optimization set-up\n"
        s += "===========================\n"
        s += "Leadfield file: {0}\n".format(self.leadfield_hdf)
        s += "Max. total current: {0} (A)\n".format(self.max_total_current)
        s += "Max. individual current: {0} (A)\n".format(self.max_individual_current)
        s += "Max. active electrodes: {0}\n".format(self.max_active_electrodes)
        s += "Name: {0}\n".format(self.name)
        s += "----------------------\n"
        s += "N targets: {0}\n".format(len(self.target))
        s += "......................\n".join(
            [
                "Target {0}:\n{1}".format(i + 1, str(t))
                for i, t in enumerate(self.target)
            ]
        )
        s += "----------------------\n"
        s += "N avoid: {0}\n".format(len(self.avoid))
        s += "......................\n".join(
            ["Avoid {0}:\n{1}".format(i + 1, str(t)) for i, t in enumerate(self.avoid)]
        )
        return s

    def summary(self, currents):
        """Returns a string with a summary of the optimization

        Parameters
        ------------
        field: ElementData or NodeData
            Field of interest

        Returns
        ------------
        summary: str
            Summary of field
        """
        s = "Optimization Summary\n"
        s += "=============================\n"
        s += "Total current: {0:.2e} (A)\n".format(np.linalg.norm(currents, ord=1) / 2)
        s += "Maximum current: {0:.2e} (A)\n".format(np.max(np.abs(currents)))
        s += "Active electrodes: {0}\n".format(int(np.linalg.norm(currents, ord=0)))
        field = self.field(currents)
        s += "Field Summary\n"
        s += "----------------------------\n"
        s += "Peak Value (99.9 percentile): {0:.2f} ({1})\n".format(
            field.get_percentiles(99.9)[0], self.field_units
        )
        s += "Mean field magnitude: {0:.2e} ({1})\n".format(
            field.mean_field_norm(), self.field_units
        )
        if np.any(self.mesh.elm.elm_type == 4):
            v_units = "mm3"
        else:
            v_units = "mm2"

        s += f"Focality: 50 %: {field.get_focality(cuttofs=50, peak_percentile=99.9)[0]:.2e} ({v_units})\n"
        s += f"          70 %: {field.get_focality(cuttofs=70, peak_percentile=99.9)[0]:.2e} ({v_units})\n"

        for i, t in enumerate(self.target):
            s += "Target {0}\n".format(i + 1)
            s += f"    Intensity specified: {t.intensity:10.4f} ({self.field_units})\n"
            s += f"    Intensity achieved:  {t.mean_intensity(field):10.4f} ({self.field_units})\n"
            if t.max_angle is not None:
                s += (
                    "    Average angle across target: {0:.1f} "
                    "(max set to {1:.1f}) (degrees)\n".format(
                        t.mean_angle(field), t.max_angle
                    )
                )
            else:
                s += "    Average angle across target: {0:.1f} (degrees)\n".format(
                    t.mean_angle(field)
                )

        for i, a in enumerate(self.avoid):
            s += "Avoid {0}\n".format(i + 1)
            s += "    Mean field magnitude in region: {0:.2e} ({1})\n".format(
                a.mean_field_norm_in_region(field), self.field_units
            )

        return s


class TDCStarget:
    """Defines a target for TDCS optimization

    Attributes
    -------------
    positions: Nx3 ndarray
        List of target positions, in x, y, z coordinates and in the subject space. Will find the closes mesh points
    indexes: Nx1 ndarray of ints
        Indexes (1-based) of elements/nodes for optimization. Overwrites positions
    directions: Nx3 ndarray
        List of Electric field directions to be optimized for each mesh point, the string
        'normal' or None (for magnitude optimization), Default: 'normal'
    intensity: float (optional)
        Target intensity of the electric field component in V/m. Default: 0.2
    max_angle: float (optional)
        Maximum angle between electric field and target direction, in degrees. Default:
        No maximum
    radius: float (optional)
        Radius of target. All the elements/nodes within the given radies of the indexes
        will be included.
    tissues: list or None (Optional)
        Tissues included in the target. Either a list of integer with tissue tags or None
        for all tissues. Default: None

    THE ONES BELOW SHOULD NOT BE FILLED BY THE USERS IN NORMAL CIRCUNTANCES:
    mesh: simnibs.msh.mesh_io.Msh (optional)
        Mesh where the target is defined. Set by the TDCSoptimize methods
    lf_type: 'node' or 'element'
        Where the electric field values are defined

    """

    def __init__(
        self,
        positions=None,
        indexes=None,
        directions="normal",
        intensity=0.2,
        max_angle=None,
        radius=2,
        tissues=None,
        mesh=None,
        lf_type=None,
    ):

        self.lf_type = lf_type
        self.mesh = mesh
        self.radius = radius
        self.tissues = tissues
        self.positions = positions
        self.indexes = indexes
        self.intensity = intensity
        self.max_angle = max_angle
        self.directions = directions

    @property
    def directions(self):
        return self._directions

    @directions.setter
    def directions(self, value):
        if value == "normal":
            pass
        elif value == "none":
            value = None
        elif isinstance(value, str):
            raise ValueError(
                "Invalid value for directions: f{directions} "
                'valid arguments are "normal", "none" or an array'
            )
        if value is None and self.max_angle is not None:
            raise ValueError("Can't constrain angle in magnitude optimizations")
        self._directions = value

    @classmethod
    def read_mat_struct(cls, mat):
        """Reads a .mat structure

        Parameters
        -----------
        mat: dict
            Dictionary from scipy.io.loadmat

        Returns
        ----------
        t: TDCStarget
            TDCStarget structure
        """
        t = cls()
        positions = try_to_read_matlab_field(mat, "positions", list, t.positions)
        indexes = try_to_read_matlab_field(mat, "indexes", list, t.indexes)
        directions = try_to_read_matlab_field(mat, "directions", list, t.directions)
        try:
            directions[0]
        except IndexError:
            directions = "normal"
        else:
            if isinstance(directions[0], str):
                directions = "".join(directions)
            if isinstance(directions[0], bytes):
                directions = "".join([d.decode() for d in directions])
        intensity = try_to_read_matlab_field(mat, "intensity", float, t.intensity)
        max_angle = try_to_read_matlab_field(mat, "max_angle", float, t.max_angle)
        radius = try_to_read_matlab_field(mat, "radius", float, t.radius)
        tissues = try_to_read_matlab_field(mat, "tissues", list, t.tissues)
        if positions is not None and len(positions) == 0:
            positions = None
        if indexes is not None and len(indexes) == 0:
            indexes = None
        if tissues is not None and len(tissues) == 0:
            tissues = None

        is_empty = True
        is_empty *= t.positions == positions
        is_empty *= t.indexes == indexes
        is_empty *= t.directions == directions
        is_empty *= t.intensity == intensity
        is_empty *= t.max_angle == max_angle
        is_empty *= t.tissues == tissues
        if is_empty:
            return None

        return cls(
            positions, indexes, directions, intensity, max_angle, radius, tissues
        )

    def get_weights(self):
        assert self.lf_type is not None, "Please set a lf_type"

        if self.lf_type == "node":
            weights = self.mesh.nodes_volumes_or_areas().value
        elif self.lf_type == "element":
            weights = self.mesh.elements_volumes_and_areas().value
        else:
            raise ValueError(
                "Invalid lf_type: {0}, should be "
                '"element" or "node"'.format(self.lf_type)
            )

        return weights

    def get_indexes_and_directions(self):
        """Calculates the mesh indexes and directions corresponding to this target
        Returns
        ----------
        indexes: (n,) ndarray of ints
            0-based region indexes

        indexes: (n,3) ndarray of floats
            Target directions
        """
        indexes, mapping = _find_indexes(
            self.mesh,
            self.lf_type,
            positions=self.positions,
            indexes=self.indexes,
            tissues=self.tissues,
            radius=self.radius,
        )

        directions = _find_directions(
            self.mesh, self.lf_type, self.directions, indexes, mapping
        )

        return indexes - 1, directions

    def as_field(self, name="target_field"):
        """Returns the target as an ElementData or NodeData field

        Parameters
        -----------
        name: str
            Name of the field. Default: 'target_field'
        Returns
        ---------
        target: ElementData or NodeData
            A vector field with a vector pointing in the given direction in the target
        """
        if (self.positions is None) == (self.indexes is None):  # negative XOR operation
            raise ValueError("Please set either positions or indexes")

        assert self.mesh is not None, "Please set a mesh"

        if self.directions is None:
            nr_comp = 1
        else:
            nr_comp = 3

        if self.lf_type == "node":
            field = np.zeros((self.mesh.nodes.nr, nr_comp))
            field_type = mesh_io.NodeData
        elif self.lf_type == "element":
            field = np.zeros((self.mesh.elm.nr, nr_comp))
            field_type = mesh_io.ElementData
        else:
            raise ValueError(
                "lf_type must be 'node' or 'element'."
                " Got: {0} instead".format(self.lf_type)
            )

        indexes, mapping = _find_indexes(
            self.mesh,
            self.lf_type,
            positions=self.positions,
            indexes=self.indexes,
            tissues=self.tissues,
            radius=self.radius,
        )

        if self.directions is None:
            field[indexes - 1] = self.intensity
        else:
            directions = _find_directions(
                self.mesh, self.lf_type, self.directions, indexes, mapping
            )
            field[indexes - 1] = directions * self.intensity

        return field_type(field, name, mesh=self.mesh)

    def mean_intensity(self, field):
        """Calculates the mean intensity of the given field in this target

        Parameters
        -----------
        field: Nx3 NodeData or ElementData
            Electric field

        Returns
        ------------
        intensity: float
            Mean intensity in this target and in the target direction
        """
        if (self.positions is None) == (self.indexes is None):  # negative XOR operation
            raise ValueError("Please set either positions or indexes")

        assert self.mesh is not None, "Please set a mesh"
        assert field.nr_comp == 3, "Field must have 3 components"

        indexes, mapping = _find_indexes(
            self.mesh,
            self.lf_type,
            positions=self.positions,
            indexes=self.indexes,
            tissues=self.tissues,
            radius=self.radius,
        )

        f = field[indexes]
        if self.directions is None:
            components = np.linalg.norm(f, axis=1)

        else:
            directions = _find_directions(
                self.mesh, self.lf_type, self.directions, indexes, mapping
            )

            components = np.sum(f * directions, axis=1)

        if self.lf_type == "node":
            weights = self.mesh.nodes_volumes_or_areas()[indexes]
        elif self.lf_type == "element":
            weights = self.mesh.elements_volumes_and_areas()[indexes]
        else:
            raise ValueError(
                "lf_type must be 'node' or 'element'."
                " Got: {0} instead".format(self.lf_type)
            )

        return np.average(components, weights=weights)

    def mean_angle(self, field):
        """Calculates the mean angle between the field and the target

        Parameters
        -----------
        field: Nx3 NodeData or ElementData
            Electric field

        Returns
        ------------
        angle: float
            Mean angle in this target between the field and the target direction, in
            degrees
        """
        if (self.positions is None) == (self.indexes is None):  # negative XOR operation
            raise ValueError("Please set either positions or indexes")

        assert self.mesh is not None, "Please set a mesh"
        assert field.nr_comp == 3, "Field must have 3 components"
        if self.directions is None:
            return np.nan

        indexes, mapping = _find_indexes(
            self.mesh,
            self.lf_type,
            positions=self.positions,
            indexes=self.indexes,
            tissues=self.tissues,
            radius=self.radius,
        )

        directions = _find_directions(
            self.mesh, self.lf_type, self.directions, indexes, mapping
        )
        if self.intensity < 0:
            directions *= -1
        f = field[indexes]
        components = np.sum(f * directions, axis=1)
        norm = np.linalg.norm(f, axis=1)
        tangent = np.sqrt(norm**2 - components**2)
        angles = np.rad2deg(np.arctan2(tangent, components))
        if self.lf_type == "node":
            weights = self.mesh.nodes_volumes_or_areas()[indexes]
        elif self.lf_type == "element":
            weights = self.mesh.elements_volumes_and_areas()[indexes]
        else:
            raise ValueError(
                "lf_type must be 'node' or 'element'."
                " Got: {0} instead".format(self.lf_type)
            )
        weights *= norm
        return np.average(angles, weights=weights)

    def __str__(self):
        s = (
            "positions: {0}\n"
            "indexes: {1}\n"
            "directions: {2}\n"
            "radius: {3}\n"
            "intensity: {4}\n"
            "max_angle: {5}\n"
            "tissues: {6}\n".format(
                str(self.positions),
                str(self.indexes),
                str(self.directions),
                self.radius,
                self.intensity,
                str(self.max_angle),
                str(self.tissues),
            )
        )
        return s


class TDCSavoid:
    """List of positions to be avoided by optimizer

    Attributes
    -------------
    positions: Nx3 ndarray
        List of positions to be avoided, in x, y, z coordinates and in the subject space.
        Will find the closest mesh points
    indexes: Nx1 ndarray of ints
        Indexes (1-based) of elements/nodes to be avoided. Overwrites positions
    weight : float (optional)
        Weight to give to avoid region. The larger, the more we try to avoid it. Default:
        1e3
    radius: float (optional)
        Radius of region. All the elements/nodes within the given radius of the indexes
        will be included.
    tissues: list or None (Optional)
        Tissues to be included in the region. Either a list of integer with tissue tags or None
        for all tissues. Default: None


    Note
    -------
    If both positions and indexes are set to None, and a tissue is set, it will set the
    given weith to all elements/nodes in the given tissues

    THE ONES BELLOW SHOULD NOT BE FILLED BY THE USERS IN NORMAL CIRCUNTANCES:

    mesh: simnibs.msh.mesh_io.Msh (optional)
        Mesh where the target is defined. Set by the TDCSoptimize methods
    lf_type: 'node' or 'element'
        Where the electric field values are defined

    Warning
    -----------
    Changing positions constructing the class
    can cause unexpected behaviour
    """

    def __init__(
        self,
        positions=None,
        indexes=None,
        weight=1e3,
        radius=2,
        tissues=None,
        mesh=None,
        lf_type=None,
    ):
        self.lf_type = lf_type
        self.mesh = mesh
        self.radius = radius
        self.tissues = tissues
        self.positions = positions
        self.indexes = indexes
        self.weight = weight

    @classmethod
    def read_mat_struct(cls, mat):
        """Reads a .mat structure

        Parameters
        -----------
        mat: dict
            Dictionary from scipy.io.loadmat

        Returns
        ----------
        t: TDCSavoid
            TDCSavoid structure
        """
        t = cls()
        positions = try_to_read_matlab_field(mat, "positions", list, t.positions)
        indexes = try_to_read_matlab_field(mat, "indexes", list, t.indexes)
        weight = try_to_read_matlab_field(mat, "weight", float, t.weight)
        radius = try_to_read_matlab_field(mat, "radius", float, t.radius)
        tissues = try_to_read_matlab_field(mat, "tissues", list, t.tissues)
        if positions is not None and len(positions) == 0:
            positions = None
        if indexes is not None and len(indexes) == 0:
            indexes = None
        if tissues is not None and len(tissues) == 0:
            tissues = None

        is_empty = True
        is_empty *= t.positions == positions
        is_empty *= t.indexes == indexes
        is_empty *= t.weight == weight
        is_empty *= t.tissues == tissues
        if is_empty:
            return None

        return cls(positions, indexes, weight, radius, tissues)

    def _get_avoid_region(self):
        if (self.indexes is not None) or (self.positions is not None):
            indexes, _ = _find_indexes(
                self.mesh,
                self.lf_type,
                positions=self.positions,
                indexes=self.indexes,
                tissues=self.tissues,
                radius=self.radius,
            )
            return indexes
        elif self.tissues is not None:
            if self.lf_type == "element":
                return self.mesh.elm.elm_number[
                    np.isin(self.mesh.elm.tag1, self.tissues)
                ]
            elif self.lf_type == "node":
                return self.mesh.elm.nodes_with_tag(self.tissues)
        else:
            raise ValueError("Please define either indexes/positions or tissues")

    def avoid_field(self):
        """Returns a field with self.weight in the target area and
        weight=1 outside the target area

        Returns
        ------------
        w: float, >= 1
            Weight field
        """
        assert self.mesh is not None, "Please set a mesh"
        assert self.lf_type is not None, "Please set a lf_type"
        assert self.weight >= 0, "Weights must be >= 0"
        if self.lf_type == "node":
            f = np.ones(self.mesh.nodes.nr)
        elif self.lf_type == "element":
            f = np.ones(self.mesh.elm.nr)
        else:
            raise ValueError(
                "lf_type must be 'node' or 'element'."
                " Got: {0} instead".format(self.lf_type)
            )

        indexes = self._get_avoid_region()
        f[indexes - 1] = self.weight
        if len(indexes) == 0:
            raise ValueError("Empty avoid region!")

        return f

    def as_field(self, name="weights"):
        """Returns a NodeData or ElementData field with the weights

        Paramets
        ---------
        name: str (optional)
            Name for the field

        Returns
        --------
        f: NodeData or ElementData
            Field with weights
        """
        w = self.avoid_field()
        if self.lf_type == "node":
            return mesh_io.NodeData(w, name, mesh=self.mesh)
        elif self.lf_type == "element":
            return mesh_io.ElementData(w, name, mesh=self.mesh)

    def mean_field_norm_in_region(self, field):
        """Calculates the mean field magnitude in the region defined by the avoid structure

        Parameters
        -----------
        field: ElementData or NodeData
            Field for which we calculate the mean magnitude
        """
        assert self.mesh is not None, "Please set a mesh"
        assert self.lf_type is not None, "Please set a lf_type"
        indexes = self._get_avoid_region()
        v = np.linalg.norm(field[indexes], axis=1)
        if self.lf_type == "node":
            weight = self.mesh.nodes_volumes_or_areas()[indexes]
        elif self.lf_type == "element":
            weight = self.mesh.elements_volumes_and_areas()[indexes]
        else:
            raise ValueError(
                "lf_type must be 'node' or 'element'."
                " Got: {0} instead".format(self.lf_type)
            )

        return np.average(v, weights=weight)

    def __str__(self):
        s = (
            "positions: {0}\n"
            "indexes: {1}\n"
            "radius: {2}\n"
            "weight: {3:.1e}\n"
            "tissues: {4}\n".format(
                str(self.positions),
                str(self.indexes),
                self.radius,
                self.weight,
                str(self.tissues),
            )
        )
        return s


class TDCSDistributedOptimize:
    """Defines a tdcs optimization problem with distributed sources

    This function uses the problem setup from

    Ruffini et al. "Optimization of multifocal transcranial current
    stimulation for weighted cortical pattern targeting from realistic modeling of
    electric fields", NeuroImage, 2014

    And the algorithm from

    Saturnino et al. "Accessibility of cortical regions to focal TES:
    Dependence on spatial position, safety, and practical constraints."
    NeuroImage, 2019

    Parameters
    --------------
    leadfield_hdf: str (optional)
        Name of file with leadfield
    max_total_current: float (optional)
        Maximum current across all electrodes (in Amperes). Default: 2e-3
    max_individual_current: float (optional)
        Maximum current for any single electrode (in Amperes). Default: 1e-3
    max_active_electrodes: int (optional)
        Maximum number of active electrodes. Default: no maximum
    name: str (optional)
        Name of optimization problem. Default: optimization
    target_image: str or pair (array, affine)
        Image to be "reproduced" via the optimization
    mni_space: bool (optional)
        Wether the image is in MNI space. Default True
    subpath: str (optional)
        Path to the subject "m2m" folder. Needed if mni_space=True
    intensity: float
        Target field intensity
    min_img_value: float >= 0 (optional)
        minimum image (for example t value) to be considered. Corresponds to T_min in
        Ruffini et al. 2014. Default: 0
    open_in_gmsh: bool (optional)
        Whether to open the result in Gmsh after the calculations. Default: False

    Attributes
    --------------
    leadfield_hdf: str
        Name of file with leadfield
    max_total_current: float (optional)
        Maximum current across all electrodes (in Amperes). Default: 2e-3
    max_individual_current: float
        Maximum current for any single electrode (in Amperes). Default: 1e-3
    max_active_electrodes: int
        Maximum number of active electrodes. Default: no maximum
    ledfield_path: str
        Path to the leadfield in the hdf5 file. Default: '/mesh_leadfield/leadfields/tdcs_leadfield'
    mesh_path: str
        Path to the mesh in the hdf5 file. Default: '/mesh_leadfield/'

    The two above are used to define:

    mesh: simnibs.msh.mesh_io.Msh
        Mesh with problem geometry

    leadfield: np.ndarray
        Leadfield matrix (N_elec -1 x M x 3) where M is either the number of nodes or the
        number of elements in the mesh. We assume that there is a reference electrode

    Alternatively, you can set the three attributes above and not leadfield_path,
    mesh_path and leadfield_hdf

    lf_type: None, 'node' or 'element'
        Type of leadfield.

    name: str
        Name for the optimization problem. Defaults tp 'optimization'

    target_image: str or pair (array, affine)
        Image to be "reproduced" via the optimization

    mni_space: bool (optional)
        Wether the image is in MNI space. Default True

    subpath: str (optional)
        Path to the subject "m2m" folder. Needed if mni_space=True

    intensity: float
        Target field intensity

    min_img_value: float >= 0 (optional)
        minimum image (for example t value) to be considered. Corresponds to T_min in
        Ruffini et al. 2014. Default: 0

    open_in_gmsh: bool (optional)
        Whether to open the result in Gmsh after the calculations. Default: False

    Warning
    -----------
    Changing leadfield_hdf, leadfield_path and mesh_path after constructing the class
    can cause unexpected behaviour
    """

    def __init__(
        self,
        leadfield_hdf=None,
        max_total_current=2e-3,
        max_individual_current=1e-3,
        max_active_electrodes=None,
        name="optimization/tdcs",
        target_image=None,
        mni_space=True,
        subpath=None,
        intensity=0.2,
        min_img_value=0,
        open_in_gmsh=True,
    ):

        self._tdcs_opt_obj = TDCSoptimize(
            leadfield_hdf=leadfield_hdf,
            max_total_current=max_total_current,
            max_individual_current=max_individual_current,
            max_active_electrodes=max_active_electrodes,
            name=name,
            target=[],
            avoid=[],
            open_in_gmsh=open_in_gmsh,
        )
        self.max_total_current = max_total_current
        self.max_individual_current = max_individual_current
        self.max_active_electrodes = max_active_electrodes
        self.leadfield_path = "/mesh_leadfield/leadfields/tdcs_leadfield"
        self.mesh_path = "/mesh_leadfield/"
        self.target_image = target_image
        self.mni_space = mni_space
        self.open_in_gmsh = open_in_gmsh
        self.subpath = subpath
        self.name = name

        self.intensity = intensity
        self.min_img_value = min_img_value

        if min_img_value < 0:
            raise ValueError("min_img_value must be > 0")

        self._tdcs_opt_obj._assert_valid_currents(
            self.max_total_current, self.max_individual_current
        )

    @property
    def lf_type(self):
        self._tdcs_opt_obj.mesh = self.mesh
        self._tdcs_opt_obj.leadfield = self.leadfield

        return self._tdcs_opt_obj.lf_type

    @property
    def leadfield_hdf(self):
        return self._tdcs_opt_obj.leadfield_hdf

    @leadfield_hdf.setter
    def leadfield_hdf(self, leadfield_hdf):
        self._tdcs_opt_obj.leadfield_hdf = leadfield_hdf

    @property
    def leadfield_path(self):
        return self._tdcs_opt_obj.leadfield_path

    @leadfield_path.setter
    def leadfield_path(self, leadfield_path):
        self._tdcs_opt_obj.leadfield_path = leadfield_path

    @property
    def mesh_path(self):
        return self._tdcs_opt_obj.mesh_path

    @mesh_path.setter
    def mesh_path(self, mesh_path):
        self._tdcs_opt_obj.mesh_path = mesh_path

    @property
    def name(self):
        return self._tdcs_opt_obj.name

    @name.setter
    def name(self, name):
        self._tdcs_opt_obj.name = name

    @property
    def leadfield(self):
        """Reads the leadfield from the HDF5 file"""
        self._tdcs_opt_obj.leadfield_hdf = self.leadfield_hdf
        return self._tdcs_opt_obj.leadfield

    @leadfield.setter
    def leadfield(self, leadfield):
        self._tdcs_opt_obj.leadfield = leadfield

    @property
    def mesh(self):
        self._tdcs_opt_obj.leadfield_hdf = self.leadfield_hdf
        return self._tdcs_opt_obj.mesh

    @mesh.setter
    def mesh(self, mesh):
        self._tdcs_opt_obj.mesh = mesh

    @property
    def field_name(self):
        self._tdcs_opt_obj.leadfield_hdf = self.leadfield_hdf
        return self._tdcs_opt_obj._field_name

    @field_name.setter
    def field_name(self, field_name):
        self._tdcs_opt_obj._field_name = field_name

    @property
    def field_units(self):
        self._tdcs_opt_obj.leadfield_hdf = self.leadfield_hdf
        return self._tdcs_opt_obj._field_units

    def to_mat(self):
        """Makes a dictionary for saving a matlab structure with scipy.io.savemat()

        Returns
        --------------------
        dict
            Dictionaty for usage with scipy.io.savemat
        """
        mat = {}
        mat["type"] = "TDCSDistributedOptimize"
        mat["leadfield_hdf"] = remove_None(self.leadfield_hdf)
        mat["max_total_current"] = remove_None(self.max_total_current)
        mat["max_individual_current"] = remove_None(self.max_individual_current)
        mat["max_active_electrodes"] = remove_None(self.max_active_electrodes)
        mat["open_in_gmsh"] = remove_None(self.open_in_gmsh)
        mat["name"] = remove_None(self.name)
        mat["target_image"] = remove_None(self.target_image)
        mat["mni_space"] = remove_None(self.mni_space)
        mat["subpath"] = remove_None(self.subpath)
        mat["intensity"] = remove_None(self.intensity)
        mat["min_img_value"] = remove_None(self.min_img_value)

        return mat

    @classmethod
    def read_mat_struct(cls, mat):
        """Reads a .mat structure

        Parameters
        -----------
        mat: dict
            Dictionary from scipy.io.loadmat

        Returns
        ----------
        p: TDCSoptimize
            TDCSoptimize structure
        """
        t = cls()
        leadfield_hdf = try_to_read_matlab_field(
            mat, "leadfield_hdf", str, t.leadfield_hdf
        )
        max_total_current = try_to_read_matlab_field(
            mat, "max_total_current", float, t.max_total_current
        )
        max_individual_current = try_to_read_matlab_field(
            mat, "max_individual_current", float, t.max_individual_current
        )
        max_active_electrodes = try_to_read_matlab_field(
            mat, "max_active_electrodes", int, t.max_active_electrodes
        )
        open_in_gmsh = try_to_read_matlab_field(
            mat, "open_in_gmsh", bool, t.open_in_gmsh
        )
        name = try_to_read_matlab_field(mat, "name", str, t.name)
        target_image = try_to_read_matlab_field(
            mat, "target_image", str, t.target_image
        )
        mni_space = try_to_read_matlab_field(mat, "mni_space", bool, t.mni_space)
        subpath = try_to_read_matlab_field(mat, "subpath", str, t.subpath)
        intensity = try_to_read_matlab_field(mat, "intensity", float, t.intensity)
        min_img_value = try_to_read_matlab_field(
            mat, "min_img_value", float, t.min_img_value
        )

        return cls(
            leadfield_hdf=leadfield_hdf,
            max_total_current=max_total_current,
            max_individual_current=max_individual_current,
            max_active_electrodes=max_active_electrodes,
            name=name,
            target_image=target_image,
            mni_space=mni_space,
            subpath=subpath,
            intensity=intensity,
            min_img_value=min_img_value,
            open_in_gmsh=open_in_gmsh,
        )

    def _target_distribution(self):
        """Gets the y and W fields, by interpolating the target_image

        Based on Eq. 1 from
        Ruffini et al. "Optimization of multifocal transcranial current
        stimulation for weighted cortical pattern targeting from realistic modeling of
        electric fields", NeuroImage, 2014
        """
        assert self.mesh is not None, "Please set a mesh"
        assert self.min_img_value >= 0, "min_img_value must be >= 0"
        assert self.intensity is not None, "intensity not set"
        # load image
        if isinstance(self.target_image, str):
            img = nibabel.load(self.target_image)
            vol = np.array(img.dataobj)
            affine = img.affine
        else:
            vol, affine = self.target_image
        vol = vol.squeeze()  # fix when image is "4D", i.e. NxMxKx1
        if vol.ndim != 3:
            raise ValueError("Target image has to be 3D")
        vol[np.isnan(vol)] = 0.0

        # if in MNI space, tranfrom coordinates
        if self.mni_space:
            if self.subpath is None:
                raise ValueError("subpath not set!")
            nodes_mni = transformations.subject2mni_coords(
                self.mesh.nodes[:], self.subpath
            )
            orig_nodes = np.copy(self.mesh.nodes[:])
            self.mesh.nodes.node_coord = nodes_mni
        # Interpolate
        if self.lf_type == "node":
            field = mesh_io.NodeData.from_data_grid(self.mesh, vol, affine)
        elif self.lf_type == "element":
            field = mesh_io.ElementData.from_data_grid(self.mesh, vol, affine)
        field = np.float64(field[:])

        # setting values in eyes to zero
        if np.any(self.mesh.elm.tag1 == ElementTags.EYE_BALLS_TH_SURFACE):
            logger.info("setting target values in eyes to zero")
            if self.lf_type == "node":
                eye_nodes = np.unique(
                    self.mesh.elm.node_number_list[
                        self.mesh.elm.tag1 == ElementTags.EYE_BALLS_TH_SURFACE, :
                    ]
                )
                eye_nodes = eye_nodes[eye_nodes > 0]
                field[eye_nodes - 1] = 0.0  # node indices in mesh are 1-based
            elif self.lf_type == "element":
                field[self.mesh.elm.tag1 == ElementTags.EYE_BALLS_TH_SURFACE] = 0.0

        if self.mni_space:
            self.mesh.nodes.node_coord = orig_nodes

        W = np.abs(field)
        W[np.abs(field) < self.min_img_value] = self.min_img_value
        y = field[:].copy()
        y[np.abs(field) < self.min_img_value] = 0
        y *= self.intensity

        if np.all(np.abs(field) < self.min_img_value):
            raise ValueError("Target image values are below min_img_value!")
        return y, W

    def normal_directions(self):
        assert self.mesh is not None, "Please set a mesh"
        assert self.lf_type is not None, "Please set a lf_type"

        if 4 in self.mesh.elm.elm_type:
            raise ValueError("Can't define a normal direction for volumetric data!")

        if self.lf_type == "node":
            normals = self.mesh.nodes_normals()[:]
        elif self.lf_type == "element":
            normals = self.mesh.triangle_normals()[:]

        return -normals

    def field(self, currents):
        """Outputs the electric fields caused by the current combination

        Parameters
        -----------
        currents: N_elec x 1 ndarray
            Currents going through each electrode, in A. Usually from the optimize
            method. The sum should be approximately zero

        Returns
        ----------
        E: simnibs.mesh.NodeData or simnibs.mesh.ElementData
            NodeData or ElementData with the field caused by the currents
        """
        return self._tdcs_opt_obj.field(currents)

    def field_mesh(self, currents):
        """Creates showing the targets and the field
        Parameters
        -------------
        currents: N_elec x 1 ndarray
            Currents going through each electrode, in A. Usually from the optimize
            method. The sum should be approximately zero

        Returns
        ---------
        results: simnibs.msh.mesh_io.Msh
            Mesh file
        """
        e_field = self.field(currents)
        e_magn_field = e_field.norm()
        normals = self.normal_directions()
        e_normal_field = np.sum(e_field[:] * normals, axis=1)
        target_map, W = self._target_distribution()
        erni = (target_map - W * e_normal_field) ** 2 - target_map**2
        erni *= len(target_map) / np.sum(W)

        m = copy.deepcopy(self.mesh)
        if self.lf_type == "node":
            add_field = m.add_node_field
        elif self.lf_type == "element":
            add_field = m.add_element_field

        add_field(e_field, e_field.field_name)
        add_field(e_magn_field, e_magn_field.field_name)
        add_field(e_normal_field, "normal" + e_field.field_name)
        add_field(target_map, "target_map")
        add_field(erni, "ERNI")
        return m

    def optimize(self, fn_out_mesh=None, fn_out_csv=None):
        """Runs the optimization problem

        Parameters
        -------------
        fn_out_mesh: str
            If set, will write out the electric field and currents to the mesh

        fn_out_mesh: str
            If set, will write out the currents and electrode names to a CSV file


        Returns
        ------------
        currents: N_elec x 1 ndarray
            Optimized currents. The first value is the current in the reference electrode
        """
        assert self.leadfield is not None, "Leadfield not defined"
        assert self.mesh is not None, "Mesh not defined"
        if self.max_active_electrodes is not None:
            assert (
                self.max_active_electrodes > 1
            ), "The maximum number of active electrodes should be at least 2"

        self._tdcs_opt_obj._assert_valid_currents(
            self.max_total_current, self.max_individual_current
        )
        max_total_current = self.max_total_current
        max_individual_current = self.max_individual_current

        assert self.min_img_value is not None, "min_img_value not set"
        assert self.intensity is not None, "intensity not set"

        y, W = self._target_distribution()
        normals = self.normal_directions()
        weights = np.sqrt(self._tdcs_opt_obj.get_weights())

        if self.max_active_electrodes is None:
            opt_problem = TESDistributed(
                W[None, :, None] * self.leadfield,
                y[:, None] * normals,
                weights[:, None] * normals,
                max_total_current,
                max_individual_current,
            )
        else:
            opt_problem = TESDistributedElecConstrained(
                self.max_active_electrodes,
                W[None, :, None] * self.leadfield,
                y[:, None] * normals,
                weights[:, None] * normals,
                max_total_current,
                max_individual_current,
            )

        currents = opt_problem.solve()

        logger.log(25, "\n" + self.summary(currents))

        if fn_out_mesh is not None:
            fn_out_mesh = os.path.abspath(fn_out_mesh)
            m = self.field_mesh(currents)
            m.write(fn_out_mesh)
            v = m.view()
            ## Configure view
            v.Mesh.SurfaceFaces = 0
            v.View[ElementTags.GM].Visible = 1
            # Electrode geo file
            el_geo_fn = os.path.splitext(fn_out_mesh)[0] + "_el_currents.geo"
            self._tdcs_opt_obj.electrode_geo(el_geo_fn, currents)
            v.add_merge(el_geo_fn)
            max_c = np.max(np.abs(currents))
            v.add_view(
                Visible=1,
                RangeType=2,
                ColorTable=gmsh_view._coolwarm_cm(),
                CustomMax=max_c,
                CustomMin=-max_c,
            )
            v.write_opt(fn_out_mesh)
            if self.open_in_gmsh:
                mesh_io.open_in_gmsh(fn_out_mesh, True)

        if fn_out_csv is not None:
            self._tdcs_opt_obj.write_currents_csv(currents, fn_out_csv)

        return currents

    def __str__(self):
        s = "Optimization set-up\n"
        s += "===========================\n"
        s += "Leadfield file: {0}\n".format(self.leadfield_hdf)
        s += "Max. total current: {0} (A)\n".format(self.max_total_current)
        s += "Max. individual current: {0} (A)\n".format(self.max_individual_current)
        s += "Max. active electrodes: {0}\n".format(self.max_active_electrodes)
        s += "Name: {0}\n".format(self.name)
        s += "----------------------\n"
        s += "Target image: {0}\n".format(self.target_image)
        s += "MNI space: {0}\n".format(self.mni_space)
        s += "Min. image value: {0}\n".format(self.min_img_value)
        s += "Target intensity: {0}\n".format(self.intensity)
        return s

    def summary(self, currents):
        """Returns a string with a summary of the optimization

        Parameters
        ------------
        field: ElementData or NodeData
            Field of interest

        Returns
        ------------
        summary: str
            Summary of field
        """
        s = self._tdcs_opt_obj.summary(currents)
        # Calculate erri
        field = self.field(currents)[:]
        normals = self.normal_directions()
        field_normal = np.sum(field * normals, axis=1)
        y, W = self._target_distribution()
        erri = np.sum((y - field_normal * W) ** 2 - y**2)
        erri *= len(y) / np.sum(W**2)
        # add Erri to messaga
        s += f"Error Relative to Non Intervention (ERNI): {erri:.2e}\n"
        return s

    def run(self, cpus=1):
        """Interface to use with the run_simnibs function

        Parameters
        ---------------
        cpus: int (optional)
            Does not do anything, it is just here for the common interface with the
            simulation's run function
        """
        return TDCSoptimize.run(self)


def _save_TDCStarget_mat(target):
    target_dt = np.dtype(
        [
            ("type", "O"),
            ("indexes", "O"),
            ("directions", "O"),
            ("positions", "O"),
            ("intensity", "O"),
            ("max_angle", "O"),
            ("radius", "O"),
            ("tissues", "O"),
        ]
    )

    target_mat = np.empty(len(target), dtype=target_dt)

    for i, t in enumerate(target):
        target_mat[i] = np.array(
            [
                (
                    "TDCStarget",
                    remove_None(t.indexes),
                    remove_None(t.directions),
                    remove_None(t.positions),
                    remove_None(t.intensity),
                    remove_None(t.max_angle),
                    remove_None(t.radius),
                    remove_None(t.tissues),
                )
            ],
            dtype=target_dt,
        )

    return target_mat


def _save_TDCSavoid_mat(avoid):
    avoid_dt = np.dtype(
        [
            ("type", "O"),
            ("indexes", "O"),
            ("positions", "O"),
            ("weight", "O"),
            ("radius", "O"),
            ("tissues", "O"),
        ]
    )

    avoid_mat = np.empty(len(avoid), dtype=avoid_dt)

    for i, t in enumerate(avoid):
        avoid_mat[i] = np.array(
            [
                (
                    "TDCSavoid",
                    remove_None(t.indexes),
                    remove_None(t.positions),
                    remove_None(t.weight),
                    remove_None(t.radius),
                    remove_None(t.tissues),
                )
            ],
            dtype=avoid_dt,
        )

    return avoid_mat


def _find_indexes(
    mesh, lf_type, indexes=None, positions=None, tissues=None, radius=0.0
):
    """Looks into the mesh to find either
        1. nodes/elements withn a given radius of a set of points (defined as positions)
        and in the specified tissues. The fist step will be to find the closest
        node/element
        2. Specific indexes
    Returns the indices of the nodes/elements in the mesh as well as a mapping saying
    from which of the oridinal points the new points were acquired"""

    if (positions is not None) == (indexes is not None):  # negative XOR operation
        raise ValueError("Please define either positions or indexes")

    if indexes is not None:
        indexes = np.atleast_1d(indexes)
        return indexes, np.arange(len(indexes))

    if lf_type == "node":
        if tissues is not None:
            mesh_indexes = mesh.elm.nodes_with_tag(tissues)
        else:
            mesh_indexes = mesh.nodes.node_number

        mesh_pos = mesh.nodes[mesh_indexes]

    elif lf_type == "element":
        if tissues is not None:
            mesh_indexes = mesh.elm.elm_number[np.isin(mesh.elm.tag1, tissues)]
        else:
            mesh_indexes = mesh.elm.elm_number

        mesh_pos = mesh.elements_baricenters()[mesh_indexes]

    else:
        raise ValueError('lf_type must be either "node" or "element"')

    assert radius >= 0.0, "radius should be >= 0"
    assert len(mesh_pos) > 0, "Could not find any elements or nodes with given tags"
    kdtree = scipy.spatial.cKDTree(mesh_pos)
    pos_projected, indexes = kdtree.query(positions)
    indexes = np.atleast_1d(indexes)
    if radius > 1e-9:
        in_radius = kdtree.query_ball_point(mesh_pos[indexes], radius)
        original = np.concatenate([(i,) * len(ir) for i, ir in enumerate(in_radius)])
        in_radius, uq_idx = np.unique(np.concatenate(in_radius), return_index=True)
        return mesh_indexes[in_radius], original[uq_idx]
    else:
        return mesh_indexes[indexes], np.arange(len(indexes))


def _find_directions(mesh, lf_type, directions, indexes, mapping=None):
    if directions is None:
        return None
    if directions == "normal":
        if 4 in np.unique(mesh.elm.elm_type):
            raise ValueError("Can't define a normal direction for volumetric data!")
        if lf_type == "node":
            directions = -mesh.nodes_normals()[indexes]
        elif lf_type == "element":
            directions = -mesh.triangle_normals()[indexes]
        return directions
    else:
        directions = np.atleast_2d(directions)
        if directions.shape[1] != 3:
            raise ValueError("directions must be the string 'normal' or a Nx3 array")
        if mapping is None:
            if len(directions) == len(indexes):
                mapping = np.arange(len(indexes))
            else:
                raise ValueError(
                    "Different number of indexes and directions and no "
                    "mapping defined"
                )
        elif len(directions) == 1:
            mapping = np.zeros(len(indexes), dtype=int)

        directions = directions / np.linalg.norm(directions, axis=1)[:, None]
        return directions[mapping]


class TESConstraints:
    def __init__(self, n, max_total_current, max_el_current):
        self.max_total_current = max_total_current
        self.max_el_current = max_el_current
        self.n = n

    def _bound_contraints(self):
        # The bound constraints are defined in the extended system
        # format C_.dot(x_) < d_
        C_ = np.vstack([np.eye(2 * self.n), -np.eye(2 * self.n)])
        d_ = np.hstack(
            [self.max_el_current * np.ones(2 * self.n), np.zeros(2 * self.n)]
        )
        return C_, d_

    def _l1_constraint(self):
        # The L1 constraints are defined in the extended system
        # format C_.dot(x_) < d_
        C_ = np.ones((1, 2 * self.n))
        d_ = np.array(2 * self.max_total_current)
        return C_, d_

    def _kirchhoff_constraint(self):
        # Kirchhoff's law defined in the extended system
        # A_.dot(x_) = b_
        A_ = np.hstack([np.ones((1, self.n)), -np.ones((1, self.n))])
        b_ = np.array([0.0])
        return A_, b_


class TESOptimizationProblem(TESConstraints):
    """Base class for TES Optimization problems. Provide basic functionality for
    manipulating leadfields and the quadratic portion of the objective

    Parameters
    -------------
    leadfield: N_elec x N_roi x N_comp ndarray
        Leadfield

    max_total_current: float
        Maximum total current flow through all electrodes

    max_el_current: float
        Maximum current flow through each electrode

    References
    ----------
    Saturnino (2019). Accessibility of cortical regions to focal TES:
        Dependence on spatial position, safety, and practical constraints.
        https://www.sciencedirect.com/science/article/pii/S1053811919307748
    """

    def __init__(
        self, leadfield, max_total_current=1e4, max_el_current=1e4, weights=None
    ):
        # + 1 because we add the reference
        super().__init__(leadfield.shape[0] + 1, max_total_current, max_el_current)
        self.n_elec_no_ref = leadfield.shape[0]
        self.n_target = leadfield.shape[1]
        self.leadfield = leadfield
        self.weights = np.ones(self.n_target) if weights is None else np.array(weights)

        assert (
            self.weights.shape[0] == self.n_target
        ), "Please define one weight per leadfield element"

        self.P = self._calc_P_mat()
        self.Q = self._calc_Q_mat()

    def _calc_P_mat(self):
        """Compute the matrix P which adds a reference electrode to the system
        (eq. 6).
        """
        return np.linalg.pinv(
            np.vstack([-np.ones(self.n_elec_no_ref), np.eye(self.n_elec_no_ref)])
        )

    def _calc_l_vec(self, target_indices, target_direction, weights):
        """Calculates the matrix `l` (eq. 14).

        This operator calculates the mean electric field in the target when
        applied to a vector of currents.

        In the paper,

        tau = target_indices
        t = weights
        n = target_direction

        """
        target_direction = np.atleast_2d(target_direction)
        if target_direction.shape[1] != 3:
            target_direction = target_direction.T
        assert target_direction.shape[1] == 3, "A direction must have 3 dimentions"

        target_indices = np.atleast_1d(target_indices)
        assert len(target_indices) == len(
            target_direction
        ), "Please define one direction per target"

        target_direction = target_direction / np.linalg.norm(
            target_direction, axis=1, keepdims=True
        )

        weights = weights[target_indices]
        leadfield = self.leadfield[:, target_indices] * weights[:, None]

        l = np.tensordot(leadfield, target_direction, axes=2) / weights.sum()

        return l @ self.P  # add a reference channel

    def _calc_Q_mat(
        self,
        indices: None | npt.NDArray = None,
        weights: None | npt.NDArray = None
    ):
        """Calculate the energy matrix, Q, for optimization

            x.dot(Q.dot(x)) = e

        where `x` is the electrode currents and `e` is the average squared
        electric field norm.

        Optionally, restrict calculations to elements in `indices`.

        Eq. 18 or 21 depending on whether `indices` and/or `weights` are
        specified or not.

        Parameters
        ----------
        indices :
            Target indices. If None, calculate for the whole domain.
        weights :
            Target weights. If None, use the weights supplied during
            initialization.

        Returns
        -------
        Q: np.ndarray
            Quadratic component
        """
        weights = self.weights if weights is None else np.atleast_1d(weights)
        if indices is not None:
            indices = np.atleast_1d(indices)
            weights = weights[indices]
            leadfield = self.leadfield[:, indices]
        else:
            leadfield = self.leadfield

        # Contract over the last/first two axes
        Q = (
            np.tensordot(
                leadfield,
                np.transpose(leadfield * weights[:, None], (1, 2, 0)),
                axes=2,
            )
            / weights.sum()
        )

        return self.P.T @ Q @ self.P  # add a reference channel

    def extend_currents(self, x):
        """
         Returns the extended version of the input currents x
         x_ext = [x, -x]
         with x_ext > 0
         This representation in very usefull when dealing with L1 constraints

        Parameters
        ------------
        x: ndarray
            Variable to be extended (normally currents)

        Returns
        ----------
        x_: ndarray
            Extended version
        """

        x_ = np.hstack([x, -x])
        x_[x_ < 0] = 0
        return x_

    def solve(self):
        raise NotImplementedError()


class TESLinearConstrained(TESOptimizationProblem):
    """Class for solving the TES Problem with linear constraints

    This corresponds to Problem 8 in Saturnino et al., 2019
    """

    def __init__(
        self, leadfield, max_total_current=1e5, max_el_current=1e5, weights=None
    ):

        super().__init__(leadfield, max_total_current, max_el_current, weights)
        self.l = np.empty((0, self.n), dtype=float)
        self.target_means = np.empty(0, dtype=float)

    def add_linear_constraint(
        self, target_indices, target_direction, target_mean, target_weights=None
    ):
        """Add a linear constrait to the problem

        Parameters
        ------------
        target_indices: list of ints
            Indices of targets.
        target_direction: ndarray
            The electric field direction to be optimized for each target position
        target_mean: float
            Target mean electric field in region
        target_weights: ndarray (optional)
            Weights (such are areas/volumes) for calculating the mean. Defined for every
            index
        """
        target_weights = self.weights if target_weights is None else target_weights

        l = self._calc_l_vec(target_indices, target_direction, target_weights)
        l *= np.sign(target_mean)
        self.l = np.vstack([self.l, l])

        self.target_means = np.hstack([self.target_means, np.abs(target_mean)])

    def solve(self, log_level=20):
        """Solves the optimization problem

        Returns
        ----------
        x: np.array
            Optimized currents
        """
        return _linear_constrained_tes_opt(
            self.l,
            self.target_means,
            self.Q,
            self.max_el_current,
            self.max_total_current,
            log_level=log_level,
        )


class TESLinearAngleConstrained(TESOptimizationProblem):
    """Class for solving the TES Problem with linear and angle contraints

    This corresponds to Problem 6 in Saturnino et al., 2019
    """

    def __init__(
        self,
        target_indices,
        target_direction,
        target_mean,
        max_angle,
        leadfield,
        max_total_current=1e5,
        max_el_current=1e5,
        weights=None,
        target_weights=None,
    ):

        super().__init__(leadfield, max_total_current, max_el_current, weights)
        target_weights = self.weights if target_weights is None else target_weights

        self.l = np.atleast_2d(
            self._calc_l_vec(target_indices, target_direction, target_weights)
        )
        self.Qnorm = self._calc_Q_mat(target_indices)
        self.target_mean = np.atleast_1d(target_mean)
        self.max_angle = max_angle

    def solve(self, log_level=20):
        return _linear_angle_constrained_tes_opt(
            self.l,
            self.target_mean,
            self.Q,
            self.max_el_current,
            self.max_total_current,
            self.Qnorm,
            self.max_angle,
            log_level=log_level,
        )


class TESLinearElecConstrained(TESLinearConstrained):
    """Class for solving the TES Problem with linear and number of electrodes
    constraints

    This corresponds to Problem 10 in Saturnino et al., 2019
    """

    def __init__(
        self, n_elec, leadfield, max_total_current=1e5, max_el_current=1e5, weights=None
    ):

        super().__init__(leadfield, max_total_current, max_el_current, weights)
        self.n_elec = n_elec

    def _solve_reduced(self, linear, quadratic, extra_ineq=None):
        l = linear[0]
        Q = quadratic[0]
        x = _linear_constrained_tes_opt(
            l,
            self.target_means,
            Q,
            self.max_el_current,
            self.max_total_current,
            extra_ineq=extra_ineq,
            log_level=10,
        )
        if np.any(l.dot(x) < self.target_means * 0.99):
            return x, 1e20
        else:
            return x, x.dot(Q).dot(x)

    def solve(
        self, log_level=20, eps_bb=1e-1, max_bb_iter=100, init_strategy="compact"
    ):
        # Heuristically eliminate electrodes
        max_el_current = min(self.max_el_current, self.max_total_current)
        el = np.arange(self.n)
        if init_strategy == "compact":
            x = _linear_constrained_tes_opt(
                self.l,
                self.target_means,
                self.Q,
                self.max_el_current,
                self.max_total_current,
                log_level=10,
            )
            active = np.abs(x) > 1e-3 * max_el_current
            if not np.any(active):
                active = np.ones(self.n, dtype=bool)
            init = bb_state([], el[~active].tolist(), el[active].tolist())
        elif init_strategy == "full":
            init = bb_state([], [], el.tolist())
        else:
            raise ValueError("Invalid initialization strategy")

        bounds_function = functools.partial(
            _bb_bounds_tes_problem,
            max_l0=self.n_elec,
            linear=[self.l],
            quadratic=[self.Q],
            max_el_current=max_el_current,
            func=self._solve_reduced,
        )

        final_state = _branch_and_bound(
            init, bounds_function, eps_bb, max_bb_iter, log_level=log_level
        )

        return final_state.x_ub


class TESLinearAngleElecConstrained(TESLinearAngleConstrained):
    """Class for solving the TES Problem with linear, angle and number of electrodes
    constraints

    This corresponds to Problem 7 in Saturnino et al., 2019
    """

    def __init__(
        self,
        n_elec,
        target_indices,
        target_direction,
        target_mean,
        max_angle,
        leadfield,
        max_total_current=1e5,
        max_el_current=1e5,
        weights=None,
        target_weights=None,
    ):
        super().__init__(
            target_indices,
            target_direction,
            target_mean,
            max_angle,
            leadfield,
            max_total_current,
            max_el_current,
            weights,
            target_weights,
        )

        self.n_elec = n_elec
        self._feasible = True

    def _solve_reduced(self, linear, quadratic, extra_ineq=None):
        l = linear[0]
        Q = quadratic[0]
        Qnorm = quadratic[1]

        x = _linear_angle_constrained_tes_opt(
            l,
            self.target_mean,
            Q,
            self.max_el_current,
            self.max_total_current,
            Qnorm,
            self.max_angle,
            extra_ineq=extra_ineq,
            log_level=10,
        )

        field = l.dot(x)
        if not self._feasible:
            return x, -field

        elif np.any(field < self.target_mean * 0.99):
            return x, 1e20

        else:
            return x, x.dot(Q).dot(x)

    def solve(
        self, log_level=20, eps_bb=1e-1, max_bb_iter=100, init_strategy="compact"
    ):
        # Heuristically eliminate electrodes
        max_el_current = min(self.max_el_current, self.max_total_current)
        el = np.arange(self.n)
        #  first determine if problem is feasible
        x = _linear_angle_constrained_tes_opt(
            self.l,
            self.target_mean,
            self.Q,
            self.max_el_current,
            self.max_total_current,
            self.Qnorm,
            self.max_angle,
            log_level=10,
        )

        feasible = np.allclose(self.l.dot(x), self.target_mean, rtol=1e-2)
        if not feasible:
            init_strategy = "full"

        if init_strategy == "compact":
            active = np.abs(x) > 1e-3 * max_el_current
            if not np.any(active):
                active = np.ones(self.n, dtype=bool)
            init = bb_state([], el[~active].tolist(), el[active].tolist())

        elif init_strategy == "full":
            init = bb_state([], [], el.tolist())

        else:
            raise ValueError("Invalid initialization strategy")

        bounds_function = functools.partial(
            _bb_bounds_tes_problem,
            max_l0=self.n_elec,
            linear=[self.l],
            quadratic=[self.Q, self.Qnorm],
            max_el_current=max_el_current,
            func=self._solve_reduced,
        )
        self._feasible = feasible
        final_state = _branch_and_bound(
            init, bounds_function, eps_bb, max_bb_iter, log_level=log_level
        )

        return final_state.x_ub


class TESNormConstrained(TESOptimizationProblem):
    """Class for solving the TES Problem with norm-type constraints"""

    def __init__(
        self, leadfield, max_total_current=1e5, max_el_current=1e5, weights=None
    ):
        super().__init__(leadfield, max_total_current, max_el_current, weights)
        self.Qnorm = np.empty((0, self.n, self.n), dtype=float)
        self.target_means = np.empty(0, dtype=float)

    def add_norm_constraint(self, target_indices, target_mean, target_weights=None):
        """Add a norm contraint of the type
        x^T Q x = t^2
        where x^T Q x is the squared field norm in target_indices and t is taret_mean

        Parameters
        ------------
        target_indices: list of ints
            Indices of targets.
        target_mean: float
            Target mean electric field norm in region
        target_weights: ndarray (optional)
            Weights (such are areas/volumes) for calculating the mean. Defined for every
            index
        """
        target_weights = self.weights if target_weights is None else target_weights

        Qnorm = self._calc_Q_mat(target_indices, target_weights)
        self.Qnorm = np.concatenate([self.Qnorm, Qnorm[None, ...]])

        self.target_means = np.hstack([self.target_means, np.abs(target_mean)])

    def solve(self, log_level=20):
        """Solves the optimization problem

        Returns
        ----------
        x: np.array
            Optimized currents
        """
        return _norm_constrained_tes_opt(
            self.Qnorm,
            self.target_means,
            self.Q,
            self.max_el_current,
            self.max_total_current,
            log_level=log_level,
        )


class TESNormElecConstrained(TESNormConstrained):
    """Class for solving the TES Problem with norm-type and number of electrodes
    constraints
    """

    def __init__(
        self, n_elec, leadfield, max_total_current=1e5, max_el_current=1e5, weights=None
    ):

        super().__init__(leadfield, max_total_current, max_el_current, weights)
        self.n_elec = n_elec

    def _solve_reduced(self, linear, quadratic, extra_ineq=None):
        Q = quadratic[0]
        Qnorm = np.stack(quadratic[1:])  # It should be stack here

        x = _norm_constrained_tes_opt(
            Qnorm,
            self.target_means,
            Q,
            self.max_el_current,
            self.max_total_current,
            log_level=10,
        )
        if np.any(np.sqrt(x.T.dot(Qnorm).dot(x)) < self.target_means * 0.99):
            return x, 1e20
        else:
            return x, x.dot(Q).dot(x)

    def solve(
        self, log_level=20, eps_bb=1e-1, max_bb_iter=100, init_strategy="compact"
    ):
        # Heuristically eliminate electrodes
        max_el_current = min(self.max_el_current, self.max_total_current)
        el = np.arange(self.n)
        if init_strategy == "compact":
            x = _norm_constrained_tes_opt(
                self.Qnorm,
                self.target_means,
                self.Q,
                self.max_el_current,
                self.max_total_current,
                log_level=10,
            )
            active = np.abs(x) > 1e-3 * max_el_current
            if not np.any(active):
                active = np.ones(self.n, dtype=bool)
            init = bb_state([], el[~active].tolist(), el[active].tolist())
        elif init_strategy == "full":
            init = bb_state([], [], el.tolist())
        else:
            raise ValueError("Invalid initialization strategy")

        bounds_function = functools.partial(
            _bb_bounds_tes_problem,
            max_l0=self.n_elec,
            linear=[np.zeros((1, self.n))],
            quadratic=np.concatenate([self.Q[None, ...], self.Qnorm]),
            max_el_current=max_el_current,
            func=self._solve_reduced,
        )

        final_state = _branch_and_bound(
            init, bounds_function, eps_bb, max_bb_iter, log_level=log_level
        )

        return final_state.x_ub


class TESDistributed(TESConstraints):
    """Class defining TES Optimization problems with distributed sources


    minimize (W(target_field - leadfield x))^2

    Parameters
    -------------
    leadfield: N_elec x N_roi x 3 ndarray
        Leadfield

    target_field: N_roi x 3
        Target electric field

    directions: N_roi x 3
        Field directions to be used for optimization

    max_total_current: float
        Maximum total current flow through all electrodes

    max_el_current: float
        Maximum current flow through each electrode

    weights: N_roi x 1 or N_roi x 3 ndarray
        Weight for each element / field component
    """

    def __init__(
        self,
        leadfield,
        target_field,
        weights=None,
        max_total_current=1e4,
        max_el_current=1e4,
    ):
        super().__init__(leadfield.shape[0] + 1, max_total_current, max_el_current)
        weights = np.ones(leadfield.shape[1]) if weights is None else weights
        self.l, self.Q = self._calc_l_Q(leadfield, target_field, weights)

    def _calc_l_Q(self, leadfield, target_field, weights):
        """Calculates the linear and quadratic parts of the optimization problem"""
        if weights.ndim == 1:
            Q = sum(
                leadfield[..., i].dot((leadfield[..., i] * weights**2).T)
                for i in range(leadfield.shape[2])
            )
            l = -2 * np.einsum(
                "ijk, jk -> i", leadfield, target_field * weights[:, None] ** 2
            )

        elif weights.ndim == 2 and weights.shape[1] == 3:
            A = np.einsum("ijk, jk -> ij", leadfield, weights)
            Q = A.dot(A.T)
            l = -2 * np.sum(target_field * weights, axis=1).dot(A.T)

        else:
            raise ValueError("Invalid shape for weights")

        # For numerical reasons
        P = np.linalg.pinv(np.vstack([-np.ones(len(l)), np.eye(len(l))]))
        l = l.dot(P)
        Q = P.T.dot(Q).dot(P)

        l /= np.sum(weights**2)
        Q /= np.sum(weights**2)
        return l, Q

    def solve(self, log_level=20):
        """Solves the optimization problem

        Returns
        ----------
        x: np.array
            Optimized currents
        """
        return _least_squares_tes_opt(
            self.l,
            self.Q,
            self.max_el_current,
            self.max_total_current,
            log_level=log_level,
        )


class TESDistributedElecConstrained(TESDistributed):
    """Class for solving the TES Distributed Problem with number of electrodes
    constraints


    """

    def __init__(
        self,
        n_elec,
        leadfield,
        target_field,
        weights=None,
        max_total_current=1e4,
        max_el_current=1e4,
    ):

        super().__init__(
            leadfield, target_field, weights, max_total_current, max_el_current
        )
        self.n_elec = n_elec

    def _solve_reduced(self, linear, quadratic, extra_ineq=None):
        l = linear[0]
        Q = quadratic[0]
        x = _least_squares_tes_opt(
            l,
            Q,
            self.max_el_current,
            self.max_total_current,
            extra_ineq=extra_ineq,
            log_level=10,
        )
        return x, l.dot(x) + x.dot(Q).dot(x)

    def solve(
        self, log_level=20, eps_bb=1e-1, max_bb_iter=500, init_strategy="compact"
    ):
        # Heuristically eliminate electrodes
        max_el_current = min(self.max_el_current, self.max_total_current)
        el = np.arange(self.n)
        if init_strategy == "compact":
            x = _least_squares_tes_opt(
                self.l,
                self.Q,
                self.max_el_current,
                self.max_total_current,
                log_level=10,
            )
            active = np.abs(x) > 1e-3 * max_el_current
            if not np.any(active):
                active = np.ones(self.n, dtype=bool)
            init = bb_state([], el[~active].tolist(), el[active].tolist())

        elif init_strategy == "full":
            init = bb_state([], [], el.tolist())

        else:
            raise ValueError("Invalid initialization strategy")

        bounds_function = functools.partial(
            _bb_bounds_tes_problem,
            max_l0=self.n_elec,
            linear=[self.l[None, :]],
            quadratic=[self.Q],
            max_el_current=max_el_current,
            func=self._solve_reduced,
        )

        final_state = _branch_and_bound(
            init, bounds_function, eps_bb, max_bb_iter, log_level=log_level
        )

        return final_state.x_ub


def _linear_constrained_tes_opt(
    l,
    target_mean,
    Q,
    max_el_current,
    max_total_current,
    extra_ineq=None,
    extra_eq=None,
    log_level=10,
):
    assert (
        l.shape[0] == target_mean.shape[0]
    ), "Please specify one target mean per target"
    assert l.shape[1] == Q.shape[0]

    scale = max_total_current

    target_mean = target_mean / scale
    max_el_current = max_el_current / scale
    max_total_current = max_total_current / scale

    n = l.shape[1]
    tes_constraints = TESConstraints(n, max_total_current, max_el_current)
    # First solve an LP to get a feasible starting point
    l_ = np.hstack([l, -l])

    C_, d_ = tes_constraints._l1_constraint()
    A_, b_ = tes_constraints._kirchhoff_constraint()

    # Inequality constraints
    if extra_ineq is not None:
        C_ = np.vstack([C_, extra_ineq[0]])
        d_ = np.hstack([d_, extra_ineq[1] / scale])

    if extra_eq is not None:
        A_ = np.vstack([A_, extra_eq[0]])
        b_ = np.hstack([b_, extra_eq[1] / scale])

    sol = scipy.optimize.linprog(
        -np.average(l_, axis=0),
        A_ub=np.vstack([C_, l_]),
        b_ub=np.hstack([d_, target_mean]),
        # constraints due to Kirchhoff's law
        A_eq=A_,
        b_eq=b_,
        bounds=(0, max_el_current),
        method="highs",
        options=dict(
            # Sets the primal/dual feasibility tolerances for `highs-ds`;
            # `highs-ipm` uses the minimum of the two
            dual_feasibility_tolerance=1e-8,
            primal_feasibility_tolerance=1e-8,
        ),
    )
    if sol.status != 0:
        raise RuntimeError(sol.message)

    x_ = sol.x

    # Test if the objective can be reached
    f = l.dot(x_[:n] - x_[n:])

    if np.any(np.abs(f - target_mean) >= np.abs(1e-2 * target_mean)):
        logger.log(log_level, "Could not reach target intensities")
        return (x_[:n] - x_[n:]) * scale

    logger.log(log_level, "Target intensity reached, optimizing focality")

    # Do the QP
    Q_ = np.block([[Q, -Q], [-Q, Q]])
    C_b, d_b = tes_constraints._bound_contraints()
    abs_x_max = np.abs(x_).max()

    x_ = _active_set_QP(
        l=np.zeros(2 * n),
        Q=Q_,
        C=np.vstack([C_b, C_]),
        d=np.hstack([d_b, d_]),
        x0=x_,
        A=np.vstack([A_, l_]),
        b=np.hstack([b_, f]),  # I use "f"
        tol_primal=1e-3 * abs_x_max * scale,
        tol_feasibility_x0=1e-3 * abs_x_max,
        tol_zero_div=1e-9 / scale,
    )

    return (x_[:n] - x_[n:]) * scale


def _calc_angle(x, Qin, l):
    tan = np.sqrt(np.abs(x.dot(Qin).dot(x) - l.dot(x) ** 2))
    return np.abs(np.arctan2(tan, l.dot(x)))[0]


def _linear_angle_constrained_tes_opt(
    l,
    target_mean,
    Q,
    max_el_current,
    max_total_current,
    Qin,
    max_angle,
    extra_ineq=None,
    extra_eq=None,
    eps_linear=1e-5,
    eps_angle=1e-1,
    log_level=20,
):

    max_angle = np.deg2rad(max_angle)
    logger.log(log_level, "Running optimization with angle constraint")
    max_iter = 20
    # Try to find where we can find values below and above the target
    it = 0
    above = 0

    # Check if the maximal focality solution alreay fulfills the constraint
    x = _linear_constrained_tes_opt(
        l,
        target_mean,
        Q,
        max_el_current,
        max_total_current,
        extra_ineq=extra_ineq,
        extra_eq=extra_eq,
        log_level=log_level - 10,
    )

    if _calc_angle(x, Qin, l) <= max_angle:
        logger.log(log_level, "Max focality solution fullfills angle constraint")
        return x

    x_above = np.copy(x)

    # calculate the smallest angle, given a fixed intensity
    def _minimize_angle(alpha):
        x_l = _linear_constrained_tes_opt(
            l,
            alpha * target_mean,
            Qin,
            max_el_current,
            max_total_current,
            extra_ineq=extra_ineq,
            extra_eq=extra_eq,
            log_level=log_level - 10,
        )
        return x_l

    x = _minimize_angle(1.0)
    angle = _calc_angle(x, Qin, l)

    # if we cannot reduce the angle to the target while keeting l^t x at the target
    # intensity, reduce the target intensity untill it's achievable.
    if angle > max_angle:
        logger.log(log_level, "Target intensity can't be reached, reducing it")
        above = 1.0
        angle_above = angle
        below = 0
        angle_below = 0
        it = 0
        # Use the secant method
        while not (angle > max_angle * (1 - eps_angle) and angle < max_angle):
            alpha = above + (max_angle * (1 - eps_angle * 0.5) - angle_above) * (
                below - above
            ) / (angle_below - angle_above)
            x = _minimize_angle(alpha)
            angle = _calc_angle(x, Qin, l)
            logger.log(
                log_level,
                "{0} alpha: {1:.3e}, angle: {2:.2e}, max_angle: {3:.2e}".format(
                    it, alpha, angle, max_angle
                ),
            )
            if angle < max_angle:
                below = alpha
                angle_below = angle
                x_below = np.copy(x)
            else:
                above = alpha
                angle_above = angle
            it += 1
            if it > max_iter:
                if below == 0:
                    return x
                else:
                    return x_below
        return x

    # In this case, we know that by minimizing x^t Qin x while keeping l^t x = t, we can
    # achieve the bound
    # find a combination between Q and Qin that maximizes focality while keeping Qin in
    # the bound
    else:
        logger.log(
            log_level,
            "Target intensity reached, optimizing focality with angle constraint",
        )
        angle_below = angle
        below = 1.0
        x_below = np.copy(x)
        angle_above = _calc_angle(x_above, Qin, l)
        above = 0

        # Start the secant method
        it = 0
        alpha = 1
        while not (angle > max_angle * (1 - eps_angle) and angle < max_angle):
            alpha = above + (max_angle * (1 - eps_angle * 0.5) - angle_above) * (
                below - above
            ) / (angle_below - angle_above)

            x = _linear_constrained_tes_opt(
                l,
                target_mean,
                (1 - alpha) * Q + alpha * Qin,
                max_el_current,
                max_total_current,
                extra_ineq=extra_ineq,
                extra_eq=extra_eq,
                log_level=log_level - 10,
            )
            angle = _calc_angle(x, Qin, l)
            logger.log(
                log_level,
                f"{it} alpha: {alpha:.2f}, angle: {angle:.2e}, max_angle: {max_angle:.2e}",
            )

            if angle > max_angle:
                above = alpha
                angle_above = angle
            else:
                below = alpha
                angle_below = angle
                x_below = np.copy(x)

            it += 1
            if it > max_iter:
                return x_below

        return x


def _bb_bounds_tes_problem(state, max_l0, linear, quadratic, max_el_current, func):
    """Returns upper bound, lower bound a child states for a TES optimization problem

    This is meant ro be a general interface for the BB algorithm. It requireas a function
    "func" with the following call
        x, objective = func(linear, quadratic, extra_ineq)
    "linear" is a list of linear factors, quadratic is a list of quadratic factors
    """
    if len(state.active) > max_l0:
        return 1e20, 1e20, None, None

    n = linear[0].shape[1]

    # Create a list of active + unasigned electrodes
    ac = np.zeros(n, dtype=bool)
    if len(state.active) < max_l0:
        ac[state.active + state.unassigned] = True
    elif len(state.active) == max_l0:
        ac[state.active] = True

    ##  Lower bound calculation

    # Create the extra l1 constraint (PS1.6)
    if len(state.active) > 0 and len(state.active) < max_l0:
        v = np.ones(n) / max_el_current
        v[state.active] = 0
        v = np.tile(v[ac], 2)
        extra_ineq = (v, max_l0 - len(state.active))
    else:
        extra_ineq = None

    # Run problem in reduced system
    linear_ac = [l[:, ac] for l in linear]
    quadratic_ac = [Q[np.ix_(ac, ac)] for Q in quadratic]
    x_ac, objective_lb = func(linear_ac, quadratic_ac, extra_ineq=extra_ineq)
    x_lb = np.zeros(n)
    x_lb[ac] = x_ac
    ## Upper bound calculation
    # Solve PS2
    x_ac, _ = func(linear_ac, quadratic_ac)
    x_ub1 = np.zeros(n)
    x_ub1[ac] = x_ac

    # Solve problem PS3
    # Select the "l0 - active" largest unassigned electrodes
    order_unasigned = np.argsort(-np.abs(x_ub1[state.unassigned]))
    selected_unasigned = [
        state.unassigned[i] for i in order_unasigned[: max_l0 - len(state.active)]
    ]

    # Select the active electrodes plus the largest unassigned electrodes
    s = state.active + selected_unasigned
    linear_s = [l[:, s] for l in linear]
    quadratic_s = [Q[np.ix_(s, s)] for Q in quadratic]
    x_s, objective_ub = func(linear_s, quadratic_s)
    x_ub = np.zeros(n)
    x_ub[s] = x_s

    # Split by activating / deactivating the unassigned electrode with the most current
    split_var = state.unassigned[np.argmax(np.abs(x_ub[state.unassigned]))]
    child1 = state.activate(split_var)
    child2 = state.inactivate(split_var)
    state.x_ub = x_ub
    state.x_lb = x_lb

    return objective_ub, objective_lb, child1, child2


def _constrained_eigenvalue(Q):
    """
    Finds the stationary values for the constrained eigenvalue problem
    Qx = lambda x
    such that 1^Tx = 0
    Q is real and symmetric

    Golub, Gene H. "Some modified matrix eigenvalue problems." Siam Review 15.2 (1973) 318-334.
    https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.454.9868&rep=rep1&type=pdf
    """
    n = Q.shape[0]
    c = 1 / np.sqrt(n) * np.ones(n)
    P = np.eye(n) - np.outer(c, c)
    K = P.dot(Q).dot(P)
    eigvals, z = np.linalg.eigh(K)
    eigvec = P.dot(z)
    return eigvals, eigvec


def _norm_constrained_tes_opt(
    Qnorm,
    target_norm,
    Q,
    max_el_current,
    max_total_current,
    extra_ineq=None,
    extra_eq=None,
    log_level=20,
    n_start=20,
    eigval_cutoff=1e-6,
):
    """Convex-concave algorithm to solve the problem
    minimize   x^T Q x
    subject to x^T Q_{i} x = t_i^2,  i = 1, 2, ...
               |x| <= I_ind
               |x|_1 <= 2I_tot

    Based on algorithm 3.1 fro Lipp and Boyd, Optim. Eng. (2016)n{align}
    """
    if Qnorm.ndim == 2:
        Qnorm = Qnorm[None, ...]
    assert eigval_cutoff > 0 and eigval_cutoff < 1
    # Statring positions given by constrained eigenvalue problem
    eigval, eigvec = _constrained_eigenvalue(np.sum(Qnorm, axis=0))
    # select n_start positions
    n_start = min(len(eigval), n_start)
    eigval = eigval[: -(n_start + 1) : -1]
    starting_x = eigvec[:, : -(n_start + 1) : -1]
    # Use the eigenvalues to cut-off small (ofter problematic) eigenvectors
    starting_x = starting_x[:, eigval > eigval_cutoff * np.max(eigval)]
    # Run optimizations
    max_norm = 0
    x_max_norm = None
    min_energy = np.inf
    x_min_energy = None
    for i, x0 in enumerate(starting_x.T):
        x, energy, norm = _norm_opt_x0(
            x0,
            Qnorm,
            target_norm,
            Q,
            max_el_current,
            max_total_current,
            extra_ineq=extra_ineq,
            extra_eq=extra_eq,
            log_level=log_level - 10,
        )
        if np.sum(norm) > np.sum(max_norm):
            max_norm = np.sum(norm)
            x_max_norm = x
        if np.all(norm > target_norm * 0.99) and energy < min_energy:
            min_energy = energy
            x_min_energy = x

    return x_max_norm if x_min_energy is None else x_min_energy


def _norm_opt_x0(
    x0,
    Qnorm,
    target_norm,
    Q,
    max_el_current,
    max_total_current,
    extra_ineq=None,
    extra_eq=None,
    log_level=20,
):

    x = x0.copy()

    scale = max_total_current

    # not necessary to rescale Q and Qnorm
    x = x / scale
    target_norm = target_norm / scale
    max_el_current = max_el_current / scale
    max_total_current = max_total_current / scale

    n = len(x)
    n_eqs = Qnorm.shape[0]

    l1norm_x0 = np.linalg.norm(x, 1)
    if l1norm_x0 > max_total_current * 2:
        x *= (max_total_current * 2) / l1norm_x0

    linfnorm_x0 = np.linalg.norm(x, np.inf)
    if linfnorm_x0 > max_el_current:
        x *= max_el_current / linfnorm_x0

    lqnorm_x0 = np.sqrt(x.dot(Qnorm).dot(x))
    if np.any(lqnorm_x0 > target_norm):
        x *= np.min(target_norm / lqnorm_x0)

    x_init = x.copy()
    # First pass: Tries to maximize the quatratic
    # This is used just to evaluate feasibility

    tes_constraints = TESConstraints(n, max_total_current, max_el_current)

    C_, d_ = tes_constraints._l1_constraint()
    A_, b_ = tes_constraints._kirchhoff_constraint()

    # Extra constraints
    if extra_ineq is not None:
        C_ = np.vstack([C_, extra_ineq[0]])
        d_ = np.hstack([d_, extra_ineq[1] / scale])

    if extra_eq is not None:
        A_ = np.vstack([A_, extra_eq[0]])
        b_ = np.hstack([b_, extra_eq[1] / scale])

    # First of all, see if the norm can be reached
    n_iter = 0
    max_iter = 100
    cur_min = np.inf
    while n_iter < max_iter:
        l = x.dot(Qnorm)
        l_ = np.hstack([l, -l])
        sol = scipy.optimize.linprog(
            c=-np.average(l_, axis=0),
            A_ub=np.vstack([C_, l_]),
            b_ub=np.hstack([d_, target_norm**2]),
            # constraints due to Kirchhoff's law
            A_eq=A_,
            b_eq=b_,
            bounds=(0, max_el_current),
            method="highs",
            options=dict(
                # Sets the primal/dual feasibility tolerances for `highs-ds`;
                # `highs-ipm` uses the minimum of the two
                dual_feasibility_tolerance=1e-8,
                primal_feasibility_tolerance=1e-8,
            ),
        )
        if sol.status != 0:
            raise RuntimeError(sol.message)

        if (delta := cur_min - sol.fun) > 0:
            x = sol.x[:n] - sol.x[n:]
            if delta < 1e-3 * np.abs(sol.fun):
                break
            else:
                cur_min = sol.fun
                n_iter += 1
        else:
            break

    # If the target norm cannot be reached, exit
    if np.any(np.sqrt(x.dot(Qnorm).dot(x)) < target_norm * 0.9):
        logger.log(log_level, "Target norm could not be reached")
        x *= scale
        return x, x.dot(Q).dot(x), np.sqrt(x.dot(Qnorm).dot(x))

    # notice, I reset x on purpose
    x = x_init
    # Slack term penalty
    tau = 1e-2 * x.dot(Q).dot(x) / np.sum(x.dot(Qnorm).dot(x))
    # How much tau gets multiplied per iteration
    mu = 2
    # Maximum value
    max_tau = 1e4 * tau
    # Stuff to observe the iterations
    prev_energy = np.inf
    prev_norms = 0
    n_iter = 0

    # Preparations for the QPs
    Q_ = np.block([[Q, -Q], [-Q, Q]])

    # Bound constraints
    C_b, d_b = tes_constraints._bound_contraints()
    C_ = np.vstack([C_, C_b])
    d_ = np.hstack([d_, d_b])

    # Add slack variables
    x_ = np.hstack([x, -x, np.zeros(n_eqs)])
    x_[x_ < 0] = 0
    # Linear portion of the objective
    a_ = np.zeros(2 * n + n_eqs)
    # Quadratic part of the objective
    Q_ = np.hstack([Q_, np.zeros((Q_.shape[0], n_eqs))])
    Q_ = np.vstack([Q_, np.zeros((n_eqs, Q_.shape[1]))])
    # Inequality
    C_ = np.vstack([C_, np.zeros((2 * n_eqs, C_.shape[1]))])
    C_ = np.hstack([C_, np.zeros((C_.shape[0], n_eqs))])
    C_[-2 * n_eqs : -n_eqs, -n_eqs:] = -np.eye(n_eqs)
    C_[-n_eqs:, -n_eqs:] = -np.eye(n_eqs)
    d_ = np.hstack([d_, np.zeros(2 * n_eqs)])
    # Equality
    A_ = np.hstack([A_, np.zeros((A_.shape[0], n_eqs))])

    while n_iter < max_iter:
        l = -2 * x.dot(Qnorm)
        x_ = np.hstack([x, -x, np.zeros(n_eqs)])
        x_[x_ < 0] = 0
        squared_norm = x.dot(Qnorm).dot(x)
        # Add the slack variables
        x_[-n_eqs:] = (target_norm**2) - squared_norm
        # Update the penalty
        a_[-n_eqs:] = tau
        # Inequalily Cx < d
        C_[-n_eqs:, : 2 * n] = np.hstack([l, -l])
        d_[-n_eqs:] = -squared_norm - (target_norm**2)
        # Run the QP
        # Sometimes due to numerical instabilities this can fail
        abs_x_max = np.abs(x).max()
        try:
            x_ = _active_set_QP(
                l=a_,
                Q=Q_,
                C=C_,
                d=d_,
                x0=x_,
                A=A_,
                b=b_,
                tol_primal=1e-3 * abs_x_max * scale,
                tol_feasibility_x0=1e-3 * abs_x_max,
                tol_zero_div=1e-9 / scale,
            )
        except ValueError:
            return None, np.inf, 0
        x = x_[:n] - x_[n : 2 * n]
        norms = np.sqrt(x.dot(Qnorm).dot(x))
        if np.any(norms > target_norm):
            x *= np.min(target_norm / norms)
        energy = x.dot(Q).dot(x)
        norms = np.sqrt(x.dot(Qnorm).dot(x))
        # 2 Stoping criteria:
        # Energy should stop decreasing
        # Norm should stop increasing
        norms_str = np.array2string(
            norms, formatter={"float_kind": lambda x: "%.2e" % x}
        )
        logger.log(log_level, f"{n_iter}: Energy: {energy: .2e} Norms: {norms_str}")
        if (prev_energy - energy) < 1e-3 * energy and np.all(
            norms - prev_norms < 1e-3 * norms
        ):
            break
        else:
            tau = min(mu * tau, max_tau)
            prev_energy = energy
            prev_norms = norms
            n_iter += 1

    return x * scale, energy * scale**2, norms * scale


def _least_squares_tes_opt(
    l,
    Q,
    max_el_current,
    max_total_current,
    extra_ineq=None,
    extra_eq=None,
    log_level=10,
):
    n = Q.shape[1]
    tes_constraints = TESConstraints(n, max_total_current, max_el_current)
    # First solve an LP to get a feasible starting point
    l_ = np.hstack([l, -l])

    C_, d_ = tes_constraints._l1_constraint()
    A_, b_ = tes_constraints._kirchhoff_constraint()

    # Inequality constraints
    if extra_ineq is not None:
        C_ = np.vstack([C_, extra_ineq[0]])
        d_ = np.hstack([d_, extra_ineq[1]])

    if extra_eq is not None:
        A_ = np.vstack([A_, extra_eq[0]])
        b_ = np.hstack([b_, extra_eq[1]])

    Q_ = np.block([[Q, -Q], [-Q, Q]])

    x = _eq_constrained_QP(np.squeeze(l), Q, A_[:, :n], b_)
    x_ = np.hstack([x, -x])
    x_[x_ < 0] = 0
    # I leave some gap just so that I don't start with too many
    # active contraints
    if np.linalg.norm(x_, 1) > 1.8 * max_total_current:
        x_ *= 1.8 * max_total_current / np.linalg.norm(x_, 1)
    if np.linalg.norm(x_, np.inf) > 0.9 * max_el_current:
        x_ *= 0.9 * max_el_current / np.linalg.norm(x_, np.inf)
    # Do the QP
    eps = 1e-3 * min(max_total_current, max_el_current, 1e-1)
    C_b, d_b = tes_constraints._bound_contraints()

    x_ = _active_set_QP(
        np.squeeze(l_),
        2 * Q_,
        np.vstack([C_b, C_]),
        np.hstack([d_b, d_]),
        x_,
        A_,
        b_,
        tol_primal=eps,
        tol_feasibility_x0=eps,
    )

    x = x_[:n] - x_[n:]
    return x


def _active_set_QP(
        l: npt.NDArray,
        Q: npt.NDArray,
        C: npt.NDArray,
        d: npt.NDArray,
        x0: npt.NDArray,
        A: None | npt.NDArray =None,
        b: None | npt.NDArray = None,
        tol_primal: float = 1e-5,
        tol_feasibility_x0: float = 1e-5,
        tol_zero_div: float = 1e-9,
    ):
    """Solves the problem

        minimize    l^T x + 1/2 x^T Q x
        subject to  Cx <= d
                    Ax = b

    Parameters
    ----------
    l :
    Q :
    C :
        Left-hand side matrix for inequality constraints.
    d :
        Right-hand side for inequality constraints.
    x0 :
        Initial point.
    A :
        Left-hand side matrix for equality constraints.
    b :
        Right-hand side for equality constraints.
    tol_primal : float
        Tolerance for primal variables, p. If all |p| are smaller than this,
        the solution is not updated in this iteration and the duals are checked
        for termination. (The duals use -`tol_primal` as tolerance.)
    tol_feasibility_x0 : float
        Tolerance for determining the feasibility of the initial point, `x0`
        (i.e., if the constraints are violated).
    tol_zero_div : float
        Tolerance for avoiding zero division. In particular, the values for
        which C @ p < tol_zero_div are ignored in (d - C @ x) / (C @ p).

    Returns
    -------
    x :
        The solution vector.

    References
    ----------
    Gill and Murray (1978). Numerically stable methods for quadratic
        programming. https://link.springer.com/article/10.1007/BF01588976
    """
    # get the active set:
    x = np.copy(x0)
    n = len(x0)
    n_iter = 0
    active = np.abs(C.dot(x) - d) < tol_feasibility_x0

    # Check feasibility of x0
    assert np.all(C.dot(x) <= d + tol_feasibility_x0), f"Infeasible start ({C.dot(x)} > {d})"
    if A is not None and b is not None:
        assert np.allclose(A.dot(x), b, atol=tol_feasibility_x0), f"Infeasible start ({A.dot(x)} != {b})"

    # if not np.allclose(A.dot(x0), b, atol=eps):
    #     print(f"Infeasible start ({A.dot(x0)} != {b})")

    if A is None:
        A = np.empty((0, len(x0)))
    else:
        A = np.atleast_2d(A)
        if A.shape[0] > A.shape[1]:
            A = A.T

    max_iter = max(200, 2 * len(x0))
    added_iq_constraint = -1
    duals_inequality = np.zeros(C.shape[0])
    Active = np.vstack([A, C[active, :]])
    n_active = Active.shape[0]
    if n_active > 0:
        ZY, R = np.linalg.qr(Active.T, "complete")

    else:
        ZY = np.eye(n)
        R = np.zeros((n, n))

    Y = ZY[:, :n_active]
    Z = ZY[:, n_active:]
    while n_iter <= max_iter:
        l_i = l + Q.dot(x)
        Y = ZY[:, :n_active]
        Z = ZY[:, n_active:]
        if n_active >= n:
            p = np.zeros(n)
        else:
            try:
                w_z = np.linalg.solve(Z.T.dot(Q).dot(Z), -Z.T.dot(l_i))
            except np.linalg.LinAlgError:
                w_z = np.linalg.lstsq(Z.T.dot(Q).dot(Z), -Z.T.dot(l_i), rcond=None)[0]
            p = Z.dot(w_z)
            p = np.squeeze(p.T)
        # If no update, check the duals and try to terminate
        if np.all(np.abs(p) < tol_primal):
            # calculate duals
            if Active.shape[0] > 0:
                try:
                    duals = scipy.linalg.solve_triangular(
                        R[:n_active, :], -Y.T.dot(l_i)
                    )
                except (scipy.linalg.LinAlgError, ValueError, TypeError):
                    duals = np.linalg.lstsq(Active.T, -l_i, rcond=None)[0]

                duals = np.atleast_1d(np.squeeze(duals.T))
                duals_inequality = duals[A.shape[0] :]

            if np.all(duals_inequality > -tol_primal):
                break

            min_C_dual = np.argmin(duals_inequality) + A.shape[0]
            # anti-cycling
            if n_iter > 0.5 * max_iter and min_C_dual == Active.shape[0] - 1:
                break

            Active = np.delete(Active, min_C_dual, axis=0)
            ZY, R = scipy.linalg.qr_delete(ZY, R, min_C_dual, which="col")
            n_active -= 1

        else:
            den = C.dot(p)
            s = den > tol_zero_div
            alpha_violation = (d[s] - C[s, :].dot(x)) / den[s]
            # if no constraint can be violated, step size is one
            if len(alpha_violation) == 0 or np.min(alpha_violation >= 1):
                alpha = 1.0

            # if one constraint can be violated, calculate the maximum step size
            else:
                alpha = np.min(alpha_violation)
                went_to_bounds = np.argmin(alpha_violation)
                indices = np.where(s)[0]
                added_iq_constraint = indices[went_to_bounds]
                Active = np.vstack([Active, C[added_iq_constraint, :]])
                # Update the QR decomposition
                if n_active > 0:
                    ZY, R = scipy.linalg.qr_insert(
                        ZY, R, C[added_iq_constraint, :], n_active, which="col"
                    )
                    n_active += 1
                else:
                    ZY, R = np.linalg.qr(
                        np.atleast_2d(C[added_iq_constraint, :]).T, "complete"
                    )
                    n_active += 1
            x = x + alpha * p

        n_iter += 1
    if n_iter >= max_iter:
        # raise ValueError('Maximal number of iterations reached!')
        pass

    return x


def _eq_constrained_QP(l, Q, A, b):
    Q_qr, _ = np.linalg.qr(A.T, "complete")
    m, n = A.shape
    Y = Q_qr[:, np.arange(m)]
    Z = Q_qr[:, np.arange(m, n)]
    w_y = np.linalg.solve(A.dot(Y), b)
    p = Y.dot(w_y)
    w_z = np.linalg.solve(Z.T.dot(Q).dot(Z), -Z.T.dot(l + Q.dot(p)))
    x = p + Z.dot(w_z)
    x = np.squeeze(x.T)
    return x


class bb_state(object):
    """State for branch and bound algorithm. Contains a list of active, inactive and
    unassigned electrodes"""

    def __init__(self, active, inactive, unassigned):
        self.active = active
        self.inactive = inactive
        self.unassigned = unassigned
        self.x_lb = None
        self.x_ub = None

    def inactivate(self, i):
        if i not in self.unassigned:
            raise ValueError("Can only inactivate unassigned element")

        active = copy.copy(self.active)
        inactive = copy.copy(self.inactive)
        unassigned = copy.copy(self.unassigned)

        unassigned.remove(i)
        inactive.append(i)
        return bb_state(active, inactive, unassigned)

    def activate(self, i):
        if i not in self.unassigned:
            raise ValueError("Can only activate unassigned element")

        active = copy.copy(self.active)
        inactive = copy.copy(self.inactive)
        unassigned = copy.copy(self.unassigned)

        unassigned.remove(i)
        active.append(i)
        return bb_state(active, inactive, unassigned)


class bb_node(object):
    """Node for branch and bound algorithm.
    Contains the current state
    bounds_funct is a funtiom wich takes in a state and return the upper bound, lower
    bound, children1 and children2"""

    def __init__(self, state, bounds_func):
        self.state = state
        self.bounds_func = bounds_func
        self.ub_val, self.lb_val, self.child1, self.child2 = self.bounds_func(
            self.state
        )

    def split(self):
        """Returns 2 child nodes"""
        return bb_node(self.child1, self.bounds_func), bb_node(
            self.child2, self.bounds_func
        )


def _branch_and_bound(init, function, eps, max_k, log_level=20):
    """Brach and Bound Algorithm
    Parameters:
    --------
    init: object
        initial state
    function: func
        Function which takes the state and returns an upper bound, lower bound and 2
        child states
    eps: float
        Tolerance between upper and lower bound
    max_k: int
        Maximum depth
    """
    active_nodes = [bb_node(init, function)]
    k = 0
    return_val = None
    while True:
        lb = np.array([n.lb_val for n in active_nodes])
        ub = np.array([n.ub_val for n in active_nodes])
        # Prune
        keep = lb <= ub.min()
        keep[ub.argmin()] = True
        lb = lb[keep]
        ub = ub[keep]
        active_nodes = [n for i, n in enumerate(active_nodes) if keep[i]]
        logger.log(
            log_level, f"{k} Upper Bound: {ub.min():.2e}, Lower Bound: {lb.min():.2e}"
        )
        if ub.min() - lb.min() <= eps * np.abs(lb.min()) or k >= max_k:
            if ub.min() - lb.min() <= eps * np.abs(lb.min()):
                logger.log(log_level, "Tolerance reached, returning")
            else:
                logger.log(log_level, "Maximum number of iterations reached, retunning")
            return_val = active_nodes[ub.argmin()].state
            break
        q = active_nodes.pop(lb.argmin())
        c1, c2 = q.split()
        active_nodes.append(c1)
        active_nodes.append(c2)
        k += 1

    return return_val
