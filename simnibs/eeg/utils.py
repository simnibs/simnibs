from pathlib import Path
from typing import Iterable, Union

from nibabel.affines import apply_affine
import numpy as np

from simnibs.utils import csv_reader
from simnibs.utils.transformations import (
    _get_nearest_triangles_on_surface,
    _project_points_to_surface,
    coordinates_nonlinear,
)
from simnibs.utils.spatial_transform import fit_matched_points_analytical
from simnibs.utils.file_finder import templates


class Montage:
    def __init__(
        self,
        name: Union[None, Path, str] = None,
        ch_names: Union[None, list[str], tuple[str]] = None,
        ch_pos: Union[None, list, tuple, np.ndarray] = None,
        ch_types: Union[None, list[str], tuple[str], str] = None,
        landmarks: Union[None, dict] = None,
        headpoints: Union[None, list, tuple, np.ndarray] = None,
    ):
        """
        montage = Montage.from_filename("easycap_BC_TMS64_X21")
        # montage = Montage( ... )
        montage.add_landmarks()
        """
        self.name = name
        self.ch_names = [] if ch_names is None else ch_names
        self.ch_pos = np.array([]).reshape(0, 3) if ch_pos is None else np.array(ch_pos)
        self.n_channels = self.ch_pos.shape[0]
        assert self.n_channels == len(self.ch_names)

        if isinstance(ch_types, str):
            ch_types = [ch_types] * self.n_channels
        self.ch_types = [] if ch_types is None else ch_types

        self.landmarks = landmarks or {}
        if self.landmarks:
            self.landmarks = {k: np.asarray(v) for k, v in self.landmarks.items()}
        self.headpoints = (
            np.array([]).reshape(0, 3) if headpoints is None else np.array(headpoints)
        )
        self.n_headpoints = self.headpoints.shape[0]

    def add_landmarks(self, filename=None):
        """The name of the file to read fiducials from..."""
        filename = templates.fiducials or filename
        # fmt: off
        ch_types, ch_pos, _, ch_names, _, _= csv_reader.read_csv_positions(filename)
        # fmt: on
        assert all(i == "Fiducial" for i in ch_types)
        self.landmarks = dict(zip(ch_names, ch_pos))

    def get_landmark_names(self):
        return list(self.landmarks.keys())

    def get_landmark_pos(self, names: Union[Iterable, None] = None):
        names = names or self.get_landmark_names()
        return np.array([self.landmarks[n] for n in names]).reshape(-1, 3)

    def apply_trans(self, trans):
        """Apply affine transformation."""
        self.ch_pos = apply_affine(trans, self.ch_pos)
        self.headpoints = apply_affine(trans, self.headpoints)
        if self.landmarks:
            for k, v in self.landmarks.items():
                self.landmarks[k] = apply_affine(trans, v).squeeze()

    def apply_deform(self, deform):
        """Apply a nonlinear deformation field, e.g., MNI2Conform_nonl or
        Conform2MNI_nonl.

        `deform` is a 4D array where the last dimension encodes the coordinates
        corresponding to a particular voxel.
        """
        # deform affine : `from` coords -> `from` voxels
        # deform data   : `from` voxels -> `to` coords
        deform_coords = deform.get_fdata()
        self.ch_pos = coordinates_nonlinear(self.ch_pos, (deform_coords, deform.affine))
        self.headpoints = coordinates_nonlinear(
            self.headpoints, (deform_coords, deform.affine)
        )
        if self.landmarks:
            for k, v in self.landmarks.items():
                self.landmarks[k] = coordinates_nonlinear(
                    v, (deform_coords, deform.affine)
                )

    def project_to_surface(self, surf, subset=None):
        n = 2

        split_idx = np.cumsum((self.n_channels, self.n_headpoints))
        points = np.concatenate((self.ch_pos, self.headpoints, self.get_landmark_pos()))
        pttris = _get_nearest_triangles_on_surface(points, surf, n, subset)
        _, _, projs, _ = _project_points_to_surface(points, surf, pttris)

        self.ch_pos, self.headpoints, proj_landmarks = np.split(projs, split_idx)
        if self.landmarks:
            for k, v in zip(self.get_landmark_names(), proj_landmarks):
                self.landmarks[k] = v

    def fit_to(self, other, **kwargs):
        """Fit the landmarks in `self` to those in `other` in a least squares
        sense, i.e., the transformation matrix moves points from `self` to
        `other`.

        kwargs to pass to `fit_matched_points_analytical`

        """
        assert isinstance(other, Montage)
        assert self.landmarks and other.landmarks
        names = set(self.get_landmark_names()).intersection(
            set(other.get_landmark_names())
        )
        return fit_matched_points_analytical(
            self.get_landmark_pos(names), other.get_landmark_pos(names), **kwargs
        )

    def write(self, filename):
        ch_types = list(self.ch_types)
        ch_pos = self.ch_pos.tolist()
        ch_names = list(self.ch_names)
        if self.landmarks:
            n = len(self.landmarks)
            ch_types += ["Fiducial"] * n
            ch_pos += list(self.landmarks.values())
            ch_names += list(self.landmarks.keys())

        csv_reader.write_csv_positions(
            filename, ch_types, ch_pos, None, ch_names, None, None
        )

    @classmethod
    def from_filename(cls, filename: Union[Path, str]):
        """Construct montage from a filename.

        filename
            Name of one of the standard montages included in SimNIBS or full
            path to a CSV file defining an EEG montage in the same format as
            those included in simnibs.
        """
        filename = (
            templates.get_eeg_montage(filename)
            if not Path(filename).exists()
            else filename
        )
        # fmt: off
        ch_types, ch_pos, _, ch_names, _, _  = csv_reader.read_csv_positions(filename)
        # fmt: on
        lm = np.array(ch_types) == "Fiducial"
        if lm.any():
            landmarks = dict(zip(ch_names[lm], ch_pos[lm]))
            not_fids = ~lm
            ch_types = ch_types[not_fids]
            ch_pos = ch_pos[not_fids]
            ch_names = ch_names[not_fids]
        else:
            landmarks = None
        assert all(np.isin(ch_types, ("Electrode", "ReferenceElectrode")))
        return cls(filename, ch_names, ch_pos, ch_types, landmarks)
