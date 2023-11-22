from .helpers import CommandLineArgument


subsampling = CommandLineArgument(
    ["-s", "--subsampling"],
    dict(
        type=int,
        help="""Target number of sources in each hemisphere. Original
        surfaces are ~100,000 vertices (default: %(default)s).""",
    ),
)

fsaverage = CommandLineArgument(
    ["--fsaverage"],
    dict(
        type=int,
        choices=[10, 40, 160],
        help="""
        The resolution of the fsaverage template to morph to. The number
        denotes the (approximate) number of vertices per hemisphere such that
            10  ->  10,242 (fsaverage 5)
            40  ->  40,962 (fsaverage 6)
            160 -> 163,842 (fsaverage 7, full resolution)
        By default, no interpolator is computed.
        """,
    ),
)
leadfield = CommandLineArgument(
    ["leadfield"], dict(type=str, help="Name of the leadfield to process (HDF5 file).")
)

#### MNE, FieldTrip

info_mne = CommandLineArgument(
    ["info"],
    dict(
        type=str,
        help="""
        Name of an MNE-Python file (raw, epochs, evoked, or info)
        containing information about the data (e.g., electrode names and
        positions).
    """,
    ),
)

trans_mne = CommandLineArgument(
    ["trans"],
    dict(
        type=str,
        help="""
        Name of an MNE-Python trans file containing an affine transformation
        mapping between EEG head coordinates and MRI coordinates.
    """,
    ),
)

info_fieldtrip = CommandLineArgument(
    ["info"],
    dict(
        type=str,
        help="""
        Name of a MAT file containing the variable `elec` which is a struct
        describing the electrodes as obtained from `ft_read_sens` or from the
        `elec` field of a FieldTrip data set.
    """,
    ),
)

trans_fieldtrip = CommandLineArgument(
    ["-t", "--trans"],
    dict(
        type=str,
        default=None,
        help="""
        Name of a MAT file containing the variable `trans` which is a 4x4
        affine transformation matrix mapping from the coordinate system of the
        electrodes to MRI subject space. This is applied to the positions
        before writing the result. Units are assumed to be mm.
    """,
    ),
)


# TDCS leadfield

montage = CommandLineArgument(
    ["montage"],
    dict(
        type=str, help="""Name of EEG montage file defining the electrode positions."""
    ),
)

# Optional

output_dir = CommandLineArgument(
    ["--output_dir", "-o"],
    dict(
        type=str,
        default="fem_{subid}",
        help="""Directory in which to store the results of the simulation. If
        it does not exist, it will be created (default: %(default)s).""",
    ),
)

mesh_electrodes = CommandLineArgument(
    ["--mesh-electrodes"],
    dict(
        action="store_true",
        help="""Model electrodes as rings with a diameter of 10 mm and a
    thickness of 4 mm. Otherwise, electrodes are modeled as points
    (default: %(default)s).""",
    ),
)

use_pardiso = CommandLineArgument(
    ["--pardiso"],
    dict(
        action="store_true",
        help="""Use PARDISO as solver for FEM calculations. Otherwise, PETSc will be
    used.""",
    ),
)
