from dataclasses import dataclass
from pathlib import Path


@dataclass()
class CommandLineArgument:
    flags: list  # [str]
    actions: dict

    def add_to(self, parser):
        parser.add_argument(*self.flags, **self.actions)


def resolve_subject_id_path(subid):
    m2m_dir = Path(subid).resolve()
    return (
        m2m_dir
        if m2m_dir.name.startswith("m2m_")
        else m2m_dir.with_name(f"m2m_{m2m_dir.name}")
    )


subid = CommandLineArgument(
    ["subid"],
    dict(
        type=str,
        help="""Subject ID or /path/to/{subid} or /path/to/m2m_{subid}.
        The former will resolve to "m2m_{subid}" in the current working
        directory. The latter cases both resolve to /path/to/m2m_{subid}.
        """,
    ),
)
subsampling = CommandLineArgument(
    ["-s", "--subsampling"],
    dict(
        type=int,
        default=10000,
        help="""Target number of sources in each hemisphere. Original
        surfaces are ~100,000 vertices (default: %(default)s).""",
    ),
)

fsaverage = CommandLineArgument(
    ["--fsaverage"],
    dict(
        type=int,
        default=5,
        choices=[5, 6, 7],
        help="""
        Which fsaverage subdivision to morph to. The number of vertices per
        hemisphere for each subdivision factor is as follows:
            fsaverage5 (10,242),
            fsaverage6 (40,962), and
            fsaverage7 (163,842).
        fsaverage7 is the default (full resolution) fsaverage template.
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
    ["trans"],
    dict(
        type=str,
        help="""
        Name of a MAT file containing the variable `trans` which is a 4x4
        affine transformation matrix mapping from the coordinate system of the
        electrodes to MRI subject space. This is applied to the positions
        before writing the result.
    """,
    ),
)

