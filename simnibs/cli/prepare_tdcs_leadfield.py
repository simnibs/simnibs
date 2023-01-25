import argparse
from pathlib import Path
import sys

from simnibs.cli.utils import clargs

from simnibs.simulation.eeg import compute_tdcs_leadfield


def parse_args(argv):

    program = dict(
        description="""
            Convenience function to run a TDCS simulation..
            This will compute the electric field and the central gray matter
            surface. Electrodes are modeled as circular with diameter of 10 mm
            and a thickness of 4 mm. Default conductivities will be used.

            For more control, use the python function
                simnibs.simulation.eeg.compute_tdcs_leadfield_for_eeg
            or define the simulation from scratch using
                simnibs.simulation.sim_struct.TDCSLEADFIELD.
            """,
    )

    montage = dict(
        type=str, help="""Name of EEG montage file defining the electrode positions."""
    )

    # Optional

    output_dir = dict(
        type=str,
        default="fem_{subid}",
        help="""Directory in which to store the results of the simulation. If
            it does not exist, it will be created (default: %(default)s).""",
    )

    mesh_electrodes = dict(
        action="store_true",
        help="""Model electrodes as rings with a diameter of 10 mm and a
        thickness of 4 mm. Otherwise, electrodes are modeled as points
        (default: %(default)s)."""
    )
    pardiso = dict(
        action="store_true",
        help="""Use PARDISO as solver for FEM calculations. Otherwise, PETSc will be
        used.""",
    )

    parser = argparse.ArgumentParser(**program)

    clargs.subid.add_to(parser)

    parser.add_argument("montage", **montage)
    parser.add_argument("--mesh_electrodes", **mesh_electrodes)
    parser.add_argument("-o", "--output_dir", **output_dir)
    parser.add_argument("--pardiso", **pardiso)
    clargs.subsampling.add_to(parser)

    return parser.parse_args(argv[1:])

def main():
    args = parse_args(sys.argv)

    m2m_dir = clargs.resolve_subject_id_path(args.subid)
    subid = m2m_dir.stem.lstrip("m2m_")
    fem_dir = Path(args.output_dir.format(subid=subid)).resolve()
    montage = Path(args.montage)
    point_electrodes = not args.mesh_electrodes

    init_kwargs = {}
    if args.pardiso:
        init_kwargs["solver_options"] = "pardiso"

    compute_tdcs_leadfield(
        m2m_dir,
        fem_dir,
        montage,
        args.subsampling,
        point_electrodes,
        init_kwargs=init_kwargs
    )

if __name__ == "__main__":
    main()