import argparse
from pathlib import Path
import sys

from simnibs.cli.utils.helpers import add_argument, resolve_subject_id_path
from simnibs.cli.utils import args_eeg, args_general

from simnibs.eeg.forward import compute_tdcs_leadfield


def parse_args(argv):

    program = dict(
        description="""
            Convenience function to run a TDCS leadfield simulation which forms
            the basis for the EEG forward solution computed by SimNIBS.

            This will compute the electric field and interpolate it to the
            central gray matter surface. By default, electrodes are modelled as
            points (projected onto the skin surface). Alternatively, electrodes
            can be modelled as circular with diameter of 10 mm and a thickness
            of 4 mm and meshed onto the head model. Default conductivities will
            be used.

            For more control, use the python function

                simnibs.eeg.forward.compute_tdcs_leadfield

            or define the simulation from scratch using

                simnibs.simulation.sim_struct.TDCSLEADFIELD.
            """,
    )

    parser = argparse.ArgumentParser(**program)

    add_argument(parser, args_general.subid)

    add_argument(parser, args_eeg.montage)
    add_argument(parser, args_eeg.mesh_electrodes)
    add_argument(parser, args_eeg.output_dir)
    add_argument(parser, args_eeg.use_pardiso)

    args_eeg.subsampling.add_to(parser)

    return parser.parse_args(argv[1:])

def main():
    args = parse_args(sys.argv)

    m2m_dir = resolve_subject_id_path(args.subid)
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