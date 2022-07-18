import argparse
from pathlib import Path
import sys

from simnibs.cli.utils import clargs, clpars
from simnibs.simulation.eeg import prepare_for_inverse

def parse_args(argv):

    program = dict(
        description = """
            Prepare a TDCS leadfield for EEG inverse calculations.
            Additionally, write the source space definition and a "morpher"
            (i.e., a sparse matrix) mapping from subject space to fsaverage
            space.
            """,
    )

    no_average_ref = dict(
        action = "store_true",
        help = "Do not apply an average reference to the final gain matrix.",
    )

    parent_parser = argparse.ArgumentParser(add_help=False)
    clargs.subid.add_to(parent_parser)
    clargs.leadfield.add_to(parent_parser)
    clargs.fsaverage.add_to(parent_parser)
    parent_parser.add_argument("--no-average-ref", **no_average_ref)

    parser = argparse.ArgumentParser(**program, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sp_eeg = parser.add_subparsers(**clpars.write_format_eeg.kwargs)

    # MNE
    sp_mne = clpars.write_format_mne.add_to(sp_eeg, parents=[parent_parser])
    clargs.info_mne.add_to(sp_mne)
    clargs.trans_mne.add_to(sp_mne)

    # FieldTrip
    sp_fieldtrip = clpars.write_format_fieldtrip.add_to(sp_eeg, parents=[parent_parser])

    return parser.parse_args(argv[1:])

if __name__ == "__main__":
    args = parse_args(sys.argv)

    kwargs = dict(
        m2m_dir = clargs.resolve_subject_id_path(args.subid),
        fname_leadfield = Path(args.leadfield).resolve(),
        out_format = args.format,
        morph_to_fsaverage = args.fsaverage,
        apply_average_proj = not args.no_average_ref,
        write = True,
    )
    if args.format == 'mne':
        kwargs["info"] = args.info
        kwargs["trans"] = args.trans

    prepare_for_inverse(**kwargs)
