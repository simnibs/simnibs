import argparse
from pathlib import Path
import sys

from simnibs.cli.utils.helpers import add_argument, resolve_subject_id_path
from simnibs.cli.utils import args_eeg, parsers_eeg, args_general
from simnibs.eeg.forward import make_forward

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
    add_argument(parent_parser, args_general.subid)
    add_argument(parent_parser, args_eeg.leadfield)
    add_argument(parent_parser, args_eeg.fsaverage)
    parent_parser.add_argument("--no-average-ref", **no_average_ref)

    parser = argparse.ArgumentParser(**program, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sp_eeg = parser.add_subparsers(**parsers_eeg.write_format_eeg.kwargs)

    # MNE
    sp_mne = parsers_eeg.write_format_mne.add_to(sp_eeg, parents=[parent_parser])
    add_argument(sp_mne, args_eeg.info_mne)
    add_argument(sp_mne, args_eeg.trans_mne)

    # FieldTrip
    sp_fieldtrip = parsers_eeg.write_format_fieldtrip.add_to(sp_eeg, parents=[parent_parser])

    return parser.parse_args(argv[1:])

def main():
    args = parse_args(sys.argv)

    kwargs = dict(
        m2m_dir = resolve_subject_id_path(args.subid),
        fname_leadfield = Path(args.leadfield).resolve(),
        out_format = args.format,
        morph_to_fsaverage = args.fsaverage,
        apply_average_proj = not args.no_average_ref,
        write = True,
    )
    if args.format == 'mne':
        kwargs["info"] = args.info
        kwargs["trans"] = args.trans

    make_forward(**kwargs)

if __name__ == "__main__":
    main()
