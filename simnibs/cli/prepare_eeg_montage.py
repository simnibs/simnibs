import argparse
from pathlib import Path
import sys

from simnibs.cli.utils.helpers import add_argument
from simnibs.cli.utils import args_eeg, parsers_eeg


def parse_args(argv):

    program = dict(
        description = """
            Convert electrode information from MNE-Python or FieldTrip data
            structure to SimNIBS format. The resulting electrode positions
            should be in MRI space and is writting to a CSV file.
        """,
    )

    montage = dict(
        type = str,
        help = """Name of resulting EEG montage file. `csv` extension is
        forced."""
    )

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("montage", **montage)

    parser = argparse.ArgumentParser(**program)
    sp_eeg = parser.add_subparsers(**parsers_eeg.read_format_eeg.kwargs)

    # MNE
    sp_mne = parsers_eeg.read_format_mne.add_to(sp_eeg, parents=[parent_parser])
    add_argument(sp_mne, args_eeg.info_mne)
    add_argument(sp_mne, args_eeg.trans_mne)

    # FieldTrip
    sp_fieldtrip = parsers_eeg.read_format_fieldtrip.add_to(sp_eeg, parents=[parent_parser])
    add_argument(sp_fieldtrip, args_eeg.info_fieldtrip)
    add_argument(sp_fieldtrip, args_eeg.trans_fieldtrip)

    return parser.parse_args(argv[1:])

def main():
    args = parse_args(sys.argv)

    montage = Path(args.montage).with_suffix(".csv")

    if args.format == 'mne':
        from simnibs.eeg.utils_mne import prepare_montage
    elif args.format == 'fieldtrip':
        from simnibs.eeg.utils_fieldtrip import prepare_montage
    else:
        raise ValueError(f'{args.format} is not a valid format.')
    prepare_montage(montage, args.info, args.trans)

if __name__ == '__main__':
    main()