'''
    Simple wrapper script ro call run_simnibs
'''
import sys
import argparse

from simnibs import run_simnibs
from simnibs import __version__

def parseArguments(argv):
    parser = argparse.ArgumentParser(
        prog="simnibs",
        description="Prepare, run and postprocess SimNIBS problems")
    parser.add_argument("simnibs_file", help="Input .mat or file")
    parser.add_argument("--cpus", type=int,
                        help="Maximum number of CPUs to run simulations",
                        default=1)
    parser.add_argument('--version', action='version', version=__version__)

    return parser.parse_args(argv)


def main():
    args = parseArguments(sys.argv[1:])
    run_simnibs(args.simnibs_file, args.cpus)


if __name__ == '__main__':
    main()
