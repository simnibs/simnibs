import sys
import argparse

from simnibs.msh.eeg_positions import *
from simnibs.msh import surface, mesh_io
from simnibs.utils.csv_reader import read_csv_positions
from simnibs import __version__

def parseArguments(argv):
    # argument parsing exception handling
    parser = argparse.ArgumentParser(prog="eeg_positions",
                                     description= "generates a list "
                                     "of 10/10 EEG electrode positions "
                                     " from the 4 fiducials by calculating "
                                     " distances in the scalp. Outputs a '.csv'"
                                     " file with electrode positions as well as"
                                     " a '.geo' fiele for visualization in GMSH")
    parser.add_argument("-m", metavar="head_model.msh", required=True,
                        help='Head mesh where the electrode positions are to be'
                        ' calculated')
    parser.add_argument("-o", help="output base file name", metavar="10_10_positions", required=True)
    parser.add_argument("-i", help="input csv file, must define all 4 fiducials",
                        metavar="fiducial_positions.csv")
    parser.add_argument('-Nz', help="Nasion coordinates x y z", nargs=3, type=float)
    parser.add_argument('-Iz', help="Inion coordinates x y z", nargs=3, type=float)
    parser.add_argument("-LPA", help="Left preauricular point coordinates x y z", nargs=3, type=float)
    parser.add_argument("-RPA", help="Right preauricular point coordinates x y z", nargs=3, type=float)
    parser.add_argument('--NE_cap', help='Wether to calculate positions for the '
                        'Neuroelectrics cap', action='store_true')
    parser.add_argument('--version', action='version', version=__version__)

    args = parser.parse_args()
    if args.i is None and any([args.Iz is None, args.Iz is None, args.LPA is None, args.RPA is None]):
        raise ValueError('either a .csv-file has to be passed or all four fiducials have to be defined')

    if args.i is not None:
        print("Reading fiducials from .csv-file")
        type_, coordinates, _, name, _, _ = read_csv_positions(args.i)

        p = [i for i in range(len(type_)) if type_[i] == 'Fiducial']

        if len(p) < 4:
            raise ValueError('.csv-file has to contain all four fiducials')
        name = [name[i] for i in p]
        coordinates = [list(coordinates[i]) for i in p]

        args.Iz = coordinates[name.index('Iz')]
        args.Nz = coordinates[name.index('Nz')]
        args.LPA = coordinates[name.index('LPA')]
        args.RPA = coordinates[name.index('RPA')]

    return args

def main():
    args = parseArguments(sys.argv)
    print("Reading Mesh")
    mesh = mesh_io.read_msh(args.m)
    print("Reference points:")
    print("Nz:", args.Nz)
    print("Iz:", args.Iz)
    print("LPA:", args.LPA)
    print("RPA:", args.RPA)
    surf = surface.Surface(mesh, [5, 1005])
    eeg = FindPositions(args.Nz, args.Iz, args.LPA, args.RPA, surf)
    eeg = FindPositions(args.Nz, args.Iz, args.LPA, args.RPA, surf, args.NE_cap)
    print("Printing positions to file", args.o + '.csv')
    eeg.print2csv(args.o + '.csv')
    eeg.print2geo(args.o + '.geo')

if __name__ == '__main__':
    main()
