'''
    Optional post-install procedures in SimNIBS
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2020 Guilherme B Saturnino

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

'''
import sys
import argparse
from .postinstall_simnibs import download_extra_coils

def main():
    parser = argparse.ArgumentParser(
        prog="download_coils",
        description="Download extra coil files"
    )
    parser.add_argument(
        "--timeout", default=None,
        help="Timeout in seconds"
    )
    args = parser.parse_args(sys.argv[1:])
    print('Downloading extra coil files')
    download_extra_coils(timeout=args.timeout)


if __name__ == '__main__':
    main()


