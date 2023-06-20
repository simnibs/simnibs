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
import os
import sys
import shutil
import tempfile
import argparse
import zipfile
import requests
import functools

from simnibs import file_finder


# from https://stackoverflow.com/questions/16694907/download-large-file-in-python-with-requests
def download_file(url, local_filename, timeout=None):
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raw.read = functools.partial(r.raw.read, decode_content=True)
        with open(local_filename, 'wb') as f:
            shutil.copyfileobj(r.raw, f)

def download_extra_coils(timeout=None):
    version = 'master'
    url = f'https://github.com/simnibs/simnibs-coils/archive/{version}.zip'
    with tempfile.NamedTemporaryFile('wb', delete=False) as tmpf:
        tmpname = tmpf.name
    download_file(url, tmpf.name, timeout)
    with zipfile.ZipFile(tmpname) as z:
        z.extractall(file_finder.coil_models)
    os.remove(tmpname)
    src = os.path.join(file_finder.coil_models, f'simnibs-coils-{version}')
    dst = os.path.join(file_finder.coil_models, f'Deng-coils-{version}')
    shutil.move(src, dst)
    
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


