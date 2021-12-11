''' Updates SimNIBS

    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2020  Guilherme B Saturnino

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
import subprocess
import argparse
from packaging import version
import requests

from simnibs import __version__

GH_RELEASES_URL = 'https://api.github.com/repos/simnibs/simnibs/releases'

def _get_release_data(url):
    ''' Get avaliable versions and release data'''
    response = requests.get(url)
    # Raise an exception if the API call fails.
    response.raise_for_status()
    data = response.json()
    release_data = {}
    for d in data:
        if d['tag_name'][0] == 'v' and not d['prerelease']:
            release_data[d['tag_name'][1:]] = d

    return release_data


def update_simnibs(requested_version='latest'):
    if __version__ == requested_version:
        print(f'SimNIBS already at the requested version ({requested_version})')
        return

    release_data = _get_release_data(GH_RELEASES_URL)
    current_version = version.parse(__version__)
    
    # Find the absolute latest version
    latest_version = version.parse('0.0.0')
    for v in release_data.keys():
        v_parsed = version.parse(v)
        if v_parsed > latest_version:
            latest_version = v_parsed

    if (latest_version.release[0] > current_version.release[0] or
        latest_version.release[1] > current_version.release[1]):
        print(
            f'A new SimNIBS feature version ({latest_version.public}) '
            f'has been released. Please go to www.simnibs.org to download it'
        )

    # If the user requested the latest version
    if requested_version == 'latest':
        requested_version = current_version
        # Go through all the releases
        for v in release_data.keys():
            v_parsed = version.parse(v)
            # Check if the first two components of the version are the same
            # and check if the version number is greater than the current
            if (v_parsed.release[:2] == current_version.release[:2]
                and v_parsed > requested_version):

                requested_version = v_parsed
    else:
        if requested_version not in release_data.keys():
            raise ValueError(f"Invalid SimNIBS version: {requested_version}")
        requested_version = version.parse(requested_version)

    # If already at the current version
    if requested_version == current_version:
        print('SimNIBS already at the requested version')
        return

    if requested_version.release[:2] != current_version.release[:2]:
        raise ValueError(
            "Can't change SimNIBS feature version."
            "\nPlease go to www.simnibs.org to install"
            " the requested version"
            )

    print(f'Updating SimNIBS to version {requested_version.public}')
    url = release_data[requested_version.public]['html_url']
    subprocess.run([
        sys.executable,
        '-m', 'pip',
        'install', '--upgrade',
        '-f', url, '--no-index',
        f'simnibs=={requested_version.public}'
        ],
        check=True
    )


def main():
    parser = argparse.ArgumentParser(
        prog="update_simnibs",
        description="Update the SimNIBS package")
    parser.add_argument("-t", "--target-version",
                        help="Target version. "
                        "By default the latest debug version",
                        default="latest")
    parser.add_argument('--version', action='version', version=__version__)

    args = parser.parse_args(sys.argv[1:])
    update_simnibs(args.target_version)


if __name__ == '__main__':
    main()