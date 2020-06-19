''' This script is meant to fix the "shebang" on windows executables coming from installation of wheels
    pip uses the distlib package to create entry points for wheels. These entry poits have a hard-coded
    path to the python executables. Here, we replace those paths with new ones

    Guilherme Saturnino, 2020
'''

import re
import struct
import sys
import glob
import os
import argparse


def replace_pyzzer_entry_point_shebang(all_data, new_prefix):
    """Code adapted from pyzzer. This is meant to deal with entry point exe's
    created by distlib, which consist of a launcher, then a shebang, then a zip
    archive of the entry point code to run.  We need to change the shebang.

    all_data: executable in binary format
    new_prefix: path to where the python executables are located

    """
    # Copyright (c) 2013 Vinay Sajip.
    #
    # Permission is hereby granted, free of charge, to any person obtaining a copy
    # of this software and associated documentation files (the "Software"), to deal
    # in the Software without restriction, including without limitation the rights
    # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    # copies of the Software, and to permit persons to whom the Software is
    # furnished to do so, subject to the following conditions:
    #
    # The above copyright notice and this permission notice shall be included in
    # all copies or substantial portions of the Software.
    #
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    # THE SOFTWARE.

    #  GBS 2020: This was taken from the conda-pack repository and further modified

    # Find if it has these bits to ensure its a distlib-created binary
    # If not, just return
    pos = all_data.rfind(b'PK\x05\x06')
    if pos < 0:
        return all_data
    # TODO: can also use pip._vendor.distlib.scripts._DEFAULT_MANIFEST
    begin_pad_shebang = all_data.rfind(b'</assembly>', 0) + len(b'</assembly>')
    if begin_pad_shebang < 0:
        raise OSError('Could not find begining of shebang line')
    end_pad_shebang = all_data.rfind(b'\n\r\n', begin_pad_shebang)
    if end_pad_shebang < 0:
        raise OSError('Could not find end of shebang line')
    shebang_line = all_data[begin_pad_shebang:end_pad_shebang]
    # find the interpreter, usually python.exe or pythonw.exe
    if shebang_line.endswith(b'python.exe') or shebang_line.endswith(b'python.exe"'):
        interpreter = b'python.exe'
    elif shebang_line.endswith(b'pythonw.exe') or shebang_line.endswith(b'pythonw.exe"'):
        interpreter = b'pythonw.exe'
    else:
        raise OSError('Could not find python interpreter name!')
    # Create the new shebang
    if hasattr(new_prefix, 'encode'):
        new_prefix = new_prefix.encode('utf-8')
    # Add quotes for escaping spaces
    new_shebang_line = b'#!"' + os.path.join(new_prefix, interpreter) + b'"'
    if len(new_shebang_line) > len(shebang_line):
        raise OSError('Cant set new shebang. Path too long?')
    new_shebang_line = (len(shebang_line) - len(new_shebang_line)) * b'\x00' + new_shebang_line

    all_data = b"".join([
        all_data[:begin_pad_shebang],
        new_shebang_line,
        all_data[end_pad_shebang:]
    ])
    return all_data



def update_prefix(path, new_prefix):
    with open(path, 'rb+') as fh:
        original_data = fh.read()
        fh.seek(0)
        data = replace_pyzzer_entry_point_shebang(original_data, new_prefix)

        # If the before and after content is the same, skip writing
        if data != original_data:
            fh.write(data)
            fh.truncate()

def fix_all_scripts_win(scripts_folder, python_prefix):
    # Remove double quotes in case we get it from the command line
    scripts_folder = scripts_folder.replace('"', '')
    scripts_folder = os.path.abspath(scripts_folder)
    python_prefix = python_prefix.replace('"', '')
    python_prefix = os.path.normpath(os.path.abspath(python_prefix))
    for exe_file in glob.glob(os.path.join(scripts_folder, '*.exe')):
        try:
            update_prefix(exe_file, python_prefix)
        except OSError as e:
            print(e)
    return

def fix_all_scripts_posix(scripts_folder, python_prefix):
    scripts_folder = os.path.abspath(scripts_folder)
    python_prefix = os.path.normpath(os.path.abspath(python_prefix))
    for script in glob.glob(os.path.join(scripts_folder, '*')):
        try:
            with open(script, 'r+') as f:
                file_contents = f.read()
                # This will update the path to the python interpreter in the shebang, if any
                updated_file = re.sub(
                    r'(^#!)([^.]*)(/python\d*\n.*)$',
                    r'#!' + python_prefix + r'\3',
                    file_contents
                )
                if updated_file != file_contents:
                    f.write(updated_file)
        # If it's a binary, I will get unicode decode errors
        except UnicodeDecodeError:
            pass

def parse_arguments(argv):
    parser = argparse.ArgumentParser(
        prog="fix_entrypoints",
        description="Changes the path to the python interpreter in scripts (posix) or executables (windows)")
    parser.add_argument("scripts_folder",
                        help="Folder with the entrypoints")
    parser.add_argument("python_prefix",
                        help="Folder with the python executables")
    return parser.parse_args(argv)


if __name__ == '__main__':
    args = parse_arguments(sys.argv[1:])
    if sys.platform == 'win32':
        fix_all_scripts_win(args.scripts_folder,  args.python_prefix)
    else:
        fix_all_scripts_posix(args.scripts_folder, args.python_prefix)