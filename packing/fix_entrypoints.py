''' This script is meant to fix the "shebang" on windows executables coming from installation of wheels
    pip uses the distlib package to create entry points for wheels. These entry poits have a hard-coded
    path to the python executables. Here, we replace those paths with new ones

    Guilherme Saturnino, 2020
'''

import re
import struct
import sys


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

    launcher = shebang = None
    # Find if it has these bits to ensure its a distlib-created binary
    pos = all_data.rfind(b'PK\x05\x06')
    if pos < 0:
        return all_data
    #find the end of the assembly XML entry
    launcher=all
    # TODO: can also use pip._vendor.distlib.scripts._DEFAULT_MANIFEST
    begin_pad_shebang = all_data.rfind(b'</assembly>', 0) + len(b'</assembly>')
    end_pad_shebang = all_data.rfind(b'\n\r\n', begin_pad_shebang)
    shebang_line = all_data[begin_pad_shebang:end_pad_shebang]
    # find the interpreter, usually python.exe or pythonw.exe
    m = re.search(b'(\\\)python.*exe$', shebang_line)
    if m is None:
        raise OSError('Could not find python interpreter from binary file')
    else:
        interpreter = m.group(0)
    if hasattr(new_prefix, 'encode'):
        new_prefix = new_prefix.encode('utf-8')
    # TODO: Add quotes for escaping spaces?
    new_shebang_line = b'#!' + new_prefix + interpreter
    if len(new_shebang_line) > len(shebang_line):
        raise OSError('Cant set new shebang. Path too long?')
    new_shebang_line = (len(shebang_line) - len(new_shebang_line)) * b'\x00' + new_shebang_line

    all_data = b"".join([
        all_data[:begin_pad_shebang],
        new_shebang_line,
        all_data[end_pad_shebang:]
    ])
    '''
    end_cdr = all_data[pos + 12:pos + 20]
    cdr_size, cdr_offset = struct.unpack('<LL', end_cdr)
    arc_pos = pos - cdr_size - cdr_offset
    data = all_data[arc_pos:]
    if arc_pos > 0:
        pos = all_data.rfind(b'#!', 0, arc_pos)
        if pos >= 0:
            shebang = all_data[pos:arc_pos]
            if pos > 0:
                launcher = all_data[:pos]

        if data and shebang and launcher:
            if hasattr(placeholder, 'encode'):
                placeholder = placeholder.encode('utf-8')
            if hasattr(new_prefix, 'encode'):
                new_prefix = new_prefix.encode('utf-8')
            shebang = shebang.replace(placeholder, new_prefix)
            all_data = b"".join([launcher, shebang, data])
    '''
    return all_data



def update_prefix(path, new_prefix, path_out):
    with open(path, 'rb+') as fh:
        original_data = fh.read()
        fh.seek(0)

        data = replace_pyzzer_entry_point_shebang(original_data, new_prefix)

        # If the before and after content is the same, skip writing
        #if data != original_data:
        #    fh.write(data)
        #    fh.truncate()
    with open(path_out, 'wb+') as fh:
        fh.write(data)
        fh.truncate()

if __name__ == '__main__':
    update_prefix('entry_point.exe', r'C:\Users\gbsa\simnibs_dev\packing\pack\simnibs_env', 'entry_point_fixed.exe')
