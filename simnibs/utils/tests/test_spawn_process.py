import sys
import subprocess

import pytest

from .. import spawn_process


def test_spawn_process_0():
    ret = spawn_process.spawn_process([sys.executable,'-h'])
    assert ret == 0

def test_spawn_process_1():
     with pytest.raises(subprocess.CalledProcessError):
         ret = spawn_process.spawn_process([sys.executable,'yeeeeet'])