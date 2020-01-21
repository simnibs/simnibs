import sys

from simnibs.utils.spawn_process import spawn_process


def test_spawn_process_0():
    ret = spawn_process(f'{sys.executable} -h')
    assert ret == 0

def test_spawn_process_1():
    ret = spawn_process(f'{sys.executable} yeeeeet')
    assert ret != 1
