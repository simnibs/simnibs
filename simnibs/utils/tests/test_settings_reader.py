import tempfile
from simnibs.utils import settings_reader

class TestReadIni:
    def test_read_ini(self):
        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(b'[section]\n')
        f.write(b'int=1\n')
        f.write(b'bool=false\n')
        f.write(b'float=1.2\n')
        f.write(b'string="sss"\n')
        f.write(b'list=[1,2, 3]\n')
        f.write(b'dict={"a": 1, "b": 2}\n')
        f.write(b'[section2]\n')
        f.write(b'int=2\n')
        f.write(b'bool=true\n')
        f.write(b'nested_list=[1, [2, "b"]]\n')

        f.close()
        out = settings_reader.read_ini(f.name)
        assert out['section']['int'] == 1
        assert isinstance(out['section']['int'], int)
        assert not out['section']['bool']
        assert isinstance(out['section']['bool'], bool)
        assert out['section']['float'] == 1.2
        assert isinstance(out['section']['float'], float)
        assert out['section']['string'] == "sss"
        assert isinstance(out['section']['string'], str)
        assert out['section']['list'] == [1, 2, 3]
        assert out['section']['dict'] == {"a": 1, "b": 2}
        assert out['section2']['int'] == 2
        assert out['section2']['nested_list'] == [1, [2, "b"]]
