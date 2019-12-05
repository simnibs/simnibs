import tempfile
import numpy as np

from simnibs.utils import csv_reader

class TestCSV:
    def test_read_csv_generic(self):
        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(b'type,x,y,z,extra1,extra2\n')
        f.write(b'Generic,1,2,3,a1,a2\n')
        f.write(b'Generic,1.2,2.4,3.6,b1,b2\n')
        f.close()
        type_, coordinates, extra, name, extra_cols, header = csv_reader.read_csv_positions(f.name)
        assert type_ == ['Generic', 'Generic']
        assert np.allclose(coordinates, [[1, 2, 3], [1.2, 2.4, 3.6]])
        assert header == ['type','x', 'y', 'z', 'extra1', 'extra2']
        assert name == ['a1', 'b1']
        assert extra_cols == [['a2'], ['b2']]

        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(b'type,x,y,z,extra1,extra2\n')
        f.write(b'Generic,1,2,3\n')
        f.write(b'Generic,1.2,2.4,3.6\n')
        f.close()
        type_, coordinates, extra, name, extra_cols, header = csv_reader.read_csv_positions(f.name)
        assert type_ == ['Generic', 'Generic']
        assert np.allclose(coordinates, [[1, 2, 3], [1.2, 2.4, 3.6]])
        assert header == ['type','x', 'y', 'z', 'extra1', 'extra2']
        assert name == [None, None]
        assert extra_cols == [None, None]


    def test_read_csv_fiducial_electrode(self):
        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(b'Fiducial,1,2,3,fiducial,a2\n')
        f.write(b'Electrode,1.2,2.4,3.6,9,8,7,electrode,b2,b3\n')
        f.write(b'ReferenceElectrode,1.2,2.4,7.1\n')
        f.close()
        type_, coordinates, extra, name, extra_cols, header = csv_reader.read_csv_positions(f.name)
        assert type_ == ['Fiducial', 'Electrode', 'ReferenceElectrode']
        assert np.allclose(coordinates, [[1, 2, 3],
                                         [1.2, 2.4, 3.6],
                                         [1.2, 2.4, 7.1]])
        assert header == []
        assert name == ['fiducial', 'electrode', None]
        assert extra_cols == [['a2'], ['b2', 'b3'], None]
        assert extra[0] is None
        assert np.allclose(extra[1], [9, 8, 7])
        assert extra[2] is None

    def test_read_csv_coil(self):
        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(b'Fiducial,1,2,3,fiducial,a2\n')
        f.write(b'ReferenceElectrode,1.2,2.4,7.1\n')
        f.write(b'CoilPos,1,2,3,4,5,6,7,8,9,10,coil,comment\n')
        f.close()
        type_, coordinates, extra, name, extra_cols, header = csv_reader.read_csv_positions(f.name)
        assert type_ == ['Fiducial', 'ReferenceElectrode', 'CoilPos']
        assert np.allclose(coordinates, [[1, 2, 3],
                                         [1.2, 2.4, 7.1],
                                         [1, 2, 3]])
        assert header == []
        assert name == ['fiducial',  None, 'coil']
        assert extra_cols == [['a2'], None, ['comment']]
        assert extra[0] is None
        assert extra[1] is None
        assert np.allclose(extra[2], [4, 5, 6, 7, 8, 9, 10])

    def test_write_csv_generic(self):
        type_ = ['Generic', 'Generic']
        coordinates = np.array([[1, 2, 3], [1.2, 2.4, 3.6]])
        header = ['type', 'x', 'y', 'z', 'extra1', 'extra2']
        name = ['a1', 'b1']
        extra_cols = [['a2'], ['b2']]
        extra = [None, None]
        f = tempfile.NamedTemporaryFile(delete=False)
        csv_reader.write_csv_positions(f.name, type_, coordinates, extra, name, extra_cols, header)
        with open(f.name, 'r') as f:
            assert f.readline().strip() == 'type,x,y,z,extra1,extra2'
            assert f.readline().strip() == 'Generic,1.0,2.0,3.0,a1,a2'
            assert f.readline().strip() == 'Generic,1.2,2.4,3.6,b1,b2'

    def test_write_csv_electrode(self):
        type_ = ['Fiducial', 'Electrode', 'ReferenceElectrode']
        coordinates = np.array([[1, 2, 3],
                                [1.2, 2.4, 3.6],
                                [1.2, 2.4, 7.1]])
        header = []
        name = ['fiducial', 'electrode', None]
        extra_cols = [['a2'], ['b2', 'b3'], None]
        extra = [None, np.array([9, 8, 7], dtype=float), None]

        f = tempfile.NamedTemporaryFile(delete=False)
        csv_reader.write_csv_positions(f.name, type_, coordinates, extra, name, extra_cols, header)
        with open(f.name, 'r') as f:
            assert f.readline().strip() == 'Fiducial,1.0,2.0,3.0,fiducial,a2'
            assert f.readline().strip() == 'Electrode,1.2,2.4,3.6,9.0,8.0,7.0,electrode,b2,b3'
            assert f.readline().strip() == 'ReferenceElectrode,1.2,2.4,7.1'


    def test_write_csv_coil(self):
        type_ = ['Fiducial', 'ReferenceElectrode', 'CoilPos']
        coordinates = np.array([[1, 2, 3],
                                [1.2, 2.4, 7.1],
                                [1, 2, 3]])
        header = []
        name = ['fiducial',  None, 'coil']
        extra_cols = [['a2'], None, ['comment']]
        extra = [None, None, np.array([4., 5., 6., 7., 8., 9., 10], dtype=float)]

        f = tempfile.NamedTemporaryFile(delete=False)
        csv_reader.write_csv_positions(f.name, type_, coordinates, extra, name, extra_cols, header)
        with open(f.name, 'r') as f:
            assert f.readline().strip() == 'Fiducial,1.0,2.0,3.0,fiducial,a2'
            assert f.readline().strip() == 'ReferenceElectrode,1.2,2.4,7.1'
            assert f.readline().strip() == \
                'CoilPos,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,coil,comment'

