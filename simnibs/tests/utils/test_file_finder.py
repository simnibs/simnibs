import os
from simnibs.utils import file_finder


class TestTemplates:
    def test_find(self):
        templates = file_finder.templates
        for k, fn in templates.__dict__.items():
            if fn is None:
                continue
            else:
                print(fn)
                is_file = os.path.isfile(fn)
                is_dir = os.path.isdir(fn)
                assert is_file or is_dir

class TestSubjectFiles:
    def test_define_fnamehead(self):
        s = file_finder.SubjectFiles(
            os.path.join('path', 'to', 'sub.msh'))
        assert s.fnamehead == os.path.join('path', 'to', 'sub.msh')
        assert s.subid == 'sub'
        assert s.basedir == os.path.join('path', 'to')
        assert s.subpath == os.path.join('path', 'to', 'm2m_sub')

    def test_define_subpath(self):
        s = file_finder.SubjectFiles(subpath=os.path.join('path', 'to', 'm2m_sub'))
        assert s.subpath == os.path.join('path', 'to', 'm2m_sub')
        assert s.subid == 'sub'
        assert s.basedir == os.path.join('path', 'to')
        assert s.fnamehead == os.path.join('path', 'to', 'sub.msh')

    def test_define_fnamehead_subpath(self):
        s = file_finder.SubjectFiles(
            os.path.join('some', 'random', 'file.msh'),
            subpath=os.path.join('path', 'to', 'm2m_sub'))

        assert s.fnamehead == os.path.join('some', 'random', 'file.msh')
        assert s.subid == 'sub'
        assert s.basedir == os.path.join('path', 'to')
        assert s.subpath == os.path.join('path', 'to', 'm2m_sub')

    def test_tensor_file(self):
        s = file_finder.SubjectFiles(os.path.join('path', 'to', 'sub.msh'))
        assert s.tensor_file == os.path.join('path', 'to', 'd2c_sub',
                                             'dti_results_T1space',
                                             'DTI_conf_tensor.nii.gz')

    def test_cap_files(self):
        s = file_finder.SubjectFiles(os.path.join('path', 'to', 'sub.msh'))
        assert s.eeg_cap_folder == os.path.join('path', 'to', 'm2m_sub',
                                                'eeg_positions')
        assert s.eeg_cap_1010 == os.path.join('path', 'to', 'm2m_sub',
                                              'eeg_positions',
                                              'EEG10-10_UI_Jurak_2007.csv')

        assert s.get_eeg_cap('test.csv') == os.path.join(
            'path', 'to', 'm2m_sub', 'eeg_positions', 'test.csv')

