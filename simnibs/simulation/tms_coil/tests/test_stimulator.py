import pytest

from simnibs.simulation.tms_coil.tms_stimulator import TmsStimulator, TmsWaveform


class TestInitStimulator:
    def test_small_max_di_dt(self):
        stimulator = TmsStimulator(name="Test", max_di_dt=0.5)
        assert stimulator.di_dt == 0.5


class TestdIdT:
    def test_set_di_dt(self):
        with pytest.warns(UserWarning):
            TmsStimulator(name="Test", max_di_dt=21.4).di_dt = 22

        with pytest.warns(UserWarning):
            TmsStimulator(name="Test", max_di_dt=10.4).di_dt = -10.5


class TestInitWaveform:
    def test_wrong_shape_time(self):
        with pytest.raises(ValueError):
            TmsWaveform([[1], [2]], [1, 2])

        with pytest.raises(ValueError):
            TmsWaveform([[1, 2, 3]], [1, 2])

    def test_wrong_size_time(self):
        with pytest.raises(ValueError):
            TmsWaveform([], [1, 2])

        with pytest.raises(ValueError):
            TmsWaveform([1], [1, 2])

    def test_wrong_shape_signal(self):
        with pytest.raises(ValueError):
            TmsWaveform([1, 2], [[1], [2]])

        with pytest.raises(ValueError):
            TmsWaveform([1, 2], [[1, 2, 3]])

    def test_wrong_size_signal(self):
        with pytest.raises(ValueError):
            TmsWaveform([1, 2], [])

        with pytest.raises(ValueError):
            TmsWaveform([1, 2], [1])

    def test_wrong_shape_fit(self):
        with pytest.raises(ValueError):
            TmsWaveform([1, 2], [1, 2], fit=[[1], [2]])

        with pytest.raises(ValueError):
            TmsWaveform([1, 2], [1, 2], fit=[[1, 2, 3]])

    def test_wrong_size_fit(self):
        with pytest.raises(ValueError):
            TmsWaveform([1, 2], [1, 2], fit=[])

        with pytest.raises(ValueError):
            TmsWaveform([1, 2], [1, 2], fit=[1])
