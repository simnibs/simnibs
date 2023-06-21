from typing import Optional

import numpy as np
import numpy.typing as npt


class TmsWaveform:
    """A representation of a TMS waveform

    Parameters
    ----------
    name : Optional[str]
        The name of the waveform
    time : npt.NDArray[np.float_]
        The recording time stamps
    signal : npt.NDArray[np.float_]
        The recorded signal
    fit : Optional[npt.NDArray[np.float_]]
        A fitted version of the recorded signal

    Attributes
    ----------------------
    name : Optional[str]
        The name of the waveform
    time : npt.NDArray[np.float_]
        The recording time stamps
    signal : npt.NDArray[np.float_]
        The recorded signal
    fit : Optional[npt.NDArray[np.float_]]
        A fitted version of the recorded signal
    """

    def __init__(
        self,
        name: Optional[str],
        time: npt.NDArray[np.float_],
        signal: npt.NDArray[np.float_],
        fit: Optional[npt.NDArray[np.float_]],
    ):
        self.name = name
        self.time = time
        self.signal = signal
        self.fit = fit

    def to_tcd(self) -> dict:
        tcd_waveform = {}
        if self.name is not None:
            tcd_waveform["name"] = self.name
        tcd_waveform["time"] = self.time.tolist()
        tcd_waveform["signal"] = self.signal.tolist()
        if self.fit is not None and len(self.fit) > 0:
            tcd_waveform["fit"] = self.fit.tolist()
        return tcd_waveform

    @classmethod
    def from_tcd(cls, tcd_element: dict):
        name = tcd_element.get("name")
        time = np.array(tcd_element["time"])
        signal = np.array(tcd_element["signal"])
        fit = None if tcd_element.get("fit") is None else np.array(tcd_element["fit"])
        return cls(name, time, signal, fit)


class TmsStimulator:
    """A representation of a TMS stimulator

    Parameters
    ----------
    name : Optional[str]
        The name of the stimulator
    brand : Optional[str]
        the brand of the stimulator
    max_di_dt :  Optional[float]
        Maximum dI/dt values for the stimulator
    waveforms : Optional[list[TmsWaveform]]
        A list of waveforms that can be generated with this stimulator

    Attributes
    ----------------------
    name : Optional[str]
        The name of the stimulator
    brand : Optional[str]
        the brand of the stimulator
    max_di_dt :  Optional[float]
        Maximum dI/dt value in A/s for the stimulator
    waveforms : list[TmsWaveform]
        A list of waveforms that can be generated with this stimulator
    di_dt : float
        The current dI/dt setting in A/s of this stimulator, used for the A field calculation of connected coil elements
    """
    def __init__(
        self,
        name: Optional[str],
        brand: Optional[str],
        max_di_dt: Optional[float],
        waveforms: Optional[list[TmsWaveform]],
    ):
        self.name = name
        self.brand = brand

        self.max_di_dt = max_di_dt
        self.waveforms = waveforms if waveforms is not None else []

        if self.max_di_dt is not None and 1.0 > self.max_di_dt:
            self._di_dt = self.max_di_dt
        else:
            self._di_dt = 1.0

    @property
    def di_dt(self) -> float:
        return self._di_dt

    @di_dt.setter
    def di_dt(self, value):
        if self.max_di_dt is not None and (self.di_dt < -self.max_di_dt or self.di_dt > self.max_di_dt):
            raise ValueError(
                f"dIdt must be within the range ({-self.max_di_dt}, {self.max_di_dt})"
            )
        else:
            self._di_dt = value

    def to_tcd(self) -> dict:
        tcd_stimulator = {}
        if self.name is not None:
            tcd_stimulator["name"] = self.name
        if self.brand is not None:
            tcd_stimulator["brand"] = self.brand
        if self.max_di_dt is not None:
            tcd_stimulator["maxdIdt"] = self.max_di_dt
        if self.waveforms is not None and len(self.waveforms) > 0:
            tcd_waveforms = []
            for waveform in self.waveforms:
                tcd_waveforms.append(waveform.to_tcd())
            tcd_stimulator["waveformList"] = tcd_waveforms
        return tcd_stimulator

    @classmethod
    def from_tcd(cls, tcd_element: dict):
        name = tcd_element.get("name")
        brand = tcd_element.get("brand")
        max_di_dt = tcd_element.get("maxdIdt")
        waveforms = []
        for waveform in tcd_element.get("waveformList", []):
            waveforms.append(TmsWaveform.from_tcd(waveform))

        return cls(name, brand, max_di_dt, waveforms)
