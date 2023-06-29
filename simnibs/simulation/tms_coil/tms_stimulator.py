import base64
from typing import Optional

import numpy as np
import numpy.typing as npt


class TmsWaveform:
    """A representation of a TMS waveform

    Parameters
    ----------
    time : npt.ArrayLike
        The recording time stamps
    signal : npt.ArrayLike
        The recorded signal
    name : Optional[str], optional
        The name of the waveform, by default None
    fit : Optional[npt.ArrayLike], optional
        A fitted version of the recorded signal, by default None

    Attributes
    ----------------------
    name : Optional[str]
        The name of the waveform
    time : npt.NDArray[np.float64]
        The recording time stamps
    signal : npt.NDArray[np.float64]
        The recorded signal
    fit : Optional[npt.NDArray[np.float64]]
        A fitted version of the recorded signal
    """

    def __init__(
        self,
        time: npt.ArrayLike,
        signal: npt.ArrayLike,
        name: Optional[str] = None,
        fit: Optional[npt.ArrayLike] = None,
    ):
        self.name = name
        self.time = np.array(time, dtype=np.float64)
        self.signal = np.array(signal, dtype=np.float64)
        self.fit = None if fit is None else np.array(fit, dtype=np.float64)

        if self.time.ndim != 1:
            raise ValueError(
                f"Expected 'time' to have the shape (N) but shape was {self.time.shape}"
            )

        if self.signal.ndim != 1:
            raise ValueError(
                f"Expected 'signal' to have the shape (N) but shape was {self.signal.shape}"
            )

        if self.fit is not None and self.fit.ndim != 1:
            raise ValueError(
                f"Expected 'fit' to have the shape (N) but shape was {self.fit.shape}"
            )

    def to_tcd(self, ascii_mode: bool = False) -> dict:
        tcd_waveform = {}
        if self.name is not None:
            tcd_waveform["name"] = self.name

        if ascii_mode:
            tcd_waveform["time"] = self.time.tolist()
            tcd_waveform["signal"] = self.signal.tolist()
            if self.fit is not None and len(self.fit) > 0:
                tcd_waveform["fit"] = self.fit.tolist()
        else:
            tcd_waveform["time"] = base64.b64encode(self.time.tobytes()).decode("ascii")
            tcd_waveform["signal"] = base64.b64encode(self.signal.tobytes()).decode(
                "ascii"
            )
            if self.fit is not None and len(self.fit) > 0:
                tcd_waveform["fit"] = base64.b64encode(self.fit.tobytes()).decode(
                    "ascii"
                )

        return tcd_waveform

    @classmethod
    def from_tcd(cls, tcd_element: dict):
        name = tcd_element.get("name")

        if isinstance(tcd_element["time"], str):
            time = np.frombuffer(
                base64.b64decode(tcd_element["time"]), dtype=np.float64
            )
        else:
            time = np.array(tcd_element["time"])

        if isinstance(tcd_element["signal"], str):
            signal = np.frombuffer(
                base64.b64decode(tcd_element["signal"]), dtype=np.float64
            )
        else:
            signal = np.array(tcd_element["signal"])

        if tcd_element.get("fit") is None:
            fit = None
        elif isinstance(tcd_element["fit"], str):
            fit = np.frombuffer(base64.b64decode(tcd_element["fit"]), dtype=np.float64)
        else:
            fit = np.array(tcd_element["fit"])

        return cls(time, signal, name, fit)


class TmsStimulator:
    """A representation of a TMS stimulator

    Parameters
    ----------
    name : Optional[str]
        The name of the stimulator
    brand : Optional[str], optional
        the brand of the stimulator, by default None
    max_di_dt :  Optional[float], optional
        Maximum dI/dt values for the stimulator, by default None
    waveforms : Optional[list[TmsWaveform]], optional
        A list of waveforms that can be generated with this stimulator, by default None

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
        brand: Optional[str] = None,
        max_di_dt: Optional[float] = None,
        waveforms: Optional[list[TmsWaveform]] = None,
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
        if self.max_di_dt is not None and (
            self.di_dt < -self.max_di_dt or self.di_dt > self.max_di_dt
        ):
            raise ValueError(
                f"dIdt must be within the range ({-self.max_di_dt}, {self.max_di_dt})"
            )
        else:
            self._di_dt = value

    def to_tcd(self, ascii_mode: bool = False) -> dict:
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
                tcd_waveforms.append(waveform.to_tcd(ascii_mode))
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
