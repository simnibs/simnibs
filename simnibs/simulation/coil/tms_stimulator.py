from typing import Optional
import numpy as np
import numpy.typing as npt


class TMSWaveform:
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


class TMSStimulator:
    def __init__(
        self,
        name: Optional[str],
        brand: Optional[str],
        max_di_dt: float,
        waveforms: Optional[list[TMSWaveform]],
    ):
        self.name = name
        self.brand = brand

        self.max_di_dt = max_di_dt
        self.waveforms = waveforms if waveforms is not None else []

    def to_tcd(self) -> dict:
        tcd_stimulator = {}
        if self.name is not None:
            tcd_stimulator["name"] = self.name
        if self.brand is not None:
            tcd_stimulator["brand"] = self.brand
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
        max_di_dt = tcd_element["maxdIdt"]
        waveforms = []
        for waveform in tcd_element.get("waveformList", []):
            waveforms.append(TMSWaveform.from_tcd(waveform))

        return cls(name, brand, max_di_dt, waveforms)
