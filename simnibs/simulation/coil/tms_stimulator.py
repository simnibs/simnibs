class TMSStimulator:
    def __init__(self, name: str, brand: str, max_di_dt: float):
        self.name: str = name
        self.brand: str = brand

        self.max_di_dt: float = max_di_dt
        self._waveform = None
