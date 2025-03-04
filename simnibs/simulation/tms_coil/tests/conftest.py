import json
import os
from typing import Any

import numpy as np
import pytest

from simnibs.mesh_tools.mesh_io import Elements, Msh, Nodes
from simnibs.simulation.tms_coil.tms_coil import TmsCoil
from simnibs.simulation.tms_coil.tms_coil_deformation import (
    TmsCoilDeformationRange,
    TmsCoilRotation,
    TmsCoilTranslation,
)
from simnibs.simulation.tms_coil.tms_coil_element import (
    DipoleElements,
    LineSegmentElements,
    SampledGridPointElements,
)
from simnibs.simulation.tms_coil.tms_coil_model import TmsCoilModel
from simnibs.simulation.tms_coil.tms_stimulator import TmsStimulator, TmsWaveform
from simnibs.utils import file_finder

from .... import SIMNIBSDIR

@pytest.fixture(scope="module")
def minimal_tcd_coil_dict() -> dict[str, Any]:
    coil_dict = {
        "coilElementList": [
            {
                "stimulator": 0,
                "type": 1,
                "points": [[1.0, 2.0, 3.0]],
                "values": [[4.0, 5.0, 6.0]],
            }
        ],
        "stimulatorList": [{"name": "SimNIBS-Stimulator"}],
    }
    return coil_dict


@pytest.fixture(scope="module")
def minimal_tcd_coil() -> TmsCoil:
    coil = TmsCoil(
        [
            DipoleElements(
                TmsStimulator("SimNIBS-Stimulator", None, None, None),
                np.array([[1, 2, 3]]),
                np.array([[4, 5, 6]]),
            )
        ]
    )
    return coil


@pytest.fixture(scope="module")
def medium_tcd_coil_dict() -> dict[str, Any]:
    coil_dict = {
        "coilElementList": [
            {
                "stimulator": 0,
                "deformations": [0, 1, 2, 3],
                "type": 1,
                "points": [[1.0, 2.0, 3.0]],
                "values": [[0.0, 1.0, 2.0]],
            },
            {
                "stimulator": 0,
                "elementCasing": 0,
                "type": 2,
                "points": [[1.0, 2.0, 3.0]],
                "values": [[0.0, 1.0, 2.0]],
            },
            {
                "stimulator": 0,
                "type": 3,
                "data": [[[[1.0, 2.0, 3.0]]]],
                "affine": [
                    [1.0, 0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0],
                ],
            },
        ],
        "stimulatorList": [
            {
                "maxdIdt": 1,
                "waveformList": [{"time": [0.0, 100.0], "signal": [0.0, 0.99]}],
            }
        ],
        "deformList": [
            {"deformRange": 0, "type": "x"},
            {"deformRange": 1, "type": "y"},
            {"deformRange": 2, "type": "z"},
            {
                "deformRange": 3,
                "type": "rot2p",
                "point1": [0.0, 0.0, 0.0],
                "point2": [0.0, 0.0, 1.0],
            },
        ],
        "deformRangeList": [
            {"initial": -75, "range": [-100.0, -50.0]},
            {"initial": 0, "range": [-50.0, 50.0]},
            {"initial": 75, "range": [50.0, 100.0]},
            {
                "initial": 180,
                "range": [0.0, 360.0],
            },
        ],
        "coilModels": [
            {
                "points": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
                "faces": [[0, 1, 2]],
            }
        ],
    }
    return coil_dict


@pytest.fixture(scope="module")
def medium_tcd_coil() -> TmsCoil:
    coil_casings = [
        TmsCoilModel(
            Msh(
                Nodes(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])),
                Elements(triangles=np.array([[0, 1, 2]]) + 1),
            ),
            None,
        ),
    ]

    deformations = [
        TmsCoilTranslation(TmsCoilDeformationRange(-75, (-100, -50)), 0),
        TmsCoilTranslation(TmsCoilDeformationRange(0, (-50, 50)), 1),
        TmsCoilTranslation(TmsCoilDeformationRange(75, (50, 100)), 2),
        TmsCoilRotation(
            TmsCoilDeformationRange(180, (0, 360)),
            np.array([0, 0, 0]),
            np.array([0, 0, 1]),
        ),
    ]
    stimulator = TmsStimulator(
        None,
        max_di_dt=1,
        waveforms=[TmsWaveform(np.array([0, 100]), np.array([0.0, 0.99]))],
    )
    coil = TmsCoil(
        [
            DipoleElements(
                stimulator,
                np.array([[1, 2, 3]]),
                np.array([[0, 1, 2]]),
                deformations=[
                    deformations[0],
                    deformations[1],
                    deformations[2],
                    deformations[3],
                ],
            ),
            LineSegmentElements(
                stimulator,
                np.array([[1, 2, 3]]),
                np.array([[0, 1, 2]]),
                casing=coil_casings[0],
            ),
            SampledGridPointElements(
                stimulator,
                np.array([[[[1, 2, 3]]]]),
                np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
            ),
        ],
    )
    return coil


@pytest.fixture(scope="module")
def full_tcd_coil_dict() -> dict[str, Any]:
    coil_dict = {
        "name": "SimNIBS-TMS-Coil",
        "brand": "SimNIBS",
        "version": "V1",
        "limits": [[0.0, 255.0], [0.0, 255.0], [0.0, 255.0]],
        "resolution": [1.0, 1.0, 1.0],
        "coilCasing": 0,
        "coilElementList": [
            {
                "name": "Dipole Elements",
                "stimulator": 0,
                "elementCasing": 1,
                "deformations": [0],
                "type": 1,
                "points": [[1.1, 2.1, 3.1]],
                "values": [[100.1, 101.1, 102.1]],
            },
            {
                "name": "Line Elements",
                "stimulator": 0,
                "elementCasing": 2,
                "deformations": [0, 1],
                "type": 2,
                "points": [[1.2, 2.2, 3.2]],
                "values": [[100.2, 101.2, 102.2]],
            },
            {
                "name": "Digitized Grid Elements",
                "stimulator": 0,
                "elementCasing": 3,
                "deformations": [2, 3],
                "type": 3,
                "data": [[[[1.4, 2.4, 3.4]]]],
                "affine": [
                    [1.0, 0.0, 0.0, 1.4],
                    [0.0, 1.0, 0.0, 1.4],
                    [0.0, 0.0, 1.0, 1.4],
                    [0.0, 0.0, 0.0, 1.0],
                ],
            },
        ],
        "stimulatorList": [
            {
                "name": "SimNIBS-Stimulator",
                "brand": "SimNIBS",
                "maxdIdt": 1,
                "waveformList": [
                    {
                        "name": "biphasic",
                        "time": [0.0, 100.0],
                        "signal": [0.0, 0.99],
                        "fit": [0.0, 1.0],
                    }
                ],
            }
        ],
        "deformList": [
            {"deformRange": 0, "type": "x"},
            {"deformRange": 1, "type": "y"},
            {"deformRange": 2, "type": "z"},
            {
                "deformRange": 3,
                "type": "rot2p",
                "point1": [0.0, 0.0, 0.0],
                "point2": [0.0, 0.0, 1.0],
            },
        ],
        "deformRangeList": [
            {"initial": -75, "range": [-100.0, -50.0]},
            {"initial": 0, "range": [-50.0, 50.0]},
            {"initial": 75, "range": [50.0, 100.0]},
            {"initial": 180, "range": [0.0, 360.0]},
        ],
        "coilModels": [
            {
                "points": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
                "faces": [[0, 1, 2]],
                "minDistancePoints": [[0.0, -1.0, -2.0]],
            },
            {
                "points": [[10.0, 11.0, 12.0], [13.0, 14.0, 15.0], [16.0, 17.0, 18.0]],
                "faces": [[0, 1, 2]],
                "minDistancePoints": [[-3.0, -4.0, -5.0]],
            },
            {
                "points": [[19.0, 20.0, 21.0], [22.0, 23.0, 24.0], [25.0, 26.0, 27.0]],
                "faces": [[0, 1, 2]],
                "minDistancePoints": [[-6.0, -7.0, -8.0]],
            },
            {
                "points": [[28.0, 29.0, 30.0], [31.0, 32.0, 33.0], [34.0, 35.0, 36.0]],
                "faces": [[0, 1, 2]],
                "minDistancePoints": [[-9.0, -10.0, -11.0]],
            },
        ],
        "selfIntersectionTest": [[0, 1], [1, 2, 3]],
    }
    return coil_dict


@pytest.fixture(scope="module")
def full_tcd_coil() -> TmsCoil:
    coil_casings = [
        TmsCoilModel(
            Msh(
                Nodes(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])),
                Elements(triangles=np.array([[0, 1, 2]]) + 1),
            ),
            np.array([[-0, -1, -2]]),
        ),
        TmsCoilModel(
            Msh(
                Nodes(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) + 9),
                Elements(triangles=np.array([[0, 1, 2]]) + 1),
            ),
            np.array([[-3, -4, -5]]),
        ),
        TmsCoilModel(
            Msh(
                Nodes(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) + 18),
                Elements(triangles=np.array([[0, 1, 2]]) + 1),
            ),
            np.array([[-6, -7, -8]]),
        ),
        TmsCoilModel(
            Msh(
                Nodes(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) + 27),
                Elements(triangles=np.array([[0, 1, 2]]) + 1),
            ),
            np.array([[-9, -10, -11]]),
        ),
    ]
    deformations = [
        TmsCoilTranslation(TmsCoilDeformationRange(-75, (-100, -50)), 0),
        TmsCoilTranslation(TmsCoilDeformationRange(0, (-50, 50)), 1),
        TmsCoilTranslation(TmsCoilDeformationRange(75, (50, 100)), 2),
        TmsCoilRotation(
            TmsCoilDeformationRange(180, (0, 360)),
            np.array([0, 0, 0]),
            np.array([0, 0, 1]),
        ),
    ]
    stimulator = TmsStimulator(
        "SimNIBS-Stimulator",
        "SimNIBS",
        1,
        [
            TmsWaveform(
                np.array([0, 100]), np.array([0.0, 0.99]), "biphasic", np.array([0, 1])
            )
        ],
    )
    elements = [
        DipoleElements(
            stimulator,
            np.array([[1.1, 2.1, 3.1]]),
            np.array([[100.1, 101.1, 102.1]]),
            "Dipole Elements",
            coil_casings[1],
            [deformations[0]],
        ),
        LineSegmentElements(
            stimulator,
            np.array([[1.2, 2.2, 3.2]]),
            np.array([[100.2, 101.2, 102.2]]),
            "Line Elements",
            coil_casings[2],
            [deformations[0], deformations[1]],
        ),
        SampledGridPointElements(
            stimulator,
            np.array([[[[1.4, 2.4, 3.4]]]]),
            np.array([[1, 0, 0, 1.4], [0, 1, 0, 1.4], [0, 0, 1, 1.4], [0, 0, 0, 1]]),
            "Digitized Grid Elements",
            coil_casings[3],
            [deformations[2], deformations[3]],
        ),
    ]
    coil = TmsCoil(
        elements,
        "SimNIBS-TMS-Coil",
        "SimNIBS",
        "V1",
        np.array([[0, 255], [0, 255], [0, 255]]),
        np.array([1, 1, 1]),
        coil_casings[0],
        [[0, 1], [1, 2, 3]],
    )
    return coil

@pytest.fixture(scope="module")
def tcd_json_schema():
    with open(file_finder.templates.tcd_json_schema, "r", encoding="utf-8") as fid:
        return json.loads(fid.read())
