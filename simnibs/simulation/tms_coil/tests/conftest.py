import json
import os
from typing import Any

import numpy as np
import pytest

from simnibs.mesh_tools.mesh_io import Elements, Msh, Nodes
from simnibs.simulation.tms_coil.tms_coil import TmsCoil
from simnibs.simulation.tms_coil.tms_coil_deformation import (
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
def sphere3_msh():
    fn = os.path.join(SIMNIBSDIR, "_internal_resources", "testing_files", "sphere3.msh")
    return Msh(fn=fn)


@pytest.fixture(scope="module")
def minimal_tcd_coil_dict() -> dict[str, Any]:
    coil_dict = {"coilElementList": [{"stimulator":0, "type": 1, "points": [[1, 2, 3]], "values":[[4,5,6]]}], "stimulatorList": [{"name": "SimNIBS-Stimulator"}]}
    return coil_dict


@pytest.fixture(scope="module")
def minimal_tcd_coil() -> TmsCoil:
    coil = TmsCoil(
        None,
        None,
        None,
        None,
        None,
        None,
        [DipoleElements(None, None, None, np.array([[1, 2, 3]]), np.array([[4, 5, 6]]), TmsStimulator("SimNIBS-Stimulator", None, None, None))],
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
                "points": [[1, 2, 3]],
                "values": [[0, 1, 2]],
            },
            {
                "stimulator": 0,
                "elementCasing": 0,
                "type": 2,
                "points": [[1, 2, 3]],
                "values": [[0, 1, 2]],
            },
            {
                "stimulator": 0,
                "type": 3,
                "data": [[[[1, 2, 3]]]],
                "affine": [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
            },
        ],
        "stimulatorList": [
            {"maxdIdt": 1, "waveformList": [{"time": [0, 100], "signal": [0.0, 0.99]}]}
        ],
        "deformList": [
            {"initial": -25, "range": [-100, -50], "type": "x"},
            {"initial": 0, "range": [-50, 50], "type": "y"},
            {"initial": 25, "range": [50, 100], "type": "z"},
            {
                "initial": 180,
                "range": [0, 360],
                "type": "rot2p",
                "point1": [0, 0, 0],
                "point2": [0, 0, 1],
            },
        ],
        "coilModels": [
            {"points": [[1, 2, 3], [4, 5, 6], [7, 8, 9]], "faces": [[0, 1, 2]]}
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
            None,
        ),
    ]
    deformations = [
        TmsCoilTranslation(-25, (-100, -50), 0),
        TmsCoilTranslation(0, (-50, 50), 1),
        TmsCoilTranslation(25, (50, 100), 2),
        TmsCoilRotation(180, (0, 360), np.array([0, 0, 0]), np.array([0, 0, 1])),
    ]
    stimulator = TmsStimulator(
        None,
        None,
        1,
        [TmsWaveform(None, np.array([0, 100]), np.array([0.0, 0.99]), None)],
    )
    coil = TmsCoil(
        None,
        None,
        None,
        None,
        None,
        None,
        [
            DipoleElements(
                None,
                None,
                [deformations[0], deformations[1], deformations[2], deformations[3]],
                np.array([[1, 2, 3]]),
                np.array([[0, 1, 2]]),
                stimulator,
            ),
            LineSegmentElements(
                None,
                coil_casings[0],
                None,
                np.array([[1, 2, 3]]),
                np.array([[0, 1, 2]]),
                stimulator,
            ),
            SampledGridPointElements(
                None,
                None,
                None,
                np.array([[[[1, 2, 3]]]]),
                np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
                stimulator,
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
        "limits": [[0, 255], [0, 255], [0, 255]],
        "resolution": [1, 1, 1],
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
                "affine": [[1.0, 0.0, 0.0, 1.4], [0.0, 1.0, 0.0, 1.4], [0.0, 0.0, 1.0, 1.4], [0.0, 0.0, 0.0, 1.0]],
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
                        "time": [0, 100],
                        "signal": [0.0, 0.99],
                        "fit": [0, 1],
                    }
                ],
            }
        ],
        "deformList": [
            {"initial": -25, "range": [-100, -50], "type": "x"},
            {"initial": 0, "range": [-50, 50], "type": "y"},
            {"initial": 25, "range": [50, 100], "type": "z"},
            {
                "initial": 180,
                "range": [0, 360],
                "type": "rot2p",
                "point1": [0, 0, 0],
                "point2": [0, 0, 1],
            },
        ],
        "coilModels": [
            {
                "points": [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                "faces": [[0, 1, 2]],
                "minDistancePoints": [[-0, -1, -2]],
                "intersectPoints": [[1, 1, 1]],
            },
            {
                "points": [[10, 11, 12], [13, 14, 15], [16, 17, 18]],
                "faces": [[0, 1, 2]],
                "minDistancePoints": [[-3, -4, -5]],
                "intersectPoints": [[2, 2, 2]],
            },
            {
                "points": [[19, 20, 21], [22, 23, 24], [25, 26, 27]],
                "faces": [[0, 1, 2]],
                "minDistancePoints": [[-6, -7, -8]],
                "intersectPoints": [[3, 3, 3]],
            },
            {
                "points": [[28, 29, 30], [31, 32, 33], [34, 35, 36]],
                "faces": [[0, 1, 2]],
                "minDistancePoints": [[-9, -10, -11]],
                "intersectPoints": [[4, 4, 4]],
            },
        ],
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
            np.array([[1, 1, 1]]),
        ),
        TmsCoilModel(
            Msh(
                Nodes(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) + 9),
                Elements(triangles=np.array([[0, 1, 2]]) + 1),
            ),
            np.array([[-3, -4, -5]]),
            np.array([[2, 2, 2]]),
        ),
        TmsCoilModel(
            Msh(
                Nodes(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) + 18),
                Elements(triangles=np.array([[0, 1, 2]]) + 1),
            ),
            np.array([[-6, -7, -8]]),
            np.array([[3, 3, 3]]),
        ),
        TmsCoilModel(
            Msh(
                Nodes(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]) + 27),
                Elements(triangles=np.array([[0, 1, 2]]) + 1),
            ),
            np.array([[-9, -10, -11]]),
            np.array([[4, 4, 4]]),
        )
    ]
    deformations = [
        TmsCoilTranslation(-25, (-100, -50), 0),
        TmsCoilTranslation(0, (-50, 50), 1),
        TmsCoilTranslation(25, (50, 100), 2),
        TmsCoilRotation(180, (0, 360), np.array([0, 0, 0]), np.array([0, 0, 1])),
    ]
    stimulator = TmsStimulator(
        "SimNIBS-Stimulator",
        "SimNIBS",
        1,
        [
            TmsWaveform(
                "biphasic", np.array([0, 100]), np.array([0.0, 0.99]), np.array([0, 1])
            )
        ],
    )
    coil = TmsCoil(
        "SimNIBS-TMS-Coil",
        "SimNIBS",
        "V1",
        np.array([[0, 255], [0, 255], [0, 255]]),
        np.array([1, 1, 1]),
        coil_casings[0],
        [
            DipoleElements(
                "Dipole Elements",
                coil_casings[1],
                [deformations[0]],
                np.array([[1.1, 2.1, 3.1]]),
                np.array([[100.1, 101.1, 102.1]]),
                stimulator,
            ),
            LineSegmentElements(
                "Line Elements",
                coil_casings[2],
                [deformations[0], deformations[1]],
                np.array([[1.2, 2.2, 3.2]]),
                np.array([[100.2, 101.2, 102.2]]),
                stimulator,
            ),
            SampledGridPointElements(
                "Digitized Grid Elements",
                coil_casings[3],
                [deformations[2], deformations[3]],
                np.array([[[[1.4, 2.4, 3.4]]]]),
                np.array([[1, 0, 0, 1.4], [0, 1, 0, 1.4], [0, 0, 1, 1.4], [0, 0, 0, 1]]),
                stimulator,
            ),
        ],
    )
    return coil


@pytest.fixture(scope="module")
def tcd_json_schema():
    with open(file_finder.templates.tcd_json_schema, "r") as fid:
        return json.loads(fid.read())
