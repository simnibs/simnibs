from copy import  deepcopy
import json
import os
from typing import Any
import numpy as np
import pytest
from simnibs.simulation.coil.coil import Coil
from simnibs.simulation.coil.coil_deformation import CoilDeformation, CoilRotation, CoilTranslation
from simnibs.simulation.coil.coil_element import CoilDipoles, CoilLineElements, CoilLinePoints, CoilSampledGridElements
from simnibs.simulation.coil.coil_model import CoilModel
from simnibs.simulation.coil.tms_stimulator import TMSStimulator, TMSWaveform

from simnibs.utils import file_finder


from simnibs.mesh_tools.mesh_io import Elements, Msh, Nodes

from .... import SIMNIBSDIR


@pytest.fixture(scope='module')
def sphere3_msh():
    fn = os.path.join(
        SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere3.msh')
    return Msh(fn=fn)

@pytest.fixture(scope='module')
def minimal_tcd_coil_dict() -> dict[str, Any]:
    coil_dict = {
        "coilElementList":[
            {"type":3, "points": [[1, 2, 3]]}
        ]
    }
    return coil_dict

@pytest.fixture(scope='module')
def minimal_tcd_coil() -> Coil:
    coil = Coil(None, None ,None, None, None, None, [CoilLinePoints(None, None, None, np.array([[1, 2, 3]]), None, None)])
    return coil

@pytest.fixture(scope='module')
def medium_tcd_coil_dict() -> dict[str, Any]:
    coil_dict = {
        "coilElementList":[
            {"deformations":[0,1,2,3], "type":1, "points": [[1, 2, 3]], "values":[[0,1,2]]},
            {"elementCasing":0, "type":2, "points": [[1, 2, 3]], "values":[[0,1,2]]},
            {"stimulator":0, "type":3, "points": [[1, 2, 3]]},
            {"type":4, "data":[[[[1,2,3]]]], "affine":[[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]}
        ],
        "stimulatorList":[{"maxdIdt":1, "waveformList": [{"time":[0, 100], "signal":[0.0, 0.99]}]}],
        "deformList":[
            {"initial":-25, "range":[-100, -50], "type":"x"},
            {"initial":0, "range":[-50, 50], "type":"y"},
            {"initial":25, "range":[50, 100], "type":"z"},
            {"initial":180, "range":[0, 360], "type":"rot2p", "point1":[0,0,0], "point2":[0,0,1]}
        ],
        "coilModels":[
            {"points":[[1,2,3], [4,5,6], [7,8,9]], "faces":[[0,1,2]]}
        ]
    }
    return coil_dict

@pytest.fixture(scope='module')
def medium_tcd_coil() -> Coil:
    coil_casings = [
        CoilModel(Msh(Nodes(np.array([[1,2,3], [4,5,6], [7,8,9]])), Elements(triangles=np.array([[0,1,2]]) + 1)), None, None),
    ]
    deformations = [
        CoilTranslation(-25, (-100, -50), 0),
        CoilTranslation(0, (-50, 50), 1),
        CoilTranslation(25, (50, 100), 2),
        CoilRotation(180, (0, 360), np.array([0,0,0]), np.array([0,0,1]))
    ]
    stimulator = TMSStimulator(None, None, 1, [TMSWaveform(None, np.array([0, 100]), np.array([0.0, 0.99]), None)])
    coil = Coil(None, None, None, None, None, None,
                [CoilDipoles(None, None, [deformations[0], deformations[1], deformations[2], deformations[3]], np.array([[1, 2, 3]]), np.array([[0,1,2]]), None, None),
                 CoilLineElements(None, coil_casings[0], None, np.array([[1, 2, 3]]), np.array([[0,1,2]]), None, None),
                 CoilLinePoints(None, None, None, np.array([[1, 2, 3]]), stimulator, None),
                 CoilSampledGridElements(None, None, None, np.array([[[[1,2,3]]]]), np.array([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]), None, None),
                 ])
    return coil

@pytest.fixture(scope='module')
def full_tcd_coil_dict() -> dict[str, Any]:
    coil_dict = {
        "name":"SimNIBS-TMS-Coil",
        "brand":"SimNIBS",
        "version":"V1",
        "limits":[[0,255], [0,255], [0,255]],
        "resolution":[1,1,1],
        "coilCasing":0,
        "coilElementList":[
            {"name": "Dipole Elements", "stimulator":0, "elementCasing":1, "deformations":[0], "weights":[1], "type":1, "points": [[1, 2, 3]], "values":[[0,1,2]]},
            {"name": "Line Elements", "stimulator":0, "elementCasing":2, "deformations":[0, 1], "weights":[1], "type":2, "points": [[1, 2, 3]], "values":[[0,1,2]]},
            {"name": "Line Point Elements", "stimulator":0, "elementCasing":3, "deformations":[2], "weights":[1], "type":3, "points": [[1, 2, 3]]},
            {"name": "Digitized Grid Elements", "stimulator":0, "elementCasing":4, "deformations":[3], "weights":[1], "type":4, "data":[[[[1,2,3]]]], "affine":[[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]}
        ],
        "stimulatorList":[{"name":"SimNIBS-Stimulator", "brand": "SimNIBS", "maxdIdt":1, "waveformList": [{"name":"biphasic", "time":[0, 100], "signal":[0.0, 0.99], "fit":[0, 1]}]}],
        "deformList":[
            {"initial":-25, "range":[-100, -50], "type":"x"},
            {"initial":0, "range":[-50, 50], "type":"y"},
            {"initial":25, "range":[50, 100], "type":"z"},
            {"initial":180, "range":[0, 360], "type":"rot2p", "point1":[0,0,0], "point2":[0,0,1]}
        ],
        "coilModels":[
            {"points":[[1,2,3], [4,5,6], [7,8,9]], "faces":[[0,1,2]], "minDistancePoints":[[0,1,2]], "intersectPoints":[[1,1,1]]},
            {"points":[[2,3,4], [4,5,6], [7,8,9]], "faces":[[0,1,2]], "minDistancePoints":[[0,1,2]], "intersectPoints":[[1,1,1]]},
            {"points":[[3,4,5], [4,5,6], [7,8,9]], "faces":[[0,1,2]], "minDistancePoints":[[0,1,2]], "intersectPoints":[[1,1,1]]},
            {"points":[[4,5,6], [4,5,6], [7,8,9]], "faces":[[0,1,2]], "minDistancePoints":[[0,1,2]], "intersectPoints":[[1,1,1]]},
            {"points":[[5,6,7], [4,5,6], [7,8,9]], "faces":[[0,1,2]], "minDistancePoints":[[0,1,2]], "intersectPoints":[[1,1,1]]},
        ]
    }
    return coil_dict

@pytest.fixture(scope='module')
def full_tcd_coil() -> Coil:
    coil_casings = [
        CoilModel(Msh(Nodes(np.array([[1,2,3], [4,5,6], [7,8,9]])), Elements(triangles=np.array([[0,1,2]]) + 1)), np.array([[0,1,2]]), np.array([[1,1,1]])),
        CoilModel(Msh(Nodes(np.array([[2,3,4], [4,5,6], [7,8,9]])), Elements(triangles=np.array([[0,1,2]]) + 1)), np.array([[0,1,2]]), np.array([[1,1,1]])),
        CoilModel(Msh(Nodes(np.array([[3,4,5], [4,5,6], [7,8,9]])), Elements(triangles=np.array([[0,1,2]]) + 1)), np.array([[0,1,2]]), np.array([[1,1,1]])),
        CoilModel(Msh(Nodes(np.array([[4,5,6], [4,5,6], [7,8,9]])), Elements(triangles=np.array([[0,1,2]]) + 1)), np.array([[0,1,2]]), np.array([[1,1,1]])),
        CoilModel(Msh(Nodes(np.array([[5,6,7], [4,5,6], [7,8,9]])), Elements(triangles=np.array([[0,1,2]]) + 1)), np.array([[0,1,2]]), np.array([[1,1,1]]))
    ]
    deformations = [
        CoilTranslation(-25, (-100, -50), 0),
        CoilTranslation(0, (-50, 50), 1),
        CoilTranslation(25, (50, 100), 2),
        CoilRotation(180, (0, 360), np.array([0,0,0]), np.array([0,0,1]))
    ]
    stimulator = TMSStimulator("SimNIBS-Stimulator", "SimNIBS", 1, [TMSWaveform("biphasic", np.array([0, 100]), np.array([0.0, 0.99]), np.array([0, 1]))])
    coil = Coil("SimNIBS-TMS-Coil", "SimNIBS", "V1", np.array([[0,255], [0,255], [0,255]]), np.array([1,1,1]), coil_casings[0],
                [CoilDipoles("Dipole Elements", coil_casings[1], [deformations[0]], np.array([[1, 2, 3]]), np.array([[0,1,2]]), stimulator, np.array([1])),
                 CoilLineElements("Line Elements", coil_casings[2], [deformations[0], deformations[1]], np.array([[1, 2, 3]]), np.array([[0,1,2]]), stimulator, np.array([1])),
                 CoilLinePoints("Line Point Elements", coil_casings[3], [deformations[2]], np.array([[1, 2, 3]]), stimulator, np.array([1])),
                 CoilSampledGridElements("Digitized Grid Elements", coil_casings[4], [deformations[3]], np.array([[[[1,2,3]]]]), np.array([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]), stimulator, np.array([1])),
                 ])
    return coil

@pytest.fixture(scope='module')
def tcd_json_schema():
    with open(file_finder.templates.tcd_json_schema, "r") as fid:
        return json.loads(fid.read())