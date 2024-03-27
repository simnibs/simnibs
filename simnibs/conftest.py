import os

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
    LineSegmentElements,
)
from simnibs.simulation.tms_coil.tms_coil_model import TmsCoilModel
from simnibs.simulation.tms_coil.tms_stimulator import TmsStimulator

from . import SIMNIBSDIR

def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)

@pytest.fixture(scope="module")
def sphere3_msh():
    fn = os.path.join(SIMNIBSDIR, "_internal_resources", "testing_files", "sphere3.msh")
    return Msh(fn=fn)


@pytest.fixture(scope="module")
def small_functional_3_element_coil() -> TmsCoil:
    casings = [
        TmsCoilModel(
            Msh(
                Nodes(
                    np.array(
                        [
                            [-20, -20, 0],
                            [-20, 20, 0],
                            [20, -20, 0],
                            [20, 20, 0],
                            [0, 0, 20],
                        ]
                    )
                ),
                Elements(
                    triangles=np.array(
                        [
                            [0, 1, 2],
                            [3, 2, 1],
                            [0, 1, 4],
                            [1, 3, 4],
                            [2, 3, 4],
                            [2, 0, 4],
                        ]
                    )
                    + 1
                ),
            ),
            None,
        ),
        TmsCoilModel(
            Msh(
                Nodes(
                    np.array(
                        [
                            [-20, -60, 0],
                            [-20, -20, 0],
                            [20, -60, 0],
                            [20, -20, 0],
                            [0, -40, 20],
                        ]
                    )
                ),
                Elements(
                    triangles=np.array(
                        [
                            [0, 1, 2],
                            [3, 2, 1],
                            [0, 1, 4],
                            [1, 3, 4],
                            [2, 3, 4],
                            [2, 0, 4],
                        ]
                    )
                    + 1
                ),
            ),
            None,
        ),
        TmsCoilModel(
            Msh(
                Nodes(
                    np.array(
                        [
                            [-20, 20, 0],
                            [-20, 60, 0],
                            [20, 20, 0],
                            [20, 60, 0],
                            [0, 40, 20],
                        ]
                    )
                ),
                Elements(
                    triangles=np.array(
                        [
                            [0, 1, 2],
                            [3, 2, 1],
                            [0, 1, 4],
                            [1, 3, 4],
                            [2, 3, 4],
                            [2, 0, 4],
                        ]
                    )
                    + 1
                ),
            ),
            None,
        ),
    ]
    deformations = [
        TmsCoilRotation(
            TmsCoilDeformationRange(0, (0, 90)),
            np.array([0, -20, 0]),
            np.array([40, -20, 0]),
        ),
        TmsCoilRotation(
            TmsCoilDeformationRange(0, (0, 90)),
            np.array([40, 20, 0]),
            np.array([0, 20, 0]),
        ),
    ]
    stimulator = TmsStimulator("SimNIBS-Stimulator")
    coil = TmsCoil(
        [
            LineSegmentElements(
                stimulator,
                np.array([[0, -10, 0], [10, 10, 0], [-10, 10, 0]]),
                casing=casings[0],
            ),
            LineSegmentElements(
                stimulator,
                np.array([[0, -50, 0], [10, -30, 0], [-10, -30, 0]]),
                casing=casings[1],
                deformations=[deformations[0]],
            ),
            LineSegmentElements(
                stimulator,
                np.array([[0, 30, 0], [10, 50, 0], [-10, 50, 0]]),
                casing=casings[2],
                deformations=[deformations[1]],
            ),
        ],
        limits=np.array([[-200, 200], [-200, 200], [-200, 200]]),
        resolution=np.array([10, 10, 10]),
    )

    return coil


@pytest.fixture(scope="module")
def small_self_intersecting_2_element_coil() -> TmsCoil:
    casings = [
        TmsCoilModel(
            Msh(
                Nodes(
                    np.array(
                        [
                            [-20, -20, 0],
                            [-20, 20, 0],
                            [20, -20, 0],
                            [20, 20, 0],
                            [0, 0, 20],
                        ]
                    )
                ),
                Elements(
                    triangles=np.array(
                        [
                            [0, 1, 2],
                            [3, 2, 1],
                            [0, 1, 4],
                            [1, 3, 4],
                            [2, 3, 4],
                            [2, 0, 4],
                        ]
                    )
                    + 1
                ),
            ),
            None,
        ),
        TmsCoilModel(
            Msh(
                Nodes(
                    np.array(
                        [[-2, -2, 0], [-2, 2, 0], [2, -2, 0], [2, 2, 0], [0, 0, 2]]
                    )
                    * 3
                    + np.array([0, 0, 10])
                ),
                Elements(
                    triangles=np.array(
                        [
                            [0, 1, 2],
                            [3, 2, 1],
                            [0, 1, 4],
                            [1, 3, 4],
                            [2, 3, 4],
                            [2, 0, 4],
                        ]
                    )
                    + 1
                ),
            ),
            None,
        ),
    ]
    deformations = [
        TmsCoilTranslation(TmsCoilDeformationRange(0, (-50, 50)), 0),
        TmsCoilTranslation(TmsCoilDeformationRange(0, (-50, 50)), 2),
        TmsCoilTranslation(TmsCoilDeformationRange(0, (-50, 50)), 0),
        TmsCoilTranslation(TmsCoilDeformationRange(0, (-50, 50)), 2),
    ]
    stimulator = TmsStimulator("SimNIBS-Stimulator")
    coil = TmsCoil(
        [
            LineSegmentElements(
                stimulator,
                np.array([[0, -10, 0], [10, 10, 0], [-10, 10, 0]]),
                casing=casings[0],
                deformations=[deformations[0], deformations[1]],
            ),
            LineSegmentElements(
                stimulator,
                np.array([[0, -1, 10], [1, 1, 10], [-1, 1, 10]]),
                casing=casings[1],
                deformations=[deformations[2], deformations[3]],
            ),
        ],
        limits=np.array([[-200, 200], [-200, 200], [-200, 200]]),
        resolution=np.array([10, 10, 10]),
        self_intersection_test=[[1, 2]],
    )

    return coil