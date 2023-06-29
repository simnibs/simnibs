"""
Example of how to create a simple parametric TMS coil and save it in the tcd format.
The coil is constructed using line segments which reconstruct the windings of the coil.
"""

import numpy as np
from simnibs.simulation.tms_coil.tms_coil import TmsCoil

from simnibs.simulation.tms_coil.tms_coil_element import LineSegmentElements
from simnibs.simulation.tms_coil.tms_stimulator import TmsStimulator


def spiral(outer_diam, inner_diam, wire_diam, segment_count):
    phi_max = 2 * np.pi * (outer_diam / 2 - inner_diam / 2) / (wire_diam / 2)

    phi = np.linspace(np.pi / 6, phi_max, segment_count)
    radius_loop = outer_diam / 2 - wire_diam / 2 * phi / (2 * np.pi)
    path = np.array(
        [radius_loop * np.cos(phi), radius_loop * np.sin(phi), radius_loop * -0.1]
    )

    return path


def figure_of_8_wire_path(
    wire_diam,
    segment_count,
    connection_segment_count,
    outer_diam,
    inner_diam,
    element_distance,
):
    path = spiral(outer_diam, inner_diam, wire_diam, segment_count)
    spiral_1 = path + np.array((-element_distance / 2, 0, 0))[:, None]

    path[1] = -path[1]
    spiral_2 = np.fliplr(
        path * np.array((-1, 1, 1))[:, None]
        + np.array((element_distance / 2, 0, 0))[:, None]
    )

    initial_wire_path = np.linspace(
        (0.1, -outer_diam / 2, 0), (0.1, 0, 0), connection_segment_count
    ).T
    ending_wire_path = np.linspace(
        (-0.1, 0, 0), (-0.1, -outer_diam / 2, 0), connection_segment_count
    ).T

    wire_coil_2_connection = np.linspace(
        initial_wire_path[:, -1], spiral_2[:, 0], connection_segment_count
    ).T
    coil_coil_connection = np.linspace(
        spiral_2[:, -1], spiral_1[:, 0], connection_segment_count
    ).T
    coil_1_wire_connection = np.linspace(
        spiral_1[:, -1], ending_wire_path[:, 0], connection_segment_count
    ).T

    return np.concatenate(
        (
            initial_wire_path,
            wire_coil_2_connection,
            spiral_2,
            coil_coil_connection,
            spiral_1,
            coil_1_wire_connection,
            ending_wire_path,
        ),
        axis=1,
    ).T


wire_diam = 0.5
segment_count = 500
connection_segment_count = 20
outer_diam = 10
inner_diam = 1
element_distance = 10.5

wire_path = figure_of_8_wire_path(
    wire_diam,
    segment_count,
    connection_segment_count,
    outer_diam,
    inner_diam,
    element_distance,
)


limits = [
    [np.min(wire_path[:, 0]) - 50, np.max(wire_path[:, 0]) + 50],
    [np.min(wire_path[:, 1]) - 50, np.max(wire_path[:, 1]) + 50],
    [np.min(wire_path[:, 2]) - 50, np.max(wire_path[:, 2]) + 50],
]

stimulator = TmsStimulator("Example Stimulator", "Example Stimulator Brand", 122.22)
line_element = LineSegmentElements(stimulator, wire_path, name="Spirals")
tms_coil = TmsCoil(
    [line_element], "Example Coil", "Example Coil Brand", "V1.0", limits, [1, 1, 1]
)

tms_coil.write("example_coil.tcd")
