"""
Example of how to create a simple parametric TMS coil consisting of a circular and a figure of 8 coil stacked.
The TMS coil is than saved in the tcd format.
The coil is constructed using line segments which reconstruct the windings of the coil.
"""

import numpy as np
from simnibs.simulation.tms_coil.tms_coil import TmsCoil

from simnibs.simulation.tms_coil.tms_coil_element import LineSegmentElements
from simnibs.simulation.tms_coil.tms_stimulator import TmsStimulator


def spiral(outer_diam: float, inner_diam: float, wire_diam: float, segment_count: int):
    """Creates a spiral path based on the inner and outer diameter of the coil as well as
    the diameter of the wire and the segment count.

    Parameters
    ----------
    outer_diam : float
        The outer diameter of the spiral
    inner_diam : float
        The inner diameter of the spiral
    wire_diam : float
        The diameter of the wire
    segment_count : int
        The amount of segments used to represent the spiral

    Returns
    -------
        The wire path of the spiral
    """
    # Calculates the maximal spiral angle
    phi_max = 2 * np.pi * (outer_diam / 2 - inner_diam / 2) / (wire_diam / 2)

    # Evenly spaced angles between 30 degrees and phi_max
    phi = np.linspace(np.pi / 6, phi_max, segment_count)
    # Calculates the radius of every line segment in the spiral
    radius_loop = outer_diam / 2 - wire_diam / 2 * phi / (2 * np.pi)
    # Calculates the cartesian coordinates of the spiral
    path = np.array(
        [radius_loop * np.cos(phi), radius_loop * np.sin(phi), radius_loop * 0]
    )

    return path


def figure_of_8_wire_path(
    wire_diam: float,
    segment_count: int,
    connection_segment_count: int,
    outer_diam: float,
    inner_diam: float,
    element_distance: float,
    winding_casing_distance: float,
):
    """Generates the windings of a figure of 8 coil

    Parameters
    ----------
    wire_diam : float
        The diameter of the wire
    segment_count : int
       The amount of segments used to represent each of the two spirals
    connection_segment_count : int
        The amount of segments used to connect the different parts of the coil
    outer_diam : float
        The outer diameter of the spiral
    inner_diam : float
        The inner diameter of the spiral
    element_distance : float
        The center to center distance of the two spirals
    winding_casing_distance : float
        The distance of the casing to the windings of the coil

    Returns
    -------
        The windings of a figure of 8 coil
    """
    # Generate left spiral of the coil
    path = spiral(outer_diam, inner_diam, wire_diam, segment_count)
    spiral_1 = (
        path + np.array((-element_distance / 2, 0, -winding_casing_distance))[:, None]
    )

    # Generate right spiral of the coil
    path[1] = -path[1]
    spiral_2 = np.fliplr(
        path * np.array((-1, 1, 1))[:, None]
        + np.array((element_distance / 2, 0, -winding_casing_distance))[:, None]
    )

    # Generate incoming wire to the coil inside handle
    initial_wire_path = np.linspace(
        (0.1, -outer_diam / 2, -winding_casing_distance),
        (0, 0, -winding_casing_distance),
        connection_segment_count,
    ).T
    # Generate outgoing wire from the coil inside handle
    ending_wire_path = np.linspace(
        (-0.1, 0, -winding_casing_distance),
        (0, -outer_diam / 2, -winding_casing_distance),
        connection_segment_count,
    ).T

    # Generate the connection from the incoming wire to the right spiral center
    wire_coil_2_connection = np.linspace(
        initial_wire_path[:, -1], spiral_2[:, 0], connection_segment_count
    ).T
    # Generate the connection from end of the right spiral to the outside start of the left spiral
    coil_coil_connection = np.linspace(
        spiral_2[:, -1], spiral_1[:, 0], connection_segment_count
    ).T
    # Generate the connection from the inside of the left spiral to the outgoing wire
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


# Set up coil parameters
wire_diam = 4
segment_count = 1000
connection_segment_count = 20
outer_diam = 80
inner_diam = 20
circular_diam = 60
element_distance = 90
winding_casing_distance = 4
figure_of_8_to_circular_distance = 0.5

wire_path_figure_of_8 = figure_of_8_wire_path(
    wire_diam,
    segment_count,
    connection_segment_count,
    outer_diam,
    inner_diam,
    element_distance,
    winding_casing_distance,
)

# Generating the angles of the circle points
angles = np.linspace(0, 2 * np.pi, segment_count, endpoint=False)
# Generating the cartesian coordinates of the circle points
wire_path_circular = np.array(
    [
        circular_diam / 2 * np.cos(angles),
        circular_diam / 2 * np.sin(angles),
        np.full_like(angles, -winding_casing_distance - figure_of_8_to_circular_distance),
    ]
).T

# The limits of the a field of the coil, used for the transformation into nifti format
limits = [[-300.0, 300.0], [-200.0, 200.0], [-100.0, 300.0]]
# The resolution used when sampling to transform into nifti format
resolution = [1, 1, 1]

# Creating two example stimulator with a name, a brand and a maximum dI/dt
stimulator_1 = TmsStimulator("Example Stimulator 1", "Example Stimulator Brand", 122.22e6)
stimulator_2 = TmsStimulator("Example Stimulator 2", "Example Stimulator Brand", 102.45e6)

# Creating the line segments from a list of wire path points
line_element_figure_of_8 = LineSegmentElements(stimulator_1, wire_path_figure_of_8, name="Figure_of_8")
line_element_circular = LineSegmentElements(stimulator_2, wire_path_circular, name="Circular")

# Creating the TMS coil with its element, a name, a brand, a version, the limits and the resolution
tms_coil = TmsCoil(
    [line_element_figure_of_8, line_element_circular], "Example Coil", "Example Coil Brand", "V1.0", limits, resolution
)

# Generating a coil casing that has a specified distance from the coil windings
tms_coil.generate_element_casings(
    winding_casing_distance, winding_casing_distance / 5, True, combined_casing=True
)

# Write the coil to a tcd file
tms_coil.write("two_stimulator_example_coil.tcd")
