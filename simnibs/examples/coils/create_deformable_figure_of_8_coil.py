"""
Example of how to create a deformable parametric figure of 8 TMS coil and save it in the tcd format.
The coil is constructed using line segments which reconstruct the windings of the coil.
"""
import os

import numpy as np
from simnibs.simulation.tms_coil.tms_coil import TmsCoil
from simnibs.simulation.tms_coil.tms_coil_deformation import TmsCoilRotation, TmsCoilDeformationRange

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
        The point index where the two spirals are split
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
        (0, -outer_diam / 2, -winding_casing_distance),
        (0, 0, -winding_casing_distance),
        connection_segment_count,
    ).T
    # Generate outgoing wire from the coil inside handle
    ending_wire_path = np.linspace(
        (0, 0, -winding_casing_distance),
        (0, -outer_diam / 2, -winding_casing_distance),
        connection_segment_count,
    ).T

    # Generate the connection from the incoming wire to the right spiral center
    wire_coil_2_connection = np.linspace(
        initial_wire_path[:, -1], spiral_2[:, 0], connection_segment_count
    ).T
    # Generate the connection from end of the right spiral to the outside start of the left spiral
    coil_coil_connection_1 = np.linspace(
        spiral_2[:, -1], (0, 0, -winding_casing_distance), connection_segment_count
    ).T
    coil_coil_connection_2 = np.linspace(
        (0, 0, -winding_casing_distance), spiral_1[:, 0], connection_segment_count
    ).T
    # Generate the connection from the inside of the left spiral to the outgoing wire
    coil_1_wire_connection = np.linspace(
        spiral_1[:, -1], ending_wire_path[:, 0], connection_segment_count
    ).T

    return (
        np.concatenate(
            (
                initial_wire_path,
                wire_coil_2_connection,
                spiral_2,
                coil_coil_connection_1,
                coil_coil_connection_2,
                spiral_1,
                coil_1_wire_connection,
                ending_wire_path,
            ),
            axis=1,
        ).T,
        initial_wire_path.shape[1] + wire_coil_2_connection.shape[1] + spiral_2.shape[1] + coil_coil_connection_1.shape[1]
    )


# Set up coil parameters
wire_diam = 4
segment_count = 1000
connection_segment_count = 20
outer_diam = 80
inner_diam = 20
element_distance = 90
winding_casing_distance = 4

wire_path, split_idx = figure_of_8_wire_path(
    wire_diam,
    segment_count,
    connection_segment_count,
    outer_diam,
    inner_diam,
    element_distance,
    winding_casing_distance,
)

# Creating line segments from the wire path
wire_dir = np.zeros(wire_path.shape)
wire_dir[:-1] = np.diff(wire_path, axis=0)
segment_points = wire_path + wire_dir / 2

# The limits of the a field of the coil, used for the transformation into nifti format
limits = [[-300.0, 300.0], [-200.0, 200.0], [-100.0, 300.0]]
# The resolution used when sampling to transform into nifti format
resolution = [1, 1, 1]

# Creating a example stimulator with a name, a brand and a maximum dI/dt
stimulator = TmsStimulator("Example Stimulator", "Example Stimulator Brand", 122.22e6)

# Creating coil element deformation to rotate the spirals around the y-axis
rotation_1 = TmsCoilRotation(TmsCoilDeformationRange(30, [0, 60]), [0, 0, -winding_casing_distance], [0, -1, -winding_casing_distance])
rotation_2 = TmsCoilRotation(TmsCoilDeformationRange(30, [0, 60]), [0, 0, -winding_casing_distance], [0, 1, -winding_casing_distance])

# Creating the line segment elements from the line segments
line_element_1 = LineSegmentElements(stimulator, segment_points[:split_idx], wire_dir[:split_idx], name="Figure_of_8_part_1", deformations=[rotation_1])
line_element_2 = LineSegmentElements(stimulator, segment_points[split_idx:], wire_dir[split_idx:], name="Figure_of_8_part_2", deformations=[rotation_2])

# Creating the TMS coil with its element, a name, a brand, a version, the limits and the resolution
tms_coil = TmsCoil(
    [line_element_1, line_element_2], "Example Coil", "Example Coil Brand", "V1.0", limits, resolution
)

# Generating a coil casing that has a specified distance from the coil windings
tms_coil.generate_element_casings(
    winding_casing_distance, winding_casing_distance / 2, True
)

# Write the coil to a tcd file
os.mkdir("coil_example")
tms_coil.write(os.path.join("coil_example", "deformable_figure_of_8_example_coil.tcd"))
