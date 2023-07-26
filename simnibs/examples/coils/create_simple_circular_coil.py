"""
Example of how to create a simple parametric TMS coil and save it in the tcd format.
The coil is constructed using line segments which reconstruct the windings of the coil.
"""

import numpy as np
from simnibs.simulation.tms_coil.tms_coil import TmsCoil

from simnibs.simulation.tms_coil.tms_coil_element import LineSegmentElements
from simnibs.simulation.tms_coil.tms_stimulator import TmsStimulator


# Set up coil parameters
wire_diam = 4
segment_count = 1000
diam = 50
winding_casing_distance = 4

# Generating the angles of the circle points
angles = np.linspace(0, 2 * np.pi, segment_count, endpoint=False)
# Generating the cartesian coordinates of the circle points
wire_path = path = np.array(
    [
        diam / 2 * np.cos(angles),
        diam / 2 * np.sin(angles),
        np.full_like(angles, -winding_casing_distance),
    ]
).T

# The limits of the a field of the coil, used for the transformation into nifti format
limits = [[-300.0, 300.0], [-200.0, 200.0], [-100.0, 300.0]]
# The resolution used when sampling to transform into nifti format
resolution = [1, 1, 1]

# Creating a example stimulator with a name, a brand and a maximum dI/dt
stimulator = TmsStimulator("Example Stimulator", "Example Stimulator Brand", 122.22)

# Creating the line segments from a list of wire path points
line_element = LineSegmentElements(stimulator, wire_path, name="Circular")
# Creating the TMS coil with its element, a name, a brand, a version, the limits and the resolution
tms_coil = TmsCoil(
    [line_element], "Example Coil", "Example Coil Brand", "V1.0", limits, resolution
)

# Generating a coil casing that has a specified distance from the coil windings
tms_coil.generate_element_casings(
    winding_casing_distance, winding_casing_distance / 20, True
)

# Write the coil to a tcd file
tms_coil.write("circular_example_coil.tcd")
