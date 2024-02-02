from enum import IntEnum


class TmsCoilElementTag(IntEnum):
    COIL_CASING = 1
    COIL_CASING_MIN_DISTANCE_POINTS = 2
    DIPOLES = 3
    LINE_ELEMENTS = 4
    SAMPLED_GRID_ELEMENTS = 5
    INDEX_OFFSET = 100
    BOUNDING_BOX = 99999
