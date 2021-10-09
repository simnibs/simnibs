import numpy as np
import pytest

from .. import charm_main


def generate_label_arr(ndim=2):
    """Generate an array of labels."""
    size = ndim * (10,)
    arr = np.zeros(size, dtype=np.uint16)
    arr[5:] = 1
    arr[0, :5] = 2
    arr[5, 5] = 3
    arr[7:, 7:] = 5
    return arr


# all cases invoke remapping!
test_label_unassigned_elements_inputs = (
    (2, [0, 1, 2, 3, 5], 3, None, 0),
    (3, None, 3, None, 1),
    (3, None, 3, [1], 0),
    (5, None, 3, None, 1),
    (5, None, 3, [1], 5),  # should show a warning
)


# @pytest.mark.filterwarnings("ignore:Some elements could not be labeled")
@pytest.mark.parametrize(
    "label_unassign, labels, window_size, ignore_labels, expected_label",
    test_label_unassigned_elements_inputs,
)
def test_label_unassigned_elements(
    label_unassign, labels, window_size, ignore_labels, expected_label,
):
    label_arr = generate_label_arr(3)
    expected_arr = label_arr.copy()
    expected_arr[label_arr == label_unassign] = expected_label
    np.testing.assert_allclose(
        expected_arr,
        charm_main.label_unassigned_elements(
            label_arr, label_unassign, labels, window_size, ignore_labels
        ),
    )
