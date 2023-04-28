from pathlib import Path
from tempfile import TemporaryDirectory

from nibabel.affines import apply_affine
import numpy as np
import scipy.io
import pytest

from simnibs.utils import csv_reader, spatial_transform
from simnibs.eeg import utils_fieldtrip

@pytest.mark.parametrize('use_trans', [False, True])
def test_prepare_montage(use_trans):
    """Create a MAT file with some electrode data, convert to CSV, and check.
    """
    rng = np.random.default_rng(0)

    n = 10
    # data for creating a MAT file
    elec = dict(
        elecpos = rng.random((n,3)),
        chantype = np.array([["eeg"]] * n, dtype=object),
        label = np.array(list(map(str, range(n))), dtype=object, ndmin=2).T,
        unit = "mm",
    )

    trans = spatial_transform.matrix_from_params(0.1, 0.2, 0.3) if use_trans else None

    with TemporaryDirectory() as d:
        d = Path(d)
        elec_file = d / "elec.mat"
        montage_file = d / "montage.csv"

        scipy.io.savemat(elec_file, dict(elec=elec))
        if use_trans:
            trans_file = d / "trans.mat"
            scipy.io.savemat(trans_file, dict(trans=trans))
            pos = apply_affine(trans, elec["elecpos"])
        else:
            trans_file = None
            pos = elec["elecpos"]

        utils_fieldtrip.prepare_montage(montage_file, elec_file, trans_file)

        elec_type, elec_pos, _, elec_name, _, _ = csv_reader.read_csv_positions(montage_file)

        np.testing.assert_allclose(pos, elec_pos)
        np.testing.assert_array_equal(elec["label"].squeeze().astype(str), np.array(elec_name))
        np.testing.assert_array_equal(
            elec["chantype"].squeeze() == "eeg",
            np.array(elec_type) == "Electrode"
        )



