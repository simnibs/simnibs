import numpy as np
from scipy.integrate import simps


def integral_focality(e1, e2, v1, v2):
    """
    Calculates the integral focality measure from Fernandez-Corazza et al. (2020) eq. (14).

    Fernandez-Corazza, M., Turovets, S., & Muravchik, C. H. (2020).
    Unification of optimal targeting methods in transcranial electrical stimulation.
    Neuroimage, 209, 116403.

    Parameters
    ----------
    e1 : np.ndarray of float [n_ele_roi]
        Electric field magnitude in ROI
    e2 : np.ndarray of float [n_ele_non-roi]
        Electric field magnitude in non-ROI
    v1 : float
        Volume (or area) of ROI
    v2 : float
        Volume (or area) of non-ROI

    Returns
    -------
    integral_focality : float
        Integral focality measure
    """
    return (np.mean(e1) / v1) / np.sqrt(np.mean(e2) / v2)


def AUC(e1, e2):
    """
    Computes the area under curve (AUC) (or receiver operating characteristic, ROC) measure.

    Parameters
    ----------
    e1 : np.ndarray of float [n_ele_roi]
        Electric field magnitude in ROI
    e2 : np.ndarray of float [n_ele_non-roi]
        Electric field magnitude in non-ROI

    Returns
    -------
    AUC : float
        Area under curve (receiver operating characteristic)
    """
    threshold = np.linspace(0.001, 1, 500)
    sensitivity = np.zeros(len(threshold))
    specificity_inv = np.zeros(len(threshold))
    e_max = np.max((np.max(e1), np.max(e2)))
    e1_norm = e1 / e_max
    e2_norm = e2 / e_max
    n_roi_1 = len(e1)
    n_roi_2 = len(e2)

    for i, t in enumerate(threshold):
        # true positive
        sensitivity[i] = np.sum(e1_norm >= t) / n_roi_1

        # false positive (1-specificity)
        specificity_inv[i] = np.sum(e2_norm >= t) / n_roi_2

    sort_idx = np.argsort(specificity_inv)

    specificity_inv = specificity_inv[sort_idx]
    sensitivity = sensitivity[sort_idx]

    specificity_inv = np.append(specificity_inv, 1)
    sensitivity = np.append(sensitivity, 1)

    # calculate area under curve
    auc = simps(y=sensitivity, x=specificity_inv)

    # import matplotlib.pyplot as plt
    # import matplotlib
    # matplotlib.use("Qt5Agg")
    # plt.plot(specificity_inv, sensitivity, "-o")
    # plt.plot(np.linspace(0, 1, 2), np.linspace(0,1,2), "k--")
    # plt.xlabel("1-spec")
    # plt.ylabel("sens")

    return auc
