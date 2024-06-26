import numpy as np
from scipy.integrate import simpson


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
    auc = simpson(y=sensitivity, x=specificity_inv)

    # import matplotlib.pyplot as plt
    # import matplotlib
    # matplotlib.use("Qt5Agg")
    # plt.plot(specificity_inv, sensitivity, "-o")
    # plt.plot(np.linspace(0, 1, 2), np.linspace(0,1,2), "k--")
    # plt.xlabel("1-spec")
    # plt.ylabel("sens")

    return auc


def ROC(e1, e2, threshold, focal=False):
    """
    Determine distance of receiving operating characteristic (ROC) to optimal solution,
    i.e. sensitivity = specificity = 1, for a given electric field and threshold(s).

    Parameters
    ----------
    e1 : np.ndarray of float [n_ele_roi]
        Electric field magnitude in ROI
    e2 : np.ndarray of float [n_ele_non-roi]
        Electric field magnitude in non-ROI
    threshold : float or list of float [2]
        Threshold(s) the receiving operating characteristic is computed for. Either a single value or two values can be
        given. If two values are given, the first threshold may not be exceeded in the non-ROI and the second value
        should be exceeded in the ROI in the sense of a safety gap.
    focal : bool, optional, default: True
        Distance to optimal solution.
            True : quantifies focality of a solution.
                   (goal: sensitivity = specificity = 1)
            False: quantifies intensity in ROI while also increasing the field in the non ROI
                   (goal: sensitivity(ROI) = 1, sensitivity(nonROI) = 1)
    Return
    ------
    dist : float
        Distance of the ROC value to the optimal solution (optimality depends on "focal" flag)
    """
    if (type(threshold) != list and type(threshold) != np.ndarray):
        threshold = [threshold, threshold]

    if len(threshold) == 1:
        threshold = [threshold[0], threshold[0]]

    threshold_array = np.linspace(threshold[0], threshold[1], 2)
    sensitivity = np.zeros(len(threshold_array))
    specificity_inv = np.zeros(len(threshold_array))
    n_roi_1 = len(e1)
    n_roi_2 = len(e2)

    for i, t in enumerate(threshold_array):
        # true positive
        sensitivity[i] = np.sum(e1 >= t) / n_roi_1

        # false positive (1-specificity)
        specificity_inv[i] = np.sum(e2 >= t) / n_roi_2

    # sort it for thresholds
    sort_idx = np.argsort(threshold_array)
    specificity_inv = specificity_inv[sort_idx]
    sensitivity = sensitivity[sort_idx]
    threshold_array = threshold_array[sort_idx]

    # get sensitivity and specificity_inv value where threshold conditions are fulfilled
    sens_idx_tmp = np.where(threshold_array >= threshold[1])[0]
    if len(sens_idx_tmp) == 0:
        sens_thresh = 0.
    else:
        sens_idx = sens_idx_tmp[0]
        sens_thresh = sensitivity[sens_idx]

    spec_inv_idx_tmp = np.where(threshold_array <= threshold[0])[0]
    if len(spec_inv_idx_tmp) == 0:
        spec_inv_thresh = 1.
    else:
        spec_inv_idx = spec_inv_idx_tmp[-1]
        spec_inv_thresh = specificity_inv[spec_inv_idx]

    # determine distance to optimum
    if focal:
        dist = np.linalg.norm([spec_inv_thresh, sens_thresh - 1])
    else:
        dist = np.linalg.norm([spec_inv_thresh, sens_thresh])

    return dist
