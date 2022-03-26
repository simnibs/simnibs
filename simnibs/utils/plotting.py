import matplotlib as mpl
import matplotlib.cm
import numpy as np
import os
import nibabel as nib
from scipy import ndimage
from . import brainsprite_helper
import scipy.ndimage.morphology as mrph

se = ndimage.generate_binary_structure(3,3)

cmaplist_final = [(0, 0, 0, 0),
                 (230/255, 230/255, 230/255, 1.0),
                 (129/255, 129/255, 129/255, 1.0),
                 (104/255, 163/255, 255/255, 1.0),
                 (255/255, 239/255, 179/255, 1.0),
                 (255/255, 166/255, 133/255, 1.0),
                 (255/255, 240/255, 0, 1.0),
                 (255/255, 239/255, 179/255, 1.0),
                 (255/255, 138/255, 57/255, 1.0),
                 (0, 65/255, 142/255, 1.0)]

cmaplist_affine = [(0, 0, 0, 0),
                   (0, 0, 0, 0),
                   (0, 0, 0, 0),
                   (0, 0, 0, 0),
                   (0, 0, 0, 0),
                   (255/255, 166/255, 133/255, 0.3),
                   (255/255, 239/255, 179/255, 1.0)]

cmap_final = mpl.colors.LinearSegmentedColormap.from_list(
            'final_segmentation', cmaplist_final, len(cmaplist_final))

cmap_affine = mpl.colors.LinearSegmentedColormap.from_list(
            'template', cmaplist_affine, len(cmaplist_affine))


def _mosaic_images(dim, num_of_blocks=5):

    patternX = np.sin(np.arange(dim[0])/dim[0] * 2 * np.pi * (num_of_blocks / 2))
    patternY = np.sin(np.arange(dim[1])/dim[1] * 2 * np.pi * (num_of_blocks / 2))
    patternZ = np.sin(np.arange(dim[2])/dim[2] * 2 * np.pi * (num_of_blocks / 2))
    patternXY = np.outer(patternX, patternY)
    mask = np.zeros(dim)

    for d in range(dim[2]):
        mask[:,:,d] = patternXY * patternZ[d]

    return mask

def _registration_overlay(T2):

    T2_data = nib.load(T2)
    bufferT2 = T2_data.get_fdata()
    # normalize
    bufferT2 = (bufferT2 - bufferT2.min())/(bufferT2.max() - bufferT2.min())
    d = ndimage.sobel(bufferT2, axis=0, mode='constant')**2 + \
        ndimage.sobel(bufferT2, axis=1, mode='constant')**2 + \
        ndimage.sobel(bufferT2, axis=2, mode='constant')**2
    d = np.sqrt(d)
    p95 = np.percentile(d, 95)
    d[d>p95] = p95
    d /= d.max()
    #Now d is between 0 and one. Let's first substract .7 to map .7 to zero,
    #and then multiply by three to map the max (0.3) to near 1
    d = 3*(d-0.7)
    #Mask out everything below zero
    d[d<0] = np.NAN

    return nib.Nifti1Image(np.squeeze(d), T2_data.affine)


def _final_overlay(labeling):
    tissues_data = nib.load(labeling)
    tissues_buffer = tissues_data.get_fdata()
    labels = np.unique(tissues_buffer)
    labels = labels[1:]
    for l in labels:
        tissue_mask = (tissues_buffer == l)
        ero = mrph.binary_erosion(np.squeeze(tissue_mask), se, 1)
        tissues_buffer[ero] = 0

    return nib.Nifti1Image(np.uint16(np.squeeze(tissues_buffer)), tissues_data.affine)

def _cap_intensities(im):
    im_data = nib.load(im)
    im_buffer = im_data.get_fdata()
    capped_buffer = (im_buffer - im_buffer.min()) / (im_buffer.max() - im_buffer.min())
    return nib.Nifti1Image(np.squeeze(capped_buffer), im_data.affine)

def write(sub_files, templates):
    names = []
    imgs = []
    cmaps = []
    interpolation_order = []
    if os.path.exists(sub_files.reference_volume):
        imgs.append(_cap_intensities(sub_files.reference_volume))
        cmaps.append(matplotlib.cm.gray)
        interpolation_order.append(1)
        names.append('Reference volume (T1)')

    if os.path.exists(sub_files.T2_reg):
        imgs.append(_registration_overlay(sub_files.T2_reg))
        cmaps.append(matplotlib.cm.inferno)
        interpolation_order.append(1)
        names.append('Registration overlay')

    if 0:#os.path.exists(sub_files.template_coregistered):
        imgs.append(nib.load(sub_files.template_coregistered))
        cmaps.append(cmap_affine)
        interpolation_order.append(0)
        names.append('Coregistered template')

    if os.path.exists(sub_files.final_labels):
        imgs.append(_final_overlay(sub_files.final_labels))
        cmaps.append(cmap_final)
        interpolation_order.append(0)
        names.append('Tissue labels')

    brainsprite_helper.write_viewer(imgs, cmaps, interpolation_order,
                                    sub_files.viewer, templates.html_template,
                                    templates.jquery, templates.brainsprite,
                                    names=names)


