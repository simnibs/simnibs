import matplotlib as mpl
from nilearn import image
from nilearn import plotting
import numpy as np
from scipy import ndimage
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
            'Custom cmap', cmaplist_final, len(cmaplist_final))

cmap_affine = mpl.colors.LinearSegmentedColormap.from_list(
            'Custom cmap', cmaplist_affine, len(cmaplist_affine))

def viewer_registration(T1, T2, save_path):

    T1_data = image.load_img(T1)
    bufferT1 = T1_data.get_fdata()
    # normalize
    bufferT1 = (bufferT1 - bufferT1.min()) / (bufferT1.max() - bufferT1.min())
    cappedT1 = image.new_img_like(T1_data, bufferT1)

    T2_data = image.load_img(T2)
    bufferT2 = T2_data.get_fdata()
    # normalize
    bufferT2 = (bufferT2 - bufferT2.min())/(bufferT2.max() - bufferT2.min())
    d = ndimage.sobel(bufferT2, axis=0, mode='constant')**2 + \
        ndimage.sobel(bufferT2, axis=1, mode='constant')**2 + \
        ndimage.sobel(bufferT2, axis=2, mode='constant')**2
    d = np.sqrt(d)
    d = 1 / (1 + np.exp(-d)) - 0.8
    d[d<0] = 0
    cappedT2 = image.new_img_like(T2_data, d)
    html_view = plotting.view_img(cappedT2, bg_img=cappedT1, cmap='inferno', symmetric_cmap=False,
                                  colorbar=False, opacity=0.8, dim=-1)
    html_view.save_as_html(save_path)

def viewer_affine(T1, template, save_path):
    T1_data = image.load_img(T1)
    bufferT1 = T1_data.get_fdata()
    # normalize
    bufferT1 = (bufferT1 - bufferT1.min()) / (bufferT1.max() - bufferT1.min())
    cappedT1 = image.new_img_like(T1_data, bufferT1)
    template_data = image.load_img(template)
    html_view = plotting.view_img(template_data, bg_img=cappedT1,
                                  cmap=cmap_affine,
                                  symmetric_cmap=False,
                                  colorbar=False,
                                  opacity=1.0)
    html_view.save_as_html(save_path)

def viewer_final(T1, labeling, save_path):
    T1_data = image.load_img(T1)
    buffer = T1_data.get_fdata()
    upper95 = np.quantile(buffer, 0.95)
    lower5 = np.quantile(buffer, 0.05)
    buffer[buffer > upper95] = upper95
    buffer[buffer < lower5] = lower5
    cappedT1 = image.new_img_like(T1_data, buffer)
    tissues_data = image.load_img(labeling)
    tissues_buffer = tissues_data.get_fdata()
    labels = np.unique(tissues_buffer)
    labels = labels[1:]
    for l in labels:
        tissue_mask = (tissues_buffer == l)
        ero = mrph.binary_erosion(np.squeeze(tissue_mask), se, 1)
        tissues_buffer[ero] = 0

    tissue_borders = image.new_img_like(tissues_data, tissues_buffer)

    html_view = plotting.view_img(tissue_borders, bg_img=cappedT1, opacity=0.9, cmap=cmap_final, symmetric_cmap=False,
                                  dim=2, colorbar=False)

    html_view.save_as_html(save_path)

def _mosaic_images(dim, num_of_blocks=5):

    patternX = np.sin(np.arange(dim[0])/dim[0] * 2 * np.pi * (num_of_blocks / 2))
    patternY = np.sin(np.arange(dim[1])/dim[1] * 2 * np.pi * (num_of_blocks / 2))
    patternZ = np.sin(np.arange(dim[2])/dim[2] * 2 * np.pi * (num_of_blocks / 2))
    patternXY = np.outer(patternX, patternY)
    mask = np.zeros(dim)

    for d in range(dim[2]):
        mask[:,:,d] = patternXY * patternZ[d]

    return mask

