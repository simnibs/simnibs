import numpy as np
import nibabel as nib

from ..utils.file_finder import SubjectFiles
from ..utils.transformations import mni2subject_coords, create_new_connectivity_list_point_mask

class RegionOfInterest():
    """
    Helper class to catch ROI parameters before initialization
    """

    def __init__(self):
        """
        Initializes RegionOfInterest class instance
        """
        self.mesh = None
        self.center = None
        self.domains = None
        self.type = "custom"                    # "custom", "GMmidlayer" or "volume"
        self.roi_sphere_center_mni = None
        self.roi_sphere_center_subject = None
        self.roi_sphere_radius = None
        self.ff_subject = None

    def initialize(self):
        """
        Run ROI initialization and return final RegionOfInterest class instance

        Returns
        -------
        roi : RegionOfInterest class instance
            Region of Interest
        """

        if self.center is not None and type(self.center is list):
            self.center = np.array(self.center)

        if self.nodes is not None and type(self.nodes is list):
            self.nodes = np.array(self.nodes)

        if self.con is not None and type(self.con is list):
            self.con = np.array(self.con)

        if self.roi_sphere_center_mni is not None and type(self.roi_sphere_center_mni is list):
            self.roi_sphere_center_mni = np.array(self.roi_sphere_center_mni)

        if self.roi_sphere_center_subject is not None and type(self.roi_sphere_center_subject is list):
            self.roi_sphere_center_subject = np.array(self.roi_sphere_center_subject)

        self.ff_subject = SubjectFiles(fnamehead=self.mesh.fn)

        if self.type == "custom":
            pass

        elif self.type == "GMmidlayer":
            if ((self.roi_sphere_center_mni is None and self.roi_sphere_center_subject is None)
                    or self.roi_sphere_radius is None):
                raise AssertionError("Please specify roi_sphere_center_mni or roi_sphere_center_subject and "
                                     "roi_sphere_radius for 'MNI' ROI type.")

            # transform MNI to subject coordinates
            if self.roi_sphere_center_subject is None:
                self.roi_sphere_center_subject = mni2subject_coords(self.roi_sphere_center_mni, self.ff_subject.subpath)

            # load midlayer from m2m_* folder
            img_lh = nib.load(self.ff_subject.surfaces["central"]["lh"])
            img_rh = nib.load(self.ff_subject.surfaces["central"]["rh"])

            lh_nodes = img_lh.darrays[0].data
            rh_nodes = img_rh.darrays[0].data

            lh_con = img_lh.darrays[1].data
            rh_con = img_rh.darrays[1].data

            # merge two hemispheres
            lh_rh_nodes = np.vstack((lh_nodes, rh_nodes))
            lh_rh_con = np.vstack((lh_con, rh_con+lh_nodes.shape[0]))

            # mask out ROI sphere
            mask = np.linalg.norm(lh_rh_nodes - self.roi_sphere_center_subject, axis=1) <= self.roi_sphere_radius
            self.nodes, self.con = create_new_connectivity_list_point_mask(points=lh_rh_nodes,
                                                                           con=lh_rh_con,
                                                                           point_mask=mask)
            self.center = np.mean(self.nodes[self.con,], axis=1)

        elif self.type == "volume":
            if ((self.roi_sphere_center_mni is None and self.roi_sphere_center_subject is None)
                    and self.roi_sphere_radius is None):
                self.roi_sphere_center_mni = np.array([0, 0, 0])
                self.roi_sphere_radius = 1e6

            if ((self.roi_sphere_center_mni is None and self.roi_sphere_center_subject is None)
                    or self.roi_sphere_radius is None):
                raise AssertionError("Please specify roi_sphere_center or roi_sphere_center_subject and "
                                     "roi_sphere_radius for 'MNI' ROI type.")

            # transform MNI to subject coordinates
            if self.roi_sphere_center_subject is None:
                self.roi_sphere_center_subject = mni2subject_coords(self.roi_sphere_center_mni, self.ff_subject.subpath)

            # mask out spherical ROI (mask targets whole mesh.elm.node_number_list)
            ele_center = np.mean(self.mesh.nodes.node_coord[self.mesh.elm.node_number_list-1, ], axis=1)
            ele_center[self.mesh.elm.elm_type != 4] = np.inf

            if self.domains is not None:
                domains_remove = [d for d in np.unique(self.mesh.elm.tag1) if d not in self.domains]

                for d in domains_remove:
                    ele_center[self.mesh.elm.tag1 == d] = np.inf

            mask = np.linalg.norm(ele_center - self.roi_sphere_center_subject, axis=1) <= self.roi_sphere_radius
            mask_cropped = mask[self.mesh.elm.elm_type == 4]

            self.center = self.mesh.elements_baricenters()[:][mask_cropped]
        else:
            raise AssertionError("Specified ROI type not implemented ('custom', 'MNI', 'GMmidlayer')")

        return self.center


