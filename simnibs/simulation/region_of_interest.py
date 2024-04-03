import numpy as np
import nibabel as nib

from simnibs.mesh_tools.mesh_io import Elements, Msh, Nodes

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

        if self.roi_sphere_center_mni is not None and type(self.roi_sphere_center_mni is list):
            self.roi_sphere_center_mni = np.array(self.roi_sphere_center_mni)

        if self.roi_sphere_center_subject is not None and type(self.roi_sphere_center_subject is list):
            self.roi_sphere_center_subject = np.array(self.roi_sphere_center_subject)

        self.ff_subject = SubjectFiles(fnamehead=self.mesh.fn)

        if self.type == "custom":
            self._cropped_mesh = None

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
            nodes, con = create_new_connectivity_list_point_mask(points=lh_rh_nodes,
                                                                           con=lh_rh_con,
                                                                           point_mask=mask)
            self._cropped_mesh = Msh(Nodes(nodes), Elements(triangles=con+1))
            self.center = self._cropped_mesh.elements_baricenters()[:]

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

            # mask out spherical ROI 
            ele_center = self.mesh.elements_baricenters()[:]

            mask = self.mesh.elm.elm_type == 4

            if self.domains is not None:
                domain_mask = np.ones_like(mask)
                for d in self.domains:
                    domain_mask = domain_mask | self.mesh.elm.tag1 == d
                mask = mask & domain_mask

            distance_mask = np.linalg.norm(ele_center[mask] - self.roi_sphere_center_subject, axis=1) <= self.roi_sphere_radius
            self.center = ele_center[mask][distance_mask]
            self._cropped_mesh = self.mesh.crop_mesh(elements= self.mesh.elm.elm_number[mask][distance_mask])
        else:
            raise AssertionError("Specified ROI type not implemented ('custom', 'MNI', 'GMmidlayer')")

        return self.center


