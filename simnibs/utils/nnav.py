# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 20:37:50 2022

@author: axthi
"""
import io
import math
import os.path
import re
import warnings
import numpy as np
import nibabel as nib
import json
from datetime import datetime
import xml.etree.ElementTree as ET

from ..simulation import TMSLIST, POSITION
from .. import __version__ as simnibs_v
from .file_finder import SubjectFiles


def _lps2ras():
    """
    Flip matrix from LPS to RAS coordinate system and vice versa.

    Example:
        if mat_in_lps.shape is (n_pos, 4, 4)
        mat_in_ras = _lps2ras() @ mat_in_lps
        mat_in_lps = _lps2ras() @ map_in_ras

        if mat_in_lps.shape is (4, 4, n_pos):
         mat_in_ras = np.tensordot(_lps2ras(), mat_in_lps, axes=[0,0])
    """
    return np.array([
        [-1, +0, +0, +0],  # +x -> -x (R -> L)
        [+0, -1, +0, +0],  # +y -> -y (A -> P)
        [+0, +0, +1, +0],  # z -> z
        [0, 0, 0, 1]])


class localite:
    """
    I/O of localite data
    """

    def read(self, fn, markertype=None):
        """
        Imports coil positions/orientations from a Localite TMS Neuronavigator .xml file as TMSLIST().
    
        File(s) can be of type TriggerMarkers_Coil*_%date%.xml or InstrumentMarker%date%.xml.
        Stimulation intensity (didt) for each TMSLIST.pos is read from .xml file or defaults to 1Aµs.
    
        Written by Ole Numssen, numssen@cbs.mpg.de, 2022.
    
        Parameter:
        ----------
        fn: str or list of str
            Filename(s) of InstrumentMarker or TriggerMarker .xml files.
        markertype: str, optional
            Can be one of ('InstrumentMarker','TriggerMarker', None). If None, Markertype is guessd from filename.
    
        Returns:
        -------
        tmslist : simnibs.simulation.sim_struct.TMSLIST or list of simnibs.simulation.sim_struct.TMSLIST
        """
        if isinstance(fn, list):
            tms_lists = []
            for f in fn:
                tms_lists.append(self.read(f, markertype))
            return tms_lists

        assert os.path.exists(fn), f"File does not exist: {fn}"
        timestamp = ''
        empty_pos = np.array(
                [[0, 0, 1, 0],
                 [0, -1, 0, 0],
                 [1, 0, 0, 0],
                 [0, 0, 0, 1]])  # 'coil position' for untracked coils.

        if not markertype:
            # guess markertype from filename
            if 'InstrumentMarker' in fn:
                markertype = 'InstrumentMarker'

                # get date from xml filename
                # timestamp = regex.findall(r"InstrumentMarker(\d*).xml", fn)[0]
                timestamp = re.findall(r"InstrumentMarker(\d*).xml", fn)[0]
                timestamp = datetime.strptime(timestamp, "%Y%m%d%H%M%f")
                timestamp = timestamp.strftime("%Y-%m-%d %H:%M:%S")
            elif 'TriggerMarker' in fn:
                markertype = 'TriggerMarker'

                # get date from xml filename
                # timestamp = regex.findall(r"TriggerMarkers_Coil\d_(\d*).xml", fn)[0]
                timestamp = re.findall(r"TriggerMarkers_Coil\d_(\d*).xml", fn)[0]
                timestamp = datetime.strptime(timestamp, "%Y%m%d%H%M%f")
                timestamp = timestamp.strftime("%Y-%m-%d %H:%M:%S")
            else:
                raise ValueError(f"'markertype' not provided and cannot guess from filename ({fn}).")

        # build TMSLIST
        tms_list = TMSLIST()
        for i, (im, descr, onset, stim_int, didt) in \
                enumerate(zip(*self._parse_localite_xml(fn, markertype=markertype))):

            if np.all(im == empty_pos):
                print(f"Skipping untracked TriggerMarker #{i:0>3} for {fn}.")
                continue

            p = POSITION()
            p.matsimnibs = im
            p.name = descr
            if onset:
                p.date = f"{timestamp}:{onset:0>7}"
            else:
                p.date = timestamp
            p.didt = didt

            tms_list.add_position(p)
        return tms_list

    def write(self, matsimnibs, xml_fn, names=None, overwrite=False, out_coord_space='RAS'):
        """
        Writes an instrument marker .xml file in the fashion of Localite TMS Navigator.
        Input can be a single or multiple 4x4 matsimnibs matrices with coil position and orientation.
    
        Parameters
        ----------
        matsimnibs: np.ndarray or TMSLIST or POSITION
            Coil position/orientation(s) to export.
            If np.ndarray, shape is (4, 4, n_positions).
        xml_fn : str
            Output filename.
        names : str or list of str
            Optional, name for each instrument markers
        overwrite : bool (Default: False)
            Overwrite existing file.
        out_coord_space : str, one of ['RAS','LPS']. Default: 'RAS'
            Coordinate space of the T1 used for neuronavigation.
            Rule of thumb:
                DICOM -> 'LPS'
                NIFTI -> 'RAS'
    
        Written by Ole Numssen, numssen@cbs.mpg.de, 2022.
    
        Returns
        -------
        file: fn_out.xml
        """
        # unpack tmslist/position into np.ndarray
        if isinstance(matsimnibs, TMSLIST):
            matsimnibs = np.stack([p.matsimnibs for p in matsimnibs.pos], axis=2)
        elif isinstance(matsimnibs, POSITION):
            matsimnibs = matsimnibs.matsimnibs
        matsimnibs = np.atleast_3d(matsimnibs)

        if isinstance(names, str):
            names = [names]

        # check inputs
        assert matsimnibs.shape[:2] == (4, 4), 'Expecting array with shape (4, 4, N instrument marker).'
        out_coord_space = out_coord_space.upper()
        assert out_coord_space in ['RAS', 'LPS'], f'out_coord_space={out_coord_space} is not one of ["RAS","LPS"].'

        if not xml_fn.lower().endswith('.xml'):
            xml_fn += '.xml'

        assert not os.path.exists(xml_fn) or overwrite, 'File {fn} already exists. Remove or set overwrite=True.'

        with io.open(xml_fn, 'w', newline='\n') as f:  # correct windows style would be \r\n, but Localite uses \n
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            f.write(f'<InstrumentMarkerList coordinateSpace="{out_coord_space}">\n')
            f.write(f'    <!--This InstrumentMarkerList was written by SimNIBS v{simnibs_v}-->\n')
            n_ims = matsimnibs.shape[-1]
            for idx in range(n_ims):
                if names is None:
                    name = f"{idx:0>{int(math.log10(n_ims)) + 1}}"
                else:
                    name = names[idx]

                im = matsimnibs[:, :, idx]

                # bring into Localite axes definition
                im = im @ self._simnibs2localite()

                # transform into correct coordinate axes system
                if out_coord_space == 'LPS':
                    im = _lps2ras() @ im

                f.write('    ' + f'<InstrumentMarker alwaysVisible="false" index="{idx}" selected="false">\n')
                f.write('    ' * 2 + f'<Marker additionalInformation="" '
                                     f'color="#ff0000" description="{name}" set="true">\n')
                f.write('    ' * 3 + '<Matrix4D \n')
                f.write('    ' * 4 + 'data00="{:+1.17f}" data01="{:+1.17f}" '
                                     'data02="{:+1.17f}" data03="{:+1.17f}"\n'.format(im[0, 0], im[0, 1], im[0, 2],
                                                                                      im[0, 3]))
                f.write('    ' * 4 + 'data10="{:+1.17f}" data11="{:+1.17f}" '
                                     'data12="{:+1.17f}" data13="{:+1.17f}"\n'.format(im[1, 0], im[1, 1], im[1, 2],
                                                                                      im[1, 3]))
                f.write('    ' * 4 + 'data20="{:+1.17f}" data21="{:+1.17f}" '
                                     'data22="{:+1.17f}" data23="{:+1.17f}"\n'.format(im[2, 0], im[2, 1], im[2, 2],
                                                                                      im[2, 3]))
                f.write('    ' * 4 + 'data30="{:+1.17f}" data31="{:+1.17f}" '
                                     'data32="{:+1.17f}" data33="{:+1.17f}"/>\n'.format(0, 0, 0, 1))
                f.write('    ' * 2 + '</Marker>\n')
                f.write('    ' + '</InstrumentMarker>\n')

            f.write('</InstrumentMarkerList>\n')

    def _parse_localite_xml(self, im_path, markertype):
        """
            Read coil position/orientation Localite TMS Neuronavigator .xml-file.
    
            Parameters
            ----------
            im_path : str or list of str
                Path to instrument-marker-file
            markertype : str
                'InstrumentMarker', or 'TriggerMarker'
    
            Returns
            -------
            marker_arrays : np.ndarray of float [Mx4x4]
                Instrument-marker-matrices in simnibs space
            marker_descr : list of (str, None)
                Labels of the instrument-marker-descriptions
            marker_times : list of (str, None)
                Pulse onset times for TriggerMarker in ms. 0 is start of recording
            stim_intens : list of (float, None)
                Stimulator intensities for TriggerMarker in %MSO
            didt_intens : list of float
                Realized stimulation intensitties for TriggerMarkers in A/s. Defaults to 1e6 for InstrumentMarkers
            """
        marker_arrays = np.empty([0, 4, 4])
        marker_descr = []
        didt_intens, stim_intens, marker_times = [], [], []  # for TriggerMarkers

        # parse XML document
        im_tree = ET.parse(im_path)
        im_root = im_tree.getroot()
        coord_system = im_root.get('coordinateSpace').upper()
        if coord_system == 'RAS':
            lps2ras_flipmat = np.eye(4)
        elif coord_system == 'LPS':
            print("Transforming 'LPS' to 'RAS' coordinate system.")
            lps2ras_flipmat = _lps2ras()
        else:
            raise ValueError(f"Coordinate system {coord_system} not supported. Use 'RAS' or 'LPS'.")
        coil_axes_flipmat = self._simnibs2localite()

        # iterate over all 'InstrumentMarker' tags
        for marker_i in im_root.iter(markertype):
            marker_arr = np.empty([1, 4, 4])
            # get tag were the matrix is
            if markertype == 'InstrumentMarker':
                marker_object = marker_i.find('Marker')
                marker_descr.append(marker_object.get('description'))
                matrix4d = marker_object.find('Matrix4D')

                # also fill TriggerMarker arrays
                marker_times.append(None)
                stim_intens.append(None)
                didt_intens.append(1e6)

            elif markertype == 'TriggerMarker':
                matrix4d = marker_i.find('Matrix4D')
                marker_descr.append(marker_i.get('description'))
                marker_times.append(marker_i.get('recordingTime'))

                # read di/dt and stimulator intensity
                im_rv = marker_i.find('ResponseValues').findall('Value')

                for _im_rv in im_rv:
                    # di/dt
                    if _im_rv.get('key') == "valueA":
                        val = _im_rv.get('response')
                        if val in ['NaN', "0.0"]:
                            warnings.warn("Couldn't determine stimulation intensity. Defaulting to 1 A/µs")
                            val = 1e6
                        else:
                            val = float(val) * 1e6
                        didt_intens.append(val)

                    # stimulator intensity
                    elif _im_rv.get('key') == "amplitudeA":
                        val = _im_rv.get('response')
                        try:
                            val = float(val)
                        except ValueError:
                            pass
                        stim_intens.append(float(val))
            else:
                raise ValueError(f"markertype={markertype} unknown.")

            # get position and orientation values
            for im_index1 in range(4):
                for im_index2 in range(4):
                    marker_arr[0, im_index1, im_index2] = (float(matrix4d.get(f"data{str(im_index1)}{str(im_index2)}")))

            # transform to simnibs space
            marker_arr = lps2ras_flipmat @ marker_arr
            marker_arr = marker_arr @ coil_axes_flipmat
            marker_arrays = np.append(marker_arrays, marker_arr, axis=0)

        return marker_arrays, marker_descr, marker_times, stim_intens, didt_intens

    def _simnibs2localite(self):
        """
        Flip matrix for localite -> simnibs and vice versa.
    
        Example:
            simnibs_mat = localite_mat @ localite._simnibs2localite()
            localite_mat = simnibs_mat @ localite._simnibs2localite()
    
        """
        return np.array([
            [+0, +0, +1, +0],  # +z -> +x
            [+0, -1, +0, +0],  # +y -> -y
            [+1, +0, +0, +0],  # +z -> +x
            [0, 0, 0, 1]])


class softaxic:
    """
    I/O of softaxic data
    """

    def read(self, filename):
        """
        Imports coil positions/orientations from a softaxic Neuronavigator as TMSLIST().
        
        Written by Kris H. Madsen, khma@dtu.dk, 2022.
    
        Parameter:
        ----------
        fn: str 
            Filename of softaxic .stmpx file
    
        Returns:
        -------
        tmslist : simnibs.simulation.sim_struct.TMSLIST 

        """
        tmslist = TMSLIST()
        M_lst, pos_id_lst = self._parse_softaxic(filename)
        for M, pos_id in zip(M_lst, pos_id_lst):
            position = tmslist.add_position()
            position.matsimnibs = M
            position.name = pos_id
        return tmslist

    def _parse_softaxic(self, fn):
        root = ET.parse(fn).getroot()  # parse the XML-like file
        coords = ('x', 'y', 'z')  # strings representing the three coordinates
        # Rotation matrix for changing into coil coordinate system in softaxis convention
        # rotation is such that x and y are interchanged
        R = np.zeros((3, 3))
        R[0, 1] = 1
        R[1, 0] = -1
        R[2, 2] = 1
        # Placeholder for final matsimnibs affine position definitions
        M = []
        # Placeholder for direction cosine matrix (not really used except inside the loop)
        D = []
        # Placeholder for positions (not really used except inside the loop)
        pos = []
        # Placeholder for list of position ids
        pos_id = []
        # Iterate over all fmp elements
        for k, elem in enumerate(root.iter('fmp')):
            pos_id.append(elem.attrib['id'])
            attrib = elem.find('fp').attrib  # extract fp attributes
            D.append(np.zeros((3, 3)))  # init D
            pos.append(np.zeros(3))  # init pos
            for i in range(3):  # loop over first dimension af M
                for j in range(3):  # loop over second dimension of M
                    D[k][i, j] = attrib[f'm{i}{j}']  # set D elements
                pos[k][i] = attrib[f'{coords[i]}']  # set positions
            M.append(np.identity(4))  # init M
            M[k][:3, :3] = D[k] @ R  # rotate DCT part
            M[k][:3, 3] = pos[k]  # set center position
        return M, pos_id


class brainsight:
    """
    I/O of brainsight data
    """

    def read(self, fn):
        """
        Import coil positions/orientations from a BrainSight neuronavigation system.
    
        Data has to be exported from Brainsight to 'NIfTI:Aligned' space.
        One TMSLIST for 'Targets' and one for 'Samples' is returned.

        We expect:
            - the same T1 scan for nnav and SimNIBS is used
            - qform and sform to be equal (automatically set by SimNIBS now)

        Notes:
        ------
        Sean: "Brainsight supports multiple coordinate systems when exporting to its .txt file format,
        and the possible options depend on the project's anatomical dataset, specifically:

        * "Brainsight" - this option is always available. It's Brainsight's internal coordinate system,
            where the X axis goes from right to left, Y axis goes from anterior to posterior, and the Z axis goes
            from inferior to superior. The origin is the centre of the voxel at the right-anterior-inferior corner.
            The units are millimetres. SimNIBS calls this `LPS`, short for the direction each axis points towards.
        * "World" - this option is only available for projects *not* based on a NIfTI file.
            Since SimNIBS 4 no longer supports anything but NIfTI, this choice is not applicable.
        * "NIfTI:Aligned" - this option is only available for projects based on a NIfTI file.
            It corresponds to the qform/sform `NIFTI_XFORM_ALIGNED_ANAT
            <https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html>`_
            coordinate system.
        * "NIfTI:Scanner" - same but for `NIFTI_XFORM_SCANNER_ANAT
            <https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html>`_
        * "NIfTI:MNI-152" - same but for `NIFTI_XFORM_MNI_152
            <https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html>`_
        * "NIfTI:Talairach" - same but for `NIFTI_XFORM_TALAIRACH
            <https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html>`_

        Many NIfTI files have only an sform, or only a qform, or they are both the same.
        In that case, only one of the `NIfTI:*` options above will be available, and that is what you should use.

        If the sform and qform are different, each will result in a `NIfTI:*` option,
        but the Brainsight user interface does not indicate which is from the `qform` and which is from the `sform`.
        To deterime that, you'd need to use a tool like
        `nifti1_tool <https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/programs/nifti1_tool_sphx.html>`_
        and look at the `qform_code` and `sform_code`."

        Parameter:
        ----------
        fn: str
            Filename to Brainsight .txt export
    
        Returns:
        --------
        results : [TMSLIST, TMSLIST]
            [Targets, Samples] as SimNIBS TMSLIST objects
    
        Written by Ole Numssen, numssen@cbs.mpg.de; Konstantin Weise, kweise@cbs.mpg.de; 2023.
        """
        # init empty values in case nothing is found in fn
        coord_sys, encoding = '', ''
        data_targets, data_samples = [], []

        with open(fn, 'r') as f:
            # read header
            while True:
                line = f.readline().rstrip()
                if line.startswith('# Target Name') or line.startswith('# Sample Name'):
                    break
                elif line.startswith('# Coordinate system:'):
                    coord_sys = line.replace("# Coordinate system: ", "")
                elif line.startswith('# Encoding: '):
                    encoding = line.replace("# Encoding: ", "")

            if coord_sys.lower() != 'nifti:aligned':
                raise ValueError(f"Coordinate system '{coord_sys}' is not supported. "
                                 f"Export targes/samples as NIfTI:Aligned.")

            # Let's only read UTF-8
            if encoding == '':
                warnings.warn(f"Cannot read encoding from {fn}. Assuming UTF-8.")
            elif encoding != 'UTF-8':
                raise ValueError(f"Encoding '{encoding}' not supported. Use UTF-8.")

            # read data, probably 'Targets'
            if line.startswith('# Target Name'):
                # get column names
                col_names_targets = line.replace('# ', '').split('\t')

                # read all target lines
                while line:
                    line = f.readline()
                    if line.startswith('#'):
                        break
                    if line:
                        data_targets.append(line.rstrip().split('\t'))

            if line.startswith('# Sample Name'):
                # get column names
                col_names_samples = line.replace('# ', '').split('\t')

                # read all sample lines
                while line:
                    line = f.readline()
                    if line.startswith('#'):
                        break
                    data_samples.append(line.rstrip().split('\t'))

            # get matsimnibs arrays in simnibs space and axes definition
            if len(data_targets):
                names_targets, matsimnibs_targets = self._transform_brainsight(data_targets,
                                                                               col_names_targets)
            else:
                names_targets, matsimnibs_targets = [], None

            if len(data_samples):
                names_samples, matsimnibs_samples = self._transform_brainsight(data_samples,
                                                                               col_names_samples)
            else:
                names_samples, matsimnibs_samples = [], None

            if (matsimnibs_targets is None or matsimnibs_targets.size == 0) and \
                    (matsimnibs_samples is None or matsimnibs_samples.size == 0):
                raise ValueError(f"Could not find any targets in {fn}.")

            # get TMSLIST for targets
            tms_list_targets = TMSLIST()
            tms_list_targets.name = 'Targets'
            for i, name in enumerate(names_targets):
                p = POSITION()
                p.matsimnibs = matsimnibs_targets[i]
                p.name = name
                tms_list_targets.add_position(p)

            # get TMSLIST for samples
            tms_list_samples = TMSLIST()
            tms_list_samples.name = 'Samples'
            for i, name in enumerate(names_samples):
                p = POSITION()
                p.matsimnibs = matsimnibs_samples[i]
                p.name = name
                tms_list_samples.add_position(p)

            # return either a single TMSLIST or a list of TMSLIST if samples and targets are found.
            return [tms_list_targets, tms_list_samples]

    def write(self, matsimnibs, fn, names=None, overwrite=False):
        """
        Writes an .txt file that can be imported with the Brainsight neuronavition system.
    
        Input can be a single or multiple 4x4 matsimnibs matrices with coil position and orientation.
        Output space is NIfTI:aligned.
    
        matsimnibs: np.ndarray or TMSLIST or POSITION
            Coil position/orientation(s) to export.
            If np.ndarray, shape is (4, 4, n_positions).
        names : str or list of str
            Optional, name for each markers
        fn : str
            Output filename.
        overwrite : bool (Default: False)
            Overwrite existing file.

        Written by Ole Numssen, numssen@cbs.mpg.de, 2023.
        """
        # unpack tmslist/position into np.ndarray
        if isinstance(matsimnibs, TMSLIST):
            matsimnibs = np.stack([p.matsimnibs for p in matsimnibs.pos], axis=2)
        elif isinstance(matsimnibs, POSITION):
            matsimnibs = matsimnibs.matsimnibs
        matsimnibs = np.atleast_3d(matsimnibs)

        if isinstance(names, str):
            names = [names]

        # check inputs
        assert matsimnibs.shape[:2] == (4, 4), 'Expecting array with shape (4, 4, N instrument marker).'

        out_coord_space = 'NIfTI:Aligned'

        # change coil axes definition to brainsight
        matsimnibs = np.matmul(self._simnibs2brainsight(), matsimnibs)

        if not fn.lower().endswith('.txt'):
            fn += '.txt'

        assert not os.path.exists(fn) or overwrite, f'File {fn} already exists. Remove or set overwrite=True.'

        with open(fn, 'w') as f:  # correct windows style would be \r\n, but Localite uses \n
            f.write('# Version: 12\n')
            f.write(f'# Coordinate system: {out_coord_space}\n')
            f.write(f'# Created by: SimNIBS v{simnibs_v}\n')
            f.write('# Units: millimetres, degrees, milliseconds, and microvolts\n')
            f.write('# Encoding: UTF-8\n')
            f.write(
                '# Notes: Each column is delimited by a tab. Each value within a column is delimited by a semicolon.\n')
            f.write('# Target Name	'
                    'Loc. X	Loc. Y	Loc. Z	'
                    'm0n0	m0n1	m0n2	'
                    'm1n0	m1n1	m1n2	'
                    'm2n0	m2n1	m2n2\n')
            for i in range(matsimnibs.shape[2]):
                name = names[i] if names is not None else i
                f.write(f'{name:0>3}\t' +
                        f'{matsimnibs[0, 3, i]:.4f}\t{matsimnibs[1, 3, i]:.4f}\t{matsimnibs[2, 3, i]:.4f}\t' +
                        f'{matsimnibs[0, 0, i]:.4f}\t{matsimnibs[1, 0, i]:.4f}\t{matsimnibs[2, 0, i]:.4f}\t' +
                        f'{matsimnibs[0, 1, i]:.4f}\t{matsimnibs[1, 1, i]:.4f}\t{matsimnibs[2, 1, i]:.4f}\t' +
                        f'{matsimnibs[0, 2, i]:.4f}\t{matsimnibs[1, 2, i]:.4f}\t{matsimnibs[2, 2, i]:.4f}\n')

    def _transform_brainsight(self, data, col_names):
        """
        Transforms Brainsight coil position/orientation into SimNIBS matsimnibs
    
        data: list of lists with positions
        col_names: list of str
    
        Returns:
        --------
        names : list of str
            Samples/target names. len(names) = n_pos
        matsimnibs: np.ndarray
            Coil position/orietnation in SimNIBS style, coil axes definition and RAS space.
            shape: (n_pos, 4, 4)
    
        Written by Ole Numssen, numssen@cbs.mpg.de, 2022.
        """
        matsimnibs = np.zeros((len(data), 4, 4))
        pos_names = []
        for pos, i in zip(data, range(len(data))):
            p_dict = {x: pos[k] for x, k in zip(col_names, range(len(col_names)))}

            m = [[p_dict['m0n0'], p_dict['m1n0'], p_dict['m2n0'], p_dict['Loc. X']],
                 [p_dict['m0n1'], p_dict['m1n1'], p_dict['m2n1'], p_dict['Loc. Y']],
                 [p_dict['m0n2'], p_dict['m1n2'], p_dict['m2n2'], p_dict['Loc. Z']],
                 [0, 0, 0, 1]]

            matsimnibs[i] = np.array(m).astype(float)
            pos_names.append(pos[0])

        # adjust coil axes definition to simnibs style
        return pos_names, matsimnibs @ self._simnibs2brainsight()

    def _simnibs2brainsight(self):
        """
        Flip matrix for brainsight -> simnibs and vice versa.
    
        Example:
            if brainsight_mat.shape = (n_pos, 4, 4):
            simnibs_mat = brainsight_mat @ simnibs2localite()
            localite_mat = simnibs_mat @ brainsight2localite()
    
        """
        return np.array([
            [-1, +0, +0, +0],  # +x -> -x
            [+0, +1, +0, +0],  # +y -> +y
            [+0, +0, -1, +0],  # +z -> -z
            [0, 0, 0, 1]])


class ant:
    """
    I/O of ANT data
    """
    
    def _compare_imageinfos(self, ImageInfo, ImageInfo2):
        ''' compare two image info dicts'''
        keys1 = set(ImageInfo.keys())
        keys2 = set(ImageInfo2.keys())
        
        shared_keys = keys1.intersection(keys2)
        if len(shared_keys) < len(keys1):
            warnings.warn("image infos have different fields")
            return False
        
        for k in keys1:
            if len(ImageInfo[k]) != len(ImageInfo2[k]):
                warnings.warn("these info field differs in length: "+str(k))
                return False
            else:
                if not np.all(np.isclose(ImageInfo[k], ImageInfo2[k], atol=1e-05)):
                    warnings.warn('these info field has varying entries: '
                                  +str(k)+' '+str(ImageInfo[k])+' '+str(ImageInfo2[k])
                                  )
                    return False    
        return True
    
    def _imageinfo_from_nifti(self, fn):
        ''' loads the relevant infos from image header into a dict'''
        hdr = nib.load(fn).header
        ImageInfo = {}
        ImageInfo['NIfTI_quatern'] = [float(hdr['quatern_b']),
                                      float(hdr['quatern_c']),
                                      float(hdr['quatern_d'])]
                                                
        ImageInfo['NIfTI_pixdim'] = np.asarray(hdr['pixdim'][:4], dtype=float).tolist()
        ImageInfo['NIfTI_qoffset'] = [float(hdr['qoffset_x']),
                                      float(hdr['qoffset_y']),
                                      float(hdr['qoffset_z'])]
        return ImageInfo
    
    def _mrk_template(self):
        mrk_template = {'Configuration': 
                        {'SourceInfo': {'AppName': 'SimNIBS',
                                        'AppVersion': str(simnibs_v)
                                        },
                         'AssociatedImageInfo': {'NIfTI_quatern': [0.0, 0.0, 0.0],
                                                 'NIfTI_qoffset': [0.0, 0.0, 0.0],
                                                 'NIfTI_pixdim':  [0.0, 0.0, 0.0, 0.0]
                                                 }
                         },
                          'TargetMarkers': []
                        }
        return mrk_template
    
    def _target_template(self):
        target_template = {'CoilTarget': 
                           {'Name': '',
                            'CoordinateSystem': 'ASSOCIATED_IMAGE',
                            'PositionalUnit': 'mm',
                            'CoilAxesConvention': 'simnibs',
                            'CoilInfo': {'Type': '',
                                         'Id': ''
                                         },
                            'CoilPose': {'RigidTransformMatrix4x4': []}
                            }
                           }
        return target_template
    
    def read(self, fn, subpath=None, imageinfo=None, return_imageinfo=False):
        """
        Imports coil positions from a .mrk-file written by an ANT 
        neuronavigation system. 
        
        Optionally, the image information stored in the .mrk-file will be compared
        to the T1 scan of the SimNIBS headmodel to ensure that the coordinates 
        spaces match.
        
        Parameters
        ----------
        fn : str
            Filename of .mrk-file.
        subpath : str, optional
            Path to m2m-folder of the subject. If given, some header fields of 
            m2m-{subID}/T1.nii.gz (qform quaternion, qform offset, pixdim) will 
            be compared to the image information stored in the .mrk-file.
            Mismatches raise warnings. The default is None.
        imageinfo : dict, optional
            Image information to compare to the information stored in the .mrk-file.
            Can be used as alternative to providing subpath. subpath takes precedence
            in case both are set. The default is None.
        return_imageinfo : bool, optional
            Returns the image information stored in the .mrk-file. The default is False.
    
        Returns
        -------
        list of TMSLIST
            The list contains one TMSLIST for each coil ID found in the .mrk-file
        dict (optional)
            The dict contains the image info stored in .mrk-file
        
        Written by Axel Thielscher, 2023
        """
        assert os.path.exists(fn), f"File does not exist: {fn}"
        
        with open(fn) as marker_file:
          marker_str = marker_file.read()
        marker_dict = json.loads(marker_str)
        
        # determine image info from .mrk, and (optionally) compare to internal info
        try:
            conf = marker_dict.get('Configuration')
            imInfo_mrk = conf.get('AssociatedImageInfo')
        except:
            imInfo_mrk = None
            warnings.warn("no image information associated with positions")
        
        imInfo_internal = None
        if imageinfo is not None:
            imInfo_internal = imageinfo
        if subpath is not None:
            ff = SubjectFiles(subpath=subpath)
            imInfo_internal = self._imageinfo_from_nifti(ff.reference_volume)
            
        if imInfo_mrk is not None and imInfo_internal is not None:
            if not self._compare_imageinfos(imInfo_mrk, imInfo_internal):
                warnings.warn("image informations DO NOT MATCH!")       
    
        # read in positions
        tms_lists = []
        targets=marker_dict.get('TargetMarkers')
        if targets is None:
            warnings.warn("marker file does not contain positions")
        else:
            for t in targets:
                coiltarget = t.get('CoilTarget')
                if coiltarget is None:
                    warnings.warn("position is not a CoilTarget, skipping")
                else:
                    assert coiltarget.get('CoordinateSystem') == 'ASSOCIATED_IMAGE'
                    assert coiltarget.get('PositionalUnit') == 'mm'
                    assert coiltarget.get('CoilAxesConvention') == 'simnibs'
                    assert 'CoilPose' in coiltarget
                    assert 'RigidTransformMatrix4x4' in coiltarget['CoilPose']
                    
                    p = POSITION()
                    if 'Name' in coiltarget:
                        p.name = coiltarget['Name']
                    p.matsimnibs = np.asarray(coiltarget['CoilPose']['RigidTransformMatrix4x4'])
                    
                    # get coil filename
                    fn_coil = ''
                    if 'CoilInfo' in coiltarget:
                        if 'Id' in coiltarget['CoilInfo']:
                            fn_coil = coiltarget['CoilInfo']['Id']
                            
                    # get tmslist to which position should be added (match coil filename)
                    current_list = None
                    for tms_list in tms_lists:
                        if tms_list.fnamecoil == fn_coil:
                            current_list = tms_list
                            break
                    
                    # no matching tmslist found, or tms_lists still empty
                    if current_list is None:
                        tms_lists.append(TMSLIST())
                        current_list = tms_lists[-1]
                        current_list.fnamecoil = fn_coil
                    
                    current_list.add_position(p)
                    
        if return_imageinfo:
            return tms_lists, imInfo_mrk
        return tms_lists
                        
    def write(self, tms_lists, fn, subpath=None, imageinfo=None):
        """
        Exports coil positions from TMSLISTs to a .mrk-file. 
        
        Optionally, the image information stored in T1 scan of the SimNIBS headmodel 
        will be added to the .mrk-file. By that, the ANT neuronavigation can ensure
        that the coordinates spaces match when importing the .mrk-file.
    
        Parameters
        ----------
        fn : str
            Filename of .mrk-file.
        tms_lists : list of TMSLIST
            List of TMSLISTs to be exported.
        subpath : str, optional
            Path to m2m-folder of the subject. If given, some header fields of 
            m2m-{subID}/T1.nii.gz (qform quaternion, qform offset, pixdim) will 
            be added as image information to the .mrk-file. The default is None.
        imageinfo : dict, optional
            Image information to add to to the .mrk-file.
            Can be used as alternative to providing subpath. subpath takes precedence
            in case both are set. The default is None.
            
        Written by Axel Thielscher, 2023
        """
        marker_dict = self._mrk_template()
        
        imInfo_internal = None
        if imageinfo is not None:
            imInfo_internal = imageinfo
        if subpath is not None:
            ff = SubjectFiles(subpath=subpath)
            imInfo_internal = self._imageinfo_from_nifti(ff.reference_volume)
                
        if imInfo_internal is not None:
            marker_dict['Configuration']['AssociatedImageInfo'] = imInfo_internal
            
        if isinstance(tms_lists,TMSLIST):
            tms_lists = [tms_lists]
            
        poslist = marker_dict['TargetMarkers']
        for tms_list in tms_lists:
            for p in tms_list.pos:
                pos_dict = self._target_template()
                pos_dict['CoilTarget']['Name'] = p.name
                pos_dict['CoilTarget']['CoilInfo']['Id'] = tms_list.fnamecoil
                pos_dict['CoilTarget']['CoilPose']['RigidTransformMatrix4x4'] = p.matsimnibs.tolist()
                poslist.append(pos_dict)
        
        marker_str = json.dumps(marker_dict)
        with open(fn, "w") as marker_file:
             marker_file.write(marker_str)

