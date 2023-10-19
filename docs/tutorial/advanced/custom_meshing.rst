.. _meshing_tutorial:

How to create and use a custom mesh
===================================

This example demonstrates how to create a mesh from a custom label image and how to set up simulations with the mesh.


Creating an example label image
-------------------------------
To get started, let's create a nifti file that contains a two-layer sphere with tissue label 5 (corresponding to scalp) for the outer shell and tissue label 2 (gray matter) for the inner part. In addition, let's add a smaller sphere with tissue type 17 somewhere in the center. Tissue label 17 is not a standard SimNIBS label. It is used here to define a new tissue type (could be a tumor or stroke lesion, for example).

* *Python*

  .. code-block:: python
  
	import numpy as np
	import nibabel as nib

	label_img = np.zeros((101,101,101), np.uint16)
	xv, yv, zv = np.meshgrid(np.linspace(-50,50,101),
				np.linspace(-50,50,101),
				np.linspace(-50,50,101))

	# make a two-layer sphere
	r = np.sqrt(xv**2 + yv**2 + zv**2)
	label_img[r<=40] = 5 # 5 corresponds to scalp
	label_img[r<=35] = 2 # 2 corresponds to gray matter

	# add a smaller decentered sphere
	r = np.sqrt((xv-15)**2 + yv**2 + zv**2)
	label_img[r<=15] = 17 # 17 is an arbitrary custom tissue label

	# save
	affine = np.eye(4)
	img = nib.Nifti1Image(label_img, affine)
	nib.save(img,'myspheres.nii.gz')

* *MATLAB*

  .. code-block:: matlab
  
	label_img = zeros([101,101,101],'uint16');
	[xv, yv, zv] = meshgrid(-50:50,-50:50,-50:50);

	% make a two-layer sphere
	r = sqrt(xv.^2 + yv.^2 + zv.^2);
	label_img(r<=40) = 5; % 5 corresponds to scalp
	label_img(r<=35) = 2; % 2 corresponds to gray matter

	% add a smaller decentered sphere
	r = sqrt((xv-15).^2 + yv.^2 + zv.^2);
	label_img(r<=15) = 17; % 17 is an arbitrary custom tissue label
	
	% save
	niftiwrite(label_img,'myspheres','Compressed',true)


Meshing the example label image
-------------------------------
To create a tetrahedral mesh from "myspheres.nii.gz", run in the command line

.. code-block::

  meshmesh myspheres.nii.gz myspheres.msh --voxsize_meshing 0.5

\

.. note:: The parameter --voxsize_meshing controls the internal voxel size to which the label image is upsampled before meshing. To better resolve thin structures in the mesh, a rule of thumb is to supply the label image at a resolution of 0.5 mm (preferred), or use the internal upsampling in case the image resolution is lower.


Run simulations with the custom mesh
------------------------------------
Running the simulations is very similar to the standard SimNIBS case. As difference, a m2m_{subID} folder is missing that contains information about the EEG positions, transformations to MNI and fsaverage space, etc. Therefore, the corresponding postprocessing options are not available. Electrode and TMS coil positions have to be supplied as coordinates, as EEG positions are not available. The coordinates can be determined by loading myspheres.nii.gz into a nifti viewer such as freeview.

A further difference is that we decided to include a custom tissue type with label 17 into the mesh, for which we have to define the conductivity.

* *Python*

.. literalinclude:: ../../../simnibs/examples/simulations/simulation_custom_mesh.py
   :language: python


* *MATLAB*

.. literalinclude:: ../../../simnibs/examples/simulations/simulation_custom_mesh.m
   :language: matlab


Output 
------

Windows, such as the following, should appear and show the results:

.. image:: ../../images/custommesh1.png
   :align: center

By making the tetrahedra visible and clipping the volume, the field inside the volume is shown (including the custom tissue type). As the conductivity of the custom tissue type was selected higher than the surrounding, the electric field strength is weaker there:

.. image:: ../../images/custommesh2.png
   :align: center

The results can also be found in the output folder 'simu_custom_mesh'.