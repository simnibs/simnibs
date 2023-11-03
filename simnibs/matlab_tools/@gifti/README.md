# GIfTI library for MATLAB/Octave

This GIfTI library allows to handle the GIfTI Geometry file format from the Neuroimaging
Informatics Technology Initiative (NIfTI) using a MATLAB/Octave class:
  * GIfTI: https://www.nitrc.org/projects/gifti/
  * NIfTI: https://nifti.nimh.nih.gov/

It relies on external libraries:
  * [yxml](https://dev.yorhel.nl/yxml), by Yoran Heling
  * [Base64](https://stackoverflow.com/a/37109258), by polfosol
  * [miniz](https://github.com/richgel999/miniz), by Rich Geldreich
  * [mVTK](https://www.artefact.tk/software/matlab/mvtk/), by Guillaume Flandin
  * [JSONio](https://www.artefact.tk/software/matlab/jsonio/), by Guillaume Flandin

Note that these tools are already included in the GIfTI library provided
here, so you don't need to install them separately.

There are import facilities from [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki/FileFormats),
[VTK](https://vtk.org/), [Wavefront OBJ](https://www.wikipedia.org/wiki/Wavefront_.obj_file),
[Stanford PLY](https://www.wikipedia.org/wiki/PLY_%28file_format%29),
[STL](https://www.wikipedia.org/wiki/STL_%28file_format%29) and
[MZ3](https://github.com/neurolabusc/surf-ice/tree/master/mz3) file formats.

There are export facilities to [VTK](https://vtk.org/),
[Collada](https://www.khronos.org/collada/),
[IDTF](http://www.meshlab.net/),
[Wavefront OBJ](https://www.wikipedia.org/wiki/Wavefront_.obj_file) and
[JS/JSON](https://plot.ly/javascript/) (for [Plotly](https://plot.ly/javascript/)) file formats.

This library is also part of [SPM](https://www.fil.ion.ucl.ac.uk/spm/).

## INSTALLATION

[MATLAB](https://www.mathworks.com/products/matlab.html) R2007a or above is required to use most of the features of
this toolbox. [GNU Octave](https://www.octave.org/) is also supported.
 
All the code is embedded in a `@gifti` class. To install it, all you need is to 
make sure that the directory containing `@gifti` is in MATLAB path:
 
```matlab
addpath /home/login/Documents/MATLAB/gifti
```
 
The library relies on a number of C-MEX files (`zstream`, `base64`, `xml_parser`).
Compiled versions for 64 bit MATLAB on Windows, Linux and Mac are provided
but they can easily be compiled by yourself otherwise, see `@gifti/private/Makefile`.
  
## TUTORIAL
 
In the following, we use the files contained in `BV_GIFTI.tar.gz` (BrainVISA examples),
available from the [NITRC website](https://www.nitrc.org/frs/?group_id=75&release_id=123): 
   
```matlab
% Read the GIfTI surface file
g = gifti('sujet01_Lwhite.surf.gii')

% Read the GIfTI shape file
c = gifti('sujet01_Lwhite.shape.gii')

% Display mesh
figure; plot(g);
% Display mesh with curvature
figure; plot(g,c);
```
   
In a similar way, a gifti object can be created from scratch and saved to a file:
   
```matlab
load mri
D = squeeze(D);
Ds = smooth3(D);
g = gifti(isosurface(Ds,5))

h = plot(g);
daspect([1,1,.4]); view(45,30); axis tight
lightangle(45,30);
set(h,'SpecularColorReflectance',0,'SpecularExponent',50)

save(g,'mri.surf.gii','Base64Binary');
```

See also: [export to Plotly](https://gllmflndn.github.io/gifti/).
