function hdr = nifti_hdr(V)
% hdr = nifti_hdr(V)
% make header so that matrix V can be stored as nifti (using nifti_save)
%
% A. Thielscher, 2021


hdr.sizeof_hdr = 348;
hdr.endian = 'l';
hdr.data_type = '';
hdr.db_name = '';
hdr.extents = 0;
hdr.session_error = 0;
hdr.regular = 'r';
hdr.dim_info = 0;

hdr.dim = ones(8,1);
hdr.dim(1) = 4;
hdr.dim(2) = size(V,1);
hdr.dim(3) = size(V,2);
hdr.dim(4) = size(V,3);
hdr.dim(5) = size(V,4);
         
hdr.intent_p1 = 0;
hdr.intent_p2 = 0;
hdr.intent_p3 = 0;
hdr.intent_code = 0;
          
switch(class(V))
    case 'uint8',  hdr.datatype = 2; hdr.bitpix = 8; % 'uchar'
    case 'int16',  hdr.datatype = 4; hdr.bitpix = 16; % 'short'
    case 'int32',  hdr.datatype = 8; hdr.bitpix = 32; % 'int'
    case 'single', hdr.datatype = 16; hdr.bitpix = 32; % 'float'
    case 'double', hdr.datatype = 64; hdr.bitpix = 64;
    case 'uint16', hdr.datatype = 512; hdr.bitpix = 16; % 'ushort'
    case 'uint32', hdr.datatype = 768; hdr.bitpix = 32; % 'uint'
    otherwise
        error(['datatype ' class(V) ' not supported'])
end
         
hdr.slice_start = 0;
hdr.pixdim = ones(1,8);
hdr.scl_slope = 0;
hdr.scl_inter = 0;
hdr.slice_end = 0;
hdr.slice_code = 0;
hdr.xyzt_units = 0; %18; % assuming [mm] and [s]
hdr.cal_max = 0;
hdr.cal_min = 0;
hdr.slice_duration = 0;
hdr.toffset = 0;
hdr.glmax = round(double(max(V(:))));
hdr.glmin = round(double(min(V(:))));
hdr.descrip = ' ';
hdr.aux_file = ' ';
hdr.qform_code = 1;
hdr.sform_code = 1;
hdr.quatern_b = 0;
hdr.quatern_c = 0;
hdr.quatern_d = 0;
hdr.quatern_x = 0;
hdr.quatern_y = 0;
hdr.quatern_z = 0;
hdr.srow_x = [1 0 0 0];
hdr.srow_y = [0 1 0 0];
hdr.srow_z = [0 0 1 0];
hdr.intent_name = '';
hdr.magic = uint8([110    43    49     0]);
hdr.vol = V;
