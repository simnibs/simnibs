function cat_stat_diff(P,rel,glob)
% compute difference between two images or surfaces
% image2 - image1
% output name will be diff_{name of image2}
% or diff_rel{name of image2} for relative differences
%
% FORMAT cat_stat_diff(P,rel,glob)
% P     - filenames for image 1/2
% rel   - compute relative difference rather than absolute (1 for relative)
% glob  - normalize global mean of images
%__________________________________________________________________________
% Christian Gaser
% $Id: cat_stat_diff.m 1129 2017-05-09 12:58:31Z gaser $

global SWD

% get filenames
%----------------------------------------------------------------------------
if nargin < 3, glob = 0; end
if nargin < 2, rel = 0; end

if nargin == 0
    % select images for each subject
    don = 0;
    for i = 1:1000,
        P = spm_select([0 Inf],'.*(img|nii|gii)',['Image(s), subj ' num2str(i)]);
        if size(P,1) < 2, don = 1; break; end;
        try
          V{i} = spm_data_hdr_read(P);
        catch
          error('Error reading file. Ensure that you either have an image file or a surface texture file with values.');
        end
    end
    rel  = 0;
    glob = 0;
else
  try
    V{i} = spm_data_hdr_read(P);
  catch
    error('Error reading file. Ensure that you either have an image file or a surface texture file with values.');
  end
end

if isfield(V{1}(1).private,'cdata')
     surf = 1;
else surf = 0; end
  
for i = 1:length(V);
    Vi = V{i};
    
    n = length(Vi);
    % compute global means to normalize images to same value
    if glob & ~surf
       gm=zeros(n,1);
       disp('Calculating globals...');
       for j=1:n, gm(j) = spm_global(Vi(j)); end
       gm_all = mean(gm);
       for j=1:n
            Vi(j).pinfo(1,:) = gm_all*Vi(j).pinfo(1,:)/gm(j);
       end
    end
    
    for j=2:n
    
        Q = char(Vi(j).fname);
        [pth,nm,xt,vr] = spm_fileparts(Q);

        if rel
            outname = ['diffrel_' nm xt vr];
            str='relative difference %';
        else
            outname = ['diff_' nm xt vr];
            str='difference';
        end

        Q = fullfile(pth,outname);
    
        Vo = struct(    'fname',    Q,...
            'dim',      Vi(j).dim,...
            'mat',      Vi(j).mat,...
            'descrip',  str);

        Vo.dt = [spm_type('float32') spm_platform('bigend')];
        
        [pth1,nm1,xt1,vr1] = spm_fileparts(Vi(1).fname);
        [pth2,nm2,xt2,vr2] = spm_fileparts(Vi(j).fname);
        if surf
          fprintf('Calculate s2-s1: %s - %s\n',[nm2 xt2],[nm1 xt1]);
        else
          fprintf('Calculate i2-i1: %s - %s\n',[nm2 xt2],[nm1 xt1]);
        end

        % implement the computing
        %---------------------------------------------------------------------------
        if surf
          if rel, formula='200*(s2-s1)./(s1+s2+eps)';
          else, formula='s2-s1'; end
          iname{1} = Vi(1).fname;
          iname{2} = Vi(j).fname;
          spm_mesh_calc(iname,Vo.fname,formula);
        else
          if rel, formula='200*(i2-i1)./(i1+i2+eps)';
          else, formula='i2-i1'; end
          spm_imcalc(Vi([1 j]),Vo,formula,{0 0 4 1});
        end
    end

end
