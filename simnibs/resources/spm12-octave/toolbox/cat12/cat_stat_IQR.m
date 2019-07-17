function cat_stat_IQR(p)
%cat_stat_IQR to read weighted overall image quality (IQR) from xml-files
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_stat_IQR.m 1026 2016-10-03 20:44:40Z gaser $

fid = fopen(p.iqr_name,'w');

if fid < 0
	error('No write access: check file permissions or disk space.');
end

spm_progress_bar('Init',length(p.data_xml),'Load xml-files','subjects completed')
for i=1:length(p.data_xml)
    xml = cat_io_xml(deblank(p.data_xml{i})); 
    try
      iqr = xml.qualityratings.IQR;
    catch % also try to use old versions
      iqr = xml.QAM.QM.rms;
    end

    [pth,nam]     = spm_fileparts(p.data_xml{i});
        fprintf(fid,'%s\n',iqr);
        fprintf('%s\n',iqr);
    spm_progress_bar('Set',i);  
end
spm_progress_bar('Clear');


if fclose(fid)==0
	fprintf('\nValues saved in %s.\n',p.iqr_name);
end
