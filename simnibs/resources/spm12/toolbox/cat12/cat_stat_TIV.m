function varargout = cat_stat_TIV(p)
%cat_stat_TIV to read total intracranial volume (TIV) from xml-files
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_stat_TIV.m 1234 2017-12-04 10:52:18Z gaser $

if ~p.calcvol_TIV
  fprintf('%60s\t%7s\t%7s\t%7s\t%7s\t%7s\n','Name','Total','GM','WM','CSF','WMH');
end
fid = fopen(p.calcvol_name,'w');

if fid < 0
	error('No write access: check file permissions or disk space.');
end

if p.calcvol_TIV
  calcvol = zeros(length(p.data_xml),1);
else
  calcvol = zeros(length(p.data_xml),5);
end

spm_progress_bar('Init',length(p.data_xml),'Load xml-files','subjects completed')
for i=1:length(p.data_xml)
    xml = cat_io_xml(deblank(p.data_xml{i})); 
    tmp  = xml.subjectmeasures.vol_abs_CGW; 
    
    name = spm_str_manip(xml.filedata.fname,'a50');

    % only save TIV
    if p.calcvol_TIV
        fprintf(fid,'%7.2f\n',sum(tmp));
        calcvol(i) = sum(tmp);
        fprintf('%60s\t%7.2f\n',spm_str_manip(name,'l60'),sum(tmp));
    else % also save GM/WM/CSF 
        fprintf(fid,'%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n',sum(tmp),tmp(2),tmp(3),tmp(1),tmp(4));
        calcvol(i,:) = [sum(tmp),tmp(2),tmp(3),tmp(1),tmp(4)];
        fprintf('%60s\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n',spm_str_manip(name,'l60'),sum(tmp),tmp(2),tmp(3),tmp(1),tmp(4));
    end
    spm_progress_bar('Set',i);  
end
spm_progress_bar('Clear');

if fclose(fid)==0
	fprintf('\nValues saved in %s.\n',p.calcvol_name);
    if nargout == 1
	    varargout{1}.calcvol = calcvol;
    end
end
