function cat_debug
%cat_debug	print debug information for SPM12 and CAT12
%
% FORMAT cat_debug
%
%__________________________________________________________________________
% Christian Gaser
% $Id: cat_debug.m 766 2015-11-17 15:05:32Z gaser $

rev = '$Rev: 766 $';

% print last error
fprintf('\nLast error message:\n');
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------------------------------\n');
try
	er = lasterror;
	fprintf('%s\n',er.message);
	if isfield(er,'stack')
		for i=1:length(er.stack)
			fprintf('%s at line %g\n',char(er.stack(i).file),er.stack(i).line);
		end
	end
catch
	fprintf('%s\n',lasterr);
end

fprintf('-------------------------------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------------------------------\n');

fprintf('\nVersion information:\n');
fprintf('-------------------------------------------------------------------------------------\n');

ver

