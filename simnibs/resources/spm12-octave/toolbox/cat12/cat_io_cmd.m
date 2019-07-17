function varargout = cat_io_cmd(str,style,strlength,verb,stime)
% ______________________________________________________________________
% Writes a string str with a specific style or color followed by a set 
% of blanks to fit a length of strlength character. 
% The style variable can be a color that should be a 1x3 matrix with RGB
% value from 0 to 1. Otherwise, it can be a str 'error' (red), 'warning'
% (orange), or 'comment' (blue). If no style is given it is black. 
%
%   stime = cat_io_cmd(str[,style,strlength,verb,stime])
%
% Example: 
%   stime = cat_io_cmd('Testfunction','comment',63); pause(3);
%   fprintf('%3.0fs\n',etime(clock,stime));
%
%   Testfunction                                                      3s
%
%
%   stime1 = cat_io_cmd('Function 0815 (without cleanup)'); fprintf('\n');
%   stime  = cat_io_cmd('  Substep 1','g5','',1);       pause(0.74);
%   stime  = cat_io_cmd('  Substep 2','g5','',1,stime); pause(1.12);
%   stime  = cat_io_cmd('  Substep 3','g5','',1,stime); pause(1.65);
%   stime  = cat_io_cmd(' ','','',1,stime); % last Substeptime
%   stime1 = cat_io_cmd('Function 0816 (with cleanup)','','',1,stime1); fprintf('\n');
%   stime  = cat_io_cmd('  Substep A','g5','',1);       pause(2.34);
%   stime  = cat_io_cmd('  Substep B','g5','',1,stime); pause(1.12);
%   stime  = cat_io_cmd('cleanup',2,'',1,stime1); 
%   stime1 = cat_io_cmd('Function 0817','','',1); pause(2.12);
%   stime  = cat_io_cmd('','','',1,stime1);
%   fprintf('done.\n');
%
% see also: cat_io_cprintf for colored command line output.
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_io_cmd.m 1212 2017-11-10 12:39:37Z dahnke $ %

  %#ok<*NASGU> 

  if ~exist('verb','var') || isempty(verb), verb=1; end
  if ~exist('strlength','var') || isempty(strlength), strlength=65; end
  strlength2 = strlength;
  if cat_io_matlabversion<20110, strlength2 = strlength2+1; end
  
  if verb
    switch str
      case 'testthisfunction'
        stime1 = cat_io_cmd('Function 0815 (without cleanup)'); fprintf('\n');
        stime  = cat_io_cmd('  Substep 1','g5','',1);       pause(0.74);
        stime  = cat_io_cmd('  Substep 2','g5','',1,stime); pause(1.12);
        stime  = cat_io_cmd('  Substep 3','g5','',1,stime); pause(1.65);
        stime  = cat_io_cmd(' ','','',1,stime); % last Substeptime
        stime1 = cat_io_cmd('Function 0816 (with cleanup)','','',1,stime1); fprintf('\n');
        stime  = cat_io_cmd('  Substep A','g5','',1);       pause(2.34);
        stime  = cat_io_cmd('  Substep B','g5','',1,stime); pause(1.12);
        stime  = cat_io_cmd('cleanup',2,'',1,stime1); % cleanup with number of lines 
        stime1 = cat_io_cmd('Function 0817','','',1); pause(2.12);
        stime  = cat_io_cmd('','','',1,stime1);
        fprintf('done.\n');
        
% ---        
% this case does not work equal on all matlab versions, because in some 
% versions an addition space is placed by cprintf by unknown reasons
%       case 'cleanup'
%         % this works not for all matlab versions correctly
%         fprintf(sprintf('%s',repmat('\b',1,style * (strlength2+7) - 5)));
%         if exist('stime','var') && ~isempty(stime)
%           cat_io_cmd('','','',verb,stime);
%         end
% ---    
      case 'cleanup'
        if exist('stime','var') && ~isempty(stime)
          fprintf('% 5.0fs\n',etime(clock,stime));
        end        
% ---        
      otherwise
        if exist('stime','var') && ~isempty(stime)
          fprintf('% 5.0fs\n',etime(clock,stime));
        end

        if ~isempty(str) 
          if exist('style','var')
            cat_io_cprintf(style,sprintf('%s%s',str,repmat(' ',1,1+strlength-length(str)))); 
          else
            fprintf('%s:%s',str,repmat(' ',1,strlength-length(str))); 
          end
        end
    end
  end
  
  if ~strcmp(str,'testthisfunction') && nargout>0;
    varargout{1} = clock; 
  end  
end