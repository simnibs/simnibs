function S_out = expand_to_center_surround(tdcslist,subpath,varargin)
% 
% tdcslist = expand_to_center_surround(tdcslist,subpath [,'OptionName',OptionValue,...])
%
% Expands a TDCSLIST with a single electrode to a center-surround montage
%
%   tdcslist: the tdcslist to expand
%   subpath: path of the m2m-folder of the subject
%
%   Options are set using the option name followed by the value.
%       radius_surround: distance (centre-to-centre) between the centre and 
%                        surround electrodes (optional; standard: 50 mm)
%                        either a single number or an array with N entries
%                        (N: number of electrodes)             
%       pos_dir_1stsurround: position indicating the direction in which 
%                            the first surround electrode should be placed 
%                            (optional; standard: [])                        
%       N: number of surround electrodes (optional; standard: 4)
%       multichannel: when set to true: Simulation of multichannel stimulator 
%                     with each suround channel receiving 1/N-th of the
%                     center channel (optional; standard: false, i.e. all 
%                     surround electrodes connected to the same channel)
%       phis_surround: Angles in degree at which the electrodes will be placed 
%                      relative to the direction defined by pos_dir_1stsurround.
%                      (optional; standard: [], resulting in equal distances
%                       between surround electrodes)
%    Examples:
%       tdcslist = expand_to_center_surround(tdcslist,subpath,'multichannel',true)
                     
if nargin<2
    error('at least 2 inputs needed: TDCSLIST and subpath')
end

% standard settings will be filled in by simnibs
s.radius_surround=[];
s.pos_dir_1stsurround=[];
s.N=[];
s.multichannel=[];
s.phis_surround=[];
s=parse_input(s,varargin{:});

fn_in  = [tempname,'.mat'];
fn_out = [tempname,'.mat'];
save(fn_in,'-struct','tdcslist', '-v7')

% Run expand_to_center_surround
cmdstr = [simnibs_cli_call('expand_to_center_surround') ...
           ' -S ' fn_in ' -p ' subpath ' -F ' fn_out ];
fn = fieldnames(s);
for k=1:numel(fn)
    if ~isempty(s.(fn{k}))
        if islogical(s.(fn{k}))
            if s.(fn{k})
                cmdstr = [cmdstr ' --' fn{k}];
            end
        else
            hlpstr = mat2str(s.(fn{k}));
            if length(s.(fn{k})) > 1
                hlpstr=hlpstr(2:end-1);
            end
            cmdstr = [cmdstr ' --' fn{k} ' ' hlpstr];
        end
    end
end
[status,result] = system(cmdstr);

if status ~= 0
    delete(fn_in);
    delete(fn_out);
    error('There was an error running expand_to_center_surround:\n %s',result)
end
S_out = load(fn_out);
delete(fn_in);
delete(fn_out);

