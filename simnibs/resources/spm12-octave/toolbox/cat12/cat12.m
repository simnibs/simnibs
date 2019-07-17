function varargout = cat12(varargin)
% ______________________________________________________________________
% CAT12 Toolbox wrapper to start CAT with different user modes or 
% default files.  Changing the user mode requires restarting of CAT and
% SPM.  The expert user mode allows to control further parameters and  
% semi-evaluated functions, whereas the developer mode contain parameter
% for internal tests and unsafe functions.
% 
%   cat12(action)
%   
%   CAT user modes:
%     action = ['default','expert','developer'] 
%
%   CAT default files for other species (in development):
%     action = ['oldwoldmonkeys'|'greaterapes']
%
%   CAT start with own default files:
%     action = 'select' 
%     action = 'mypath/cat_defaults_mydefaults'
%
% ______________________________________________________________________
% Christian Gaser, Robert Dahnke
% $Id: cat12.m 1275 2018-02-12 22:00:05Z gaser $

% CAT12 M-file for cat12.fig
%      CAT12, by itself, creates a new CAT12 or raises the existing
%      singleton*.
%
%      H = CAT12 returns the handle to a new CAT12 or the handle to
%      the existing singleton*.
%
%      CAT12('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAT12.M with the given input arguments.
%
%      CAT12('Property','Value',...) creates a new CAT12 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DEM_demo_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cat12_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above title to modify the response to help cat12

% Last Modified by GUIDE v2.5 23-Jan-2018 22:36:22

if nargin==0 
  spm_cat12;
  return;
elseif nargin==1 && ~strcmp(varargin{1},'fig')
  spm_cat12(varargin{1});
  return;
elseif nargin==2 && ~strcmp(varargin{1},'fig')
  spm_cat12(varargin{1},varargin{2});
  return;
end

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cat12_OpeningFcn, ...
                   'gui_OutputFcn',  @cat12_OutputFcn, ...
                   'gui_LayoutFcn',  @cat12_LayoutFcn, ...
                   'gui_Callback',   []);
                   
if nargin && ~strcmp(varargin{1},'fig') && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before cat12 is made visible.
function cat12_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cat12 (see VARARGIN)

% Choose default command line output for cat12
handles.output = hObject; 

% Update handles structure
guidata(hObject, handles);

% enable/disable different menus if TFCE is installed or not
if exist(fullfile(spm('dir'),'toolbox','TFCE'))
    set(handles.popupmenu052,'String',{'Treshold-Free Cluster Enhancement...','Call TFCE Toolbox'});
else
    set(handles.popupmenu052,'String',{'Treshold-Free Cluster Enhancement...','Install TFCE Toolbox'});
end

% --- Outputs from this function are returned to the command line.
function varargout = cat12_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spm_clf('Interactive'); 

% Get default command line output from handles structure
varargout{1} = handles.output;

expert  = cat_get_defaults('extopts.expertgui'); 
species = cat_get_defaults('extopts.species'); 
switch expert
  case 1, set(handles.CAT,'color', [0.85 0.85 0.85]);
  case 2, set(handles.CAT,'color', [0.93 0.93 0.93]); 
end

FS = spm('FontSizes');

% This creates the 'background' image
handles.ha = axes('units','normalized','position',[0 0.87 1 0.13]);
I = imread(fullfile(spm('dir'),'toolbox','cat12','html','images','contact.jpg'));
imagesc(I);
axis off; 
text(80,140,'Computational Anatomy Toolbox','Color',[1 1 1],'Fontsize',FS(14),'Fontweight','bold');
switch species
  case 'human',           speciesdisp = ''; 
  case 'ape_greater',     speciesdisp = ' (greater apes)';
  case 'ape_lesser',      speciesdisp = ' (lesser apes)';
  case 'monkey_oldworld', speciesdisp = ' (oldworld monkeys)'; 
  case 'monkey_newworld', speciesdisp = ' (newworld monkeys)'; 
  case 'chimpanzee',      speciesdisp = ' (chimpanzee)'; 
  case 'dog',             speciesdisp = ' (dogs)'; 
  otherwise               speciesdisp = ''; 
end
switch expert
  case 1, text(80,90,['Expert Mode'    speciesdisp],'Color',[0.1 0.7 1.0],'Fontsize',FS(10),'Fontweight','bold'); 
  case 2, text(80,90,['Developer Mode' speciesdisp],'Color',[1.0 0.0 0.0],'Fontsize',FS(10),'Fontweight','bold');
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function CAT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


%-------------------------------------------------------------------

% --- Executes on button press in pushbutton011.
function pushbutton011_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.cat.estwrite');

% --- Executes on button press in pushbutton012.
function pushbutton012_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.cat.tools.long');

% --- Executes on button press in pushbutton021.
function pushbutton021_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.cat.tools.showslice');

% --- Executes on button press in pushbutton051.
function pushbutton051_Callback(hObject, eventdata, handles)
spm_jobman('interactive','cat_stat_factorial_design.m');

% --- Executes on button press in pushbutton053.
function pushbutton053_Callback(hObject, eventdata, handles)
cat_stat_spm;

% --- Executes on button press in pushbutton032.
function pushbutton032_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.cat.stools.surfresamp');

% --- Executes on button press in pushbutton034.
function pushbutton034_Callback(hObject, eventdata, handles)
cat_surf_display;

% --- Executes on button press in pushbutton074.
function pushbutton074_Callback(hObject, eventdata, handles)
F = spm_figure('FindWin','Menu');
% close SPM windows, if no Menu window exist
if isempty(F)
  spm('Quit')
end
close(gcf);


% --- Executes on selection change in popupmenu033.
function popupmenu033_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu033 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu033 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu033


% --- Executes during object creation, after setting all properties.
function popupmenu033_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu033 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton054.
function pushbutton054_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.cat.tools.calcvol');


% --- Executes on button press in pushbutton073.
function pushbutton073_Callback(hObject, eventdata, handles)
cat_io_senderrormail;


% --- Executes on selection change in popupmenu052.
function popupmenu052_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu052 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu052 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu052
% --- Executes during object creation, after setting all properties.

% Determine the selected data set.
if get(hObject,'Value') == 2
    if exist(fullfile(spm('dir'),'toolbox','TFCE'))
        % call TFCE toolbox 
        spm_TFCE;
    else % install TFCE toolbox
        d0 = spm('Dir');
        d = fullfile(spm('Dir'),'toolbox'); 
        s = unzip('http://www.neuro.uni-jena.de/tfce/tfce_latest.zip', d);
        addpath(d0);
        rehash
        rehash toolboxcache;
        toolbox_path_cache
        eval(['spm fmri;clear cat_version;spm_cat12']);
    end
end


% --- Executes when CAT is resized.
function CAT_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to CAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu041.
function popupmenu041_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu041 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu041 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu041


% --- Executes during object creation, after setting all properties.
function popupmenu041_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu041 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu062.
function popupmenu062_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu062 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu062 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu062


% --- Executes during object creation, after setting all properties.
function popupmenu062_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu062 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton042.
function pushbutton042_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton042 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cat_stat_analyze_ROIs;


% --- Executes on selection change in popupmenu071.
function popupmenu071_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu071 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu071 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu071


% --- Executes during object creation, after setting all properties.
function popupmenu071_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu071 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu072.
function popupmenu072_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu072 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu072 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu072


% --- Executes during object creation, after setting all properties.
function popupmenu072_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu072 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu061.
function popupmenu061_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu061 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu061 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu061


% --- Executes during object creation, after setting all properties.
function popupmenu061_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu061 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu031.
function popupmenu031_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu031 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu031 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu031


% --- Executes during object creation, after setting all properties.
function popupmenu031_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu031 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popupmenu052_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu052 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu022.
function popupmenu022_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu022 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu022 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu022


% --- Executes during object creation, after setting all properties.
function popupmenu022_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu022 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Creates and returns a handle to the GUI figure. 
function h000 = cat12_LayoutFcn(policy)
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

persistent hsingleton;
if strcmpi(policy, 'reuse') & ishandle(hsingleton)
    h000 = hsingleton;
    return;
end


appdata = [];
appdata.GUIDEOptions = struct(...
    'active_h', [], ...
    'taginfo', struct(...
    'figure', 2, ...
    'uipanel', 38, ...
    'text', 13, ...
    'axes', 7, ...
    'pushbutton', 175, ...
    'popupmenu', 18, ...
    'togglebutton', 2, ...
    'activex', 2), ...
    'override', 1, ...
    'release', 13, ...
    'resize', 'simple', ...
    'accessibility', 'callback', ...
    'mfile', 1, ...
    'callbacks', 1, ...
    'singleton', 1, ...
    'syscolorfig', 0, ...
    'blocking', 0);
appdata.lastValidTag = 'CAT';
appdata.GUIDELayoutEditor = [];
appdata.initTags = struct(...
    'handle', [], ...
    'tag', 'CAT');

FS = spm('FontSizes');

%-------------------------------------------------------------------
% gui positions
x2  = [0.07 0.53];   % x-pos two columns
x2w = 0.39;          % width button/popupmenu for two columns
y1t = 0.77;          % y-pos one row text
y1b = 0.39;          % y-pos one row button
y1p = 0.34;          % y-pos one row popupmenu
y2b = [0.592 0.215]; % y-pos two rows button
y2p = [0.562 0.185]; % y-pos two rows popupmenu
y1h  = 0.45;         % height button/popupmenu for one row
yph  = 0.28;         % height button/popupmenu for two rows
%-------------------------------------------------------------------

h000 = figure(...
'Units','characters',...
'PaperUnits','centimeters',...
'Color',get(0,'defaultfigureColor'),...
'Colormap',gray(64),...
'DockControls','off',...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','cat12',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'PaperSize',[21 29.7],...
'PaperType','A4',...
'Position',[36 44 80 54],...
'ResizeFcn',@(hObject,eventdata)cat12('CAT_ResizeFcn',hObject,eventdata,guidata(hObject)),...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)cat12('CAT_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'UserData',[],...
'Tag','CAT',...
'Visible','on');

%-------------------------------------------------------------------
appdata = [];
appdata.lastValidTag = 'uipanel010';

h010 = uipanel(...
'Parent',h000,...
'FontSize',FS(7),...
'Title','Preprocessing',...
'Tag','uipanel010',...
'UserData',[],...
'Clipping','on',...
'BackgroundColor',[0.8 0.8 0.8],...
'Position',[0.069 0.76658163265306 0.853 0.095],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton011';

h011 = uicontrol(...
'Parent',h010,...
'Units','normalized',...
'Callback',@(hObject,eventdata)cat12('pushbutton011_Callback',hObject,eventdata,guidata(hObject)),...
'FontSize',FS(6),...
'Position',[x2(1) y1b x2w y1h],...
'String','Segment Data',...
'TooltipString','Segment (cross-sectional) data',...
'Tag','pushbutton011',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton012';

h012 = uicontrol(...
'Parent',h010,...
'Units','normalized',...
'Callback',@(hObject,eventdata)cat12('pushbutton012_Callback',hObject,eventdata,guidata(hObject)),...
'FontSize',FS(6),...
'Position',[x2(2) y1b x2w y1h],...
'String','Segment Longitudinal Data',...
'TooltipString','Segment longitudinal data that are characterized by small changes and short times between time points',...
'UserData',[],...
'Tag','pushbutton012',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

%-------------------------------------------------------------------
appdata = [];
appdata.lastValidTag = 'uipanel020';

h020 = uipanel(...
'Parent',h000,...
'FontSize',FS(7),...
'Title','Check Data Quality',...
'Tag','uipanel020',...
'Clipping','on',...
'BackgroundColor',[0.8 0.8 0.8],...
'Position',[0.069 0.660714285714285 0.853 0.095],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton021';

h021 = uicontrol(...
'Parent',h020,...
'Units','normalized',...
'Callback',@(hObject,eventdata)cat12('pushbutton021_Callback',hObject,eventdata,guidata(hObject)),...
'FontSize',FS(6),...
'Position',[x2(1) y1b x2w y1h],...
'String','Display One Slice For All Images',...
'TooltipString','Display a selected slice for all images',...
'UserData',[],...
'Tag','pushbutton021',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupmenu022';

h022 = uicontrol(...
'Parent',h020,...
'Units','normalized',...
'Callback','spm(''PopUpCB'',gcbo)',...
'FontSize',FS(6),...
'ListboxTop',0,...
'Position',[x2(2) y1p x2w y1h],...
'String',{  'Check Sample Homogeneity...'; 'VBM Data'; 'Surface Data'; 'SPM Design' },...
'Style','popupmenu',...
'TooltipString','Check sample homogeneity using correlation across sample',...
'Value',1,...
'UserData',{  'spm_jobman(''interactive'','''',''spm.tools.cat.tools.check_cov'');' 'spm_jobman(''interactive'','''',''spm.tools.cat.stools.check_mesh_cov'');' 'cat_stat_check_SPM;'},...
'Tag','popupmenu022',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

%-------------------------------------------------------------------
appdata = [];
appdata.lastValidTag = 'uipanel030';

h030 = uibuttongroup(...
'Parent',h000,...
'FontSize',FS(7),...
'Title','Surface Tools',...
'Tag','uipanel030',...
'Clipping','on',...
'BackgroundColor',[0.8 0.8 0.8],...
'Position',[0.069 0.505102040816326 0.853 0.14],...
'SelectedObject',[],...
'SelectionChangeFcn',[],...
'OldSelectedObject',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupmenu031';

h031 = uicontrol(...
'Parent',h030,...
'Units','normalized',...
'Callback','spm(''PopUpCB'',gcbo)',...
'FontSize',FS(6),...
'ListboxTop',0,...
'Position',[x2(1) y2p(1) x2w yph],...
'String',{  'Extract & Map Surface Data...'; 'Extract Additional Surface Parameters '; 'Map Volume (Native Space) to Individual Surface'; 'Nonlinear Co-Register and Map Volume (Native Space) to Individual Surface'; 'Map Normalized Volume (Template Space) to Template Surface' },...
'Style','popupmenu',...
'TooltipString','Choose a method to extract data from volumes or surfaces such as gyrification or sulcus depth',...
'Value',1,...
'UserData',{  'spm_jobman(''interactive'','''',''spm.tools.cat.stools.surfextract'')' 'spm_jobman(''interactive'','''',''spm.tools.cat.stools.vol2surf'')' 'spm_jobman(''interactive'',''cat_surf_coregvol2surf.m'')' 'spm_jobman(''interactive'','''',''spm.tools.cat.stools.vol2surftemp'')' },...
'Tag','popupmenu031',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton032';

h032 = uicontrol(...
'Parent',h030,...
'Units','normalized',...
'Callback',@(hObject,eventdata)cat12('pushbutton032_Callback',hObject,eventdata,guidata(hObject)),...
'FontSize',FS(6),...
'Position',[x2(2) y2b(1) x2w yph],...
'String','Resample & Smooth Surfaces',...
'TooltipString','Resample and smooth surface data',...
'Tag','pushbutton032',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupmenu033';

h033 = uicontrol(...
'Parent',h030,...
'Units','normalized',...
'Callback','spm(''PopUpCB'',gcbo)',...
'FontSize',FS(6),...
'ListboxTop',0,...
'Position',[x2(1) y2p(2) x2w yph],...
'String',{  'Surface Calculator...'; 'Surface Calculator '; 'Surface Calculator (Subject-wise)' },...
'Style','popupmenu',...
'TooltipString','Choose a method to extract data from volumes or surfaces such as gyrification or sulcus depth',...
'Value',1,...
'UserData',{  'spm_jobman(''interactive'','''',''spm.tools.cat.stools.surfcalc'')' 'spm_jobman(''interactive'','''',''spm.tools.cat.stools.surfcalcsub'')' },...
'Tag','popupmenu033',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton034';

h034 = uicontrol(...
'Parent',h030,...
'Units','normalized',...
'Callback',@(hObject,eventdata)cat12('pushbutton034_Callback',hObject,eventdata,guidata(hObject)),...
'FontSize',FS(6),...
'Position',[x2(2) y2b(2) x2w yph],...
'String','Display Surfaces',...
'TooltipString','Display surfaces',...
'Tag','pushbutton034',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

%-------------------------------------------------------------------
appdata = [];
appdata.lastValidTag = 'uipanel040';

h040 = uipanel(...
'Parent',h000,...
'FontSize',FS(7),...
'Title','Regions of Interest Tools',...
'Tag','uipanel040',...
'Clipping','on',...
'BackgroundColor',[0.8 0.8 0.8],...
'Position',[0.069 0.392857142857143 0.853 0.095],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupmenu041';

h041 = uicontrol(...
'Parent',h040,...
'Units','normalized',...
'Callback','spm(''PopUpCB'',gcbo)',...
'FontSize',FS(6),...
'ListboxTop',0,...
'Position',[x2(1) y1p x2w y1h],...
'String',{  'Extract ROI Data...'; 'Extract ROI-based Surface Values'; 'Estimate Mean Values inside ROI for external Analysis' },...
'Style','popupmenu',...
'TooltipString','Extract ROI data using Atlases',...
'Value',1,...
'UserData',{  'spm_jobman(''interactive'','''',''spm.tools.cat.stools.surf2roi'');' 'spm_jobman(''interactive'','''',''spm.tools.cat.tools.calcroi'');' },...
'Tag','popupmenu041',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton042';

h042 = uicontrol(...
'Parent',h040,...
'Units','normalized',...
'Callback',@(hObject,eventdata)cat12('pushbutton042_Callback',hObject,eventdata,guidata(hObject)),...
'FontSize',FS(6),...
'Position',[x2(2) y1b x2w y1h],...
'String','Analyze ROIs',...
'TooltipString','Analyze ROIs using existing SPM design',...
'Tag','pushbutton042',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

%-------------------------------------------------------------------
appdata = [];
appdata.lastValidTag = 'uipanel050';

h050 = uibuttongroup(...
'Parent',h000,...
'FontSize',FS(7),...
'Title','Statistical Analysis',...
'Tag','uipanel050',...
'Clipping','on',...
'BackgroundColor',[0.8 0.8 0.8],...
'Position',[0.069 0.239795918367347 0.853 0.14],...
'SelectedObject',[],...
'SelectionChangeFcn',[],...
'OldSelectedObject',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton051';

h051 = uicontrol(...
'Parent',h050,...
'Units','normalized',...
'Callback',@(hObject,eventdata)cat12('pushbutton051_Callback',hObject,eventdata,guidata(hObject)),...
'FontSize',FS(6),...
'Position',[x2(1) y2b(1) x2w yph],...
'String','Basic Models',...
'TooltipString','Basic staitistical models (2nd-level)',...
'Tag','pushbutton051',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupmenu052';

h052 = uicontrol(...
'Parent',h050,...
'Units','normalized',...
'Callback',@(hObject,eventdata)cat12('popupmenu052_Callback',hObject,eventdata,guidata(hObject)),...
'FontSize',FS(6),...
'ListboxTop',0,...
'Position',[x2(2) y2p(1) x2w yph],...
'String',{  'Treshold-Free Cluster Enhancement...'; 'Estimate'; 'Results' },...
'Style','popupmenu',...
'TooltipString','Call threshold-free cluster enhancement (TFCE) Toolbox',...
'Value',1,...
'UserData',[],...
'Tag','popupmenu052',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton053';

h053 = uicontrol(...
'Parent',h050,...
'Units','normalized',...
'Callback',@(hObject,eventdata)cat12('pushbutton053_Callback',hObject,eventdata,guidata(hObject)),...
'FontSize',FS(6),...
'Position',[x2(1) y2b(2) x2w yph],...
'String','Estimate Surface Models',...
'TooltipString','Model estimation for surface data',...
'Tag','pushbutton053',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton054';

h054 = uicontrol(...
'Parent',h050,...
'Units','normalized',...
'Callback',@(hObject,eventdata)cat12('pushbutton054_Callback',hObject,eventdata,guidata(hObject)),...
'FontSize',FS(6),...
'Position',[x2(2) y2b(2) x2w yph],...
'String','Estimate TIV',...
'TooltipString','Estimate Total Intracranial Volume (TIV)',...
'Tag','pushbutton054',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

%-------------------------------------------------------------------
appdata = [];
appdata.lastValidTag = 'uipanel060';

h060 = uibuttongroup(...
'Parent',h000,...
'FontSize',FS(7),...
'Title','Data Presentation',...
'Tag','uipanel060',...
'Clipping','on',...
'BackgroundColor',[0.8 0.8 0.8],...
'Position',[0.069 0.13265306122449 0.853 0.095],...
'SelectedObject',[],...
'SelectionChangeFcn',[],...
'OldSelectedObject',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupmenu061';

h061 = uicontrol(...
'Parent',h060,...
'Units','normalized',...
'Callback','spm(''PopUpCB'',gcbo)',...
'FontSize',FS(6),...
'ListboxTop',0,...
'Position',[x2(1) y1p x2w y1h],...
'String',{  'Transform SPM-maps...'; 'spmT images'; 'spmF images'; 'spmT surfaces'; 'spmF surfaces' },...
'Style','popupmenu',...
'TooltipString','Threshold and transform SPM-maps to (log-scaled) p-maps or correlation maps',...
'Value',1,...
'UserData',{  'spm_jobman(''interactive'','''',''spm.tools.cat.tools.T2x'');' 'spm_jobman(''interactive'','''',''spm.tools.cat.tools.F2x'');' 'spm_jobman(''interactive'','''',''spm.tools.cat.tools.T2x_surf'');' 'spm_jobman(''interactive'','''',''spm.tools.cat.tools.F2x_surf'');' },...
'Tag','popupmenu061',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupmenu062';

h062 = uicontrol(...
'Parent',h060,...
'Units','normalized',...
'Callback','spm(''PopUpCB'',gcbo)',...
'FontSize',FS(6),...
'ListboxTop',0,...
'Position',[x2(2) y1p x2w y1h],...
'String',{  'Display Results...'; 'Slice Overlay'; 'Display Surface Results' },...
'Style','popupmenu',...
'TooltipString','Display results as overlay or render view',...
'Value',1,...
'UserData',{  'cat_vol_slice_overlay' 'y = cat_surf_results;' },...
'Tag','popupmenu062',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

%-------------------------------------------------------------------
appdata = [];
appdata.lastValidTag = 'uipanel070';

h070 = uibuttongroup(...
'Parent',h000,...
'FontSize',FS(7),...
'Title','Tools',...
'Tag','uipanel070',...
'UserData',[],...
'Clipping','on',...
'BackgroundColor',[0.8 0.8 0.8],...
'Position',[0.069 0.0267857142857143 0.853 0.095],...
'SelectedObject',[],...
'SelectionChangeFcn',[],...
'OldSelectedObject',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupmenu071';

h071 = uicontrol(...
'Parent',h070,...
'Units','normalized',...
'Callback','spm(''PopUpCB'',gcbo)',...
'FontSize',FS(6),...
'ListboxTop',0,...
'Position',[x2(1) y1p 0.65*x2w y1h],...
'String',{  'Internet/Help...     '; 'Interactive Help'; 'Manual'; 'Update'; 'VBM Website' },...
'Style','popupmenu',...
'TooltipString','Call website/manual or update VBM12',...
'Value',1,...
'UserData',{  'web(fullfile(spm(''Dir''),''toolbox'',''cat12'',''html'',''cat.html''))' 'try,open(fullfile(spm(''dir''),''toolbox'',''cat12'',''CAT12-Manual.pdf''));end' 'spm(''alert'',evalc(''cat_update(1)''),''CAT Update'');' 'try,web(''http://dbm.neuro.uni-jena.de/vbm'',''-browser'');end' },...
'Tag','popupmenu071',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'popupmenu072';

h072 = uicontrol(...
'Parent',h070,...
'Units','normalized',...
'Callback','spm(''PopUpCB'',gcbo)',...
'FontSize',FS(6),...
'ListboxTop',0,...
'Position',[0.33 y1p 0.65*x2w y1h],...
'String',{  'Utils...      '; 'Denoising Filter (SANLM)'; 'Apply Deformations (Many Subjects)'; 'Apply Deformations (Many Images)'; 'Calculate Differences between Images or Surfaces'; 'Install CAT Atlases to SPM' },...
'Style','popupmenu',...
'TooltipString','Additional tools',...
'Value',1,...
'UserData',{  'spm_jobman(''interactive'','''',''spm.tools.cat.tools.sanlm'');' 'spm_jobman(''interactive'','''',''spm.tools.cat.tools.defs2'');' 'spm_jobman(''interactive'','''',''spm.tools.cat.tools.defs'');' 'cat_stat_diff' 'cat_install_atlases' },...
'Tag','popupmenu072',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton073';

h073 = uicontrol(...
'Parent',h070,...
'Units','normalized',...
'Callback',@(hObject,eventdata)cat12('pushbutton073_Callback',hObject,eventdata,guidata(hObject)),...
'FontSize',FS(6),...
'Position',[0.6 y1b 0.38*x2w y1h],...
'String','Report Error',...
'TooltipString','Close',...
'Tag','pushbutton073',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pushbutton074';

h074 = uicontrol(...
'Parent',h070,...
'Units','normalized',...
'Callback',@(hObject,eventdata)cat12('pushbutton074_Callback',hObject,eventdata,guidata(hObject)),...
'FontSize',FS(6),...
'ForegroundColor',[1 0 0],...
'Position',[0.77 y1b 0.38*x2w y1h],...
'String','Close',...
'TooltipString','Close',...
'Tag','pushbutton074',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

%-------------------------------------------------------------------
appdata = [];
appdata.lastValidTag = 'text080';

h080 = uicontrol(...
'Parent',h000,...
'Units','normalized',...
'FontUnits','pixels',...
'BackgroundColor',[0.95 0.95 0.95],...
'FontSize',FS(9),...
'ForegroundColor',[0.6 0.6 0.6],...
'ListboxTop',0,...
'Position',[-0.00170068027210884 -0.00127551020408163 1 0.0216836734693878],...
'String','Copyright (c) Structural Brain Mapping Group',...
'Style','text',...
'Tag','text080',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );


hsingleton = h000;


% --- Set application data first then calling the CreateFcn. 
function local_CreateFcn(hObject, eventdata, createfcn, appdata)

if ~isempty(appdata)
   names = fieldnames(appdata);
   for i=1:length(names)
       name = char(names(i));
       setappdata(hObject, name, getfield(appdata,name));
   end
end

if ~isempty(createfcn)
   if isa(createfcn,'function_handle')
       createfcn(hObject, eventdata);
   else
       eval(createfcn);
   end
end


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)

gui_StateFields =  {'gui_Name'
    'gui_Singleton'
    'gui_OpeningFcn'
    'gui_OutputFcn'
    'gui_LayoutFcn'
    'gui_Callback'};
gui_Mfile = '';
for i=1:length(gui_StateFields)
    if ~isfield(gui_State, gui_StateFields{i})
        error('MATLAB:gui_mainfcn:FieldNotFound', 'Could not find field %s in the gui_State struct in GUI M-file %s', gui_StateFields{i}, gui_Mfile);
    elseif isequal(gui_StateFields{i}, 'gui_Name')
        gui_Mfile = [gui_State.(gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % CAT12
    % create the GUI only if we are not in the process of loading it
    % already
    gui_Create = true;
elseif local_isInvokeActiveXCallback(gui_State, varargin{:})
    % CAT12(ACTIVEX,...)
    vin{1} = gui_State.gui_Name;
    vin{2} = [get(varargin{1}.Peer, 'Tag'), '_', varargin{end}];
    vin{3} = varargin{1};
    vin{4} = varargin{end-1};
    vin{5} = guidata(varargin{1}.Peer);
    feval(vin{:});
    return;
elseif local_isInvokeHGCallback(gui_State, varargin{:})
    % CAT12('CALLBACK',hObject,eventData,handles,...)
    gui_Create = false;
else
    % CAT12(...)
    % create the GUI and hand varargin to the openingfcn
    gui_Create = true;
end

if ~gui_Create
    % In design time, we need to mark all components possibly created in
    % the coming callback evaluation as non-serializable. This way, they
    % will not be brought into GUIDE and not be saved in the figure file
    % when running/saving the GUI from GUIDE.
    designEval = false;
    if (numargin>1 && ishghandle(varargin{2}))
        fig = varargin{2};
        while ~isempty(fig) && ~isa(handle(fig),'figure')
            fig = get(fig,'parent');
        end
        designEval = isappdata(0,'CreatingGUIDEFigure');
    end
        
    if designEval
        beforeChildren = findall(fig);
    end
    
    % evaluate the callback now
    varargin{1} = gui_State.gui_Callback;
    if nargout
        [varargout{1:nargout}] = feval(varargin{:});
    else       
        feval(varargin{:});
    end
    
    % Set serializable of objects created in the above callback to off in
    % design time. Need to check whether figure handle is still valid in
    % case the figure is deleted during the callback dispatching.
    if designEval && ishandle(fig)
        set(setdiff(findall(fig),beforeChildren), 'Serializable','off');
    end
else
    if gui_State.gui_Singleton
        gui_SingletonOpt = 'reuse';
    else
        gui_SingletonOpt = 'new';
    end

    % Check user passing 'visible' P/V pair first so that its value can be
    % used by oepnfig to prevent flickering
    gui_Visible = 'auto';
    gui_VisibleInput = '';
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end

        % Recognize 'visible' P/V pair
        len1 = min(length('visible'),length(varargin{index}));
        len2 = min(length('off'),length(varargin{index+1}));
        if ischar(varargin{index+1}) && strncmpi(varargin{index},'visible',len1) && len2 > 1
            if strncmpi(varargin{index+1},'off',len2)
                gui_Visible = 'invisible';
                gui_VisibleInput = 'off';
            elseif strncmpi(varargin{index+1},'on',len2)
                gui_Visible = 'visible';
                gui_VisibleInput = 'on';
            end
        end
    end
    
    % Open fig file with stored settings.  Note: This executes all component
    % specific CreateFunctions with an empty HANDLES structure.

    
    % Do feval on layout code in m-file if it exists
    gui_Exported = ~isempty(gui_State.gui_LayoutFcn);
    % this application data is used to indicate the running mode of a GUIDE
    % GUI to distinguish it from the design mode of the GUI in GUIDE. it is
    % only used by actxproxy at this time.   
    setappdata(0,genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]),1);
    if gui_Exported
        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);

        % make figure invisible here so that the visibility of figure is
        % consistent in OpeningFcn in the exported GUI case
        if isempty(gui_VisibleInput)
            gui_VisibleInput = get(gui_hFigure,'Visible');
        end
        set(gui_hFigure,'Visible','off')

        % openfig (called by local_openfig below) does this for guis without
        % the LayoutFcn. Be sure to do it here so guis show up on screen.
        movegui(gui_hFigure,'onscreen');
    else
        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt, gui_Visible);
        % If the figure has InGUIInitialization it was not completely created
        % on the last pass.  Delete this handle and try again.
        if isappdata(gui_hFigure, 'InGUIInitialization')
            delete(gui_hFigure);
            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt, gui_Visible);
        end
    end
    if isappdata(0, genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]))
        rmappdata(0,genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]));
    end

    % Set flag to indicate starting GUI initialization
    setappdata(gui_hFigure,'InGUIInitialization',1);

    % Fetch GUIDE Application options
    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');
    % Singleton setting in the GUI M-file takes priority if different
    gui_Options.singleton = gui_State.gui_Singleton;

    if ~isappdata(gui_hFigure,'GUIOnScreen')
        % Adjust background color
        if gui_Options.syscolorfig
            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));
        end

        % Generate HANDLES structure and store with GUIDATA. If there is
        % user set GUI data already, keep that also.
        data = guidata(gui_hFigure);
        handles = guihandles(gui_hFigure);
        if ~isempty(handles)
            if isempty(data)
                data = handles;
            else
                names = fieldnames(handles);
                for k=1:length(names)
                    data.(char(names(k)))=handles.(char(names(k)));
                end
            end
        end
        guidata(gui_hFigure, data);
    end

    % Apply input P/V pairs other than 'visible'
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end

        len1 = min(length('visible'),length(varargin{index}));
        if ~strncmpi(varargin{index},'visible',len1)
            try set(gui_hFigure, varargin{index}, varargin{index+1}), catch break, end
        end
    end

    % If handle visibility is set to 'callback', turn it on until finished
    % with OpeningFcn
    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
    if strcmp(gui_HandleVisibility, 'callback')
        set(gui_hFigure,'HandleVisibility', 'on');
    end

    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});

    if isscalar(gui_hFigure) && ishandle(gui_hFigure)
        % Handle the default callbacks of predefined toolbar tools in this
        % GUI, if any
        guidemfile('restoreToolbarToolPredefinedCallback',gui_hFigure); 
        
        % Update handle visibility
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);

        % Call openfig again to pick up the saved visibility or apply the
        % one passed in from the P/V pairs
        if ~gui_Exported
            gui_hFigure = local_openfig(gui_State.gui_Name, 'reuse',gui_Visible);
        elseif ~isempty(gui_VisibleInput)
            set(gui_hFigure,'Visible',gui_VisibleInput);
        end
        if strcmpi(get(gui_hFigure, 'Visible'), 'on')
            figure(gui_hFigure);
            
            if gui_Options.singleton
                setappdata(gui_hFigure,'GUIOnScreen', 1);
            end
        end

        % Done with GUI initialization
        if isappdata(gui_hFigure,'InGUIInitialization')
            rmappdata(gui_hFigure,'InGUIInitialization');
        end

        % If handle visibility is set to 'callback', turn it on until
        % finished with OutputFcn
        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
        if strcmp(gui_HandleVisibility, 'callback')
            set(gui_hFigure,'HandleVisibility', 'on');
        end
        gui_Handles = guidata(gui_hFigure);
    else
        gui_Handles = [];
    end

    if nargout
        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    else
        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    end

    if isscalar(gui_hFigure) && ishandle(gui_hFigure)
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
    end
end

function gui_hFigure = local_openfig(name, singleton, visible)

% openfig with three arguments was new from R13. Try to call that first, if
% failed, try the old openfig.
if nargin('openfig') == 2
    % OPENFIG did not accept 3rd input argument until R13,
    % toggle default figure visible to prevent the figure
    % from showing up too soon.
    gui_OldDefaultVisible = get(0,'defaultFigureVisible');
    set(0,'defaultFigureVisible','off');
    gui_hFigure = openfig(name, singleton);
    set(0,'defaultFigureVisible',gui_OldDefaultVisible);
else
    gui_hFigure = openfig(name, singleton, visible);
end

function result = local_isInvokeActiveXCallback(gui_State, varargin)

try
    result = ispc && iscom(varargin{1}) ...
             && isequal(varargin{1},gcbo);
catch
    result = false;
end

function result = local_isInvokeHGCallback(gui_State, varargin)

try
    fhandle = functions(gui_State.gui_Callback);
    result = ~isempty(findstr(gui_State.gui_Name,fhandle.file)) || ...
             (ischar(varargin{1}) ...
             && isequal(ishandle(varargin{2}), 1) ...
             && (~isempty(strfind(varargin{1},[get(varargin{2}, 'Tag'), '_'])) || ...
                ~isempty(strfind(varargin{1}, '_CreateFcn'))) );
catch
    result = false;
end


