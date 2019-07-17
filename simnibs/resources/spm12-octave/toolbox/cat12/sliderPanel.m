function [sliderHandle,panelHandle,editHandle] = sliderPanel(parent,PanelPVs,SliderPVs,EditPVs,LabelPVs,numFormat,varargin)
% [sliderHandle,panelHandle,editHandle] = sliderPanel(parent,PanelPVs,SliderPVs,EditPVs,LabelPVs,numFormat)
%
% Creates a slider in a separate uipanel, with an associated
% interactive EditBox, and left and right labels showing the
% minimum and maximum values of the slider, respectively.
% Moving the slider automatically updates the textbox, and vice
% versa. Both slider movement and text edits will trigger
% (non-recursively) the callback of the slider.
%
% The EditBox automatically disallows the entry of
% non-numeric values, or of values outside of [min,max].
% Attempts to enter disallowed values will be ignored.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTAXES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1) The main syntax for sliderPanel allows full control of
% all elements of the uitool. (See INPUT ARGUMENTS below for
% details).
%
% [sliderHandle,panelHandle,editHandle] =
%     sliderPanel(parent,PanelPVs,SliderPVs,EditPVs,LabelPVs,numFormat);
%
% 2)The following syntax captures a small subset of the
% sliderPanel functionality;
%
% hFig = gcf;
% sliderPanel(...
% 		'Parent'  , hFig, ...
% 		'Title'   , 'Slider Panel', ...
% 		'Position', [0.3 0.5 0.4 0.2], ...
%       'Backgroundcolor', 'r',...
% 		'Min'     , 0, ...
% 		'Max'     , 100, ...
% 		'Value'   , 50, ...
% 		'FontName', 'Verdana', ...
% 		'Callback', @myCallback);
%
% Note: This simplified syntax supports a small subset of
% the available PV pairs that can be controlled via the
% primary syntax. Supported parameters, and the object to
% which they assigned, are shown below:
%
% Parameter       Object(s) affected
%____________________________________
% Parent            uipanel
% Title             uipanel
% Position          uipanel
% Backgroundcolor   uipanel
% Bordertype        uipanel
% Tag               uipanel
% Fontname          uipanel, edit box, labels
% Fontweight        uipanel, edit box, labels
% Fontsize          uipanel, edit box, labels
% Min               slider
% Max               slider
% Value             slider
% Sliderstep        slider
% Callback          slider (and, by extension, edit box)
% Units             uipanel, slider, edit box, labels
% Visible           uipanel (as a parent)
% Numberformat      edit box
%
% (The code is easily modifiable to add new PV support.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT ARGUMENTS (ALL OPTIONAL):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    PARENT:    the handle of the parent object for the
%               uipanel. Default is the current figure.
%
%    PANELPVS:  a cell array of any parameter-value pairs
%               valid for the uipanel object, or a structure
%               of same. (Defaults are those used by
%               UIPANEL.)
%
%    SLIDERPVS: a cell array of any parameter-value pairs
%               valid for the slider object, or a structure
%               of same. (Defaults are those used by
%               UICONTROL('STYLE','SLIDER').)
%
%    EDITPVS:   a cell array of any parameter-value pairs
%               valid for the UICONTROL EDIT object, or a
%               structure of same. (Defaults are those used
%               by UICONTROL('STYLE','EDIT').)
%
%    LABELPVS:  a cell array of any parameter-value pairs
%               valid for the UICONTROL TEXT objects used to
%               label the slider, or a structure of same.
%               (Defaults are those used by
%               UICONTROL('STYLE','TEXT').)
%               NOTE: The following are valid syntaxes for LABELPVS:
%               1) Applied to both (l/r) labels:
%                  {param1, val1, param2, val2,...}
%               2) First array applied to left label, second
%                  array applied to right label:
%                  {{P-V array},{P-V array}}
%               3) First array applied to left label, second
%                  array applied to right label, third array
%                  applied to both (l/r) labels:
%                  {{P-V array},{P-V array},{P-V array}}
%
%    NUMFORMAT: a format string accepted by SPRINTF, which
%               controls the formatting of the EditBox when
%               the slider is dragged. (Typing directly in
%               the EditBox does not trigger the formatting
%               constraint.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT ARGUMENTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    SLDRHNDL:  handles of slider object.
%
%    PNLHNDL:   handle of uipanel object.
%
%    EDITHNDL:  handle of edit box.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Create a single-slider sliderPanel using cell-array inputs.
%
% [sldr,pnl,edt] = ...
%                   sliderPanel(gcf,...
%                   {'title','Threshold','pos',[0.1 0.2 0.8 0.15],'fontweight','b','units','pixels'},...
%                   {'max',80,'value',60,'callback','disp(''Slid'')'},...
%                   {},...
%                   {{'string','Low','foregroundcolor','b'},...
%                      {'string','High','foregroundcolor','r'},...
%                      {'fontweight','b'}},...
%                   '%0.1f');
%
% 2) Create a single-slider sliderPanel as a child of a
% uipanel, using struct and cell-array inputs. Moving the
% slider or updating the edit box will immediately refresh
% the value of variable sliderVal in the base workspace.
%
% hFig = figure;
% hPanel = uipanel(hFig,'title','MASTER','pos',[0.3 0.1 0.4 0.6]);
% PnlOpt.title = 'Threshold';
% PnlOpt.position = [0.05 0.2 0.9 0.4];
% PnlOpt.fontsize = 10;
% SldrOpt.min = 10;
% SldrOpt.max = 100;
% SldrOpt.value = 50;
% SldrOpt.callback = 'assignin(''base'',''sliderVal'',get(gcbo,''value''));';
% sliderPanel(hPanel,PnlOpt,SldrOpt,{'fontsize',12},{},'%0.0f')
%
% 3) Create multiple sliderPanels as children of a UIPANEL.
%
% figure;
% h = uipanel(gcf,'title','MAIN','units','normalized','pos',[0.2 0.1 0.6 0.8]);
% PnlOpt.title = 'Parameter Tool';
% PnlOpt.bordertype = 'none';
% PnlOpt.titleposition = 'centertop';
% PnlOpt.fontweight = 'bold';
% SldrOpt.min = 0;
% SldrOpt.max = 255;
% SldrOpt.value = 50;
% EditOpts = {'fontsize',10};
% LabelOpts = {'fontsize',9,'fontweight','b'};
% numFormat = '%0.0f';
% titleStrings = {'Slider 1','Slider 2', 'Slider 3', 'Slider 4'};
% startPos = {[0.1 0.05 0.8 0.21];
% 	        [0.1 0.28 0.8 0.21];
% 	        [0.1 0.51 0.8 0.21];
% 	        [0.1 0.74 0.8 0.21]};
% 		sldrCallbacks = {'disp(''Slider 1 moved'')';
% 			'disp(''Slider 2 moved'')';
% 			'disp(''Slider 3 moved'')';
% 			'disp(''Slider 4 moved'')'};
% for ii = 1:4
% 	PnlOpt.position = startPos{ii};
% 	PnlOpt.title = titleStrings{ii};
% 	SldrOpt.callback = sldrCallbacks{ii};
% 	sliderPanel(h,PnlOpt,SldrOpt,EditOpts,LabelOpts,numFormat);
% end
%
% 4) Demonstrate the use of sliderPanel in a function that
% interactively thresholds an image. (To try this example,
% paste the following code (both TEST and THRESH) into a new
% mfile, save it, and run it.)
%
% function test
% figure;
% ax1 = axes('units','normalized','pos',[0.1 0.25 0.8 0.7]);
% myimg = im2double(imread('cameraman.tif'));
% imobj = imshow(myimg);
% sliderPanel(gcf,{'pos',[0.1 0.05 0.8 0.15]},{'callback',{@thresh,myimg,imobj}},{},{},'%0.1f')
%
% % SUBFUNCTION:
% function thresh(varargin)
% sldr = varargin{1};
% newval = get(sldr,'value');
% myimg = varargin{3};
% imobj = varargin{4};
% set(imobj,'cdata',myimg > newval);
%
% 5) SIMPLIFIED SYNTAX
%
% hFig = gcf;
% [a,b]=sliderPanel(...
% 'Title'   , 'Slider Panel', ...
%   'Position', [0.3, 0.5, 0.4, 0.2], ...
%   'Min'     , 0, ...
%   'Max'     , 100, ...
%   'Value'   , 50, ...
%   'String'  , 'Slider 1', ...
%   'FontName', 'Verdana', ...
%   'units','pixels',...
%   'numformat','%0.0f',...
%   'Callback', 'disp(''slid'')')

% REVISIONS:
% 06/07/2010
% Modified simple syntax to specify label colors the same as background
% colors, unless otherwise set.
% 
% 10/10/2012
% Implemented right-click resetting to default (initial) value. Right-click
% anywhere on the slider itself, and the slider and text reset
% automatically.
%
% 1/18/2013 Right-clicking to reset default now (appropriately) triggers
% slider's callback.

% Written by Brett Shoelson, PhD
% brett.shoelson@mathworks.com
% 01/21/07
% Copyright 2007 - 2012 MathWorks, Inc.


if nargin < 6, numFormat = []; end
if nargin < 5, LabelPVs = {}; end
if nargin < 4, EditPVs = {}; end
if nargin < 3, SliderPVs = {}; end
if nargin < 2, PanelPVs = {}; end
if nargin == 0, parent = gcf; end

% OPTIONAL CALLING SYNTAX:
% Subset of functionality available for this calling syntax
if nargin > 1 && ~iscell(PanelPVs) && ~isa(PanelPVs,'struct')
    allargs = {parent,PanelPVs,SliderPVs,EditPVs,LabelPVs,numFormat,varargin{:}};
    PanelPVs = {};SliderPVs = {};EditPVs = {};LabelPVs = {};numFormat = [];
    if ishandle(allargs{1}),
        parent = allargs{1};
    end
    loc = cellfind(allargs,'parent');
    if ~isempty(loc)
        parent = allargs{loc+1};
    end
    if ~ishandle(parent)
        parent = gcf;
    end
    loc = cellfind(allargs,'numformat');
    if ~isempty(loc)
        numFormat = allargs{loc+1};
    end
    loc = cellfind(allargs,'numberformat');
    if ~isempty(loc)
        numFormat = allargs{loc+1};
    end

    PanelPVs = validate(allargs,PanelPVs,...
        {'title','pos','position','fontname','fontname','fontsize',...
        'fontweight','backgroundcolor','bordertype','tag','units','visible','userdata'});
    SliderPVs = validate(allargs,SliderPVs,...
        {'min','max','value','callback','sliderstep','units'});
    EditPVs = validate(allargs,EditPVs,...
        {'fontname','fontsize','fontweight'});
    LabelPVs = validate(allargs,LabelPVs,{'fontname','fontsize',...
        'fontweight','units'});
    if isfield(PanelPVs,'backgroundcolor') && ~isfield(LabelPVs,'backgroundcolor')
        LabelPVs.backgroundcolor = PanelPVs.backgroundcolor;
    end
end

% CREATE PANEL (DEFAULT)
% Uipanels can be children of figures, uipanels, or
% uibuttongroups; the latter two are of 'type' 'uipanel'

% CREATE UIPANEL AS PARENT
panelHandle = uipanel('parent',parent);

% CREATE SLIDER (DEFAULT)
sliderHandle = uicontrol(panelHandle,'style','slider','units','normalized',...
    'pos',[0.05 0.5 0.9 0.40]);
% CREATE EDIT BOX (DEFAULT)
editHandle = uicontrol(panelHandle,'style','edit','units','normalized',...
    'pos',[0.35 0.05 0.3 0.35]);
% CREATE LABELS (DEFAULT)
labelHandle(1) = uicontrol(panelHandle,'style','text','units','normalized',...
    'pos',[0.05 0.05 0.25 0.25],'horizontalalignment','l',...
    'Backgroundcolor',get(panelHandle,'Backgroundcolor'));
labelHandle(2) = uicontrol(panelHandle,'style','text','units','normalized',...
    'pos',[0.7 0.05 0.25 0.25],'horizontalalignment','r',...
    'Backgroundcolor',get(panelHandle,'Backgroundcolor'));

% CUSTOMIZE PER USER-DEFINED PV PAIRS
applyPVs(panelHandle,PanelPVs);
applyPVs(sliderHandle,SliderPVs);

% EXTRACT SLIDER PARAMETERS
sldr.minval = get(sliderHandle,'min');
sldr.maxval = get(sliderHandle,'max');
sldr.value = get(sliderHandle,'value');
%Ensure that slider's value is in  acceptable range of
%[min,max]
if sldr.value < sldr.minval || sldr.value > sldr.maxval
    disp('Out-of-range value ignored for slider object.');
    sldr.value = sldr.minval;
    set(sliderHandle,'value',sldr.value);
end
sldr.callback = get(sliderHandle,'callback');

set(labelHandle(1),'string',sldr.minval);
set(labelHandle(2),'string',sldr.maxval);
set(editHandle,'string',sldr.value);

applyPVs(editHandle,EditPVs);
% If LabelPVs is a 1x2 array of cells, then cell 1 is
% applied to label 1, and cell 2 is applied to label 2.
% Otherwise,
if numel(LabelPVs)>1 && iscell(LabelPVs{1}) && iscell(LabelPVs{2})
    applyPVs(labelHandle(1),LabelPVs{1});
    applyPVs(labelHandle(2),LabelPVs{2});
else
    applyPVs(labelHandle,LabelPVs);
end
% If a third array of PVs is provided, use it for both
% labels.
if numel(LabelPVs)>2 && iscell(LabelPVs{3})
    applyPVs(labelHandle,LabelPVs{3});
end

%USERDATA IS A CHAR
set(editHandle,'userdata',get(editHandle,'string'));

%GIVE ADDHNDLEVENT/EVALHNDLEVENT A TRY
addHndlEvent(sliderHandle,'callback',@updateText);
addHndlEvent(editHandle,'callback',@updateSlider);

    function applyPVs(obj,pvarray)
        if isstruct(pvarray)
            set(obj,pvarray);
        else %Cell
            if ~isempty(pvarray)
                for ii = 1:2:numel(pvarray)
                    set(obj,pvarray{ii},pvarray{ii+1});
                end
                %set(obj,pvarray{:});
            end
        end
        % NEW: 10/10/2012: implement right-click resetting of default
        %      initial value
        if strcmp(get(obj,'type'),'uicontrol')
            if strcmp(get(obj,'style'),'slider')
            val = get(obj,'value');
            set(obj,'buttondownfcn',{@resetDefault,gcbo,val});
            end
        end
    end

    function resetDefault(varargin)
        set(varargin{1},'value',varargin{4});
        updateText;
        if ~isempty(sldr.callback)
            feval(sldr.callback,varargin{1});
        end
    end

    function updateText(varargin)
        %Triggered by slider move
        newVal = get(sliderHandle,'value');
        if ~isempty(numFormat)
            newVal = sprintf(numFormat,newVal);
        end
        set(editHandle,'string',newVal,'userdata',newVal);
        drawnow;
    end

    function updateSlider(varargin)
        %Triggered by text change
        newVal = get(editHandle,'string');
        if isnan(str2double(newVal)) || str2double(newVal) < get(sliderHandle,'min') || str2double(newVal) > get(sliderHandle,'max')
            set(editHandle,'string',get(editHandle,'userdata'));
            return
        end
        set(sliderHandle,'value',str2double(newVal));
        set(editHandle,'userdata',newVal);
        set(editHandle,'value',str2double(newVal));%BUG FIX, 7/1/10
        drawnow;
        currFun = sldr.callback;
        if isempty(currFun)
            %Do nothing
        elseif isa(currFun,'function_handle')
            %REPLACE OBJECT HANDLE WITH THAT OF SLIDER
            varargin{1} = sliderHandle;
            currFun(varargin{:});
            %feval(sldr.callback,varargin{:});
        elseif isa(currFun,'char')
            %Just in case the callback specifies GCBO (which
            %will initially point to the edit box, rather
            %than the slider:
            currFun = strrep(currFun,'gcbo','sliderHandle');
            eval(currFun);
        elseif isa(currFun,'cell');
            %REPLACE OBJECT HANDLE WITH THAT OF SLIDER
            varargin{1} = sliderHandle;
            currFun{1}(varargin{:},currFun{2:end});
        else
            %Shouldn't ever get here...but if you do, I
            %would like to know about it.
            error('Unrecognized event registered in eventid %d.',ii);
        end
    end

    function addHndlEvent(hndl,eventtype,newevent)
        hndlevent = getappdata(hndl,[eventtype 'hndlevent']);
        if isempty(hndlevent)
            %Initialize to original event comand
            hndlevent.cmdset = get(hndl,eventtype);
        end
        if isempty(hndlevent.cmdset)
            numevents = 0;
        else
            numevents = numel(hndlevent);
        end
        hndlevent(numevents+1).cmdset = newevent;

        set(hndl,eventtype,@evalHndlEvent);
        setappdata(hndl,[eventtype 'hndlevent'],hndlevent);

        function evalHndlEvent(varargin)
            hndlevent = getappdata(hndl,[eventtype 'hndlevent']);
            if nargin < 4
                eventList = 1:numel(hndlevent);
            end
            for ii = eventList
                currFun = hndlevent(ii).cmdset;
                if isempty(currFun)
                    continue
                elseif ischar(currFun)
                    eval(currFun);
                elseif iscell(currFun)
                    currFun{1}(varargin{:},currFun{2:end});
                else
                    currFun(varargin{:});
                end
            end
        end
    end

    function PVarray = validate(allargs, PVarray, parameterStrings)
        for ii = 1:numel(parameterStrings)
            parameter = parameterStrings{ii};
            loc = cellfind(allargs, parameter);
            if ~isempty(loc)
                eval(['PVarray.' parameter ' = allargs{loc(1)+1};']);
            end
        end
    end

    function posns = cellfind(cellarray, searchval)
        posns = [];
        if ischar(searchval)
            searchval = lower(searchval);
            for ii = 1:numel(cellarray)
                tmp = cellarray{ii};
                if ischar(tmp)
                    tmp = lower(tmp);
                end
                if isequal(searchval,tmp)
                    posns = [posns;ii];
                end
            end
        end
    end

end