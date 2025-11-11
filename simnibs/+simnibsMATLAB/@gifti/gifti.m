classdef gifti
% GIfTI Geometry file format class
% Geometry format under the Neuroimaging Informatics Technology Initiative
% (NIfTI):
%                 http://www.nitrc.org/projects/gifti/
%                      http://nifti.nimh.nih.gov/
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

    properties

        metadata (1,1) struct = struct % TODO: short documentation here.

        label (1,1) struct = struct % TODO: short documentation here.

        data (1,:) cell = struct % TODO: short documentation here.

        mat % TODO: add size, class and default value restrictions here.

        faces % TODO: add size, class and default value restrictions here.

        vertices % TODO: add size, class and default value restrictions here.

        cdata (1,1) struct = struct % TODO: short documentation here.

        indices % TODO: add size, class and default value restrictions here.

    end % properties

    methods

        function this = gifti(varargin)
        % GIfTI Geometry file format class
        % Geometry format under the Neuroimaging Informatics Technology Initiative
        % (NIfTI):
        %                 http://www.nitrc.org/projects/gifti/
        %                      http://nifti.nimh.nih.gov/
        %__________________________________________________________________________
        % Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

        % Guillaume Flandin
        % $Id: gifti.m 7621 2019-06-20 16:58:59Z guillaume $

            switch nargin

                case 0

                    this = simnibsMATLAB.gifti(giftistruct) ;

                case 1
                    if isa(varargin{1},'simnibsMATLAB.gifti')
                        this = varargin{1};

                    elseif isstruct(varargin{1})
                        f       = {'faces', 'face', 'tri' 'vertices', 'vert', 'pnt', 'cdata', 'indices'};
                        ff      = {'faces', 'faces', 'faces', 'vertices', 'vertices', 'vertices', 'cdata', 'indices'};
                        [c, ia] = intersect(f,fieldnames(varargin{1}));
                        if ~isempty(c)

                            for i=1:length(c)
                                fieldName = ff{ia(i)} ;
                                this.(fieldName) = varargin{1}.(c{i}) ;
                            end
                            if isfield(varargin{1},'mat')
                                this.mat = varargin{1}.mat ;
                            end
                        elseif isempty(setxor(fieldnames(varargin{1}), {'metadata','label','data'}))
                            this = this ;
                        else
                            error('[GIFTI] Invalid structure.');
                        end

                    elseif ishandle(varargin{1})

                        this = struct('vertices',get(varargin{1},'Vertices'), ...
                                      'faces',   get(varargin{1},'Faces'));
                        if ~isempty(get(varargin{1},'FaceVertexCData'))
                              this.cdata = get(varargin{1},'FaceVertexCData');
                        end
                        this = simnibsMATLAB.gifti(this);

                    elseif isnumeric(varargin{1})

                        this = simnibsMATLAB.gifti;
                        this.cdata = varargin{1};

                    elseif iscell(varargin{1}) && numel(varargin{1}) == 1 && isnumeric(varargin{1}{1})
                        this = simnibsMATLAB.gifti;
                        for i=1:size(varargin{1}{1},2)
                            this.data{i}.metadata = struct([]);
                            this.data{i}.space    = [];
                            this.data{i}.attributes.Intent = 'NIFTI_INTENT_NONE';
                            this.data{i}.attributes.DataType = 'NIFTI_TYPE_FLOAT32';
                            this.data{i}.attributes.Dim = size(varargin{1}{1},1);
                            this.data{i}.data     = single(varargin{1}{1}(:,i));
                        end

                    elseif ischar(varargin{1})

                        if size(varargin{1},1)>1
                            this = simnibsMATLAB.gifti(cellstr(varargin{1}));
                            return;
                        end
                        [~,~,e] = fileparts(varargin{1});
                        if strcmpi(e,'.mat')
                            try
                                this = simnibsMATLAB.gifti(load(varargin{1}));
                            catch
                                error('[GIFTI] Loading of file %s failed.', varargin{1});
                            end
                        elseif ismember(lower(e),{'.asc','.srf','.mgh','.mgz','.pial',...
                                '.white','.inflated','.nofix','.orig','.smoothwm',...
                                '.sphere','.reg','.surf','.curv','.area','.sulc','.annot'})
                            this = freesurfer_read(varargin{1});
                            this = simnibsMATLAB.gifti(this);
                        elseif strcmpi(e,'.vtk')
                            this = mvtk_read(varargin{1});
                            this = simnibsMATLAB.gifti(this);
                        elseif strcmpi(e,'.obj')
                            this = obj_read(varargin{1});
                            this = simnibsMATLAB.gifti(this);
                        elseif strcmpi(e,'.ply')
                            this = ply_read(varargin{1});
                            this = simnibsMATLAB.gifti(this);
                        elseif strcmpi(e,'.off')
                            this = off_read(varargin{1});
                            this = simnibsMATLAB.gifti(this);
                        elseif strcmpi(e,'.stl')
                            this = stl_read(varargin{1});
                            this = simnibsMATLAB.gifti(this);
                        elseif strcmpi(e,'.mz3')
                            this = mz3_read(varargin{1});
                            this = simnibsMATLAB.gifti(this);
                        else
                            this = gifti_read(varargin{1},giftistruct);
                            this = simnibsMATLAB.gifti(this);
                        end

                    elseif iscellstr(varargin{1})
                        fnames = varargin{1};
                        this(numel(fnames)) = simnibsMATLAB.gifti(giftistruct);
                        for i=1:numel(fnames)
                            this(i) = simnibsMATLAB.gifti(fnames{i});
                        end

                    else
                        error('[GIFTI] Invalid object construction.');
                    end

                otherwise
                    error('[GIFTI] Invalid object construction.');
            end % switch
        end % function
    end % methods
end % classdef

%==========================================================================
function s = giftistruct
s = struct(...
    'metadata', ...
        struct(...
            'name',       {}, ...
            'value',      {} ...
        ), ...
    'label', ...
        struct(...
            'name',       {}, ...
            'index',      {} ...
        ), ...
    'data', ...
        struct(...
            'attributes', {}, ...
            'metadata',   struct('name',{}, 'value',{}), ...
            'space',      {}, ...
            'data',       {} ...
        ) ...
    );

end % function
