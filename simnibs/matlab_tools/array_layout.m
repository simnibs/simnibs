function S = array_layout(type)
%
% create an empty data structure
% to set up a simnibs simulation
%
% S = array_layout(type)
% 
% K. Weise, 2023

%    This program is part of the SimNIBS package.
%    Please check on www.simnibs.org how to cite our work in publications.
%
%    Copyright (C) 2013-2023 Axel Thielscher, Andre Antunes, Guilherme Saturnino
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.




validtypes={'ElectrodeInitializer'};

if ~any(strcmp(type,validtypes))
    disp(validtypes)
    error('structure type has to be one of the above')
end

S.type=type;

switch S.type
       
    case 'ElectrodeInitializer'
        % general properties
        S.type = '';
        S.current = -1;

        % CircularArray properties
        S.radius_inner = -1;
        S.radius_inner_bounds = -1
        S.radius_outer = -1
        S.radius_outer_bounds = -1
        S.distance = -1
        S.distance_bounds = -1
        S.n_outer = -1
        S.n_outer_bounds = -1

        % ElectrodeArrayPair properties
        S._center = -1
        S._radius = -1
        S.radius_bounds = -1
        S._length_x = -1
        S.length_x_bounds = -1
        S._length_y = -1
        S.length_y_bounds = -1

        % simulation properties
        S.dirichlet_correction = true
        S.dirichlet_correction_detailed = false
        S.current_estimator_method = -1

    otherwise
        error('unknown structure type')
end
