function result = run_simnibs(session)

%Runs simnibs, given a session from sim_struct

% Guilherme Saturnino, 2018

%    This program is part of the SimNIBS package.
%    Please check on www.simnibs.org how to cite our work in publications.
%
%    Copyright (C) 2013-2018 Axel Thielscher, Andre Antunes, Guilherme Saturnino
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

name_mat = [tempname,'.mat'];
save('-v7', name_mat, 'session');
[status,result] = system([simnibs_cli_call('run_simnibs') ' "' name_mat '"']);

if status ~= 0
    error('There was an error running SimNIBS:\n %s',result)
end
