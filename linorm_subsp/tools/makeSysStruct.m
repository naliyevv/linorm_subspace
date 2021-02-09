function sys = makeSysStruct( E, A, B, C, D )
%
% PURPOSE:
%
% To write the realization matrices defining a linear (descriptor) system
% into a struct that can be used as input to linorm_subsp.
%
% ARGUMENTS:
%
% Inputs:
%
% E, A, B, 
% C, D    : Matrices defining a linear decriptor system
%             Ex' = Ax + Bu, y = Cx + Du.
%
% Outputs:
%   
% sys     : Struct containing the input data to be useable in linorm_subsp.
%    sys.E                : = E.
%    sys.A                : = A.
%    sys.B                : = B.
%    sys.C                : = C.
%    sys.D                : = D.
%
% REFERENCES:
%
% None.
%
% AUTHORS:
%
% Paul Schwerdtner and Matthias Voigt, Technische Universitaet Berlin,
% Institut fuer Mathematik, Berlin, Germany.
%
% 24/04/2017.
%
% REVISIONS:
%
% Matthias Voigt, 09/2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LINORM_SUBSP 1.2 Copyright (C) 2018 Nicat Aliyev, Emre Mengi,
%                                      Paul Schwerdtner, Matthias Voigt 
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% COPY THE DATA.
%
sys.A = A;
sys.B = B;
sys.C = C;
sys.D = D;
sys.E = E;
%
return