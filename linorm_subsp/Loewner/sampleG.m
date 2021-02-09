function data = sampleG( sys, smpl, maxDir )
%
% PURPOSE:
%
% To sample the values and first derivatives of the function G(s) at given 
% sample points, where G(s) is holomorphic in an area around the sample 
% points.
%
% ARGUMENTS:
%
% Inputs:
%            
% sys     : Function handle containing the function G(s).
% smpl    : vector in which the complex interpolation points are given.
% maxDir  : specifies whether the singular vectors of all sampled function
%           values corresponding to the maximum singular values should be 
%           computed as follows:
%           = 0: no singular vectors are computed;
%           = 1: the left and right singular vectors of the sampled
%                function values corresponding to the maximum singular
%                values are computed.
%
% Output:
%
% data   : Struct containing information about the evaluated function.
%    data.smpl : array in which the interpolation points are stored.
%    data.G    : cell array in which the function values of G(s) at the
%                sample points are stored as 
%                  data.G = { G(smpl(1)), ..., G(smpl(length(smpl))) }.
%    data.dG   : cell array in which the derivative of G(s) evaluated at
%                the sample points are stored as 
%                  data.dG = { dG/ds(smpl(1)), ...,
%                              dG/ds(smpl(length(smpl))) }.
%    data.b, 
%    data.c    : arrays containing the left and right singular vectors
%                associated with the maximum singular value of the
%                function values of G(s) at the sample points,
%                respectively. Hereby, data.b( :, k ) and data.c( :, k )
%                contain the corresponding singular vectors of G(smpl(k)).
%                These fields are only computed, if maxDir == 1.
%
% REFERENCES:
%
% [1] C. Beattie and S. Gugercin. Realization-independent 
%     H_2-approximation. Proc. 51st IEEE Conference on Decision and 
%     Control, Maui, HI, USA, 2012.
%
% AUTHORS: 
%
% Paul Schwerdtner and Matthias Voigt, Technische Universitaet Berlin,
% Institut fuer Mathematik, Berlin, Germany.
%
% 27/09/2017.
%
% REVISIONS:
%
% Matthias Voigt, 12/2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LINORM_SUBSP 1.2 Copyright ( C ) 2018 Paul Schwerdtner, Matthias Voigt
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  ( at your option ) any later version.
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
%% COMPUTATIONS.
%
% Initialize the data struct.
%
data.smpl = smpl;
data.G = cell( length( smpl ), 1 );
data.dG = cell( length( smpl ), 1 );
h = 1e-8;
%
for k = 1:length( smpl )
    %
    % Evaluate the function G(s) at the given sample points.
    %
    data.G{ k }  = sys( smpl( k ) );
    data.dG{ k } = ( sys( smpl( k ) + h ) - sys( smpl( k ) - h ) )/( 2*h );
    if maxDir == 1
        [ C, ~, B ] = svd( data.Hw{ k } );
        data.c( :, k ) = C( :, 1 );
        data.b( :, k ) = B( :, 1 );
    end
end
%
return