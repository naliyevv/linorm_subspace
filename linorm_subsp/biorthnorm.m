function [ Uout, Vout ] = biorthnorm( Uin, Vin, tol )
%
% PURPOSE:
%
% To compute biorthonormal bases for the subspaces spanned by the columns
% of Uin and Vin, i.e., Uout'*Vout = eye(k) for some k. It is assumend that
% Uin and Vin are already both spanned by orthonormal columns.
%
% ARGUMENTS:
%
% Inputs:
%
% Uin, Vin : two matrices with orthonormal columns.
% tol      : relative truncation tolerance for the singular values of the
%            product Uin'*Vin.
%
% Outputs:
%
% Uout, 
% Vout     : two matrices with biorthonormal columns. Uout and Vout may
%            have less columns than Uin and Vin if Uin'*Vin has singular
%            values which are smaller than tol*max(svd(Uin'*Vin)). In this
%            case, the corresponding columns in Uin and Vin are truncated.
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
% 20/04/2017
%
% REVISIONS:
%
% -
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
%% COMPUTATIONS.
%
[ L, S, R ] = svd( Uin'*Vin );
%
% Scale the singular values such that the maximal one is equal to 1.
%
Srel = diag( S )./S( 1 );
%
% Determine the parts of Uin and Vin that correspond to values in Srel
% which are larger than the truncation tolerance tol.
%
[ m, ~ ] = size( Srel );
Srel = [ Srel; 0 ];
dim = 0;
j = 1;
while( Srel( j ) > tol && j <= m )
    dim = dim + 1;
    j = j + 1;
end
S = S( 1:dim, 1:dim );
L = L( :, 1:dim );
R = R( :, 1:dim );
%
% Compute the biorthonormal bases.
%
Uout = Uin*L*diag( 1./sqrt( diag( S ) ) );
Vout = Vin*R*diag( 1./sqrt( diag( S ) ) );
%
return