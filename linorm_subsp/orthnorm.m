function [ Uout1, Uout2 ] = orthnorm( Uin, tol )
%
% PURPOSE:
%
% To compute an orthonormal basis for the subspace spanned by the columns
% of Uin, i.e., Uout'*Uout = eye(k) for some k.
%
% ARGUMENTS:
%
% Inputs:
%
% Uin      : a matrix to be orthonormalized.
% tol      : relative truncation tolerance for the singular values of Uin.
%            
% Outputs:
%
% Uout1    : a matrix with orthonormal columns approximately spanning the
%            image of Uin. Uout may have less columns than Uin if Uin has
%            singular values which are smaller than tol*max(svd(Uin)). 
%            In this case, the corresponding parts Uin are truncated.
% Uout2    : a matrix with orthonormal columns containing the truncated
%            parts of the matrix Uin correponding the singular values 
%            smaller than tol*max(svd(Uin)).
%
% REFERENCES:
%
% None.
%
% AUTHORS:
%
% Matthias Voigt, Technische Universitaet Berlin, Institut fuer Mathematik,
% Berlin, Germany.
%
% 20/04/2017
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
%% COMPUTATIONS.
%
[ Uout, S ] = svd( Uin, 'econ' );
%
% Scale the singular values such that the maximal one is equal to 1.
%
m = size( Uout, 2 );
S = diag( S );
S = S./S( 1 );
%
% Determine the parts of Uin that correspond to values in S which are 
% larger than the truncation tolerance tol.
%
dim = min( length( S( S > tol ) ), m );
%
% Determine the final orthonormal bases.
%
Uout1 = Uout( :, 1:dim );
Uout2 = Uout( :, dim+1:m );
%
return