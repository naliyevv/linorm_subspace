function data = sampleFun( sys, smpl, maxDir )
%
%% PURPOSE:
%
% To sample the values and first derivatives of the function 
%
%   G(s) = C(s)*(s*E(s) - A(s))^(-1)*B(s) + D(s)                        (1)
%
% at given sample points, where G(s) is holomorphic in an area around the
% sample points and A(s), B(s), C(s), D(s), and E(s) are functions that are
% meromorphic in area around the sample points. These functions are of the
% form
%
%   A(s) = a_1(s)*A_1 + ... + a_{k_a}(s)*A_{k_a},
%   B(s) = b_1(s)*B_1 + ... + b_{k_b}(s)*B_{k_b},
%   C(s) = c_1(s)*C_1 + ... + c_{k_c}(s)*C_{k_c},                       (2)
%   D(s) = d_1(s)*D_1 + ... + d_{k_d}(s)*D_{k_d},
%   E(s) = e_1(s)*E_1 + ... + e_{k_e}(s)*E_{k_e},
%
% where all a_j, b_j, c_j, d_j, and e_j are scalar-valued functions which
% are meromorphic in an area enclosing the imaginary axis, and all A_j,
% B_j, C_j, D_j, and E_j are fixed matrices. Here, all A_j and E_j are 
% assumed to be sparse.
%
% ARGUMENTS:
%
% Inputs:
%            
% sys     : Struct containing the structural information of the function
%           G(s).
%    sys.A, sys.B, sys.C,
%    sys.D, sys.E         : arrays containing the matrices A_j, B_j, C_j,
%                           D_j, and E_j in (2), respectively. If k in k_x
%                           is larger than 1, i.e., the number of summands
%                           within the function X(s) is larger than 1 and
%                           multiple constant matrices X_1 to X_{k_x} must 
%                           be given, they must be stored in a cell array 
%                           as sys.X = { X_1, ..., X_{k_x} }, where x is 
%                           either a, b, c, d, or e and X is either A, B,
%                           C, D, or E.  
%    sys.fct.type         : specifies the type of the functions a_j, b_j,
%                           c_j, d_j, e_j in (2) as follows:
%                           = 'd': a_j, b_j, c_j, d_j, e_j are all delay 
%                                  functions, i.e., of the form 
%                                     x_j(s) = exp(-i*s*tau_x_j),       (3)
%                                  where x is either a, b, c, d, or e;
%                           = 'l': a_j, b_j, c_j, d_j, e_j are general 
%                                  functions.
%    sys.fct.a, sys.fct.b,
%    sys.fct.c, sys.fct.d,
%    sys.fct.e            : vectors specifying the functions a_j, b_j, c_j,
%                           d_j, and e_j, respectively.
%                           If sys.fct.type == 'd', then they are
%                           represented as sys.fct.x(j) = tau_x_j, where
%                           tau_x_j is as in (3) and x is either a, b, c,
%                           d, or e.
%                           If sys.fct.type == 'l', then sys.fct.x is a
%                           function handle that returns the evaluation of
%                           all x_j(s) as in (2), where x is either a, b,
%                           c, d, or e. That is 
%                           sys.fct.x(s) = [ x_1(s); ...; x_{k_x}(s) ].
%    sys.des              : specifies whether the function G(s) can be
%                           realized by a linear (descriptor) system as
%                           follows:
%                           = 0: the function is not realized by a linear
%                                system;
%                           = 1: the function is realized by a linear
%                                system.
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
% -
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
%
for k = 1:length( smpl )
    %
    % Evaluate the G(s) in (1) at the given sample points.
    %
    [ data.G{ k }, data.dG{ k } ] = evalFun( sys, smpl( k ) );
    if maxDir == 1
        [ C, ~, B ] = svd( data.Hw{ k } );
        data.c( :, k ) = C( :, 1 );
        data.b( :, k ) = B( :, 1 );
    end
end
%
return