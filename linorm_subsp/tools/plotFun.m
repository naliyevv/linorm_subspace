function funvals = plotFun( sys, wmin, wmax, nPoints )
%
% PURPOSE:
%
% To evaluate the maximum singular values of the function 
%
%   G(s) = C(s)*(s*E(s) - A(s))^(-1)*B(s) + D(s)                        (1)
%
% on nPoints equidistantly spaced points on the imaginary interval 
% I := [i*wmin, i*wmax]. Here we assume that G(s) is holomorphic in an open
% domain around I and that that A(s), B(s), C(s), D(s), and E(s) are 
% functions that are meromorphic in an area enclosing I. These function
% s are of the form
%
%   A(s) = a_1(s)*A_1 + ... + a_{k_a}(s)*A_{k_a},
%   B(s) = b_1(s)*B_1 + ... + b_{k_b}(s)*B_{k_b},
%   C(s) = c_1(s)*C_1 + ... + c_{k_c}(s)*C_{k_c},                       (2)
%   D(s) = d_1(s)*D_1 + ... + d_{k_d}(s)*D_{k_d},
%   E(s) = e_1(s)*E_1 + ... + e_{k_e}(s)*E_{k_e},
%
% where all a_j, b_j, c_j, d_j, and e_j are scalar-valued functions which
% are meromorphic in an area enclosing the imaginary axis, and all A_j,
% B_j, C_j, D_j, and E_j are fixed matrices.
%
% Arguments:
%
% Inputs:
%
% sys     : Struct containing the structural information of the the
%           function G(s).
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
%                                     x_j(s) = exp(-s*tau_j),         (3)
%                                  where x is either a, b, c, d, or e;
%                           = 'l': a_j, b_j, c_j, d_j, e_j are general 
%                                  functions.
%    sys.fct.a, sys.fct.b,
%    sys.fct.c, sys.fct.d,
%    sys.fct.e            : vectors specifying the functions a_j, b_j, c_j,
%                           d_j, and e_j, respectively.
%                           If sys.fct.type == 'd', then they are
%                           represented as sys.fct.x(j) = tau_j, where
%                           tau_j is as in (3) and x is either a, b, c, d,
%                           or e.
%                           If sys.fct.type == 'l', then sys.fct.x is a
%                           function handle that returns the evaluation of
%                           all x_j(s) as in (2), where x is either a, b,
%                           c, d, or e. That is 
%                           sys.fct.x(s) = [ x_1(s); ...; x_{k_x}(s) ].
%    sys.des              : specifies whether the function G(s) can is
%                           realized by a linear (descriptor) system as
%                           follows:
%                           = 0: the function is not realized by a linear
%                                system;
%                           = 1: the function is realized by a linear
%                                system.        
% wmin    : lower bound on the interval in which the function is evaluated.
% wmax    : upper bound on the interval in which the function is evaluated.
% nPoints : number of evaluation points.
%
% Output:
% 
% funvals : array containing the evaluations of the maximum singular values
%           of the function G(s) in (1) on the imaginary interval
%           I := [i*wmin, i*wmax].
%
% REFERENCES:
%
% -
%
% AUTHORS:
%
% Paul Schwerdtner and Matthias Voigt, Technische Universitaet Berlin,
% Institut fuer Mathematik, Berlin, Germany.
%
% 29/09/2017.
%
% REVISIONS:
%
% Matthias Voigt, 06/2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LINORM_SUBSP 1.2 Copyright (C) 2018 Paul Schwerdtner, Matthias Voigt     
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
%% CHECK ARGUMENTS.
% 
% Since the function evalFUN checks the arguments, they are not checked
% here. Just for convenience we add a zero matrix for D, if it is not
% given in sys.
%
if not( isfield( sys, 'D' ) )
    %
    % The matrix D is not defined in sys.
    %
    sys.D = zeros( size( sys.C, 1 ), size( sys.B,2 ) );
end
%
%% COMPUTATIONS.
%
w = linspace( wmin, wmax, nPoints );
funvals = zeros( length( w ), 1 );
for iPlot = 1:nPoints
    %
    % Compute evaluations of the transfer function.
    %
    funvals( iPlot ) = max( svd( evalFun( sys, 1i*w( iPlot ) ) ) );
end
%
% Plot the function values over the interval [wmin, wmax].
%
figure;
plot( w, funvals, '-b', 'LineWidth', 2 );
title( 'Maximum singular value plot' );
xlabel( 'frequency' );
ylabel( 'maximum singular value' );
%
return
