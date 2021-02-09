function redsys = reduceSystem( sys, U, V, opt )
%
% PURPOSE:
%
% To construct the matrices representing the reduced function Gr(s) 
% obtained from a function G(s) by projecting on the spaces represented by
% the matrices U and V.    
% The function G(s) is given by
%
%   G(s) = C(s)*(s*E(s) - A(s))^(-1)*B(s),
%
% where A(s), B(s), C(s), and E(s) are matrix-valued functions. These 
% functions are of the form
%
%   A(s) = a_1(s)*A_1 + ... + a_{k_a}(s)*A_{k_a},
%   B(s) = b_1(s)*B_1 + ... + b_{k_b}(s)*B_{k_b},                       (1)
%   C(s) = c_1(s)*C_1 + ... + c_{k_c}(s)*C_{k_c},
%   E(s) = e_1(s)*E_1 + ... + e_{k_e}(s)*E_{k_e},
%
% where all a_j, b_j, c_j, and e_j are scalar-valued functions, and all 
% A_j, B_j, C_j, and E_j are fixed matrices. Here, all A_j and E_j are 
% assumed to be sparse. The reduced function Gr(s) is given by
%
%   Gr(s) = Cr(s)*(s*Er(s) - Ar(s))^(-1)*Br(s),   
%
% where 
%
%   Ar(s) = a_1(s)*U'*A_1*V + ... + a_{k_a}(s)*U'*A_{k_a}*V,
%   Br(s) = b_1(s)*U'*B_1   + ... + b_{k_b}(s)*U'*B_{k_b},              (2)
%   Cr(s) = c_1(s)*   C_1*V + ... + c_{k_c}(s)*   C_{k_c}*V,                     
%   Er(s) = e_1(s)*U'*E_1*V + ... + e_{k_e}(s)*U'*E_{k_e}*V.
%
% ARGUMENTS:
%
% Inputs:
%
% sys     : Struct containing the structural information of the function
%           G(s)
%    sys.A, sys.B, sys.C,
%    sys.E                : arrays containing the matrices A_j, B_j, C_j,
%                           and E_j in (1), respectively. If k in k_x is
%                           larger than 1, i.e., the number of summands 
%                           within the function X(s) is larger than 1 and 
%                           multiple constant matrices X_1 to X_{k_x} need 
%                           to be given, they must be stored in a cell 
%                           array as sys.X = { X_1, ..., X_{k_x} }, where
%                           x is either a, b, c, e and X is either A, B, C,
%                           or E.  
%    sys.fct.type         : specifies the type of the functions a_j, b_j,
%                           c_j, e_j in (1) as follows:
%                           = 'd': a_j, b_j, c_j, e_j are all delay 
%                                  functions, i.e., of the form 
%                                     x_j(s) = exp(-i*s*tau_x_j),       (2)
%                                  where x is either a, b, c, or e;
%                           = 'l': a_j, b_j, c_j, e_j are general 
%                                  functions.
%    sys.fct.a, sys.fct.b,
%    sys.fct.c, sys.fct.e : vectors specifying the functions a_j, b_j, c_j,
%                           and e_j, respectively.
%                           If sys.fct.type == 'd', then they are
%                           represented as sys.fct.x(j) = tau_x_j, where
%                           tau_x_j is as in (2) and x is either a, b, c,
%                           or e.
%                           If sys.fct.type == 'l', then sys.fct.x is a
%                           function handle that returns the evaluation of
%                           all x_j(s) as in (1), where x is either a, b,
%                           c, or e. That is 
%                           sys.fct.x(s) = [ x_1(s); ...; x_{k_x}(s) ].
%    sys.des              : specifies whether the function G(s) can be
%                           realized by a linear (descriptor) system as
%                           follows:
%                           = 0: the function is not realized by a linear
%                                system;
%                           = 1: the function is realized by a linear
%                                system.
%    sys.d2               : specifies whether each of the functions A(s),
%                           B(s), C(s), and E(s) is represented by only one
%                           matrix as follows:
%                           = 0: one of the functions is represented by
%                                more than one matrix;
%                           = 1: all functions are represented by only one
%                                matrix.
% U, V    : Matrices of the same dimensions representing projection spaces.
% opt     : Struct containing options (use default values if empty). The
%           following options can be specified:
%    opt.makeRegular      : specifies whether a regularization should be
%                           performed as follows:
%                           = 0: no regularization is performed;
%                           = 1: regularization is performed;
%                           (default = 0).
%                           This option is only used, if G(s) is rational
%                           and can be realized by a linear (descriptor)
%                           system. Then s*E(s) - A(s) = s*E - A for some
%                           fixed matrices E and A. This algorithm removes
%                           small singular values of [ E, A ] and [ E; A ]
%                           according to the tolerance specified in
%                           opt.trunctol.
%    opt.trunctol         : truncation tolerance for removing small
%                           singular values in the regularization
%                           procedure (default = 1e-12). This option is
%                           only used if opt.makeRegular == 1.
%
% Outputs:
%            
% redsys  : Struct containing the structural information of the the reduced
%           function Gr(s). Only the matrices are updated, all the other
%           information is copied from sys.
%    redsys.A, redsys.B,
%    redsys.C, redsys.D   : arrays containing the reduced matrices 
%                           U'*A_j*V, U'*B_j, C_j*V, and U'*E_j*V in (2),
%                           respectively. If k in k_x is larger than 1, 
%                           i.e., the number of summands within the 
%                           function Xr(s) is larger than 1, then they are
%                           stored in a cell array as 
%                           redsys.X = { X_1, ..., X_{k_x} }, where x is
%                           either a, b, c, e and X is either A, B, C, or 
%                           E.
%    redsys.fct.type      : = sys.fct.type.
%    redsys.fct.a         : = sys.fct.a.
%    redsys.fct.b         : = sys.fct.b.
%    redsys.fct.c         : = sys.fct.c.
%    redsys.fct.e         : = sys.fct.e.
%    redsys.des           : = sys.des.
%    redsys.d2            : = sys.d2.
%
% REFERENCES:
%
% [1] N. Aliyev, P. Benner, E. Mengi, P. Schwerdtner, and M. Voigt.
%     Large-scale computation of L-infinity norms by a greedy subspace
%     method, SIAM J. Matrix Anal. Appl., 38(4):1496-1516, 2017.
%
% AUTHORS:
%
% Nicat Aliyev and Emre Mengi, Koc University, Istanbul, Turkey.
%
% Paul Schwerdtner and Matthias Voigt, Technische Universitaet Berlin,
% Institut fuer Mathematik, Berlin, Germany.
%
% 21/04/2017.
%
% REVISIONS:
%
% Paul Schwerdtner, Matthias Voigt, 09/2017, 12/2017.
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
%% SET OPTIONS.
%
% opt.makeRegular: specifies whether a regularization should be performed. 
%
if isfield( opt, 'makeRegular' ) 
    makeRegular = opt.makeRegular;
else % default
    makeRegular = 0;
end
%
% opt.trunctol: truncation tolerance for removing small singular values in 
%               the regularization procedure.
%
if isfield( opt, 'trunctol' ) 
    trunctol = opt.trunctol;
else % default
    trunctol = 1e-12;
end
%
%% COMPUTATIONS. 
%
% The reduced function gets the same properties that were specified in the 
% original function and which will not be changed by the reduction of 
% dimension.
%
redsys.D = sys.D;
redsys.fct = sys.fct;
redsys.des = sys.des;
redsys.d2 = sys.d2;
%
% Check whether all defining functions are constructed by only one matrix.
%
if sys.d2
    redsys.A = U'*sys.A*V;
    redsys.E = U'*sys.E*V;
    redsys.B = U'*sys.B;
    redsys.C = sys.C*V;
else
    if iscell( sys.A ) == 1 % if ( k_a > 1 ) all summands are reduced.
        for i = 1:length( sys.A )
            redsys.A{ i } = U'*sys.A{ i }*V;
        end
    else
        redsys.A = U'*sys.A*V;  
    end
    %
    if iscell( sys.B ) == 1 % if ( k_b > 1 ) all summands are reduced.
        for i = 1:length( sys.B )
            redsys.B{ i } = U'*sys.B{ i };
        end
    else
        redsys.B = U'*sys.B;  
    end
    %
    if iscell( sys.C ) == 1 % if ( k_c > 1 ) all summands are reduced.
        for i = 1:length( sys.C )
            redsys.C{ i } = sys.C{ i }*V;
        end
    else
        redsys.C = sys.C*V;  
    end
    %
    if iscell( sys.E ) == 1 % if ( k_e > 1 ) all summands are reduced.
        for i = 1:length( sys.E )
            redsys.E{ i } = U'*sys.E{ i }*V;
        end
    else
        redsys.E = U'*sys.E*V;  
    end
end
%
if makeRegular
    if redsys.des == 1
        %
        % Check whether the pencil sE-A is regular and regularize it if
        % necessary.
        %
        nr = size( redsys.E, 1 );
        [ Y, Sl, ~ ] = svd( [ redsys.E, redsys.A ] );
        [ ~, Sr, Z ] = svd( [ redsys.E; redsys.A ] );
        Sl = diag( Sl )./Sl( 1 );
        Sr = diag( Sr )./Sr( 1 );
        Sl = [ Sl; 0 ];
        Sr = [ Sr; 0 ];
        dim = 0;
        j = 1;
        %
        while Sl( j ) > trunctol && Sr( j ) > trunctol && j <= nr
            dim = dim + 1;
            j = j + 1;
        end
        %
        % Truncate the small singular values.
        %
        U = U*Y( :, 1:dim );
        V = V*Z( :, 1:dim );
        redsys.A = U'*sys.A*V;
        redsys.E = U'*sys.E*V;
        redsys.B = U'*sys.B;
        redsys.C = sys.C*V;
    end
end
%
return