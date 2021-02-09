function [ V, U ] = getSubspace( sys, w, maxSing )
%
% PURPOSE:
%
% To obtain matrices V and U such that the maximum singular value of the 
% function 
%
%   G(s) = C(s)*(s*E(s) - A(s))^(-1)*B(s)
%
% and its derivative are equal to the maximum singular value and its 
% derivative of the reduced function
%
%   Gr(s) = C(s)*V*(s*U'*E(s)*V - U'*A(s)*V)^(-1)*U'*B(s)
%
% coincide at the point s = i*w. The method is based on Hermite
% interpolation of G(s) at the point s = i*w. The matrix-valued functions
% are of the form
%
%   A(s) = a_1(s)*A_1 + ... + a_{k_a}(s)*A_{k_a},
%   B(s) = b_1(s)*B_1 + ... + b_{k_b}(s)*B_{k_b},                       (1)
%   C(s) = c_1(s)*C_1 + ... + c_{k_c}(s)*C_{k_c},
%   E(s) = e_1(s)*E_1 + ... + e_{k_e}(s)*E_{k_e},
%
% where x_1, ..., x_{k_x} are a scalar-valued functions assumed to be
% meromorphic in an area ecnlosing the imaginary axis and X_1, ..., X_{k_x}
% is a constand matrix. Here, x is either a, b, c, or e and X is either A,
% B, C, or E.
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
%                                     x_j(s) = exp(-i*s*tau_j),         (2)
%                                  where x is either a, b, c, or e;
%                           = 'l': a_j, b_j, c_j, e_j are general 
%                                  functions.
%    sys.fct.a, sys.fct.b,
%    sys.fct.c, sys.fct.e : vectors specifying the functions a_j, b_j, c_j,
%                           and e_j, respectively.
%                           If sys.fct.type == 'd', then they are
%                           represented as sys.fct.x(j) = tau_j, where
%                           tau_j is as in (2) and x is either a, b, c, or 
%                           e.
%                           If sys.fct.type == 'l', then sys.fct.x is a
%                           function handle that returns the evaluation of
%                           all x_j(s) as in (1), where x is either a, b,
%                           c, or e. That is 
%                           sys.fct.x(s) = [ x_1(s); ...; x_{k_x}(s) ].
%    sys.des              : specifies whether the function G(s) can is
%                           realized by a linear (descriptor) system as
%                           follows:
%                           = 0: the function is not realized by a linear
%                                system,
%                           = 1: the function is realized by a linear
%                                system.
% w       : Frequency at which the function G(s) is to be interpolated.
% maxSing : specifies how the projection spaces are constructed as follows:
%           = 0: all the singular vectors of G(i*w) at an interpolation 
%                point i*w are included in the projection spaces; 
%           = 1: only the singular vectors corresponding to the largest 
%                singular value of G(i*w) at an interpolation point i*w are
%                included in the projection spaces.
%
% Outputs:
%
% V, U    : Subspaces used to construct the reduced function Gr(s).
%
% REFERENCES:
%
% [1] N. Aliyev, P. Benner, E. Mengi and M. Voigt. Large-scale computation
%     of H-infinity norms by a greedy subspace method.
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
% -
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LINORM_SUBSP 1.0 Copyright (C) 2017 Nicat Aliyev, Emre Mengi,
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
% Check the type of function.
%
if sys.des  
    %
    % The function G(s) is realized by a linear system.
    %
    [ ~, m ] = size( sys.B );
    [ p, ~ ] = size( sys.C );
    As = sys.A;
    Bs = sys.B;
    Cs = sys.C;
    Es = sys.E;
else
    % 
    % The function G(s) is not realized by a linear system, the 
    % matrix-valued functions given in (1) have to be evaluated.
    %
    if ( iscell( sys.A ) )
        [ n, ~ ] = size( sys.A{ 1 } );
        As = sparse( n, n );
    else
        [ n, ~ ] = size( sys.A );
        As = sparse( n, n );  
    end
    %
    if ( iscell( sys.B ) )
        [ ~, m ] = size( sys.B{ 1 } );
        Bs = sparse( n, m );
    else
        [ ~, m ] = size( sys.B );
        Bs = sparse( n, m );  
    end
    %
    if ( iscell( sys.C ) )  
        [ p, ~ ] = size( sys.C{ 1 } );
        Cs = sparse( p, n );
    else
        [ p, ~ ] = size( sys.C );
        Cs = sparse( p, n );  
    end
    %
    Es = sparse( n, n );  
    %
    switch sys.fct.type
        %
        % For all matrices the delay function e^{-1i*w*tau_k} is evaluated
        % at the given value. 
        %
        case 'd'
        if ( iscell( sys.A ) )    
            for i = 1:length( sys.A )
                As = As + exp( -1i.*w.*sys.fct.a( i ) )*sys.A{ i };
            end
        else
            As = sys.A*exp( -1i.*w.*sys.fct.a );  
        end
        %
        if ( iscell( sys.B ) )    
            for i = 1:length( sys.B )
                Bs = Bs + exp( -1i.*w.*sys.fct.b( i ) )*sys.B{ i };
            end
        else
            Bs = sys.B*exp( -1i.*w.*sys.fct.b );  
        end
        %
        if ( iscell( sys.C ) )    
            for i = 1:length( sys.C )
                Cs = Cs + exp( -1i.*w.*sys.fct.c( i ) )*sys.C{ i };
            end
        else
            Cs = sys.C*exp( -1i.*w.*sys.fct.c );  
        end
        %
        if ( iscell( sys.E ) )    
            for i = 1:length( sys.E )
                Es = Es + exp( -1i.*w.*sys.fct.e( i ) )*sys.E{ i };
            end
        else
            Es = sys.E*exp( -1i.*w.*sys.fct.e );  
        end
        %
        case 'l'
        %
        % For all matrices the corresponding function is evaluated at the
        % given value.
        %
        if ( iscell( sys.A ) )
            f = sys.fct.a( 1i*w );
            for i = 1:length( sys.A )
                As = As + f( i )*sys.A{ i };
            end
        else
            As = sys.A*sys.fct.a( 1i*w );  
        end
        %
        if ( iscell( sys.B ) )
            f = sys.fct.b( 1i*w );
            for i = 1:length( sys.B )
                Bs = Bs + f( i )*sys.B{ i };
            end
        else
            Bs = sys.B*sys.fct.b( 1i*w );  
        end
        %
        if ( iscell( sys.C ) )
            f = sys.fct.c( 1i*w );
            for i = 1:length( sys.C )
                Cs = Cs + f( i )*sys.C{ i };
            end
        else
            Cs = sys.C*sys.fct.c( 1i*w );  
        end
        %
        if ( iscell( sys.E ) )
            f = sys.fct.e( 1i*w );
            for i = 1:length( sys.E )
                Es = Es + f( i )*sys.E{ i };
            end
        else
            Es = sys.E*sys.fct.e( 1i*w );  
        end
    end 
end
%
if maxSing
    %
    % Compute the maximum singular value of G(i*w) and construct a subspace
    % from the corresponding singular vectors.
    %
    [ U1, ~, V1 ] = svd( full( Cs*( ( 1i.*w.*Es-As )\Bs ) ) );
    V = ( 1i.*w.*Es-As )\( Bs*V1( :, 1 ) );
    U = ( 1i.*w.*Es-As )'\( Cs'*U1( :, 1 ) );
else
    % 
    % Solve the linear systems to construct an initial subspace with the
    % desired interpolation property.
    %
    X = ( 1i.*w.*Es-As )\Bs;
    Y = ( 1i.*w.*Es-As )'\Cs';
    %
    % Compare the input dimension m and output dimension p to make the 
    % subspaces U and V have the same dimension.
    %
    if ( m == p )
        V = X;
        U = Y;
    elseif ( m < p )                  
        V = X;
        U = Y*Cs*X;
    else
        V = X*Bs'*Y;
        U = Y;
    end
end
%
if sys.des
    %
    % linorm_h only takes real valued matrices, so transform U and V to
    % real matrices.
    %
    V = [ real( V ), imag( V ) ];
    U = [ real( U ), imag( U ) ];
end
%
return