function [ value, der ] = evalFun( sys, s0 )
%
% PURPOSE:
%
% To evaluate the function 
%
%   G(s) = C(s)*(s*E(s) - A(s))^(-1)*B(s) + D(s)                        (1)
%
% and its derivative at a complex point s0. Here we assume that G(s) is 
% holomorphic in an open domain around s0 and that that A(s), B(s), C(s),
% D(s), and E(s) are functions that are meromorphic in an area enclosing
% the point s0. These functions are of the form
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
% ARGUMENTS:
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
%                                     x_j(s) = exp(-s*tau_x_j);         (3)
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
%                           all x_j(s) as in (3), where x is either a, b,
%                           c, d, or e. That is 
%                           sys.fct.x(s) = [ x_1(s); ...; x_{k_x}(s) ].
%    sys.des              : specifies whether the function G(s) can is
%                           realized by a linear (descriptor) system as
%                           follows:
%                           = 0: the function is not realized by a linear
%                                system;
%                           = 1: the function is realized by a linear
%                                system.        
% s0      : Complex point at which the function G(s) and its derivatives 
%           are evaluated.
% 
% Outputs:
%            
% value   : value of the functions G(s) in (1) at s0.
% der     : first derivative of the function G(s) in (1) at s0.
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
% 27/09/2017.
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
%% SET OPTIONS AND CHECK ARGUMENTS.
%
% Check the dimensions of the inputs.
%
if iscell( sys.A )
    [ na1, na2 ] = size( sys.A{ 1 } );
    dA = length( sys.A );
else
    [ na1, na2, dA ] = size( sys.A );
end
%
if iscell( sys.B )
    [ n1, m1 ] = size( sys.B{ 1 } );
    dB = length( sys.B );
else
    [ n1, m1, dB ] = size( sys.B );
end
%
if iscell( sys.C )
    [ p1, n2 ] = size( sys.C{ 1 } );
    dC = length( sys.C );
else
    [ p1, n2, dC ] = size( sys.C );
end
%
if iscell( sys.D )
    [ p2, m2 ] = size( sys.D{ 1 } );
    dD = length( sys.D );
else
    [ p2, m2, dD ] = size( sys.D );
end
%
if iscell( sys.E )
    [ ne1, ne2 ] = size( sys.E{ 1 } );
    dE = length( sys.E );
else
    [ ne1, ne2, dE ] = size( sys.E );
end
%
if n1 ~= n2 || p1 ~= p2 || m1 ~= m2
    error( 'The dimensions of B, C, and D are not compatible.' );
end
%
if ne1 ~= ne2
    error( 'The matrix E must be square.' );
elseif na1 ~= na2
    error( 'The matrix A must be square.' );
elseif ne1 ~= na1
    error( 'The dimensions of E and A are not compatible.' );
elseif ne1 ~= n1
    error( [ 'The dimensions of E and A are not compatible with B, ', ... 
        'C, and D.' ] );
end
%
% Check whether sys.A, sys.B, sys.C, sys.D, and sys.E are all consisting of
% only one matrix.
%
sys.d2 = 0;
%
if dA == 1 && dB == 1 && dC == 1 && dD == 1 && dE == 1
    sys.d2 = 1;
end
%
% Check the type of the function G(s) in (1).
%
sys.des = 0;
%
if not( isfield( sys, 'fct' ) ) && sys.d2 == 1
    %
    % The function G(s) corresponds to a linear system. This is indicated
    % by setting sys.des = 1.
    %
    sys.fct.type = 'd';
    sys.des = 1;
    %
end
%
switch sys.fct.type
    case 'd'
        %
        % If no delay is defined then the delay parameter is set to zero.
        %
        if not( iscell( sys.A ) ) && not( isfield( sys.fct, 'a' ) )
            sys.fct.a = 0;
        end
        %
        if not( iscell( sys.B ) ) && not( isfield( sys.fct, 'b' ) )
            sys.fct.b = 0;
        end
        %
        if not( iscell( sys.C ) ) && not( isfield( sys.fct, 'c' ) )
            sys.fct.c = 0;
        end
        %
        if not( iscell( sys.D ) ) && not( isfield( sys.fct, 'd' ) )
            sys.fct.d = 0;
        end
        %
        if not( iscell( sys.E ) ) && not( isfield( sys.fct, 'e' ) )
            sys.fct.e = 0;
        end
        %
        if dB ~= length( sys.fct.b ) || dC ~= length( sys.fct.c ) || ...
                dA ~= length( sys.fct.a ) || dE ~= length( sys.fct.e ) ...
                || dD ~= length( sys.fct.d )
             error( [ 'The number of the delay entries in sys.fct ', ...
                 'does not match the dimensions of sys.A, sys.B, ', ...
                 'sys.C, sys.D, or sys.E.' ] );
        end
    case 'l'
        %
        % If no function sys.fct.x is defined then sys.fct.x = 0.
        %
        if dA == 1 && not( isfield( sys.fct, 'a' ) )
            sys.fct.a = @( x ) 1;
        end
        %
        if dB == 1 && not( isfield( sys.fct, 'b' ) )
            sys.fct.b = @( x ) 1;
        end
        %
        if dC == 1 && not( isfield( sys.fct, 'c' ) )
            sys.fct.c = @( x ) 1;
        end
        %
        if dD == 1 && not( isfield( sys.fct, 'd' ) )
            sys.fct.d = @( x ) 1;
        end
        %
        if dE == 1 && not( isfield( sys.fct, 'e' ) )
            sys.fct.e = @( x ) 1;
        end
        %
        if dB ~= length( sys.fct.b( 0 ) ) || ...
                dC ~= length( sys.fct.c( 0 ) ) || ... 
                dA ~= length( sys.fct.a( 0 ) ) || ...
                dE ~= length( sys.fct.e( 0 ) ) || ...
                dD ~= length( sys.fct.d( 0 ) )
            error( [ 'The dimensions of the functions in sys.fct do ', ...
                'not match the dimensions of sys.A, sys.B, sys.C, ', ...
                'sys.D, or sys.E.' ] );
        end
end
%
%% COMPUTATIONS.
%
% Check whether the derivative is needed.
%
if nargout == 1
    %
    % The derivative is not needed.
    %
    h = 0;
    nloop = 1;
elseif nargout == 2
    %
    % The derivative is needed. Choose step size for the computation of the
    % derivatives.
    %
    h = 1e-6;
    s0 = s0 - h;
    nloop = 3;
end
%
f = cell( 1, nloop );
%
% Evaluate the factors in G(s). 
%
for j = 1:nloop
    if sys.des
        %
        % The function G(s) is realized by a linear (descriptor) system.
        %
        As = sys.A;
        Bs = sys.B;
        Cs = sys.C;
        Ds = sys.D;
        Es = sys.E;
    else
        %
        % The function G(s) is not realized by a linear (descriptor)
        % system.
        %
        if iscell( sys.A )
            [ n, ~ ] = size( sys.A{ 1 } );
            As = sparse( n, n );
        else
            [ n, ~ ] = size( sys.A );
            As = sparse( n, n );  
        end
        %
        if iscell( sys.B )
            [ ~, m ] = size( sys.B{ 1 } );
            Bs = sparse( n, m );
        else
            [ ~, m ] = size( sys.B );
            Bs = sparse( n, m );  
        end
        %
        if iscell( sys.C )  
            [ p, ~ ] = size( sys.C{ 1 } );
            Cs = sparse( p, n );
        else
            [ p, ~ ] = size( sys.C );
            Cs = sparse( p, n );  
        end
        %
        Ds = sparse( p, m );  
        Es = sparse( n, n );  
        %
        switch sys.fct.type
            case 'd'
                %
                % All functions stored in sys.fct are delays.
                %
                if iscell( sys.A )    
                    for ind = 1:length( sys.A )
                        As = As + exp( -s0.*sys.fct.a( ind ) )* ...
                            sys.A{ ind };
                    end
                else
                    As = sys.A*exp(- s0.*sys.fct.a );  
                end
                %
                if iscell( sys.B )    
                    for ind = 1:length( sys.B )
                        Bs = Bs + exp(- s0.*sys.fct.b( ind ) )* ...
                            sys.B{ ind };
                    end
                else
                    Bs = sys.B*exp(- s0.*sys.fct.b );  
                end
                %
                if iscell( sys.C )    
                    for ind = 1:length( sys.C )
                        Cs = Cs + exp( -s0.*sys.fct.c( ind ) )* ...
                            sys.C{ ind };
                    end
                else
                    Cs = sys.C*exp( -s0.*sys.fct.c );  
                end
                %
                if iscell( sys.D )    
                    for ind = 1:length( sys.D )
                        Ds = Ds + exp( -s0.*sys.fct.d( ind ) )* ...
                            sys.D{ ind };
                    end
                else
                    Ds = sys.D*exp( -s0.*sys.fct.d );  
                end
                %
                if iscell( sys.E )    
                    for ind = 1:length( sys.E )
                        Es = Es + exp( -s0.*sys.fct.e( ind ) )* ...
                            sys.E{ ind };
                    end
                else
                    Es = sys.E*exp( -s0.*sys.fct.e );  
                end
                %
            case 'l'
                %
                % The functions stored in sys.fct can be arbitrary.
                %
                if iscell( sys.A )
                    feval = sys.fct.a( s0 );
                    for ind = 1:length( sys.A )
                        As = As + feval( ind )*sys.A{ ind };
                    end
                else
                    As = sys.A*sys.fct.a( s0 );  
                end
                %
                if iscell( sys.B )
                    feval = sys.fct.b( s0 );
                    for ind = 1:length( sys.B )
                        Bs = Bs + feval( ind )*sys.B{ ind };
                    end
                else
                    Bs = sys.B*sys.fct.b( s0 );  
                end
                %
                if iscell( sys.C )
                    feval = sys.fct.c( s0 );
                    for ind = 1:length( sys.C )
                        Cs = Cs + feval( ind )*sys.C{ ind };
                    end
                else
                    Cs = sys.C*sys.fct.c( s0 );  
                end
                %
                if iscell( sys.D )
                    feval = sys.fct.d( s0 );
                    for ind = 1:length( sys.D )
                        Ds = Ds + feval( ind )*sys.D{ ind };
                    end
                else
                    Ds = sys.D*sys.fct.d( s0 );  
                end
                %
                if iscell( sys.E )
                    feval = sys.fct.e( s0 );
                    for ind = 1:length( sys.E )
                        Es = Es + feval( ind )*sys.E{ ind };
                    end
                else
                    Es = sys.E*sys.fct.e( s0 );  
                end
        end
    end
    %
    % The function G(s) is evaluated at s0.
    %
    H = Cs*( ( s0*Es - As )\Bs ) + Ds;
    f{ j } = H;
    s0 = s0 + h;
end
%
% Assign the outputs.
%
if nargout == 1
    %
    % The function value is assigned.
    %
    value = f{ 1 };
elseif nargout == 2
    %
    % The function value is assigned.
    %
    value = f{ 2 };
    %
    % The complex derivative values are computed using a central difference
    % scheme.
    %
    der = ( f{ 3 } - f{ 1 } )/( 2*h );
end
%
return
