function [ value, der, der2 ] = Linf( w, sys )
%
% PURPOSE:
%
% To evaluate the function 
%
%   s(w) = sigma_max(G(i*w)),
%
% and its derivative and second derivative at a given point w, where
%
%   G(s) = C(s)*(s*E(s) - A(s))^(-1)*B(s) + D(s),
%
% and sigma_max denotes the maximum singular value.
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
%                                     x_j(s) = exp(-i*s*tau_j),         (3)
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
% w       : Frequency at which the function s(w) and its derivatives should
%           be evaluated.
% 
% Outputs:
%            
% value   : value of the functions s(w).
% der     : first derivative of s(w).
% der2    : second derivative of s(w).
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
% Matthias Voigt, 09/2017, 06/2018.
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
%% Check input arguments.
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
    error( 'The dimensions of E and A are not compatible with B, C, ', ...
        'and D.' );
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
if not( isfield( sys, 'fct' ) ) && sys.d2 == 1
    %
    % The function G(s) corresponds to a linear system. This is indicated
    % by setting sys.des = 1. In this case, the Boyd-Balakrishnan algorithm
    % to compute the L-infinity norm is called. 
    %
    sys.fct.type = 'd';
    sys.des = 1;
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
             error( 'The number of the delay entries in sys.fct does ', ...
                 'not match the dimensions of sys.A, sys.B, sys.C, ', ... 
                 'sys.D, or sys.E.' );
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
            error( 'The dimensions of the functions in sys.fct do ', ...
                'not match the dimensions of sys.A, sys.B, sys.C, ', ...
                'sys.D, or sys.E.' );
        end
end
%
%% COMPUTATIONS.
%
% Choose step size for the computation of the derivatives.
%
h = 1e-6;
w = w - h;
f = zeros( 1, 3 );
%
% Evaluate the factors in G(i*w). 
%
for j = 1:3
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
                    for i = 1:length( sys.A )
                        As = As + exp( -1i.*w.*sys.fct.a( i ) )*sys.A{ i };
                    end
                else
                    As = sys.A*exp( -1i.*w.*sys.fct.a );  
                end
                %
                if iscell( sys.B )
                    for i = 1:length( sys.B )
                        Bs = Bs + exp( -1i.*w.*sys.fct.b( i ) )*sys.B{ i };
                    end
                else
                    Bs = sys.B*exp( -1i .* w .* sys.fct.b );  
                end
                %
                if iscell( sys.C )
                    for i = 1:length( sys.C )
                        Cs = Cs + exp( -1i.*w.*sys.fct.c( i ) )*sys.C{ i };
                    end
                else
                    Cs = sys.C*exp( -1i.*w.*sys.fct.c );  
                end
                %
                if iscell( sys.D )
                    for i = 1:length( sys.D )
                        Ds = Ds + exp( -1i.*w.*sys.fct.d( i ) )*sys.D{ i };
                    end
                else
                    Ds = sys.D*exp( -1i.*w.*sys.fct.d );  
                end
                %
                if iscell( sys.E )
                    for i = 1:length( sys.E )
                        Es = Es + exp( -1i.*w.*sys.fct.e( i ) )*sys.E{ i };
                    end
                else
                    Es = sys.E*exp( -1i.*w.*sys.fct.e );  
                end
                %
            case 'l'
                %
                % The functions stored in sys.fct can be arbitrary.
                %
                if iscell( sys.A )
                    feval = sys.fct.a( 1i*w );
                    for i = 1:length( sys.A )
                        As = As + feval( i )*sys.A{ i };
                    end
                else
                    As = sys.A*sys.fct.a( 1i*w );  
                end
                %
                if iscell( sys.B )
                    feval = sys.fct.b( 1i*w );
                    for i = 1:length( sys.B )
                        Bs = Bs + feval( i )*sys.B{ i };
                    end
                else
                    Bs = sys.B*sys.fct.b( 1i*w );  
                end
                %
                if iscell( sys.C )
                    feval = sys.fct.c( 1i*w );
                    for i = 1:length( sys.C )
                        Cs = Cs + feval( i )*sys.C{ i };
                    end
                else
                    Cs = sys.C*sys.fct.c( 1i*w );  
                end
                %
                if iscell( sys.D )
                    feval = sys.fct.d( 1i*w );
                    for i = 1:length( sys.D )
                        Ds = Ds + feval( i )*sys.D{ i };
                    end
                else
                    Ds = sys.D*sys.fct.d( 1i*w );  
                end
                %
                if iscell( sys.E )
                    feval = sys.fct.e( 1i*w );
                    for i = 1:length( sys.E )
                        Es = Es + feval( i )*sys.E{ i };
                    end
                else
                    Es = sys.E*sys.fct.e( 1i*w );  
                end
        end
    end
    %
    % The function G(i*w) is evaluated.
    %
    H = Cs*( ( 1i*w*Es - As )\Bs ) + Ds;
    H = full( H );
    smax = svd( H );
    f( j ) = smax( 1, 1 );
    w = w + h;
end
%
% The function value is assigned.
%
value = f( 2 );
%
% The derivative values are computed using a central difference scheme.
%
der = ( f( 3 ) - f( 1 ) )/( 2*h );
der2 = ( f( 3 ) + f( 1 ) - 2*f( 2 ) )/( h^2 );
%
return