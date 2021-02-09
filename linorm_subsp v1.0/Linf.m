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
%                                system,
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
% -
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
% Choose step size for the computation of the derivatives.
%
h = 1e-6;
w = w - h;
fval = zeros( 1, 3 );
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
        Ds = sparse( p, m );  
        Es = sparse( n, n );  
    end
    %
    switch sys.fct.type
        case 'd'
            %
            % All functions stored in sys.fct are delays.
            %
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
                Bs = sys.B*exp( -1i .* w .* sys.fct.b );  
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
            if ( iscell( sys.D ) )    
                for i = 1:length( sys.D )
                    Ds = Ds + exp( -1i.*w.*sys.fct.d( i ) )*sys.D{ i };
                end
            else
                Ds = sys.D*exp( -1i.*w.*sys.fct.d );  
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
            % The functions stored in sys.fct can be arbitrary.
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
            if ( iscell( sys.D ) )
                f = sys.fct.d( 1i*w );
                for i = 1:length( sys.D )
                    Ds = Ds + f( i )*sys.D{ i };
                end
            else
                Ds = sys.D*sys.fct.d( 1i*w );  
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
    %
    % The function G(i*w) is evaluated.
    %
    H = Cs*( ( 1i*w*Es - As )\Bs ) + Ds;
    H = full( H );
    smax = svd( H );
    fval( j ) = smax( 1, 1 );
    w = w + h;
end
%
% The function value is assigned.
%
value = fval( 2 );
%
% The derivative values are computed using a central difference scheme.
%
der = ( fval( 3 ) - fval( 1 ) )/( 2*h );
der2 = ( fval( 3 ) + fval( 1 ) - 2*fval( 2 ) )/( h^2 );
%
return
