function [ f, z, info ] = linorm_subsp( sys, opt )
%
%% PURPOSE:
%
% To compute the L-infinity norm of a function G(s) via a function handle
% or a function of the form
%
%   G(s) = C(s)*(s*E(s) - A(s))^(-1)*B(s) + D(s),                       (1)
%
% where A(s), B(s), C(s), D(s), and E(s) are functions that are meromorphic
% in an area enclosing the imaginary axis. These functions are of the form
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
% assumed to be sparse. It is also assumed that G(s) is an L-infinity-
% function meaning that it is bounded on the imaginary axis.
%
% The L-infinity norm of G is defined as  
%
%   ||G||_Linf := sup_{w real} ||G(iw)||_2 
%               = sup_{w real} sigma_max(G(iw)),
%
% where sigma_max denotes the maximum singular value of a matrix.
%
% ARGUMENTS:
%
% Inputs:
%
% sys     : A function handle that evaluates the function G(s), if
%           opt.fnHandle == 0; or otherwise a struct containg the structual 
%           information of the function G(s) specified as follows:
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
% opt     : Struct containing options (use default values if empty). The
%           following options can be specified:
%    opt.tolz             : relative tolerance on the change of the optimal
%                           frequencies between two consecutive iterations.
%                           If the computed optimal frequencies between two
%                           consecutive iterations have relative distance
%                           less than opt.tolz, then the algorithm is
%                           assumed to be converged (default = 1e-6).
%    opt.tolf             : relative tolerance on the change of the
%                           computed L-infinity norms between two 
%                           consecutive iterations. If the computed 
%                           L-infinity norms between two consecutive
%                           iterations have relative distance less than
%                           opt.tolf, then the algorithm is assumed to be
%                           converged (default = 1e-6).
%    opt.maxit 		      : maximum number of iterations allowed until 
%                           termination of the algorithm (default = 30).
%    opt.initialPoints	  : initial frequencies of the initial 
%                           interpolation points on the imaginary axis.
%    opt.prtlevel         : specifies the print level as follows:
%                           = 0: return no information;
%                           = 1: return minimal information;
%                           = 2: return full information;
%                           (default = 0).
%    opt.boydbalak        : specifies which implementation of the
%                           Boyd-Balakrishnan algorithm is used for the
%                           computation of the L-infinity-norm of the
%                           reduced functions as follows:
%                           = 0: use the FORTRAN routine AB13HD.f which is
%                                called by the gateway function generated
%                                by linorm_h.F;
%                           = 1: use the function 'norm' of the MATLAB
%                                Control System Toolbox (WARNING: This may
%                                not work, if the reduced descriptor matrix
%                                pars.E is singular.);
%                           (default = 0).
%    opt.eigopt.bounds    : a vector [ lowerBound,  upperBound ] specifying
%                           the interval in which eigopt shall optimize the
%                           maximum singular value of the function G(s) in 
%                           (1). This option is only used if 
%                           opt.doLoewner == 0.
%    opt.eigopt.gamma     : a lower bound for the second derivative of the
%                           the function -sigma_max(G(i*w)) where w is 
%                           within opt.eigopt.bounds. This option is only
%                           used if opt.doLoewner == 0.
%    opt.keepSubspaces    : specifies whether the intermediate subspaces
%                           obtained during the iteration are kept as
%                           follows:
%                           = 0: only the subspaces of the initial reduced
%                                function and the past two subspaces
%                                obtained during the iteration are kept;
%                           = 1: all the interpolation subspaces are kept 
%                                during the iteration. This will let the 
%                                dimensions of projection spaces grow in
%                                every iteration but is more robust in many
%                                cases;
%                           (default = 1).
%    opt.biorth           : specifies which orthonormalization scheme is
%                           used as follows:
%                           = 0: the intermediate projection spaces
%                                contained in U and V are orthonormalized
%                                separately, i.e., U'*U = V'*V = eye(k) for
%                                some k;
%                           = 1: the intermediate projection spaces are
%                                bi-orthonormalized, i.e., U'*V = eye(k)
%                                for some k;
%                           (default = 0).
%    opt.orthtol          : relative truncation tolerance in svd for the 
%                           determination of an orthonormal basis of a
%                           subspace or bi-orthonormal bases of two
%                           subspaces (default = 1e-12).
%    opt.maxSing          : specifies how the projection spaces are
%                           constructed as follows:
%                           = 0: all the singular vectors of G(i*w) at an 
%                                interpolation point i*w are included in
%                                the updated projection spaces; 
%                           = 1: only the singular vectors corresponding to
%                                the largest singular value of G(i*w) at an 
%                                interpolation point i*w are included in 
%                                the updated projection spaces;
%                           (default = 0).
%    opt.doLoewner        : specifies whether the Loewner approach is
%                           used to calculate the L-infinity norm as
%                           follows:
%                           = 0: the Loewner framework is not used; 
%                           = 1: the Loewner framework is used;
%                           (default = 0).
%    opt.fnHandle         : specifies whether a function handle is passed
%                           to linorm_subsp. In this case, the Loewner
%                           approach is used to construct a rational
%                           function that matches the function evaluations.
%                            = 0: no function handle, but a function of the
%                                 previously defined format is passed;
%                            = 1: a function handle is passed
%                           (default = 0).
%
% Outputs:
%   
% f       : L-infinity norm of the function G(s) in (1).
% w       : Frequency at which the norm value is attained.
% info    : Struct which contains the following information:
%    info.time             : time needed to compute the result.
%    info.iterations       : number of iterations at termination of the
%                            algorithm.
%    info.finaltolz        : relative distance of the variable z between
%                            the last two iterations before termination.
%    info.finaltolf        : relative distance of the variable f between
%                            the last two iterations before termination.
%    info.termcrit         : specifies which termination criteria are
%                            satisfied as follows:
%                            = 0: both the optimal frequency z and the 
%                                 L-infinity norm f have converged;
%                            = 1: only the optimal frequency z has
%                                 converged, but the L-infinity norm may be
%                                 inaccurate;
%                            = 2: only the L-infinity norm f has converged,
%                                 but the optimal frequency may be
%                                 inaccurate;
%                            = 3: the algorithm has terminated by reaching
%                                 the maximum number of iterations.
%    info.error            : contains an error indicator as follows:
%                            = 0: return without an error;
%                            = 1: the maximum number of iterations 
%                                 specified in opt.maxit has been exceeded.
%
% REFERENCES:
%
% [1] N. Aliyev, P. Benner, E. Mengi, P. Schwerdtner, and M. Voigt.
%     Large-scale computation of L-infinity norms by a greedy subspace
%     method, SIAM J. Matrix Anal. Appl., 38(4):1496-1516, 2017.
%
% [2] P. Schwerdtner and M. Voigt. Computation of the L-infinity-norm
%     using rational interpolation, FAC-PapersOnLine, 2018. Joint 9th IFAC 
%     Symposium on Robust Control Design and 2nd IFAC Workshop on Linear
%     Parameter Varying Systems, Florianópolis, Brazil, 2018. Accepted for
%     publication.
%
% AUTHORS:
%
% Nicat Aliyev and Emre Mengi, Koc University, Istanbul, Turkey.
%
% Paul Schwerdtner and Matthias Voigt, Technische Universitaet Berlin,
% Institut fuer Mathematik, Berlin, Germany.
%
% 24/04/2017.
%
% REVISIONS:
%
% Paul Schwerdtner, Matthias Voigt, 09/2017, 12/2017.
% Matthias Voigt, 06/2018
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
%% SET OPTIONS AND CHECK ARGUMENTS.
%
% Start timing.
%
tic;
warning( 'off', 'all' );
%
info.error = 0;
%
% tolz: relative tolerance on the change of the optimal frequencies between
%       two consecutive iterations. If the computed optimal frequencies
%       between two consecutive iterations have relative distance less than
%       opt.tolz, then the algorithm is assumed to be converged.
%
if isfield( opt, 'tolz' )
    tolz = opt.tolz;
else % default
    tolz = 1e-6;
end
%
% tolf: relative tolerance on the change of the computed L-infinity norms
%       between two consecutive iterations. If the computed 
%       L-infinity norms between two consecutive iterations have relative
%       distance less than opt.tolf, then the algorithm is assumed to be
%       converged.
%
if isfield( opt, 'tolf' )
    tolf = opt.tolf;
else % default
    tolf = 1e-6;
end
%
% maxit: maximum number of iterations allowed until termination of the 
%        algorithm.
%
if isfield( opt, 'maxit' )
    maxit = opt.maxit;
else % default
    maxit = 30;
end
%
% boydbalak: specifies which implementation of the Boyd-Balakrishnan 
%            algorithm is used for the computation of the L-infinity-norm
%            of the reduced functions.
%
if isfield( opt, 'boydbalak' )
    boydbalak = opt.boydbalak;
else % default
    boydbalak = 0;
end
%
% doLoewner: specifies whether the Loewner approach is used to calculate
%            the L-infinity-norm. 
%
if isfield( opt, 'doLoewner' )
    doLoewner = opt.doLoewner;
else % default
    doLoewner = 0;
end
%
% opt.fnHandle: specifies whether a function handle is passed to 
%               linorm_subsp. In this case, the Loewner approach is used to
%               construct a rational function that matches the function
%               evaluations.
%
if isfield( opt, 'fnHandle' )
    fnHandle = opt.fnHandle;
else % default
    fnHandle = 0;
end
%
% opt.initialPoints: initial frequencies of the initial interpolation 
%                    points on the imaginary axis.
%
if isfield( opt, 'initialPoints' )
    w = opt.initialPoints;
else
    w = [ ];
    error( 'The argument opt.initialPoints must be specified.' );
end
%
% keepSubspaces: specifies whether the intermediate subspaces obtained
%                during the iteration are kept.
%
if isfield( opt, 'keepSubspaces' )
    if ( opt.keepSubspaces == 0 || opt.keepSubspaces == 1 )
        keepSubspaces = opt.keepSubspaces;
    else
        error( 'The argument opt.keepSubspaces must be 0 or 1.' );
    end
else % default
    keepSubspaces = 1;
end
%
% maxSing: specifies how the projection spaces are constructed.
%
if isfield( opt, 'maxSing' )
    if ( opt.maxSing == 0 || opt.maxSing == 1 )
        maxSing = opt.maxSing;
    else
        error( 'The argument opt.maxSing must be 0 or 1.' );
    end
else % default
    maxSing = 0;
end
%
% biorth: specifies the orthonormalization scheme.
%
if isfield( opt, 'biorth' )
    if ( opt.biorth == 0 || opt.biorth == 1 )
        biorth = opt.biorth;
    else
        error( 'The argument opt.biorth must be 0 or 1.' );
    end
else % default
    biorth = 0;
end
%
% orthtol: relative truncation tolerance in svd for the determination of an
%          orthonormal basis of a subspace or biorthonormal bases of two 
%          subspaces.
%
if isfield( opt, 'orthtol' )
    orthtol = opt.orthtol;
else % default
    orthtol = 1e-12;
end
%
% prtlevel: specifies the print level.
%
if isfield( opt, 'prtlevel' )
    prtlevel = opt.prtlevel;
    if ( prtlevel ~= 0 && prtlevel ~= 1 && prtlevel ~= 2 )
        error( 'The argument opt.prtlevel must be 0, 1, or 2.' );
    end
else % default
    prtlevel = 0;
end
%
% Check the dimensions of the inputs.
%
if fnHandle == 0
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
        error( [ 'The dimensions of E and A are not compatible with ', ...
            'B, C, and D.' ] );
    end
    %
    % Check whether sys.A, sys.B, sys.C, sys.D, and sys.E are all 
    % consisting of only one matrix.
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
    algo = 'e'; % use eigopt to find the L-infinity norm in the main loop.
    %
    if not( isfield( sys, 'fct' ) ) && sys.d2 == 1
        %
        % The function G(s) corresponds to a linear system. This is
        % indicated by setting sys.des = 1. In this case, the 
        % Boyd-Balakrishnan algorithm to compute the L-infinity norm is
        % called. 
        %
        sys.fct.type = 'd';
        algo = 'b';
        sys.des = 1;
        %
        if prtlevel > 0
            %
            % Print info.
            %
            fprintf( [ 'The function corresponds to a linear system. ', ... 
                '\nThe Boyd-Balakrishnan algorithm is called to ', ...
                'compute the L-infinity norm in the subspace ', ... 
                'iterations. \n' ] );
        end
    end
    %
    if sys.des == 0 && doLoewner ~= 1
        if isfield( opt, 'eigopt' )
            %
            % eigopt.bounds: lower and upper bounds for eigopt.
            %
            if isfield( opt.eigopt, 'bounds' )
                bounds.lb = opt.eigopt.bounds( 1 );
                bounds.ub = opt.eigopt.bounds( 2 );
            else
                error( 'The variable eigopt.bounds is not specified.' );
            end
            %
            % eigopt.gamma: a lower bound for the second derivative of the  
            %               function -sigma_max(G(i*w)) where w is within 
            %               opt.eigopt.bounds.
            %
            if isfield( opt.eigopt, 'gamma' )
                gamma = opt.eigopt.gamma;
            else
                error( 'The variable eigopt.gamma is not specified.' );
            end
        else
            error( 'The variable eigopt is not specified.' );
        end
    end
    %
    switch sys.fct.type
        case 'd'
            %
            % If no delay is defined then the delay parameter is set to
            % zero.
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
            if dB ~= length( sys.fct.b ) || dC ~= length( sys.fct.c ) ...
                    || dA ~= length( sys.fct.a ) ... 
                    || dE ~= length( sys.fct.e ) ...
                    || dD ~= length( sys.fct.d )
                 error( [ 'The number of the delay entries in ', ...
                     'sys.fct does not match the dimensions of ', ...
                     'sys.A, sys.B, sys.C, sys.D, or sys.E.' ] );
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
                error( [ 'The dimensions of the functions in sys.fct ', ...
                    'do not match the dimensions of sys.A, sys.B, ', ...
                    'sys.C, sys.D, or sys.E.' ] );
            end
    end
end
%
%% COMPUTATIONS.
%
if doLoewner
    %
    % Determine the initial interpolation points.
    %
    if isempty( w ) 
        error( 'The variable opt.initialPoints is not specified.' );
    else
        smpl = 1i*w( : );
    end
    %
    % With the frequency points obtained above, a rational approximation is
    % obtained using the Loewner framework.
    %
    if maxSing
        %
        % Run tangential Loewner. Take the singular vectors associated with
        % the largest singular values as tangential directions.
        %
        maxDir = 1;
        optionsLoewner.tan = 1;
        optionsLoewner.makeReal = 1;
        optionsLoewner.trunctol = orthtol;
    else
        %
        % Run block Loewner.
        %
        maxDir = 0;
        optionsLoewner.tan = 0;
        optionsLoewner.makeReal = 1;
        optionsLoewner.trunctol = orthtol;
    end
    %
    if fnHandle == 1
        data = sampleG( sys, smpl, maxDir );
    else
        data = sampleFun( sys, smpl, maxDir );
    end
    [ pars, cells ] = getLoewnerSystem( data, optionsLoewner );
    zOld = -1e16;
    fOld = -1e16;
    %
    % The L-infinity norm and a new interpolation point are computed.
    %
    if boydbalak == 0
        %
        % Use AB13HD.f/linorm_h.F.
        %
        [ f, z ] = linorm_h( full( pars.A ), full( pars.E ), ...
            full( pars.B ), full( pars.C ), full ( pars.D ), 1, 0, 1, 0 );
        %
        if f( 2 ) == 0 || z( 2 ) == 0
            error( [ 'linorm_h: The reduced function is not in ', ...
                'L-infinity.' ] );
        else
            f = f( 1 )/f( 2 );
            z = z( 1 )/z( 2 );
        end
    else
        %
        % Use the MATLAB function 'norm'.
        %
        parssys = dss( pars.A, pars.B, pars.C, pars.D, pars.E );
        [ f, z ] = norm( parssys, Inf );
    end
    %
    if prtlevel > 0
        %
        % Print info.
        %
        subSize = size( pars.A, 1 );
        fprintf( [ 'Current iteration: 0 \t Current interpolation ', ... 
            'frequency: %.6f \t Model order: %i \n' ], z, subSize );
    end
    %
    if keepSubspaces
        optionsLoewner.add = 1;
    else
        optionsLoewner.add = 0;
    end
    %
    % Main loop: Run the algorithm until the relative accuracy or the
    % maximum number of iterations is reached.
    %
    iter = 0;
    %
    while ( abs( z-zOld ) > tolz*( 0.5*( abs( z )+abs( zOld ) ) ) && ...
            ( ( abs( f-fOld ) > tolf*( 0.5*( abs( f )+abs( fOld ) ) ) ) ...
            || abs( z-zOld ) > 1e3*tolz*( 0.5*( abs( z )+ ...
            abs( zOld ) ) ) ) && iter < maxit )
        %
        % The interpolation order is changed if interpolation point is
        % zero, since all real interpolation points are treated in a
        % different way than the complex ones in the Loewner matrix
        % framework.
        %
        if z == 0
            optionsLoewner.add = 0;
        end
        iter = iter + 1;
        zOld = z;
        fOld = f;
        %
        % Sample at the critical frequency and add to the measured data.
        %
        if fnHandle == 1
            dataAdd = sampleG( sys, 1i*z, maxDir );
        else
            dataAdd = sampleFun( sys, 1i*z, maxDir );
        end
        data.smpl = [ data.smpl; dataAdd.smpl ];
        data.G = [ data.G; dataAdd.G ]; 
        data.dG = [ data.dG; dataAdd.dG ];
        %
        if maxSing
            data.c = [ data.c, dataAdd.c ];
            data.b = [ data.b, dataAdd.b ];
        end
        %
        % Obtain the updated rational approximation.
        %
        data.Prev = cells;
        [ pars, cells ]  = getLoewnerSystem( data, optionsLoewner );
        %
        % Determine the L-infinity norm of the reduced function and a new
        % interpolation point.
        %
        if boydbalak == 0
            %
            % Use AB13HD.f/linorm_h.F.
            %
            [ f, z ] = linorm_h( full( pars.A ), full( pars.E ), ...
                full( pars.B ), full( pars.C ), full ( pars.D ), 1, 0, ...
                1, 0 );
            %
            if f( 2 ) == 0 || z( 2 ) == 0
                error( [ 'linorm_h: The reduced function is not in ', ...
                    'L-infinity.' ] );
            else
                f = f( 1 )/f( 2 );
                z = z( 1 )/z( 2 );
            end
        else
            %
            % Use the MATLAB function 'norm'.
            %
            parssys = dss( pars.A, pars.B, pars.C, pars.D, pars.E );
            [ f, z ] = norm( parssys, Inf );
        end
        %
        if prtlevel > 0
            %
            % Print info.
            %
            subSize = size( pars.A, 1 );
            fprintf( [ 'Current iteration: %i \t Current ', ...
                'interpolation frequency: %.6f \t Model order: %i ', ...
                '\n' ], iter, z, subSize );
        end
        %
        if keepSubspaces == 0
            %
            % If only the last two iterations should be used, older data is
            % deleted.
            %
            data.smpl = 1i*zOld;
            data.G = data.G{ end };
            data.dG = data.dG{ end };
        end
    end        
else
    %
    % Projection is used to compute the L-infinity norm (as in v1.0).
    %
    % Get the initial projection subspaces.
    %
    U_initial = cell( 1, length( w ) );
    V_initial = cell( 1, length( w ) );
    %
    for i = 1:length( w )
        if maxSing
            [ V_initial{ i }, U_initial{ i } ] = getSubspace( sys, ...
                w( i ), 1 );
        else
            [ V_initial{ i }, U_initial{ i } ] = getSubspace( sys, ...
                w( i ), 0 );
        end
    end
    %
    if prtlevel > 1
        %
        % Print info.
        %
        fprintf( [ 'The initial subspaces for interpolation have ', ...
            'been determined.\n' ] );       
    end
    %
    % Convert U and V into full matrices.
    %
    for i = 1:length( w )
        U_initial{ i } = full( U_initial{ i } );
        V_initial{ i } = full( V_initial{ i } );
    end
    %
    U = [ U_initial{ : } ];
    V = [ V_initial{ : } ];
    %
    % Construct orthonormal or biorthonormal bases for the subspaces given
    % by U and V.
    %
    [ U, U2 ] = orthnorm( U, orthtol );
    [ V, V2 ] = orthnorm( V, orthtol );
    nu = size( U, 2 );
    nv = size( V, 2 );
    %
    if ( nu > nv )
        d = nu - nv;
        V = [ V, V2( :, 1:d ) ];
    elseif ( nv > nu )
        d = nv - nu;
        U = [ U, U2( :, 1:d ) ];
    end
    %
    if biorth
        %
        % Construct biorthonormal bases. 
        %
        [ U, V ] = biorthnorm( U, V, orthtol );
    end
    %
    % Initially, the left and right subspaces are Col( U ), Col( V ).
    %
    U_init = U;
    V_init = V;
    Vprev = V;
    Uprev = U;
    V_iteration = [];
    U_iteration = [];
    %
    % Reduce the function with the subspaces U and V.
    %
    optionsRedSys.makeRegular = 1;
    optionsRedSys.trunctol = orthtol;
    pars = reduceSystem( sys, U, V, optionsRedSys );
    %
    % Initialization: zOld keeps the optimal frequency of the previous
    % reduced function.
    %
    zOld = -1e16;
    fOld = -1e16;
    %
    % Choose the algorithm.
    %
    switch algo
        case 'e'
            %
            % Use eigopt.
            %
            % Set the parameters for eigopt.
            %
            pars.gamma = gamma;
            pars.minmax = 1;
            pars.isprint = 0;
            %
            % Call eigopt.
            %
            [ f, z ] = eigopt( 'Linf', bounds, pars );
        case 'b'
            %
            % Use the Boyd-Balakrishnan algorithm.
            %
            if boydbalak == 0
                %
                % Use AB13HD.f/linorm_h.F.
                %
                [ f, z ] = linorm_h( full( pars.A ), full( pars.E ), ...
                    full( pars.B ), full( pars.C ), full ( pars.D ), 1, ...
                    0, 1, 0 );
                %
                if f( 2 ) == 0 || z( 2 ) == 0
                    error( [ 'linorm_h: The reduced function is not ', ...
                        'in L-infinity.' ] );
                else
                    f = f( 1 )/f( 2 );
                    z = z( 1 )/z( 2 );
                end
            else
                %
                % Use the MATLAB function 'norm'.
                %
                parssys = dss( pars.A, pars.B, pars.C, pars.D, pars.E );
                [ f, z ] = norm( parssys, Inf );
            end
    end
    %
    if prtlevel > 0
        %
        % Print info.
        %
        [ ~, subSize ] = size( V );
        fprintf( [ 'Current iteration: 0 \t Current interpolation ', ...
            'frequency: %.6f \t Subspace dimension: %i \n' ], z, subSize );
    end
    %
    % Main loop: Run the algorithm until the relative accuracy or the
    % maximum number of iterations is reached.
    %
    iter = 0;
    %
    while ( abs( z-zOld ) > tolz*( 0.5*( abs( z )+abs( zOld ) ) ) && ...
            ( ( abs( f-fOld ) > tolf*( 0.5*( abs( f )+abs( fOld ) ) ) ) ... 
            || abs( z-zOld ) > 1e3*tolz*( 0.5*( abs( z )+ ...
            abs( zOld ) ) ) ) && iter < maxit )
        %
        zOld = z;
        fOld = f;
        %
        % Compute the subspace at the new interpolation point.
        %
        if maxSing
            [ v, u ] = getSubspace( sys, z, 1 ); 
        else
            [ v, u ] = getSubspace( sys, z, 0 );
        end
        %
        % Convert U and V into full matrices.
        %
        u = full( u );
        v = full( v );
        %
        % The interpolation points and subspaces are kept as specified in 
        % keepSubspaces.
        %
        U_iteration = [ U_iteration u ];
        V_iteration = [ V_iteration v ];
        %
        if keepSubspaces
            %
            % All subspaces are kept for the next iteration.
            %
            U = [ U_init, U_iteration ];
            V = [ V_init, V_iteration ];   
        else
            %
            % Only the current subspace is kept for the next iteration.
            %
            U = [ Uprev, u ];
            V = [ Vprev, v ];
            Uprev = u;
            Vprev = v;
        end
        % 
        % Construct orthonormal or biorthonormal bases for the subspaces
        % given by U and V.
        %
        [ U, U2 ] = orthnorm( U, orthtol );
        [ V, V2 ] = orthnorm( V, orthtol );
        nu = size( U, 2 );
        nv = size( V, 2 );
        %
        if ( nu > nv )
            d = nu - nv;
            V = [ V, V2( :, 1:d ) ];
        elseif ( nv > nu )
            d = nv - nu;
            U = [ U, U2( :, 1:d ) ];
        end
        %
        if biorth
            %
            % Construct biorthonormal bases. 
            %
            [ U, V ] = biorthnorm( U, V, orthtol );
        end
        %
        % Determine the reduced function.
        %
        pars = reduceSystem( sys, U, V, optionsRedSys );
        switch algo
            case 'e'
                %   
                % Call eigopt to calculate the L-infinity norm of the
                % reduced function.
                %
                % Set the eigopt parameters.
                %
                pars.gamma = gamma;
                pars.minmax = 1;
                %
                [ f, z ] = eigopt( 'Linf', bounds, pars );
            case 'b'
                %
                % Use the Boyd-Balakrishnan algorithm.
                %
                if boydbalak == 0
                    %
                    % Use AB13HD.f/linorm_h.F.
                    %
                    [ f, z ] = linorm_h( full( pars.A ), ...
                        full( pars.E ), full( pars.B ), full( pars.C ), ...
                        full ( pars.D ), 1, 0, 1, 0 );
                    %
                    if f( 2 ) == 0 || z( 2 ) == 0
                        error( [ 'linorm_h: The reduced function is ', ...
                            'not in L-infinity.' ] );
                    else
                        f = f( 1 )/f( 2 );
                        z = z( 1 )/z( 2 );
                    end
                else
                    %
                    % Use the MATLAB function 'norm'.
                    %
                    parssys = dss( pars.A, pars.B, pars.C, pars.D, pars.E );
                    [ f, z ] = norm( parssys, Inf );
                end
        end
        %
        % Update the number of iterations.
        %
        iter = iter + 1;
        %
        if prtlevel > 0
            %
            % Print info.
            %
            [ ~, subSize ] = size( V );
            fprintf( [ 'Current iteration: %i \t Current ', ...
                'interpolation frequency: %g \t Subspace dimension: ', ...
                '%i \n' ], iter, z, subSize );      
        end
    end
end
%
%% OUTPUT INFORMATION.
%
info.time = toc;
info.iterations = iter;
info.reducedSystem = pars;
info.finaltolz = abs( z-zOld )/( 0.5*( abs( z )+abs( zOld ) ) );
info.finaltolf = abs( f-fOld )/( 0.5*( abs( f )+abs( fOld ) ) );
%
if iter >= maxit    
    info.error = 1;
    if ( prtlevel > 0 )
        %
        % Print info.
        %
        warning( [ 'The algorithm has terminated by reaching the ', ...
            'maximum number of iterations.' ] );
    end
end
%
% Check which of the termination criteria are satisfied.
%
if abs( z-zOld ) < tolz*( 0.5*( abs( z )+abs( zOld ) ) ) && ...
        abs( f-fOld ) < tolf*( 0.5*( abs( f )+abs( fOld ) ) )
    %
    % Both the optimal frequency z and the L-infinity norm f have 
    % converged.
    %
    info.termcrit = 0;
    %
    if prtlevel > 0
        %
        % Print info.
        %
        fprintf( [ 'Both the optimal frequency and the L-infinity ',... 
            'norm have converged.\n' ] );
    end
elseif abs( z-zOld ) < tolz*( 0.5*( abs( z )+abs( zOld ) ) )
    %
    % Only the optimal frequency z has converged, but the L-infinity norm
    % may be inaccurate.
    %
    info.termcrit = 1;
    %
    if prtlevel > 0
        %
        % Print info.
        %
        fprintf( [ 'The optimal frequency has converged, but the ', ...
            'L-infinity norm may be inaccurate.\n' ] );
    end
elseif ( ( abs( f-fOld ) < tolf*( 0.5*( abs( f )+abs( fOld ) ) ) ) ... 
        && abs( z-zOld ) < 1e3*tolz*( 0.5*( abs( z )+abs( zOld ) ) ) )
    %
    % Only the L-infinity norm f has converged, but the optimal frequency 
    % may be inaccurate.
    %
    info.termcrit = 2;
    %
    if prtlevel > 0
        %
        % Print info.
        %
        fprintf( [ 'The L-infinity norm has converged, but the ', ...
            'optimal frequency may be inaccurate.\n' ] );
    end
else
    %
    % The algorithm has terminated by reaching the maximum number of
    % iterations.
    %
    info.termcrit = 3;
end
%
if prtlevel > 0
    %
    % Print info.
    %
    fprintf( [ 'The L-infinity norm is attained for fopt = %g. The ', ...
        'norm value is %g.\n' ], z, f );
    fprintf( 'The runtime is %g seconds.\n', info.time );
end
%
return
