function [ sys, cells ] = getLoewnerSystem( data, opt )
%
%% PURPOSE:
%
% To construct a state-space realization of a linear system with a transfer
% function matrix G(s) = C*(s*E - A)^(-1)*B that matches some given data.
% If opt.tan == 0, then the data (smpl, G, dG) and G(s) fulfill the
% interpolation conditions
%
%   G(smpl(k)) = G{k} and G'(smpl(k)) = dG{k} for k = 1, ..., length(smpl);
%
% if opt.tan == 1, then the data (smpl, G, dG, b, c) and G(s) fulfill the
% tangential interpolation conditions
%
%   G(smpl(k))*b(:,k) = G{k}*b(:,k), c(:,k)'*G(smpl(k)) = c(:,k)'*G{k}, and
%   
%   c(:,k)'*G'(smpl(k))*b(:,k) = c(:,k)'*dG{k}*b(:,k) 
%
%   for k = 1, ..., length(smpl).
%
% Here smpl is a set of interpolation points, G{k} and dG{k} are function
% value and derivative information of the function G(s) at the point 
% smpl(k) and b(:,k) and c(:,k) are the k-th right and left tangential
% directions, respectively.
%
% The realization is constructed using the Loewner framework, see [1].
%
% ARGUMENTS:
%
% Inputs:
%
% data    : Struct containing the following information:
%    data.smpl            : array containing the interpolation points.
%    data.G               : cell array containing the function value
%                           information such that G(smpl(k)) = data.G{k}.
%    data.dG              : cell array containing the derivative 
%                           information such that G'(smpl(k)) = data.dG{k}.
%    data.Prev.A, 
%    data.Prev.E          : cell arrays containing the realization matrices
%                           A and E of the function G(s) in (1) if only the
%                           interpolation points smpl(1), ..., 
%                           smpl(length(smpl))-1 are used. The matrices 
%                           E{ i, j } and A{ i, j } contain the ij-th block
%                           of the Loewner and the shifted Loewner matrix, 
%                           respectively, see [1] for details. This data 
%                           can be made use of, if only one interpolation
%                           point is added and avoids the complete
%                           recomputation of the realization matrices. This
%                           data is only used if opt.add == 0.
%    data.b, data.c       : the right and left tangential interpolation 
%                           directions. Here data.b(:,k) and data.c(:,k) 
%                           are the k-th right and left tangential
%                           directions, respectively. This data is only
%                           used if opt.tan == 1.
% opt     : Struct containing options (use default values if empty). The
%           following options can be specified:
%    opt.makeReal         : specifies whether the realization matrices of
%                           G(s) should be made real as follows:
%                           = 0: the realization matrices should not be
%                                made real;
%                           = 1: the realization matrices should be made
%                                real. This increases the state-space
%                                dimension compared to a complex
%                                realization;
%                           (default = 0).
%    opt.add              : specifies whether only one additional
%                           interpolation point is added as follows:
%                           = 0: all realization matrices are computed from
%                                scratch;
%                           = 1: only one additional interpolation point is
%                                added. The Loewner matrices contained in
%                                data.Prev are extended. This should only
%                                be used if opt.makeReal == 1;
%                           (default = 0).
%    opt.tan              : specifies whether tangential interpolation
%                           should be performed using the tangential
%                           directions provided in data.b and data.c as
%                           follows:
%                           = 0: the full interpolation should be 
%                                performed;
%                           = 1: tangential interpolation should be
%                                performed;
%                           (default = 0).
%    opt.trunctol         : truncation tolerance for removing small
%                           singular values when making the realization
%                           matrices E and A square (default = 1e-12).
%
% Outputs:
%            
% sys     : Struct containing the structural information of the
%           function G(s).
%    sys.A, sys.B, sys.C,
%    sys.E                : the realization matrices A, B, C, and E of the 
%                           function G(s). 
%    sys.D                : = 0.
%    sys.des              : = 1 (G(s) is realized by a descriptor system).
% cells   : Struct containing the data of the Loewner matrix in sys.E and 
%           the shifted Loewner matrix in sys.A to reuse them in future
%           computations.           
%    cells.A, cells.E     : cell arrays containing the realization matrices
%                           A and E of the function G(s) in (1). The
%                           matrices E{ i, j } and A{ i, j } contain the
%                           ij-th block of the Loewner and the shifted
%                           Loewner matrix, respectively, see [1] for
%                           details.
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
% 28/09/2017.
%
% REVISIONS:
%
% Matthias Voigt, 06/2018.
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
%% SET OPTIONS AND CHECK ARGUMENTS.
%
% smpl, G and dG must always be provided.
%
if isfield( data, 'smpl' ) && isfield( data, 'G' ) && isfield( data, 'dG' )
    smpl = data.smpl;
    G = data.G;
    dG = data.dG;
else
    error( [ 'For interpolation with getLoewnerSystem, the ', ... 
        'interpolation points, the function values and the ', ...
        'derivative values must be provided.' ] );
end
%
% opt.makeReal: specifies whether the realization matrices of G(s) should
%               be made real
%
if isfield( opt, 'makeReal' )
    makeReal = opt.makeReal;
else % default
    makeReal = 0;
end
%
% opt.add: specifies whether only one additional interpolation point is 
%          added.
%
if isfield( opt, 'add' )
    add = opt.add;
else % default
    add = 0;
end
%
% opt.tan: specifies whether tangential interpolation should be performed
%          using the tangential directions provided in data.b and data.c.
%
if isfield( opt, 'tan' )
    tan = opt.tan;
else % default
    tan = 0;
end
%
% opt.trunctol: truncation tolerance for removing small singular values
%               when making the realization matrices E and A square.
%
if isfield( opt, 'trunctol' ) 
    trunctol = opt.trunctol;
else % default
    trunctol = 1e-12;
end
%
% Check wether the option specific data is provided.
%
if tan == 1
    if isfield( data, 'b' ) && isfield( data, 'c' )
        b = data.b;
        c = data.c;
    else
        error( [ 'If tangential interpolation is performed ', ... 
            '( opt.tan = 1 ), then the tangential directions in ', ... 
            'data.b and data.c must be provided.' ] );
    end
end
%
if add == 1
    if isfield( data.Prev, 'E' ) && isfield( data.Prev, 'A' )
        E = data.Prev.E;
        A = data.Prev.A;
        if makeReal == 0
            error( [ 'This function is not designed for the case ', ...
                'opt.add == 1 with opt.makeReal == 0.' ] );
        end
    else
        error( 'If opt.add == 1, then data.Prev must be specified.' );
    end
end
%
%% COMPUTATIONS.
%
smpl = smpl( : );
nsmpl = length( smpl );
%
% Check for repetitions in interpolation set and remove them if necessary.
%
if length( unique( smpl ) ) < nsmpl
    [ smpl, ia, ~ ] = unique( smpl, 'stable' );
    G = G( ia );
    add = 0;
end
%
[ p, m ] = size( G{ 1 } );
%
if makeReal == 1
    %
    % Sort for real and complex interpolation points and sort the data in G
    % and dG accordingly.
    %
    smplCmpl = smpl( imag( smpl ) ~= 0 );
    smplCmpl = [ smplCmpl, conj( smplCmpl ) ].';
    smplCmpl = smplCmpl( : );
    %
    % For the complex sample points, add the complex conjugates of G and dG
    % to the interpolation data and sort the cell arrays such that the
    % interpolation data fits the interpolation points as given in smpl.
    %
    GCmpl = G( imag( smpl ) ~= 0 );
    GCmplConj = cellfun( @conj, GCmpl, 'un', 0 );
    GCmpl = [ GCmpl, GCmplConj ].';
    GCmpl = GCmpl( : );
    G = [ G( imag( smpl ) == 0 ); GCmpl ];
    dGCmpl = dG( imag( smpl ) ~= 0 );
    dGCmplConj = cellfun( @conj, dGCmpl, 'un', 0 );
    dGCmpl = [ dGCmpl, dGCmplConj ].';
    dGCmpl = dGCmpl( : );
    dG = [ dG( imag( smpl ) == 0 ); dGCmpl ];
    %
    if tan == 1
        %
        % Also sort the tangential directions c and b.
        %
        ccmpl = c( :, imag( smpl ) ~= 0 );
        ccmpl = ccmpl.';
        ccmpl = [ ccmpl, conj( ccmpl ) ];
        ccmpl = ccmpl.';
        ccmpl = ccmpl( : );
        ccmpl = reshape( ccmpl, p, length( ccmpl )/p );
        c = [ c( :, imag( smpl ) == 0 ), ccmpl ];
        %
        bcmpl = b( :, imag( smpl ) ~= 0 );
        bcmpl = bcmpl.';
        bcmpl = [ bcmpl, conj( bcmpl ) ];
        bcmpl = bcmpl.';
        bcmpl = bcmpl( : );
        bcmpl = reshape( bcmpl, m, length( bcmpl )/m );
        b = [ b( :, imag( smpl ) == 0 ), bcmpl ];
    end
    smpl = [ smpl( imag( smpl ) == 0 ); smplCmpl ];
    nsmpl = length( smpl );
end
%
if add == 0
    %
    % Initialize E and A.
    %
    E = cell( nsmpl );
    A = cell( nsmpl );
end
%
% Contruct off-diagonal blocks of the Loewner and the shifted Loewner 
% matrix as in [1].
%
if add == 1 && tan == 1
    %
    % Only add the last rows and columns.
    %
    for i = 1:nsmpl
        for j = [ nsmpl-1, nsmpl ]
            if j ~= i
                E{ i, j } = c( :, i )'*( G{ i }-G{ j } )*b( :, j )/ ...
                    ( smpl( i )-smpl( j ) );
                A{ i, j } = c( :, i )'*( smpl( i )*G{ i }- ...
                    w( j )*G{ j } )*b( :, j )/( smpl( i )-smpl( j ) );
                E{ j, i } = c( :, j )'*( G{ j }-G{ i } )*b( :, i )/ ...
                    ( smpl( j )-smpl( i ) );
                A{ j, i } = c( :, j )'*( smpl( j )*G{ j }- ...
                    smpl( i )*G{ i } )*b( :, i )/( smpl( j )-smpl( i ) );
            end
        end
    end
elseif add == 1
    %
    % Only add the last block column.
    %
    for i = 1:( nsmpl )
        for j = [ nsmpl-1, nsmpl ]
            if j > i
                E{ i, j } = ( G{ i }-G{ j } )/( smpl( i )-smpl( j ) );
                A{ i, j } = ( smpl( i )*G{ i }-smpl( j )*G{ j } )/ ...
                    ( smpl( i )-smpl( j ) );
            end
        end
    end
    %
    % Make the cells square.
    %
    E = [ E; cell( 1, size( E, 2 ) ) ];
    A = [ A; cell( 1, size( A, 2 ) ) ];
elseif tan == 1
    %
    % Contruct the matrices from scratch. Perform tangential interpolation.
    %
    for i = 1:nsmpl
        for j = 1:nsmpl
            if j ~= i
                E{ i, j } = c( :, i )'*( G{ i }-G{ j } )*b( :, j )/ ...
                    ( smpl( i )-smpl( j ) );
                A{ i, j } = c( :, i )'*(  smpl( i )*G{ i }- ...
                    smpl( j )* G{ j } )*b( :, j )/( smpl( i )-smpl( j ) );
            end
        end
    end
else
    %
    % Construct the matrices from scratch.
    %
    for i = 1:nsmpl
        for j = 1:nsmpl
            if j > i
                E{ i, j } =  ( G{ i }-G{ j } )/( smpl( i )-smpl( j ) );
                A{ i, j } =  ( smpl( i )*G{ i }-smpl( j )*G{ j } )/ ...
                    ( smpl( i )-smpl( j ) );
            end
        end
    end
end
%
cells.E = E;
cells.A = A;
%
% Fill empty the block matrices with zeros and convert to a matrix.
%
tfE = cellfun( 'isempty', E );
tfA = cellfun( 'isempty', A );
%
if tan == 1
    E( tfE ) = { 0 }; 
    A( tfA ) = { 0 };
else
    E( tfE ) = { zeros( p, m ) }; 
    A( tfA ) = { zeros( p, m ) };
end
%
if tan == 1
    sys.E = cell2mat( E );
    sys.A = cell2mat( A );
else
    %
    % In the nontangential case, E and A are symmetric.
    %
    sys.E = cell2mat( E' ) + cell2mat( E );
    sys.A = cell2mat( A' ) + cell2mat( A );
end
%
% Construct the diagonal blocks of the Loewner and the shifted Loewner 
% matrix as in [1].
%
if tan == 1
    for i = 1:nsmpl
        sys.E( i, i ) = c( :, i )'*dG{ i }*b( :, i );
        sys.A( i, i ) = c( :, i )'*( G{ i }+smpl( i )*dG{ i } )*b( :, i );
    end
else
    for i = 1:nsmpl
        sys.E( ( ( i-1 )*p+1 ):( p*( i-1 )+p ), ...
            ( ( i-1 )*m+1 ):( m*( i-1 )+m ) ) = dG{ i };
        sys.A( ( ( i-1 )*p+1 ):( p*( i-1 )+p ), ...
            ( ( i-1 )*m+1 ):( m*( i-1 )+m ) ) = G{ i }+smpl( i )*dG{ i };
    end
end
%
% Assign input and output matrices B and C.
%
if tan == 1
    sys.C = zeros( p, nsmpl );
    sys.B = zeros( nsmpl, m );
    for i = 1:nsmpl
        sys.C( :, i ) = G{ i }*b( :, i );
        sys.B( i, : ) = c( :, i )'*G{ i };
    end
else    
    sys.C = cell2mat( G.' );
    sys.B = cell2mat( G );
end
%
if makeReal == 1
    %
    % Make the system real. 
    %
    if tan == 1
        R = blkdiag( eye( length( smpl( imag( smpl ) == 0 ) ) ), ...
            kron( eye( length( smplCmpl )/2 ), 1/sqrt( 2 )* ...
            [ 1, 1i; 1, -1i ] ) );
        L = blkdiag( eye( length( smpl( imag( smpl ) == 0 ) ) ), ...
            kron( eye( length( smplCmpl )/2 ), 1/sqrt( 2 )* ...
            [ 1, 1i; 1, -1i ] ) );
    else
        R = blkdiag( eye( m*length( smpl( imag( smpl ) == 0 ) ) ), ...
            kron( eye( length( smplCmpl )/2 ), ...
            kron( [ 1, 1i; 1, -1i ], eye( m ) ) ) );
        L = blkdiag( eye( p*length( smpl( imag( smpl ) == 0 ) ) ), ...
            kron( eye( length( smplCmpl )/2 ), ...
            kron( [ 1, 1i; 1, -1i ], eye( p ) ) ) );
    end
    % 
    % The system matrices become real by the transformation, the real 
    % function is only called to make up for rounding errors. 
    % 
    sys.E = real( L'*sys.E*R );
    sys.A = real( L'*sys.A*R );
    sys.C = real( sys.C*R );
    sys.B = real( L'*sys.B );
end
%
% Make A, E square.
%
[ Y, Sig, X ] = svd( 1*sys.E-sys.A, 'econ' ); 
%
% Truncate small singular values and remove the corresponding singular
% vectors.
%

X = X';
Sig = diag( Sig )/Sig( 1 );
k = length( Sig( Sig > trunctol ) );
X = X( 1:k, : );
Y = Y( :, 1:k );
%
% Obtain a system of right dimension.
%
sys.E = -Y'*sys.E*X';
sys.A = -Y'*sys.A*X';
sys.B = Y'*sys.B;
sys.C = sys.C*X';
%
% Indicate that sys is a descriptor system.
%
sys.des = 1;
sys.D = zeros( p, m );
%
return