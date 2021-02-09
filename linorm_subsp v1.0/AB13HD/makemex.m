function makemex
%
% Function for generating gateway function from MEX file
%
%% Set flags depending on machine architecture.
%
flags = '';
is64 = ( ~isempty( strfind( computer, '64' ) ) );
if is64
    %
    % 64-bit MATLAB
    %
    flags = 'FFLAGS="$FFLAGS -fPIC -fno-omit-frame-pointer -fdefault-integer-8" -largeArrayDims';
end
%
%% Set location of the libraries. You might want to modify these. 
%
libslicot = '../../libs/slicot64.a';
%
%% Mexing.
%
fprintf( 'mex %s linorm_h.F %s %s -lmwlapack -lmwblas\n', flags, ...
    'AB13HD.o', libslicot );
eval( sprintf( 'mex %s linorm_h.F %s %s -lmwlapack -lmwblas\n', flags, ...
    'AB13HD.o', libslicot ) );
