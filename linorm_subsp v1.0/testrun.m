%
% Test script to run a set of example systems with LINORM_SUBSP.
%
% The example data can be downloaded from
%
%    https://sites.google.com/site/rommes/software; and
%    http://www.slicot.org/index.php?site=benchmodred.
%
% REFERENCES:
%
% [1] N. Aliyev, P. Benner, E. Mengi and M. Voigt. Large-scale computation
%     of H-infinity norms by a greedy subspace method.
%
% [2] P. Benner and M. Voigt. L-infinity-norm computation for 
%     continuous-time descriptor systems using structured matrix pencils, 
%     IEEE Trans. Automat. Control, 57(1):233-238, 2012.
%
% AUTHORS:
%
% Paul Schwerdtner and Matthias Voigt, Technische Universitaet Berlin,
% Institut fuer Mathematik, Berlin, Germany.
%
% 25/04/2017.
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
delete diary.dia
diary( 'diary.dia' );
fprintf( '===========================================================\n' );
fprintf( 'TEST SERIES FOR linorm_subsp.m\n' );
fprintf( '===========================================================\n' );
fprintf( 'Date: %s.\n', date );
fprintf( 'Version: %s.\n', version );
fprintf( 'Computer: %s.\n', computer );
%% EXAMPLE 1
clear
load testexample/build
options.prtlevel = 2;
options.initialPoints = linspace( 0, 100, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 1: build\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( speye( 48 ), A, B, C, 0 );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 2
clear
load testexample/pde
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 2: pde\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct(speye( 84 ), A, B, C, 0 );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 3
clear
load testexample/CDplayer
options.prtlevel = 2;
options.initialPoints = linspace( 0, 100, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 3: CDplayer\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( speye( 120 ), A, B, C, zeros( 2 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 4
clear
load testexample/iss
options.prtlevel = 2;
options.initialPoints = linspace( 0, 100, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 4: iss\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( speye( 270 ), A, B, C, zeros( 3 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 5
clear
load testexample/beam
options.prtlevel = 2;
options.initialPoints = linspace( 0, 1, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 5: beam\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( speye( 348 ), A, B, C, 0 );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 6
clear
load testexample/S10PI_n1
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 6: S10PI_n1\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, 0 );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 7
clear
load testexample/S20PI_n1
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 7: S20PI_n1\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, 0 );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 8
clear
load testexample/S40PI_n1
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10  );
fprintf( '===========================================================\n' );
fprintf( 'Example 8: S40PI_n1\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, 0 );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 9
clear
load testexample/S80PI_n1
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 9: S80PI_n1\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, 0 );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 10
clear
load testexample/M10PI_n1
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 10: M10PI_n1\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, zeros( 3 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 11
clear
load testexample/M20PI_n1
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 11: M20PI_n1\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, zeros( 3 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 12
clear
load testexample/M40PI_n1
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 12: M40PI_n1\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, zeros( 3 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 13
clear
load testexample/M80PI_n1
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 13: M80PI_n1\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, zeros( 3 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 14
clear
load testexample/peec
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10, 80 );
fprintf( '===========================================================\n' );
fprintf( 'Example 14: peec\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( E, A, B, C, 0 );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 15
clear
load testexample/S10PI_n
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 15: S10PI_n\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, zeros( 1 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 16
clear
load testexample/S20PI_n
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 16: S20PI_n\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, zeros( 1 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 17
clear
load testexample/S40PI_n
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 17: S40PI_n\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, zeros( 1 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 18
clear
load testexample/S80PI_n
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 18: S80PI_n\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, zeros( 1 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 19
clear
load testexample/M10PI_n
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
options.keepSubspaces = 0;
fprintf( '===========================================================\n' );
fprintf( 'Example 19: M10PI_n\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, zeros( 3 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 20
clear
load testexample/M20PI_n
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 20: M20PI_n\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, zeros( 3 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 21
clear
load testexample/M40PI_n
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 21: M40PI_n\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, zeros( 3 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 22
clear
load testexample/M80PI_n
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 22: M80PI_n\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( T, A, B, C, zeros( 3 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 23
clear
load testexample/bips98_606
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 23: bips98_606\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( E, A, b, c, zeros( 4 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 24
clear
load testexample/bips98_1142
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 24: bips98_1142\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( E, A, b, c, zeros( 4 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 25
clear
load testexample/bips98_1450
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 25: bips98_1450\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( E, A, b, c, zeros( 4 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 26
clear
load testexample/bips07_1693
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 26: bips07_1693\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( E, A, b, c, zeros( 4 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 27
clear
load testexample/bips07_1998
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 1e-4, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 27: bips07_1998\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( E, A, b, c, zeros( 4 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 28
clear
load testexample/bips07_2476
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 1e-4, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 28: bips07_2476\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( E, A, b, c, zeros( 4 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 29
clear
load testexample/bips07_3078
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 1e-4, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 29: bips07_3078\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( E, A, b, c, zeros( 4 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 30
clear
load testexample/xingo_afonso_itaipu
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 30: xingo_afonso_itaipu\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( E, A, B, C', zeros( 1 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 31
clear
load testexample/mimo8x8_system
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
options.maxSing = 1;
fprintf( '===========================================================\n' );
fprintf( 'Example 31: mimo8x8_system\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( E, A, b, c', zeros( 8 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 32
clear
load testexample/mimo28x28_system
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
options.maxSing = 1;
fprintf( '===========================================================\n' );
fprintf( 'Example 32: mimo28x28_system\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( E, A, b, c', zeros( 28 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 33
clear
load testexample/mimo46x46_system
% Set options
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
options.maxSing = 1;
fprintf( '===========================================================\n' );
fprintf( 'Example 33: mimo46x46_system\n' );
fprintf( '===========================================================\n' );
sys = makeSysStruct( E, A, B, C', zeros( 46 ) );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 34
% Delay Model
%
% Add path for eigopt.
%
addpath( 'eigopt' );
%
clear
addpath( 'testexample' );
sys = delay_model( 500, 5, 0.01, 1 );
%
options.prtlevel = 2;
options.eigopt.bounds = [ 0, 50 ];
options.eigopt.gamma = -100;
options.initialPoints = linspace( 0, 50, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 34: Delay Model\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%
diary off
