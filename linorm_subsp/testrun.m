%
% Test script to run a set of example systems with LINORM_SUBSP v1.1.
%
% The example data can be downloaded from
%
%    https://sites.google.com/site/rommes/software;
%    http://slicot.org/20-site/126-benchmark-examples-for-model-reduction;
%    http://www.maths.manchester.ac.uk/our-research/research-groups/
%       numerical-analysis-and-scientific-computing/numerical-analysis/
%       software/nlevp/;
%    http://twr.cs.kuleuven.be/research/software/delay-control/hinf/.
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
% Paul Schwerdtner and Matthias Voigt, Technische Universitaet Berlin,
% Institut fuer Mathematik, Berlin, Germany.
%
% 25/04/2017.
%
% REVISIONS:
%
% Paul Schwerdtner, Matthias Voigt, 09/2017, 12/2017.
% Matthias Voigt, 06/2018.
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
delete diary.dia
diary( 'diary.dia' );
fprintf( '===========================================================\n' );
fprintf( 'TEST SERIES FOR linorm_subsp.m\n' );
fprintf( '===========================================================\n' );
fprintf( 'Date: %s.\n', date );
fprintf( 'Version: %s.\n', version );
fprintf( 'Computer: %s.\n', computer );
addpath( 'eigopt' );
addpath( 'Loewner' );
%% EXAMPLE 1
clear
load ../testexamples/rational/build
options.prtlevel = 2;
options.initialPoints = linspace( 0, 100, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 1: build\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 2
clear
load ../testexamples/rational/pde
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 2: pde\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 3
clear
load ../testexamples/rational/CDplayer
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 3: CDplayer\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 4
clear
load ../testexamples/rational/iss
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 4: iss\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 5
clear
load ../testexamples/rational/beam
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 5: beam\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 6
clear
load ../testexamples/rational/S10PI_n1
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 6: S10PI_n1\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 7
clear
load ../testexamples/rational/S20PI_n1
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 7: S20PI_n1\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 8
clear
load ../testexamples/rational/S40PI_n1
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 8: S40PI_n1\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 9
clear
load ../testexamples/rational/S80PI_n1
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 9: S80PI_n1\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 10
clear
load ../testexamples/rational/M10PI_n1
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 10: M10PI_n1\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 11
clear
load ../testexamples/rational/M20PI_n1
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 11: M20PI_n1\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 12
clear
load ../testexamples/rational/M40PI_n1
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 12: M40PI_n1\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 13
clear
load ../testexamples/rational/M80PI_n1
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 13: M80PI_n1\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 14
clear
load ../testexamples/rational/peec
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10, 80 );
fprintf( '===========================================================\n' );
fprintf( 'Example 14: peec\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 15
clear
load ../testexamples/rational/S10PI_n
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 15: S10PI_n\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 16
clear
load ../testexamples/rational/S20PI_n
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 16: S20PI_n\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 17
clear
load ../testexamples/rational/S40PI_n
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 17: S40PI_n\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 18
clear
load ../testexamples/rational/S80PI_n
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 18: S80PI_n\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 19
clear
load ../testexamples/rational/M10PI_n
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 19: M10PI_n\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 20
clear
load ../testexamples/rational/M20PI_n
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 20: M20PI_n\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 21
clear
load ../testexamples/rational/M40PI_n
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 21: M40PI_n\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 22
clear
load ../testexamples/rational/M80PI_n
options.prtlevel = 2;
options.initialPoints = linspace( 0.1, 10000, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 22: M80PI_n\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 23
clear
load ../testexamples/rational/bips98_606
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 23: bips98_606\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 24
clear
load ../testexamples/rational/bips98_1142
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 24: bips98_1142\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 25
clear
load ../testexamples/rational/bips98_1450
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 25: bips98_1450\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 26
clear
load ../testexamples/rational/bips07_1693
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 26: bips07_1693\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 27
clear
load ../testexamples/rational/bips07_1998
options.prtlevel = 2;
options.initialPoints = linspace( 1e-4, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 27: bips07_1998\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 28
clear
load ../testexamples/rational/bips07_2476
options.prtlevel = 2;
options.initialPoints = linspace( 1e-4, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 28: bips07_2476\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 29
clear
load ../testexamples/rational/bips07_3078
options.prtlevel = 2;
options.initialPoints = linspace( 1e-4, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 29: bips07_3078\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 30
clear
load ../testexamples/rational/xingo_afonso_itaipu
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 30: xingo_afonso_itaipu\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 31
clear
load ../testexamples/rational/mimo8x8_system
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
options.maxSing = 1;
fprintf( '===========================================================\n' );
fprintf( 'Example 31: mimo8x8_system\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 32
clear
load ../testexamples/rational/mimo28x28_system
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
options.maxSing = 1;
fprintf( '===========================================================\n' );
fprintf( 'Example 32: mimo28x28_system\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 33
clear
load ../testexamples/rational/mimo46x46_system
options.prtlevel = 2;
options.initialPoints = linspace( 0, 10, 10 );
options.maxSing = 1;
fprintf( '===========================================================\n' );
fprintf( 'Example 33: mimo46x46_system\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 34a
clear
addpath( '../testexamples/nonrational' );
sys = delay_model( 500, 5, 0.01, 1 );
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 50, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 34: Delay Model (via Loewner)\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 34b
addpath( '../testexamples/nonrational' );
sys = delay_model( 500, 5, 0.01, 1 );
options.doLoewner = 0;
options.prtlevel = 2;
options.eigopt.bounds = [ 0, 50 ];
options.eigopt.gamma = -100;
options.initialPoints = linspace( 0, 50, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 34: Delay Model (via eigopt)\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 35
clear
load ../testexamples/nonrational/butterfly
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 35: butterfly (with Loewner)\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 36
clear
load ../testexamples/nonrational/dirac
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 36: dirac (with Loewner)\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info] = linorm_subsp( sys, options );
%% EXAMPLE 37
clear
load ../testexamples/nonrational/gen_hyper2
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 100, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 37: gen_hyper2 (with Loewner)\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 38
clear
load ../testexamples/nonrational/gen_tantipal2
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 38: gen_tantipal2 (with Loewner)\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 39
clear
load ../testexamples/nonrational/gen_tpal2
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 39: gen_tpal2 (with Loewner)\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 40
clear
load ../testexamples/nonrational/hadeler
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 40: hadeler (with Loewner)\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 41
clear
load ../testexamples/nonrational/loaded_string
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 100, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 41: loaded_string (with Loewner)\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 42
clear
load ../testexamples/nonrational/sleeper
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 42: sleeper (with Loewner)\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info] = linorm_subsp( sys, options );
%% EXAMPLE 43
clear
load ../testexamples/nonrational/spring
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 43: spring (with Loewner)\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 44
clear
load ../testexamples/nonrational/spring_dashpot
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 0.01, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 44: spring_dashpot (with Loewner)\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 45
clear
load ../testexamples/nonrational/wiresaw2
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 100, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 45: wiresaw2 (with Loewner)\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 46
clear
load ../testexamples/nonrational/hinfn_ex1
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 1, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 46: hinfn_ex1\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 47
clear
load ../testexamples/nonrational/hinfn_ex2
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 1, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 47: hinfn_ex2\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 48
clear
load ../testexamples/nonrational/hinfn_ex3
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 48: hinfn_ex3\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 49
clear
load ../testexamples/nonrational/hinfn_ex4
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 10, 20 );
fprintf( '===========================================================\n' );
fprintf( 'Example 49: hinfn_ex4\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 50
clear
load ../testexamples/nonrational/hinfn_ex5
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 100, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 50: hinfn_ex5\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 51
clear
load ../testexamples/nonrational/hinfn_ex6
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 100, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 51: hinfn_ex6\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 52
clear
load ../testexamples/nonrational/hinfn_ex7
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 100, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 52: hinfn_ex7\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 53
clear
load ../testexamples/nonrational/hinfn_ex8
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 100, 20 );
fprintf( '===========================================================\n' );
fprintf( 'Example 53: hinfn_ex8\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 54
clear
load ../testexamples/nonrational/hinfn_ex9
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 100, 20 );
fprintf( '===========================================================\n' );
fprintf( 'Example 54: hinfn_ex9\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 55
clear
load ../testexamples/nonrational/hinfn_ex10
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 10, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 55: hinfn_ex10\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 56
clear
load ../testexamples/nonrational/hinfn_ex11
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 100, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 56: hinfn_ex11\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%% EXAMPLE 57
clear
load ../testexamples/nonrational/hinfn_ex12
options.prtlevel = 2;
options.doLoewner = 1;
options.initialPoints = linspace( 0, 100, 10 );
fprintf( '===========================================================\n' );
fprintf( 'Example 57: hinfn_ex12\n' );
fprintf( '===========================================================\n' );
[ normval, fopt, info ] = linorm_subsp( sys, options );
%
diary off