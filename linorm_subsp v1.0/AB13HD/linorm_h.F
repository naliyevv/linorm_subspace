#include "fintrf.h"
C LINORM_H.F - Gateway function for computation of the L_infinity norm
C              of a continuous-time or discrete-time system, in standard
C              or descriptor form, using SLICOT routine AB13HD.
C
C RELEASE 2.0 of SLICOT Robust Control Toolbox.
C Based on SLICOT RELEASE 5.0, Copyright (c) 2002-2012 NICONET e.V.
C
C Matlab call:
C   [gpeak(,fpeak)] = linorm_h(A,E,B,C,D(,dico,systype,ranke,
C                                    equil(,ckprop,reduce),fpeak0,tol))
C
C Purpose:
C   To compute the L-infinity norm of the continuous-time or
C   discrete-time system, either standard or in the descriptor form,
C
C                                     -1
C        G(lambda) = C*( lambda*E - A ) *B + D .
C
C Input parameters:
C   A       - the n-by-n system state matrix A.
C   E       - the n-by-n descriptor matrix E of the system, or an empty
C             matrix (i.e., with 0 rows and/or columns), in which case
C             E is taken as an identity matrix of order n (standard
C             system). The contents of E are ignored if it is non-empty
C             and systype = 1, but its order is then checked out.
C   B       - the n-by-m system input matrix B.
C   C       - the p-by-n system output matrix C.
C   D       - the p-by-m system matrix D.
C   dico    - (optional) specifies the type of the system:
C             = 1 : continuous-time system;
C             = 2 : discrete-time system.
C             Default: dico = 1.
C   systype - (optional) specifies whether or not the system is of
C             descriptor type:
C             = -1 : descriptor system in compressed form (E = [ T  0 ]
C                                                              [ 0  0 ]
C                    with a square full-rank matrix T);
C             =  0 : descriptor system;
C             =  1 : standard system (E = I).
C             Default: systype = 0.
C   ranke   - (optional) if systype = -1, specifies the rank of the
C             matrix E.
C   equil   - (optional) specifies whether the user wishes to
C             preliminarily equilibrate the system (A,E,B,C) or (A,B,C):
C             = 1 : do not perform equilibration;
C             = 2 : perform equilibration (scaling).
C             Default: equil = 1.
C   ckprop  - (optional) if dico = 1 and systype <> 1, specifies whether
C             the user wishes to check the properness of the transfer
C             function of the descriptor system:
C             = 0 : do not check the properness;
C             = 1 : check the properness.
C             If dico = 2 or systype = 1, this option is ineffective.
C             Default: ckprop = 1.
C   reduce  - (optional) if ckprop = 1, specifies whether the user
C             wishes to reduce the system order, by removing all
C             uncontrollable and unobservable poles before computing
C             the norm:
C             = 0: compute the norm without reducing the order;
C             = 1: reduce the system order.
C             If ckprop = 0, this option is ineffective.
C             Default: reduce = 1.
C
C   fpeak0  - (optional) array of length 2 containing an estimate of the
C             frequency where the gain of the frequency response would
C             achieve its peak value. If fpeak0(2) = 0, the frequency is
C             infinite.
C             Default: fpeak0 = [0; 1].
C   tol     - (optional) tolerance used to set the accuracy in
C             determining the norm.
C             Default: sqrt(epsilon_machine) where epsilon_machine is
C             the relative machine precision.
C
C Output parameters:
C   gpeak   - array of length 2 containing the value of the L_infinity
C             norm of the system. If gpeak(2) = 0, the norm is infinite.
C   fpeak   - (optional) array of length 2 containing the frequency
C             for which the frequency response achieves its peak
C             value gpeak. If fpeak(2) = 0, the frequency is infinite.
C
C Contributor:
C   V. Sima, Research Institute for Informatics, Bucharest, Dec. 2011.
C
C Revisions:
C   M. Voigt, TU Berlin, Apr. 2017, May 2017.
C
C **********************************************************************
C
C
      SUBROUTINE MEXFUNCTION( NLHS, PLHS, NRHS, PRHS )
C
C .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO, HUNDRD
      PARAMETER         ( ONE = 1.0D0, ZERO = 0.0D0, HUNDRD = 1.0D2 )
C
C .. Mex-file interface parameters ..
      mwPointer         PLHS( * ), PRHS( * )
      INTEGER*4         NLHS, NRHS
C
C .. Mex-file integer functions ..
      mwPointer         mxCreateDoubleMatrix, mxGetPr
      INTEGER*4         mxGetM, mxGetN, mxIsNumeric, mxIsComplex
C
C .. Scalar parameters used by SLICOT subroutines ..
      CHARACTER         CKPROP, DICO, EQUIL, JOBD, JOBE, REDUCE
      INTEGER           INFO, LDA, LDB, LDC, LDD, LDE, LDWORK, LZWORK,
     $                  M, N, P, RANKE
C
C .. Allocatable arrays ..
C !Fortran 90/95 (Fixed dimensions should be used with Fortran 77.)
      LOGICAL,          ALLOCATABLE :: BWORK(:)
      INTEGER,          ALLOCATABLE :: IWORK(:)
      mwPointer         A, B, C, D, DWORK, E, ZWORK
C
C .. Local variables and constant dimension arrays ..
      LOGICAL           NODYN, UNITE, WITHD
      CHARACTER*120     TEXT
      INTEGER           ICK, ICOL, IP, ISIZE, ISYS, ISYSD, ITM, ITMP,
     $                  LBWORK, LIWORK, MAXPM, MINPM, NBLK, NN
      DOUBLE PRECISION  BNORM, CNORM, TEMP
      DOUBLE PRECISION  FPEAK(2), FPK(2), GPEAK(2), TOL(4)
C
C .. External functions ..
      LOGICAL           LSAME, MA02HD
      DOUBLE PRECISION  DLAMCH, DLANGE
      EXTERNAL          DLAMCH, DLANGE, LSAME, MA02HD
C
C .. External subroutines ..
      EXTERNAL          AB13HD
C
C .. Intrinsic functions..
      INTRINSIC         MAX, MIN, SQRT
C
C Check for proper number of arguments.
C
      IF ( NRHS.LT.5 ) THEN
         CALL mexErrMsgTxt
     $        ( 'LINORM_H requires at least 5 input arguments' )
      ELSE IF ( NLHS.LT.1 ) THEN
         CALL mexErrMsgTxt
     $        ( 'LINORM_H requires at least 1 output argument' )
      END IF
C
C Check dimensions of input parameters and read/set scalar parameters.
C
C   A(nxn), E(nxn), B(nxm), C(pxn), D(pxm)
C
      N = mxGetM( PRHS(1) )
      M = mxGetN( PRHS(3) )
      P = mxGetM( PRHS(4) )
C
      IF ( mxGetN( PRHS(1) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'A must be a square matrix' )
      END IF
      IF ( mxIsNumeric( PRHS(1) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(1) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'A must be a numeric matrix' )
      END IF
      IF ( mxGetM( PRHS(2) )*mxGetN( PRHS(2) ).NE.0 ) THEN
         IF ( mxGetM( PRHS(2) ).NE.N .OR. 
     $        mxGetN( PRHS(2) ).NE.N ) THEN
            CALL mexErrMsgTxt
     $           ( 'E must be a square matrix of the same order as A' )
         END IF
         IF ( mxIsNumeric( PRHS(2) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(2) ).EQ.1 ) THEN
            CALL mexErrMsgTxt( 'E must be a numeric matrix' )
         END IF
         UNITE = .FALSE.
      ELSE
         UNITE = .TRUE.
      END IF
      IF ( mxGetM( PRHS(3) ).NE.N ) THEN
         CALL mexErrMsgTxt( 'B must have the same row dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(3) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(3) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'B must be a numeric matrix' )
      END IF
      IF ( mxGetN( PRHS(4) ).NE.N ) THEN
         CALL mexErrMsgTxt
     $       ('C must have the same column dimension as A' )
      END IF
      IF ( mxIsNumeric( PRHS(4) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(4) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'C must be a numeric matrix' )
      END IF
      IF ( mxGetM( PRHS(5) ).NE.P ) THEN
         CALL mexErrMsgTxt( 'D must have the same row dimension as C' )
      END IF
      IF ( mxGetN( PRHS(5) ).NE.M ) THEN
         CALL mexErrMsgTxt
     $       ( 'D must have the same column dimension as B' )
      END IF
      IF ( mxIsNumeric( PRHS(5) ).EQ.0 .OR.
     $     mxIsComplex( PRHS(5) ).EQ.1 ) THEN
         CALL mexErrMsgTxt( 'D must be a numeric matrix' )
      END IF
C
C   dico, systype, ranke, equil, ckprop, reduce, fpeak0, tol
C
      IP = 6
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'DICO must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'DICO must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ISYS = TEMP
         IF ( ISYS.NE.1 .AND. ISYS.NE.2 ) THEN
              CALL mexErrMsgTxt( 'DICO must be 1 or 2' )
         END IF
         IF ( ISYS.EQ.2 ) THEN
            DICO = 'D'
         ELSE
            DICO = 'C'
         END IF
         IP = IP + 1
      ELSE
         ISYS = 1
         DICO = 'C'
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'SYSTYPE must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'SYSTYPE must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ISYSD = TEMP
         IF ( ISYSD.LT.-1 .OR. ISYSD.GT.1 ) THEN
              CALL mexErrMsgTxt( 'SYSTYPE must be -1, 0, or 1' )
         END IF
         IF ( ISYSD.EQ.1 .OR. UNITE ) THEN
            ISYSD = 1
            JOBE = 'I'
         ELSE IF ( ISYSD.EQ.-1 ) THEN
            JOBE = 'C'
         ELSE
            JOBE = 'G'
         END IF
         IP = IP + 1
      ELSE
         IF ( UNITE ) THEN
            ISYSD = 1
            JOBE  = 'I'
         ELSE
            ISYSD = 0
            JOBE  = 'G'
         END IF
      END IF
C
      IF ( NRHS.GE.IP .AND. ISYSD.LT.0 ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'RANKE must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'RANKE must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         RANKE = TEMP
         IF ( RANKE.LT.0 .OR. RANKE.GT.N ) THEN
              CALL mexErrMsgTxt
     $             ( 'RANKE must be at least 0 and at most order of A' )
         END IF
         IP = IP + 1
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'EQUIL must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'EQUIL must be an integer scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
         ITMP = TEMP
         IF ( ITMP.NE.1 .AND. ITMP.NE.2 ) THEN
              CALL mexErrMsgTxt( 'EQUIL must be 1 or 2' )
         END IF
         IF ( ITMP.EQ.2 ) THEN
            EQUIL = 'S'
         ELSE
            EQUIL = 'N'
         END IF
         IP = IP + 1
      ELSE
         ITMP  = 1
         EQUIL = 'N'
      END IF
C
      IF ( ISYS.EQ.1 .AND. ISYSD.LE.0 ) THEN
         IF ( NRHS.GE.IP ) THEN
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
                 CALL mexErrMsgTxt( 'CKPROP must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
                 CALL mexErrMsgTxt( 'CKPROP must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
            ICK = TEMP
            IF ( ICK.NE.0 .AND. ICK.NE.1 ) THEN
                 CALL mexErrMsgTxt( 'CKPROP must be 0 or 1' )
            END IF
            IF ( ICK.EQ.1 ) THEN
               CKPROP = 'C'
            ELSE
               CKPROP = 'N'
            END IF
            IP = IP + 1
         ELSE
            ICK    = 1
            CKPROP = 'C'
         END IF
C
         IF ( ICK.EQ.1 .AND. NRHS.GE.IP ) THEN
            IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $           mxGetN( PRHS(IP) ).NE.1 ) THEN
                 CALL mexErrMsgTxt( 'REDUCE must be a scalar' )
            END IF
            IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $           mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
                 CALL mexErrMsgTxt( 'REDUCE must be an integer scalar' )
            END IF
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TEMP, 1 )
            ITM = TEMP
            IF ( ITM.NE.0 .AND. ITM.NE.1 ) THEN
                 CALL mexErrMsgTxt( 'REDUCE must be O or 1' )
            END IF
            IF ( ITM.EQ.1 ) THEN
               REDUCE = 'R'
            ELSE
               REDUCE = 'N'
            END IF
            IP = IP + 1
         ELSE
            REDUCE = 'R'
         END IF
      END IF
C
      FPEAK(1) = ZERO
      FPEAK(2) = ONE
      IF ( NRHS.GE.IP ) THEN
         ISIZE = mxGetM( PRHS(IP) ) * mxGetN( PRHS(IP) )
         IF ( ISIZE.GT.2 ) THEN
            CALL mexErrMsgTxt
     $         ( 'FPEAK0 must be a vector with at most 2 elements' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'FPEAK0 must be a real vector' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), FPK, ISIZE )
         IF ( ISIZE.GT.0 ) FPEAK(1) = FPK(1)
         IF ( ISIZE.GT.1 ) FPEAK(2) = FPK(2)
         IP = IP + 1
      END IF
C
      IF ( NRHS.GE.IP ) THEN
         IF ( mxGetM( PRHS(IP) ).NE.1 .OR.
     $        mxGetN( PRHS(IP) ).NE.1 ) THEN
              CALL mexErrMsgTxt( 'TOL must be a scalar' )
         END IF
         IF ( mxIsNumeric( PRHS(IP) ).EQ.0 .OR.
     $        mxIsComplex( PRHS(IP) ).EQ.1 ) THEN
              CALL mexErrMsgTxt( 'TOL must be a numeric scalar' )
         END IF
         CALL mxCopyPtrToReal8( mxGetPr( PRHS(IP) ), TOL(1), 1 )
         IF ( TOL(1).LE.ZERO ) THEN
              CALL mexErrMsgTxt( 'TOL must be a positive scalar' )
         END IF
         TOL(1) = MAX( HUNDRD*DLAMCH( 'Eps' ), TOL(1) )
         TOL(2) = ZERO
         TOL(3) = ZERO
         TOL(4) = ZERO
      ELSE
         TOL(1) = SQRT( DLAMCH( 'Eps' ) )
         TOL(2) = ZERO
         TOL(3) = ZERO
         TOL(4) = ZERO
      END IF
C
C Determine the lenghts of working arrays.
C Use a larger value for LDWORK for enabling calls of block algorithms
C in DGGES, etc.
C
      NN  = N*N
      LDA = MAX( 1, N )
      LDB = LDA
      LDC = MAX( 1, P )
      LDD = LDC
      IF ( ISYSD.EQ.0 ) THEN
         LDE  = LDA
         ICOL = N
      ELSE IF ( ISYSD.LT.0 ) THEN
         LDE  = MAX( 1, RANKE )
         ICOL = RANKE
      ELSE
         LDE  = 1
         ICOL = 1
      END IF
C
      MAXPM = MAX( P, M )
      IF ( ISYSD.EQ.0 .AND. ISYS.EQ.1 ) THEN
         LBWORK = N
      ELSE
         LBWORK = 0
      END IF
C
      NBLK = N + MAXPM
      LIWORK = MAX( 2*NBLK+12, 4*NBLK + 3 )
C
C Allocate variable dimension local arrays.
C !Fortran 90/95
C
      ALLOCATE ( BWORK( LBWORK ), IWORK( LIWORK ) )
      A = mxCalloc( LDA*N,    8 )
      B = mxCalloc( LDB*M,    8 )
      C = mxCalloc( LDC*N,    8 )
      D = mxCalloc( LDD*M,    8 )
      E = mxCalloc( LDE*ICOL, 8 )
C
C Copy inputs from MATLAB workspace to locally allocated arrays.
C
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(1) ), %VAL( A ), NN )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(3) ), %VAL( B ), N*M )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(4) ), %VAL( C ), P*N )
      CALL mxCopyPtrToReal8( mxGetPr( PRHS(5) ), %VAL( D ), P*M )
C
      IF ( .NOT.UNITE ) THEN
         IF ( ISYSD.EQ.0 ) THEN
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL( E ),
     $                             NN )
C
C           Test whether E is the unit matrix.
C
            IF ( MA02HD( 'Full', N, N, ONE, %VAL( E ), LDE ) ) THEN
               JOBE  = 'I'
               ISYSD = 1
               CALL mxFree( E )
               LDE  = 1
               ICOL = 1
               E    = mxCalloc( LDE*ICOL, 8 )
            END IF
         ELSE IF ( ISYSD.LT.0 ) THEN
            CALL mxCopyPtrToReal8( mxGetPr( PRHS(2) ), %VAL( E ),
     $                             RANKE*RANKE )
         END IF
      END IF
C
C Test whether D is the zero matrix.
C
      IF ( MA02HD( 'Full', P, M, ZERO, %VAL( D ), LDD ) ) THEN
         JOBD = 'Z'
      ELSE
         JOBD = 'D'
      END IF
      WITHD = LSAME( JOBD, 'D' )
C
C Test whether B is the zero matrix.
C
      BNORM = DLANGE( '1-norm', N, M, %VAL( B ), LDB, %VAL( A ) )
C
C Test whether C is the zero matrix.
C
      CNORM = DLANGE( '1-norm', P, N, %VAL( C ), LDC, %VAL( A ) )
C
C Compute LDWORK and LZWORK and allocate DWORK and ZWORK.
C
      MINPM = MIN( P, M )
      NODYN = N.EQ.0 .OR. MIN( BNORM, CNORM ).EQ.ZERO
C
      IF( NODYN ) THEN
         IF( WITHD ) THEN
            LDWORK = P*M + MAX( 4*MINPM + MAXPM, 6*MINPM )
         ELSE
            LDWORK = 1
         END IF
         LZWORK = 1
      ELSE
         IF( UNITE .AND. ISYS.EQ.1 .AND. WITHD ) THEN
            IA = N*( M + P ) + MINPM
         ELSE
            IA = 0
         END IF
         IF( ISYSD.EQ.0 )
     $      IA = IA + NN
         LDWORK = IA + NN + N*( M + P ) + 11*NBLK*NBLK + 5*NBLK +
     $            MAX( 8*NBLK + 4, 36 )
         LZWORK = ( N + M )*( N + P ) + 2*MINPM + MAXPM
      END IF
C
      LDWORK = 2*LDWORK + 128*NBLK
      DWORK = mxCalloc( LDWORK,  8 )
      ZWORK = mxCalloc( LZWORK, 16 )
C
C Do the actual computations.
C
      CALL AB13HD( DICO, JOBE, EQUIL, JOBD, CKPROP, REDUCE, N, M, P,
     $             RANKE, FPEAK, %VAL( A ), LDA, %VAL( E ), LDE,
     $             %VAL( B ), LDB, %VAL( C ), LDC, %VAL( D ), LDD,
     $             GPEAK, TOL, IWORK, %VAL( DWORK ), LDWORK,
     $             %VAL( ZWORK ), LZWORK, BWORK, INFO )
C
C Copy output to MATLAB workspace.
C
      PLHS(1) = mxCreateDoubleMatrix( 2, 1, 0 )
      CALL mxCopyReal8ToPtr( GPEAK, mxGetPr( PLHS(1) ), 2 )
      IF( NLHS.GT.1 ) THEN
         PLHS(2) = mxCreateDoubleMatrix( 2, 1, 0 )
         CALL mxCopyReal8ToPtr( FPEAK, mxGetPr( PLHS(2) ), 2 )
      END IF
C
C Deallocate local arrays.
C !Fortran 90/95
C
      DEALLOCATE( BWORK, IWORK )
      CALL mxFree( A     )
      CALL mxFree( B     )
      CALL mxFree( C     )
      CALL mxFree( D     )
      CALL mxFree( E     )
      CALL mxFree( DWORK )
      CALL mxFree( ZWORK )
C
C Error and warning handling.
C
      IF ( INFO.NE.0 .AND. INFO.LT.7 ) THEN
         WRITE( TEXT, '('' INFO = '',I4,'' ON EXIT FROM AB13HD'')'
     $        ) INFO
         CALL mexErrMsgTxt( TEXT )
      ELSE IF ( INFO.EQ.7 ) THEN
         WRITE( TEXT, '('' The pencil is singular.'')' )
         CALL mexPrintf( TEXT )
      ELSE IF ( INFO.EQ.8 ) THEN
         WRITE( TEXT, '('' The system is improper.'')' )
         CALL mexPrintf( TEXT )
      END IF
C
      RETURN
C *** Last line of LINORM_H ***
      END
