      SUBROUTINE AB13HD( DICO, JOBE, EQUIL, JOBD, CKPROP, REDUCE, N, M,
     $                   P, RANKE, FPEAK, A, LDA, E, LDE, B, LDB, C,
     $                   LDC, D, LDD, GPEAK, TOL, IWORK, DWORK, LDWORK,
     $                   ZWORK, LZWORK, BWORK, INFO )
C
C     This is a preliminarily version of the function AB13HD to be
C     included in SLICOT.
C
C     This program is free software: you can redistribute it and/or
C     modify it under the terms of the GNU General Public License as
C     published by the Free Software Foundation, either version 3 of
C     the License, or (at your option) any later version.
C
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with this program.  If not, see
C     <http://www.gnu.org/licenses/>.
C
C     PURPOSE
C
C     To compute the L-infinity norm of a proper continuous-time or
C     causal discrete-time system, either standard or in the descriptor
C     form,
C
C                                     -1
C        G(lambda) = C*( lambda*E - A ) *B + D .
C
C     The norm is finite if and only if the matrix pair (A,E) has no
C     finite eigenvalue on the boundary of the stability domain, i.e.,
C     the imaginary axis, or the unit circle, respectively.
C
C     ARGUMENTS
C
C     Mode Parameters
C
C     DICO    CHARACTER*1
C             Specifies the type of the system, as follows:
C             = 'C':  continuous-time system;
C             = 'D':  discrete-time system.
C
C     JOBE    CHARACTER*1
C             Specifies whether E is a general square or an identity
C             matrix, as follows:
C             = 'I':  E is the identity matrix;
C             = 'G':  E is a general matrix;
C             = 'C':  E is in compressed form, i.e. E = [ T  0 ]
C                                                       [ 0  0 ]
C                     with a square full-rank matrix T.
C
C     EQUIL   CHARACTER*1
C             Specifies whether the user wishes to preliminarily
C             equilibrate the system (A,E,B,C) or (A,B,C), as follows:
C             = 'S':  perform equilibration (scaling);
C             = 'N':  do not perform equilibration.
C
C     JOBD    CHARACTER*1
C             Specifies whether or not a non-zero matrix D appears in
C             the given state space model:
C             = 'D':  D is present;
C             = 'Z':  D is assumed a zero matrix.
C
C     CKPROP  CHARACTER*1
C             If DICO = 'C' and JOBE <> 'I', specifies whether the user
C             wishes to check the properness of the transfer function of
C             the descriptor system, as follows: 
C             = 'C':  check the properness;
C             = 'N':  do not check the properness.
C             If the test is requested and the system is found improper
C             then GPEAK and FPEAK are both set to infinity, i.e., their
C             second component is zero; in addition, INFO is set to 8.
C             If the test is not requested, but the system is improper,
C             the resulted GPEAK and FPEAK may be wrong.
C             If DICO = 'D' or JOBE = 'I', this option is ineffective.
C
C     REDUCE  CHARACTER*1
C             If CKPROP = 'C', specifies whether the user wishes to
C             reduce the system order, by removing all uncontrollable
C             and unobservable poles before computing the norm, as
C             follows:
C             = 'R': reduce the system order;
C             = 'N': compute the norm without reducing the order.
C             If CKPROP = 'N', this option is ineffective.
C
C     Input/Output Parameters
C
C     N       (input) INTEGER
C             The order of the system.  N >= 0.
C
C     M       (input) INTEGER
C             The column size of the matrix B.  M >= 0.
C
C     P       (input) INTEGER
C             The row size of the matrix C.  P >= 0.
C
C     RANKE   (input) INTEGER
C             If JOBE = 'C', RANKE denotes the rank of the descriptor
C             matrix E or the size of the full-rank block T.
C             0 <= RANKE <= N.
C
C     FPEAK   (input/output) DOUBLE PRECISION array, dimension (2)
C             On entry, this parameter must contain an estimate of the
C             frequency where the gain of the frequency response would
C             achieve its peak value. Setting FPEAK(2) = 0 indicates an
C             infinite frequency. An accurate estimate could reduce the
C             number of iterations of the iterative algorithm. If no
C             estimate is available, set FPEAK(1) = 0, and FPEAK(2) = 1.
C             FPEAK(1) >= 0, FPEAK(2) >= 0.
C             On exit, if INFO = 0, this array contains the frequency
C             OMEGA, where the gain of the frequency response achieves
C             its peak value GPEAK, i.e.,
C
C                 || G ( j*OMEGA ) || = GPEAK ,  if DICO = 'C', or
C
C                         j*OMEGA
C                 || G ( e       ) || = GPEAK ,  if DICO = 'D',
C
C             where OMEGA = FPEAK(1), if FPEAK(2) > 0, and OMEGA is
C             infinite, if FPEAK(2) = 0.
C
C     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
C             The leading N-by-N part of this array must contain the
C             state dynamics matrix A.
C
C     LDA     INTEGER
C             The leading dimension of the array A.  LDA >= max(1,N).
C
C     E       (input) DOUBLE PRECISION array, dimension (LDE,K), where
C             K is N, RANKE, or 0, if JOBE = 'G', 'C', or 'I',
C             respectively.
C             If JOBE = 'G', the leading N-by-N part of this array must
C             contain the descriptor matrix E of the system.
C             If JOBE = 'C', the leading RANKE-by-RANKE part of this
C             array must contain the full-rank block T of the descriptor
C             matrix E.
C             If JOBE = 'I', then E is assumed to be the identity
C             matrix and is not referenced.
C
C     LDE     INTEGER
C             The leading dimension of the array E.
C             LDE >= MAX(1,N),     if JOBE = 'G';
C             LDE >= MAX(1,RANKE), if JOBE = 'C';
C             LDE >= 1,            if JOBE = 'I'.
C
C     B       (input) DOUBLE PRECISION array, dimension (LDB,M)
C             The leading N-by-M part of this array must contain the
C             system input matrix B.
C
C     LDB     INTEGER
C             The leading dimension of the array B.  LDB >= max(1,N).
C
C     C       (input) DOUBLE PRECISION array, dimension (LDC,N)
C             The leading P-by-N part of this array must contain the
C             system output matrix C.
C
C     LDC     INTEGER
C             The leading dimension of the array C.  LDC >= max(1,P).
C
C     D       (input) DOUBLE PRECISION array, dimension (LDD,M)
C             If JOBD = 'D', the leading P-by-M part of this array must
C             contain the direct transmission matrix D.
C             The array D is not referenced if JOBD = 'Z'.
C
C     LDD     INTEGER
C             The leading dimension of array D.
C             LDD >= MAX(1,P), if JOBD = 'D';
C             LDD >= 1,        if JOBD = 'Z'.
C
C     GPEAK   (output) DOUBLE PRECISION array, dimension (2)
C             The L-infinity norm of the system, i.e., the peak gain
C             of the frequency response (as measured by the largest
C             singular value in the MIMO case), coded in the same way
C             as FPEAK.
C
C     Tolerances
C
C     TOL     DOUBLE PRECISION array, dimension K, where K = 2, if
C             CKPROP = 'N' and ( DICO = 'D' or JOBE = 'I' ), and K = 4,
C             otherwise.
C             TOL(1) is the tolerance used to set the accuracy in
C             determining the norm.  0 <= TOL(1) < 1.
C             TOL(2) is the threshold value for magnitude of the matrix
C             elements, if EQUIL = 'S': elements with magnitude less
C             than or equal to TOL(2) are ignored for scaling. If the
C             user sets TOL(2) >= 0, then the given value of TOL(2) is
C             used. If the user sets TOL(2) < 0, then an implicitly
C             computed, default threshold, defined by THRESH = c*EPS, is
C             used instead, where c = MAX(norm_1(A,E,B,C)) and EPS is
C             the machine precision (see LAPACK Library routine DLAMCH).
C             TOL(2) = 0 is not always a good choice.  TOL(2) < 1.
C             TOL(2) is not used if EQUIL = 'N'.
C             TOL(3) is the tolerance to be used in rank determinations
C             when transforming (lambda*E-A,B,C), if CKPROP = 'C'. If
C             the user sets TOL(3) > 0, then the given value of TOL(3)
C             is used as a lower bound for reciprocal condition numbers
C             in rank determinations; a (sub)matrix whose estimated
C             condition number is less than 1/TOL(3) is considered to be
C             of full rank.  If the user sets TOL(3) <= 0, then an
C             implicitly computed, default tolerance, defined by
C             TOLDEF1 = N*N*EPS, is used instead.  TOL(3) < 1.
C             TOL(4) is the tolerance to be used for checking pencil
C             singularity when CKPROP = 'N', or singularity of the
C             matrices A and E when CKPROP = 'C'. If the user sets
C             TOL(4) > 0, then the given value of TOL(4) is used.
C             If the user sets TOL(4) <= 0, then an implicitly
C             computed, default tolerance, defined by  TOLDEF2 = 10*EPS,
C             is used instead.  TOL(4) < 1.
C
C     Workspace
C
C     IWORK   INTEGER array, dimension (LIWORK)
C             LIWORK >= MAX(N+MAX(M,P)+12,2*N+MAX(M,P)+7),
C                          if CKPROP = 'C', DICO = 'C' and JOBE = 'G';
C             LIWORK >= MAX(N+MAX(M,P)+12,2*N), otherwise.
C
C     DWORK   DOUBLE PRECISION array, dimension (LDWORK)
C             On exit, if INFO = 0, DWORK(1) contains the optimal value
C             of LDWORK.
C
C     LDWORK  INTEGER
C             The dimension of the array DWORK.
C             LDWORK >= 1, if N = 0, or B = 0, or C = 0; and JOBD = 'Z';
C             LDWORK >= P*M + MAX( 4*MIN(P,M) + MAX(P,M), 6*MIN(P,M) ),
C                       if  N = 0, or B = 0, or C = 0; and JOBD = 'D';
C             LDWORK >= a + c*N*N + N*(M+P) + 11*NBLK*NBLK + 5*NBLK +
C                       MAX(L,36) ) where NBLK = N+MAX(P,M), and 
C                       L = 8*NBLK + 4, if NBLK is even, and L = 8*NBLK,
C                       if NBLK is odd otherwise, where
C                       a = N*(M+P) + MIN(P,M), if DICO = 'C' and 
C                                   JOBE = 'I' and JOBD = 'D';
C                       a = 0, otherwise;
C                       c = 1, if JOBE = 'I', and c = 2, otherwise.
C                         + MAX(N*N+4*N+4,2*(MAX(M,P)+N-1)),
C                                         if CKPROP = 'C', REDUCE = 'R',
C                         + MAX(4*N+4,2*(MAX(M,P)+N-1),N*N+4*N)+2*N*N,
C                                         if CKPROP = 'C', REDUCE = 'N'.
C                         (from IWRK at call AB13ID).
C
C             For good performance, LDWORK must generally be larger.
C
C     ZWORK   COMPLEX*16 array, dimension (LZWORK)
C             On exit, if INFO = 0, ZWORK(1) contains the optimal
C             LZWORK.
C
C     LZWORK  INTEGER
C             The dimension of the array ZWORK.
C             LZWORK >= 1,  if N = 0, or B = 0, or C = 0;
C             LZWORK >= MAX(1, (N+M)*(N+P) + 2*MIN(P,M) + MAX(P,M)),
C                           otherwise.
C             For good performance, LZWORK must generally be larger.
C
C     BWORK   LOGICAL array, dimension (N)
C             Not referenced, if JOBE <> 'G' or DICO = 'D'.
C
C     Error Indicator
C
C     INFO    INTEGER
C             = 0:  successful exit;
C             < 0:  if INFO = -i, the i-th argument had an illegal
C                   value;
C             = 1:  a matrix used for computing the largest singular
C                   value of G(infinity) is (numerically) singular;
C             = 2:  the (periodic) QR (or QZ) algorithm for computing
C                   eigenvalues did not converge;
C             = 3:  the SVD algorithm for computing singular values did
C                   not converge;
C             = 4:  the tolerance is too small and the algorithm did
C                   not converge;
C             = 5:  other computations than QZ iteration, or reordering
C                   of eigenvalues failed in the LAPACK Library routine
C                   DHGEQZ;
C             = 6:  the numbers of "finite" eigenvalues before and after
C                   reordering differ; the threshold used might be
C                   unsuitable.
C             = 7:  the descriptor system is singular. GPEAK(1) and
C                   GPEAK(2) are set to 1 and 0, respectively,
C                   corresponding to infinity. FPEAK(1) and FPEAK(2) are
C                   set similarly. This is a warning, not error.
C             = 8:  the descriptor system is improper. GPEAK(1) and
C                   GPEAK(2) are set to 1 and 0, respectively,
C                   corresponding to infinity. FPEAK(1) and FPEAK(2) are
C                   set similarly. This is a warning, not error. It can
C                   appear only if CKPROP = 'C'.
C
C     METHOD
C
C     The routine implements the method presented in [2] which is an
C     extension of the method in [1] for descriptor systems. There are
C     several improvements and refinements to increase numerical
C     robustness, accuracy and efficiency such as the usage of
C     structure-preserving eigenvalue computations for skew-Hamiltonian/
C     Hamiltonian eigenvalue problems in the iterative method in [2].
C
C     REFERENCES
C
C     [1] Bruinsma, N.A. and Steinbuch, M.
C         A fast algorithm to compute the H-infinity-norm of a transfer
C         function matrix.
C         Systems & Control Letters, vol. 14, pp. 287-293, 1990.
C
C     [2] Voigt, M.
C         L-infinity-Norm Computation for Descriptor Systems.
C         Diploma Thesis, Fakultaet fuer Mathematik, TU Chemnitz,
C         http://nbn-resolving.de/urn:nbn:de:bsz:ch1-201001050.
C
C     NUMERICAL ASPECTS
C
C     If the algorithm does not converge in MAXIT = 30 iterations
C     (INFO = 4), the tolerance must be increased, or the system is im-
C     proper.
C
C     CONTRIBUTORS
C
C     M. Voigt, Max Planck Institute for Dynamics of Complex Technical 
C     Systems, March 2011.
C     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2011.
C     Partly based on the SLICOT Library routine AB13DD by D. Sima,
C     V. Sima, D.W. Gu, and M.M. Konstantinov.
C
C     REVISIONS
C
C     V. Sima, Mar. 2012, Apr. 2012, May 2012, June 2012.
C     M. Voigt, Apr. 2017, Sep. 2017.
C
C     KEYWORDS
C
C     H-infinity optimal control, robust control, system norm.
C
C     ******************************************************************
C
C     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
      DOUBLE PRECISION   ZERO, ONE, TWO, FOUR, P25, HND2
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0,  TWO  = 2.0D+0,
     $                     FOUR = 4.0D+0, P25 = 0.25D+0, HND2 = 1.0D2 )
      DOUBLE PRECISION   TEN, HUNDRD, THOUSD
      PARAMETER          ( TEN    = 1.0D+1, HUNDRD = 1.0D+2,
     $                     THOUSD = 1.0D+3 )
C     ..
C     .. Scalar Arguments ..
      CHARACTER          CKPROP, DICO, EQUIL, JOBD, JOBE, REDUCE
      INTEGER            INFO, LDA, LDB, LDC, LDD, LDE, LDWORK, LZWORK,
     $                   M, N, P, RANKE
      DOUBLE PRECISION   TOL( 4 )
C     ..
C     .. Array Arguments ..
      COMPLEX*16         ZWORK(  * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   D( LDD, * ), DWORK(  * ), E( LDE, * ),
     $                   FPEAK(  2 ), GPEAK(  2 )
      INTEGER            IWORK(  * )
      LOGICAL            BWORK(  * )
C     ..
C     .. Local Scalars ..
      CHARACTER          JOB, JOBEIG, JOBSYS, RESTOR, UPDATE, VECT
      LOGICAL            CMPRE, DISCR, FULLE, GENE, ILASCL, ILESCL,
     $                   IND1, ISPROP, LEQUIL, NCMPRE, NODYN, NSRT,
     $                   UNITE, USEPEN, WCKPRP, WITHD, WREDUC
      LOGICAL            LINF
      INTEGER            I, IA, IAR, IAS, IAW, IB, IBI, IBS, IBT, IBV,
     $                   IBW, IC, ICI, ICS, ICU, ICW, ID, IDI, IE, IERR,
     $                   IES, IEW, IH, IH12, IH22, IHC, IHI, II, IJ,
     $                   IJ12, ILFT, ILO, IM, IMIN, IN, IQ, IR, IRHT,
     $                   IS, ISB, ISC, ISL, IT, IT12, ITAU, ITER, IU,
     $                   IV, IWARN, IWRK, IZ, J, K, LIW, LW, M0, MAXCWK,
     $                   MAXPM, MAXWRK, MINCWR, MINPM, MINWRK, NBLK,
     $                   NCC, NE, NEI, NINF, NN, NR, NWS, P0, PM, Q, R,
     $                   RNKE, SDIM, SDIM1
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNORM, BOUND, CNORM, CND,
     $                   DIF, ENRM, ENRMTO, EPS, FPEAKI, FPEAKS, GAMMA,
     $                   GAMMAL, GAMMAS, MAXRED, OMEGA, PI, RAT, RCOND,
     $                   RTOL, SAFMAX, SAFMIN, SCL, SMLNUM, SV1, SVP,
     $                   THRESH, TM, TMP, TOLDEF, TOLE, TOLER, TOLN,
     $                   TZER, WMAX, WRMIN
      COMPLEX*16         CTMP1, CTMP2
C     ..
C     .. Local Arrays ..
      INTEGER            IDUM( 1 )
      DOUBLE PRECISION   TEMP( 1 ), TOLI( 2 )
C     ..
C     .. External Functions ..
      LOGICAL            AB13ID, LSAME
      DOUBLE PRECISION   AB13DX, DLAMCH, DLANGE, DLANHS, DLANTR, DLAPY2
      EXTERNAL           AB13DX, DLAMCH, DLANGE, DLANHS, DLANTR, DLAPY2,
     $                   LSAME
C     ..
C     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEBAL, DGECON, DGEHRD, DGEMM, DGESV,
     $                   DGESVD, DGGBAL, DGGBAK, DHGEQZ, DHSEQR, DLABAD,
     $                   DLACPY, DLASCL, DLASET, DLASRT, DORMHR, DSCAL,
     $                   DSWAP, DSYRK, DTGSEN, DTRCON, DTRSM, MA02AD,
     $                   MB01SD, MB02RD, MB02SD, MB02TD, MB03XD, MB04BD,
     $                   SB04OD, TB01ID, TG01AD, TG01BD, XERBLA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, ATAN, ATAN2, COS, DBLE, DCMPLX, INT, LOG,
     $                   MAX, MIN, SIN, SQRT
C     ..
C     .. Executable Statements ..
C
C     Test the input scalar parameters.
C
      NN     = N*N
      MINPM  = MIN( P, M )
      MAXPM  = MAX( P, M )
      INFO   = 0
      DISCR  = LSAME( DICO,   'D' )
      UNITE  = LSAME( JOBE,   'I' )
      GENE   = LSAME( JOBE,   'G' )
      CMPRE  = LSAME( JOBE,   'C' )
      LEQUIL = LSAME( EQUIL,  'S' )
      WITHD  = LSAME( JOBD,   'D' )
      WCKPRP = LSAME( CKPROP, 'C' )
      WREDUC = LSAME( REDUCE, 'R' )
      FULLE  = GENE .OR. CMPRE
C
      IF( .NOT. ( DISCR .OR. LSAME( DICO, 'C' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT. ( FULLE  .OR. UNITE ) ) THEN
         INFO = -2
      ELSE IF( .NOT. ( LEQUIL .OR. LSAME( EQUIL,  'N' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT. ( WITHD  .OR. LSAME( JOBD,   'Z' ) ) ) THEN
         INFO = -4
      ELSE IF( .NOT. ( WCKPRP .OR. LSAME( CKPROP, 'N' ) ) ) THEN
         IF( .NOT.( DISCR .OR. UNITE ) )
     $      INFO = -5
      ELSE IF( .NOT. ( WREDUC .OR. LSAME( REDUCE, 'N' ) ) ) THEN
         IF( WCKPRP )
     $      INFO = -6
      ELSE IF( N.LT.0 ) THEN
         INFO = -7
      ELSE IF( M.LT.0 ) THEN
         INFO = -8
      ELSE IF( P.LT.0 ) THEN
         INFO = -9
      ELSE IF( CMPRE .AND. ( RANKE.LT.0 .OR. RANKE.GT.N ) ) THEN
         INFO = -10
      ELSE IF( MIN( FPEAK( 1 ), FPEAK( 2 ) ).LT.ZERO ) THEN
         INFO = -11
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -13
      ELSE IF( LDE.LT.1 .OR. ( GENE  .AND. LDE.LT.N ) .OR.
     $                       ( CMPRE .AND. LDE.LT.RANKE ) ) THEN
         INFO = -15
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -17
      ELSE IF( LDC.LT.MAX( 1, P ) ) THEN
         INFO = -19
      ELSE IF( LDD.LT.1 .OR. ( WITHD .AND. LDD.LT.P ) ) THEN
         INFO = -21
      ELSE IF( TOL( 1 ).LT.ZERO .OR. TOL( 1 ).GE.ONE ) THEN
         INFO = -23
      ELSE IF( TOL( 2 ).LT.ZERO .OR. TOL( 2 ).GE.ONE ) THEN
         INFO = -23
      ELSE IF( ( WCKPRP .OR. ( .NOT.DISCR .AND. FULLE ) ) .AND.
     $         ( TOL( 3 ).LT.ZERO .OR. TOL( 3 ).GE.ONE ) ) THEN
         INFO = -23
      ELSE IF( ( .NOT.WCKPRP .AND. ( DISCR .OR. UNITE ) ) .AND.
     $         ( TOL( 4 ).LT.ZERO .OR. TOL( 4 ).GE.ONE ) ) THEN
         INFO = -23
      ELSE
         BNORM  = DLANGE( '1-norm', N, M, B, LDB, DWORK )
         CNORM  = DLANGE( '1-norm', P, N, C, LDC, DWORK )
         NODYN  = N.EQ.0 .OR. MIN( BNORM, CNORM ).EQ.ZERO
         USEPEN = FULLE .OR. DISCR
C
         IF ( CMPRE ) THEN
            NCMPRE = RANKE.EQ.N
            CMPRE  = RANKE.LT.N
         ELSE
            NCMPRE = .FALSE.
         END IF
C
         IF( DISCR .OR. UNITE )
     $      WCKPRP = .FALSE.
C
C        Compute workspace.
C
C        Note: The workspace requirement is dominated by the solution of
C              the generalized skew-Hamiltonian/Hamiltonian eigenvalue
C              problem.
C
         IF( NODYN ) THEN
            IF( WITHD ) THEN
               MINWRK = P*M + MAX( 4*MINPM + MAXPM, 6*MINPM )
            ELSE
               MINWRK = 1
            END IF
         ELSE
            R    = MOD( M+P, 2 )
            NBLK = N + ( M + P + R )/2
            IF( UNITE .AND. .NOT.DISCR .AND. WITHD ) THEN
               IA = N*( M + P ) + MINPM
            ELSE
               IA = 0
            END IF
            IF( FULLE )
     $         IA = IA + NN
            MINWRK = IA + NN + N*( M + P ) + 11*NBLK*NBLK + 5*NBLK +
     $               MAX( 2*NBLK, 32 )
         END IF
         IF( LDWORK.LT.MINWRK ) THEN
            INFO = -26
         ELSE
            IF ( NODYN ) THEN
               MINCWR = 1
            ELSE
               MINCWR = MAX( 1, ( N + M )*( N + P ) + 2*MINPM + MAXPM )
            END IF
            IF( LZWORK.LT.MINCWR )
     $         INFO = -28
         END IF
      END IF 
C
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'AB13HD', -INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( MINPM.EQ.0 ) THEN
         GPEAK( 1 ) = ZERO
         FPEAK( 1 ) = ZERO
         GPEAK( 2 ) = ONE
         FPEAK( 2 ) = ONE
         DWORK( 1 ) = ONE
         ZWORK( 1 ) = ONE
         RETURN
      END IF
C
C     Determine the maximum singular value of G(infinity) = D,
C     if B = 0 or C = 0; or DICO = 'C', and (JOBE = 'I' or JOBE = 'C'
C     and RANKE = N) and JOBD = 'D'.
C
C     If JOBE = 'I' and DICO = 'C', the full SVD of D, D = U*S*V', is
C     computed and saved for later use.
C
C     (Note: Comments in the code beginning "Workspace:" describe the
C     minimal amount of real workspace needed at that point in the
C     code, as well as the preferred amount for good performance.
C     NB refers to the optimal block size for the immediately
C     following subroutine, as returned by ILAENV.)
C
      IF ( WITHD ) THEN
         IF ( NODYN .OR. ( ( UNITE .OR. NCMPRE ) .AND. .NOT.DISCR ) )
     $         THEN
            IF ( .NOT.NODYN .AND. UNITE ) THEN
               IBV  = 1
               ICU  = IBV + N*M
               IS   = ICU + P*N
               IU   = IS  + MINPM
               IV   = IU  + P*P
               ID   = IV  + M*M
               VECT = 'A'
            ELSE
               IU   = 1
               IV   = 1
               IS   = 1
               ID   = 1 + MINPM
               VECT = 'N'
            END IF
            IWRK = ID + P*M
C
C           Workspace: need   P*M + MIN(P,M) + V +
C                             MAX( 3*MIN(P,M) + MAX(P,M), 5*MIN(P,M) ),
C                             where V = N*(M+P) + P*P + M*M,
C                                       if JOBE = 'I', DICO = 'C',
C                                       JOBD = 'D', and B <> 0, C <> 0,
C                                   V = 0,  otherwise;
C                      prefer larger.
C
            CALL DLACPY( 'Full', P, M, D, LDD, DWORK( ID ), P )
            CALL DGESVD( VECT, VECT, P, M, DWORK( ID ), P, DWORK( IS ),
     $                   DWORK( IU ), P, DWORK( IV ), M, DWORK( IWRK ),
     $                   LDWORK-IWRK+1, IERR )
            IF( IERR.GT.0 ) THEN
               INFO = 3
               RETURN
            END IF
            GAMMAL = DWORK( IS )
            MAXWRK = INT( DWORK( IWRK ) ) + IWRK - 1
C
C           Save needed singular values.
C
            SV1 = GAMMAL
            SVP = DWORK( IS+MINPM-1 )
         ELSE
            GAMMAL = ZERO
            MAXWRK = 1
         END IF
      ELSE
         GAMMAL = ZERO
         MAXWRK = 1
      END IF
C
C     Quick return if possible.
C
      IF( NODYN ) THEN
         GPEAK( 1 ) = GAMMAL
         FPEAK( 1 ) = ZERO
         GPEAK( 2 ) = ONE
         FPEAK( 2 ) = ONE
         DWORK( 1 ) = MAXWRK
         ZWORK( 1 ) = ONE
         RETURN
      END IF
C
      IF ( .NOT.USEPEN .AND. WITHD ) THEN
C
C        Standard continuous-time case, D <> 0: Compute B*V and C'*U .
C
         CALL DGEMM( 'No Transpose', 'Transpose', N, M, M, ONE, B, LDB,
     $               DWORK( IV ), M, ZERO, DWORK( IBV ), N )
         CALL DGEMM( 'Transpose', 'No Transpose', N, P, P, ONE, C,
     $               LDC, DWORK( IU ), P, ZERO, DWORK( ICU ), N )
C
C        U and V are no longer needed: free their memory space.
C        Total workspace here: need   N*(M+P) + MIN(P,M)
C        (JOBE = 'I', DICO = 'C', JOBD = 'D').
C
         IA = IU
      ELSE
         IA = 1
      END IF
C
C     Get machine constants.
C
      EPS    = DLAMCH( 'Precision' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFMAX = ONE / SAFMIN
      CALL DLABAD( SAFMIN, SAFMAX )
      SMLNUM = SQRT( SAFMIN ) / EPS
      BIGNUM = ONE / SMLNUM
      TOLER  = SQRT( EPS )
C
C     Initiate the transformation of the system to an equivalent one,
C     to be used for eigenvalue computations.
C
C     Additional workspace: need   N*N + N*d + 2*N, if JOBE = 'I';
C                                2*N*N + N*d + 2*N, otherwise,
C                           where  d = M + P,       if CKPROP = 'N';
C                                  d = 2*MAX(M,P),  otherwise.
C
      IE = IA + NN
      IF ( FULLE ) THEN
         IB = IE + NN
      ELSE
         IB = IE
      END IF
      IF ( WCKPRP ) THEN
         IC = IB + N*MAXPM
         IR = IC + MAXPM*N
      ELSE
         IC = IB + N*M
         IR = IC + P*N
      END IF
      II  = IR + N
      IBT = II + N
C
      CALL DLACPY( 'Full', N, N, A, LDA, DWORK( IA ), N )
      CALL DLACPY( 'Full', N, M, B, LDB, DWORK( IB ), N )
      CALL DLACPY( 'Full', P, N, C, LDC, DWORK( IC ), P )
C
C     Scale A if maximum element is outside the range [SMLNUM,BIGNUM].
C
      ANRM   = DLANGE( 'Max', N, N, DWORK( IA ), N, DWORK )
      ILASCL = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      END IF
      IF( ILASCL )
     $   CALL DLASCL( 'General', 0, 0, ANRM, ANRMTO, N, N, DWORK( IA ),
     $                N, IERR )
C
      BOUND = EPS*THOUSD
      TOLN = TOL( 1 )
C
      NR = N
      IF ( FULLE ) THEN
C
C        Descriptor system.
C
         IWRK = IBT

         IF( CMPRE ) THEN
C
C           Set up the full matrix E, if JOBE = 'C'.
C
            K = RANKE
            CALL DLACPY( 'Full', K, K, E, LDE, DWORK( IE ), N )
            CALL DLASET( 'Full', N-K, K, ZERO, ZERO, DWORK( IE+K ), N )
            CALL DLASET( 'Full', N, N-K, ZERO, ZERO, DWORK( IE+N*K ),
     $                   N )
         ELSE
            K = N
            CALL DLACPY( 'Full', N, N, E, LDE, DWORK( IE ), N )
         END IF
C
C        Scale E if maximum element is outside the range
C        [SMLNUM,BIGNUM].
C
         ENRM   = DLANGE( 'Max', K, K, DWORK( IE ), N, DWORK )
         ILESCL = .FALSE.
         IF( ENRM.GT.ZERO .AND. ENRM.LT.SMLNUM ) THEN
            ENRMTO = SMLNUM
            ILESCL = .TRUE.
         ELSE IF( ENRM.GT.BIGNUM ) THEN
            ENRMTO = BIGNUM
            ILESCL = .TRUE.
         END IF
         IF( ILESCL )
     $      CALL DLASCL( 'General', 0, 0, ENRM, ENRMTO, K, K,
     $                   DWORK( IE ), N, IERR )
C
C        Set the tolerances.
C
         IF( LEQUIL ) THEN
            THRESH = TOL( 2 )
            IF( THRESH.LT.ZERO )
     $          THRESH = MAX( DLANGE( '1-norm', N, N, A, LDA, DWORK ),
     $                        DLANGE( '1-norm', N, N, E, LDE, DWORK ),
     $                        BNORM, CNORM )*EPS
         END IF
         TOLDEF = TOL( 3 )
         TZER   = TOL( 4 )
         IF( TOLDEF.LE.ZERO )
     $       TOLDEF = N*N*EPS
         IF( TZER.LE.ZERO )
     $       TZER = N*EPS
         TOLI( 1 ) = TOLDEF
         TOLI( 2 ) = TZER
C
C        Equilibrate the system, if required.
C
C        Additional workspace: need   6*N.
C
         IF( LEQUIL )
     $      CALL TG01AD( 'All', N, N, M, P, THRESH, DWORK( IA ), N,
     $                   DWORK( IE ), N, DWORK( IB ), N, DWORK( IC ), P,
     $                   DWORK( II ), DWORK( IR ), DWORK( IWRK ),
     $                   IERR )
         IF( WCKPRP ) THEN
C
C           Check properness on the transformed system.
C
C           Addit. workspace: need   MAX(w+4*N+4,x)+z,
C                             where  w = N*N,           if REDUCE = 'R',
C                                    w = 0,             if REDUCE = 'N',
C                                    x = MAX(2*(MAX(M,P)+N-1),N*N+4*N),
C                                                       if RESTOR = 'N',
C                                    z = 0,             if REDUCE = 'R',
C                                    z = 2*N*N,         if REDUCE = 'N';
C                             equivalently,
C                                    MAX(N*N+4*N+4,2*(MAX(M,P)+N-1)),
C                                                       if REDUCE = 'R',
C                                    MAX(4*N+4,2*(MAX(M,P)+N-1),N*N+4*N)
C                                             +2*N*N,   otherwise;
C                             prefer larger.
C
C           Integer workspace: need  2*N+MAX(M,P)+7, if JOBE = 'G';
C                                    N,              if JOBE = 'C'.
C
            RESTOR = 'N'
C
C           RESTOR = 'R', for saving and restoring.
C
            IF( CMPRE .OR. NCMPRE ) THEN
               JOBSYS = 'N'
            ELSE
               JOBSYS = 'R'
            END IF
C
            IF( WREDUC ) THEN
               JOBEIG = 'A'
               UPDATE = 'U'
               IAW    = IA
               IEW    = IE
               IBW    = IB
               ICW    = IC
               IF( P.LT.M ) THEN
                  DO 6 J = N, 1, -1
                     CALL DCOPY( P, DWORK( IC+(J-1)*P     ), -1,
     $                              DWORK( IC+(J-1)*MAXPM ), -1 )
    6             CONTINUE
               END IF
            ELSE
               JOBEIG = 'I'
               UPDATE = 'N'
               IAW    = IWRK
               IEW    = IAW + NN
               IBW    = IEW + NN
               ICW    = IBW + N*MAXPM
               IWRK   = ICW + N*MAXPM
               CALL DLACPY( 'Full', N, N, DWORK( IA ), N, DWORK( IAW ),
     $                      N )
               CALL DLACPY( 'Full', N, N, DWORK( IE ), N, DWORK( IEW ),
     $                      N )
               CALL DLACPY( 'Full', N, M, DWORK( IB ), N, DWORK( IBW ),
     $                      N )
               CALL DLACPY( 'Full', P, N, DWORK( IC ), P, DWORK( ICW ),
     $                      MAXPM )
            END IF
C
            ISPROP = AB13ID( JOBSYS, JOBEIG, 'No equil', 'Not ck sing',
     $                       RESTOR, UPDATE, N, M, P, DWORK( IAW ), N,
     $                       DWORK( IEW ), N, DWORK( IBW ), N,
     $                       DWORK( ICW ), MAXPM, NR, RNKE, TOLI, IWORK,
     $                       DWORK( IWRK ), LDWORK-IWRK+1, IWARN, INFO )
            IF( .NOT.ISPROP ) THEN
               INFO = 8
               GPEAK( 1 ) = ONE
               GPEAK( 2 ) = ZERO
               FPEAK( 1 ) = ONE
               FPEAK( 2 ) = ZERO
               RETURN
            END IF
            IF( P.LT.M .AND. WREDUC )
     $         CALL DLACPY( 'Full', P, NR, DWORK( IC ), MAXPM,
     $                      DWORK( IC ), P )
            IF( .NOT.WREDUC ) THEN
               IWRK = IAW
               NR   = N
            END IF
         END IF
C
         IF( CMPRE .AND. .NOT.DISCR ) THEN
C
C           Determine the largest singular value of G(infinity).
C
C           Additional workspace if DICO = 'C' and RANKE < N: need  P*M.
C
            ID   = IWRK
            IS   = ID + P*M
            NINF = N - RANKE
            IF( NINF.GT.0 ) THEN
C
C              Compute G(infinity) =
C                 D-C(1:P,RANKE+1:N)*inv(A(RANKE+1:N,RANKE+1:N))*
C                   B(RANKE+1,1:M).
C
C              Additional workspace: need   NINF*M + NINF*NINF + 4*NINF,
C                                           NINF = N-RANKE.
C              Integer    workspace: need   NINF.
C
               IBS  = IS
               IAS  = IBS + NINF*M
               IWRK = IAS + NINF*NINF
               CALL DLACPY( 'Full', NINF, M, DWORK( IB+RANKE ), N,
     $                      DWORK( IBS ), NINF )
               CALL DLACPY( 'Full', NINF, NINF, A( RANKE+1, RANKE+1 ),
     $                      LDA, DWORK( IAS ), NINF )
               TMP = DLANGE( '1-norm', NINF, NINF, DWORK( IAS ), NINF,
     $                       DWORK )
               CALL DGESV(  NINF, M, DWORK( IAS ), NINF, IWORK,
     $                      DWORK( IBS ), NINF, IERR )
               CALL DGECON( '1-norm', NINF, DWORK( IAS ), NINF, TMP,
     $                      RCOND, DWORK( IWRK ), IWORK, IERR )
               IF( RCOND.LE.TEN*DBLE( NINF )*EPS ) THEN
C
C                 The matrix A(RANKE+1:N,RANKE+1:N) is numerically
C                 singular, no safe computation of G(infinity) is
C                 possible.
C
                  INFO = 1
                  RETURN
               END IF
C
               IF( WITHD ) THEN
C
C                 Copy D to DWORK.
C
                  CALL DLACPY( 'Full', P, M, D, LDD, DWORK( ID ), P )
                  CALL DGEMM( 'No Transpose', 'No Transpose', P, M,
     $                        NINF, -ONE, DWORK( IC+P*NINF ), P,
     $                        DWORK( IBS ), NINF, ONE, DWORK( ID ), P )
               ELSE
                  CALL DGEMM( 'No Transpose', 'No Transpose', P, M,
     $                        NINF, -ONE, DWORK( IC+P*NINF ), P,
     $                        DWORK( IBS ), NINF, ZERO, DWORK( ID ), P )
               END IF
            ELSE IF( WITHD ) THEN
C
C              Copy D to DWORK.
C
               CALL DLACPY( 'Full', P, M, D, LDD, DWORK( ID ), P )
            ELSE
               IS = ID
            END IF
C
C           Compute the maximum singular value of G(infinity).
C
C           Additional workspace: need   MAX( 4*MIN(P,M) + MAX(P,M),
C                                             6*MIN(P,M) );
C                                 prefer larger.
C
            IWRK = IS + MINPM
            CALL DGESVD( 'No Computation', 'No Computation', P, M,
     $                   DWORK( ID ), P, DWORK( IS ), DWORK, 1,  DWORK,
     $                   1, DWORK( IWRK ), LDWORK-IWRK+1, IERR )
            IF( IERR.GT.0 ) THEN
C
C              The SVD algorithm did not converge, no computation of
C              the largest singular value of G(infinity) is possible.
C
               INFO = 3
               RETURN
            END IF
            MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
            GAMMAL = DWORK( IS )
            IWRK   = IBT
         END IF
C
         IES = IBT + N
         IAS = IES + NN
         IQ  = IES
         IZ  = IAS
C
C        For efficiency of later calculations, the system (A,E,B,C)
C        is reduced to an equivalent one with the state matrix A in
C        Hessenberg form, and E upper triangular.
C        First, permute (A,E) to make it more nearly triangular.
C
         NSRT = CMPRE .OR. DISCR
         IF( NSRT ) THEN
            ILFT = II
            IRHT = IR
         ELSE
            ILFT = IZ   + NN
            IRHT = ILFT + N
         END IF
C
         CALL DGGBAL( 'Permute', NR, DWORK( IA ), N, DWORK( IE ), N,
     $                ILO, IHI, DWORK( ILFT ), DWORK( IRHT ),
     $                DWORK( IWRK ), IERR )
         IF( NSRT ) THEN
            VECT = 'N'
            M0   = M
            P0   = P
C
C           Apply the permutations to (the copies of) B and C.
C
            DO 10 I = NR, IHI + 1, -1
               K = DWORK( II+I-1 )
               IF( K.NE.I )
     $            CALL DSWAP( M, DWORK( IB+I-1 ), N,
     $                           DWORK( IB+K-1 ), N )
               K = DWORK( IR+I-1 )
               IF( K.NE.I )
     $            CALL DSWAP( P, DWORK( IC+(I-1)*P ), 1,
     $                           DWORK( IC+(K-1)*P ), 1 )
   10       CONTINUE
C
            DO 20 I = 1, ILO - 1
               K = DWORK( II+I-1 )
               IF( K.NE.I )
     $            CALL DSWAP( M, DWORK( IB+I-1 ), N,
     $                           DWORK( IB+K-1 ), N )
               K = DWORK( IR+I-1 )
               IF( K.NE.I )
     $            CALL DSWAP( P, DWORK( IC+(I-1)*P ), 1,
     $                           DWORK( IC+(K-1)*P ), 1 )
   20       CONTINUE
C
         ELSE
            VECT = 'I'
            M0   = 0
            P0   = 0
            IWRK = IRHT + N
         END IF
C
C        Reduce (A,E) to generalized Hessenberg form. Apply the
C        transformations to B and C, if NSRT = .TRUE., i.e., (JOB = 'C'
C        and RANKE = N) or DICO = 'D'.
C        Additional workspace: need   N + MAX(N,M);
C                              prefer N + MAX(N,M)*NB.
C
         CALL TG01BD( 'General', VECT, VECT, NR, M0, P0, ILO, IHI,
     $                DWORK( IA ), N, DWORK( IE ), N, DWORK( IB ), N,
     $                DWORK( IC ), P, DWORK( IQ ), N, DWORK( IZ ), N,
     $                DWORK( IWRK ), LDWORK-IWRK+1, IERR )
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
C        Perform QZ algorithm, computing eigenvalues.
C        Additional workspace: need   2*N*N + 2*N, if NSRT = .TRUE.;
C                                     2*N*N + 4*N, otherwise.
C                              prefer larger.
C
         IF( NSRT ) THEN
C
C           The generalized Hessenberg form is saved for later use.
C
            CALL DLACPY( 'Full', NR, NR, DWORK( IA ), N, DWORK( IAS ),
     $                   N )
            CALL DLACPY( 'Full', NR, NR, DWORK( IE ), N, DWORK( IES ),
     $                   N )
            CALL DHGEQZ( 'Eigenvalues', 'No Vectors', 'No Vectors', NR,
     $                   ILO, IHI, DWORK( IAS ), N, DWORK( IES ), N,
     $                   DWORK( IR ), DWORK( II ), DWORK( IBT ), DWORK,
     $                   N, DWORK, N, DWORK( IWRK ), LDWORK-IWRK+1, IERR
     $                    )
         ELSE
C
C           The generalized Schur form will be used.
C
            CALL DHGEQZ( 'Schur', 'Vectors', 'Vectors', NR, ILO, IHI,
     $                   DWORK( IA ), N, DWORK( IE ), N, DWORK( IR ),
     $                   DWORK( II ), DWORK( IBT ), DWORK( IQ ), N,
     $                   DWORK( IZ ), N, DWORK( IWRK ), LDWORK-IWRK+1,
     $                   IERR )
         END IF
C
         IF( IERR.NE.0 ) THEN
            INFO = 2
            RETURN
         ELSE IF( IERR.EQ.NR+1 ) THEN
            INFO = 5
            RETURN
         END IF
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
C        Check pencil regularity.
C
         DO 3 I = 0, NR - 1
            IF( ABS( DWORK( IBT+I ) ).LE.TZER ) THEN
               IF( DLAPY2( DWORK( IR+I ), DWORK( II+I ) ).LE.TZER ) THEN
                  INFO = 7
                  GPEAK( 1 ) = ONE
                  GPEAK( 2 ) = ZERO
                  FPEAK( 1 ) = ONE
                  FPEAK( 2 ) = ZERO
                  RETURN
               END IF
            END IF
    3    CONTINUE
C
         IF( .NOT.NSRT ) THEN
C
C           Reorder finite eigenvalues to the top and infinite
C           eigenvalues to the bottom and update B and C.
C
C           Additional workspace: need  MAX(8*N+16,N*MAX(M,P));
C                                 prefer larger.
C
C           Undo scaling on eigenvalues before selecting them.
C
            IF( ILASCL ) THEN
               CALL DLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, DWORK( IR ),
     $                      N, IERR )
               CALL DLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, DWORK( II ),
     $                      N, IERR )
            END IF
            IF( ILESCL )
     $         CALL DLASCL( 'G', 0, 0, ENRMTO, ENRM, N, 1, DWORK( IBT ),
     $                      N, IERR )
C
C           Select eigenvalues.
C
            IM    = IWRK
            WMAX  = ZERO
            WRMIN = SAFMAX
            DO 1 I = 0, NR - 1
               BWORK( I+1 ) = DWORK( IBT+I ).NE.ZERO
               IF( BWORK( I+1 ) ) THEN
                  IF( MIN( ABS( DWORK( IR+I ) ), ABS( DWORK( II+I ) ) )
     $                  .EQ.ZERO ) THEN
                     TMP = MAX( ABS( DWORK( IR+I ) ),
     $                          ABS( DWORK( II+I ) ) )
                  ELSE
                     TMP = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
                  END IF
                  DWORK( IM+I ) = TMP
                  IF( DWORK( IBT+I ).GT.ONE ) THEN
                     TMP = TMP/DWORK( IBT+I )
                  ELSE IF( TMP.LE.ONE ) THEN
                     TMP = TMP/MAX( DWORK( IBT+I ), SAFMIN )
                  ELSE IF( TMP/SAFMAX.LT.DWORK( IBT+I ) ) THEN
                     TMP = TMP/DWORK( IBT+I )
                  ELSE
                     BWORK( I+1 ) = .FALSE.
                  END IF
C
                  IF( BWORK( I+1 ) ) THEN
                     WMAX  = MAX( WMAX,  TMP )
                     WRMIN = MIN( WRMIN, TMP )
                  END IF
               ELSE
                  DWORK( IM+I ) = SAFMAX
               END IF
    1       CONTINUE
C
            IF( WRMIN.GT.WMAX*SAFMIN ) THEN
               RAT = WMAX/WRMIN
            ELSE
               RAT = SAFMAX
            END IF
            IF( RAT.LT.ONE/SQRT( EPS ) ) THEN
               THRESH = EPS
            ELSE
               THRESH = THOUSD*DBLE( NR )*EPS
            END IF
C
            SDIM1  = 0
            DO 2 I = 0, NR - 1
               IF( BWORK( I+1 ) ) THEN
                  IF( DWORK( IBT+I ).LE.THRESH*DWORK( IM+I ) )
     $               BWORK( I+1 ) = .FALSE.
                  IF( BWORK( I+1 ) )
     $               SDIM1 = SDIM1 + 1
               END IF
    2       CONTINUE
C
            IM = IES
C
            IF( SDIM1.LT.NR ) THEN
C
C              Reorder eigenvalues.
C              Additional workspace: need   4*N+16;
C                                    prefer larger.
C
               CALL DTGSEN( 0, .TRUE., .TRUE., BWORK, NR, DWORK( IA ),
     $                      N, DWORK( IE ), N, DWORK( IR ), DWORK( II ),
     $                      DWORK( IBT ), DWORK( IQ ), N, DWORK( IZ ), 
     $                      N, SDIM, CND, CND, DIF, DWORK( IWRK ),
     $                      LDWORK-IWRK+1, IDUM, 1, IERR )
               IF( IERR.EQ.1 ) THEN
C
C                 The eigenvalue computation succeeded, but the reordering
C                 failed.  Check the pencil singularity.
C
                  INFO = 5
                  DO 111 I = 0, NR - 1
                     IF( ABS( DWORK( IBT+I ) ).LE.TZER ) THEN
                        IF( DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
     $                        .LE.TZER ) THEN
                           INFO = 7
                           RETURN
                        END IF
                     END IF
  111             CONTINUE
                  RETURN
               END IF
               IF( SDIM.NE.SDIM1 ) THEN
                  INFO = 6
                  RETURN
               END IF
               MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
            ELSE
               SDIM = NR
            END IF
C
C           Apply back-permutation to Q and Z
C
            CALL DGGBAK( 'Permute', 'Left',  NR, ILO, IHI,
     $                   DWORK( ILFT ), DWORK( IRHT ), N, DWORK( IQ ),
     $                   N, IERR )
C
            CALL DGGBAK( 'Permute', 'Right', NR, ILO, IHI,
     $                   DWORK( ILFT ), DWORK( IRHT ), N, DWORK( IZ ),
     $                   N, IERR )
C
C           Update B and C.
C
            CALL DGEMM( 'Transpose', 'No Transpose', NR, M, NR, ONE,
     $                  DWORK( IQ ), N, DWORK( IB ), N, ZERO,
     $                  DWORK( IWRK ), N )
            CALL DLACPY( 'Full', NR, M, DWORK( IWRK ), N, DWORK( IB ),
     $                   N )
C
            CALL DGEMM( 'No Transpose', 'No Transpose', P, NR, NR, ONE,
     $                  DWORK( IC ), P, DWORK( IZ ), N, ZERO,
     $                  DWORK( IWRK ), P )
            CALL DLACPY( 'Full', P, NR, DWORK( IWRK ), P, DWORK( IC ),
     $                   P )
C
C           Determine the largest singular value of G(infinity).
C
            NINF  = NR - SDIM
            SDIM1 = MAX( 1, SDIM )
            IF( NINF.NE.0 ) THEN
C
C              Save E(1:SDIM,SDIM+1:N) and use the 1-norms of
C              E(SDIM+1:N,SDIM+1:N) and A(SDIM+1:N,SDIM+1:N) to decide
C              if the algebraic index of the system is 1 or larger.
C
C              Additional workspace: need   SDIM*NINF (from IES).
C
               CALL DLACPY( 'Full', SDIM, NINF, DWORK( IE+SDIM*N ), N,
     $                      DWORK( IES ), SDIM1 )
               TM   = DLANTR( '1-norm', 'Upper', 'Non-unit', NINF, NINF,
     $                        DWORK( IE+SDIM*( N+1 ) ), N, DWORK )
               TMP  = DLANHS( '1-norm', NINF, DWORK( IA+SDIM*( N+1 ) ),
     $                        N, DWORK )
               IND1 = TM.LT.DBLE( MAX( SDIM, NINF ) )*EPS*TMP
               IF( MAX( TM, TMP ).EQ.ZERO ) THEN
C
C                 The descriptor system is singular.
C
                  INFO = 7
                  RETURN
               ELSE IF( IND1 ) THEN
C
C                 The system has algebraic index one.
C                 Solve NINF linear systems of equations
C                    E(1:SDIM,1:SDIM)*Y = -E(1:SDIM,SDIM+1:N).
C                 Check first whether E(1:SDIM,1:SDIM) is nonsingular.
C
C                 Additional workspace: need   3*SDIM (from IWRK).
C                 Integer    workspace: need   SDIM.
C
                  IWRK = IES + SDIM*NINF
                  CALL DTRCON( '1-norm', 'Upper', 'Non Unit', SDIM,
     $                         DWORK( IE ), N, RCOND, DWORK( IWRK ),
     $                         IWORK, IERR )
                  IF( RCOND.LE.TEN*DBLE( SDIM )*EPS ) THEN
C
C                    The matrix E(1:SDIM,1:SDIM) is numerically singular,
C                    no safe computation of G(infinity) is possible.
C
                     INFO = 1
                     RETURN
                  END IF
                  CALL DTRSM( 'Left', 'Upper', 'No Transpose',
     $                        'Non-unit', SDIM, NINF, ONE, DWORK( IE ),
     $                        N, DWORK( IES ), SDIM1 )
                  CALL DSCAL( SDIM*NINF, -ONE, DWORK( IES ), 1 )
               ELSE
C
C                 The system has higher algebraic index.
C                 Solve the generalized Sylvester equation
C
C                    A(1:SDIM,1:SDIM)*Y + Z*A(SDIM+1:N,SDIM+1:N) +
C                       A(1:SDIM,SDIM+1:N) = 0,
C                    E(1:SDIM,1:SDIM)*Y + Z*E(SDIM+1:N,SDIM+1:N) +
C                       E(1:SDIM,SDIM+1:N) = 0,
C
C                 in order to decouple the system into its slow and fast
C                 part.
C
C                 Additional workspace: need   4*SDIM*NINF (from IES);
C                                       prefer larger.
C                 Integer    workspace: need   N+6.
C
                  IAS  = IES + SDIM*NINF
                  IWRK = IAS + SDIM*NINF
                  CALL DLACPY( 'Full', SDIM, NINF, DWORK( IA+SDIM*N ),
     $                         N, DWORK( IAS ), SDIM1 )
                  CALL DSCAL( 2*SDIM*NINF, -ONE, DWORK( IES ), 1 )
C
C                 Solve the generalized Sylvester equation.
C
                  CALL SB04OD( 'No Reduction', 'No Transpose',
     $                         'Frobenius', SDIM, NINF, DWORK( IE ), N,
     $                         DWORK( IE+SDIM*( N+1 ) ), N,
     $                         DWORK( IES ), SDIM1, DWORK( IA ), N,
     $                         DWORK( IA+SDIM*( N+1 ) ), N,
     $                         DWORK( IAS ), SDIM1, SCL, DIF, DWORK, 1,
     $                         DWORK, 1, DWORK, 1, DWORK, 1, IWORK,
     $                         DWORK( IWRK ), LDWORK-IWRK+1, IERR )
                  IF( IERR.GT.0 .OR. DIF.GT.ONE/BOUND ) THEN
C
C                    The generalized Sylvester equation is very
C                    ill-conditioned, no safe computation of G(infinity)
C                    is possible.
C
                     INFO = 1
                     RETURN
                  END IF
                  MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1,
     $                          MAXWRK )
C
C                 Estimate the condition of the transformation matrix.
C
                  CND = DLANGE( '1-norm', SDIM, NINF, DWORK( IES ),
     $                          SDIM1, DWORK )
                  IF( CND.GT.ONE/TOLER ) THEN
C
C                    The (right) transformation matrix is very
C                    ill-conditioned, no safe computation of G(infinity)
C                    is possible.
C
                     INFO = 1
                     RETURN
                  END IF
C
               END IF
C
C              Update C.
C
C              Additional workspace: need   P*NINF (from IWRK)
C
               ICS = IWRK
               IBS = ICS + P*NINF
               CALL DLACPY( 'Full', P, NINF, DWORK( IC+P*SDIM ), P,
     $                      DWORK( ICS ), P )
               CALL DGEMM( 'No Transpose', 'No Transpose', P, NINF,
     $                     SDIM,  ONE, DWORK( IC ),  P, DWORK( IES ),
     $                     SDIM1, ONE, DWORK( ICS ), P )
C
C              Compute G(infinity) =
C                 D-C(1:P,SDIM+1:N)*inv(A(SDIM+1:N,SDIM+1:N))*
C                   B(SDIM+1:N,1:M).
C
C              Additional workspace: need   NINF*( NINF+M+3 ).
C              Integer    workspace: need   2*NINF.
C
               IAS  = IBS + NINF*M
               IWRK = IAS + NINF*NINF
               CALL DLACPY( 'Full', NINF, NINF, DWORK( IA+SDIM*( N+1 ) )
     $                      , N, DWORK( IAS ), NINF )
               CALL DLACPY( 'Full', NINF, M, DWORK( IB+SDIM ), N,
     $                      DWORK( IBS ), NINF )
               TMP = DLANHS( '1-norm', NINF, DWORK( IAS ), NINF, DWORK )
               CALL MB02SD( NINF, DWORK( IAS ), NINF, IWORK, INFO )
               CALL MB02TD( '1-norm', NINF, TMP, DWORK( IAS ), NINF,
     $                      IWORK, RCOND, IWORK(NINF+1), DWORK( IWRK ),
     $                      INFO )
               IF( RCOND.LE.TEN*DBLE( NINF )*EPS ) THEN
C
C                 The matrix A(SDIM+1:N,SDIM+1:N) is numerically singular,
C                 no safe computation of G(infinity) is possible.
C
                  INFO = 1
                  RETURN
               END IF
               CALL MB02RD( 'No Transpose', NINF, M, DWORK( IAS ),
     $                      NINF, IWORK, DWORK( IBS ), NINF, IERR )
C
               ID = IAS
               IS = ID + P*M
               IF( WITHD ) THEN
C
C                 Copy D to DWORK.
C
                  CALL DLACPY( 'Full', P, M, D, LDD, DWORK( ID ), P )
                  CALL DGEMM( 'No Transpose', 'No Transpose', P, M,
     $                        NINF, -ONE, DWORK( ICS ), P, DWORK( IBS ),
     $                        NINF, ONE, DWORK( ID ), P )
               ELSE
                  CALL DGEMM( 'No Transpose', 'No Transpose', P, M,
     $                        NINF, -ONE, DWORK( ICS ), P, DWORK( IBS ),
     $                        NINF, ZERO, DWORK( ID ), P )
               END IF
            ELSE IF( WITHD ) THEN
C
C                 Copy D to DWORK.
C
                  ID = IES
                  IS = ID + P*M
               CALL DLACPY( 'Full', P, M, D, LDD, DWORK( ID ), P )
            END IF
C
            IF( NINF.NE.0 .OR. WITHD ) THEN
C
C              Compute the maximum singular value of G(infinity).
C
C              Additional workspace: need   MAX( 3*MIN(P,M) + MAX(P,M),
C              (from IWRK)                       5*MIN(P,M) );
C                                    prefer larger.
C
               IWRK = IS + MINPM
               CALL DGESVD( 'No Computation', 'No Computation', P, M,
     $                      DWORK( ID ), P, DWORK( IS ), DWORK, 1,
     $                      DWORK, 1, DWORK( IWRK ), LDWORK-IWRK+1,
     $                      IERR )
               IF( IERR.GT.0 ) THEN
C
C                 The SVD algorithm did not converge, computation of the
C                 largest singular value of G(infinity) is not possible.
C
                  INFO = 3
                  RETURN
               END IF
               MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
               GAMMAL = DWORK( IS )
            END IF
C
         END IF
C
C        Check if unscaling would cause over/underflow; if so, rescale
C        eigenvalues (DWORK( IR+I-1 ),DWORK( II+I-1 ),DWORK( IBT+I-1 ))
C        so that DWORK( IBT+I-1 ) is on the order of E(I,I) and
C        DWORK( IR+I-1 ) and DWORK( II+I-1 ) are on the order of A(I,I).
C
         IF( ILASCL ) THEN
C
            DO 30 I = 0, NR - 1
               IF( DWORK( IR+I ).NE.ZERO ) THEN
                  IF( ( DWORK( IR+I ) / SAFMAX ).GT.( ANRMTO / ANRM )
     $                                                              .OR.
     $                ( SAFMIN / DWORK( IR+I ) ).GT.( ANRM / ANRMTO )
     $              ) THEN
                     TM = ABS( DWORK( IA+I*(N+1)+1 ) / DWORK( IR+I ) )
                     DWORK( IBT+I ) = DWORK( IBT+I )*TM
                     DWORK(  IR+I ) = DWORK(  IR+I )*TM
                     DWORK(  II+I ) = DWORK(  II+I )*TM
                  ELSE IF( ( DWORK( II+I ) / SAFMAX ).GT.
     $                     ( ANRMTO / ANRM ) .OR. DWORK( II+I ).NE.ZERO
     $                                                             .AND.
     $               ( SAFMIN / DWORK( II+I ) ).GT.( ANRM / ANRMTO ) )
     $                     THEN
                     TM = ABS( DWORK( IA+(I+1)*(N+1) ) / DWORK( II+I ) )
                     DWORK( IBT+I ) = DWORK( IBT+I )*TM
                     DWORK(  IR+I ) = DWORK(  IR+I )*TM
                     DWORK(  II+I ) = DWORK(  II+I )*TM
                  END IF
               END IF
   30       CONTINUE
C
         END IF
C
         IF( ILESCL ) THEN
C
            DO 40 I = 0, NR - 1
               IF( DWORK( IBT+I ).NE.ZERO ) THEN
                  IF( ( DWORK( IBT+I ) / SAFMAX ).GT.( ENRMTO / ENRM )
     $                                                              .OR.
     $                ( SAFMIN / DWORK( IBT+I ) ).GT.( ENRM / ENRMTO )
     $              ) THEN
                     TM = ABS( DWORK( IE+I*(N+1)+1 ) / DWORK( IBT+I ) )
                     DWORK( IBT+I ) = DWORK( IBT+I )*TM
                     DWORK(  IR+I ) = DWORK(  IR+I )*TM
                     DWORK(  II+I ) = DWORK(  II+I )*TM
                  END IF
               END IF
   40       CONTINUE
C
         END IF
C
C        Undo scaling.
C
         IF( ILASCL ) THEN 
            CALL DLASCL( 'Hessenberg', 0, 0, ANRMTO, ANRM, N, N,
     $                   DWORK( IA ), N, IERR )
            CALL DLASCL( 'General', 0, 0, ANRMTO, ANRM, N, 1,
     $                   DWORK( IR ), N, IERR )
            CALL DLASCL( 'General', 0, 0, ANRMTO, ANRM, N, 1,
     $                   DWORK( II ), N, IERR )
         END IF
C
         IF( ILESCL ) THEN
            CALL DLASCL( 'Upper', 0, 0, ENRMTO, ENRM, N, N,
     $                   DWORK( IE ), N, IERR )
            CALL DLASCL( 'General', 0, 0, ENRMTO, ENRM, N, 1,
     $                   DWORK( IBT ), N, IERR )
         END IF
C
      ELSE
C
C        Standard state-space system.
C
         IF( LEQUIL ) THEN
C
C           Equilibrate the system.
C
            MAXRED = HUNDRD
            CALL TB01ID( 'All', N, M, P, MAXRED, DWORK( IA ), N,
     $                   DWORK( IB ), N,  DWORK( IC ), P, DWORK( II ),
     $                   IERR )
         END IF
C
C        For efficiency of later calculations, the system (A,B,C) is
C        reduced to a similar one with the state matrix in Hessenberg
C        form.
C
C        First, permute the matrix A to make it more nearly triangular
C        and apply the permutations to B and C.
C
         CALL DGEBAL( 'Permute', N, DWORK( IA ), N, ILO, IHI,
     $                DWORK( IR ), IERR )
C
         DO 50 I = N, IHI + 1, -1
            K = DWORK( IR+I-1 )
            IF( K.NE.I ) THEN
               CALL DSWAP( M, DWORK( IB+I-1 ), N,
     $                        DWORK( IB+K-1 ), N )
               CALL DSWAP( P, DWORK( IC+(I-1)*P ), 1,
     $                        DWORK( IC+(K-1)*P ), 1 )
            END IF
   50    CONTINUE
C
         DO 60 I = 1, ILO - 1
            K = DWORK( IR+I-1 )
            IF( K.NE.I ) THEN
               CALL DSWAP( M, DWORK( IB+I-1 ), N,
     $                        DWORK( IB+K-1 ), N )
               CALL DSWAP( P, DWORK( IC+(I-1)*P ), 1,
     $                        DWORK( IC+(K-1)*P ), 1 )
            END IF
   60    CONTINUE
C
C        Reduce A to upper Hessenberg form and apply the transformations
C        to B and C.
C        Additional workspace: need   N;   (from II)
C                              prefer N*NB.
C
         ITAU = IR
         IWRK = II
         CALL DGEHRD( N, ILO, IHI, DWORK( IA ), N, DWORK( ITAU ),
     $                DWORK( IWRK ), LDWORK-IWRK+1, IERR )
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
C        Additional workspace: need   M;
C                              prefer M*NB.
C
         CALL DORMHR( 'Left', 'Transpose', N, M, ILO, IHI, DWORK( IA ),
     $                N, DWORK( ITAU ), DWORK( IB ), N, DWORK( IWRK ),
     $                LDWORK-IWRK+1, IERR )
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
C        Additional workspace: need   P;
C                              prefer P*NB.
C
         CALL DORMHR( 'Right', 'NoTranspose', P, N, ILO, IHI,
     $                DWORK( IA ), N, DWORK( ITAU ), DWORK( IC ), P,
     $                DWORK( IWRK ), LDWORK-IWRK+1, IERR )
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
C        Compute the eigenvalues. The Hessenberg form is saved for
C        later use.
C        Additional workspace:  need   N*N + N;   (from IBT)
C                               prefer larger.
C
         IAS  = IBT
         IWRK = IAS + NN
         IM   = IAS
         CALL DLACPY( 'Full', N, N, DWORK( IA ), N, DWORK( IAS ), N )
         CALL DHSEQR( 'Eigenvalues', 'No Vectors', N, ILO, IHI,
     $                DWORK( IAS ), N, DWORK( IR ), DWORK( II ), DWORK,
     $                N, DWORK( IWRK ), LDWORK-IWRK+1, IERR )
         IF( IERR.GT.0 ) THEN
            INFO = 2
            RETURN
         END IF
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
C
C        Annihilate the lower part of the Hessenberg matrix.
C
         IF( N.GT.2 )
     $      CALL DLASET( 'Lower', N-2, N-2, ZERO, ZERO, DWORK( IA+2 ),
     $                   N )
C
         IF( ILASCL ) THEN
C
C           Undo scaling for the Hessenberg form of A and eigenvalues.
C
            CALL DLASCL( 'Hessenberg', 0, 0, ANRMTO, ANRM, N, N,
     $                   DWORK( IA ), N, IERR )
            CALL DLASCL( 'General', 0, 0, ANRMTO, ANRM, N, 1,
     $                   DWORK( IR ), N, IERR )
            CALL DLASCL( 'General', 0, 0, ANRMTO, ANRM, N, 1,
     $                   DWORK( II ), N, IERR )
         END IF
C
      END IF
C
C     Look for (generalized) eigenvalues on the boundary of the
C     stability domain. (Their existence implies an infinite norm.)
C     Additional workspace:  need   2*N.   (from IM)
C
      IAR   = IM + N
      IMIN  = II
      WRMIN = SAFMAX
C
C     NEI defines the number of finite eigenvalues. An eigenvalue is 
C     declared finite if it is less than ONE/BOUND.
C
      NEI = 0
C
      IF( DISCR ) THEN
C
C        For discrete-time case, compute the logarithm of the finite 
C        non-zero eigenvalues and save their moduli and absolute real
C        parts. (The logarithms are overwritten on the eigenvalues.)
C        Also, find the minimum distance to the unit circle.
C
         IF( FULLE ) THEN
            DO 70 I = 0, NR - 1
               TM = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
               IF( TM.LT.ONE/BOUND*DWORK( IBT+I ) ) THEN
C
C                 Finite eigenvalues.
C
                  TM = TM / DWORK( IBT+I )
                  IF( TM.NE.ZERO ) THEN
                     DWORK( II+NEI ) = ATAN2( DWORK( II+I ),
     $                                 DWORK( IR+I ) )
                     DWORK( IR+NEI ) = LOG( TM )
                  ELSE
                     DWORK( II+NEI ) = DWORK( II+I )
                     DWORK( IR+NEI ) = DWORK( IR+I )
                  END IF
                  DWORK( IM ) = DLAPY2( DWORK( IR+NEI ),
     $                          DWORK( II+NEI ) )
                  TM = ABS( ONE - TM )
                  IF( TM.LT.WRMIN ) THEN
                     IMIN  = II + NEI
                     WRMIN = TM
                  END IF
                  IM = IM + 1
                  DWORK( IBT+NEI ) = DWORK( IBT+I )               
                  DWORK( IAR+NEI ) = ABS( DWORK( IR+NEI ) )
                  NEI = NEI + 1
               END IF
   70       CONTINUE
C
         ELSE
C
            DO 80 I = 0, NR - 1
               TM = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
               IF( TM.NE.ZERO ) THEN
                  DWORK( II+I ) = ATAN2( DWORK( II+I ), DWORK( IR+I ) )
                  DWORK( IR+I ) = LOG( TM )
               END IF
               DWORK( IM ) = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
               TM = ABS( ONE - TM )
               IF( TM.LT.WRMIN ) THEN
                  IMIN  = II + I
                  WRMIN = TM
               END IF
               IM = IM + 1
               DWORK( IAR+I ) = ABS( DWORK( IR+I ) )
   80       CONTINUE
C
            NEI = NR
C
         END IF
C
      ELSE
C
C        For continuous-time case, save moduli and absolute real parts
C        of finite eigenvalues and find the maximum modulus and minimum
C        absolute real part.
C
         LINF = .FALSE.
         WMAX = ZERO
C
         IF( FULLE ) THEN
C
            DO 90 I = 0, NR - 1
               TM = ABS( DWORK( IR+I ) )
               DWORK( IM ) = DLAPY2( TM, DWORK( II+I ) )
               IF( THRESH*DWORK( IM ).LT.DWORK( IBT+I ) ) THEN
                  TM = TM / DWORK( IBT+I )
                  DWORK( IM ) = DWORK( IM ) / DWORK( IBT+I )
               ELSE
C
C                 The pencil has infinite eigenvalues. SAFMAX is used.
C
                  TM = SAFMAX
                  DWORK( IM ) = SAFMAX
               END IF
               IF( TM.LT.WRMIN ) THEN
                  IMIN  = II + NEI
                  WRMIN = TM
                  IF( WRMIN.LE.TEN*EPS*DWORK( IM ) .OR.
     $              ( WRMIN.LE.TEN*EPS .AND. DWORK( IMIN ).EQ.ZERO ) )
     $               LINF = .TRUE.
               END IF
               DWORK( IAR+NEI ) = TM
               IF( DWORK( IM ).NE.SAFMAX .AND. DWORK( IM ).GT.WMAX )
     $            WMAX = DWORK( IM )
               IM  = IM  + 1
               NEI = NEI + 1
   90       CONTINUE
C
         ELSE
C
            DO 100 I = 0, NR - 1
               TM = ABS( DWORK( IR+I ) )
               DWORK( IM ) = DLAPY2( TM, DWORK( II+I ) )
               IF( TM.LT.WRMIN .OR. DWORK( IM ).EQ.ZERO ) THEN
                  IMIN  = II + I
                  WRMIN = TM
                  IF( WRMIN.LE.TEN*EPS*DWORK( IM ) .OR.
     $              ( WRMIN.LE.TEN*EPS .AND. DWORK( IMIN ).EQ.ZERO ) )
     $               LINF = .TRUE.
               END IF
               IF( DWORK( IM ).GT.WMAX )
     $            WMAX = DWORK( IM )
               IM = IM + 1
               DWORK( IAR+I ) = TM
  100       CONTINUE
C
            NEI = NR
         END IF
C
         BOUND = BOUND + EPS*WMAX
C
      END IF
C
      IM = IM - NEI
C
      IF( LINF ) THEN
C
C        The L-infinity norm was found as infinite.
C
         IF ( DISCR )
     $      THRESH = BOUND
         GPEAK( 1 ) = ONE
         GPEAK( 2 ) = ZERO
         TM = ABS( DWORK( IMIN ) )
         IF ( FULLE ) THEN
            IF ( TM*THRESH.LT.DWORK( IBT+IMIN-II ) ) THEN
               TM = TM/DWORK( IBT+IMIN-II )
            ELSE
               TM = SAFMAX
            END IF
         END IF
         IF ( DISCR )
     $      TM = ABS( ATAN2( SIN( TM ), COS( TM ) ) )
         FPEAK( 1 ) = TM
         IF ( TM.LT.SAFMAX ) THEN
            FPEAK( 2 ) = ONE
         ELSE
            FPEAK( 2 ) = ZERO
         END IF
C
         DWORK( 1 ) = MAXWRK
         ZWORK( 1 ) = ONE
         RETURN
      END IF
C
C     Determine the maximum singular value of
C        G(lambda) = C*inv(lambda*E - A)*B + D,
C     over a selected set of frequencies. Besides the frequencies w = 0,
C     w = pi (if DICO = 'D'), and the given value FPEAK, this test set
C     contains the peak frequency for each mode (or an approximation
C     of it). The (generalized) Hessenberg form of the system is used.
C
C     First, determine the maximum singular value of G(0) and set FPEAK
C     accordingly.
C     Additional workspace:
C           complex: need   1, if DICO = 'C';
C                           (N+M)*(N+P)+2*MIN(P,M)+MAX(P,M)), otherwise;
C                    prefer larger;
C           real:    need   LDW1+LDW2+LDW3+N, where
C                           LDW1 = N*N+N*M+P*M, if DICO = 'C';
C                           LDW1 = 0,           if DICO = 'D';
C                           LDW2 = MIN(P,M)+MAX(3*MIN(P,M)+MAX(P,M),
C                                               5*MIN(P,M)),
C                                              if DICO = 'C';
C                           LDW2 = 6*MIN(P,M), if DICO = 'D';
C                           LDW3 = P*M,   if DICO = 'D', JOBD = 'D';
C                           LDW3 = 0,     otherwise.
C                    prefer larger.
C     Integer    workspace: need   N.
C
      IF ( FULLE ) THEN
         JOB = 'G'
      ELSE
         JOB = 'I'
      END IF
C
      IF ( DISCR ) THEN
         IAS = IA
         IBS = IB
         ID  = IAR + N
         IF( WITHD ) THEN
            IWRK = ID  + P*M
            CALL DLACPY( 'Full', P, M, D, LDD, DWORK( ID ), P )
         ELSE
            IWRK = ID
         END IF
      ELSE
         RTOL = HUNDRD*TOLER
         IAS  = IAR + N
         IBS  = IAS + NN
         IF( WITHD ) THEN
            ID   = IBS + N*M
            IWRK = ID  + P*M
            CALL DLACPY( 'Full', P, M, D, LDD, DWORK( ID ), P )
         ELSE
            ID   = IBS
            IWRK = IBS + N*M
         END IF
         CALL DLACPY( 'Upper', NR, NR, DWORK( IA ), N, DWORK( IAS ), N )
         CALL DCOPY(  NR-1, DWORK( IA+1 ), N+1, DWORK( IAS+1 ), N+1 )
         CALL DLACPY( 'Full', NR, M, DWORK( IB ), N, DWORK( IBS ), N )
      END IF
C
      GAMMA = AB13DX( DICO, JOB, JOBD, NR, M, P, ZERO, DWORK( IAS ), N,
     $                DWORK( IE ), N, DWORK( IBS ), N, DWORK( IC ), P,
     $                DWORK( ID ), P, IWORK, DWORK( IWRK ),
     $                LDWORK-IWRK+1, ZWORK, LZWORK, IERR )
      MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
      IF( IERR.GE.1 .AND. IERR.LE.NR ) THEN
         GPEAK( 1 ) = ONE
         FPEAK( 1 ) = ZERO
         GPEAK( 2 ) = ZERO
         FPEAK( 2 ) = ONE
         GO TO 240
      ELSE IF( IERR.EQ.NR+1 ) THEN
         INFO = 3
         RETURN
      END IF
C
      FPEAKS = FPEAK( 1 )
      FPEAKI = FPEAK( 2 )
      IF( GAMMAL.LT.GAMMA ) THEN
         GAMMAL     = GAMMA
         FPEAK( 1 ) = ZERO
         FPEAK( 2 ) = ONE
      ELSE IF( .NOT.DISCR ) THEN
         FPEAK( 1 ) = ONE
         FPEAK( 2 ) = ZERO
      END IF
C
      MAXCWK = INT( ZWORK( 1 ) )
C
      IF( DISCR ) THEN
C
C        Try the frequency w = pi.
C
         PI    = FOUR*ATAN( ONE )
         GAMMA = AB13DX( DICO, JOB, JOBD, NR, M, P, PI, DWORK( IA ),
     $                   N, DWORK( IE ), N, DWORK( IB ), N, DWORK( IC ),
     $                   P, D, LDD, IWORK, DWORK( IWRK ), LDWORK-IWRK+1,
     $                   ZWORK, LZWORK, IERR )
         MAXCWK = MAX( INT( ZWORK( 1 ) ), MAXCWK )
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
         IF( IERR.GE.1 .AND. IERR.LE.NR ) THEN
            GPEAK( 1 ) = ONE
            FPEAK( 1 ) = PI
            GPEAK( 2 ) = ZERO
            FPEAK( 2 ) = ONE
            GO TO 240
         ELSE IF( IERR.EQ.NR+1 ) THEN
            INFO = 3
            RETURN
         END IF
C
         IF( GAMMAL.LT.GAMMA ) THEN
            GAMMAL     = GAMMA
            FPEAK( 1 ) = PI
            FPEAK( 2 ) = ONE
         END IF
C
      ELSE
         IWRK = IAS
      END IF
C
C     Build the remaining set of frequencies.
C     Complex workspace:  need   (N+M)*(N+P)+2*MIN(P,M)+MAX(P,M));
C                         prefer larger.
C     Real workspace:     need   LDW2 = 6*MIN(P,M) (from IWRK);
C                         prefer larger.
C
      IF ( MIN( FPEAKS, FPEAKI ).NE.ZERO ) THEN
C
C        Compute also the norm at the given (finite) frequency.
C
         GAMMA = AB13DX( DICO, JOB, JOBD, NR, M, P, FPEAKS, DWORK( IA ),
     $                   N, DWORK( IE ), N, DWORK( IB ), N, DWORK( IC ),
     $                   P, D, LDD, IWORK, DWORK( IWRK ), LDWORK-IWRK+1,
     $                   ZWORK, LZWORK, IERR )
         MAXCWK = MAX( INT( ZWORK( 1 ) ), MAXCWK )
         MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
         IF ( DISCR ) THEN
            TM = ABS( ATAN2( SIN( FPEAKS ), COS( FPEAKS ) ) )
         ELSE
            TM = FPEAKS
         END IF
         IF( IERR.GE.1 .AND. IERR.LE.NR ) THEN
            GPEAK( 1 ) = ONE
            FPEAK( 1 ) = TM
            GPEAK( 2 ) = ZERO
            FPEAK( 2 ) = ONE
            GO TO 240
         ELSE IF( IERR.EQ.NR+1 ) THEN
            INFO = 3
            RETURN
         END IF
C
         IF( GAMMAL.LT.GAMMA ) THEN
            GAMMAL     = GAMMA
            FPEAK( 1 ) = TM
            FPEAK( 2 ) = ONE
         END IF
C
      END IF
C
      DO 110 I = 0, NEI - 1
         IF( DWORK( II+I ).GE.ZERO .AND. DWORK( IM+I ).GT.ZERO .AND.
     $       DWORK( IM+I )*THRESH.LE.ONE ) THEN
            IF ( ( DWORK( IM+I ).GE.ONE ) .OR. ( DWORK( IM+I ).LT.ONE
     $            .AND. DWORK( IAR+I ).LT.SAFMAX*DWORK( IM+I ) ) ) THEN
               RAT = DWORK( IAR+I ) / DWORK( IM+I )
            ELSE
               RAT = ONE
            END IF
C
            OMEGA = DWORK( IM+I )*SQRT( MAX( P25, ONE - TWO*RAT**2 ) )
C
            GAMMA = AB13DX( DICO, JOB, JOBD, NR, M, P, OMEGA,
     $                      DWORK( IA ), N, DWORK( IE ), N, DWORK( IB ),
     $                      N, DWORK( IC ), P, D, LDD, IWORK,
     $                      DWORK( IWRK ), LDWORK-IWRK+1, ZWORK, LZWORK,
     $                      IERR )
            MAXCWK = MAX( INT( ZWORK( 1 ) ), MAXCWK )
            MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
            IF ( DISCR ) THEN
               TM = ABS( ATAN2( SIN( OMEGA ), COS( OMEGA ) ) )
            ELSE
               TM = OMEGA
            END IF
            IF( IERR.GE.1 .AND. IERR.LE.NR ) THEN
               GPEAK( 1 ) = ONE
               FPEAK( 1 ) = TM
               GPEAK( 2 ) = ZERO
               FPEAK( 2 ) = ONE
               GO TO 240
            ELSE IF( IERR.EQ.NR+1 ) THEN
               INFO = 3
               RETURN
            END IF
C
            IF( GAMMAL.LT.GAMMA ) THEN
               GAMMAL     = GAMMA
               FPEAK( 1 ) = TM
               FPEAK( 2 ) = ONE
            END IF
C
         END IF
  110 CONTINUE
C
C     Return if the lower bound is zero.
C
      IF( GAMMAL.EQ.ZERO ) THEN
         GPEAK( 1 ) = ZERO
         FPEAK( 1 ) = ZERO
         GPEAK( 2 ) = ONE
         FPEAK( 2 ) = ONE
         GO TO 240
      END IF
C
C     Start the modified gamma iteration for the Bruinsma-Steinbuch
C     algorithm.
C
      ITER = 0
C
C     Use the structure-preserving embedding method on a
C     skew-Hamiltonian/Hamiltonian pencil.
C
C     Add at most one auxiliary variable.
C
      K = ( M+P+R )/2
      Q = NR - K
C
C     WHILE ( Iteration may continue ) DO
C    
  120 CONTINUE
C
         ITER  = ITER + 1
         GAMMA = ( ONE + TOLN )*GAMMAL
         USEPEN = FULLE .OR. DISCR
         IF ( .NOT.USEPEN ) THEN
            TMP = MAX( BNORM, CNORM )*SQRT( GAMMA )
            IF ( WITHD ) THEN
C
C              Check whether one can use an explicit Hamiltonian matrix:
C              compute
C              min(rcond(GAMMA**2*Im - S'*S), rcond(GAMMA**2*Ip - S*S')).
C              If P = M = 1, then GAMMA**2 - S(1)**2 is used instead.
C
               IF ( M.NE.P ) THEN
                  RCOND = ONE - ( SV1 / GAMMA )**2
               ELSE IF ( MINPM.GT.1 ) THEN
                  RCOND = ( GAMMA**2 - SV1**2 ) / ( GAMMA**2 - SVP**2 )
               ELSE
                  RCOND = GAMMA**2 - SV1**2
               END IF
C
               USEPEN = MAX( RCOND, TMP ).LT.RTOL
            ELSE
               USEPEN = TMP.LT.RTOL
            END IF
         END IF
C
         IF ( .NOT.USEPEN ) THEN
            IF ( .NOT.WITHD ) THEN
C
C              Standard continuous-time case with D = 0.
C              Form the needed part of the Hamiltonian matrix explicitly:
C                 H = H11 - H12*inv(H22)*H21/g.
C              Additional workspace: need   2*N*N+N.   (from IBT)
C
               IH   = IBT
               IH12 = IH   + NN
               ISL  = IH12 + NN + N
               CALL DLACPY( 'Full', N, N, A, LDA, DWORK( IH ), N )
C
C              Compute triangles of -C'*C/GAMMA and B*B'/GAMMA.
C
               CALL DSYRK( 'Lower', 'Transpose', N, P, -ONE/GAMMA, C,
     $                     LDC, ZERO, DWORK( IH12 ), N )
               CALL DSYRK( 'Upper', 'No Transpose', N, M, ONE/GAMMA, B,
     $                     LDB, ZERO, DWORK( IH12+N ), N )
C
            ELSE
C
C              Standard continuous-time case with D <> 0 and the SVD of
C              D can be used. Compute explicitly the needed part of the
C              Hamiltonian matrix:
C
C              H =
C               (A+B1*S'*inv(g^2*Ip-S*S')*C1' g*B1*inv(g^2*Im-S'*S)*B1')
C               (                                                      )
C               (  -g*C1*inv(g^2*Ip-S*S')*C1'            -H11'         )
C
C              where g = GAMMA, B1 = B*V, C1 = C'*U, and H11 is the first
C              block of H.
C              Primary additional workspace: need   2*N*N+N   (from IBT)
C              (for building the relevant part of the Hamiltonian matrix).
C
C              Compute C1*sqrt(inv(g^2*Ip-S*S')) .
C              Additional workspace: need   MAX(M,P)+N*P.
C
               IH   = IBT
               IH12 = IH   + NN
               ISL  = IH12 + NN + N
C
               DO 255 I = 0, MINPM - 1
                  DWORK( ISL+I ) = ONE/SQRT( GAMMA**2 -
     $                                       DWORK( IS+I )**2 )
  255          CONTINUE
C
               IF ( M.LT.P ) THEN
                  DWORK( ISL+M ) = ONE / GAMMA
                  CALL DCOPY( P-M-1, DWORK( ISL+M ), 0,
     $                        DWORK( ISL+M+1 ), 1 )
               END IF
               ISC = ISL + MAX( M, P )
               CALL DLACPY( 'Full', N, P, DWORK( ICU ), N, DWORK( ISC ),
     $                      N )
               CALL MB01SD( 'Column', N, P, DWORK( ISC ), N, DWORK,
     $                      DWORK( ISL ) )
C
C              Compute B1*S' .
C              Additional workspace: need   N*M.
C
               ISB = ISC + P*N
               CALL DLACPY( 'Full', N, M, DWORK( IBV ), N, DWORK( ISB ),
     $                      N )
               CALL MB01SD( 'Column', N, MINPM, DWORK( ISB ), N, DWORK,
     $                      DWORK( IS ) )
C
C              Compute B1*S'*sqrt(inv(g^2*Ip-S*S')) .
C
               CALL MB01SD( 'Column', N, MINPM, DWORK( ISB ), N, DWORK,
     $                      DWORK( ISL ) )
C
C              Compute H11 .
C
               CALL DLACPY( 'Full', N, N, A, LDA, DWORK( IH ), N )
               CALL DGEMM( 'No Transpose', 'Transpose', N, N, MINPM,
     $                     ONE, DWORK( ISB ), N, DWORK( ISC ), N, ONE,
     $                     DWORK( IH ), N )
C
C              Compute B1*sqrt(inv(g^2*Im-S'*S)) .
C
               IF ( P.LT.M ) THEN
                  DWORK( ISL+P ) = ONE / GAMMA
                  CALL DCOPY( M-P-1, DWORK( ISL+P ), 0,
     $                        DWORK( ISL+P+1 ), 1 )
               END IF
               CALL DLACPY( 'Full', N, M, DWORK( IBV ), N, DWORK( ISB ),
     $                      N )
               CALL MB01SD( 'Column', N, M, DWORK( ISB ), N, DWORK,
     $                      DWORK( ISL ) )
C
C              Compute the lower triangle of H21 and the upper triangle
C              of H12.
C
               CALL DSYRK( 'Lower', 'No Transpose', N, P, -GAMMA,
     $                     DWORK( ISC ), N, ZERO, DWORK( IH12 ), N )
               CALL DSYRK( 'Upper', 'No Transpose', N, M, GAMMA,
     $                     DWORK( ISB ), N, ZERO, DWORK( IH12+N ), N )
            END IF
C
C           Compute the eigenvalues of the Hamiltonian matrix by the
C           symplectic URV and the periodic Schur decompositions.
C           Additional workspace: need   (2*N+8)*N;
C                                 prefer larger.
C
            IWRK = ISL + NN
            CALL MB03XD( 'Both',  'Eigenvalues', 'No vectors',
     $                   'No vectors', N, DWORK( IH ), N, DWORK( IH12 ),
     $                   N, DWORK( ISL ), N, TEMP, 1, TEMP, 1, TEMP, 1,
     $                   TEMP, 1, DWORK( IR ), DWORK( II ), ILO,
     $                   DWORK( IWRK ), DWORK( IWRK+N ),
     $                   LDWORK-IWRK-N+1, IERR )
            IF( IERR.GT.0 ) THEN
               INFO = 2
               RETURN
            END IF
            MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK + N - 1, MAXWRK )
C
         ELSE
C
            IF( .NOT.DISCR .AND. CMPRE ) THEN
               NE = RANKE
            ELSE
               NE = NR
            END IF
C
            NBLK = NR   + ( M + P + R )/2
            II   = IR   + NBLK
            IBT  = II   + NBLK
            IH   = IBT  + NBLK
            IH12 = IH   + NBLK*NBLK
            IJ   = IH12 + NBLK*( NBLK + 1 )
            IJ12 = IJ   + NBLK*NBLK
            IT   = IJ12 + NBLK*( NBLK + 1 )
            IT12 = IT   + NBLK*NBLK
            IH22 = IT12 + NBLK*NBLK
            IWRK = IH22 + NBLK*NBLK
C
C           Initialize the pencil by zero!
C
            CALL DLASET( 'Full', 4*NBLK*NBLK+2*NBLK, 1, ZERO, ZERO,
     $                   DWORK( IH ), 1 )
C
            IF( DISCR ) THEN          
C
C              Set up the needed parts of the skew-Hamiltonian/
C              Hamiltonian pencil (J,H),
C
C                  ( H11  H12 )         ( J11  J12 )
C              H = (          ),    J = (          ),
C                  ( H21  H22 )         ( J21  J22 )
C
C              with
C
C                    ( -A+E  0  )          (  0   B  )
C              H11 = (          ),   H12 = (         ),
C                    (   0   D' )          ( B' -g*I )
C
C                    ( 0   C' )
C              H21 = (        ),     H22 = -H11',   
C                    ( C  g*I )
C
C                    ( A+E  0 )            
C              J11 = (        ),     J12 = 0, 
C                    (  0   0 )
C
C                    (  0  C' )
C              J21 = (        ),     J22 = J11',
C                    ( -C  0  )
C
C              where B, C, and D are extended such that the number of
C              inputs and outputs are equal (= MAX(P,M)), and g = GAMMA.
C
C              Additional workspace: need   4*NBLK*NBLK+5*NBLK (from IR)
C
               IF( FULLE ) THEN
                  DO 140 I = 1, Q
                     IN = ( I - 1 )*NBLK
                     DO 130 J = 0, NR - 1
                        DWORK( IH+IN+J ) =  DWORK( IA+(K+I-1)*N+J ) -
     $                                      DWORK( IE+(K+I-1)*N+J )
                        DWORK( IJ+IN+J ) = -DWORK( IA+(K+I-1)*N+J ) -
     $                                      DWORK( IE+(K+I-1)*N+J )
  130                CONTINUE
  140             CONTINUE
               ELSE
                  CALL DLACPY( 'Full', NR, Q, DWORK( IA+K*N ), N,
     $                         DWORK( IH ), NBLK )
                  DO 160 I = 1, Q
                     IN = ( I - 1 )*NBLK
                     DO 150 J = 0, NR - 1
                        DWORK( IJ+IN+J ) = -DWORK( IA+(K+I-1)*N+J )
  150                CONTINUE
  160             CONTINUE
                  DO 170 I = 1, Q
                     IN = ( I - 1 )*NBLK + K + I - 1
                     DWORK( IH+IN ) = DWORK( IH+IN ) - ONE
                     DWORK( IJ+IN ) = DWORK( IJ+IN ) - ONE
  170             CONTINUE
               END IF
C
C              Construct the rest of H11.
C
               CALL MA02AD( 'Full', P, K, DWORK( IC ), P,
     $                      DWORK( IH+Q*NBLK+NR ), NBLK )
               DO 300 I = 0, M - 1
                  DO 180 J = 0, NR - 1
                     DWORK( IH+(Q+P+I)*NBLK+J ) = -DWORK( IB+I*N+J )
  180             CONTINUE
  300          CONTINUE
C
C              Construct the rest of J11.
C
               CALL MA02AD( 'Full', P, K, DWORK( IC ), P,
     $                      DWORK( IJ+Q*NBLK+NR ), NBLK )
C
C              Construct the lower triangular part of H12.
C
               CALL DLACPY( 'Full', P, Q, DWORK( IC+K*P ), P,
     $                      DWORK( IH12+Q ), NBLK )
               CALL DLASET( 'Lower', P+M, P+M, ZERO, GAMMA,
     $                      DWORK( IH12+Q*NBLK+Q ), NBLK )
               IF( WITHD ) THEN
                  DO 320 I = 0, P - 1
                     DO 310 J = 0, M - 1
                        DWORK( IH12+(Q+I)*NBLK+Q+P+J ) = -D( I, J )
  310                CONTINUE
  320             CONTINUE
               END IF
               IF( R.EQ.1 )
     $            DWORK( IH12+NBLK*NBLK-1 ) = ONE
C
C              Construct the lower triangular part of J12.
C
               DO 340 I = 1, Q
                  IN = ( I - 1)*NBLK + Q
                  DO 330 J = 0, P - 1
                     DWORK( IJ12+IN+J ) = -DWORK( IC+(K+I-1)*P+J )
  330             CONTINUE
  340          CONTINUE
C
C              Construct the upper triangular parts of H12 and J12.
C
               IF( FULLE ) THEN
                  DO 360 I = 0, K - 1
                     DO 350 J = 0, NR - 1
                        DWORK( IH12+(NR+I)*NBLK+J ) = -DWORK( IA+I*N+J )
     $                                               + DWORK( IE+I*N+J )
                        DWORK( IJ12+(NR+I)*NBLK+J ) =  DWORK( IA+I*N+J )
     $                                               + DWORK( IE+I*N+J )
  350                CONTINUE
  360             CONTINUE
               ELSE
                  DO 380 I = 0, K - 1
                     DO 370 J = 0, NR - 1
                        DWORK( IH12+(NR+I)*NBLK+J ) = -DWORK( IA+I*N+J )
  370                CONTINUE
  380             CONTINUE
                  CALL DLACPY( 'Full', NR, K, DWORK( IA ), N,
     $                         DWORK( IJ12+(NR+1)*NBLK ), NBLK )
                  DO 390 I = 0, K - 1
                     DWORK( IH12+(NR+I)*NBLK+I )
     $                  = DWORK( IH12+(NR+I)*NBLK+I ) + ONE
                     DWORK( IJ12+(NR+I)*NBLK+I )
     $                  = DWORK( IJ12+(NR+I)*NBLK+I ) + ONE
  390             CONTINUE
               END IF
C
            ELSE
C
C              Set up the needed parts of the skew-Hamiltonian/
C              Hamiltonian pencil (H,J),
C
C                  ( H11  H12 )        ( S11   0  )        (  0   I )
C              H = (          ) ,  S = (          ) ,  J = (        ) ,
C                  ( H21  H22 )        (  0   S22 )        ( -I   0 )
C
C              with
C
C                    ( A  B )            ( 0   0  )            ( E  0 )
C              H11 = (      ),     H12 = (        ),     S11 = (      ),
C                    ( C  D )            ( 0 -g*I )            ( 0  0 )
C
C              H21 = -H12,         H22 = -H11',          S22 = S11',
C
C              where B and D are extended with one zero column if M+P is
C              odd, and g = GAMMA.
C
C              Additional workspace: need   4*NBLK*NBLK+5*NBLK (from IR)
C
C              Construct H11.
C
               IF( Q.GE.0 ) THEN
                  ICI = IC
                  IHC = IH + Q*NBLK + NR
                  NCC = P
               ELSE IF( K.LE.NR+P ) THEN
                  ICI = IC + K - NR
                  IHC = IH + NR
                  NCC = Q  + P
               ELSE
                  IBI = IB + ( K - NR - P )*NR 
                  IDI = K - NR - P + 1
               END IF
               IF( Q.GT.0 )
     $            CALL DLACPY( 'Full', NR, Q, DWORK( IA+K*N ), N,
     $                         DWORK( IH ), NBLK )
               IF( K.LE.NR+P ) THEN
                  CALL MA02AD( 'Full', NCC, MIN( NR, K ), DWORK( ICI ),
     $                         P, DWORK( IHC ), NBLK )
C
                  DO 430 I = 0, M - 1
                     DO 420 J = 0, NR - 1
                        DWORK( IH+(Q+P+I)*NBLK+J ) = -DWORK( IB+I*N+J )
  420                CONTINUE
  430             CONTINUE
                  IF( WITHD .AND. Q.LT.0 )
     $               CALL DLACPY( 'Full', -Q, M, D, LDD,
     $                            DWORK( IH+(Q+P)*NBLK+2*NR ), NBLK )
               ELSE
                  CALL DLACPY( 'Full', NR, NR+K-R, DWORK( IBI ), N,
     $                         DWORK( IH ), NBLK )
C
                  IF( WITHD )
     $               CALL DLACPY( 'Full', P, NR+K-R, D( 1, IDI ), LDD,
     $                            DWORK( IH+2*NR ), NBLK )
               END IF
C
C              Construct J11.
C
               IF( Q.GT.0 ) THEN
                  IF( FULLE ) THEN
                     IF( NE.GT.K )
     $                  CALL DLACPY( 'Full', NE, NE-K, DWORK( IE+K*N ),
     $                               N, DWORK( IJ ), NBLK )
                  ELSE
                     TEMP( 1 ) = ONE
                     CALL DCOPY( Q, TEMP, 0, DWORK( IJ ), NBLK+1 )
                  END IF
               END IF
C
C              Construct the lower triangular part of H12.
C
               IF( Q.GT.0 ) THEN
                  DO 436 I = 0, Q - 1
                     DO 435 J = 0, P - 1
                        DWORK( IH12+Q+I*NBLK+J ) =
     $                        -DWORK( IC+(K+I)*P+J )
  435                CONTINUE
  436             CONTINUE
                  IQ = Q*NBLK + Q
                  PM = P + M
               ELSE
                  IQ = 0
                  PM = P + M + Q
               END IF
               TEMP( 1 ) = GAMMA
               CALL DCOPY( PM, TEMP, 0, DWORK( IH12+IQ ), NBLK+1 )
               IF( WITHD ) THEN
                  IF( Q.GE.0 ) THEN
                     CALL MA02AD( 'Full', P, M, D, LDD,
     $                            DWORK( IH12+Q*NBLK+Q+P ), NBLK )
                  ELSE IF( K.LE.NR+P ) THEN
                     CALL MA02AD( 'Full', P+Q, M, D( 1-Q, 1 ), LDD,
     $                            DWORK( IH12+Q+P ), NBLK )
                  END IF
               END IF
               IF( R.EQ.1 )
     $            DWORK( IH12+NBLK*NBLK-1 ) = ONE
C
C              Construct the upper triangular parts of H12 and J12.
C
               CALL DLACPY( 'Full', NR, MIN( NR, K ), DWORK( IA ), N,
     $                      DWORK( IH12+(NR+1)*NBLK ), NBLK )
               IF( Q.LT.0 ) THEN
                  TEMP( 1 ) = -GAMMA
                  CALL DCOPY( -Q, TEMP, 0,
     $                        DWORK( IH12+(2*NR+1)*NBLK+2*NR ), NBLK+1 )
                  IF( K.LE.NR+P ) THEN
                     DO 454 I = 0, K - NR - 1
                        DO 453 J = 0, NR - 1
                           DWORK( IH12+(2*NR+1+I)*NBLK+NR+J ) =
     $                        -DWORK( IC+J*P+I )
  453                   CONTINUE
  454                CONTINUE
                  ELSE
                     CALL MA02AD( 'Full', P, NR, DWORK( IC ), P,
     $                            DWORK( IH12+(2*NR+1)*NBLK+NR ), NBLK )
                  END IF
               END IF
               IF( K.GT.NR+P ) THEN
                  CALL DLACPY( 'Full', NR, -Q-P, DWORK( IB ), N,
     $                         DWORK( IH12+(2*NR+1+P)*NBLK ), NBLK )
                  CALL DLACPY( 'Full', P, -Q-P, D, LDD,
     $                         DWORK( IH12+(2*NR+1+P)*NBLK+2*NR ),
     $                         NBLK )
               END IF
C
               IF( FULLE ) THEN
                  CALL DLACPY( 'Full', NE, MIN( K, NE ), DWORK( IE ), N,
     $                         DWORK( IJ12+(NR+1)*NBLK ), NBLK )
               ELSE
                  TEMP( 1 ) = ONE
                  CALL DCOPY( MIN( NR, K ), TEMP, 0,
     $                        DWORK( IJ12+(NR+1)*NBLK ), NBLK+1 )
               END IF
            END IF
C
C           Compute the generalized eigenvalues using the structure-
C           preserving method for skew-Hamiltonian/Hamiltonian pencils.
C
C           Additional workspace: need    7*NBLK*NBLK + MAX(L,36), where
C                                   L = 8*NBLK + 4, if NBLK is even, and
C                                   L = 8*NBLK,     if NBLK is odd. 
C                                 (from IT)
C                                 prefer larger.
C           Integer    workspace: need   NBLK + 12.
C
            LIW = MAX( NBLK+12, 4*NBLK+3 )
            CALL MB04BD( 'Eigenvalues', 'No Computation',
     $                   'No Computation', 2*NBLK, DWORK( IJ ), NBLK,
     $                   DWORK( IJ12 ), NBLK, DWORK( IH ), NBLK,
     $                   DWORK( IH12 ), NBLK, DWORK, 1, DWORK, 1,
     $                   DWORK( IT ), NBLK, DWORK( IT12 ), NBLK,
     $                   DWORK( IH22 ), NBLK, DWORK( IR ), DWORK( II ),
     $                   DWORK( IBT ), IWORK, LIW, DWORK( IWRK ),
     $                   LDWORK-IWRK+1, IERR )
            IF( IERR.EQ.1 .OR. IERR.EQ.2 ) THEN
               INFO = 2
               RETURN
            END IF
            MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
         END IF
C
C        Detect eigenvalues on the boundary of the stability domain,
C        if any. The test is based on a round-off level of eps*rho(H)
C        (after balancing) resulting in worst-case perturbations of
C        order sqrt(eps*rho(H)), for continuous-time systems, on the
C        real part of poles of multiplicity two (typical as GAMMA
C        approaches the infinity norm). Similarly, in the discrete-time
C        case. Above, rho(H) is the maximum modulus of eigenvalues
C        (continuous-time case).
C
C        Compute maximum eigenvalue modulus and check the absolute real
C        parts (if DICO = 'C'), or moduli (if DICO = 'D').
C
         WMAX = ZERO
C
C        Additional workspace: need   NBLK, if DICO = 'D';
C                                           (from IBT+NBLK)
C                                     0,    if DICO = 'C'.
C
         IF ( USEPEN ) THEN
         TOLE = TEN*DBLE( NBLK*SDIM )*EPS
         DO 190 I = 0, NBLK - 1
C
C              The pencil has infinite eigenvalues. ZERO is used.
C
               TM = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
               IF ( ( DWORK( IBT+I ).GE.ONE ) .OR.
     $              ( DWORK( IBT+I ).LT.ONE  .AND.
     $               TM.LT.DWORK( IBT+I )/TOLE ) ) THEN
                  TM = TM / DWORK( IBT+I )
               ELSE
C
C                 The pencil has too large eigenvalues. ZERO is used.
C
                  TM = ZERO
               END IF
               WMAX = MAX( WMAX, TM )
  190    CONTINUE
         ELSE
C
            DO 270 I = 0, NR - 1
               TM   = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
               WMAX = MAX( WMAX, TM )
  270       CONTINUE
C
         END IF
C
         NEI = 0
C
         IF ( USEPEN ) THEN
            DO 200 I = 0, NBLK - 1
C
C                    The pencil has too large eigenvalues.
C                    SAFMAX is used.
C
               TM  = ABS( DWORK( IR+I ) )
               TMP = DLAPY2( DWORK( IR+I ), DWORK( II+I ) )
               IF ( TMP.LT.DWORK( IBT+I )/TOLE ) THEN
                  TM = TM / DWORK( IBT+I )
               ELSE
C
C                 The pencil has infinite eigenvalues. SAFMAX is used.
C
                  TM = SAFMAX
               END IF
               IF ( TM.LE.TOLER*SQRT( HUNDRD + WMAX ) ) THEN
                  DWORK( IR+NEI ) = DWORK( IR+I ) / DWORK( IBT+I )
                  DWORK( II+NEI ) = DWORK( II+I ) / DWORK( IBT+I )
                  NEI = NEI + 1
               END IF
  200       CONTINUE
         ELSE
C
            DO 290 I = 0, NR - 1
               TM = ABS( DWORK( IR+I ) )
               IF ( TM.LE.TOLER*SQRT( HUNDRD + WMAX ) ) THEN
                  DWORK( IR+NEI ) = DWORK( IR+I )
                  DWORK( II+NEI ) = DWORK( II+I )
                  NEI = NEI + 1
               END IF
  290       CONTINUE
C
         END IF
C
         IF( NEI.EQ.0 ) THEN
C
C           There is no eigenvalue on the boundary of the stability
C           domain for G = ( ONE + TOLN )*GAMMAL. The norm was found.
C
            GPEAK( 1 ) = GAMMAL
            GPEAK( 2 ) = ONE
            GO TO 240
         END IF
C
C        Compute the frequencies where the gain G is attained and
C        generate new test frequencies.
C
         NWS = 0
C
         IF ( DISCR ) THEN
C
            DO 210 I = 0, NEI - 1
C
C              Back transformation of eigenvalues. Switch to complex
C              arithmetic to avoid over-/underflow.
C
               CTMP1 = DCMPLX( ONE+DWORK( II+I ), ZERO ) /
     $                 DCMPLX( ONE, DWORK( II+I ) )
               CTMP2 = DCMPLX( ONE-DWORK( II+I ), ZERO ) /
     $                 DCMPLX( ONE, -DWORK( II+I ) )
               DWORK( IR+I ) = DBLE( CTMP1*CTMP2 )
               CTMP1 = DCMPLX( 2*DWORK( II+I ), ZERO ) / 
     $                 DCMPLX( ONE, DWORK( II+I ) )
               CTMP2 = DCMPLX( ONE, ZERO ) / 
     $                 DCMPLX( ONE, -DWORK( II+I ) )
               DWORK( II+I ) = DBLE( CTMP1*CTMP2 )
C
               TM = ATAN2( DWORK( II+I ), DWORK( IR+I ) )
               DWORK( IR+I ) = MAX( EPS, TM )
               NWS = NWS + 1
  210       CONTINUE
C
         ELSE
C
            J = 0
C
            DO 220 I = 0, NEI - 1
               IF ( DWORK( II+I ).GT.EPS ) THEN
                  DWORK( IR+NWS ) = DWORK( II+I )
                  NWS = NWS + 1
               ELSE IF ( DWORK( II+I ).EQ.EPS ) THEN
                  J = J + 1
                  IF ( J.EQ.1 ) THEN
                     DWORK( IR+NWS ) = EPS
                     NWS = NWS + 1
                  END IF
               END IF
  220       CONTINUE
C
         END IF
C
         CALL DLASRT( 'Increasing', NWS, DWORK( IR ), IERR )
         LW = 1
C
         DO 230 I = 0, NWS - 1
            IF ( DWORK( IR+LW-1 ).NE.DWORK( IR+I ) ) THEN
               DWORK( IR+LW ) = DWORK( IR+I )
               LW = LW + 1
            END IF
  230    CONTINUE
C
         IF ( LW.EQ.1 ) THEN
            IF ( ITER.EQ.1 .AND. NWS.GE.1 ) THEN
C
C              Duplicate the frequency trying to force iteration.
C
               DWORK( IR+1 ) = DWORK( IR )
               LW = LW + 1
            ELSE
C
C              The norm was found.
C
               GPEAK( 1 ) = GAMMAL
               GPEAK( 2 ) = ONE
               GO TO 240
            END IF
         END IF
C
C        Form the vector of mid-points and compute the gain at new test
C        frequencies. Save the current lower bound.
C
         IWRK   = IR + LW
         GAMMAS = GAMMAL
C
         DO 250 I = 0, LW - 2
            IF ( DISCR ) THEN
               OMEGA = ( DWORK( IR+I ) + DWORK( IR+I+1 ) ) / TWO
            ELSE
               OMEGA = SQRT( DWORK( IR+I )*DWORK( IR+I+1 ) )
            END IF
C
C           Additional workspace:  need   LDW2, see above (from IWRK);
C                                  prefer larger.
C
            GAMMA = AB13DX( DICO, JOB, JOBD, NR, M, P, OMEGA,
     $                      DWORK( IA ), N, DWORK( IE ), N, DWORK( IB ),
     $                      N, DWORK( IC ), P, D, LDD, IWORK,
     $                      DWORK( IWRK ), LDWORK-IWRK+1, ZWORK, LZWORK,
     $                      IERR )
            MAXCWK = MAX( INT( ZWORK( 1 ) ), MAXCWK )
            MAXWRK = MAX( INT( DWORK( IWRK ) ) + IWRK - 1, MAXWRK )
            IF ( DISCR ) THEN
               TM = ABS( ATAN2( SIN( OMEGA ), COS( OMEGA ) ) )
            ELSE
               TM = OMEGA
            END IF
            IF( IERR.GE.1 .AND. IERR.LE.NR ) THEN
               GPEAK( 1 ) = ONE
               FPEAK( 1 ) = TM
               GPEAK( 2 ) = ZERO
               FPEAK( 2 ) = ONE
               GO TO 240
            ELSE IF( IERR.EQ.NR+1 ) THEN
               INFO = 3
               RETURN
            END IF
C
            IF( GAMMAL.LT.GAMMA ) THEN
               GAMMAL     = GAMMA
               FPEAK( 1 ) = TM
               FPEAK( 2 ) = ONE
            END IF
  250    CONTINUE
C
C        If the lower bound has not been improved, return. (This is a
C        safeguard against undetected modes of Hamiltonian matrix on the
C        boundary of the stability domain.)
C
         IF ( GAMMAL.LT.GAMMAS*( ONE + TOLN/TEN ) ) THEN
            GPEAK( 1 ) = GAMMAL
            GPEAK( 2 ) = ONE
            GO TO 240
         END IF
C
C     END WHILE
C
      IF ( ITER.LE.MAXIT ) THEN
         GO TO 120
      ELSE
         INFO = 4
         RETURN
      END IF
C
  240 CONTINUE
      DWORK( 1 ) = MAXWRK
      ZWORK( 1 ) = MAXCWK
C
      RETURN
C *** Last line of AB13HD ***
      END
