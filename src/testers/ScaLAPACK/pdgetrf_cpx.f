      SUBROUTINE PDGETRF_CPX(M, N,
     $                      A, IA, JA, DESCA,
     $                      ACP, IACP, JACP, DESCACP,
     $                      IPIV, IPIVCP, MIPIV,
     $                      CP_INTERVAL, FAULTNUM, FAULTLIST, JFAULT,
     $                      RCVR, ICTXTALL, CALCPROCS, INFO)
*
*       Modified to perform X>1 checkpoints to different X processes
*           (restore not working)
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
*
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, M, N, MIPIV, ICTXTALL,
     $                   CALCPROCS, EL, FAULTNUM,
     $                   IACP, JACP, JFAULT, CP_INTERVAL, RCVR
*     ..
*     .. Array Arguments ..
      INTEGER          DESCA( * ), IPIV( * ),
     $                 DESCACP( 9,* ), IPIVCP( * ), FAULTLIST( * )
      DOUBLE PRECISION A( * ), ACP( * )
*     ..
*
*  Purpose
*  =======
*
*  PDGETRF computes an LU factorization of a general M-by-N distributed
*  matrix sub( A ) = (IA:IA+M-1,JA:JA+N-1) using partial pivoting with
*  row interchanges.
*
*  The factorization has the form sub( A ) = P * L * U, where P is a
*  permutation matrix, L is lower triangular with unit diagonal ele-
*  ments (lower trapezoidal if m > n), and U is upper triangular
*  (upper trapezoidal if m < n). L and U are stored in sub( A ).
*
*  This is the right-looking Parallel Level 3 BLAS version of the
*  algorithm.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  This routine requires square block decomposition ( MB_A = NB_A ).
*
*  Arguments
*  =========
*
*  M       (global input) INTEGER
*          The number of rows to be operated on, i.e. the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) DOUBLE PRECISION pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the M-by-N
*          distributed matrix sub( A ) to be factored. On exit, this
*          array contains the local pieces of the factors L and U from
*          the factorization sub( A ) = P*L*U; the unit diagonal ele-
*          ments of L are not stored.
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  IPIV    (local output) INTEGER array, dimension ( LOCr(M_A)+MB_A )
*          This array contains the pivoting information.
*          IPIV(i) -> The global row local row i was swapped with.
*          This array is tied to the distributed matrix A.
*
*  INFO    (global output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  If INFO = K, U(IA+K-1,JA+K-1) is exactly zero.
*                The factorization has been completed, but the factor U
*                is exactly singular, and division by zero will occur if
*                it is used to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          COLBTOP, COLCTOP, ROWBTOP
      INTEGER            I, ICOFF, ICTXT, IINFO, IN, IROFF, J, JB, JN,
     $                   MN, MYCOL, MYROW, NPCOL, NPROW, P,
     $                   NPROCS, MYPNUM, IERR,
     $                   CP_COUNTDOWN, CP_DONE, CP_AFTER_FAULT, LASTCP,
     $                   FAULT_OCCURRED, FAULT_RECOVERABLE,
     $                   COMM_SPARE, COMM_CALC_FIRST_SPARE,
     $                   GROUPING, FIRSTSPAREPNUM
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 1 ), IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, CHK1MAT, IGAMN2D, PCHK1MAT,
     $                   PB_TOPGET, PB_TOPSET, PDGEMM, PDGETF2,
     $                   PDLASWP, PDTRSM, PXERBLA,
     $                   BLACS_PINFO, BLACS_BARRIER
*     ..
*     .. External Functions ..
      INTEGER            ICEIL
      EXTERNAL           ICEIL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. my Executable Statements ..
*
      INCLUDE "mpif.h"
      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, MYPNUM, IERR )
      CALL BLACS_PINFO ( MYPNUM, NPROCS )

*     mpi rank of the first spare process
      FIRSTSPAREPNUM=CALCPROCS

*     make proper communicators:
*       all calc procs with the first spare
      IF ( MYPNUM.LE.FIRSTSPAREPNUM ) THEN
        GROUPING=1
      ELSE
        GROUPING=0
      END IF
      CALL MPI_COMM_SPLIT ( MPI_COMM_WORLD, GROUPING, MYPNUM,
     $                       COMM_CALC_FIRST_SPARE, IERR )
*       only the spare procs
      IF ( MYPNUM.GE.FIRSTSPAREPNUM ) THEN
        GROUPING=1
      ELSE
        GROUPING=0
      END IF
      CALL MPI_COMM_SPLIT ( MPI_COMM_WORLD, GROUPING, MYPNUM,
     $                       COMM_SPARE, IERR )
*
*     init counters for checkpointing
      IF ( CP_INTERVAL.LT.0 ) THEN
*       checkpointing disabled
        CP_COUNTDOWN=-1
      ELSE
*       checkpointing enabled
        CP_COUNTDOWN = CP_INTERVAL
      END IF
      CP_DONE=0
*
      FAULT_RECOVERABLE=1
      FAULT_OCCURRED=0
*      continue (-1) checkpointing after first fault or not (1)
*      CP_AFTER_FAULT=1
      CP_AFTER_FAULT=-1
*
      IF ( MYPNUM.GE.FIRSTSPAREPNUM ) THEN
*        PRINT*, "I'm a spare proc, doing"
        MN = MIN( M, N )
        IN = MIN( ICEIL( IA, DESCA( MB_ ) )*DESCA( MB_ ), IA+M-1 )
        JN = MIN( ICEIL( JA, DESCA( NB_ ) )*DESCA( NB_ ), JA+MN-1 )
        JB = JN - JA + 1
*        do while in F77
        J = JN+1
  102   IF (J.LE.JA+MN-1) THEN
          JB = MIN( MN-J+JA, DESCA( NB_ ) )
          I = IA + J - JA
          IF ( (FAULT_OCCURRED.NE.CP_AFTER_FAULT)
     $                   .AND.(CP_COUNTDOWN.EQ.0) ) THEN
            CALL BLACS_BARRIER ( ICTXTALL, 'A' )
*            do checkpointing
*            iterate on spare procs (one checkpoint per spare proc)
            P=1
  202       IF ( MYPNUM.EQ.(P+CALCPROCS-1) ) THEN
*              PRINT '(A16,I0.5,A9,I0.5,A6,I0.5,A21)',
*     $              "## pgetrf_cpx: (", MYPNUM, ") @ iter ",
*     $              J," with ", FAULT_OCCURRED," fault(s): checkpoint"
            END IF
*            checkpoint matrix
            CALL PDGEMR2D(M, N,
     $                    A, IA, JA, DESCA,
     $                    ACP, IACP, JACP, DESCACP(:,P),
     $                    ICTXTALL)
*            barrier to flush the communication buffer
            CALL BLACS_BARRIER ( ICTXTALL, 'A' )
            P=P+1
            IF ( P.LE.(NPROCS-CALCPROCS) ) THEN
              GO TO 202
            END IF
*            checkpoint work space
            IF ( MYPNUM.EQ.FIRSTSPAREPNUM ) THEN
*              checkpoint on first spare proc
              CALL MPI_GATHER(IPIV, MIPIV, MPI_INTEGER, IPIVCP,
     $                        MIPIV, MPI_INTEGER, FIRSTSPAREPNUM,
     $                        COMM_CALC_FIRST_SPARE, IERR)
*              and broadcasts to the other spare procs
              CALL MPI_BCAST(IPIVCP, MIPIV*CALCPROCS, MPI_INTEGER,
     $                       FIRSTSPAREPNUM-CALCPROCS, COMM_SPARE, IERR)
            ELSE
              CALL MPI_BCAST(IPIVCP, MIPIV*CALCPROCS, MPI_INTEGER,
     $                       FIRSTSPAREPNUM-CALCPROCS, COMM_SPARE, IERR)
            END IF
*            save checkpoint instant
            LASTCP=J
*            reset counter for checkpointing
            CP_COUNTDOWN=CP_INTERVAL
            CP_DONE=1
            CALL BLACS_BARRIER ( ICTXTALL, 'A' )
          ELSE
*            skip checkpointing
*            PRINT*, J,"iter, with",FAULT_OCCURRED,"faults"
*            decrease counter for checkpointing
            CP_COUNTDOWN=CP_COUNTDOWN-1
          END IF
*
          IF ( FAULTNUM.GT.0 ) THEN
            IF ( FAULT_OCCURRED.EQ.0 ) THEN
              IF ( (J.LE.JFAULT).AND.(JFAULT.LT.(J+DESCA(NB_))) ) THEN
*                fault!
*                PRINT*, "##(spare) pgetrf_cp:",
*     $                  J,"iter, fault! at",JFAULT,"<",J+DESCA(NB_)
*                PRINT '(A16,I0.5,A9,I0.5,A6,I0.5,A21)',
*     $                "## pgetrf_cpx: (", MYPNUM, ") @ iter ",
*     $                J," with ", FAULT_OCCURRED," fault(s): fault"
                FAULT_OCCURRED=1
                IF ( RCVR.NE.0 ) THEN
                  IF ( CP_DONE.EQ.0 ) THEN
*                     can't recover, exit
                    J=JA+MN
*                    PRINT*, "## pgetrf_cp: ..unrecoverable! exiting.."
                    FAULT_RECOVERABLE=0
                    CALL BLACS_BARRIER ( ICTXTALL, 'A' )
                  ELSE
*                    PRINT*, "## pgetrf_cp: recovering.."
                    CALL BLACS_BARRIER ( ICTXTALL, 'A' )
                    CALL PDGEMR2D ( M, N, ACP, IACP, JACP, DESCACP,
     $                              A, IA, JA, DESCA, ICTXTALL )
                    CALL MPI_SCATTER(IPIVCP, MIPIV, MPI_INTEGER, IPIV,
     $                               MIPIV, MPI_INTEGER, NPROCS-1,
     $                               MPI_COMM_WORLD, IERR)
                    J=LASTCP
                    CP_COUNTDOWN=CP_INTERVAL
                    CP_DONE=0
                    CALL BLACS_BARRIER ( ICTXTALL, 'A' )
*                    PRINT*, "## pgetrf_cp: ..recovered"
                  END IF
                END IF
              END IF
            END IF
          END IF
*
          J=J+DESCA( NB_ )
*
          GO TO 102
        END IF
*        end do while
*
      ELSE
*     .. Executable Statements ..
*
*     Get grid parameters
*
        ICTXT = DESCA( CTXT_ )
        CALL BLACS_GRIDINFO ( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*       Test the input parameters
*
        INFO = 0
        IF ( NPROW.EQ.-1 ) THEN
          INFO = -(600+CTXT_)
        ELSE
          CALL CHK1MAT ( M, 1, N, 2, IA, JA, DESCA, 6, INFO )
          IF ( INFO.EQ.0 ) THEN
            IROFF = MOD( IA-1, DESCA( MB_ ) )
            ICOFF = MOD( JA-1, DESCA( NB_ ) )
            IF ( IROFF.NE.0 ) THEN
               INFO = -4
            ELSE IF ( ICOFF.NE.0 ) THEN
               INFO = -5
            ELSE IF ( DESCA( MB_ ).NE.DESCA( NB_ ) ) THEN
               INFO = -(600+NB_)
            END IF
          END IF
          CALL PCHK1MAT ( M, 1, N, 2, IA, JA, DESCA, 6, 0, IDUM1,
     $                  IDUM2, INFO )
        END IF
*
        IF ( INFO.NE.0 ) THEN
          CALL PXERBLA ( ICTXT, 'PDGETRF', -INFO )
          RETURN
        END IF
*
*       Quick return if possible
*
        IF ( DESCA( M_ ).EQ.1 ) THEN
          IPIV( 1 ) = 1
          RETURN
        ELSE IF ( M.EQ.0 .OR. N.EQ.0 ) THEN
          RETURN
        END IF
*
*       Split-ring topology for the communication along process rows
*
        CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
        CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
        CALL PB_TOPGET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
        CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', 'S-ring' )
        CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', ' ' )
        CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', ' ' )
*
*       Handle the first block of columns separately
*
        MN = MIN( M, N )
        IN = MIN( ICEIL( IA, DESCA( MB_ ) )*DESCA( MB_ ), IA+M-1 )
        JN = MIN( ICEIL( JA, DESCA( NB_ ) )*DESCA( NB_ ), JA+MN-1 )
        JB = JN - JA + 1
*
*       Factor diagonal and subdiagonal blocks and test for exact
*       singularity.
*
        CALL PDGETF2( M, JB, A, IA, JA, DESCA, IPIV, INFO )
*
        IF ( JB+1.LE.N ) THEN
*
*         Apply interchanges to columns JN+1:JA+N-1.
*
          CALL PDLASWP( 'Forward', 'Rows', N-JB, A, IA, JN+1, DESCA,
     $                 IA, IN, IPIV )
*
*         Compute block row of U.
*
          CALL PDTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                N-JB, ONE, A, IA, JA, DESCA, A, IA, JN+1, DESCA )
*
          IF ( JB+1.LE.M ) THEN
*
*           Update trailing submatrix.
*
            CALL PDGEMM( 'No transpose', 'No transpose', M-JB, N-JB, JB,
     $                  -ONE, A, IN+1, JA, DESCA, A, IA, JN+1, DESCA,
     $                   ONE, A, IN+1, JN+1, DESCA )
*
          END IF
        END IF
*
*       Loop over the remaining blocks of columns.
*
*       do while in F77
        J = JN+1
  100   IF ( J.LE.JA+MN-1 ) THEN
          JB = MIN( MN-J+JA, DESCA( NB_ ) )
          I = IA + J - JA
*
*         Factor diagonal and subdiagonal blocks and test for exact
*         singularity.
*
          CALL PDGETF2( M-J+JA, JB, A, I, J, DESCA, IPIV, IINFO )
*
          IF ( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + J - JA
*
*         Apply interchanges to columns JA:J-JA.
*
          CALL PDLASWP( 'Forward', 'Rowwise', J-JA, A, IA, JA, DESCA,
     $                 I, I+JB-1, IPIV )
*
          IF ( J-JA+JB+1.LE.N ) THEN
*
*           Apply interchanges to columns J+JB:JA+N-1.
*
            CALL PDLASWP( 'Forward', 'Rowwise', N-J-JB+JA, A, IA, J+JB,
     $                   DESCA, I, I+JB-1, IPIV )
*
*           Compute block row of U.
*
            CALL PDTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                  N-J-JB+JA, ONE, A, I, J, DESCA, A, I, J+JB,
     $                  DESCA )
*
            IF ( J-JA+JB+1.LE.M ) THEN
*
*             Update trailing submatrix.
*
              CALL PDGEMM( 'No transpose', 'No transpose', M-J-JB+JA,
     $                     N-J-JB+JA, JB, -ONE, A, I+JB, J, DESCA, A,
     $                     I, J+JB, DESCA, ONE, A, I+JB, J+JB, DESCA )
*
            END IF
          END IF
*
          IF ( (FAULT_OCCURRED.NE.CP_AFTER_FAULT)
     $         .AND.(CP_COUNTDOWN.EQ.0) )           THEN
            CALL BLACS_BARRIER ( ICTXTALL, 'A' )
*            do checkpointing
            IF ( MYPNUM.EQ.0 ) THEN
              PRINT '(A31,I0.5,A14)',
     $             "## pgetrf_cpx: ( all ) @ iter ", J, "checkpointing"
            END IF
*            iterate on spare procs (one checkpoint per spare proc)
            P=1
*            checkpoint matrix
  200       CALL PDGEMR2D(M, N,
     $                    A, IA, JA, DESCA,
     $                    ACP, IACP, JACP, DESCACP(:,P),
     $                    ICTXTALL)
*            barrier to flush the communication buffer
            CALL BLACS_BARRIER ( ICTXTALL, 'A' )
            P=P+1
            IF ( P.LE.(NPROCS-CALCPROCS) ) THEN
              GO TO 200
            END IF
*            checkpoint work space (to first spare, only)
            CALL MPI_GATHER(IPIV, MIPIV, MPI_INTEGER, IPIVCP,
     $                      MIPIV, MPI_INTEGER, FIRSTSPAREPNUM,
     $                      COMM_CALC_FIRST_SPARE, IERR)
*            save checkpoint instant
            LASTCP=J
*            reset counter for checkpointing
            CP_COUNTDOWN=CP_INTERVAL
            CP_DONE=1
            CALL BLACS_BARRIER ( ICTXTALL, 'A' )
          ELSE
*            skip checkpointing
*            decrease counter for checkpointing
            CP_COUNTDOWN=CP_COUNTDOWN-1
          END IF
*
          IF ( FAULTNUM.GT.0 ) THEN
           IF ( FAULT_OCCURRED.EQ.0 ) THEN
            IF ( (J.LE.JFAULT).AND.(JFAULT.LT.(J+DESCA(NB_))) ) THEN
*              fault!
              FAULT_OCCURRED=1

                DO 400 P = 1, FAULTNUM
                  IF ( MYPNUM.EQ.FAULTLIST(P) ) THEN
                   PRINT '(A17,I0.5,A9,I0.5,A12,I0.5,A14,I0.5,A3,I0.5)',
     $                    "## pgetrf_cpx: (", MYPNUM,
     $                    ") @ iter ", J,
     $                    " with fault ", FAULT_OCCURRED,
     $                    " occurred at ", JFAULT," < ",J+DESCA(NB_)
                    LOCM = NUMROC(DESCA( M_ ), DESCA( MB_ ),
     $                            MYROW, DESCA( RSRC_ ), NPROW)
                    LOCN = NUMROC(DESCA( N_ ), DESCA( NB_ ),
     $                            MYROW, DESCA( CSRC_ ), NPCOL)
                    DO 300 EL = 1, LOCM*LOCN
                      A(EL)=-99.0
300                 CONTINUE
                  END IF
400             CONTINUE
*              END IF
*
              IF ( RCVR.NE.0 ) THEN
                IF ( CP_DONE.EQ.0 ) THEN
*                  can't recover, exit
                  J=JA+MN
                  IF ( MYPNUM.EQ.0 ) THEN
                    PRINT*, " ## ..unrecoverable! exiting.."
                  END IF
                  FAULT_RECOVERABLE=0
                  CALL BLACS_BARRIER ( ICTXTALL, 'A' )
                ELSE
*                  IF (MYPNUM.EQ.0) THEN
*                    PRINT*, " ## recovering.."
*                  END IF
                  CALL BLACS_BARRIER ( ICTXTALL, 'A' )
                  CALL PDGEMR2D(M, N, ACP, IACP, JACP, DESCACP,
     $                          A, IA, JA, DESCA, ICTXTALL)
                  CALL MPI_SCATTER(IPIVCP, MIPIV, MPI_INTEGER, IPIV,
     $                             MIPIV, MPI_INTEGER, NPROCS-1,
     $                             MPI_COMM_WORLD, IERR)
                  IF ( FAULTNUM.GT.0 ) THEN
                    DO 500 P = 1, FAULTNUM
                      IF ( MYPNUM.EQ.FAULTLIST(P) ) THEN
                        PRINT '(A17,I0.5,A9,I0.5,A12,I0.5,A14,I0.5)',
     $                        "## pgetrf_cpx: (", MYPNUM, ") @ iter ",
     $                        J," with fault ", FAULT_OCCURRED,
     $                        " recovered to ", LASTCP
                      END IF
500                 CONTINUE
                  J=LASTCP
                  CP_COUNTDOWN=CP_INTERVAL
                  CP_DONE=0
                  CALL BLACS_BARRIER ( ICTXTALL, 'A' )
                  END IF
                END IF
              END IF
            END IF
          END IF
         END IF
*
          J=J+DESCA( NB_ )
*
          GO TO 100
        END IF
*        end do while
*
*
        IF ( INFO.EQ.0 ) INFO = MN + 1
*
        CALL IGAMN2D( ICTXT, 'Rowwise', ' ', 1, 1, INFO, 1,
     $               IDUM1, IDUM2, -1, -1, MYCOL )
*
        IF( INFO.EQ.MN+1 ) INFO = 0
*
        CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP )
        CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP )
        CALL PB_TOPSET( ICTXT, 'Combine', 'Columnwise', COLCTOP )
*
*       if not recoverable, set error code
        IF ( FAULT_RECOVERABLE.EQ.0 ) THEN
          INFO = -99
        END IF
*
*     end of main if for CP
      END IF
*
      RETURN
*
*     End of PDGETRF
*
      END
