      MODULE UPCGMODULE
        DOUBLEPRECISION, PARAMETER :: PI = 3.141592653589793
        DOUBLEPRECISION, DIMENSION(5) :: SCAL = 
     2    (/ 1.0D0, 0.1D0, 0.1D0, 0.1D0, 1.0D0 /)
        TYPE TGLSPOLY
          INTEGER :: NDEGREE
          INTEGER :: NLANSTEP
          INTEGER :: NLAN2
          INTEGER :: NINTV = 1
          INTEGER :: IEIGCALC
          DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: D_LANCZOS
          DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: ALPHA
          DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: BETA
          DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: GAMMA
          DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: D_V1
          DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: D_V0
          DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: D_V
          DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: D_E
          DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: INTV
          DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: P0
          DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: PPOL
          DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: APPOL
          DOUBLEPRECISION, DIMENSION(:), ALLOCATABLE :: QPOL
        END TYPE TGLSPOLY
        INTEGER,SAVE,POINTER  :: ITER1C,NPC,NOPT,NTRD,NTRDV
        INTEGER,SAVE,POINTER  :: NITERC,NNZC,NIAC
        INTEGER,SAVE,POINTER  :: NIAPC,NIWC,NPOL,NEIG
        REAL   ,SAVE,POINTER  :: HCLOSEUPCG,RCLOSEUPCG
        DOUBLE PRECISION, SAVE, POINTER :: UPCGTOTT, UPCGFMAT
        DOUBLE PRECISION, SAVE, POINTER :: UPCGPCUT, UPCGPCAT
        DOUBLE PRECISION, SAVE, POINTER :: UPCGDPT, UPCGMVT
        DOUBLE PRECISION, SAVE, POINTER :: UPCGAXPYT,UPCGVVPT,UPCGMISCT
        DOUBLE PRECISION, SAVE, POINTER :: UPCGGPUTT
        INTEGER,SAVE,POINTER  :: IUPCGO,IUPCGI
        INTEGER,          SAVE, POINTER, DIMENSION(:,:,:) :: NODEC
        DOUBLE PRECISION, SAVE, POINTER, DIMENSION(:)     :: BC
        DOUBLE PRECISION, SAVE, POINTER, DIMENSION(:)     :: XC
        DOUBLE PRECISION, SAVE, POINTER, DIMENSION(:)     :: AC
        DOUBLE PRECISION, SAVE, POINTER, DIMENSION(:)     :: APC
        INTEGER,          SAVE, POINTER, DIMENSION(:)     :: IAC
        INTEGER,          SAVE, POINTER, DIMENSION(:)     :: JAC
        INTEGER,          SAVE, POINTER, DIMENSION(:)     :: IUC
        INTEGER,          SAVE, POINTER, DIMENSION(:)     :: IXMAP
C         WORKING ARRAYS        
        INTEGER,         SAVE, POINTER, DIMENSION(:)      :: IWC
        DOUBLEPRECISION, SAVE, POINTER, DIMENSION(:)      :: DC
        DOUBLEPRECISION, SAVE, POINTER, DIMENSION(:)      :: PC
        DOUBLEPRECISION, SAVE, POINTER, DIMENSION(:)      :: QC
        DOUBLEPRECISION, SAVE, POINTER, DIMENSION(:)      :: ZC
C         DIAGONAL SCALING VECTOR
        INTEGER,SAVE,POINTER  :: ISCL
        DOUBLEPRECISION, SAVE, POINTER, DIMENSION(:)      :: SCL
        DOUBLEPRECISION, SAVE, POINTER, DIMENSION(:)      :: SCLI
C         POLYNOMIAL PRECONDITIONER        
        TYPE (TGLSPOLY),SAVE, POINTER :: GLSPOLY
C         GPU POINTERS 
        INTEGER(KIND=8), SAVE, POINTER  :: CU_HDL,CU_STAT,CU_DES
        INTEGER(KIND=8), SAVE, POINTER  :: CU_IAC,CU_JAC
        INTEGER(KIND=8), SAVE, POINTER  :: CU_AC,CU_APC,CU_XC
        INTEGER(KIND=8), SAVE, POINTER  :: CU_DC,CU_ZC,CU_PC,CU_QC
        INTEGER(KIND=8), SAVE, POINTER  :: CU_SCL,CU_SCLI
        INTEGER(KIND=8), SAVE, POINTER  :: CU_V,CU_V0,CU_V1
        INTEGER(KIND=8), SAVE, POINTER  :: PL_DC,PL_ZC

      TYPE UPCGTYPE
        INTEGER,POINTER  :: ITER1C,NPC,NOPT,NTRD,NTRDV
        INTEGER,POINTER  :: NITERC,NNZC,NIAC
        INTEGER,POINTER  :: NIAPC,NIWC,NPOL,NEIG
        REAL   ,POINTER  :: HCLOSEUPCG,RCLOSEUPCG
        DOUBLE PRECISION, POINTER :: UPCGTOTT, UPCGFMAT
        DOUBLE PRECISION, POINTER :: UPCGPCUT, UPCGPCAT
        DOUBLE PRECISION, POINTER :: UPCGDPT, UPCGMVT
        DOUBLE PRECISION, POINTER :: UPCGAXPYT,UPCGVVPT,UPCGMISCT
        DOUBLE PRECISION, POINTER :: UPCGGPUTT
        INTEGER,POINTER  :: IUPCGO,IUPCGI
        INTEGER,          POINTER, DIMENSION(:,:,:) :: NODEC
        DOUBLE PRECISION, POINTER, DIMENSION(:)     :: BC
        DOUBLE PRECISION, POINTER, DIMENSION(:)     :: XC
        DOUBLE PRECISION, POINTER, DIMENSION(:)     :: AC
        DOUBLE PRECISION, POINTER, DIMENSION(:)     :: APC
        INTEGER,          POINTER, DIMENSION(:)     :: IAC
        INTEGER,          POINTER, DIMENSION(:)     :: JAC
        INTEGER,          POINTER, DIMENSION(:)     :: IUC
        INTEGER,          POINTER, DIMENSION(:)     :: IXMAP
C         WORKING ARRAYS        
        INTEGER,         POINTER, DIMENSION(:)      :: IWC
        DOUBLEPRECISION, POINTER, DIMENSION(:)      :: DC
        DOUBLEPRECISION, POINTER, DIMENSION(:)      :: PC
        DOUBLEPRECISION, POINTER, DIMENSION(:)      :: QC
        DOUBLEPRECISION, POINTER, DIMENSION(:)      :: ZC
C         DIAGONAL SCALING VECTOR
        INTEGER,POINTER  :: ISCL
        DOUBLEPRECISION, POINTER, DIMENSION(:)      :: SCL
        DOUBLEPRECISION, POINTER, DIMENSION(:)      :: SCLI
C         POLYNOMIAL PRECONDITIONER        
        TYPE (TGLSPOLY),POINTER :: GLSPOLY
C         GPU POINTERS 
        INTEGER(KIND=8),POINTER  :: CU_HDL,CU_STAT,CU_DES
        INTEGER(KIND=8),POINTER  :: CU_IAC,CU_JAC
        INTEGER(KIND=8),POINTER  :: CU_AC,CU_APC,CU_XC
        INTEGER(KIND=8),POINTER  :: CU_DC,CU_ZC,CU_PC,CU_QC
        INTEGER(KIND=8),POINTER  :: CU_SCL,CU_SCLI
        INTEGER(KIND=8),POINTER  :: CU_V,CU_V0,CU_V1
        INTEGER(KIND=8),POINTER  :: PL_DC,PL_ZC
      END TYPE
      TYPE(UPCGTYPE), SAVE ::UPCGDAT(10)
      DOUBLEPRECISION, POINTER :: DTDP
      DOUBLEPRECISION, POINTER :: DTMV
      DOUBLEPRECISION, POINTER :: DTAXPY
      DOUBLEPRECISION, POINTER :: DTVVP
      DOUBLEPRECISION, POINTER :: DTMISC
      END MODULE UPCGMODULE


      SUBROUTINE UPCG7AR(IN,MXITER,IGRID)
C     ******************************************************************
C     ALLOCATE STORAGE FOR UPCG ARRAYS AND READ UPCG DATA
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY,IBOUND
      USE UPCGMODULE
      USE OMP_LIB 
C
      CHARACTER*200 LINE
      INTEGER IN,MXITER
      INTEGER         :: ntc, ntdp, ntmv
      DOUBLEPRECISION :: tserialc,  tompc,  ttc,  suc
      DOUBLEPRECISION :: tserialdp, tompdp, ttdp, sudp
      DOUBLEPRECISION :: tserialmv, tompmv, ttmv, sumv
C
C-------FUNCTIONS      
      DOUBLEPRECISION :: SUPCGDP
C     ------------------------------------------------------------------
      ALLOCATE( ITER1C,NPC,NOPT,NTRD,NTRDV )
      ALLOCATE( NITERC,NNZC,NIAC,NIAPC,NIWC )
      ALLOCATE( NPOL,NEIG )
      ALLOCATE( HCLOSEUPCG, RCLOSEUPCG )
      ALLOCATE( UPCGTOTT, UPCGFMAT )
      ALLOCATE( UPCGPCUT, UPCGPCAT, UPCGDPT, UPCGMVT )
      ALLOCATE( UPCGAXPYT, UPCGVVPT, UPCGMISCT )
      ALLOCATE( UPCGGPUTT )

      ALLOCATE( IUPCGO, IUPCGI )
      ALLOCATE( GLSPOLY )
      ALLOCATE( DTDP, DTMV, DTAXPY, DTVVP, DTMISC )
      
      ALLOCATE( ISCL )

      ALLOCATE( NODEC(NCOL,NROW,NLAY) )
      NODESC = NCOL*NROW*NLAY
      NODEC = 0
      
      UPCGTOTT   = 0.0D0
      UPCGFMAT   = 0.0D0
      UPCGPCUT   = 0.0D0
      UPCGPCAT   = 0.0D0
      UPCGDPT    = 0.0D0
      UPCGMVT    = 0.0D0
      UPCGAXPYT  = 0.0D0
      UPCGVVPT   = 0.0D0
      UPCGMISCT  = 0.0D0
      UPCGGPUTT  = 0.0D0
      IUPCGO     = 0
      IUPCGI     = 0
      DTDP       = 0.0D0
      DTMV       = 0.0D0
      DTAXPY     = 0.0D0
      DTVVP      = 0.0D0
      DTMISC     = 0.0D0
      
      ISCL       = 0
C
C-------PRINT A MESSAGE IDENTIFYING UPCG PACKAGE
      WRITE (IOUT,500)
  500 FORMAT (1X,/1X,'UPCG -- UNSTRUCTURED CONJUGATE-GRADIENT SOLUTION',
     &        ' PACKAGE, VERSION 7.01, 02/09/2012',
     &        /1X,8X,'INCLUDES CPU, CPU-OPENMP, AND GPU-CUDA SUPPORT')
C
C-------READ AND PRINT COMMENTS, MXITER,ITER1 AND NPCOND
      CALL URDCOM(IN,IOUT,LINE)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXITER,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ITER1C,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPC,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NOPT,R,IOUT,IN)
      IF ( NPC.LT.0 ) THEN
        NPC  = ABS( NPC )
        ISCL = 1
      END IF
      IF ( NPC.LT.0 .OR. NPC.GT.4 ) THEN
          WRITE (IOUT,'(//,A)') 'UPCG7AR: NPC MUST BE >= 0 AND < 5'
          CALL USTOP('UPCG7AR: NPC MUST BE >= 0 AND < 5')
      END IF
      IF ( NOPT.LT.1 .OR. NOPT.GT.3 ) THEN
          WRITE (IOUT,'(//,A)') 'UPCG7AR: NOPT MUST BE > 0 AND < 4'
          CALL USTOP('UPCG7AR: NOPT MUST BE > 0 AND < 4')
      END IF
      IF ( NPC.EQ.4 ) THEN
        ISCL = 1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,N,R,IOUT,IN)
        GLSPOLY%NDEGREE  = N
        IF ( N.LT.1 ) THEN
          WRITE (IOUT,'(//,A)') 'UPCG7AR: POLYNOMIAL DEGREE MUST BE > 0'
          CALL USTOP('UPCG7AR: POLYNOMIAL DEGREE MUST BE > 0')
        END IF
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,N,R,IOUT,IN)
        GLSPOLY%NLANSTEP = N
        GLSPOLY%IEIGCALC = 1
        IF ( N.EQ.-2 ) THEN
          N = ABS( N )
          GLSPOLY%IEIGCALC = 0
          GLSPOLY%NLANSTEP = 0
        END IF
        IF ( N.LT.0 ) THEN
          WRITE (IOUT,'(//,A)') 'UPCG7AR: NLANSTEP MUST BE > 0 OR = -2'
          CALL USTOP('UPCG7AR: POLYNOMIAL NLANSTEP MUST BE > 0 OR = -2')
        END IF
      END IF
      NTRD  = 1
      NTRDV = 1
      IF ( NOPT.EQ.2 ) THEN
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NTRD,R,-IOUT,IN)
        IF ( NTRD.LT.0 ) THEN
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NTRDV,R,-IOUT,IN)
          IF ( NTRDV.LT.1 ) THEN
              WRITE (IOUT,'(//,A)') 'UPCG7AR: NTRDV MUST BE > 0'
              CALL USTOP('UPCG7AR: NTRDV MUST BE > 0')
          END IF
        END IF
      END IF
C
C-------READ HCLOSEPCG,RCLOSEPCG,RELAXPCG,NBPOL,IPRPCG,MUTPCG
      CALL URDCOM(IN,IOUT,LINE)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,HCLOSEUPCG,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,RCLOSEUPCG,IOUT,IN)
      IF ( HCLOSEUPCG.LE.0.0 ) THEN
          WRITE (IOUT,'(//,A)') 'UPCG7AR: HCLOSE MUST BE > 0.0'
          CALL USTOP('UPCG7AR: HCLOSE MUST BE > 0.0')
      END IF
      IF ( RCLOSEUPCG.LE.0.0 ) THEN
          WRITE (IOUT,'(//,A)') 'UPCG7AR: RCLOSE MUST BE > 0.0'
          CALL USTOP('UPCG7AR: RCLOSE MUST BE > 0.0')
      END IF
C
C-------PRINT MXITER,ITER1C,NPC,HCLOSEPCG,RCLOSEPCG,NPC,NOPT
C-------MUTPCG,DAMPPCG
        WRITE (IOUT,505)
  505   FORMAT (1X,/,18X,'SOLUTION BY THE CONJUGATE-GRADIENT METHOD',
     &        /,1X,75('-'))
      WRITE (IOUT,510) MXITER
  510 FORMAT (1X,1X,'MAXIMUM NUMBER OF CALLS TO PCG ROUTINE =',I9)
      WRITE (IOUT,515) ITER1C
  515 FORMAT (1X,5X,'MAXIMUM ITERATIONS PER CALL TO PCG =',I9)
      WRITE (IOUT,520) NPC
  520 FORMAT (1X,12X,'MATRIX PRECONDITIONING TYPE =',I9)
      WRITE (IOUT,535) HCLOSEUPCG
  535 FORMAT (1X,6X,'HEAD CHANGE CRITERION FOR CLOSURE =',E15.5)
      WRITE (IOUT,540) RCLOSEUPCG
  540 FORMAT (1X,2X,'RESIDUAL CHANGE CRITERION FOR CLOSURE =',E15.5)
C
      WRITE (IOUT,545) NPC, NOPT
  545 FORMAT (/1X,75('-'),/,
     &        19X,' MATRIX PRECONDITIONING TYPE :',I5,/,
     &        19X,'   NONE      : NPC = 0',/,
     &        19X,'   JACOBI    : NPC = 1',/,
     &        19X,'   ILU0      : NPC = 2',/,
     &        19X,'   MILU0     : NPC = 3',/,
     &        19X,'   GLS POLY. : NPC = 4',/,
     &        19X,' HARDWARE SOLUTION OPTION    :',I5,/,
     &        19X,'   CPU       : NOPT = 1',/,
     &        19X,'   CPU - OMP : NOPT = 2',/,
     &        19X,'   GPU - CUDA: NOPT = 3',/,
     &        1X,75('-'))
C
C-------INITIALIZE NITERC  
      NITERC = 0
C
C-------CALCULATE NUMBER OF NON-ZERO ENTRIES IN MODEL GRID
      NNZC  = 0
      NIAC  = 0
      NIAPC = 0
      NIWC  = 0
      NPOL  = 0
      NEIG  = 0
      IC    = 0
      ieq   = 0
      NRC = NROW * NCOL
      DO K = 1, NLAY
        DO I = 1, NROW
          DO J = 1, NCOL
            IC = IC + 1
            IF ( IBOUND(J,I,K).GT.0 ) THEN
              NIAC = NIAC + 1
              NNZC = NNZC + 1
              ieq  = ieq  + 1
              NODEC(J,I,K) = ieq
C               TOP FACE
              IF ( K.GT.1 ) THEN
                IF ( IBOUND(J,I,K-1).GT.0 ) NNZC = NNZC + 1
              END IF
C               UPPER FACE
              IF ( I.GT.1 ) THEN
                IF ( IBOUND(J,I-1,K).GT.0 ) NNZC = NNZC + 1
              END IF
C               LEFT FACE
              IF ( J.GT.1 ) THEN
                IF ( IBOUND(J-1,I,K).GT.0 ) NNZC = NNZC + 1
              END IF
C               RIGHT FACE
              IF ( J.LT.NCOL ) THEN
                IF ( IBOUND(J+1,I,K).GT.0 ) NNZC = NNZC + 1
              END IF
C               LOWER FACE
              IF ( I.LT.NROW ) THEN
                IF ( IBOUND(J,I+1,K).GT.0 ) NNZC = NNZC + 1
              END IF
C               BOTTOM FACE
              IF ( K.LT.NLAY ) THEN
                IF ( IBOUND(J,I,K+1).GT.0 ) NNZC = NNZC + 1
              END IF
            END IF
          END DO
        END DO
      END DO
C
C-------ALLOCATE AND INITIALIZE COMPRESSED ROW STORAGE VECTORS
C       COEFFICIENT MATRIX AND PRECONDITIONER MATRIX
      ALLOCATE(AC(NNZC))
      NIAPC = NNZC
      IF ( NPC.EQ.0 ) THEN
        NIAPC = 1
      ELSE IF ( NPC.EQ.1 ) THEN
        NIAPC = NIAC
      ELSE IF ( NPC.EQ.4 ) THEN
        NIAPC = 1
      END IF
      ALLOCATE(APC(NIAPC))
      ALLOCATE(IAC(NIAC+1),JAC(NNZC),IUC(NIAC),IXMAP(NIAC))
C       ALLOCATE WORKING VECTORS FOR UPCG SOLVER      
      ALLOCATE(BC(NIAC),XC(NIAC))
      ALLOCATE(DC(NIAC),PC(NIAC))
      ALLOCATE(QC(NIAC),ZC(NIAC))
C       INITIALIZE PCG WORKING ARRAYS
      DO n = 1, NNZC
        AC(n)  = 0.0D0
        JAC(n) = 0
      END DO
      DO n = 1, NIAPC
        APC(n) = 0.0D0
      END DO
      DO n = 1, NIAC+1
        IAC(n) = 0
      END DO
      DO n = 1, NIAC
        IUC(n)   = 0
        IXMAP(n) = 0
        BC(n)    = 0.0D0
        XC(n)    = 0.0D0
C         WORKING ARRAYS
        DC(n)    = 0.0D0
        PC(n)    = 0.0D0
        QC(n)    = 0.0D0
        ZC(n)    = 0.0D0
      END DO
C       ALLOCATE SPACE FOR ILU0 AND MILU0 NON-ZERO ROW ENTRY VECTOR
      IF ( NPC.EQ.2 .OR. NPC.EQ.3 ) THEN
        NIWC = NIAC
      ELSE
        NIWC = 1
      END IF
      ALLOCATE(IWC(NIWC))
C       INITIALIZE ILU0 AND MILU0 PRECONDITIONER WORKING VECTOR      
      DO n = 1, NIWC
        IWC(n)   = 0
      END DO
C       ALLOCATE DIAGONAL SCALING VECTOR
      ALLOCATE( SCL(NIAC)  )
      ALLOCATE( SCLI(NIAC) )
C       INITIALIZE DIAGONAL SCALING VECTOR
      DO n = 1, NIAC
        SCL(n)               = 1.0D0
        SCLI(n)              = 1.0D0
      END DO
C       ALLOCATE WORKING VECTORS FOR POLYNOMIAL PRECONDITIONER
      IF ( NPC.EQ.4 ) THEN
        NPOL = NIAC
        IF ( GLSPOLY%NLANSTEP.GT.NIAC ) THEN
          GLSPOLY%NLANSTEP = NIAC
        END IF
!        m    = 512 / 8
!        GLSPOLY%NLAN2 = ( NIAC + m - 1 ) / ( m * m )
        GLSPOLY%NLAN2 = NIAC
        NEIG = NPOL
        IF ( GLSPOLY%IEIGCALC.EQ.0 ) THEN
          NEIG             = 1
          GLSPOLY%NLANSTEP = 0
          GLSPOLY%NLAN2    = 1
        END IF
      ELSE
        NPOL             = 1
        NEIG             = 1
        GLSPOLY%NDEGREE  = 1
        GLSPOLY%NLANSTEP = 0
        GLSPOLY%NLAN2    = 1
      END IF
      nlansize = ( GLSPOLY%NLANSTEP + 1 ) * GLSPOLY%NLAN2
      ALLOCATE( GLSPOLY%D_LANCZOS(nlansize)                 )
      ALLOCATE( GLSPOLY%ALPHA(GLSPOLY%NDEGREE)              )
      ALLOCATE( GLSPOLY%BETA(GLSPOLY%NDEGREE+1)             )
      ALLOCATE( GLSPOLY%GAMMA(GLSPOLY%NDEGREE+1)            )
      ALLOCATE( GLSPOLY%D_V1(NPOL)                          )
      ALLOCATE( GLSPOLY%D_V0(NPOL)                          )
      ALLOCATE( GLSPOLY%D_V(NPOL)                           )
      ALLOCATE( GLSPOLY%D_E(NEIG)                           )
      ALLOCATE( GLSPOLY%INTV(2)                             )
      ALLOCATE( GLSPOLY%P0(2)                               )
      ALLOCATE( GLSPOLY%PPOL(GLSPOLY%NDEGREE+3)             )
      ALLOCATE( GLSPOLY%APPOL(GLSPOLY%NDEGREE+3)            )
      ALLOCATE( GLSPOLY%QPOL(GLSPOLY%NDEGREE+2)             )
C       INITIALIZE POLYNOMIAL PRECONDITIONER WORKING ARRAYS
      DO n = 1, nlansize
        GLSPOLY%D_LANCZOS(n) = 0.0D0
      END DO
      DO n = 1, NPOL
        GLSPOLY%D_V1(n)      = 0.0D0
        GLSPOLY%D_V0(n)      = 0.0D0
        GLSPOLY%D_V(n)       = 0.0D0
      END DO
      DO n = 1, NEIG
        GLSPOLY%D_E(n)       = 0.0D0
      END DO
      DO n = 1, 2
        IF ( GLSPOLY%IEIGCALC.NE.0 ) THEN
          GLSPOLY%INTV(n)      = 0.0D0
        ELSE
          IF ( n.EQ.1 ) THEN
            GLSPOLY%INTV(n) = 0.0D+00
          ELSE
            GLSPOLY%INTV(n) = 2.0D+00
          END IF
        END IF
        GLSPOLY%P0(n)       = 0.0D0
      END DO
      DO n = 1, GLSPOLY%NDEGREE+3
        IF ( n.LE.GLSPOLY%NDEGREE ) THEN
          GLSPOLY%ALPHA(n)   = 0.0D0
        END IF
        IF ( n.LE.GLSPOLY%NDEGREE+1 ) THEN
          GLSPOLY%BETA(n)    = 0.0D0
          GLSPOLY%GAMMA(n)   = 0.0D0
        END IF
        IF ( n.LE.GLSPOLY%NDEGREE+2 ) THEN
          GLSPOLY%QPOL(n)    = 0.0D0
        END IF
        GLSPOLY%PPOL(n)      = 0.0D0
        GLSPOLY%APPOL(n)     = 0.0D0
      END DO
C
C-------PRINT POLYNOMIAL PRECONDITIONER INFORMATION
      IF ( NPC.EQ.4 ) THEN
        WRITE (IOUT,550) GLSPOLY%NDEGREE
        IF ( GLSPOLY%IEIGCALC.NE.0 ) THEN
          WRITE (IOUT,555) GLSPOLY%NLANSTEP
        ELSE
          WRITE (IOUT,560)
        END IF
        WRITE (IOUT,565)
      END IF  
  550 FORMAT (/1X,75('-'),/,
     &        13X,' GENERAL LEAST-SQUARES POLYNOMIAL PRECONDITIONER',/,
     &        13X,'   DEGREE POLYNOMIAL = ',I5,/,
     &        13X,'   DIAGONAL RESCALING APPLIED')
  555 FORMAT (13X,'   LANCZOS STEPS     = ',I5,/,
     &        13X,'     MAY HAVE BEEN REDUCED TO PROBLEM SIZE')
  560 FORMAT (13X,'   MAX. AND MIN. EIGENVALUES ASSUMED TO BE',/,
     &        13X,'     2.0 AND 0.0, RESPECTIVELY.')
  565 FORMAT (1X,75('-'))
C
C-------FILL IA AND JA
      IND = 0
      IC  = 0
      ieq = 0
      DO K = 1, NLAY
        DO I = 1, NROW
          DO J = 1, NCOL
            IND = IND + 1
            IF ( IBOUND(J,I,K).GT.0 ) THEN
              IC = IC + 1
              ieq = ieq + 1
              IAC(ieq) = IC
              IXMAP(ieq) = IND
              JAC(IC) = NODEC(J,I,K)
C               TOP FACE
              IF ( K.GT.1 ) THEN
                IF ( IBOUND(J,I,K-1).GT.0 ) THEN
                  IC = IC + 1
                  JAC(IC) = NODEC(J,I,K-1)
                END IF
              END IF
C               UPPER FACE
              IF ( I.GT.1 ) THEN
                IF ( IBOUND(J,I-1,K).GT.0 ) THEN
                  IC = IC + 1
                  JAC(IC) = NODEC(J,I-1,K)
                END IF
              END IF
C               LEFT FACE
              IF ( J.GT.1 ) THEN
                IF ( IBOUND(J-1,I,K).GT.0 ) THEN
                  IC = IC + 1
                  JAC(IC)  = NODEC(J-1,I,K)
                END IF
              END IF
C               RIGHT FACE
              IF ( J.LT.NCOL ) THEN
                IF ( IBOUND(J+1,I,K).GT.0 ) THEN
                  IC = IC + 1
                  JAC(IC) = NODEC(J+1,I,K)
                END IF
              END IF
C               LOWER FACE
              IF ( I.LT.NROW ) THEN
                IF ( IBOUND(J,I+1,K).GT.0 ) THEN
                  IC = IC + 1
                  JAC(IC) = NODEC(J,I+1,K)
                END IF
              END IF
C               BOTTOM FACE
              IF ( K.LT.NLAY ) THEN
                IF ( IBOUND(J,I,K+1).GT.0 ) THEN
                  IC = IC + 1
                  JAC(IC) = NODEC(J,I,K+1)
                END IF
              END IF
              
            END IF
          END DO
        END DO
      END DO
C
C-------SET LAST POSITION IN IAC
      IAC(NIAC+1) = IC + 1
C
C-------SET IUC
      DO N = 1, NIAC
        DO I = IAC(N), IAC(N+1)-1
          IF ( JAC(I).GT.N .AND. IUC(N).EQ.0 ) IUC(N) = I
        END DO
C         SET IUC TO FIRST ELEMENT OF THE NEXT EQUATION IF DIAGONAL
C         IS LARGEST NODE NUMBER IN AC FOR CURRENT EQUATION              
        IF ( IUC(N).EQ.0 ) IUC(N) = IAC(N+1)
      END DO
C
C-------TIME SERIAL AND OPEN MP VECTOR OPERATIONS
      TSTOMP: IF ( NOPT.EQ.2 ) THEN
C         DETERMINE IF SPECIFIED NUMBER OF THREADS 
C         EXCEEDS THE MAXIMUM NUMBER OF THREADS AVAILABLE
C         MINUS ONE. IF SO, REDUCE TO THE MAXIMUM.
        NTRDMAX = OMP_GET_MAX_THREADS() - 1
        IF ( NTRD.NE.0 ) THEN
          IF ( ABS( NTRD ).GT.NTRDMAX ) THEN
            NTRD = ( NTRD / ABS( NTRD ) ) * NTRDMAX
          END IF
          !CALL OMP_SET_NUM_THREADS( ABS( NTRD ) )
        ELSE
          NTRD = NTRDMAX
        END IF
C         
C         TIME OPEN MP OPERATIONS IF A POSITIVE NTRD VALUE
C         IS SPECIFIED (WHICH MEANS NTRDV IS NOT SPECIFIED)
        TIMEOMP: IF ( NTRD.GT.0 ) THEN
          NTRDV = NTRD
          nop   = 100
C           VECTOR COPY OPERATIONS
          tserialc  = 0.0D0
          tompc     = 0.0D0
C           SERIAL
          DTMISC = 0.0D0
          DO i = 1, nop
            CALL SUPCGDCOPY(1,NTRDV,NIAC,DC,ZC)
          END DO
          tserialc = DTMISC
C           OPENMP
          ntc    = 1
          tompc  = tserialc
          DO np = 2, NTRDMAX
            !CALL OMP_SET_NUM_THREADS( np )
            DTMISC = 0.0D0
            ttc    = 0.0D0
            DO i = 1, nop
              CALL SUPCGDCOPY(NOPT,np,NIAC,DC,ZC)
            END DO
            ttc    = DTMISC
            IF ( ttc.LT.tserialc .AND.
     2           ttc.LT.tompc ) THEN
              ntc    = np
              tompc  = ttc
            END IF       
          END DO
          tserialc  = tserialc  * 1.0D+3 / REAL( nop, 8 )
          tompc     = tompc     * 1.0D+3 / REAL( nop, 8 )
          suc       = 0.0D0
          IF ( tompc.GT.0.0D0 ) suc = tserialc / tompc
C           DOT PRODUCT        
          tserialdp = 0.0D0
          tompdp    = 0.0D0
C           SERIAL
          DTDP = 0.0D0
          DO i = 1, nop
            v = SUPCGDP(1,NTRDV,NIAC,DC,ZC)
          END DO
          tserialdp = DTDP
C           OPENMP
          ntdp   = 1
          tompdp = tserialdp
          DO np = 2, NTRDMAX
           !CALL OMP_SET_NUM_THREADS( np )
           DTDP = 0.0D0
            ttdp = 0.0D0
            DO i = 1, nop
              v = SUPCGDP(NOPT,np,NIAC,DC,ZC)
            END DO
            ttdp    = DTDP
            IF ( ttdp.LT.tserialdp .AND.
     2           ttdp.LT.tompdp ) THEN
              ntdp   = np
              tompdp = ttdp
            END IF       
          END DO
          tserialdp = tserialdp * 1.0D+3 / REAL( nop, 8 )
          tompdp    = tompdp    * 1.0D+3 / REAL( nop, 8 )
          sudp      = 0.0D0
          IF ( tompdp.GT.0.0D0 ) sudp = tserialdp / tompdp
C           SPARSE MATRIX VECTOR PRODUCT        
          tserialmv = 0.0D0
          tompmv    = 0.0D0
C           SERIAL
          DTMV = 0.0D0
          DO i = 1, nop
            CALL SUPCGMV(1,NTRD,NNZC,NIAC,AC,DC,ZC,IAC,JAC)
          END DO
          tserialmv = DTMV
C           OPENMP
          ntmv   = 1
          tompmv = tserialmv
          DO np = 2, NTRDMAX
            !CALL OMP_SET_NUM_THREADS( np )
            DTMV = 0.0D0
            ttmv = 0.0D0
              DO i = 1, nop
                CALL SUPCGMV(NOPT,np,NNZC,NIAC,AC,DC,ZC,IAC,JAC)
              END DO
            ttmv    = DTMV
            IF ( ttmv.LT.tserialmv .AND.
     2           ttmv.LT.tompmv ) THEN
              ntmv   = np
              tompmv = ttmv
            END IF       
          END DO
          tserialmv = tserialmv * 1.0D+3 / REAL( nop, 8 )
          tompmv    = tompmv    * 1.0D+3 / REAL( nop, 8 )
          sumv      = 0.0D0
          IF ( tompmv.GT.0.0D0 ) sumv = tserialmv / tompmv
C           USE OPENMP THREADS THAT RESULT IN THE SMALLEST
C           EXECUTION TIMES
          NTRD  = ntmv
          NTRDV = ntdp
C
C         WRITE SUMMARY OF THE AVERAGE TIME TO COMPLETE OPENMP OPERATIONS
          WRITE (IOUT,570) nop,  
     2                     tserialc,  tompc,  suc,
     3                     tserialdp, tompdp, sudp,
     4                     tserialmv, tompmv, sumv,
     5                     NIAC, NNZC
        ELSE
          NTRD = ABS( NTRD )
        END IF TIMEOMP
C
C---------WRITE SUMMARY OF THREADS BEING USED FOR OPEN MP 
C         SMV AND VECTOR OPERATIONS
        WRITE (IOUT,575) NTRDV, NTRD
C
C---------SET MAXIMUM NUMBER OF THREADS
        NTRDMAX = MAX( NTRD, NTRDV )
        !CALL OMP_SET_NUM_THREADS( NTRDMAX )
      END IF TSTOMP
C      
  570 FORMAT (/1X,75('-'),
     &        /34X,'OPEN MP',
     &        /23X,'AVERAGE TIME FOR',1X,I3,1X,'OPERATIONS',
     &        /47X,'TIME',
     &        /12X,'OPERATION',22X,'MILLISECONDS',6X,'SPEEDUP',
     &        /1X,75('-')
     &        /12X,' SERIAL VECTOR COPY TIME      :',G15.7,
     &        /12X,' OPENMP VECTOR COPY TIME      :',G15.7,1X,F10.3,
     &        /12X,' SERIAL DOT PRODUCT TIME      :',G15.7,
     &        /12X,' OPENMP DOT PRODUCT TIME      :',G15.7,1X,F10.3,
     &        /12X,' SERIAL SMV PRODUCT TIME      :',G15.7,
     &        /12X,' OPENMP SMV PRODUCT TIME      :',G15.7,1X,F10.3,
     &        /1X,75('.'),
     &        /12X,' NUMBER OF ACTIVE CELLS       :',I10,
     &        /12X,' NUMBER OF NON-ZERO ENTRIES   :',I10,
     &        /1X,75('-'))
  575 FORMAT (/1X,75('-'),
     &        /19X,' NUMBER OF OPENMP   V THREADS :',I5,
     &        /19X,' NUMBER OF OPENMP SMV THREADS :',I5,
     &        /1X,75('-'))

C
C-------ALLOCATE GPU POINTERS
      ALLOCATE(CU_HDL,CU_STAT,CU_DES)
      ALLOCATE(CU_IAC,CU_JAC)
      ALLOCATE(CU_AC,CU_APC,CU_XC)
      ALLOCATE(CU_DC,CU_ZC,CU_PC,CU_QC)
      ALLOCATE(CU_SCL,CU_SCLI,CU_V,CU_V0,CU_V1)
      ALLOCATE(PL_DC,PL_ZC)
C
C-------INITIALIZE GPU MEMORY
      IF ( NOPT.EQ.3 ) THEN
C         INITIALIZE MEMORY ON THE GPU 
        CALL UPCGC7_INIT(CU_HDL,CU_STAT,CU_DES,
     &                   NNZC,NIAC,NIAPC,
     &                   NPC,GLSPOLY%NDEGREE,
     &                   CU_IAC,IAC,CU_JAC,JAC,
     &                   CU_AC,CU_APC,CU_XC,
     &                   CU_DC,CU_ZC,CU_PC,CU_QC,
     &                   CU_SCL,CU_SCLI,CU_V,CU_V0,CU_V1,
     &                   PL_DC,PL_ZC)      
      END IF
C
C-------SET POINTERS FOR GRID
      CALL UPCG7PSV(IGRID)
C
C-------RETURN
      RETURN
      END SUBROUTINE UPCG7AR
      
      SUBROUTINE UPCG7AP(HNEW,IBOUND,CR,CC,CV,HCOF,RHS,
     &                 ICNVG,KSTP,KPER,MXITER,KITER,
     &                 NCOL,NROW,NLAY,NODES,HNOFLO,IOUT,
     &                 NPC,NOPT,NTRD,NTRDV,NITER,ITER1,NNZC,NIAC,
     &                 NIAPC,NIWC,NPOL,NEIG,
     &                 HCLOSE,RCLOSE,
     &                 UPCGTOTT,UPCGFMAT,
     &                 UPCGPCUT,UPCGPCAT,UPCGDPT,UPCGMVT,
     &                 UPCGAXPYT,UPCGVVPT,UPCGMISCT,UPCGGPUTT,
     &                 IUPCGO,IUPCGI,
     &                 NODEC,BC,XC,AC,APC,IAC,JAC,IUC,IXMAP,IWC,
     &                 DC,ZC,PC,QC,ISCL,SCL,SCLI,GLSPOLY,
     &                 CU_HDL,CU_STAT,CU_DES,CU_JAC,CU_IAC,
     &                 CU_AC,CU_APC,CU_XC,
     &                 CU_DC,CU_ZC,CU_PC,CU_QC,
     &                 CU_SCL,CU_SCLI,CU_V,CU_V0,CU_V1,
     &                 PL_DC,PL_ZC)
C
C     ******************************************************************
C     SOLUTION BY THE CONJUGATE GRADIENT METHOD -
C                                          UP TO ITER1 ITERATIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE UPCGMODULE, ONLY: TGLSPOLY, 
     2                      DTDP, DTMV, DTAXPY, DTVVP, DTMISC
      IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
      DOUBLEPRECISION, DIMENSION(NODES), INTENT(INOUT) :: HNEW
      INTEGER, DIMENSION(NODES), INTENT(INOUT)         :: IBOUND
      REAL, DIMENSION(NODES), INTENT(IN)               :: CR
      REAL, DIMENSION(NODES), INTENT(IN)               :: CC
      REAL, DIMENSION(NODES), INTENT(IN)               :: CV
      REAL, DIMENSION(NODES), INTENT(IN)               :: HCOF
      REAL, DIMENSION(NODES), INTENT(IN)               :: RHS
      INTEGER, INTENT(INOUT)                           :: ICNVG
      INTEGER, INTENT(IN)                              :: KSTP
      INTEGER, INTENT(IN)                              :: KPER
      INTEGER, INTENT(IN)                              :: MXITER
      INTEGER, INTENT(IN)                              :: KITER
      INTEGER, INTENT(IN)                              :: NCOL
      INTEGER, INTENT(IN)                              :: NROW
      INTEGER, INTENT(IN)                              :: NLAY
      INTEGER, INTENT(IN)                              :: NODES
      REAL, INTENT(IN)                                 :: HNOFLO
      INTEGER, INTENT(IN)                              :: IOUT
      INTEGER, INTENT(IN)                              :: NPC
      INTEGER, INTENT(IN)                              :: NOPT
      INTEGER, INTENT(IN)                              :: NTRD
      INTEGER, INTENT(IN)                              :: NTRDV
      INTEGER, INTENT(IN)                              :: NITER
      INTEGER, INTENT(INOUT)                           :: ITER1
      INTEGER, INTENT(IN)                              :: NNZC
      INTEGER, INTENT(IN)                              :: NIAC
      INTEGER, INTENT(IN)                              :: NIAPC
      INTEGER, INTENT(IN)                              :: NIWC
      INTEGER, INTENT(IN)                              :: NPOL
      INTEGER, INTENT(IN)                              :: NEIG
      REAL, INTENT(IN)                                 :: HCLOSE
      REAL, INTENT(IN)                                 :: RCLOSE
      DOUBLEPRECISION, INTENT(INOUT)                   :: UPCGTOTT
      DOUBLEPRECISION, INTENT(INOUT)                   :: UPCGFMAT
      DOUBLEPRECISION, INTENT(INOUT)                   :: UPCGPCUT
      DOUBLEPRECISION, INTENT(INOUT)                   :: UPCGPCAT
      DOUBLEPRECISION, INTENT(INOUT)                   :: UPCGDPT
      DOUBLEPRECISION, INTENT(INOUT)                   :: UPCGMVT
      DOUBLEPRECISION, INTENT(INOUT)                   :: UPCGAXPYT
      DOUBLEPRECISION, INTENT(INOUT)                   :: UPCGGPUTT
      DOUBLEPRECISION, INTENT(INOUT)                   :: UPCGVVPT
      DOUBLEPRECISION, INTENT(INOUT)                   :: UPCGMISCT
      INTEGER, INTENT(INOUT)                           :: IUPCGO
      INTEGER, INTENT(INOUT)                           :: IUPCGI
      INTEGER, DIMENSION(NCOL,NROW,NLAY)               :: NODEC
      DOUBLEPRECISION, DIMENSION(NIAC), INTENT(INOUT)  :: BC
      DOUBLEPRECISION, DIMENSION(NIAC), INTENT(INOUT)  :: XC
      DOUBLEPRECISION, DIMENSION(NNZC), INTENT(INOUT)  :: AC
      DOUBLEPRECISION, DIMENSION(NIAPC), INTENT(INOUT) :: APC
      INTEGER, DIMENSION(NIAC+1), INTENT(IN)           :: IAC
      INTEGER, DIMENSION(NNZC), INTENT(IN)             :: JAC
      INTEGER, DIMENSION(NIAC), INTENT(IN)             :: IUC
      INTEGER, DIMENSION(NIAC), INTENT(IN)             :: IXMAP
C       WORKING ARRAYS
      INTEGER, DIMENSION(NIWC), INTENT(INOUT)          :: IWC
      DOUBLEPRECISION, DIMENSION(NIAC), INTENT(INOUT)  :: DC
      DOUBLEPRECISION, DIMENSION(NIAC), INTENT(INOUT)  :: PC
      DOUBLEPRECISION, DIMENSION(NIAC), INTENT(INOUT)  :: QC
      DOUBLEPRECISION, DIMENSION(NIAC), INTENT(INOUT)  :: ZC
C       DIAGONAL SCALING VECTOR
      INTEGER, INTENT(IN) :: ISCL
      DOUBLEPRECISION, DIMENSION(NIAC), INTENT(INOUT)  :: SCL
      DOUBLEPRECISION, DIMENSION(NIAC), INTENT(INOUT)  :: SCLI
C       POLYNOMIAL PRECONDITIONER
      TYPE (TGLSPOLY), INTENT(INOUT)                   :: GLSPOLY
C     GPU VARIABLES      
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_HDL
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_STAT
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_DES
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_JAC
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_IAC
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_AC
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_APC      
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_XC
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_DC
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_ZC
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_PC
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_QC
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_SCL    
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_SCLI    
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_V      
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_V0      
      INTEGER(KIND=8),                    INTENT(IN)   :: CU_V1
      INTEGER(KIND=8),                    INTENT(IN)   :: PL_DC
      INTEGER(KIND=8),                    INTENT(IN)   :: PL_ZC
C     + + + LOCAL DEFINITIONS + + +
      DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
      DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
      INTEGER :: i, j, k, n
      INTEGER :: ii
      INTEGER :: iapos
      INTEGER :: ieq
      INTEGER :: nrc
      INTEGER :: iinactive
      INTEGER :: iadiag
      INTEGER :: nrn, nrl, ncn, ncl, nln, nll
      INTEGER :: ncf, ncd, nrb, nrh, nls, nlz
      INTEGER :: ncount
      INTEGER :: iiter
      INTEGER :: irc
      INTEGER :: iicnvg
      INTEGER :: nx, nr

      DOUBLE PRECISION :: dhclose, drclose
      DOUBLE PRECISION :: rrhs, hhcof, rsq, fbar, fmax
      DOUBLE PRECISION :: z, b, d, e, f, h, s
      DOUBLE PRECISION :: zhnew, bhnew, dhnew, fhnew, hhnew, shnew

      DOUBLEPRECISION :: t, ta
      DOUBLEPRECISION :: deltax
      DOUBLEPRECISION :: rmax
      DOUBLEPRECISION :: alpha, beta
      DOUBLEPRECISION :: rho, rho0
      
      DOUBLEPRECISION :: ttf, tt1, tt2
      DOUBLEPRECISION :: tpcu1, tpcu2
      DOUBLEPRECISION :: tpca1, tpca2
      DOUBLEPRECISION :: tdp1, tdp2
      DOUBLEPRECISION :: tmv1, tmv2
      
C       GPU LOCALS DEFINITIONS
      DOUBLEPRECISION :: cu_rho, cu_rho0        

C       FUNCTIONS      
      DOUBLEPRECISION :: SUPCGDP
C
C-------CODE
C
C-------START UPCG TOTAL TIMER
      CALL SUPCGTIMER(0,tt1,UPCGTOTT)
C
C       SET TEMPORARY TIMER VARIABLES TO CURRENT VALUE OF 
C       DUMMY ARGUMENT TIMER VARIABLES
      DTDP   = UPCGDPT
      DTMV   = UPCGMVT
      DTAXPY = UPCGAXPYT
      DTVVP  = UPCGVVPT
      DTMISC = UPCGMISCT
C
C-------START FORMULATE TIMER
      CALL SUPCGTIMER(0,ttf,UPCGFMAT)
C
C       SET LOCAL VARIABLES
      nrc     = NROW * NCOL
      iapos   = 0
      ieq     = 0
      fmax    = dzero
      ncount  = 0
      dhclose = REAL( hclose, 8 )
      drclose = REAL( rclose, 8 )
      DO n = 1, NNZC
        AC(n) = DZERO
      END DO
C
C-------LOOP THROUGH ALL NODES IN THE GRID AND SET UP MATRIX EQUATIONS.
C-------NOTE THAT THE FORMULATION OF THESE EQUATIONS IS OPPOSITE IN SIGN
C-------FROM WHAT IS GIVEN IN THE MODFLOW MANUAL SO THAT THE DIAGONAL
C-------AND RHS ARE BOTH POSITIVE (LHS AND RHS ARE MULTIPLIED BY -1)
C-------THIS LOOP STRUCTURE AND INDEXING IS IDENTICAL TO THAT OF PCG2 
C-------AND IS BLATANTLY COPIED FROM HILL, 1990.
      LFILL: DO k = 1, NLAY
        RFILL: DO i = 1, NROW
          CFILL: DO j = 1, NCOL
C
C---------------CALCULATE 1 DIMENSIONAL SUBSCRIPT OF CURRENT CELL AND
C---------------INITIALIZE MATRIX COEFFICIENTS TO ZERO. CHECK IF CELL IS ACTIVE 
C---------------SKIP COEFFICIENT CALCULATIONS IF CELL IS INACTIVE
            n = j + (i-1) * NCOL + (k-1) * nrc
            e = dzero
            z = dzero
            b = dzero
            d = dzero
            f = dzero
            h = dzero
            s = dzero
            iinactive = 1
            IF( IBOUND(n).GT.0 ) THEN
              iinactive = 0
              ieq       = ieq + 1
C
C---------------CALCULATE 1 DIMENSIONAL SUBSCRIPTS FOR LOCATING THE 6
C---------------SURROUNDING CELLS
              nrn = n + NCOL
              nrl = n - NCOL
              ncn = n + 1
              ncl = n - 1
              nln = n + nrc
              nll = n - nrc
C
C---------------CALCULATE 1 DIMENSIONAL SUBSCRIPTS FOR CONDUCTANCE TO THE 6
C---------------SURROUNDING CELLS.
              ncf = n
              ncd = n - 1
              nrb = n - NCOL
              nrh = n
              nls = n
              nlz = n - nrc
C
C---------------STORE DOUBLE PRECISION VALUE OF RHS FOR CALCULATION OF RESIDUALS
              rrhs    =  RHS(n)
              BC(ieq) = -rrhs
C
C---------------GET CONDUCTANCES TO NEIGHBORING CELLS.  
C---------------ACCUMULATE CONTRIBUTIONS TO DIAGONAL COEFFICIENT. IF NEIGHBOR IS 
C---------------CONSTANT HEAD, MODIFY RHS AND SET OFF-DIAGONAL COEFFICIENT TO 0
C
C
!              iapos  = iapos + 1
              iadiag = IAC(ieq)
              iapos  = iadiag
C
C               TOP FACE
C---------------NEIGHBOR IS 1 LAYER BEHIND
              zhnew = dzero
              IF ( k.NE.1 ) THEN
                z = CV(nlz)
                e = e + z 
                zhnew = z*(HNEW(nll) - HNEW(n))
                IF ( IBOUND(nll).GT.0 ) THEN
                  iapos = iapos + 1
                  AC(iapos) = -z
                END IF
                IF( IBOUND(nll).LT.0 ) THEN
                  BC(ieq) = BC(ieq) + z*HNEW(nll) 
                  z = DZERO
                END IF
              END IF
C
C               UPPER FACE
C---------------NEIGHBOR IS 1 ROW BACK
              bhnew = dzero
              IF ( i.NE.1 ) THEN
                b = CC(nrb)
                e = e + b
                bhnew = b*(HNEW(nrl) - HNEW(n))
                IF( IBOUND(nrl).GT.0 ) THEN
                  iapos = iapos + 1
                  AC(iapos) = -b
                END IF 
                IF( IBOUND(NRL).LT.0 ) THEN
                  BC(ieq) = BC(ieq) + b*HNEW(nrl) 
                  b = dzero
                END IF 
              END IF
C
C               LEFT FACE
C---------------NEIGHBOR IS 1 COLUMN BACK
              dhnew = dzero
              IF ( j.NE.1 ) THEN
                d = CR(ncd)
                e = e + d
                dhnew = d*(HNEW(ncl) - HNEW(n))
                IF( IBOUND(ncl).GT.0 ) THEN
                  iapos = iapos + 1
                  AC(iapos) = -d
                END IF
                IF( IBOUND(ncl).LT.0 ) THEN
                  BC(ieq) = BC(ieq) + d*HNEW(ncl) 
                  d = dzero
                END IF 
              END IF
C
C               RIGHT FACE
C---------------NEIGHBOR IS 1 COLUMN AHEAD
              fhnew = dzero
              IF ( j.NE.NCOL ) THEN
                f = CR(ncf)
                e = e + f
                fhnew = f*(HNEW(ncn) - HNEW(n))
                IF( IBOUND(ncn).GT.0 ) THEN
                  iapos = iapos + 1
                  AC(iapos) = -f
                END IF
                IF( IBOUND(ncn).LT.0 ) THEN
                  BC(ieq) = BC(ieq) + f*HNEW(ncn) 
                  f = dzero 
                END IF
              END IF
C
C               LOWER FACE
C---------------NEIGHBOR IS 1 ROW AHEAD
              hhnew = dzero
              IF ( i.NE.NROW ) THEN
                h = CC(nrh)
                e = e + h
                hhnew = h*(HNEW(nrn) - HNEW(n))
                IF( IBOUND(nrn).GT.0 ) THEN
                  iapos = iapos + 1
                  AC(iapos) = -h
                END IF
                IF( IBOUND(nrn).LT.0 ) THEN
                  BC(ieq) = BC(ieq) + h*HNEW(nrn) 
                  h = dzero
                END IF 
              END IF
C
C               BOTTOM FACE
C---------------NEIGHBOR IS 1 LAYER AHEAD
              shnew = dzero
              IF ( k.NE.NLAY ) THEN
                s = CV(nls)
                e = e + s
                shnew = s*(HNEW(nln) - HNEW(n))
                IF( IBOUND(nln).GT.0 ) THEN
                  iapos = iapos + 1
                  AC(iapos) = -s
                END IF
                IF( IBOUND(nln).LT.0 ) THEN
                  BC(ieq) = BC(ieq) + s*HNEW(nln) 
                  s = dzero
                END IF
              END IF
C    
C---------------CHECK IF SURROUNDING CELLS ARE ACTIVE (E > 0).  IF SO, CALCULATE 
C---------------L2 NORM.  ACCUMULATE THE AVERAGE ABSOLUTE VALUE  OF THE RHS 
C---------------VECTOR FOR ALL ACTIVE CELLS.  THIS IS USED TO SCALE THE THE 
C---------------CLOSURE CRITERIA. 
C---------------IF SURROUNDING CELLS ARE INACTIVE BUT CURRENT CELL IS ACTIVE,
C---------------SET HNEW TO HNOFLO, IBOUND TO 0, AND CHANGE INACTIVE FLAG TO 1
              IF ( e.GT.dzero ) THEN
                hhcof = HNEW(n)*HCOF(n)
                rsq = rsq + (rrhs - zhnew - bhnew - dhnew - hhcof - 
     &             fhnew - hhnew - shnew)**2
                e = e - HCOF(n)
                fbar = fbar + ABS( BC(ieq) )
                ncount = ncount + 1
              ELSE
                HNEW(n)   = HNOFLO
                IBOUND(n) = 0
                iinactive = 1
C                 IF INACTIVE OR CONSTANT HEAD, SET DIAGONAL TO 1.0, AND ADJUST RHS ACCORDINGLY.  
                e = done
                BC(ieq) = HNEW(n)
              END IF
C
C-------------FIND THE MAXIMUM VALUE OF THE RHS VECTOR FOR ALL CELLS (ACTIVE
C-------------AND INACTIVE) FOR CLOSURE SCALING USED BY THE UNSTRUCTURED PCG SOLVER
              fmax = MAX( fmax, ABS( BC(ieq) ) )
C
C---------------STORE THE COEFFICENTS OF THE DIAGONAL IN A
              AC(iadiag) = e
C---------------STORE INITIAL GUESS OF HEADS
              XC(ieq) = HNEW(n)     
C
C-------------END IBOUND(N) .GT. 0
            END IF
          END DO CFILL
        END DO RFILL
      END DO LFILL
C
C-------END FORMULATE TIMER
      CALL SUPCGTIMER(1,ttf,UPCGFMAT)
C
C-------UPDATE UPCG OUTER ITERATION COUNTER
      IUPCGO = IUPCGO + 1
C
C-------UPDATE PRECONDITIONER
      CALL SUPCGTIMER(0,tpcu1,UPCGPCUT)
C       SCALE MATRIX FOR POLYNOMIAL PRECONDITIONER
!      IF ( NPC.EQ.4 ) THEN
      IF ( ISCL.NE.0 ) THEN
        CALL SUPCGSCL(1,NNZC,NIAC,AC,XC,BC,SCL,SCLI,IAC,JAC)
      END IF
      CALL SUPCGPCU(NOPT,NTRD,NTRDV,NNZC,NIAC,NIAPC,NIWC,NPC,
     2              AC,APC,IAC,JAC,IUC,IWC,
     3              GLSPOLY)
      CALL SUPCGTIMER(1,tpcu1,UPCGPCUT)
C
C-------WRITE SUMMARY OF INITIAL EIGENVALUES FOR
C       THE POLYNOMIAL PRECONDITIONER
      IF ( NPC.EQ.4 ) THEN
        IF ( KPER.EQ.1 ) THEN
          IF ( KSTP.EQ.1 ) THEN
            IF ( KITER.EQ.1 ) THEN
              WRITE (IOUT,2000) GLSPOLY%INTV(1), GLSPOLY%INTV(2)
            END IF
          END IF
        END IF
      END IF
2000  FORMAT(//1X,'INITIAL EIGENVALUE RANGE FOR STRESS PERIOD 1, ',
     2            'TIME STEP 1, AND ITERATION 1',/1X,80('-'),
     3        /1X,'MINIMUM EIGENVALUE:',G15.7,            
     4        /1X,'MAXIMUM EIGENVALUE:',G15.7,/1X,80('-'))            
C-------INITIALIZE SOLUTION VARIABLE AND ARRAYS
      iiter = 0
      IF ( KITER.EQ.1 ) ITER1 = 0
      irc    = 1
      ICNVG  = 0
      iicnvg = 0
      alpha  = dzero
      beta   = dzero
      rho    = dzero
      rho0   = dzero
      DO n = 1, NIAC
        DC(n) = DZERO
        PC(n) = DZERO
        QC(n) = DZERO
        ZC(n) = DZERO
      END DO
C-------CALCULATE INITIAL RESIDUAL
      CALL SUPCGMV(NOPT,NTRD,NNZC,NIAC,AC,XC,DC,IAC,JAC)
      DO n = 1, NIAC
        t     = DC(n)
        DC(n) = BC(n) - t
      END DO
C       CPU
      CPUGPUT: IF ( NOPT.EQ.1 .OR. NOPT.EQ.2 ) THEN
C---------INNER ITERATION          
        INNER: DO iiter = 1, NITER
           IUPCGI = IUPCGI + 1
           ITER1  = ITER1  + 1
C-----------APPLY PRECONDITIONER
          CALL SUPCGTIMER(0,tpca1,UPCGPCAT)
          SELECT CASE (NPC)
C             NO PRECONDITIONER
            CASE (0)
              DO n = 1, NIAC
                ZC(n) = DC(n)
              END DO
C             JACOBI PRECONDITIONER
            CASE (1)
              CALL SUPCGJACA(NOPT,NTRDV,NIAC,APC,DC,ZC)
C             ILU0 AND MILU0 PRECONDITIONERS
            CASE (2,3)
              CALL SUPCGILU0A(NNZC,NIAC,NIAPC,
     2                        APC,IAC,JAC,IUC,DC,ZC)
C             POLYNOMIAL PRECONDITIONER
            CASE (4)
              CALL SUPCGPOLYA(NOPT,NTRD,NTRDV,NNZC,NIAC,
     2                        AC,IAC,JAC,GLSPOLY,DC,ZC)
          END SELECT
          CALL SUPCGTIMER(1,tpca1,UPCGPCAT)

          rho = SUPCGDP( NOPT,NTRDV,NIAC,DC,ZC )
C-----------COMPUTE DIRECTIONAL VECTORS
          IF (iiter.EQ.1) THEN
            DO n = 1, NIAC
              PC(n) = ZC(n)
            END DO
          ELSE
            beta = rho / rho0
            DO n = 1, NIAC
              PC(n) = ZC(n) + beta * PC(n)
            END DO
          END IF
C-----------COMPUTE ITERATES
C           UPDATE qc
          CALL SUPCGMV(NOPT,NTRD,NNZC,NIAC,AC,PC,QC,IAC,JAC)

          alpha = rho / SUPCGDP( NOPT,NTRDV,NIAC,PC,QC)
C-----------UPDATE X AND RESIDUAL
          deltax = DZERO
          rmax   = DZERO
C-----------UNSCALE HEAD CHANGE AND RESIDUAL FOR POLYNOMIAL PRECONDITION
          DO n = 1, NIAC
            t      = alpha * PC(n)
            XC(n)  = XC(n) + t
            ta     = t * SCL(n)
            IF ( ABS( ta ).GT.deltax ) THEN
              nx     = n
              deltax = ABS( ta )
            END IF
            t      = DC(n)
            t      = t - alpha * QC(n)
            DC(n)  = t
            ta     = t * SCLI(n)
            IF ( ABS( ta ).GT.rmax ) THEN
              nr   = n
              rmax = ABS( ta )
            END IF
          END DO
          IF ( deltax.LE.dhclose .AND. rmax.LE.drclose ) THEN
            iicnvg = 1
          END IF
          IF ( MXITER.EQ.1 ) THEN
            IF ( iicnvg.EQ.1 ) ICNVG = 1
          ELSE
            IF ( iiter.EQ.1 .AND. iicnvg.EQ.1 ) ICNVG = 1
          ENDIF
          IF ( iicnvg.EQ.1 ) EXIT INNER
C-----------SAVE CURRENT INNER ITERATES
          rho0 = rho
        END DO INNER
C         GPU
      ELSE
        CALL UPCGC7(iiter,CU_HDL,CU_STAT,CU_DES,
     &              NNZC,NIAC,NIAPC,
     &              NPC,GLSPOLY%NDEGREE,
     &              CU_IAC,IAC,CU_JAC,JAC,IUC,
     &              CU_AC,AC,CU_APC,APC,CU_XC,XC,
     &              CU_DC,DC,CU_ZC,ZC,
     &              CU_PC,CU_QC,
     &              GLSPOLY%ALPHA,GLSPOLY%BETA,GLSPOLY%GAMMA,
     &              CU_SCL,SCL,CU_SCLI,SCLI,
     &              CU_V,GLSPOLY%D_V,
     &              CU_V0,GLSPOLY%D_V0,
     &              CU_V1,GLSPOLY%D_V1,
     &              PL_DC,PL_ZC,
     &              cu_rho0,cu_rho,
     &              MXITER,ICNVG,NITER,dhclose,drclose,
     &              deltax,rmax,
     &              UPCGPCAT,UPCGDPT,UPCGMVT,
     &              UPCGAXPYT,UPCGVVPT,UPCGMISCT,UPCGGPUTT)
         IUPCGI = IUPCGI + iiter
         ITER1  = ITER1  + iiter

      END IF CPUGPUT
C
C       UNSCALE XC
!      IF ( NPC.EQ.4 ) THEN
      IF ( ISCL.NE.0 ) THEN
        CALL SUPCGSCL(0,NNZC,NIAC,AC,XC,BC,SCL,SCLI,IAC,JAC)
      END IF
C
C-------FILL HNEW WITH NEW ESTIMATE
      DO n = 1, NIAC
        HNEW(IXMAP(n)) = XC(n)
      END DO
C
C-------IF END OF TIME STEP, PRINT # OF ITERATIONS THIS STEP
      IF ( ICNVG.NE.0 .OR. KITER.EQ.MXITER ) THEN
        WRITE (IOUT,510)
  510   FORMAT (1X,/1X)
        WRITE (IOUT,515) KITER, KSTP, KPER, ITER1
  515   FORMAT (I6,' CALLS TO PCG ROUTINE FOR TIME STEP',I4,
     &          ' IN STRESS PERIOD ',I4,/,I6,' TOTAL ITERATIONS')
        ITER1 = 0
      ENDIF
C
C       SET DUMMY ARGUMENT TIMER VARIABLES TO CURRENT VALUE OF
C       TEMPORARY TIMER VARIABLES
      UPCGDPT     = DTDP
      UPCGMVT     = DTMV
      UPCGAXPYT   = DTAXPY
      UPCGVVPT    = DTVVP
      UPCGMISCT   = DTMISC
C
C-------END UPCG TIMER
      CALL SUPCGTIMER(1,tt1,UPCGTOTT)
C
C-------RETURN
      RETURN
C
      END SUBROUTINE UPCG7AP

      SUBROUTINE UPCG7OT(IGRID)
C     ******************************************************************
C     OUTPUT UPCG TIMER RESULTS - FLOATING POINT OPERATIONS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT
      USE UPCGMODULE
      IMPLICIT NONE
C
C     + + + DUMMY ARGUMENTS + + +
      INTEGER, INTENT(IN) :: IGRID
C     + + + LOCAL DEFINITIONS + + +
      DOUBLEPRECISION :: tupcg
!      INTEGER :: i
C     + + + FUNCTIONS + + +
C     + + + FORMATS + + +
2000  FORMAT(//,1X,'SUMMARY OF UPCG EXECUTION TIME',
     &        /,1X,'TOTAL NUMBER OF OUTER UPCG ITERATIONS: ',I10,
     &        /,1X,'TOTAL NUMBER OF INNER UPCG ITERATIONS: ',I10,
     &       //,1X,'TIMER ITEM',32X,'TIME (SEC.)',4X,' PERCENTAGE',
     &        /,1X,69('-'),
     &        /,1X,'TOTAL UPCG EXECUTION TIME:             ',G16.7,
     &          1X,F10.3,
     &        /,1X,'TOTAL UPCG A MATRIX FORMULATE TIME:    ',G16.7,
     &          1X,F10.3,
     &        /,1X,'TOTAL UPCG PRECOND. FORMULATE TIME:    ',G16.7,
     &          1X,F10.3,
     &        /,1X,'TOTAL UPCG PRECOND. SOLUTION TIME:     ',G16.7,
     &          1X,F10.3,
     &        /,1X,'TOTAL UPCG DOT PROD. EXECUTION TIME:   ',G16.7,
     &          1X,F10.3,
     &        /,1X,'TOTAL UPCG MV PROD. EXECUTION TIME:    ',G16.7,
     &          1X,F10.3,
     &        /,1X,'TOTAL UPCG AXPY EXECUTION TIME:        ',G16.7,
     &          1X,F10.3,
     &        /,1X,'TOTAL UPCG VV PROD. EXECUTION TIME:    ',G16.7,
     &          1X,F10.3,
     &        /,1X,'TOTAL UPCG MISC. EXECUTION TIME:       ',G16.7,
     &          1X,F10.3,
     &        /,1X,'TOTAL UPCG GPU MEMORY TRANSFER TIME:   ',G16.7,
     &          1X,F10.3,
     &        /,1X,69('-'))
C     + + + CODE + + +
C
C         SET POINTERS
        CALL UPCG7PNT(IGRID)
C
        tupcg = 1.0D0
        IF ( UPCGTOTT.GT.0.0D0 ) tupcg = tupcg / UPCGTOTT
        WRITE (IOUT,2000) IUPCGO, IUPCGI, 
     &                    UPCGTOTT,  100.0D0 * UPCGTOTT * tupcg,
     &                    UPCGFMAT,  100.0D0 * UPCGFMAT * tupcg,
     &                    UPCGPCUT,  100.0D0 * UPCGPCUT * tupcg, 
     &                    UPCGPCAT,  100.0D0 * UPCGPCAT * tupcg, 
     &                    UPCGDPT,   100.0D0 * UPCGDPT  * tupcg,
     &                    UPCGMVT,   100.0D0 * UPCGMVT  * tupcg,
     &                    UPCGAXPYT, 100.0D0 * UPCGAXPYT* tupcg, 
     &                    UPCGVVPT,  100.0D0 * UPCGVVPT * tupcg,
     &                    UPCGMISCT, 100.0D0 * UPCGMISCT* tupcg,
     &                    UPCGGPUTT, 100.0D0 * UPCGGPUTT* tupcg
     
!C
!C-------SAVE FINAL A MATRIX IN FULL N X N FORM
!      OPEN(UNIT=99,FILE='AMATRIX.BIN',FORM='BINARY',STATUS='REPLACE')
!      WRITE(99) NIAC, NNZC, IAC, JAC, AC
!      CLOSE( 99 )
!C
!C-------SAVE FINAL A MATRIX IN FULL N X N FORM
!      OPEN(UNIT=99,FILE='APCMATRIX.DAT',STATUS='REPLACE')
!      DO i = 1, NIAPC
!        WRITE(99,*) APC(i)
!      END DO
!      CLOSE( 99 )
C
C-------RETURN
      RETURN
C
      END SUBROUTINE UPCG7OT
C
C
      SUBROUTINE UPCG7DA(IGRID)
C  Deallocate UPCG DATA
        USE UPCGMODULE
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: IGRID
C     + + + LOCAL DEFINITIONS + + +
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C
C         SET POINTERS
        CALL UPCG7PNT(IGRID)
C
C         FREE GPU MEMORY
        IF ( NOPT.EQ.3 ) THEN
          CALL UPCGC7_FINAL(NPC,CU_HDL,
     &                      CU_JAC,CU_IAC,
     &                      CU_AC,CU_APC,CU_XC,
     &                      CU_DC,CU_ZC,CU_PC,CU_QC,
     &                      PL_DC,PL_ZC)
        END IF
C
C         DEALLOCATE UPCG MEMORY
        DEALLOCATE(ITER1C,NPC,NOPT,NTRD,NTRDV)
        DEALLOCATE(NITERC,NNZC,NIAC)
        DEALLOCATE(NIAPC,NIWC,NPOL,NEIG)
        DEALLOCATE(HCLOSEUPCG,RCLOSEUPCG)
        DEALLOCATE(UPCGTOTT,UPCGFMAT)
        DEALLOCATE(UPCGPCUT,UPCGPCAT,UPCGDPT,UPCGMVT)
        DEALLOCATE(UPCGAXPYT,UPCGVVPT,UPCGMISCT)
        DEALLOCATE(UPCGGPUTT)
        DEALLOCATE(DTDP,DTMV,DTAXPY,DTVVP,DTMISC)
        DEALLOCATE(IUPCGO,IUPCGI)
        DEALLOCATE(NODEC)
        DEALLOCATE(BC)
        DEALLOCATE(XC)
        DEALLOCATE(AC)
        DEALLOCATE(APC)
        DEALLOCATE(IAC)
        DEALLOCATE(JAC)
        DEALLOCATE(IUC)
        DEALLOCATE(IXMAP)
C       WORKING ARRAYS
        DEALLOCATE(IWC)
        DEALLOCATE(DC)
        DEALLOCATE(PC)
        DEALLOCATE(QC)
        DEALLOCATE(ZC)
        DEALLOCATE(ISCL)
        DEALLOCATE(SCL)
        DEALLOCATE(SCLI)
        DEALLOCATE(GLSPOLY)
C       GPU POINTERS
        DEALLOCATE(CU_HDL,CU_STAT,CU_DES)
        DEALLOCATE(CU_IAC,CU_JAC)
        DEALLOCATE(CU_AC,CU_APC,CU_XC)
        DEALLOCATE(CU_DC,CU_ZC,CU_PC,CU_QC)
        DEALLOCATE(CU_SCL,CU_SCLI,CU_V,CU_V0,CU_V1)
        DEALLOCATE(PL_DC,PL_ZC)        
C
C---------RETURN
      RETURN
      END SUBROUTINE UPCG7DA
      
      SUBROUTINE UPCG7PNT(IGRID)
C  Set pointers to UPCG data for a grid
      USE UPCGMODULE
C
      ITER1C=>UPCGDAT(IGRID)%ITER1C
      NPC=>UPCGDAT(IGRID)%NPC
      NOPT=>UPCGDAT(IGRID)%NOPT
      NTRD=>UPCGDAT(IGRID)%NTRD
      NTRDV=>UPCGDAT(IGRID)%NTRDV
      NITERC=>UPCGDAT(IGRID)%NITERC
      NNZC=>UPCGDAT(IGRID)%NNZC
      NIAC=>UPCGDAT(IGRID)%NIAC
      NIAPC=>UPCGDAT(IGRID)%NIAPC
      NIWC=>UPCGDAT(IGRID)%NIWC
      NPOL=>UPCGDAT(IGRID)%NPOL
      NEIG=>UPCGDAT(IGRID)%NEIG
      HCLOSEUPCG=>UPCGDAT(IGRID)%HCLOSEUPCG
      RCLOSEUPCG=>UPCGDAT(IGRID)%RCLOSEUPCG
      UPCGTOTT=>UPCGDAT(IGRID)%UPCGTOTT
      UPCGFMAT=>UPCGDAT(IGRID)%UPCGFMAT
      UPCGPCUT=>UPCGDAT(IGRID)%UPCGPCUT
      UPCGPCAT=>UPCGDAT(IGRID)%UPCGPCAT
      UPCGDPT=>UPCGDAT(IGRID)%UPCGDPT
      UPCGMVT=>UPCGDAT(IGRID)%UPCGMVT
      UPCGAXPYT=>UPCGDAT(IGRID)%UPCGAXPYT
      UPCGVVPT=>UPCGDAT(IGRID)%UPCGVVPT
      UPCGMISCT=>UPCGDAT(IGRID)%UPCGMISCT
      UPCGGPUTT=>UPCGDAT(IGRID)%UPCGGPUTT
      IUPCGO=>UPCGDAT(IGRID)%IUPCGO
      IUPCGI=>UPCGDAT(IGRID)%IUPCGI
      NODEC=>UPCGDAT(IGRID)%NODEC
      BC=>UPCGDAT(IGRID)%BC
      XC=>UPCGDAT(IGRID)%XC
      AC=>UPCGDAT(IGRID)%AC
!      TAPC=>UPCGDAT(IGRID)%TAPC
      APC=>UPCGDAT(IGRID)%APC
      IAC=>UPCGDAT(IGRID)%IAC
      JAC=>UPCGDAT(IGRID)%JAC
      IUC=>UPCGDAT(IGRID)%IUC
      IXMAP=>UPCGDAT(IGRID)%IXMAP
C       WORKING ARRAYS
      IWC=>UPCGDAT(IGRID)%IWC
      DC=>UPCGDAT(IGRID)%DC
      PC=>UPCGDAT(IGRID)%PC
      QC=>UPCGDAT(IGRID)%QC
      ZC=>UPCGDAT(IGRID)%ZC
C       POLYNOMIAL PRECONDITIONER
      ISCL=>UPCGDAT(IGRID)%ISCL
      SCL=>UPCGDAT(IGRID)%SCL
      SCLI=>UPCGDAT(IGRID)%SCLI
      GLSPOLY=>UPCGDAT(IGRID)%GLSPOLY
C       GPU POINTERS
      CU_HDL=>UPCGDAT(IGRID)%CU_HDL
      CU_STAT=>UPCGDAT(IGRID)%CU_STAT
      CU_DES=>UPCGDAT(IGRID)%CU_DES
      CU_IAC=>UPCGDAT(IGRID)%CU_IAC
      CU_JAC=>UPCGDAT(IGRID)%CU_JAC
      CU_AC=>UPCGDAT(IGRID)%CU_AC
      CU_APC=>UPCGDAT(IGRID)%CU_APC
      CU_XC=>UPCGDAT(IGRID)%CU_XC
      CU_DC=>UPCGDAT(IGRID)%CU_DC
      CU_ZC=>UPCGDAT(IGRID)%CU_ZC
      CU_PC=>UPCGDAT(IGRID)%CU_PC
      CU_QC=>UPCGDAT(IGRID)%CU_QC
      CU_SCL=>UPCGDAT(IGRID)%CU_SCL
      CU_SCLI=>UPCGDAT(IGRID)%CU_SCLI
      CU_V=>UPCGDAT(IGRID)%CU_V
      CU_V0=>UPCGDAT(IGRID)%CU_V0
      CU_V1=>UPCGDAT(IGRID)%CU_V1
      PL_DC=>UPCGDAT(IGRID)%PL_DC
      PL_ZC=>UPCGDAT(IGRID)%PL_ZC
C
      RETURN
      END SUBROUTINE UPCG7PNT

      SUBROUTINE UPCG7PSV(IGRID)
C  Save pointers to UPCG data
      USE UPCGMODULE
C
      UPCGDAT(IGRID)%ITER1C=>ITER1C
      UPCGDAT(IGRID)%NPC=>NPC
      UPCGDAT(IGRID)%NOPT=>NOPT
      UPCGDAT(IGRID)%NTRD=>NTRD
      UPCGDAT(IGRID)%NTRDV=>NTRDV
      UPCGDAT(IGRID)%NITERC=>NITERC
      UPCGDAT(IGRID)%NNZC=>NNZC
      UPCGDAT(IGRID)%NIAC=>NIAC
      UPCGDAT(IGRID)%NIAPC=>NIAPC
      UPCGDAT(IGRID)%NIWC=>NIWC
      UPCGDAT(IGRID)%NPOL=>NPOL
      UPCGDAT(IGRID)%NEIG=>NEIG
      UPCGDAT(IGRID)%HCLOSEUPCG=>HCLOSEUPCG
      UPCGDAT(IGRID)%RCLOSEUPCG=>RCLOSEUPCG
      UPCGDAT(IGRID)%UPCGTOTT=>UPCGTOTT
      UPCGDAT(IGRID)%UPCGFMAT=>UPCGFMAT
      UPCGDAT(IGRID)%UPCGPCUT=>UPCGPCUT
      UPCGDAT(IGRID)%UPCGPCAT=>UPCGPCAT
      UPCGDAT(IGRID)%UPCGDPT=>UPCGDPT
      UPCGDAT(IGRID)%UPCGMVT=>UPCGMVT
      UPCGDAT(IGRID)%UPCGAXPYT=>UPCGAXPYT
      UPCGDAT(IGRID)%UPCGVVPT=>UPCGVVPT
      UPCGDAT(IGRID)%UPCGMISCT=>UPCGMISCT
      UPCGDAT(IGRID)%UPCGGPUTT=>UPCGGPUTT
      UPCGDAT(IGRID)%IUPCGO=>IUPCGO
      UPCGDAT(IGRID)%IUPCGI=>IUPCGI
      UPCGDAT(IGRID)%NODEC=>NODEC
      UPCGDAT(IGRID)%BC=>BC
      UPCGDAT(IGRID)%XC=>XC
      UPCGDAT(IGRID)%AC=>AC
!      UPCGDAT(IGRID)%TAPC=>TAPC
      UPCGDAT(IGRID)%APC=>APC
      UPCGDAT(IGRID)%IAC=>IAC
      UPCGDAT(IGRID)%JAC=>JAC
      UPCGDAT(IGRID)%IUC=>IUC
      UPCGDAT(IGRID)%IXMAP=>IXMAP
C       WORKING ARRAYS
      UPCGDAT(IGRID)%IWC=>IWC
      UPCGDAT(IGRID)%DC=>DC
      UPCGDAT(IGRID)%PC=>PC
      UPCGDAT(IGRID)%QC=>QC
      UPCGDAT(IGRID)%ZC=>ZC
C       POLYNOMIAL PRECONDITIONER
      UPCGDAT(IGRID)%ISCL=>ISCL
      UPCGDAT(IGRID)%SCL=>SCL
      UPCGDAT(IGRID)%SCLI=>SCLI
      UPCGDAT(IGRID)%GLSPOLY=>GLSPOLY
C       GPU POINTERS
      UPCGDAT(IGRID)%CU_HDL=>CU_HDL
      UPCGDAT(IGRID)%CU_STAT=>CU_STAT
      UPCGDAT(IGRID)%CU_DES=>CU_DES
      UPCGDAT(IGRID)%CU_IAC=>CU_IAC
      UPCGDAT(IGRID)%CU_JAC=>CU_JAC
      UPCGDAT(IGRID)%CU_AC=>CU_AC
      UPCGDAT(IGRID)%CU_APC=>CU_APC
      UPCGDAT(IGRID)%CU_XC=>CU_XC
      UPCGDAT(IGRID)%CU_DC=>CU_DC
      UPCGDAT(IGRID)%CU_ZC=>CU_ZC
      UPCGDAT(IGRID)%CU_PC=>CU_PC
      UPCGDAT(IGRID)%CU_QC=>CU_QC
      UPCGDAT(IGRID)%CU_SCL=>CU_SCL
      UPCGDAT(IGRID)%CU_SCLI=>CU_SCLI
      UPCGDAT(IGRID)%CU_V=>CU_V
      UPCGDAT(IGRID)%CU_V0=>CU_V0
      UPCGDAT(IGRID)%CU_V1=>CU_V1
      UPCGDAT(IGRID)%PL_DC=>PL_DC
      UPCGDAT(IGRID)%PL_ZC=>PL_ZC
C
      RETURN
      END SUBROUTINE UPCG7PSV
C
C-------ROUTINE TO SCALE THE COEFFICIENT MATRIX
      SUBROUTINE SUPCGSCL(ISCALE,NNZC,NIAC,AC,XC,BC,SCL,SCLI,IAC,JAC)
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: ISCALE
        INTEGER, INTENT(IN) :: NNZC
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NNZC),  INTENT(INOUT)  :: AC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT)  :: XC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT)  :: BC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT)  :: SCL
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT)  :: SCLI
        INTEGER, DIMENSION(NIAC+1), INTENT(IN)   :: IAC
        INTEGER, DIMENSION(NNZC), INTENT(IN)     :: JAC
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n, i
        INTEGER :: ic
        INTEGER :: i0, i1
        DOUBLEPRECISION :: c1, c2, v
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C
C---------SCALE AC, XC, AND BC
        IF ( ISCALE.NE.0 ) THEN
          DO n = 1, NIAC
            ic = IAC(n)
            SCL(n)  = 1.0D0 / SQRT( ABS( AC(ic) ) )
            SCLI(n) = SQRT( ABS( AC(ic) ) )
          END DO
          DO n = 1, NIAC
            c1 = SCL(n)
            i0 = IAC(n)
            i1 = IAC(n+1) - 1
            DO i = i0, i1
              ic    = JAC(i)
              c2    = SCL(ic)
              v     = c1 * AC(i) * c2
              AC(i) = v
            END DO
          END DO
C-----------SCALE XC AND BC
          DO n = 1, NIAC
            c1     = SCL(n)
            XC(n)  = XC(n) / c1
            BC(n)  = BC(n) * c1
          END DO
C---------UNSCALE XC -- NO NEED TO UNSCALE AC AND BC BECAUSE THEY ARE NOT REUSED
        ELSE
          DO n = 1, NIAC
            c1 = SCL(n)
!            c2 = SCLI(n)
            i0 = IAC(n)
            i1 = IAC(n+1) - 1
!C             UNSCALE AC
!            DO i = i0, i1
!              jc = JAC(i)
!              c2 = SCL(jc)
!              AC(i) = ( 1.0D0 / c1 ) * AC(i) * ( 1.0D0 / c2 ) 
!            END DO
C             UNSCALE XC
            XC(n) = XC(n) * c1
!            BC(n) = BC(n) / c1
!            BC(n) = BC(n) * c2
          END DO     
        END IF
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGSCL
C
C-------ROUTINE TO UPDATE THE PRECONDITIONER
      SUBROUTINE SUPCGPCU(NOPT,NTRD,NTRDV,NNZC,NIAC,NIAPC,NIWC,NPC,
     2                    AC,APC,IAC,JAC,IUC,IWC,
     3                    GLSPOLY)
        USE UPCGMODULE, ONLY: TGLSPOLY
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NOPT
        INTEGER, INTENT(IN) :: NTRD
        INTEGER, INTENT(IN) :: NTRDV
        INTEGER, INTENT(IN) :: NNZC
        INTEGER, INTENT(IN) :: NIAC
        INTEGER, INTENT(IN) :: NIAPC
        INTEGER, INTENT(IN) :: NIWC
        INTEGER, INTENT(IN) :: NPC
        DOUBLEPRECISION, DIMENSION(NNZC),  INTENT(IN)     :: AC
!        DOUBLEPRECISION, DIMENSION(NIAPC), INTENT(INOUT)  :: TAPC
        DOUBLEPRECISION, DIMENSION(NIAPC), INTENT(INOUT)  :: APC
        INTEGER, DIMENSION(NIAC+1), INTENT(IN)   :: IAC
        INTEGER, DIMENSION(NNZC), INTENT(IN)     :: JAC
        INTEGER, DIMENSION(NIAC), INTENT(IN)     :: IUC
        INTEGER, DIMENSION(NIWC), INTENT(INOUT)  :: IWC
        TYPE (TGLSPOLY), INTENT(INOUT)           :: GLSPOLY
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        SELECT CASE(NPC)
C           NO PRE-CONDITIONER
          CASE (0)
C           JACOBI PRE-CONDITIONER
          CASE (1)
            CALL SUPCGPCJ(NNZC,NIAC,AC,APC,IAC)
C           ILU0
          CASE (2,3)
            CALL SUPCGPCILU0(NPC,NNZC,NIAC,NIAPC,NIWC,
     2                       AC,APC,IAC,JAC,IUC,IWC)
C           NEUMAN POLYNOMIAL
          CASE (4)
            CALL SUPCGGLSPOL(NOPT,NTRD,NTRDV,NNZC,NIAC,AC,IAC,JAC,
     2                       GLSPOLY)
C           ADDITIONAL PRECONDITIONERS - ILU, etc.
        END SELECT
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGPCU
C
C-------JACOBI PRECONDITIONER - INVERSE OF DIAGONAL 
      SUBROUTINE SUPCGPCJ(NNZC,NIAC,AC,APC,IAC)
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NNZC
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NNZC),  INTENT(IN)      :: AC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT)   :: APC
        INTEGER, DIMENSION(NIAC+1), INTENT(IN) :: IAC
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: i, n
        INTEGER :: ic0, ic1
        INTEGER :: id
        DOUBLEPRECISION :: t
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
        DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        DO n = 1, NIAC
            id = IAC(n)
            t  = AC(id)
            IF ( ABS( t ).EQ.DZERO ) THEN
              CALL USTOP('SUPCGPCJ ERROR: ABS(AC)=0.0')
            END IF
            APC(n) = DONE / t
        END DO
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGPCJ

      SUBROUTINE SUPCGJACA(NOPT,NTRDV,NIAC,A,D1,D2)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NOPT
        INTEGER, INTENT(IN) :: NTRDV
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: A
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: D1
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT) :: D2
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        DOUBLEPRECISION :: t, djac
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
        DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        CALL SUPCGVVP(NOPT,NTRDV,NIAC,A,D1,D2)
!        SELECT CASE ( NOPT )
!C           CPU
!          CASE (1)
!            DO n = 1, NIAC
!              t     = A(n) * D1(n)
!              D2(n) = t
!            END DO
!C           CPU - OPEN MP          
!          CASE (2)
!!$OMP  PARALLEL
!!$OMP& NUM_THREADS(NTRDV)
!!$OMP& DEFAULT(SHARED)
!!$OMP& PRIVATE(n, t)
!!$OMP  DO
!            DO n = 1, NIAC
!C               MULTIPLY INVERSE OF DIAGONAL AND D1
!              t     = A(n) * D1(n)
!              D2(n) = t
!            END DO
!!$OMP  END DO
!!$OMP  END PARALLEL
!        END SELECT
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGJACA

      SUBROUTINE SUPCGPCILU0(NPC,NNZC,NIAC,NIAPC,NIWC,
     2                       AC,APC,IAC,JAC,IUC,IWC)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NPC
        INTEGER, INTENT(IN) :: NNZC
        INTEGER, INTENT(IN) :: NIAC
        INTEGER, INTENT(IN) :: NIAPC
        INTEGER, INTENT(IN) :: NIWC
        DOUBLEPRECISION, DIMENSION(NNZC),  INTENT(IN)     :: AC
        DOUBLEPRECISION, DIMENSION(NIAPC), INTENT(INOUT)  :: APC
        INTEGER, DIMENSION(NIAC+1), INTENT(IN)   :: IAC
        INTEGER, DIMENSION(NNZC), INTENT(IN)     :: JAC
        INTEGER, DIMENSION(NIAC), INTENT(IN)     :: IUC
        INTEGER, DIMENSION(NIWC), INTENT(INOUT)  :: IWC
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: ic0, ic1, id0, iu1
        INTEGER :: iic0, iic1
        INTEGER :: j, n
        INTEGER :: jj, nn
        INTEGER :: jpos, jcol, jw
        INTEGER :: id
        INTEGER :: izero
        DOUBLEPRECISION :: tl
        DOUBLEPRECISION :: t
        DOUBLEPRECISION :: rs
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
        DOUBLEPRECISION, PARAMETER :: DONE  = 1.0D0
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        izero = 0
        DO n = 1, NIAPC
          APC(n) = AC(n)
        END DO
        DO n = 1, NIAC
          IWC(n)  = 0
        END DO
        MAIN: DO n = 1, NIAC
          ic0 = IAC(n)
          ic1 = IAC(n+1)-1
          DO j = ic0, ic1
            jcol = JAC(j)
            IWC(jcol) = j
          END DO
          rs    = DZERO
          jpos = IAC(n)
          iu1 = IUC(n) - 1
          LOWER: DO j = ic0+1, iu1
            jpos = j
            jcol = JAC(j)
            id0  = IAC(jcol)
            tl   = APC(j) * APC(id0)
            APC(j) = tl
            iic0 = IUC(jcol)
            iic1 = IAC(jcol+1) - 1
            DO jj = iic0, iic1
              jw = IWC(JAC(jj))
              IF ( jw.NE.0 ) THEN
                APC(jw) = APC(jw) - tl * APC(jj)
              ELSE
                IF ( NPC.EQ.3 ) rs = rs + tl * APC(jj)
                !rs = rs + tl * APC(jj)
              END IF
            END DO
          END DO LOWER
C           DIAGONAL - CALCULATE INVERSE OF DIAGONAL FOR SOLUTION
          id0 = IAC(n)
          tl  = APC(id0) - rs
          IF ( tl.GT.DZERO ) THEN
            APC(id0) = DONE / tl
          ELSE
            CALL USTOP('SUPCGPCILU0: tl <= 0.0')
            !izero = 1
            !EXIT MAIN
            !APC(id0) = 1.0D+20
          END IF
C           RESET POINTER FOR IW TO ZERO
          DO j = ic0, ic1
            jcol = JAC(j)
            IWC(jcol) = 0
          END DO
        END DO MAIN
C---------REVERT TO A IF ZERO ON DIAGONAL ENCOUNTERED
        IF ( izero.NE.0 ) THEN
          DO n = 1, NIAPC
            APC(n) = AC(n)
          END DO
        END IF
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGPCILU0

      SUBROUTINE SUPCGILU0A(NNZC,NIAC,NIAPC,
     2                      APC,IAC,JAC,IUC,R,D)
!        USE OMP_LIB
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NNZC
        INTEGER, INTENT(IN) :: NIAC
        INTEGER, INTENT(IN) :: NIAPC
        DOUBLEPRECISION, DIMENSION(NIAPC),  INTENT(INOUT)  :: APC
        INTEGER, DIMENSION(NIAC+1), INTENT(IN) :: IAC
        INTEGER, DIMENSION(NNZC), INTENT(IN)   :: JAC
        INTEGER, DIMENSION(NIAC), INTENT(IN)   :: IUC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)     :: R
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT)  :: D
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: ic0, ic1, id0
        INTEGER :: jcol
        INTEGER :: j, n
        DOUBLEPRECISION :: t
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C         FORWARD SOLVE - APC * D = R
        FORWARD: DO n = 1, NIAC
          t   = R(n)
          ic0 = IAC(n) + 1
          ic1 = IUC(n) - 1
          LOWER: DO j = ic0, ic1
            jcol = JAC(j)
            t    = t - APC(j) * D(jcol)
          END DO LOWER
          D(n) = t
        END DO FORWARD
C         BACKWARD SOLVE - D = D / U
        BACKWARD: DO n = NIAC, 1, -1
          id0 = IAC(n)
          ic0 = IUC(n)
          ic1 = IAC(n+1)-1
          t   = D(n)
          UPPER: DO j = ic0, ic1
            jcol = JAC(j)
            t    = t - APC(j) * D(jcol)
          END DO UPPER
C           COMPUTE D FOR DIAGONAL - D = D / U
          D(n) = APC(id0) * t
        END DO BACKWARD
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGILU0A

C-------POLYNOMIAL PRECONDITIONER
      SUBROUTINE SUPCGGLSPOL(NOPT,NTRD,NTRDV,NNZC,NIAC,AC,IAC,JAC,
     2                       GLSPOLY)
        USE UPCGMODULE, ONLY: TGLSPOLY
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NOPT
        INTEGER, INTENT(IN) :: NTRD
        INTEGER, INTENT(IN) :: NTRDV
        INTEGER, INTENT(IN) :: NNZC
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NNZC),  INTENT(IN)     :: AC
        INTEGER, DIMENSION(NIAC+1), INTENT(IN)   :: IAC
        INTEGER, DIMENSION(NNZC), INTENT(IN)     :: JAC
        TYPE (TGLSPOLY), INTENT(INOUT)           :: GLSPOLY 
C     + + + LOCAL DEFINITIONS + + +
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C
C---------ESTIMATE MAXIMUM AND MINIMUM EIGENVALUES
C         ACTUAL MAXIMUM AND MINIMUM EIGENVALUES ARE CALCULATED IF
C         GLSPOLY%NLANSTEP EQUALS NIAC
        IF ( GLSPOLY%IEIGCALC.NE.0 ) THEN
          CALL UPCGLANCZOS(NOPT,NTRD,NTRDV,NIAC,NNZC,IAC,JAC,AC,
     2      GLSPOLY%NLANSTEP,GLSPOLY%NLAN2,GLSPOLY%D_LANCZOS,
     3      GLSPOLY%D_V1,GLSPOLY%D_V0,GLSPOLY%D_V,
     4      GLSPOLY%D_E,GLSPOLY%INTV)
        END IF
C
C---------CALCULATE ALPHA, BETA, AND GAMMA VALUES FOR POLYNOMIAL
C         PRECONDITIONER
        CALL UPCGUPDPOLY(GLSPOLY)        
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGGLSPOL


      SUBROUTINE SUPCGMV(NOPT,NTRD,NNZC,NIAC,A,D1,D2,IAC,JAC)
!        USE OMP_LIB
        USE UPCGMODULE, ONLY: DTMV
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NOPT
        INTEGER, INTENT(IN) :: NTRD
        INTEGER, INTENT(IN) :: NNZC
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NNZC),  INTENT(IN)    :: A
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: D1
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT) :: D2
        INTEGER, DIMENSION(NIAC+1), INTENT(IN) :: IAC
        INTEGER, DIMENSION(NNZC), INTENT(IN)   :: JAC
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: i, j
        INTEGER :: istart, iend
        INTEGER :: jstart, jend, jlen
        INTEGER :: jcol
        INTEGER :: n0
        INTEGER :: iblksize
        DOUBLEPRECISION :: tv
        DOUBLEPRECISION :: t
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C
C---------START TIMER
        CALL SUPCGTIMER(0,tv,DTMV)
C
        SELECT CASE ( NOPT )
          CASE (1,3)
            DO i = 1, NIAC
C               ADD DIAGONAL AND OFF-DIAGONAL TERMS
              t   = DZERO
              jstart = IAC(i)
              jend   = IAC(i+1)-1
              DO j = jstart, jend
                jcol = JAC(j) 
                t  = t + A(j) * D1(jcol)
              END DO
              D2(i) = t
            END DO
          CASE (2)
            n0 = ( ( NIAC - 1 ) / NTRD ) + 1
!$OMP  PARALLEL 
!$OMP& SHARED(NIAC,A,D1,D2)
!$OMP& PRIVATE(iblksize,istart,iend,jstart,jend,jlen)
!$OMP& NUM_THREADS(NTRD)
!$OMP  DO SCHEDULE(STATIC) 
            do i = 1, NTRD
              iblksize = min( n0, NIAC - ( i - 1 ) * n0 )
              istart   = ( i - 1 ) * n0 + 1
              iend     = istart + iblksize
              jstart   = IAC(istart)
              jend     = IAC(iend)
              jlen     = jend - jstart + 1
              IF ( iblksize.GT.0 ) THEN 
                CALL SUPCGSGEMV(iblksize,NIAC,jlen,jstart,
     2                          IAC(istart),JAC(jstart),
     3                          A(jstart),D1,D2(istart))
              END IF
            END DO
!$OMP  END DO NOWAIT
!$OMP  END PARALLEL
!!$OMP  PARALLEL DO
!!$OMP& NUM_THREADS(NTRD)
!!$OMP& DEFAULT(SHARED)
!!$OMP& PRIVATE(i, ic0, ic1, j, t)
!!!$OMP& PRIVATE(i, ic0, ic1, j)
!!!$OMP& REDUCTION(+:t)
!            DO i = 1, NIAC
!C               ADD DIAGONAL AND OFF-DIAGONAL TERMS
!              t   = 0.0D0
!              ic0 = IAC(i)
!              ic1 = IAC(i+1) - 1
!              DO j = ic0, ic1
!!                WRITE (*,1070) OMP_GET_THREAD_NUM(), i, JAC(j), 
!!     2                            j, ic0, ic1, t
!                t = t + A(j) * D1( JAC(j) )
!              END DO
!              D2(i) = t
!            END DO
!!$OMP  END PARALLEL DO
!!            WRITE (*,*)
        END SELECT
1070    FORMAT('TID:',I5,1X,'Row:',I5,1X,'Col:',I5,1X,
     2         'j:',I5,1X,'ic0:',I5,1X,'ic1',I5,1X't:'G15.7)
C
C         END TIMER
        CALL SUPCGTIMER(1,tv,DTMV)
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGMV

      SUBROUTINE SUPCGSGEMV(IBLKSIZE,NIAC,JLEN,JSTART,IA,JA,A,D1,D2)
        IMPLICIT NONE
C         + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: IBLKSIZE
        INTEGER, INTENT(IN) :: NIAC
        INTEGER, INTENT(IN) :: JLEN
        INTEGER, INTENT(IN) :: JSTART
        INTEGER, DIMENSION(IBLKSIZE+1),       INTENT(IN)    :: IA
        INTEGER, DIMENSION(JLEN),             INTENT(IN)    :: JA
        DOUBLEPRECISION, DIMENSION(JLEN),     INTENT(IN)    :: A
        DOUBLEPRECISION, DIMENSION(NIAC),     INTENT(IN)    :: D1
        DOUBLEPRECISION, DIMENSION(IBLKSIZE), INTENT(INOUT) :: D2
C         + + + LOCAL DEFINITIONS + + +
        INTEGER :: i, j
        INTEGER :: j0, j1, jcol
        DOUBLEPRECISION, PARAMETER :: dzero = 0.0d0
C         + + + FUNCTIONS + + +
C         + + + CODE + + +
        DO i = 1, IBLKSIZE
          j0 = IA(i)   - JSTART + 1
          j1 = IA(i+1) - JSTART
          D2(i) = dzero
          DO j = j0, j1
            jcol = JA(j)
            D2(i) = D2(i) + A(j) * D1(jcol)
          END DO
        END DO
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGSGEMV

      SUBROUTINE SUPCGAXPY(NOPT,NTRDV,NIAC,C,D1,D2)
        USE UPCGMODULE, ONLY: DTAXPY
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NOPT
        INTEGER, INTENT(IN) :: NTRDV
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION,  INTENT(IN)                     :: C
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: D1
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT) :: D2
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        DOUBLEPRECISION :: tv
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C---------START TIMER
        CALL SUPCGTIMER(0,tv,DTAXPY)
C
        SELECT CASE ( NOPT )
          CASE ( 1 )
            DO n = 1, NIAC
              D2(n) = D2(n) + C * D1(n)
            END DO
          CASE ( 2 )
!$OMP  PARALLEL
!$OMP& NUM_THREADS(NTRDV)
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(n)
!$OMP  DO SCHEDULE(STATIC)
            DO n = 1, NIAC
              D2(n) = D2(n) + C * D1(n)
            END DO
!$OMP  END DO NOWAIT
!$OMP  END PARALLEL
        END SELECT
C---------END TIMER
        CALL SUPCGTIMER(1,tv,DTAXPY)
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGAXPY

      SUBROUTINE SUPCGSETX(NOPT,NTRDV,NIAC,D1,C)
        USE UPCGMODULE, ONLY: DTMISC
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NOPT
        INTEGER, INTENT(IN) :: NTRDV
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT) :: D1
        DOUBLEPRECISION,  INTENT(IN)                     :: C
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        DOUBLEPRECISION :: tv
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C---------START TIMER
        CALL SUPCGTIMER(0,tv,DTMISC)
C
        SELECT CASE ( NOPT )
          CASE ( 1 )
            DO n = 1, NIAC
              D1(n) = C
            END DO
          CASE ( 2 )
!$OMP  PARALLEL
!$OMP& NUM_THREADS(NTRDV)
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(n)
!$OMP  DO SCHEDULE(STATIC)
            DO n = 1, NIAC
              D1(n) = C
            END DO
!$OMP  END DO NOWAIT
!$OMP  END PARALLEL
        END SELECT
C---------END TIMER
        CALL SUPCGTIMER(1,tv,DTMISC)
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGSETX

      SUBROUTINE SUPCGDCOPY(NOPT,NTRDV,NIAC,D1,D2)
        USE UPCGMODULE, ONLY: DTMISC
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NOPT
        INTEGER, INTENT(IN) :: NTRDV
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: D1
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT) :: D2
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        DOUBLEPRECISION :: tv
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C---------START TIMER
        CALL SUPCGTIMER(0,tv,DTMISC)
C
        SELECT CASE ( NOPT )
          CASE ( 1 )
            DO n = 1, NIAC
              D2(n) = D1(n)
            END DO
          CASE ( 2 )
!$OMP  PARALLEL
!$OMP& NUM_THREADS(NTRDV)
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(n)
!$OMP  DO SCHEDULE(STATIC)
            DO n = 1, NIAC
              D2(n) = D1(n)
            END DO
!$OMP  END DO NOWAIT
!$OMP  END PARALLEL
        END SELECT
C---------END TIMER
        CALL SUPCGTIMER(1,tv,DTMISC)
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGDCOPY

      SUBROUTINE SUPCGDSCAL(NOPT,NTRDV,NIAC,C,D1)
        USE UPCGMODULE, ONLY: DTMISC
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NOPT
        INTEGER, INTENT(IN) :: NTRDV
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, INTENT(IN) :: C
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT) :: D1
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        DOUBLEPRECISION :: tv
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C---------START TIMER
        CALL SUPCGTIMER(0,tv,DTMISC)
C
        SELECT CASE ( NOPT )
          CASE ( 1 )
            DO n = 1, NIAC
              D1(n) = C * D1(n)
            END DO
          CASE ( 2 )
!$OMP  PARALLEL
!$OMP& NUM_THREADS(NTRDV)
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(n)
!$OMP  DO SCHEDULE(STATIC)
           DO n = 1, NIAC
              D1(n) = C * D1(n)
            END DO
!$OMP  END DO NOWAIT
!$OMP  END PARALLEL
        END SELECT
C---------END TIMER
        CALL SUPCGTIMER(1,tv,DTMISC)
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGDSCAL

      SUBROUTINE SUPCGVVP(NOPT,NTRDV,NIAC,D1,D2,D3)
        USE UPCGMODULE, ONLY: DTVVP
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NOPT
        INTEGER, INTENT(IN) :: NTRDV
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: D1
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: D2
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(INOUT) :: D3
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        DOUBLEPRECISION :: tv
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C---------START TIMER
        CALL SUPCGTIMER(0,tv,DTVVP)
C
        SELECT CASE ( NOPT )
          CASE ( 1 )
            DO n = 1, NIAC
              D3(n) = D1(n) * D2(n)
            END DO
          CASE ( 2 )
!$OMP  PARALLEL
!$OMP& NUM_THREADS(NTRDV)
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(n)
!$OMP  DO SCHEDULE(STATIC)
            DO n = 1, NIAC
              D3(n) = D1(n) * D2(n)
            END DO
!$OMP  END DO NOWAIT
!$OMP  END PARALLEL
        END SELECT
C---------END TIMER
        CALL SUPCGTIMER(1,tv,DTVVP)
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGVVP

      DOUBLEPRECISION FUNCTION SUPCGDP(NOPT,NTRDV,NIAC,A,B) RESULT(C)
        USE UPCGMODULE, ONLY: DTDP
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NOPT
        INTEGER, INTENT(IN) :: NTRDV
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: A
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: B
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        DOUBLEPRECISION :: tv
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C---------START TIMER
        CALL SUPCGTIMER(0,tv,DTDP)
C
        C = DZERO
        SELECT CASE ( NOPT )
          CASE ( 1 )
            DO n = 1, NIAC
              C = C + A(n) * B(n)
            END DO
          CASE ( 2 )
!$OMP  PARALLEL
!$OMP& NUM_THREADS(NTRDV)
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(n)
!$OMP& REDUCTION(+: C)
!$OMP  DO SCHEDULE(STATIC)
            DO n = 1, NIAC
              C = C + A(n) * B(n)
            END DO
!$OMP  END DO NOWAIT
!$OMP  END PARALLEL
        END SELECT
C---------END TIMER
        CALL SUPCGTIMER(1,tv,DTDP)
C---------RETURN
        RETURN
      END FUNCTION SUPCGDP
C
C-------INFINITY NORM
      DOUBLEPRECISION FUNCTION SUPCGLINFNORM(NIAC,A) RESULT(B)
        USE UPCGMODULE, ONLY: DTMISC
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: NIAC
        DOUBLEPRECISION, DIMENSION(NIAC),  INTENT(IN)    :: A
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: n
        DOUBLEPRECISION :: tv
        DOUBLEPRECISION :: babs
        DOUBLEPRECISION, PARAMETER :: DZERO = 0.0D0
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
C---------START TIMER
        CALL SUPCGTIMER(0,tv,DTMISC)
C
        B    = DZERO
        babs = DZERO
        DO n = 1, NIAC
          IF ( ABS( A(n) ).GT.babs ) THEN
            B    = A(n)
            babs = ABS( B )
          END IF
        END DO
C---------END TIMER
        CALL SUPCGTIMER(1,tv,DTMISC)
C---------RETURN
        RETURN
      END FUNCTION SUPCGLINFNORM
C
C-------TIMER FOR UPCG CALCULATIONS
      SUBROUTINE SUPCGTIMER(It,T1,Ts)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: It
        DOUBLEPRECISION, INTENT(INOUT) :: T1
        DOUBLEPRECISION, INTENT(INOUT) :: Ts
C     + + + LOCAL DEFINITIONS + + +
        DOUBLEPRECISION :: tt
        DOUBLEPRECISION :: dt
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        IF ( It.EQ.0 ) THEN
          T1 = SECNDS(0.0)
        ELSE
          dt = SECNDS(T1)
          Ts = Ts + dt
        END IF
C---------RETURN
        RETURN
      END SUBROUTINE SUPCGTIMER
      
      
