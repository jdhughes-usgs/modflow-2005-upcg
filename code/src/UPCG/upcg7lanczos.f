      SUBROUTINE UPCGLANCZOS(Nopt,Ntrd,Ntrdv,N,Nnz,Ia,Ja,AC,
     2                       Nsteps,N2,Vv,D_w,D_u,D,E,Ev)
        USE GLOBAL,   ONLY:IOUT
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: Nopt
        INTEGER, INTENT(IN) :: Ntrd
        INTEGER, INTENT(IN) :: Ntrdv
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: Nnz
        INTEGER, INTENT(IN), DIMENSION(N+1) :: Ia
        INTEGER, INTENT(IN), DIMENSION(NNZ) :: Ja
        DOUBLEPRECISION, INTENT(IN), DIMENSION(NNZ) :: AC
        INTEGER, INTENT(IN) :: Nsteps
        INTEGER, INTENT(IN) :: N2
        DOUBLEPRECISION, INTENT(INOUT), DIMENSION((Nsteps+1)*N2) :: Vv
        DOUBLEPRECISION, INTENT(INOUT), DIMENSION(N) :: D_w
        DOUBLEPRECISION, INTENT(INOUT), DIMENSION(N) :: D_u
        DOUBLEPRECISION, INTENT(INOUT), DIMENSION(N) :: D
        DOUBLEPRECISION, INTENT(INOUT), DIMENSION(N) :: E
        DOUBLEPRECISION, INTENT(INOUT), DIMENSION(2) :: Ev
C     + + + LOCAL DEFINITIONS + + +
        CHARACTER (LEN=1) :: compz
        INTEGER :: i
        INTEGER :: ipos, ipos2
        INTEGER :: info
        INTEGER, DIMENSION(1) :: seed
        DOUBLEPRECISION :: harvest
        DOUBLEPRECISION :: t
        DOUBLEPRECISION :: alpha
        DOUBLEPRECISION :: beta
        DOUBLEPRECISION :: orthtol
        DOUBLEPRECISION :: wn
        DOUBLEPRECISION :: sg
        DOUBLEPRECISION, PARAMETER :: dzero=0.0D0
        DOUBLEPRECISION, PARAMETER :: done =1.0D0
C     + + + FUNCTIONS + + +
        DOUBLEPRECISION :: DNRM2
        DOUBLEPRECISION :: SUPCGDP
C     + + + CODE + + +
C---------INITIALIZE Vv, D_w, D_u, D, E       
        CALL SUPCGSETX(Nopt,Ntrdv,(Nsteps+1)*N2,Vv,dzero)
        CALL SUPCGSETX(Nopt,Ntrdv,N,D_w,dzero)
        CALL SUPCGSETX(Nopt,Ntrdv,N,D_u,dzero)
        CALL SUPCGSETX(Nopt,Ntrdv,N,D,dzero)
        CALL SUPCGSETX(Nopt,Ntrdv,N,E,dzero)
        
        seed(1) = 100
        CALL RANDOM_SEED (PUT = seed)
        DO i = 1, N
          CALL RANDOM_NUMBER (harvest)
          Vv(i) = harvest
        END DO

C---------CALCULATE EUCLIDIAN NORM  
        t = DNRM2(N, VV(1), 1)
!        CALL DSCAL(N, 1.0/t, Vv(1), 1)
        CALL SUPCGDSCAL(Nopt, Ntrdv, N, 1.0/t, Vv(1))
        beta    = dzero
        orthtol = 1.0D-8
        wn      = dzero
C
C---------TRIDIAGONAL MATRIX
C         SAVED IN D(DIAG), E(OFF-DIAG)
C---------MAIN LOOP
!        WRITE (*,*) 'LANCZOS ALG. BEGINS...'
        LANCZOSL: DO i = 1, Nsteps
          ipos = (i-1)*N2 + 1
          CALL SUPCGMV(Nopt,Ntrd,Nnz,N,AC,Vv(ipos),D_w,Ia,Ja)
          IF ( i.EQ.1 ) THEN
            CALL SUPCGAXPY(Nopt,Ntrdv,N,-beta,Vv(ipos),D_w)
          ELSE
            CALL SUPCGAXPY(Nopt,Ntrdv,N,-beta,Vv(ipos-N),D_w)
          END IF
          alpha = SUPCGDP(Nopt,Ntrdv,N,D_w,Vv(ipos))
          wn    = wn + alpha * alpha
          D(i)  = alpha
          CALL SUPCGAXPY(Nopt,Ntrdv,N,-alpha,Vv(ipos),D_w)
          ! U = V'*W
          CALL DGEMV('T',N,i,done,Vv,N2,D_w,1,dzero,D_u,1)
          ! W = V*U - W
          CALL DGEMV('N',N,i,-done,Vv,N2,D_u,1,done,D_w,1)

          beta = SUPCGDP(Nopt,Ntrdv,N,D_w,D_w)
          IF ( beta*REAL(i,8).LT.orthtol*wn ) THEN
            EXIT LANCZOSL
          END IF
          wn = wn + 2.0 * beta
          beta = SQRT(beta);
          
          ipos2 = i*N2 + 1
          CALL SUPCGAXPY(Nopt,Ntrdv,N,done/beta,D_w,Vv(ipos2))

          IF (i.LT.Nsteps) THEN
            E(i) = beta
          END IF

        END DO LANCZOSL
!        WRITE (*,*) 'LANCZOS ALG. ENDS.'

C-------COMPUTE EIGENVALUE OF T BY LAPACK
!        WRITE (*,*) 'LAPACK STEQF BEGINS ...'
        CALL DSTERF( Nsteps, D, E, info )
        IF ( INFO.NE.0 ) THEN
          WRITE (IOUT,'(//,A)') 'LAPACK: FAILED TO FIND EIGENVALUES !!!'
          CALL USTOP('LAPACK: FAILED TO FIND EIGENVALUES !!!')
        END IF
!        WRITE (*,*) 'LAPACK STERF ENDS.'
C
        Ev(2) = D(Nsteps)
        Ev(1) = D(1)

!        WRITE (*,2000) Ev(1), Ev(2)
!2000  FORMAT('EIGEN-VALUE:[',E10.4,',',E10.4,']')
C
C---------RETURN
        RETURN
      END SUBROUTINE UPCGLANCZOS