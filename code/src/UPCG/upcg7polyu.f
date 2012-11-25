      SUBROUTINE UPCGUPDPOLY(GLSPOLY)
        USE GLOBAL,   ONLY:IOUT
        USE UPCGMODULE, ONLY: TGLSPOLY
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        TYPE (TGLSPOLY), INTENT(INOUT)           :: GLSPOLY 
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: i,j,k,m
        INTEGER :: ilen
        DOUBLEPRECISION :: bet1, bet, alp
        DOUBLEPRECISION, PARAMETER :: dzero = 0.0D0
        DOUBLEPRECISION, PARAMETER :: done  = 1.0D0
C     + + + FUNCTIONS + + +
        DOUBLEPRECISION :: DOTP
C     + + + CODE + + +
        m = 2
C---------INITIALIZE P0
        DO i = 1, GLSPOLY%NINTV
          GLSPOLY%P0(i) = done
          GLSPOLY%P0(GLSPOLY%NINTV+i) = dzero
        END DO
!C---------INITIALIZE PPOL, APPOL, AND QPOL
C        
        CALL XMUL(   m, GLSPOLY%NINTV, GLSPOLY%P0, 
     2                  GLSPOLY%NINTV, GLSPOLY%INTV, 
     3             m+1, GLSPOLY%NINTV, GLSPOLY%PPOL )
        
        bet1 = DOTP( m+1,GLSPOLY%NINTV,GLSPOLY%PPOL,
     2               m+1,GLSPOLY%NINTV,GLSPOLY%PPOL )
        bet1 = SQRT( bet1 )
        
        DO i = 1, ((m+1)*GLSPOLY%NINTV)
          GLSPOLY%PPOL(i) = GLSPOLY%PPOL(i) / bet1
        END DO
        
        GLSPOLY%GAMMA(1) = DOTP( m+1,GLSPOLY%NINTV,GLSPOLY%PPOL,
     2                             m,GLSPOLY%NINTV,GLSPOLY%P0 )
        bet = dzero
        GLSPOLY%BETA(1) = bet1
C        
C---------CALCULATE ALPHA, BETA, AND GAMMA FOR EACH
C         POLYNOMIAL DEGREE
        j = m
        DO i = 1, GLSPOLY%NDEGREE
          call XMUL( j+1,GLSPOLY%NINTV,GLSPOLY%PPOL,
     2                   GLSPOLY%NINTV,GLSPOLY%INTV,
     3               j+2,GLSPOLY%NINTV,GLSPOLY%APPOL )
          alp = DOTP( j+2,GLSPOLY%NINTV,GLSPOLY%APPOL,
     2                j+1,GLSPOLY%NINTV,GLSPOLY%PPOL )
          GLSPOLY%ALPHA(i) = alp
          CALL POLSUM( j+2,GLSPOLY%NINTV,GLSPOLY%APPOL,
     2                 j+1,GLSPOLY%NINTV,GLSPOLY%PPOL,-alp )
          
          IF ( i.GT.1 ) THEN
            CALL POLSUM( j+2,GLSPOLY%NINTV,GLSPOLY%APPOL,
     2                     j,GLSPOLY%NINTV,GLSPOLY%QPOL,-bet )
          END IF
          
          bet = DOTP( j+2,GLSPOLY%NINTV,GLSPOLY%APPOL,
     2                j+2,GLSPOLY%NINTV,GLSPOLY%APPOL )
          bet = SQRT( bet )
          
          GLSPOLY%BETA(i+1) = bet
          
          IF ( bet.EQ.dzero ) THEN
            WRITE (IOUT,'(//,A)') 'UPCGUPDPOLY: bet = 0.0'
            CALL USTOP('UPCGUPDPOLY: bet = 0.0')
          END IF
          
C-----------COPY PPOL TO QPOL
          ilen = ( j+1 ) * GLSPOLY%NINTV
          CALL DCOPY( ilen, GLSPOLY%PPOL, 1, GLSPOLY%QPOL, 1 )
          
          DO k = 1, ( (j+2) * GLSPOLY%NINTV )
            GLSPOLY%PPOL(k) = GLSPOLY%APPOL(k) / bet
          END DO
          
          j = j + 1
          GLSPOLY%GAMMA(i+1) = DOTP( j+1,GLSPOLY%NINTV,GLSPOLY%PPOL,
     2                                 m,GLSPOLY%NINTV,GLSPOLY%P0 )
          
        END DO
C---------RETURN
        RETURN
      END SUBROUTINE UPCGUPDPOLY


      SUBROUTINE XMUL(P1, P2, P, NINTV, INTV, Q1, Q2, Q)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: P1
        INTEGER, INTENT(IN) :: P2
        DOUBLEPRECISION, DIMENSION(P1*P2), INTENT(IN) :: P
        INTEGER, INTENT(IN) :: NINTV
        DOUBLEPRECISION, DIMENSION(2), INTENT(IN) :: INTV
        INTEGER, INTENT(IN) :: Q1
        INTEGER, INTENT(IN) :: Q2
        DOUBLEPRECISION, DIMENSION(Q1*Q2), INTENT(INOUT) :: Q
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: i, j
        INTEGER :: jq, jp
        DOUBLEPRECISION :: c, h, h2
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        IF ( P1.EQ.0 ) GOTO 9999
        
        DO i = 1, NINTV
          c = 0.5D0 * ( INTV(i+1) + INTV(i) )
          h = 0.5D0 * ( INTV(i+1) - INTV(i) )
          
          DO j = 1, P1
            jq = (j - 1) * Q2 + i
            jp = (j - 1) * P2 + i
            Q(jq) = c * P(jp)
          END DO
          
          Q(P1*Q2+i) = 0.0D0
          Q(Q2+i) = Q(Q2+i) + h * P(i)
          h2 = 0.5D0 * h
          
          DO j = 2, P1
            jq = j * Q2 + i
            jp = ( j - 1 ) * P2  + i
            Q(jq) = Q(jq) + h2 * P(jp)
            jq = ( j - 2 ) * Q2 + i
            Q(jq) = Q(jq) + h2 * P(jp)
          END DO
          
        END DO
        
09999   RETURN  
      END SUBROUTINE XMUL

      DOUBLEPRECISION FUNCTION DOTP(P1,P2,P,Q1,Q2,Q) RESULT(R)
        USE UPCGMODULE, ONLY: PI, SCAL
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: P1
        INTEGER, INTENT(IN) :: P2
        DOUBLEPRECISION, DIMENSION(P1*P2), INTENT(IN) :: P
        INTEGER, INTENT(IN) :: Q1
        INTEGER, INTENT(IN) :: Q2
        DOUBLEPRECISION, DIMENSION(Q1*Q2), INTENT(IN) :: Q
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: i, j, m
        INTEGER :: nintv
        INTEGER :: jq, jp
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        R = 0.0D0
        m = MIN( Q1, P1 )
        nintv = Q2
        
        IF ( m.NE.0 ) THEN
          DO i = 1, nintv
            R = R + SCAL(i) * Q(i) * P(i)
          END DO
        END IF
        
        DO j = 1, m
          DO i = 1, nintv
            jq = (j - 1) * Q2 + i
            jp = (j - 1) * P2 + i
            R = R + SCAL(i) * Q(jq) * P(jp)
          END DO
        END DO

        R = R * PI / 2.0D0
C---------RETURN
09999    RETURN  
      END FUNCTION DOTP

      SUBROUTINE POLSUM(P1,P2,P,Q1,Q2,Q,GAM)
        IMPLICIT NONE
C     + + + DUMMY ARGUMENTS + + +
        INTEGER, INTENT(IN) :: P1
        INTEGER, INTENT(IN) :: P2
        DOUBLEPRECISION, DIMENSION(P1*P2), INTENT(INOUT) :: P
        INTEGER, INTENT(IN) :: Q1
        INTEGER, INTENT(IN) :: Q2
        DOUBLEPRECISION, DIMENSION(Q1*Q2), INTENT(IN) :: Q
        DOUBLEPRECISION, INTENT(IN) :: GAM
C     + + + LOCAL DEFINITIONS + + +
        INTEGER :: i, j
        INTEGER :: jq, jp
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
        DO j = 1, Q1
          DO i = 1, P2
            jp = (j - 1) * P2 + i
            jq = (j - 1) * Q2 + i
            P(jp) = P(jp) + GAM * Q(jq)
          END DO
        END DO
C---------RETURN
09999   RETURN  
      END SUBROUTINE POLSUM

C     + + + DUMMY ARGUMENTS + + +
C     + + + LOCAL DEFINITIONS + + +
C     + + + FUNCTIONS + + +
C     + + + CODE + + +
