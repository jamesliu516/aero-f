C
C
C
      SUBROUTINE ERIEMANNGW(DL,UL,PL,DR,UR,PR,PM,UM,RIL,RIR,ALPHA,BETA,
     &                      PREF,gam)
*
*----------------------------------------------------------------------*
*                                                                      *
C     Exact Riemann Solver for the Time-Dependent                      *
C     One Dimensional Euler Equations                                  *
*                                                                      *
*----------------------------------------------------------------------*
*
      IMPLICIT NONE
*
C     Declaration of variables:
*
      INTEGER I, NRITER
*
      REAL*8  gam, G1, G2, G3, G4, G5, G6, G7, G8, G9,
     &        DL, UL, PL, CL, DR, UR, PR, CR,
     &        DIAPH1, DOMLEN, DS, DX, PM, PSCALE, PS, S,
     &        TIMEOU, UM, US, XPOS, ALPHA, BETA, RIL, RIR
      REAL*8  QUSER, CUV, PPV, PMIN, PMAX, PRATIO
      REAL*8  QMAX, CUP
      REAL*8  PQ, PTL, PTR, GEL, GER
      REAL*8  POLD, UDIFF, FL, FLD, AK, BK, QRT
      REAL*8  FR, FRD, P1, P2, CHANGE, PREF
      REAL*8  DP1, DP2, PR1, PR2
      REAL*8  DRP, F1, F2, F3
      REAL*8  TOLPRE
*
      TOLPRE  = 1.E-6
      NRITER  = 100

      IF (PL.LT.0.0 .OR. PR.LT.0.0) WRITE(*,*) 'NEGATIVE PRESSURES',
     &                              DL,UL,PL,DR,UR,PR
*
C     Compute gamma related constants
*
C     gam  = 1.4
      G1 = (gam - 1.0)/(2.0*gam)
      G2 = (gam + 1.0)/(2.0*gam)
      G3 = 2.0*gam/(gam - 1.0)
      G4 = 2.0/(gam - 1.0)
      G5 = 2.0/(gam + 1.0)
      G6 = (gam - 1.0)/(gam + 1.0)
      G7 = (gam - 1.0)/2.0
      G8 = gam - 1.0
      G9 = gam/(gam  -1.0)
*
C     CONSTANTS FOR BAROTROPIC LIQUID
*
C     ALPHA  = 1.3986E-01
C     BETA   = 7.15
C     PREF   = 3.14E-02
      PREF   =-PREF
*
C     Compute sound speeds
*
      CL = DSQRT(gam*PL/DL)
      CR = DSQRT(ALPHA*BETA*DR**(BETA-1.0))
*
      QUSER = 2.0
*
C     Compute guess pressure from PVRS Riemann solver
*
      CUP  = 0.25*(DL + DR)*(CL + CR)
      PPV  = 0.5*(PL + PR) + 0.5*(UL - UR)*CUP
      PPV  = MAX(0.0, PPV)
      PMIN = MIN(PL,  PR)
      PMAX = MAX(PL,  PR)
      QMAX = PMAX/PMIN
*
      IF(QMAX.LE.QUSER.AND.
     & (PMIN.LE.PPV.AND.PPV.LE.PMAX))THEN
*
C        Select PVRS Riemann solver
*
         PM = PPV
      ELSE
         IF(PPV.LT.PMIN)THEN
*
C           Select Two-Rarefaction Riemann solver
*
            PQ  = (PL/PR)**G1
            US  = (PQ*UL/CL + UR/CR +
     &            G4*(PQ - 1.0))/(PQ/CL + 1.0/CR)
            PTL = 1.0 + G7*(UL - US)/CL
            PTR = 1.0 + G7*(US - UR)/CR
            PM  = 0.5*(PL*PTL**G3 + PR*PTR**G3)
         ELSE
*
C           Select Two-Shock Riemann solver with
C           PVRS as estimate
*
            GEL = SQRT((G5/DL)/(G6*PL + PPV))
            GER = SQRT((G5/DR)/(G6*PR + PPV))
            PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER)
         ENDIF
      ENDIF
*

*
      PM    = 0.5*(PR  +PL)
      POLD  = PM
      UDIFF = UR - UL
*
*
      DO 10 I = 1, NRITER
*
C        GAS MEDIUM RELATIONS
*
         IF(POLD.LE.PL)THEN
*
C        Rarefaction wave
*
           PRATIO = POLD/PL
           FL     = G4*CL*(PRATIO**G1 - 1.0)
           FLD    = (1.0/(DL*CL))*PRATIO**(-G2)
         ELSE
*
C        Shock wave
*
           AK   = G5/DL
           BK   = G6*PL
           QRT  = SQRT(AK/(BK + POLD))
           FL   = (POLD - PL)*QRT
           FLD  = (1.0 - 0.5*(POLD - PL)/(BK + POLD))*QRT
           
         ENDIF
*
C        WATER MEDIUM RELATIONS
*
         IF(POLD.LE.PR)THEN
*
C        Rarefaction wave
*
           PRATIO = (POLD  +PREF)/(PR  +PREF)
           FR   = 2.0*CR*(PRATIO**((BETA  -1.0)/(2.0*BETA)) -1.0)/
     &            (BETA  -1.0)
           PR1  = (POLD  +PREF)**(-(BETA +1)/(2*BETA))
           PR2  = (PR    +PREF)**( (BETA -1)/(2*BETA))
           FRD  = CR*PR1/(PR2 *BETA)
         ELSE
*
C        Shock wave
*
           P1  = ((POLD  +PREF)/ALPHA)**(1/BETA)
           P2  = ((PR    +PREF)/ALPHA)**(1/BETA)
           F1  = ((POLD  -PR)*(P1  -P2)/(DR*P1))
           FR  = F1**0.5
           DRP = (((POLD  +PREF)/ALPHA)**(1/BETA  -1))/(ALPHA*BETA)
           FRD = (P1*((POLD  -PR)*DRP  + (P1  -P2))  -
     &                (POLD  -PR)*(P1  -P2)*DRP)/P1**2
           FRD = FRD*0.5/SQRT(DR)
         ENDIF

         PM     = POLD - (FL + FR + UDIFF)/(FLD + FRD)
         CHANGE = 2.0*ABS((PM- POLD)/(PM+ POLD))
         IF(CHANGE.LE.TOLPRE)GOTO 20
         IF(PM.LT.0.0)PM = TOLPRE
         POLD  = PM
*
 10   CONTINUE
*
       WRITE(6,*)'Divergence in Newton-Raphson iteration', PM, PL, PR
*
 20   CONTINUE

      UM = 0.5*(UL + UR + FR - FL)

      IF (PM.LE.PL) THEN
         RIL  = DL*(PM/PL)**(1.0/gam)
      ELSE
         RIL  = DL*(G9*PM  - 0.5*(PM  -PL))/(G9*PL  +0.5*(PM -PL))
      ENDIF

      IF (PM.LE.PR) THEN
         RIR  = DR*((PM  +PREF)/(PR  +PREF))**(1/BETA)
      ELSE
         RIR  = DR*((PM  +PREF)/(PR  +PREF))**(1/BETA)
C        RIR  = DR*UR/UM
      ENDIF
*
*
      END
