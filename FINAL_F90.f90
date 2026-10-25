!===============================================================================
!
! PROGRAM FINAL -- FORTRAN 90 VERSION
!
! SOURCE: J. C. NEWMAN, JR., "FINITE-ELEMENT ANALYSIS OF FATIGUE
! CRACK PROPAGATION -- INCLUDING THE EFFECTS OF CRACK CLOSURE",
! PH.D. THESIS, VIRGINIA POLYTECHNIC INSTITUTE AND STATE UNIVERSITY,
! 1974.  APPENDIX E, COMPUTER PROGRAM LISTING (PP. 102-121).
!
! 32-BIT REAL IS NOT ADEQUATE FOR THIS PROGRAM (THE RIGID-SPRING
! STIFFNESS SK2 = 1.0E7*CE IS ADDED TO AND LATER REMOVED FROM THE
! FACTORED MATRIX), SO DOUBLE PRECISION IS USED.
!
!===============================================================================

!*** DIMENSIONS ENLARGED (PUBLISHED: 400 ELEMENTS, 250 NODES, 500 DOF,
!*** 12000 STIFFNESS WORDS, 6 BOUNDARIES, 50 NODES PER BOUNDARY,
!*** 8 ELEMENTS PER NODE).
!***
!***   MAXEL  - ELEMENTS (AFTER SUBDIVISION)
!***   MAXNOD - NODES
!***   MAXDOF - DEGREES OF FREEDOM (2*MAXNOD)
!***   MAXA   - STIFFNESS MATRIX STORAGE (VARIABLE-BANDWIDTH PROFILE)
!***   MAXBND - BOUNDARIES
!***   MAXBN  - NODES ON ONE BOUNDARY
!***   MAXEPN - ELEMENTS CONNECTED TO ONE NODE

!===============================================================================

MODULE C_MAIN

  USE ISO_FORTRAN_ENV

  IMPLICIT NONE

  INTEGER,PARAMETER :: MAXEL  = 4000
  INTEGER,PARAMETER :: MAXNOD = 2500
  INTEGER,PARAMETER :: MAXDOF = 2*MAXNOD
  INTEGER,PARAMETER :: MAXA   = 2000000
  INTEGER,PARAMETER :: MAXBND = 10
  INTEGER,PARAMETER :: MAXBN  = 500
  INTEGER,PARAMETER :: MAXEPN = 12

  REAL(REAL64), DIMENSION(MAXA)          :: AA
  REAL(REAL64)                           :: AM,CE,CMU,CRACK,PMAX,SIGYS,SK2,SNET,SYLD,WIDTH
  REAL(REAL64), DIMENSION(MAXEL)         :: AREA,SIGE,SIGX,SIGY,STRX,STRXY,STRY,TAU
  REAL(REAL64), DIMENSION(MAXDOF,1)      :: BB
  INTEGER,      DIMENSION(MAXNOD)        :: IIMAX
  INTEGER                                :: KLU,KMAX,KNEW,LMAX,LPRIT,NBC,NDOF,NELE,NEP,NODE
  INTEGER,      DIMENSION(MAXBND,MAXBN)  :: LSAME
  INTEGER,      DIMENSION(MAXDOF)        :: MB,MPTAB,MSUM
  INTEGER,      DIMENSION(MAXEL)         :: MPLAS
  INTEGER,      DIMENSION(MAXEL,3)       :: MSAME
  INTEGER,      DIMENSION(MAXBND)        :: NOD
  INTEGER,      DIMENSION(MAXNOD,MAXEPN) :: NSAME
  REAL(REAL64), DIMENSION(MAXDOF)        :: Q,R,T1,T2,T3,Z
  REAL(REAL64), DIMENSION(MAXNOD)        :: U,UOLD,V,V2,VOLD,XR,YR
  REAL(REAL64), DIMENSION(MAXEL,3)       :: X,Y

END MODULE C_MAIN

!===============================================================================

    PROGRAM FINAL

    USE C_MAIN
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE

    ! LOCAL VARIABLE DECLARATIONS

    INTEGER                         :: I,ICARD,IERR,IFAC,II,IM,IMAX,IU,IV,IX
    INTEGER                         :: J,JJ,JU,JV
    INTEGER                         :: K,KCHK,KK,KKK,KU,KV
    INTEGER                         :: L,LDOF,LL,LO,LX
    INTEGER                         :: M,M1,MAX,MAXDIF,MM,MM1,MN,MPS,MX
    INTEGER                         :: N,N1,N2,N3,NB,NBK,NC,ND,NDIFF,NDIV,NTRI
    INTEGER                         :: NU,NUIU,NUIV,NUJU,NUJV,NUKU,NUKV,NUM,NUNU
    INTEGER                         :: NV,NVIU,NVIV,NVJU,NVJV,NVKU,NVKV,NVNV,NX,NY
    INTEGER,DIMENSION(MAXBND)       :: NBX,NBY
    REAL(REAL64)                    :: ALP,ALP2,ALP3,ALP5,ALP6,BETA1,BETA2,BETA3
    REAL(REAL64)                    :: CRACK1,CX,CY,G,QQ,RES,RES1,RS,RX,S,S1,S2,S3,SBK
    REAL(REAL64),DIMENSION(MAXBND)  :: B,C,D,E,SBX,SBY,XMAX,XMIN,YMAX,YMIN
    REAL(REAL64),DIMENSION(MAXEL,3) :: XINT,YINT

    EXTERNAL DIMERR,PLAS,SYMBAN

    ! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !
    ! CONFIGURATION AND LOADING MUST BE SYMMETRIC ABOUT X-AXIS
    !
    ! MINIMUM ELEMENT SIZE IS 0.001 UNITS
    ! CRACK = INITIAL CRACK LENGTH (MUST LIE ON X-AXIS AND BE MEASURED FROM X=0)
    ! WIDTH = WIDTH OF PLATE ALONG X-AXIS
    ! THICKNESS - PLATE HAS UNIT THICKNESS
    ! CE  = MODULUS OF ELASTICITY
    ! CMU = POISSON RATIO
    ! MPS = 1  PLANE STRAIN , = 0  PLANE STRESS
    ! IIMAX(N) = NUMBER OF ELEMENTS CONNECTED TO NODE N
    ! NSAME(N,I) = ELEMENT NUMBER CONNECTED TO NODE N (I=1,IIMAX(N))
    ! MSAME(M,I) = NODE NUMBERS FOR ELEMENT M  (I=1,3)
    ! LMAX = TOTAL NUMBER OF NODES
    ! IMAX = TOTAL NUMBER OF INITIAL ELEMENTS
    ! NOD(I) = NUMBER OF NODES ON BOUNDARY I
    ! LSAME(I,J) = NODE NUMBERS ON BOUNDARY I (J=1,NOD(I))
    ! AREA(I) = AREA OF ELEMENT I
    ! XR(N),YR(N) = COORDINATES FOR NODE N
    ! X(I,J),Y(I,J) = NODAL COORDINATES FOR ELEMENT I (J=1,3)
    ! LPRIT = 0  NO INTERNAL PRINTOUT, = 1  INTERNAL PRINTOUT
    ! KMAX = TOTAL NUMBER OF ELEMENTS AFTER SUBDIVISIONS
    ! NDIV = NUMBER OF SUBDIVISIONS OF THE ELEMENT MESH
    ! NB = TOTAL NUMBER OF BOUNDARIES (MUST BE DESCRIBED IN ANALYTICAL FORM
    !      EITHER STRAIGHT OR CIRCULAR BOUNDARIES)
    ! XMAX(I),XMIN(I),YMAX(I),YMIN(I) = MAX. AND MIN. VALUE OF X AND Y
    !                                   COORDINATES FOR BOUNDARY I
    ! B(I),C(I),D(I),E(I) ARE THE COEFFICIENTS OF THE ANALYTICAL
    !                     EQUATION FOR BOUNDARY I
    !    STRAIGHT BOUNDARY -  Y=B*X+C  WITH D=0 AND E=0
    !    FOR STRAIGHT BOUNDARY PARALLEL TO Y-AXIS USE X=C WITH D=0 AND E=1
    !    CIRCULAR BOUNDARY - (X-B)**2+(Y-C)**2=D**2   WITH E=0
    ! NBC = BOUNDARY NUMBER CONTAINING CRACK
    ! BOUNDARY CONDITIONS FOR BOUNDARY I ARE INPUT IN NBX(I) AND NBY(I)
    !       NBX(I) OR NBY(I) = 0   FREE
    !                        = 1   FIXED
    !                        = 2   LOADED
    ! LOADED BOUNDARIES HAVE UNIFORM STRESS
    !    SBX(I) OR SBY(I) = UNIT APPLIED STRESS ON BOUNDARY I
    !                       ACTING IN X- OR Y-DIRECTION
    ! TOTAL LOAD ACTING IN Y-DIRECTION (NORMAL TO CRACK PLANE) IS
    ! CONSIDERED THE PRIMARY LOADING AND MUST NOT BE ZERO.
    ! NBX(I) = BOUNDARY CONDITION IN X-DIRECTION FOR BOUNDARY I
    ! NBY(I) = BOUNDARY CONDITION IN Y-DIRECTION FOR BOUNDARY I
    ! FOR BOUNDARY ON X-AXIS CONTAINING THE CRACK SET NBX(I)=0 AND NBY(I)=1
    !     (PROGRAM WILL FREE NODES ALONG SPECIFIED CRACK LENGTH)
    ! NEP = 0  ELASTIC ANALYSIS,  NEP GREATER THAN 0 INDICATES PLASTIC
    !          ANALYSIS WITH PRINT OUT AT NEP INCREMENTS
    ! NODAL COORDINATES X(I,J) AND Y(I,J) OF ELEMENT I ARE INPUT IN COUNTER
    ! CLOCKWISE DIRECTION STARTING WITH ANY NODE (J=1,2,3) FOR THAT ELEMENT.

    ! INPUT DATA

    READ(5,*) LPRIT,IMAX,NDIV,NB,NBC,NEP
    IF(NB.GT.MAXBND) CALL DIMERR('BOUNDARIES (MAXBND)',NB,MAXBND)
    KCHK = IMAX*4**NDIV
    IF(KCHK.GT.MAXEL) CALL DIMERR('ELEMENTS (MAXEL)',KCHK,MAXEL)
    READ(5,*) CRACK,WIDTH
    DO I = 1,NB
      READ(5,*) B(I),C(I),D(I),E(I)
      READ(5,*) XMAX(I),XMIN(I),YMAX(I),YMIN(I)
      READ(5,*) NBX(I),NBY(I)
      READ(5,*) SBX(I),SBY(I)
    ENDDO
    READ(5,*) CE,CMU
    READ(5,*) MPS
    DO ICARD = 1,IMAX
      READ(5,*) I,X(I,1),X(I,2),X(I,3),Y(I,1),Y(I,2),Y(I,3)
      IF(LPRIT.NE.0) WRITE(6,'(I5,6F10.4)') I,X(I,1),X(I,2),X(I,3),Y(I,1),Y(I,2),Y(I,3)
    ENDDO

    ! SUBDIVISION
    !
    ! SUBDIVIDES MESH USING THE INITIAL ELEMENT NUMBERING SCHEME IN RENUMBERING
    ! ELEMENTS IN ORDER TO MAINTAIN REASONABLE BANDWIDTHS IN STIFFNESS MATRIX
    ! (THE INITIAL ELEMENT NUMBERING SCHEME SHOULD BE NUMBERED TO HAVE
    ! MINIMUM BANDWIDTHS)

    DO M = 1,NB
      XMAX(M) = XMAX(M)+0.0001
      XMIN(M) = XMIN(M)-0.0001
      YMAX(M) = YMAX(M)+0.0001
      YMIN(M) = YMIN(M)-0.0001
    ENDDO
    ND = -1
10  CONTINUE
    ND = ND+1
    IF(LPRIT.NE.0) WRITE(6,'(A10,A7,2A16)') 'TRIANGLE','NODE','X','Y'
    I = 0
    KMAX = IMAX*(4**ND)
    IF(ND.EQ.0) GOTO 35
    K = IMAX*(4**(ND-1))
    KKK = K
20  CONTINUE
    I = I+1
    J = 0
25  CONTINUE
    J = J+1
    K = K+1
    II = I
    JJ = J+1
    KK = K+1
    IF(J.EQ.3)II = K-2
    IF(J.EQ.3)JJ = J-2
    IF(J.EQ.3)KK = K-2
    X(K,J) = X(I,J)
    X(I,J) = 0.5*(X(I,J)+X(II,JJ))
    X(K,JJ) = X(I,J)
    X(KK,J) = X(I,J)
    Y(K,J) = Y(I,J)
    Y(I,J) = 0.5*(Y(I,J)+Y(II,JJ))
    Y(K,JJ) = Y(I,J)
    Y(KK,J) = Y(I,J)
    M = 0
30  CONTINUE
    M = M+1
    IF(M.GT.NB) THEN
      IF(J.EQ.3) THEN
        IF(I.EQ.KKK) THEN
          ! RENUMBER NODES AND CALCULATION OF AREAS
          IX = IMAX*(4**(ND-1))
          M1 = IX
          DO M = 1,IX
            MM = 4*M-3
            DO I = 1,3
              XINT(MM,I) = X(M,I)
              YINT(MM,I) = Y(M,I)
            ENDDO
            DO K = 1,3
              MM1 = MM+K
              M1 = M1+1
              DO I = 1,3
                XINT(MM1,I) = X(M1,I)
                YINT(MM1,I) = Y(M1,I)
              ENDDO
            ENDDO
          ENDDO
          DO M = 1,KMAX
            DO I = 1,3
              X(M,I) = XINT(M,I)
              Y(M,I) = YINT(M,I)
            ENDDO
          ENDDO
          GOTO 35
        ENDIF
        GOTO 20
      ENDIF
      GOTO 25
    ENDIF
    IF(D(M).LT.0.0001) GOTO 30
    RES = SQRT(((X(K,J)-B(M))**2+(Y(K,J)-C(M))**2-(D(M))**2)**2)
    IF(RES.GT.0.0001) GOTO 30
    RES1 = SQRT(((X(II,JJ)-B(M))**2+(Y(II,JJ)-C(M))**2-(D(M))**2)**2)
    IF(RES1.GT.0.0001) GOTO 30
    ! CHECK IF COORDINATES ARE OUTSIDE THE BOUNDS
    IF(X(K,J)  .LT.XMIN(M) .OR. X(K,J)  .GT.XMAX(M)) GOTO 30
    IF(Y(K,J)  .LT.YMIN(M) .OR. Y(K,J)  .GT.YMAX(M)) GOTO 30
    IF(X(II,JJ).LT.XMIN(M) .OR. X(II,JJ).GT.XMAX(M)) GOTO 30
    IF(Y(II,JJ).LT.YMIN(M) .OR. Y(II,JJ).GT.YMAX(M)) GOTO 30
    RS = SQRT((X(I,J)-B(M))**2+(Y(I,J)-C(M))**2)
    X(I,J) = B(M)+((X(I,J)-B(M))*D(M)/RS)
    Y(I,J) = C(M)+((Y(I,J)-C(M))*D(M)/RS)
    X(K,JJ) = X(I,J)
    X(KK,J) = X(I,J)
    Y(K,JJ) = Y(I,J)
    Y(KK,J) = Y(I,J)
    GOTO 30
35  CONTINUE
    IF(LPRIT.NE.0) THEN
      DO I = 1,KMAX
        DO J = 1,3
          WRITE(6,'(I10,I7,2E16.7)') I,J,X(I,J),Y(I,J)
        ENDDO
      ENDDO
    ENDIF
    I = 0
    L = 0
45  CONTINUE
    I = I+1
    IF(I.GT.KMAX) THEN
      LMAX = L
      IF(LPRIT.EQ.0) GOTO 50
      WRITE(6,'(A8,2A16)') 'NODE','X','Y'
      DO N = 1,LMAX
        WRITE(6,'(I8,2E16.7)')N,XR(N),YR(N)
      ENDDO
      WRITE(6,'(A10,A16)') 'TRIANGLE','AREA'
      DO I = 1,KMAX
        WRITE(6,'(I10,E16.7)') I,AREA(I)
      ENDDO
     ! BOUNDARY-NODE IDENTIFICATION
50    CONTINUE
      N = 0
      LMAX = L
      GOTO 70
    ENDIF
    J = 0
55  CONTINUE
    L = L+1
    IF(L.GT.MAXNOD) CALL DIMERR('NODES (MAXNOD)',L,MAXNOD)
    J = J+1
    XR(L) = X(I,J)
    YR(L) = Y(I,J)
    IF(I.EQ.1) GOTO 65
    LL = L-1
    DO K = 1,LL
      QQ = ABS(XR(L)-XR(K))
      IF(QQ.LT.0.0001) THEN
        QQ = ABS(YR(L)-YR(K))
        IF(QQ.LT.0.0001) GOTO 60
      ENDIF
    ENDDO
    GOTO 65
60  CONTINUE
    L = L-1
65  CONTINUE
    IF(J.EQ.3) THEN
      S1 = SQRT((X(I,J-2)-X(I,J-1))**2+(Y(I,J-2)-Y(I,J-1))**2)
      S2 = SQRT((X(I,J-1)-X(I,J))**2+(Y(I,J-1)-Y(I,J))**2)
      S3 = SQRT((X(I,J)-X(I,J-2))**2+(Y(I,J)-Y(I,J-2))**2)
      S = 0.5*(S1+S2+S3)
      AREA(I) = SQRT(S*(S-S1)*(S-S2)*(S-S3))
    ENDIF
    IF(J.EQ.3) GOTO 45
    GOTO 55
70  CONTINUE
    N = N+1
    IF(N.GT.NB) THEN
      IF(ND.LT.NDIV) GOTO 10
      ! ELEMENT-NODE IDENTIFICATION
      DO NODE = 1,LMAX
        I = 0
        DO NTRI = 1,KMAX
          J = 0
75        CONTINUE
          J = J+1
          CX = ABS(XR(NODE)-X(NTRI,J))
          CY = ABS(YR(NODE)-Y(NTRI,J))
          IF(CX.LT.0.001) THEN
            IF(CY.LT.0.001) THEN
              I = I+1
              IF(I.GT.MAXEPN) CALL DIMERR('ELEMENTS AT A NODE (MAXEPN)',I,MAXEPN)
              NSAME(NODE,I) = NTRI
              IIMAX(NODE) = I
              MSAME(NTRI,J) = NODE
              CYCLE
            ENDIF
            IF(J.EQ.3) CYCLE
            GOTO 75
          ENDIF
          IF(J.EQ.3) CYCLE
          GOTO 75
        ENDDO
      ENDDO
      IF(LPRIT.NE.0) THEN
        WRITE(6,'(A10,A7,A9)') 'NODE','I','NSAME'
        DO NODE = 1,LMAX
          II = IIMAX(NODE)
          DO I = 1,II
            WRITE(6,'(I10,I7,I9)') NODE,I,NSAME(NODE,I)
          ENDDO
        ENDDO
        WRITE(6,'(A10,A7,A9)') 'NTRI','J','MSAME'
        DO NTRI = 1,KMAX
          DO J = 1,3
            WRITE(6,'(I10,I7,I9)') NTRI,J,MSAME(NTRI,J)
          ENDDO
        ENDDO
      ENDIF

      ! MATERIAL CONSTANTS

80    CONTINUE
      IF(MPS.EQ.0) THEN
        ! MINIMIZATION OF POTENTIAL ENERGY - PLANE STRESS
        G = CE/(1.0-CMU**2)
        BETA1 = 1.0
        !*** CORRECTED: EQUATIONS FOR BETA2 AND BETA3 HAVE BEEN
        !*** SWITCHED IN ORDER TO ENSURE FORMULATION IS CORRECT
        !*** FOR VALUES OF CMU .NE. 1/3.
        BETA2 = CMU
        BETA3 = (1.0-CMU)/2.0
      ELSE
        ! MINIMIZATION OF POTENTIAL ENERGY - PLANE STRAIN
        G = CE/(2.0*(1.0+CMU)*(1.0-2.0*CMU))
        BETA1 = 2.0*(1.0-CMU)
        BETA2 = 1.0-2.0*CMU
        BETA3 = 2.0*CMU
      ENDIF

      ! DETERMINE BAND-WIDTHS FOR STIFFNESS MATRIX

85    CONTINUE
      DO NODE = 1,LMAX
        MAXDIF = 0
        IM = IIMAX(NODE)
        DO M = 1,IM
          NTRI = NSAME(NODE,M)
          DO L = 1,3
            NUM = MSAME(NTRI,L)
            NDIFF = 2*IABS(NUM-NODE)
            IF(NDIFF.GT.MAXDIF) MAXDIF = NDIFF
          ENDDO
        ENDDO
        NU = 2*NODE-1
        NV = NU+1
        MB(NU) = MAXDIF+1
        IF(MB(NU).GT.NU)MB(NU) = NU
        MB(NV) = MB(NU)+1
        IF(LPRIT.NE.0) WRITE(6,'(10X,5I5)') NODE,NU,MB(NU),NV,MB(NV)
      ENDDO
      MSUM(1) = 0
      MSUM(2) = 1
      LO = 2*LMAX
      DO I = 3,LO
        II = I-1
        MSUM(I) = MSUM(II)+MB(II)
      ENDDO
      LDOF = MSUM(LO)+MB(LO)
      IF(LDOF.GT.MAXA) CALL DIMERR('STIFFNESS STORAGE (MAXA)',LDOF,MAXA)
      WRITE(6,'(A,I0)') 'STORAGE REQUIREMENT FOR STIFFNESS MATRIX = ',LDOF

      ! ASSEMBLE STIFFNESS MATRIX

      LL = 2*LMAX
      DO I = 1,LL
        MPTAB(I) = 0
        R(I) = 0.0
        BB(I,1) = 0.0
      ENDDO
      DO L = 1,LDOF
        AA(L) = 0.0
      ENDDO
      DO NODE = 1,LMAX
        IM = IIMAX(NODE)
        M = 0
90      CONTINUE
        M = M+1
        IF(M.GT.IM)CYCLE
        NTRI = NSAME(NODE,M)
        L = 0
95      CONTINUE
        L = L+1
        CX = ABS(XR(NODE)-X(NTRI,L))
        CY = ABS(YR(NODE)-Y(NTRI,L))
        IF(CX.LT.0.001) THEN
          IF(CY.LT.0.001) THEN
            IF(L.EQ.1) I = MSAME(NTRI,1)
            IF(L.EQ.1) J = MSAME(NTRI,2)
            IF(L.EQ.1) K = MSAME(NTRI,3)
            IF(L.EQ.2) I = MSAME(NTRI,2)
            IF(L.EQ.2) J = MSAME(NTRI,3)
            IF(L.EQ.2) K = MSAME(NTRI,1)
            IF(L.EQ.3) I = MSAME(NTRI,3)
            IF(L.EQ.3) J = MSAME(NTRI,1)
            IF(L.EQ.3) K = MSAME(NTRI,2)
            GOTO 100
          ENDIF
          GOTO 95
        ENDIF
        IF(L.EQ.3) GOTO 100
        GOTO 95
100     CONTINUE
        N = NODE
        IF(LPRIT.NE.0) WRITE(6,'(8(3X,I5))') N,IM,M,NTRI,L,I,J,K
        NU = 2*N-1
        NV = NU+1
        IU = 2*I-1
        IV = IU+1
        JU = 2*J-1
        JV = JU+1
        KU = 2*K-1
        KV = KU+1
        NUIU = MSUM(NU)-NU+MB(NU)+IU
        IF(IU.GT.NU) GOTO 110
        AA(NUIU) = ((BETA1*((YR(J)-YR(K))**2)+BETA2*((XR(K)-XR(J))**2))*G/4.0*1.0/AREA(NTRI))+AA(NUIU)
110     CONTINUE
        NUJU = MSUM(NU)-NU+MB(NU)+JU
        IF(JU.GT.NU) GOTO 115
        AA(NUJU) = ((BETA1*(YR(K)-YR(I))*(YR(J)-YR(K))+BETA2*(XR(I)-XR(K))*(XR(K)-XR(J)))*G/4.0*1.0/AREA(NTRI))+AA(NUJU)
115     CONTINUE
        NUKU = MSUM(NU)-NU+MB(NU)+KU
        IF(KU.GT.NU) GOTO 120
        AA(NUKU) = ((BETA1*(YR(I)-YR(J))*(YR(J)-YR(K))+BETA2*(XR(J)-XR(I))*(XR(K)-XR(J)))*G/4.0*1.0/AREA(NTRI))+AA(NUKU)
120     CONTINUE
        NVIU = MSUM(NV)-NV+MB(NV)+IU
        IF(IU.GT.NV) GOTO 125
        AA(NVIU) = ((BETA3*(XR(K)-XR(J))*(YR(J)-YR(K))+BETA2*(YR(J)-YR(K))*(XR(K)-XR(J)))*G/4.0*1.0/AREA(NTRI))+AA(NVIU)
125     CONTINUE
        NVJU = MSUM(NV)-NV+MB(NV)+JU
        IF(JU.GT.NV) GOTO 130
        AA(NVJU) = ((BETA3*(XR(I)-XR(K))*(YR(J)-YR(K))+BETA2*(YR(K)-YR(I))*(XR(K)-XR(J)))*G/4.0*1.0/AREA(NTRI))+AA(NVJU)
130     CONTINUE
        NVKU = MSUM(NV)-NV+MB(NV)+KU
        IF(KU.GT.NV) GOTO 135
        AA(NVKU) = ((BETA3*(XR(J)-XR(I))*(YR(J)-YR(K))+BETA2*(YR(I)-YR(J))*(XR(K)-XR(J)))*G/4.0*1.0/AREA(NTRI))+AA(NVKU)
135     CONTINUE
        NVIV = MSUM(NV)-NV+MB(NV)+IV
        IF(IV.GT.NV) GOTO 140
        AA(NVIV) = ((BETA1*((XR(K)-XR(J))**2)+BETA2*((YR(J)-YR(K))**2))*G/4.0*1.0/AREA(NTRI))+AA(NVIV)
140     CONTINUE
        NVJV = MSUM(NV)-NV+MB(NV)+JV
        IF(JV.GT.NV) GOTO 145
        AA(NVJV) = ((BETA1*(XR(I)-XR(K))*(XR(K)-XR(J))+BETA2*(YR(K)-YR(I))*(YR(J)-YR(K)))*G/4.0*1.0/AREA(NTRI))+AA(NVJV)
145     CONTINUE
        NVKV = MSUM(NV)-NV+MB(NV)+KV
        IF(KV.GT.NV) GOTO 150
        AA(NVKV) = ((BETA1*(XR(J)-XR(I))*(XR(K)-XR(J))+BETA2*(YR(I)-YR(J))*(YR(J)-YR(K)))*G/4.0*1.0/AREA(NTRI))+AA(NVKV)
150     CONTINUE
        NUIV = MSUM(NU)-NU+MB(NU)+IV
        IF(IV.GT.NU) GOTO 155
        AA(NUIV) = ((BETA3*(YR(J)-YR(K))*(XR(K)-XR(J))+BETA2*(XR(K)-XR(J))*(YR(J)-YR(K)))*G/4.0*1.0/AREA(NTRI))+AA(NUIV)
155     CONTINUE
        NUJV = MSUM(NU)-NU+MB(NU)+JV
        IF(JV.GT.NU) GOTO 160
        AA(NUJV) = ((BETA3*(YR(K)-YR(I))*(XR(K)-XR(J))+BETA2*(XR(I)-XR(K))*(YR(J)-YR(K)))*G/4.0*1.0/AREA(NTRI))+AA(NUJV)
160     CONTINUE
        NUKV = MSUM(NU)-NU+MB(NU)+KV
        IF(KV.GT.NU) GOTO 90
        AA(NUKV) = ((BETA3*(YR(I)-YR(J))*(XR(K)-XR(J))+BETA2*(XR(J)-XR(I))*(YR(J)-YR(K)))*G/4.0*1.0/AREA(NTRI))+AA(NUKV)
        GOTO 90
      ENDDO
      ! PRINT OUT STIFFNESS MATRIX IF REQUIRED
      IF(LPRIT.NE.0) THEN
        DO L = 1,LDOF
          WRITE(6,'(5X,I8,4X,E16.7)') L,AA(L)
        ENDDO
      ENDIF
      GOTO 185
    ENDIF
    NOD(N) = 0
    L = 0
165 CONTINUE
    L = L+1
    IF(D(N).LT.0.0001) THEN
      IF(E(N).GT.0.0001) THEN
        RES = SQRT((XR(L)-B(N)*YR(L)-C(N))**2)
        GOTO 170
      ENDIF
      RES = SQRT((YR(L)-B(N)*XR(L)-C(N))**2)
170   CONTINUE
      IF(RES.LT.0.0001) GOTO 175
      GOTO 180
    ENDIF
    RES1 = SQRT(((XR(L)-B(N))**2+(YR(L)-C(N))**2-(D(N))**2)**2)
    IF(RES1.LT.0.0001) GOTO 175
    GOTO 180
175 CONTINUE
    IF(XR(L).LT.XMIN(N)) GOTO 180
    IF(XR(L).GT.XMAX(N)) GOTO 180
    IF(YR(L).LT.YMIN(N)) GOTO 180
    IF(YR(L).GT.YMAX(N)) GOTO 180
    NOD(N) = NOD(N)+1
    IF(NOD(N).GT.MAXBN) CALL DIMERR('NODES ON A BOUNDARY (MAXBN)',NOD(N),MAXBN)
    I = NOD(N)
    LSAME(N,I) = L
    IF(LPRIT.NE.0) WRITE(6,'(5X,I3,4X,I5,4X,I5)') N,I,LSAME(N,I)
180 CONTINUE
    IF(L.EQ.LMAX) GOTO 70
    GOTO 165

    ! LOADING (UNIFORM STRESS ACTING IN X- OR Y-DIRECTION ON BOUNDARY I)

185 CONTINUE
    DO N = 1,2
      DO K = 1,NB
        IF(N.EQ.1 .AND. NBX(K).NE.2)CYCLE
        IF(N.EQ.2 .AND. NBY(K).NE.2)CYCLE
        SBK = SBX(K)
        IF(N.EQ.2)SBK = SBY(K)
        NBK = NOD(K)
        DO I = 1,NBK
          L = LSAME(K,I)
          NX = 2*L-1
          IF(N.EQ.2)NX = NX+1
          IX = IIMAX(L)
          DO M = 1,IX
            MN = NSAME(L,M)
            DO J = 1,3
              MX = MSAME(MN,J)
              DO II = 1,NBK
                LX = LSAME(K,II)
                IF(MX.NE.LX)CYCLE
                RX = SQRT((XR(L)-XR(MX))**2+(YR(L)-YR(MX))**2)
                BB(NX,1) = BB(NX,1)+SBK*RX/2.0
                R(NX) = BB(NX,1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    ! PRINT OUT APPLIED LOAD VECTOR IF REQUIRED

    IF(LPRIT.NE.0) THEN
      WRITE(6,'(7X,"I",8X,"BB(I,1)")')
      DO I = 1,LL
        WRITE(6,'(5X,I5,E16.7)') I,BB(I,1)
      ENDDO
    ENDIF

    ! FIXED BOUNDARY CONDITIONS ARE SATISFIED BY RIGID SPRINGS

    SK2 = 1.0E+07*CE
    DO K = 1,NB
      IF(NBX(K).NE.1)CYCLE
      NX = NOD(K)
      DO I = 1,NX
        L = LSAME(K,I)
        NU = 2*L-1
        NUNU = MSUM(NU)+MB(NU)
        MPTAB(NU) = 1
        AA(NUNU) = AA(NUNU)+SK2
      ENDDO
    ENDDO
    CRACK1 = CRACK-0.001
    DO K = 1,NB
      IF(NBY(K).NE.1) CYCLE
      NY = NOD(K)
      DO I = 1,NY
        L = LSAME(K,I)
        IF(K.EQ.NBC .AND. XR(L).LT.CRACK1) CYCLE
        NV = 2*L
        NVNV = MSUM(NV)+MB(NV)
        MPTAB(NV) = 1
        AA(NVNV) = AA(NVNV)+SK2
      ENDDO
    ENDDO

    ! MATRIX SOLUTION ROUTINE

    IFAC = 0
    CALL SYMBAN(MAXA,LL,MB,MSUM,AA,1,BB,IFAC,T1,IERR,ALP,Z,T2,T3,NC)
    WRITE(6,'(A,I0)') 'SYMBAN IERR = ',IERR

    ! PRINT OUT UNIT LOAD DISPLACEMENTS AND STRESSES

    DO N = 1,LMAX
      NU = 2*N-1
      NV = NU+1
      U(N) = BB(NU,1)
      V(N) = BB(NV,1)
    ENDDO
    WRITE(6,'(/,A)') 'UNIT LOAD DISPLACEMENTS AND STRESSES'
    WRITE(6,'(/,A10,4A16)') 'NODE','X','Y','U','V'
    DO I = 1,LMAX
      WRITE(6,'(I10,1P,4E16.7)') I,XR(I),YR(I),U(I),V(I)
    ENDDO
    DO NTRI = 1,KMAX
      N1 = MSAME(NTRI,1)
      N2 = MSAME(NTRI,2)
      N3 = MSAME(NTRI,3)
      ALP2 = (U(N1)*(YR(N2)-YR(N3))+U(N2)*(YR(N3)-YR(N1))+U(N3)*(YR(N1)-YR(N2)))/(2.0*AREA(NTRI))
      ALP5 = (V(N1)*(YR(N2)-YR(N3))+V(N2)*(YR(N3)-YR(N1))+V(N3)*(YR(N1)-YR(N2)))/(2.0*AREA(NTRI))
      ALP3 = (U(N1)*(XR(N3)-XR(N2))+U(N2)*(XR(N1)-XR(N3))+U(N3)*(XR(N2)-XR(N1)))/(2.0*AREA(NTRI))
      ALP6 = (V(N1)*(XR(N3)-XR(N2))+V(N2)*(XR(N1)-XR(N3))+V(N3)*(XR(N2)-XR(N1)))/(2.0*AREA(NTRI))
      STRX(NTRI) = ALP2
      STRY(NTRI) = ALP6
      STRXY(NTRI) = ALP3+ALP5
      SIGX(NTRI) = (CE/(1.0-CMU**2))*(ALP2+CMU*ALP6)
      SIGY(NTRI) = (CE/(1.0-CMU**2))*(ALP6+CMU*ALP2)
      TAU(NTRI) = (CE/(2.0*(1.0+CMU)))*(ALP3+ALP5)
    ENDDO
    WRITE(6,'(/,A10,3A16,3A8)') 'TRIANGLE','SIGX','SIGY','TAU','I','J','K'
    SYLD = 0.0
    MAX = 0
    DO M = 1,KMAX
      I = MSAME(M,1)
      J = MSAME(M,2)
      K = MSAME(M,3)
      SIGE(M) = SQRT(SIGX(M)**2+SIGY(M)**2-SIGX(M)*SIGY(M)+3*TAU(M)**2)
      IF(SIGE(M).GT.SYLD) MAX = M
      IF(SIGE(M).GT.SYLD) SYLD = SIGE(M)
      WRITE(6,'(I10,1P,3E16.7,3I8)') M,SIGX(M),SIGY(M),TAU(M),I,J,K
    ENDDO
    WRITE(6,'()')
    WRITE(6,'("MODULUS OF ELASTICITY           =",1PE11.4," KSI")') CE
    WRITE(6,'("POISSON''S RATIO                 =",1PE11.4)'      ) CMU
    IF(NEP.EQ.0) STOP
    READ(5,*) SIGYS,AM
    WRITE(6,'("TENSILE YIELD STRESS            =",1PE11.4," KSI")') SIGYS
    WRITE(6,'("STRAIN HARDENING COEFFICIENT    =",1PE11.4)'       ) AM
    SYLD = SIGYS/11.4
    WRITE(6,'("APPLIED STRESS AT INITIAL YIELD =",1PE11.4," KSI AT ELEMENT ",I0/)') SYLD,MAX
    IF(MPS.EQ.1 .AND. NEP.GT.0) THEN
      WRITE(6,'(/,"PROGRAM CANNOT EXECUTE PLANE STRAIN PLASTICITY")')
      STOP
    ENDIF

    CALL PLAS

    END PROGRAM FINAL

!===============================================================================

    SUBROUTINE PLAS

    ! PERFORMS AN INCREMENTAL PLASTICITY ANALYSIS FOR AN ELASTIC-
    ! PERFECTLY PLASTIC MATERIAL

    USE C_MAIN
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE

    ! LOCAL VARIABLE DECLARATIONS

    REAL(REAL64)      :: A,A11,A12,A13,A21,A22,A23,A31,A32,A33
    REAL(REAL64)      :: ALP,ALP2,ALP3,ALP5,ALP6,AMIN,AR,ARD
    REAL(REAL64)      :: B,B1,B2,B3,C,C1,C2,C3
    REAL(REAL64)      :: DELP,DISCR,DSIGX,DSIGY,DSTRX,DSTRXY,DSTRY,DTAU
    REAL(REAL64)      :: ERIT,FX,FXY,FY,G,PCT,PONE,PT,PTWO,PX,PYLD,RE,REAS,RO
    REAL(REAL64)      :: S,S1,S2,SIGO,SIGXP,SIGYP,SX,SY,TAUP,TXY
    CHARACTER(4),SAVE :: GROW,HALT
    INTEGER           :: I,ICON,IERR,IFAC,ITER,M,MAXIT,MP
    INTEGER           :: N,N1,N2,N3,NC,NELE1,NELE2,NODE1,NODE2
    INTEGER           :: NPL,NU,NU1,NU2,NU3,NV,NV1,NV2,NV3
    CHARACTER(4)      :: WORD
    LOGICAL           :: CONVERGED
    INTEGER           :: NRHS

    EXTERNAL BREAK,CONTACT,PLOUT,SYMBAN

    DATA HALT/'HALT'/
    DATA GROW/'GROW'/

    ! SIGYS = TENSILE YIELD STRESS
    ! AM    = STRAIN HARDENING COEFFICIENT
    !         (SLOPE OF EFFECTIVE STRESS-STRAIN CURVE IN PLASTIC RANGE
    !         DIVIDED BY THE MODULUS OF ELASTICITY)
    !         MUST BE SET EQUAL TO 0 IN PRESENT PROGRAM
    ! PT    = TOTAL LOAD ACTING IN Y-DIRECTION
    ! WIDTH = WIDTH OF PLATE ALONG X-AXIS
    ! S     = DESIRED NOMINAL STRESS
    ! SNET  = NOMINAL STRESS (= PT/WIDTH)
    ! PCT   = PERCENTAGE USED TO COMPUTE INCREMENTAL LOAD
    ! MAXIT = MAXIMUM NUMBER OF ITERATIONS ALLOWED
    ! ERIT  = ERROR CHECK ON CONVERGENCE OF DISPLACEMENTS
    ! G     = RELAXATION PARAMETER (IMPROVES CONVERGENCE)
    !         G = 1.0  NORMAL CONVERGENCE

    ! INITIALIZATIONS

    NODE = LMAX
    NDOF = 2*NODE
    NELE = KMAX
    SIGO = SIGYS
    PT = 0.0
    DO I = 1,NDOF
      Q(I) = 0.0
      Z(I) = 0.0
      R(I) = SYLD*R(I)
    ENDDO
    DO I = 1,NODE
      NV = 2*I
      PT = PT+R(NV)
      PYLD = PT
      UOLD(I) = U(I)*SYLD
      VOLD(I) = V(I)*SYLD
    ENDDO
    DO M = 1,NELE
      MPLAS(M) = 0
      SIGX(M) = SIGX(M)*SYLD
      SIGY(M) = SIGY(M)*SYLD
      TAU(M) = TAU(M)*SYLD
      STRX(M) = STRX(M)*SYLD
      STRY(M) = STRY(M)*SYLD
      STRXY(M) = STRXY(M)*SYLD
    ENDDO
    READ(5,*) G
    READ(5,*) PCT,ERIT,MAXIT,NODE1,NODE2,NELE1,NELE2
    READ(5,'(E10.3,1X,A4)') S,WORD
    WRITE(6,'(A,F9.3,A,3X,A,1P,E10.3,3X,A,E10.3,3X,A,I0)') &
      'S = ',S,' KSI','PCT = ',PCT,'ERIT = ',ERIT,'MAXIT = ',MAXIT
    DELP = PCT*PT

    ! INCREMENT LOAD

    NPL = 0
    MP = 0
10  CONTINUE
    PMAX = PT
    PT = PT+DELP
    NPL = NPL+1
    KLU = 0
    ICON = 0
    DO M = 1,NELE
      MPLAS(M) = 0
    ENDDO
    DO I = 1,NODE
      V2(I) = VOLD(I)
    ENDDO

    ! HOLDING APPLIED LOAD CONSTANT - ITERATE UNTIL DISPLACEMENTS CONVERGE

15  CONTINUE
    ITERATE_LOOP: DO ITER = 1,MAXIT
      DO I = 1,NDOF
        BB(I,1) = R(I)*PT/PYLD+Q(I)
        IF(LPRIT.EQ.0)CYCLE
        WRITE(6,'(3E10.3)') R(I),Q(I),BB(I,1)
      ENDDO

      ! SOLVE FOR NEW DISPLACEMENTS

      IFAC = 1
      NRHS = 1
      CALL SYMBAN(MAXA,NDOF,MB,MSUM,AA,NRHS,BB,IFAC,T1,IERR,ALP,Z,T2,T3,NC)

      DO N = 1,NODE
        NU = 2*N-1
        NV = NU+1
        U(N) = BB(NU,1)
        V(N) = BB(NV,1)
        U(N) = U(N)-UOLD(N)
        V(N) = V(N)-VOLD(N)
        IF(LPRIT.NE.0) WRITE(6,'(10X,I5,2E10.3)') ITER,U(N),V(N)
        UOLD(N) = BB(NU,1)
        VOLD(N) = BB(NV,1)
      ENDDO

      ! COMPUTE TOTAL STRAIN INCREMENTS FROM DISPLACEMENT INCREMENTS
      ! COMPUTE ELASTIC STRESS INCREMENTS AND ADD TO CURRENT STRESSES
      ! CHECK YIELD CONDITION FOR PLASTIC ELEMENTS

      DO M = 1,NELE
        N1 = MSAME(M,1)
        N2 = MSAME(M,2)
        N3 = MSAME(M,3)
        ARD = 0.5/AREA(M)
        B1 = YR(N2)-YR(N3)
        C1 = XR(N3)-XR(N2)
        B2 = YR(N3)-YR(N1)
        C2 = XR(N1)-XR(N3)
        B3 = YR(N1)-YR(N2)
        C3 = XR(N2)-XR(N1)
        ALP2 = (U(N1)*B1+U(N2)*B2+U(N3)*B3)*ARD
        ALP5 = (V(N1)*B1+V(N2)*B2+V(N3)*B3)*ARD
        ALP3 = (U(N1)*C1+U(N2)*C2+U(N3)*C3)*ARD
        ALP6 = (V(N1)*C1+V(N2)*C2+V(N3)*C3)*ARD
        STRX(M) = STRX(M)+ALP2
        STRY(M) = STRY(M)+ALP6
        STRXY(M) = STRXY(M)+ALP3+ALP5
        DSTRX = ALP2
        DSTRY = ALP6
        DSTRXY = ALP3+ALP5
        DSIGX = (CE/(1.0-CMU**2))*(ALP2+CMU*ALP6)
        DSIGY = (CE/(1.0-CMU**2))*(ALP6+CMU*ALP2)
        DTAU = (CE/(2.0*(1.0+CMU)))*(ALP3+ALP5)
        SX = SIGX(M)
        SY = SIGY(M)
        TXY = TAU(M)
        S1 = SQRT(SX**2+SY**2-SX*SY+3.0*TXY**2)
        SIGX(M) = SX+DSIGX
        SIGY(M) = SY+DSIGY
        TAU(M) = TXY+DTAU
        IF(ICON.EQ.1) CYCLE
        S2 = SQRT(SIGX(M)**2+SIGY(M)**2-SIGX(M)*SIGY(M)+3.0*TAU(M)**2)
        IF(ABS(S2).LT.ABS(S1) .OR. ABS(S2).LE.SIGO) CYCLE
        MPLAS(M) = M
        MP = 1

        ! COMPUTE THE PLASTIC LOAD VECTOR FOR ELEMENT M

        A = DSIGX**2+DSIGY**2-DSIGX*DSIGY+3.0*DTAU**2
        B = SX*(2.0*DSIGX-DSIGY)+SY*(2.0*DSIGY-DSIGX)+6.0*TXY*DTAU
        C = S1**2-SIGO**2
        AMIN = 1.0E-10
        IF(A.LT.AMIN) CYCLE
        DISCR = B**2-4.0*A*C
        IF(DISCR.LT.0.0) THEN
          WRITE(6,'("NEGATIVE DISCRIMINANT",5X,"A =",1P,E10.3,5X,"B =",E10.3,5X,"C =",E10.3)') &
            A,B,C
          PX = 0.0
        ELSE
          DISCR = SQRT(DISCR)
          PONE = (-B+DISCR)/(2.0*A)
          PTWO = (-B-DISCR)/(2.0*A)
          PX = PONE
          IF(ABS(PONE).GT.ABS(PTWO)) PX = PTWO
          IF(LPRIT.NE.0) WRITE(6,'(10X,"PONE = ",1P,E10.3,5X,"PTWO = ",E10.3)') PONE,PTWO
        ENDIF
        SX = SX+PX*DSIGX
        SY = SY+PX*DSIGY
        TXY = TXY+PX*DTAU
        S1 = SQRT(SX**2+SY**2-SX*SY+3.0*TXY**2)
        DSTRX = (1.0-PX)*DSTRX
        DSTRY = (1.0-PX)*DSTRY
        DSTRXY = (1.0-PX)*DSTRXY
        ! CONVERGENCE IS IMPROVED BY USE OF THE RELAXATION PARAMETER G
        DSTRX = G*DSTRX
        DSTRY = G*DSTRY
        DSTRXY = G*DSTRXY
        FX = (SX-0.5*SY)/SIGO
        FY = (SY-0.5*SX)/SIGO
        FXY = 3.0*TXY/SIGO
        RO = FX**2+2.0*CMU*FX*FY+FY**2+(1.-CMU)*0.5*FXY**2+(1.0-CMU**2)*AM
        A11 = FX**2+2.0*CMU*FX*FY+CMU**2*FY**2
        A12 = FX*FY+CMU*(FX**2+FY**2)
        A21 = A12
        A13 = (1.0-CMU)*0.5*FXY*(FX+CMU*FY)
        A31 = A13
        A22 = FY**2+2.0*CMU*FX*FY+CMU**2*FX**2
        A23 = (1.0-CMU)*0.5*FXY*(FY+CMU*FX)
        A32 = A23
        A33 = ((1.0-CMU)*0.5*FXY)**2
        RE = CE/(RO*(1.0-CMU**2))
        SIGXP = RE*(A11*DSTRX+A12*DSTRY+A13*DSTRXY)
        SIGYP = RE*(A21*DSTRX+A22*DSTRY+A23*DSTRXY)
        TAUP = RE*(A31*DSTRX+A32*DSTRY+A33*DSTRXY)
        SIGX(M) = SIGX(M)-SIGXP
        SIGY(M) = SIGY(M)-SIGYP
        TAU(M) = TAU(M)-TAUP
        SIGE(M) = SQRT(SIGX(M)**2+SIGY(M)**2-SIGX(M)*SIGY(M)+3.*TAU(M)**2)
        NU1 = 2*N1-1
        NV1 = NU1+1
        NU2 = 2*N2-1
        NV2 = NU2+1
        NU3 = 2*N3-1
        NV3 = NU3+1
        AR = 0.5
        Q(NU1) = Q(NU1)+AR*(B1*SIGXP+C1*TAUP)
        Q(NU2) = Q(NU2)+AR*(B2*SIGXP+C2*TAUP)
        Q(NU3) = Q(NU3)+AR*(B3*SIGXP+C3*TAUP)
        Q(NV1) = Q(NV1)+AR*(C1*SIGYP+B1*TAUP)
        Q(NV2) = Q(NV2)+AR*(C2*SIGYP+B2*TAUP)
        Q(NV3) = Q(NV3)+AR*(C3*SIGYP+B3*TAUP)
        IF(LPRIT.NE.0) THEN
          WRITE(6,'(5X,I5,3X,1P,3E10.3,3X,2E10.3,3X,2E10.3,3X,2E10.3)') &
            M,S1,S2,SIGE(M),Q(NU1),Q(NV1),Q(NU2),Q(NV2),Q(NU3),Q(NV3)
        ENDIF
      ENDDO

      ! CHECK CONVERGENCE OF DISPLACEMENTS

      IF(ICON.EQ.1 .OR. MP.EQ.0) EXIT ITERATE_LOOP
      CONVERGED = ALL(ABS(U(1:NODE)).LE.ERIT*ABS(UOLD(1:NODE))) .AND. &
                  ALL(ABS(V(1:NODE)).LE.ERIT*ABS(VOLD(1:NODE)))
      IF(CONVERGED) THEN
        WRITE(6,'()')
        WRITE(6,'("DISPLACEMENTS CONVERGED IN ",I6," ITERATIONS")') ITER
        ICON = 1
      ENDIF
      IF(ITER.EQ.MAXIT) THEN
        WRITE(6,'()')
        WRITE(6,'("NO CONVERGENCE IN ",I6," ITERATIONS")') ITER
        CALL PLOUT(NODE1,NODE2,NELE1,NELE2)
        RETURN
      ENDIF

    ENDDO ITERATE_LOOP

    MP = 0
    ICON = 0
    SNET = PT/WIDTH

    ! CHECK CRACK DISPLACEMENTS FOR CLOSURE OR OPENING

    CALL CONTACT
    IF(KNEW.EQ.1) GOTO 15

    IF(SNET.GE.S .AND. DELP.GT.0.0) CALL PLOUT(NODE1,NODE2,NELE1,NELE2)
    IF(SNET.LE.S .AND. DELP.LT.0.0) CALL PLOUT(NODE1,NODE2,NELE1,NELE2)
    IF(NPL.EQ.NEP) CALL PLOUT(NODE1,NODE2,NELE1,NELE2)
    IF(NPL.EQ.NEP) NPL = 0
    IF(KLU.EQ.1) GOTO 35
    IF(SNET.LT.S .AND. DELP.GT.0.0) GOTO 10
    IF(SNET.GT.S .AND. DELP.LT.0.0) GOTO 10

    ! CHECK FOR MAXIMUM LOAD AND GROW CRACK (CHANGE SPRING STIFFNESS)

    IF(WORD.EQ.GROW) CALL BREAK
    IF(KLU.EQ.1) GOTO 15

35  CONTINUE
    DELP = -DELP
    READ(5,'(E10.3,1X,A4)') S,WORD
    IF(WORD.EQ.HALT) RETURN
    GOTO 10

    END SUBROUTINE PLAS

!===============================================================================

    SUBROUTINE CONTACT

    ! CHANGES SPRING STIFFNESS IF CRACK CLOSES OR OPENS

    USE C_MAIN
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE

    REAL(REAL64) :: ALP,CX,SMAX,SN
    INTEGER      :: I,IERR,IFAC,L,MPTX,NB,NC,NV,NX

    EXTERNAL SYMBAN

    NB = NBC
    NX = NOD(NB)
    DO I = 1,NX
      KNEW = 0
      L = LSAME(NB,I)
      CX = XR(L)-CRACK
      IF(CX.LT.-0.001) THEN
        NV = 2*L
        MPTX = MPTAB(NV)
        IF(VOLD(L).LE.0.0) MPTAB(NV) = 1
        IF(VOLD(L).GT.0.0) MPTAB(NV) = 0
        IF(MPTX.NE.MPTAB(NV)) KNEW = 1
        IF(KNEW.EQ.0) CYCLE
        SMAX = PMAX/WIDTH
        SN = SNET-((SNET-SMAX)/(VOLD(L)-V2(L)))*VOLD(L)
        Z(NV) = 1.0
        NC = NV
        ALP = SK2
        IF(MPTAB(NV).EQ.0) ALP = -SK2
        IFAC = 3
        CALL SYMBAN(MAXA,NDOF,MB,MSUM,AA,1,BB,IFAC,T1,IERR,ALP,Z,T2,T3,NC)
        IF(MPTAB(NV).EQ.0) WRITE(6,'("NODE ",I5," OPENED AT ",F9.3," KSI")') L,SN
        IF(MPTAB(NV).EQ.1) WRITE(6,'("NODE ",I5," CLOSED AT ",F9.3," KSI")') L,SN
      ENDIF
    ENDDO

    END SUBROUTINE CONTACT

!===============================================================================

    SUBROUTINE BREAK

    ! GROWS CRACK AT MAXIMUM LOAD (CHANGE SPRING STIFFNESS)

    USE C_MAIN
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE

    REAL(REAL64) :: ALP,CMIN,CX
    INTEGER      :: I,IERR,IFAC,L,N,NB,NC,NV,NX,NRHS

    EXTERNAL SYMBAN

    NB = NBC
    NX = NOD(NB)
    DO I = 1,NX
      L = LSAME(NB,I)
      CX = ABS(XR(L)-CRACK)
      IF(CX.LT.0.001) EXIT
    ENDDO
    NV = 2*L
    MPTAB(NV) = 0
    KLU = 1
    ALP = -SK2
    Z(NV) = 1.0
    NC = NV
    IFAC = 3
    NRHS = 1
    CALL SYMBAN(MAXA,NDOF,MB,MSUM,AA,NRHS,BB,IFAC,T1,IERR,ALP,Z,T2,T3,NC)
    WRITE(6,'("NODE ",I0," BROKE AT ",F9.3," KSI")') L,SNET
    CMIN = 1.0E+10
    DO I = 1,NX
      N = LSAME(NB,I)
      CX = XR(L)-XR(N)
      IF(CX.GT.0.001 .AND. CX.LT.CMIN) CMIN = CX
    ENDDO
    CRACK = CRACK+CMIN

    END SUBROUTINE BREAK

!===============================================================================

    SUBROUTINE PLOUT(NODE1,NODE2,NELE1,NELE2)

    ! PRINTS OUT APPLIED STRESS, CRACK LENGTH, DISPLACEMENTS  AT
    ! NODES NODE1 TO NODE2, AND STRESSES AND STRAINS FOR ELEMENTS
    ! NELE1 TO NELE2

    USE C_MAIN
    USE ISO_FORTRAN_ENV
    IMPLICIT NONE

    ! DUMMY ARGUMENT DECLARATIONS

    INTEGER,INTENT(IN) :: NODE1
    INTEGER,INTENT(IN) :: NODE2
    INTEGER,INTENT(IN) :: NELE1
    INTEGER,INTENT(IN) :: NELE2

    ! LOCAL VARIABLE DECLARATIONS

    REAL(REAL64) :: CRACK1
    INTEGER      :: I,L,M,N,NB,NX

    WRITE(6,'(/,"APPLIED GROSS STRESS = ",F9.3," KSI",5X,"CRACK = ",F9.3,5X,"WIDTH = ",F9.3)') &
      SNET,CRACK,WIDTH
    WRITE(6,'(/,A8,A12,A13)') 'NODE','X','COD'
    CRACK1 = CRACK+0.001
    NB = NBC
    NX = NOD(NB)
    DO I = 1,NX
      L = LSAME(NB,I)
      IF(YR(L).GT.0.001)CYCLE
      IF(XR(L).GT.CRACK1)CYCLE
      WRITE(6,'(I8,F12.3,1PE13.4)') L,XR(L),VOLD(L)
    ENDDO
    WRITE(6,'(/,"DISPLACEMENTS",/)')
    WRITE(6,'(A8,2A15)') 'NODE','U','V'
    DO N = NODE1,NODE2
      WRITE(6,'(I8,1P,2E15.6)') N,UOLD(N),VOLD(N)
    ENDDO
    WRITE(6,'(/,"STRESSES AND STRAINS",/)')
    WRITE(6,'("* DENOTES PLASTIC ELEMENTS",/)')
    WRITE(6,'(A7,A1,3A12,3A12)') 'ELEMENT',' ','SIGX','SIGY','TAU','STRX','STRY','STRXY'
    DO M = NELE1,NELE2
      IF(MPLAS(M).NE.0) THEN
        WRITE(6,'(I7,A1,3F12.3,1P,3E12.3)') &
          MPLAS(M),'*',SIGX(M),SIGY(M),TAU(M),STRX(M),STRY(M),STRXY(M)
      ELSE
        WRITE(6,'(I7,A1,3F12.3,1P,3E12.3)') &
          M,       ' ',SIGX(M),SIGY(M),TAU(M),STRX(M),STRY(M),STRXY(M)
      ENDIF
    ENDDO
    WRITE(6,'()')

    END SUBROUTINE PLOUT

!================================================================================

    SUBROUTINE SYMBAN(MAXN,N,M,MSUM,A,NRHS,B,IFAC,P,IERR,ALP,Z,W,D,NC)

    ! SOLVE MATRIX EQUATION AX=B WHERE A IS SYMMETRIC POSITIVE DEFINITE
    ! AND B IS A MATRIX OF CONSTANT VECTORS.
    !
    ! SYMMETRY AND BANDWIDTH OF MATRIX A IS ACCOUNTED FOR BY STORING A(I)
    ! IN LOWER TRIANGULAR MATRIX (INCLUDING DIAGONAL) BY ROWS
    !
    ! MAXN - MAXIMUM DIMENSION OF MATRIX A
    !
    ! N - NUMBER OF ROWS OF MATRIX A
    !
    ! M(I) - BAND WIDTHS (LARGEST NUMBER OF NON-ZERO COLUMNS TO THE LEFT
    !        AND INCLUDING THE DIAGONAL OF ROW I IN MATRIX A
    !        M(I) MUST BE GREATER THAN OR EQUAL TO 2 WITH M(1) = 1
    !
    ! MSUM(I) - AN ARRAY COMPUTED IN SYMBAN
    !
    ! NRHS - NUMBER OF RIGHT-HAND SIDES OF COLUMN VECTOR B
    !
    ! IFAC - INPUT INTEGER SPECIFYING WHETHER OR NOT A CHOLESKY
    !        DECOMPOSITION OF MATRIX A IS TO BE COMPUTED
    !        WHERE A = L*D*L**T
    !
    !      = 0   CHOLESKY DECOMPOSITION IS COMPUTED. IFAC SET TO 1.
    !
    !      = 1   THE CHOLESKY DECOMPOSED FORM OF MATRIX A IS INPUT.
    !            SOLUTION IS RETURNED IN B.
    !
    !      = 2   THE CHOLESKY DECOMPOSITION OF MATRIX A IS MODIFIED.
    !            MODIFICATION IS OF THE FORM A = A + ALP*Z*Z**T.
    !            NC IS THE FIRST NONZERO ROW IN COLUMN VECTOR Z.
    !            IFAC IS SET TO 1 AND Z SET TO 0 UPON RETURN.
    !            SOLUTION IS RETURNED IN B.
    !
    !      = 3   ONLY THE CHOLESKY DECOMPOSITION OF MATRIX A IS MODIFIED.
    !            IFAC IS SET TO 1 AND Z SET TO 0 UPON RETURN.

    USE ISO_FORTRAN_ENV
    IMPLICIT NONE

    ! DUMMY ARGUMENT DECLARATIONS

    INTEGER,INTENT(IN)         :: MAXN,N,NRHS
    INTEGER,INTENT(IN)         :: M(N)
    INTEGER,INTENT(INOUT)      :: MSUM(N)
    REAL(REAL64),INTENT(INOUT) :: A(MAXN)
    REAL(REAL64),INTENT(INOUT) :: B(N,NRHS)
    INTEGER,INTENT(INOUT)      :: IFAC
    REAL(REAL64),INTENT(INOUT) :: P(N)
    INTEGER,INTENT(OUT)        :: IERR
    REAL(REAL64),INTENT(INOUT) :: ALP
    REAL(REAL64),INTENT(INOUT) :: Z(N)
    REAL(REAL64),INTENT(INOUT) :: W(N)
    REAL(REAL64),INTENT(INOUT) :: D(N)
    INTEGER,INTENT(IN)         :: NC

    ! LOCAL VARIABLE DECLARATIONS

    REAL(REAL64) :: BETA,X,Y
    INTEGER      :: I,II,IM1,J,JM1,JO,K,L,L1,L2,LO,MN,MO,NM1,NN

    NM1 = N-1

    IF(NM1.GT.1) THEN
      IF(IFAC.GE.2) THEN
        !**************   IFAC = 2 OR 3   *************
        IF(NC.NE.N) THEN
          NN = N-1
          DO I = 1,N
            W(I) = Z(I)
          ENDDO
          DO J = NC,NN
            Z(J) = W(J)
            D(J) = D(J)+ALP*Z(J)**2
            BETA = ALP*Z(J)/D(J)
            ALP = ALP/(D(J)*P(J))
            P(J) = 1.0/D(J)
            DO I = J+1,N
              LO = I-M(I)
              MO = LO+1
              IF(MO.LE.J) THEN
                L = MSUM(I)+J-LO
                W(I) = W(I)-Z(J)*A(L)
                A(L) = A(L)+BETA*W(I)
              ENDIF
            ENDDO
          ENDDO
        ELSE
          D(NC) = D(NC)+ALP*Z(NC)**2
          P(NC) = 1.0/D(NC)
        ENDIF
        DO I = 1,N
          Z(I) = 0.0
        ENDDO
        IF(IFAC.EQ.3) THEN
          IFAC = 1
          RETURN
        ENDIF
        IFAC = 1
      ENDIF
      IF(IFAC.NE.1) THEN
        ! IFAC = 0
        MSUM(1) = 0
        MSUM(2) = 1
        DO I = 3,N
          II = I-1
          MSUM(I) = MSUM(II)+M(II)
        ENDDO
        IFAC = 1
        DO I = 1,N
          IM1 = I-1
          MN = I-M(I)+1
          LO = MSUM(I)-I+M(I)
          DO J = MN,I
            L = LO+J
            X = A(L)
            IF(I.LE.J) THEN
              IF(IM1.GT.0 .AND. MN.LE.IM1) THEN
                DO K = MN,IM1
                  L1 = LO+K
                  Y = A(L1)
                  A(L1) = Y*P(K)
                  X = X-Y*A(L1)
                ENDDO
              ENDIF
              IF(X.LE.0) THEN
                IERR = 1
                RETURN
              ENDIF
              D(I) = X
              P(I) = 1.0/X
            ENDIF
            JM1 = J-1
            IF(JM1.GT.0 .AND. MN.LE.JM1) THEN
              MO = J-M(J)+1
              JO = MSUM(J)-J+M(J)
              DO K = MN,JM1
                IF(K.GE.MO) THEN
                  L1 = LO+K
                  L2 = JO+K
                  X = X-A(L1)*A(L2)
                ENDIF
              ENDDO
            ENDIF
            A(L) = X
          ENDDO
        ENDDO
      ENDIF
      ! IFAC = 1
      DO I = 2,N
        IM1 = I-1
        MN = I-M(I)+1
        LO = MSUM(I)-I+M(I)
        DO J = 1,NRHS
          X = B(I,J)
          DO K = MN,IM1
            L = LO+K
            X = X-A(L)*B(K,J)
          ENDDO
          B(I,J) = X
        ENDDO
      ENDDO
      Y = P(N)
      DO J = 1,NRHS
        B(N,J) = B(N,J)*Y
      ENDDO
      DO WHILE(NM1.GT.0)
        I = NM1+1
        Y = P(NM1)
        DO J = 1,NRHS
          X = B(NM1,J)*Y
          DO K = I,N
            LO = K-M(K)
            MO = LO+1
            IF(MO.LE.NM1) THEN
              L = MSUM(K)+NM1-LO
              X = X-A(L)*B(K,J)
            ENDIF
          ENDDO
          B(NM1,J) = X
        ENDDO
        NM1 = NM1-1
      ENDDO
      IERR = 0
    ELSE
      I = 1
      IF(A(I).LE.0.0) THEN
        IERR = 1
        RETURN
      ENDIF
      P(I) = 1.0/A(I)
      DO J = 1,NRHS
        B(I,J) = B(I,J)*P(I)
      ENDDO
      IERR = 0
    ENDIF

    END SUBROUTINE SYMBAN

!================================================================================

    SUBROUTINE DIMERR(WHAT,NEED,LIMIT)

    ! STOPS WITH A MESSAGE WHEN A PROBLEM EXCEEDS ONE OF THE
    ! ARRAY DIMENSIONS SET IN PARAMETER.

    IMPLICIT NONE

    ! DUMMY ARGUMENT DECLARATIONS

    CHARACTER(*),INTENT(IN) :: WHAT
    INTEGER,INTENT(IN)      :: NEED
    INTEGER,INTENT(IN)      :: LIMIT

    WRITE(6,'()')
    WRITE(6,'("*** DIMENSION EXCEEDED: ",A," NEEDS ",I0,", LIMIT IS ",I0)') WHAT,NEED,LIMIT
    WRITE(6,'("*** INCREASE THE PARAMETER VALUE AND RECOMPILE")')
    WRITE(6,'("*** PROGRAM TERMINATED")')
    STOP

    END SUBROUTINE DIMERR
