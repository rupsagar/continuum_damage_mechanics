      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
	  
	  INCLUDE 'ABA_PARAM.INC'
      DIMENSION TIME(2)
	  
	  call MutexInit( 1 )
	  
	  RETURN
	  END
	    
	  
	  SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
c
      INCLUDE 'ABA_PARAM.INC'
	  INCLUDE 'SMAAspUserArrays.hdr'
      INCLUDE 'SMAAspUserSubroutines.hdr'
      INCLUDE 'SMAAspUserUtilities.hdr'
c
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
c     
c	  ! USER DEFINITION BEGINS
c	  ! USER DEFINED VARIABLES
	  DOUBLE PRECISION NU,C,BETA,ALPHA,KAPPA,E11,E22,E33,E12,E13,E23
	  DOUBLE PRECISION I1,J2,DMG_PREV,G1,G2,DMG_CUR,EPS_EQ
	  DOUBLE PRECISION G11,G12,G21,G22,CYC_JUMP1, X, MAX_G 
      INTEGER CYC_REAL,CYC_PSEUDO,CONST_JUMP,DMG_INCR_FLAG,CYC_JUMP,IO
	  INTEGER ISOPEN
	  
c	  
c	  ! USER DEFINED PARAMETERS
	  NU=0.3
      ALPHA=7.0
	  C= 1.56/ALPHA
	  BETA= 0.7
	  KAPPA = 174.0/70000.0
	  G11 = 0.0
	  G12 = 0.0
	  G21 = 0.0
	  G22 = 0.0
	  G1 = 0.0
	  G2 = 0.0
	  DMG_PREV = STATEV(1)
	  EA1 = STATEV(2)
	  EA2 = STATEV(3)
	  CYC_REAL = STATEV(4)
	  CYC_PSEUDO = STATEV(5)
	  !CYC_JUMP = 3000
	  IO = 100
	  ISOPEN = 1
c	  
! c	  ! STRAIN VALUES
	  CALL GETVRM('E',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1 MATLAYO,LACCFLA)
	  IF (JRCD .EQ. 0) THEN
		 E11=ARRAY(1)
		 E22=ARRAY(2)
		 E33=ARRAY(3)
		 E12=ARRAY(4)
		 E13=ARRAY(5)
		 E23=ARRAY(6)
	  END IF
	  
	  I1=E11+E22+E33
	  J2=I1**2.0/6-0.5*(E11**2+E22**2+E33**2+2*(E12**2+E13**2+E23**2))
	  EPS_EQ=1.0/(1.0+NU)*(-3.0*J2)**0.5
	  
 	  IF ((TIME(2)+DTIME-CYC_REAL) .GE. 1) THEN 
		 DMG_INCR_FLAG=1
	  ELSE
		 DMG_INCR_FLAG=0
	  ENDIF
	  
	  IF (DMG_INCR_FLAG .NE. 1) THEN
	     IF (EPS_EQ .GT. EA1) THEN
		    EA1 = EPS_EQ
	     ELSE
		    EA1 = EA1
	     ENDIF
	     IF (EPS_EQ .GT. EA2) THEN
		    EA2 = EPS_EQ
	     ELSE
		    EA2 = EA2
	     ENDIF
	  END IF
	  
	  IF (EA1 .GT. KAPPA) THEN
	     !X = EA1**(BETA+1) - KAPPA**(BETA+1)
		 X = (EA1 - KAPPA)**(BETA+1)
	     G11 = C/(BETA+1)*EXP(ALPHA*DMG_PREV)*(X)
	  ELSE
	     G11 =0.0
	  END IF
	  IF (EA2 .GT. KAPPA) THEN
	     !X = EA2**(BETA+1) - KAPPA**(BETA+1)
		 X = (EA2 - KAPPA)**(BETA+1)
	     G12 = C/(BETA+1)*EXP(ALPHA*DMG_PREV)*(X)
	  ELSE
	     G12 = 0.0
	  END IF
	  G1 = 2*G11
	  
	  IF ((TIME(2)-CYC_REAL) .GT. 0.25) THEN 
	     call MutexLock(1)
	     OPEN(15,file='/Structure_HPC5/Vivek_Agarwal/01_Abaqus/08_Coupon/txt1.dat',status='OLD',iostat=IO)
	     READ(15,*) MAX_G
         CLOSE(15)
	     IF (G1 .GT. MAX_G) THEN
		    MAX_G = G1
		    OPEN(15,file='//Structure_HPC5/Vivek_Agarwal/01_Abaqus/08_Coupon/txt1.dat',status='REPLACE')
	        WRITE(15,*) MAX_G
            CLOSE(15)
	     ELSE
	        MAX_G = MAX_G
         END IF
         call MutexUnlock(1)
      END IF
	  
	
	  IF (DMG_INCR_FLAG .EQ. 1) THEN	  
	     IF (DMG_PREV .LT. 1.0) THEN
			CYC_JUMP = 0.3/ALPHA/MAX_G
			IF (CYC_JUMP .GT. 10000) THEN
			   CYC_JUMP = 10000
			ELSE
			   CYC_JUMP = CYC_JUMP
			END IF
			IF (CYC_JUMP .LT. 1000) THEN
			   CYC_JUMP = 1000
			ELSE
			   CYC_JUMP = CYC_JUMP
			END IF
			CYC_REAL=CYC_REAL+1
			CYC_PSEUDO=CYC_PSEUDO+CYC_JUMP
		    DMG_P = DMG_PREV + G1*CYC_JUMP
			IF (EA1 .GT. KAPPA) THEN
			   !X = EA1**(BETA+1) - KAPPA**(BETA+1)
			   X = (EA1 - KAPPA)**(BETA+1)
	           G21 = C/(BETA+1)*EXP(ALPHA*DMG_P)*(X)
			ELSE
	           G21 =0
	        END IF
		    IF (EA2 .GT. KAPPA) THEN
			   !X = EA2**(BETA+1) - KAPPA**(BETA+1)
			   X = (EA2 - KAPPA)**(BETA+1)
	           G22 = C/(BETA+1)*EXP(ALPHA*DMG_P)*(X)
			ELSE
	           G22 = 0
	        END IF
		    G2 = 2*G21
		    DMG_CUR = DMG_PREV + 1.0/2.0*(G1+G2)*CYC_JUMP
	     ELSE
	        DMG_CUR = DMG_PREV
	     END IF
		 IF (DMG_CUR .GE. 1) THEN
		    DMG_CUR = 1.0
		 END IF
		 EA1 = 0
		 EA2 = 0
	  ELSE
	     DMG_CUR = DMG_PREV
	  ENDIF
	  
	  IF (DMG_INCR_FLAG .EQ. 1) THEN
	     IF (EPS_EQ .GT. EA1) THEN
		    EA1 = EPS_EQ
	     ELSE
		    EA1 = EA1
	     ENDIF
	     IF (EPS_EQ .GT. EA2) THEN
		    EA2 = EPS_EQ
	     ELSE
		    EA2 = EA2
	     ENDIF
      END IF	  

	  FIELD(1) = DMG_CUR
	  STATEV(1) = DMG_CUR
	  STATEV(2) = EA1
	  STATEV(3) = EA2
	  STATEV(4) = CYC_REAL
	  STATEV(5) = CYC_PSEUDO
      RETURN
      END	  
