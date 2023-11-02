!DIR$ FREEFORM

! MODULE FOR SHARED PARAMETERS
MODULE SHARED_DATA
	IMPLICIT NONE
	! CHARACTER(*), PARAMETER :: PATH='/Structure_HPC5/Rupsagar/06_Continuum_Damage_Mechanics/02_WIP/02_Peerlings/v4_v3_with_r_0_01'
	CHARACTER(*), PARAMETER :: PATH='D:/Programming/fortran/Continuum_Damage_Mechanics/02_Peerlings'
	
	CHARACTER(*), PARAMETER :: RUNTIME_FILE='Runtime_Parameters.txt'
	CHARACTER(*), PARAMETER :: OUTPUT_FILE='Analysis_Output.txt'
	
	CHARACTER(*), PARAMETER :: PATH_RUNTIME_FILE=PATH//'/'//RUNTIME_FILE
	CHARACTER(*), PARAMETER :: PATH_OUTPUT_FILE=PATH//'/'//OUTPUT_FILE
	
	INTEGER NUMEL,NINTP,NCYC_REAL,NCYC_LIFE,NCYC_JUMP,INCR_FLAG,EXIT_FLAG
	DOUBLE PRECISION CUTOFF,NU,ALPHA,BETA,C,KAPPA,ETA
	DOUBLE PRECISION, ALLOCATABLE :: DMG(:,:),DDOT(:,:)
	SAVE
END

! SUBROUTINE UEXTERNALDB
SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
	USE SHARED_DATA
	INCLUDE 'ABA_PARAM.INC'
	DIMENSION TIME(2)
	
	INTEGER IO_STAT,I,J
	
	! INITIALIZE MUTEX FOR THREAD SAFETY
	CALL MUTEXINIT(1)
	CALL MUTEXINIT(2)
	
	IF (LOP.EQ.0) THEN		! START OF ANALYSIS
		! INITIALIZE VARIABLES AND READ DAMAGE PARAMETERS OF MATERIAL
		IO_STAT=-1
		OPEN(10,FILE=PATH_RUNTIME_FILE,STATUS='UNKNOWN',IOSTAT=IO_STAT)
		IF (IO_STAT.EQ.0) THEN
			READ(10,*) NUMEL
			READ(10,*) NINTP
			READ(10,*) CUTOFF
			READ(10,*) ETA
			READ(10,*) NU
			READ(10,*) KAPPA
			READ(10,*) ALPHA
			READ(10,*) BETA
			READ(10,*) C
			WRITE(6,'(/,A,I8)') 'RUNTIME FILE OPENED. IOSTAT=',IO_STAT
			WRITE(6,'(2(A,I8),7(A,F8.6),/)') 'NUMEL=',NUMEL,', NINTP=',NINTP,', CUTOFF=',CUTOFF,', ETA=',ETA,', NU=',NU, ', KAPPA=',KAPPA,', ALPHA=',ALPHA,', BETA=',BETA,', C=',C
		ELSE
			WRITE(6,'(/,A,I8,/)') 'RUNTIME FILE CANNOT BE OPENED. IOSTAT=',IO_STAT
		ENDIF
		CLOSE(10)
		ALLOCATE(DMG(NUMEL,NINTP),DDOT(NUMEL,NINTP))
		DO I=1,NUMEL
			DO J=1,NINTP
				IF (DDOT(I,J).GT.DDOT_MAX) THEN
					DDOT(I,J)=0
					DMG(I,J)=0
				ENDIF
			ENDDO
		ENDDO
	ELSEIF (LOP.EQ.1) THEN		! START OF CURRENT ANALYSIS INCREMENT
		IF ((EXIT_FLAG.EQ.1).AND.(INCR_FLAG.EQ.0)) THEN
			CALL XIT
		ENDIF
	ELSEIF (LOP.EQ.2) THEN		! END OF CURRENT ANALYSIS INCREMENT
		CALL CHECK_JUMP(TIME(1),TIME(2),DTIME)
	ENDIF

	RETURN
END

! SUBROUTINE USDFLD
SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,&
	TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,&
	KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
	
	USE SHARED_DATA
	INCLUDE 'ABA_PARAM.INC'
	
	CHARACTER*80 CMNAME,ORNAME
	CHARACTER*3  FLGRAY(15)
	DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),T(3,3),TIME(2)
	DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
	
	! HEADER FILES FOR MUTEX
	! INCLUDE 'SMAAspUserArrays.hdr'
	! INCLUDE 'SMAAspUserSubroutines.hdr'
    ! INCLUDE 'SMAAspUserUtilities.hdr'
	
	! USER DEFINED VARIABLES
	DOUBLE PRECISION DMG_PREV,DMG_CUR,EPS_A1,EPS_A2,YOUNGS_MOD,SIGMA_12,&
		EPS_11,EPS_22,EPS_33,EPS_12,EPS_13,EPS_23,I1,J2,EPS_EQ,G1,G2,TEMP1,TEMP2

	! SDV RECALL 
	DMG_PREV=STATEV(1)
	YOUNGS_MOD=STATEV(5)
	EPS_A1=STATEV(6)
	EPS_A2=STATEV(7)

	! SHEAR STRESS COMPONENT
	CALL GETVRM('S',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,LACCFLA)
	IF(JRCD.EQ.0) THEN
		SIGMA_12=ARRAY(4)
	ENDIF
	
	! STRAIN TENSOR
	CALL GETVRM('E',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,LACCFLA)
	IF(JRCD.EQ.0) THEN
	  EPS_11=ARRAY(1)
	  EPS_22=ARRAY(2)
	  EPS_33=ARRAY(3)
	  EPS_12=ARRAY(4)/2.0
	  EPS_13=ARRAY(5)/2.0
	  EPS_23=ARRAY(6)/2.0
	ENDIF
	
	! YOUNG'S MODULUS
	IF (((TIME(1)-NCYC_REAL).GT.0).AND.((TIME(1)-NCYC_REAL).LE.0.25)) THEN
		YOUNGS_MOD=(1+NU)*SIGMA_12/EPS_12
	ENDIF
	
	! EQUIVALENT STRAIN AMPLITUDE
	I1=EPS_11+EPS_22+EPS_33
	J2=(I1**2)/6.0-0.5*(EPS_11**2+EPS_22**2+EPS_33**2+2*(EPS_12**2+EPS_13**2+EPS_23**2))
	EPS_EQ=((-3*J2)**0.5)/(1+NU)
	IF (((TIME(1)-NCYC_REAL).GT.0).AND.((TIME(1)-NCYC_REAL).LT.0.5)) THEN
		IF (EPS_EQ.GT.EPS_A1) THEN
			EPS_A1=EPS_EQ
		ENDIF	
	ELSEIF (((TIME(1)-NCYC_REAL).GT.0.5).AND.((TIME(1)-NCYC_REAL).LT.1.0)) THEN
		IF (EPS_EQ.GT.EPS_A2) THEN
			EPS_A2=EPS_EQ
		ENDIF
	ENDIF
	
	! DAMAGE RATE OF PRESENT CYCLE PRIOR TO DAMAGE UPDATE
	IF ((TIME(1)+DTIME-NCYC_REAL).EQ.1) THEN
		IF (EPS_A1.GT.KAPPA) THEN
			TEMP1=(EPS_A1-KAPPA)**(1+BETA)
		ELSE
			TEMP1=0
		ENDIF
		IF (EPS_A2.GT.KAPPA) THEN
			TEMP2=(EPS_A2-KAPPA)**(1+BETA)
		ELSE
			TEMP2=0
		ENDIF
		G1=C*EXP(ALPHA*DMG(NOEL,NPT))*(TEMP1)/(1+BETA)
		G2=C*EXP(ALPHA*DMG(NOEL,NPT))*(TEMP2)/(1+BETA)
		! THREAD SAFETY
		CALL MUTEXLOCK(1)
			DDOT(NOEL,NPT)=G1+G2
		CALL MUTEXUNLOCK(1)
	ENDIF
	
	! DAMAGE UPDATE 
	IF (INCR_FLAG.EQ.1) THEN
		DMG_CUR=DMG_PREV+DDOT(NOEL,NPT)*NCYC_JUMP
		IF (DMG_CUR.GT.CUTOFF) THEN
			DMG_CUR=DMG_PREV
		ENDIF
		! THREAD SAFETY
		CALL MUTEXLOCK(2)
			DMG(NOEL,NPT)=DMG_CUR
		CALL MUTEXUNLOCK(2)
	ENDIF

	! FIELD VARIABLE TO UPDATE STIFFNESS
	FIELD(1)=DMG(NOEL,NPT)

	! SDV UPDATE 
	STATEV(1)=DMG(NOEL,NPT)
	STATEV(3)=NCYC_LIFE
	STATEV(4)=NCYC_REAL
	STATEV(5)=YOUNGS_MOD
	STATEV(6)=EPS_A1
	STATEV(7)=EPS_A2
	
	RETURN
END

! SUBROUTINE FOR CYCLE JUMP
SUBROUTINE CHECK_JUMP(TIME1,TIME2,DTIME)
	USE SHARED_DATA
	DOUBLE PRECISION, INTENT(IN) :: TIME1,TIME2,DTIME
	DOUBLE PRECISION DMG_MAX,DDOT_MAX
	INTEGER IO_STAT,I,J,ITEMP,JTEMP
	LOGICAL EXIST_FLAG
	
	IF ((TIME1-NCYC_REAL).EQ.1) THEN
		! ADAPTIVE JUMP
		DO I=1,NUMEL
			DO J=1,NINTP
				IF (DDOT(I,J).GT.DDOT_MAX) THEN
					DDOT_MAX=DDOT(I,J)
					DMG_MAX=DMG(I,J)
					ITEMP=I
					JTEMP=J
				ENDIF
			ENDDO
		ENDDO
		! NCYC_JUMP=ETA/(DDOT_MAX*ALPHA)
		NCYC_JUMP=ETA/(DDOT_MAX)
		WRITE(6,'(/,2(A,F8.4),2(A,F8.6),4(A,I8))') 'AT ANALYSIS TIME(1)=',TIME1,', ETA=',ETA,', DDOT_MAX=',DDOT_MAX,', DMG_MAX=',DMG_MAX,', NOEL=',ITEMP,', NPT=',JTEMP,', NCYC_LIFE=',NCYC_LIFE,', NCYC_JUMP=',NCYC_JUMP
		
		! JUMP SCALED TO CUTOFF
		IF ((DMG_MAX+DDOT_MAX*NCYC_JUMP).GT.CUTOFF) THEN
			WRITE(6,'(A,F8.4,A,I8,A,F8.6)') 'AT ANALYSIS TIME(1)=',TIME1,', NCYC_JUMP=',NCYC_JUMP,', DMG_NEW=',(DMG_MAX+DDOT_MAX*NCYC_JUMP)
			NCYC_JUMP=(CUTOFF-DMG_MAX)/DDOT_MAX
			WRITE(6,'(A,F8.4,A,I8,A,F8.6)') 'AT ANALYSIS TIME(1)=',TIME1,', NCYC_JUMP=',NCYC_JUMP,', DMG_NEW=',(DMG_MAX+DDOT_MAX*NCYC_JUMP)
			EXIT_FLAG=1
			WRITE(6,'(2(A,F8.4),2(A,F8.6),4(A,I8))') 'AT ANALYSIS TIME(1)=',TIME1,', ETA=',ETA,', DDOT_MAX=',DDOT_MAX,', DMG_MAX=',DMG_MAX,', NOEL=',ITEMP,', NPT=',JTEMP,', NCYC_LIFE=',NCYC_LIFE,', NCYC_JUMP=',NCYC_JUMP
		ENDIF
		
		IF (NCYC_JUMP.GT.0) THEN
			NCYC_LIFE=NCYC_LIFE+NCYC_JUMP
			NCYC_REAL=NCYC_REAL+1
			INCR_FLAG=1
		ENDIF
	ELSE
		! PRINT OUTPUT FOR LATEST JUMP
		IF (INCR_FLAG.EQ.1) THEN
			IO_STAT=-1
			INQUIRE(FILE=PATH_OUTPUT_FILE, EXIST=EXIST_FLAG)
			IF (EXIST_FLAG.AND.(NCYC_REAL.EQ.1)) THEN
				OPEN(30, FILE=PATH_OUTPUT_FILE, STATUS='OLD', ACTION='WRITE',IOSTAT=IO_STAT)
			ELSEIF (EXIST_FLAG.AND.(NCYC_REAL.NE.1)) THEN
				OPEN(30, FILE=PATH_OUTPUT_FILE, STATUS='OLD', POSITION='APPEND', ACTION='WRITE',IOSTAT=IO_STAT)
			ELSE
				OPEN(30, FILE=PATH_OUTPUT_FILE, STATUS='NEW', ACTION='WRITE',IOSTAT=IO_STAT)
			END IF
			
			IF (IO_STAT.EQ.0) THEN
				WRITE(30,'(/,/,A,F8.4,2(A,I8),/,/)') 'ANALYSIS OUTPUT AT TIME(1)=',TIME1,', NCYC_REAL=',NCYC_REAL,', NCYC_LIFE=',NCYC_LIFE
				DO I=1,NUMEL
					DO J=1,NINTP
						WRITE(30,'(2(A,I8),2(A,F8.6))') 'NOEL=',I,', NPT',J,', DMG=',DMG(I,J),', DDOT=',DDOT(I,J)
					ENDDO
				ENDDO
				
				WRITE(6,'(A,F8.4,A,I8)') 'ANALYSIS OUTPUT WRITTEN TO FILE. TIME(1)=',TIME1,', IOSTAT=',IO_STAT
			ELSE
				WRITE(6,'(A,F8.4,A,I8)') 'ANALYSIS OUTPUT CANNOT BE WRITTEN TO FILE. TIME(1)=',TIME1,', IOSTAT=',IO_STAT
			ENDIF
			
			CLOSE(30)
			
		ENDIF
		
		INCR_FLAG=0
	ENDIF

	RETURN
END
