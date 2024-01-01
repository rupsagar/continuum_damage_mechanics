!DIR$ FREEFORM

! MODULE FOR SHARED PARAMETERS
MODULE SHARED_DATA
	IMPLICIT NONE
	
	CHARACTER(*), PARAMETER :: RUNTIME_FILE='Data_Input_Runtime.txt'
	CHARACTER(*), PARAMETER :: OUTPUT_FILE='Data_Output_Analysis.txt'
	CHARACTER(*), PARAMETER :: PROC_OUTPUT_FILE='Data_Output_Rank'
	CHARACTER(LEN=256) PATH_SRC,PATH_RUNTIME_FILE,PATH_OUTPUT_FILE
	CHARACTER(LEN=50) INT_OPTIONS(2),INT_METHOD
	
	INTEGER, PARAMETER :: NUM_LINE_OUTPUT=9
	INTEGER LEN_PATH,NUM_ELEM,NUM_INTP,NCYC_REAL,NCYC_LIFE,NCYC_JUMP,CYC_END_FLAG,OUTPUT_FLAG,EXIT_FLAG
	
	INTEGER NUM_PROC,PROC_RANK
	
	DOUBLE PRECISION CUTOFF,DEL_DMG,NU,KAPPA,ALPHA,BETA,C
	DOUBLE PRECISION, ALLOCATABLE :: DMG(:,:),DDOT(:,:)
	
	DOUBLE PRECISION DMG_MAX_TEMP,DDOT_MAX_TEMP,DMG_TEMP
	INTEGER NOEL_TEMP,NPT_TEMP,NCYC_JUMP_TEMP
	
	SAVE
END

! SUBROUTINE UEXTERNALDB
SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
	USE SHARED_DATA
	INCLUDE 'ABA_PARAM.INC'
	DIMENSION TIME(2)

	! INITIALIZE MUTEX FOR THREAD SAFETY
	CALL MUTEXINIT(1)
	CALL MUTEXINIT(2)
	
	! START OF ANALYSIS
	IF (LOP.EQ.0) THEN
		CALL GETOUTDIR(PATH_SRC,LEN_PATH)
		CALL GETNUMCPUS(NUM_PROC)
		CALL GETRANK(PROC_RANK)
		CALL DATA_RUNTIME
		CALL INIT_VAR
	! START OF CURRENT ANALYSIS INCREMENT
	ELSEIF (LOP.EQ.1) THEN
		IF (CYC_END_FLAG.EQ.1) THEN
			CALL CALC_GLOBAL_CYC_JUMP(TIME(1))
		ELSE
			IF (EXIT_FLAG.EQ.1) THEN
				CALL XIT
			ENDIF
		ENDIF
	! END OF CURRENT ANALYSIS INCREMENT
	ELSEIF (LOP.EQ.2) THEN
		CALL COND_CYC_END(TIME(1))
		IF (CYC_END_FLAG.EQ.1) THEN
			CALL CALC_PROC_CYC_JUMP(TIME(1))
		ENDIF
	ENDIF
	
	RETURN
END

! READ PARAMETERS
SUBROUTINE DATA_RUNTIME()
	USE SHARED_DATA
	IMPLICIT NONE
	INTEGER IO_STAT
	INTEGER, PARAMETER :: FILE_UNIT=10
	CHARACTER(LEN=50) PARAM_NAME
	PATH_RUNTIME_FILE=TRIM(PATH_SRC)//'/'//RUNTIME_FILE
	OPEN(FILE_UNIT,FILE=TRIM(PATH_RUNTIME_FILE),IOSTAT=IO_STAT,STATUS='OLD')
	IF (IO_STAT.EQ.0) THEN
		READ(FILE_UNIT,*) PARAM_NAME,NUM_ELEM
		READ(FILE_UNIT,*) PARAM_NAME,NUM_INTP
		READ(FILE_UNIT,*) PARAM_NAME,OUTPUT_FLAG
		READ(FILE_UNIT,*) PARAM_NAME,INT_METHOD
		READ(FILE_UNIT,*) PARAM_NAME,CUTOFF
		READ(FILE_UNIT,*) PARAM_NAME,DEL_DMG
		READ(FILE_UNIT,*) PARAM_NAME,NU
		READ(FILE_UNIT,*) PARAM_NAME,KAPPA
		READ(FILE_UNIT,*) PARAM_NAME,ALPHA
		READ(FILE_UNIT,*) PARAM_NAME,BETA
		READ(FILE_UNIT,*) PARAM_NAME,C
		WRITE(*,'(/,A,I8,2(A,I2),/,3(A,I8),1(A,A),/,6(A,F8.6),A,D12.4,/)') 'RUNTIME FILE OPENED. IOSTAT=',IO_STAT,', NUM_PROC=',NUM_PROC,', PROC_RANK=',PROC_RANK,&
			'NUM_ELEM=',NUM_ELEM,', NUM_INTP=',NUM_INTP,', OUTPUT_FLAG=',OUTPUT_FLAG,', INT_METHOD=',INT_METHOD,&
			'CUTOFF=',CUTOFF,', DEL_DMG=',DEL_DMG,', NU=',NU, ', KAPPA=',KAPPA,', ALPHA=',ALPHA,', BETA=',BETA,', C=',C
	ELSE
		WRITE(*,'(/,A,I8,2(A,I2),/)') 'RUNTIME FILE CANNOT BE OPENED. IOSTAT=',IO_STAT,', NUM_PROC=',NUM_PROC,', PROC_RANK=',PROC_RANK
	ENDIF
	CLOSE(FILE_UNIT)
	
	RETURN
END

! INITIALIZE VARIABLES
SUBROUTINE INIT_VAR()
	USE SHARED_DATA
	IMPLICIT NONE
	INTEGER I,J
	INT_OPTIONS=['EULER_EXPLICIT','EULER_EXPLICIT_PREDICTOR']
	NCYC_REAL=0
	NCYC_LIFE=0
	NCYC_JUMP=0
	CYC_END_FLAG=0
	EXIT_FLAG=0
	ALLOCATE(DMG(NUM_ELEM,NUM_INTP),DDOT(NUM_ELEM,NUM_INTP))
	DO I=1,NUM_ELEM
		DO J=1,NUM_INTP
			DMG(I,J)=0.0
			DDOT(I,J)=0.0
		ENDDO
	ENDDO
	RETURN
END

! WRITE DATA FOR DISTRIBUTED COMPUTING
SUBROUTINE MPI_WRITE(TIME1)
	USE SHARED_DATA
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: TIME1
	INTEGER, PARAMETER :: FILE_UNIT=20
	CHARACTER(LEN=256) PATH_PROC_OUTPUT_FILE
	INTEGER IO_STAT
	WRITE(PATH_PROC_OUTPUT_FILE,'(I2)') PROC_RANK
	PATH_PROC_OUTPUT_FILE=TRIM(PATH_SRC)//'/'//PROC_OUTPUT_FILE//'_'//TRIM(ADJUSTL(PATH_PROC_OUTPUT_FILE))//'.txt'
	IF (NCYC_REAL.EQ.0) THEN
		OPEN(FILE_UNIT,FILE=TRIM(PATH_PROC_OUTPUT_FILE),IOSTAT=IO_STAT,STATUS='UNKNOWN',ACTION='WRITE')
	ELSE
		OPEN(FILE_UNIT,FILE=TRIM(PATH_PROC_OUTPUT_FILE),IOSTAT=IO_STAT,STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
	ENDIF
	IF (IO_STAT.EQ.0) THEN
		WRITE(FILE_UNIT,'(A)') 'MPI_WRITE_UEXTERNALDB:'
		WRITE(FILE_UNIT,'(A,F8.4)') 'TIME=',TIME1
		WRITE(FILE_UNIT,'(A,I2)') 'NUM_PROC=',NUM_PROC
		WRITE(FILE_UNIT,'(A,I2)') 'PROC_RANK=',PROC_RANK
		WRITE(FILE_UNIT,'(A,F8.6)') 'DMG_MAX=',DMG_MAX_TEMP
		WRITE(FILE_UNIT,'(A,F8.6)') 'DDOT_MAX=',DDOT_MAX_TEMP
		WRITE(FILE_UNIT,'(A,I8)') 'NOEL_TEMP=',NOEL_TEMP
		WRITE(FILE_UNIT,'(A,I8)') 'NPT_TEMP=',NPT_TEMP
		WRITE(FILE_UNIT,'(A,I8)') 'NCYC_JUMP_TEMP=',NCYC_JUMP_TEMP
	ELSE
		WRITE(*,'(1(A,I8))') 'DATA OUTPUT RANK FILE NOT OPENED. IOSTAT=',IO_STAT
	ENDIF
	CLOSE(FILE_UNIT)
	RETURN
END

! READ DATA FOR DISTRIBUTED COMPUTING
SUBROUTINE MPI_READ()
	USE SHARED_DATA
	IMPLICIT NONE
	INTEGER, PARAMETER :: FILE_UNIT=30
	CHARACTER(LEN=256) PATH_PROC_OUTPUT_FILE
	CHARACTER(LEN=50) PARAM_NAME
	INTEGER I,J,IO_STAT,NCYC_JUMP_READ
	NCYC_JUMP_READ=-2
	DO I=0,(NUM_PROC-1)
		WRITE(PATH_PROC_OUTPUT_FILE,'(I2)') I
		PATH_PROC_OUTPUT_FILE=TRIM(PATH_SRC)//'/'//PROC_OUTPUT_FILE//'_'//TRIM(ADJUSTL(PATH_PROC_OUTPUT_FILE))//'.txt'
		OPEN(FILE_UNIT,FILE=TRIM(PATH_PROC_OUTPUT_FILE),IOSTAT=IO_STAT,STATUS='OLD')
		IF (IO_STAT.EQ.0) THEN
			DO J=1,((NCYC_REAL+1)*NUM_LINE_OUTPUT-1)
				READ(FILE_UNIT,*)
			ENDDO
			READ(FILE_UNIT,*) PARAM_NAME,NCYC_JUMP_READ
			IF (NCYC_JUMP_READ.GE.0) THEN
				IF (NCYC_JUMP_TEMP.LT.0) THEN
					NCYC_JUMP_TEMP=NCYC_JUMP_READ
				ELSE
					IF (NCYC_JUMP_READ.LT.NCYC_JUMP_TEMP) THEN
						NCYC_JUMP_TEMP=NCYC_JUMP_READ
					ENDIF
				ENDIF
			ENDIF
		ENDIF
		CLOSE(FILE_UNIT)
	ENDDO
	NCYC_JUMP=NCYC_JUMP_TEMP
	RETURN
END

! SUBROUTINE TO CHECK CYCLE END CONDITION
SUBROUTINE COND_CYC_END(TIME1)
	USE SHARED_DATA
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: TIME1
	IF ((TIME1-NCYC_REAL).EQ.1) THEN
		CYC_END_FLAG=1
	ELSE
		CYC_END_FLAG=0
	ENDIF
	RETURN
END

! SUBROUTINE TO FIND GLOBAL CYCLE JUMP MAGNITUDE
SUBROUTINE CALC_GLOBAL_CYC_JUMP(TIME1)
	USE SHARED_DATA
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: TIME1
	
	IF (NUM_PROC.GT.1) THEN
		CALL MPI_READ
	ELSE
		NCYC_JUMP=NCYC_JUMP_TEMP
	ENDIF
	
	IF (NCYC_JUMP.GT.0) THEN
		NCYC_LIFE=NCYC_LIFE+NCYC_JUMP
		NCYC_REAL=NCYC_REAL+1
		IF (OUTPUT_FLAG.EQ.1) THEN
			CALL WRITE_OUTPUT(TIME1)
		ENDIF
	ENDIF
	
	WRITE(*,'(A,F8.4,2(A,I2),2(A,F8.6),4(A,I8))') 'TIME=',TIME1,', NUM_PROC=',NUM_PROC, ', PROC_RANK=',PROC_RANK,', DMG_MAX_TEMP=',DMG_MAX_TEMP,&
		', DDOT_MAX_TEMP=',DDOT_MAX_TEMP,', NOEL_TEMP=',NOEL_TEMP,', NPT_TEMP=',NPT_TEMP,', NCYC_JUMP=     ',NCYC_JUMP,', NCYC_LIFE=',NCYC_LIFE
	
	RETURN
END

! SUBROUTINE TO FIND CYCLE JUMP MAGNITUDE OF PROCESS
SUBROUTINE CALC_PROC_CYC_JUMP(TIME1)
	USE SHARED_DATA
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: TIME1
	INTEGER I,J,NCYC_JUMP_TEMP_CO
	
	DMG_MAX_TEMP=0.0
	DDOT_MAX_TEMP=0.0
	DMG_TEMP=0.0
	NOEL_TEMP=-1
	NPT_TEMP=-1
	NCYC_JUMP_TEMP=-1
	DO I=1,NUM_ELEM
		DO J=1,NUM_INTP
			IF ((NOEL_TEMP.EQ.(-1)).AND.(NPT_TEMP.EQ.(-1))) THEN
				IF (DDOT(I,J).GE.DDOT_MAX_TEMP) THEN
					DMG_MAX_TEMP=DMG(I,J)
					DDOT_MAX_TEMP=DDOT(I,J)
					NOEL_TEMP=I
					NPT_TEMP=J
				ENDIF
			ELSE
				IF (DDOT(I,J).GT.DDOT_MAX_TEMP) THEN
					DMG_MAX_TEMP=DMG(I,J)
					DDOT_MAX_TEMP=DDOT(I,J)
					NOEL_TEMP=I
					NPT_TEMP=J
				ENDIF
			ENDIF
		ENDDO
	ENDDO
	IF (DDOT_MAX_TEMP.GT.0) THEN
		NCYC_JUMP_TEMP=CEILING(DEL_DMG/DDOT_MAX_TEMP)
	ENDIF
	
	! JUMP SCALED TO CUTOFF
	CALL CALC_DMG(DMG_TEMP,DMG_MAX_TEMP,DDOT_MAX_TEMP,NCYC_JUMP_TEMP)

	IF (DMG_TEMP.LE.CUTOFF) THEN
		WRITE(*,'(A,F8.4,2(A,I2),2(A,F8.6),3(A,I8))') 'TIME=',TIME1,', NUM_PROC=',NUM_PROC, ', PROC_RANK=',PROC_RANK,&
			', DMG_MAX_TEMP=',DMG_MAX_TEMP,', DDOT_MAX_TEMP=',DDOT_MAX_TEMP,', NOEL_TEMP=',NOEL_TEMP,', NPT_TEMP=',NPT_TEMP,', NCYC_JUMP_TEMP=',NCYC_JUMP_TEMP
	ELSE
		EXIT_FLAG=1
		NCYC_JUMP_TEMP_CO=NCYC_JUMP_TEMP
		NCYC_JUMP_TEMP=(CUTOFF-DMG_MAX_TEMP)/DDOT_MAX_TEMP
		WRITE(*,'(A,F8.4,2(A,I2),2(A,F8.6),3(A,I8),1(A,I8),2(A,F8.6))') 'TIME=',TIME1,', NUM_PROC=',NUM_PROC, ', PROC_RANK=',PROC_RANK,', DMG_MAX_TEMP=',DMG_MAX_TEMP,', DDOT_MAX_TEMP=',DDOT_MAX_TEMP,&
		', NOEL_TEMP=',NOEL_TEMP,', NPT_TEMP=',NPT_TEMP,', NCYC_JUMP_TEMP=',NCYC_JUMP_TEMP,', NCYC_JUMP_TEMP_CO=',NCYC_JUMP_TEMP_CO,', DMG_NEW=',(DMG_MAX_TEMP+DDOT_MAX_TEMP*NCYC_JUMP_TEMP),', DMG_NEW_CO=',(DMG_MAX_TEMP+DDOT_MAX_TEMP*NCYC_JUMP_TEMP_CO)
	ENDIF
	
	IF (NUM_PROC.GT.1) THEN
		CALL MPI_WRITE(TIME1)
	END IF
	
	RETURN
END

SUBROUTINE CALC_DMG(DMG_NEW,DMG_PREV,DDOT_CUR,DEL_NCYC_LIFE)
	USE SHARED_DATA
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(OUT) :: DMG_NEW
	DOUBLE PRECISION, INTENT(IN) :: DMG_PREV,DDOT_CUR
	INTEGER, INTENT(IN) :: DEL_NCYC_LIFE
	
	IF (TRIM(INT_METHOD).EQ.INT_OPTIONS(1)) THEN		! EULER FORWARD INTEGRATION (EXPLICIT)
		DMG_NEW=DMG_PREV+DDOT_CUR*DEL_NCYC_LIFE
	ELSEIF (TRIM(INT_METHOD).EQ.INT_OPTIONS(2))	THEN	! EULER FORWARD INTEGRATION WITH PREDICTOR (EXPLICIT)
		DMG_NEW=DMG_PREV+0.5*(DDOT_CUR+DDOT_CUR*EXP(DDOT_CUR*DEL_NCYC_LIFE))*DEL_NCYC_LIFE
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
	
	! USER DEFINED VARIABLES
	DOUBLE PRECISION DMG_PREV,DMG_CUR,YOUNGS_MOD,EPS_A1,EPS_A2,SIGMA_12,&
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
	IF (CYC_END_FLAG.EQ.1) THEN
		CALL CALC_DMG(DMG_CUR,DMG_PREV,DDOT(NOEL,NPT),NCYC_JUMP)
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

! PRINT OUTPUT FOR LATEST JUMP
SUBROUTINE WRITE_OUTPUT(TIME1)
	USE SHARED_DATA
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: TIME1
	INTEGER, PARAMETER :: FILE_UNIT=40
	INTEGER IO_STAT,I,J
	LOGICAL EXIST_FLAG
	PATH_OUTPUT_FILE=TRIM(PATH_SRC)//'/'//OUTPUT_FILE
	
	INQUIRE(FILE=TRIM(PATH_OUTPUT_FILE),EXIST=EXIST_FLAG)
	IF (EXIST_FLAG.AND.(NCYC_REAL.EQ.1)) THEN
		OPEN(FILE_UNIT,FILE=TRIM(PATH_OUTPUT_FILE),IOSTAT=IO_STAT,STATUS='OLD',ACTION='WRITE')
	ELSEIF (EXIST_FLAG.AND.(NCYC_REAL.NE.1)) THEN
		OPEN(FILE_UNIT,FILE=TRIM(PATH_OUTPUT_FILE),IOSTAT=IO_STAT,STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
	ELSE
		OPEN(FILE_UNIT,FILE=TRIM(PATH_OUTPUT_FILE),IOSTAT=IO_STAT,STATUS='NEW',ACTION='WRITE')
	END IF
	
	IF (IO_STAT.EQ.0) THEN
		WRITE(FILE_UNIT,'(/,/,A,F8.4,2(A,I8),/,/)') 'ANALYSIS OUTPUT AT TIME=',TIME1,', NCYC_REAL=',NCYC_REAL,', NCYC_LIFE=',NCYC_LIFE
		DO I=1,NUM_ELEM
			DO J=1,NUM_INTP
				WRITE(FILE_UNIT,'(2(A,I8),2(A,F8.6))') 'NOEL=',I,', NPT',J,', DMG=',DMG(I,J),', DDOT=',DDOT(I,J)
			ENDDO
		ENDDO
		WRITE(*,'(A,F8.4,A,I8)') 'ANALYSIS OUTPUT WRITTEN TO FILE. TIME=',TIME1,', IOSTAT=',IO_STAT
	ELSE
		WRITE(*,'(A,F8.4,A,I8)') 'ANALYSIS OUTPUT CANNOT BE WRITTEN TO FILE. TIME=',TIME1,', IOSTAT=',IO_STAT
	ENDIF
	
	CLOSE(FILE_UNIT)

	RETURN
END
