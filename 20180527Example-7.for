CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE USER
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This is the user part of the program. In this part, the users need
C  to specify all of the aspects related with the individual physical 
C  problems, including:
C  (1) the governing equations to be solved;
C  (2) the calculation domain, and the corresponding mesh generation;
C  (3) the physical properties of the materials;
C  (4) the boundary and/or initial conditions (if any);
C  (5) the source terms;
C  (6) the time interval if any, and the other parameters for the
C      iterations, i.e., the relaxation factors, the criteria for
C      ending the internal/external iterations.
***********************************************************************
      PARAMETER (NX=520,NY=56,NXY=520,NUM=16) !添加或修改齿数NUM
      CHARACTER (LEN=20) TITLE
      LOGICAL LSOLVE,LPRINT,LBLK,LSTOP,LEND
      COMMON F(NX,NY,13),FOLD(NX,NY,13),P(NX,NY),RHO(NX,NY),GAM(NX,NY),
     & CON(NX,NY),AIP(NX,NY),AIM(NX,NY),AJP(NX,NY),AJM(NX,NY),AP(NX,NY),
     & X(NXY),XU(NXY),XDIF(NXY),XCV(NXY),XCVS(NXY),
     & Y(NXY),YV(NXY),YDIF(NXY),YCV(NXY),YCVS(NXY),
     & YCVR(NXY),YCVRS(NXY),ARX(NXY),ARXJ(NXY),ARXJP(NXY),
     & R(NXY),RMN(NXY),SX(NXY),SXMN(NXY),XCVI(NXY),XCVIP(NXY)
      COMMON DU(NX,NY),DV(NX,NY),FV(NXY),FVP(NXY),
     & FX(NXY),FXM(NXY),FY(NXY),FYM(NXY),PT(NXY),QT(NXY),
     & BL(NUM,NY),RB(NUM,NY) !添加全局数组，存储密封齿边界位置
      COMMON/INDX/NF,NFMAX,NP,NRHO,NGAM,L1,L2,L3,M1,M2,M3,
     &  IST,JST,ITER,LAST,TITLE(13),RELAX(13),TIME,DT,XL,YL,
     &  IPREF,JPREF,LSOLVE(10),LPRINT(13),LBLK(10),MODE,NTIMES(10),
     &  RHOCON,SPACE,RATIO,SWIDTH, BWIDTH,THT !齿数,间隙,隙齿比,宽度
      REAL BRIDGEXL,BRIDGEXR !临时变量
      COMMON/CNTL/LSTOP,LEND
      COMMON ITIME
      COMMON/SORC/SMAX,SSUM,RSMAX
      COMMON/COEF/FLOW,DIFF,ACOF
      DIMENSION U(NX,NY),V(NX,NY),PC(NX,NY)
      EQUIVALENCE (F(1,1,1),U(1,1)),(F(1,1,2),V(1,1)),(F(1,1,3),PC(1,1))
      DIMENSION TH(NXY),THU(NXY),THDIF(NXY),THCV(NXY),THCVS(NXY)
      EQUIVALENCE(X,TH),(XU,THU),(XDIF,THDIF),(XCV,THCV),
     &  (XCVS,THCVS),(XL,THL)
***********************************************************************
*                                                                     *
*        Example-5: IMPINGING FLOW WITH A ROTATING PLATE              *
*                                                                     *
***********************************************************************
      DIMENSION WR(NX,NY)
      EQUIVALENCE(F(1,1,4),WR(1,1))
      DATA UIN,OMEGA/20.,300./ !数值过大会使解变成NAN原因未知
      
      ENTRY GRID
        XL=75.E-3
        YL=5.4E-3
        L1=520
        M1=56
        SPACE=1.0E-3 !间隙<YL
        RATIO=2.     !宽窄通道的比值>1
        THT=6.       !梯形齿倾斜角的正切值，不要取太小>5，那样齿会很矮
	  CALL UGRID
      RETURN
C----------------------------------------------------------------------
      ENTRY START
        MODE=2
        R(1)=200.E-3
        
        LAST=3000
        RSMAX=1.E-11
	          
        TSTART=0.
	  DT=1.E10
	  MTIME=1
	  ITIME=1
	  TIME=TSTART
	  LEND=.FALSE.

        DO I=1,4
          LSOLVE(I)=.TRUE.
          LPRINT(I)=.TRUE.
          LBLK(I)=.FALSE.
        END DO
	  RELAX(1)=0.5
	  RELAX(2)=0.5	  
	  RELAX(4)=0.5
	  RELAX(11)=0.8	  
	  TITLE(1)='.VEL-U.'
	  TITLE(2)='.VEL-V.'
	  TITLE(3)='.STRF.'	  
	  TITLE(4)='.VEL-WR.'

        DO I=1,L1
          DO J=1,M1
            U(I,J)=0.
            V(I,J)=0.
            WR(I,J)=0.
	    END DO
	  END DO
	    
	  DO K=1,NFMAX
	    DO I=1,L1
	      DO J=1,M1
	        FOLD(I,J,K)=F(I,J,K)
	      END DO
	    END DO
	  END DO
        TIME=TSTART+DT
      RETURN
C----------------------------------------------------------------------
      ENTRY DENSE
        DO I=1,L1
          DO J=1,M1
            RHO(I,J)=1.176
          END DO
        END DO
      RETURN
C----------------------------------------------------------------------
      ENTRY BOUND

	  DO J=1,M1
          IF(R(J) >= (R(1)+YL-SPACE)) THEN
	      U(2,J)=UIN ! 速度入口
            U(L1,J)=U(L2,J) !速度出口
            V(1,J)=0.
            V(L1,J)=0.
            WR(1,J)=0.
            WR(L1,J)=WR(L2,J)
          ELSE
            U(2,J)=0.
            U(L1,J)=0.
            V(1,J)=0.
            V(L1,J)=0.
            WR(1,J)=OMEGA*(R(J)**2.)
            WR(L1,J)=OMEGA*(R(J)**2.)
          END IF
	  END DO

	  DO I=1,L1
	    U(I,1)=0
	    V(I,2)=0.
	    WR(I,1)=0
	    U(I,M1)=0.
	    V(I,M1)=0.
	    WR(I,M1)=0         
	  END DO
       
        SWIDTH=XL/(RATIO*(NUM-1)+NUM)
        BWIDTH=RATIO*SWIDTH
        DY=YL/FLOAT(M1-2) !Y方向网格尺度
        DX=DY/THT ! 梯形齿单位高度X方向偏移

        ! 保证第一齿和末尾齿是直尺
        BRIDGEXL=0
        DO K=1,NUM
          BRIDGEXR=BRIDGEXL+SWIDTH
            DO J=1,M1
              BL(K,J)=BRIDGEXL
              RB(K,J)=BRIDGEXR
            END DO
          BRIDGEXL=BRIDGEXR+BWIDTH
        END DO

        ! 中间的齿是梯形齿，覆盖直齿，去掉此段，则全为直齿
        BRIDGEXL=SWIDTH+BWIDTH
        DO K=2,NUM-1
          BRIDGEXR=BRIDGEXL+SWIDTH
          BRIDGEXL1=BRIDGEXL
          BRIDGEXR1=BRIDGEXR
          DO J=1,M1
            BL(K,J)=BRIDGEXL1
            RB(K,J)=BRIDGEXR1
            BRIDGEXL1=BRIDGEXL1+DX
            BRIDGEXR1=BRIDGEXR1-DX
            ! PRINT*,DX,DY,BL(K,J),RB(K,J)
          END DO
          BRIDGEXL=BRIDGEXR+BWIDTH
        END DO

        ! 密封齿内全为边界，并设置旋转速度
        DO K=1,NUM
        DO I=1,L1                              !轴封内速度为0
	  DO J=1,M1
        IF(X(I)>=BL(K,J).AND.X(I)<RB(K,J).AND.R(J)<(R(1)+YL-SPACE))THEN
	          U(I,J)=0.
	          V(I,J)=0.
                WR(I,J)=OMEGA*(R(J)**2.)
        END IF
        END DO
	  END DO
        END DO   

      RETURN
C----------------------------------------------------------------------
      ENTRY OUTPUT
        IF(ITER.GE.2000) THEN
        FLOWIN=0.
        FLOWOUT=0.
        DO J=1,6
          FLOWIN=FLOWIN+U(2,J)*ARX(J)
        END DO          
        DO I=47,L1
          FLOWOUT=FLOWOUT+V(I,M1)*XCV(I)*R(M1)
        END DO
        FACTOR=FLOWOUT/FLOWIN
        DO I=47,L1
          V(I,M1)=V(I,M1)/FACTOR
        END DO
        END IF
        
        IF(ITER.EQ.1) WRITE(8,400)
        IF(MOD(ITER,20).EQ.0) THEN
          WRITE(8,403) ITER,P(1,52),P(L2,52),WR(25,25),SMAX
          WRITE(*,403) ITER,P(1,52),P(L2,52),WR(25,25),SMAX          
        END IF

  400   FORMAT(13X,'ITER',12X,'U(25,25)',12X,'V(25,25)',12X,
     &                       'WR(25,25)',12X,'SMAX')
  403   FORMAT(1X,I8,4F20.10)  
  
        IF(LSTOP.OR.(ITER.EQ.LAST)) CALL PRINT_RESULT        
      
        IF(LSTOP.OR.(ITER.EQ.LAST)) THEN
          OPEN(1,FILE='STREAM_FUNCTION.dat')
          WRITE(1,*) 'TITLE="SF"'
          WRITE(1,*) 'VARIABLES="X","Y","SF"'
          WRITE(1,*) 'ZONE I=',L1,',J=',2*M1-1
          DO J=M1,1,-1  
           DO I=1,L1
             WRITE(1,*) X(I)*1.E3,R(J)*1.E3,F(I,J,3)
            END DO
          END DO
          DO J=2,M1  
           DO I=1,L1
             WRITE(1,*) X(I)*1.E3,-R(J)*1.E3,F(I,J,3)
            END DO
          END DO
          CLOSE(1)
        END IF
    
         IF(LSTOP.OR.(ITER.EQ.LAST)) THEN  
          OPEN(1,FILE='VELOCITY.dat')            !速度大小
          WRITE(1,*) 'TITLE="VEL"'
          WRITE(1,*) 'VARIABLES="X","Y","U"'
          WRITE(1,*) 'ZONE I=',L1,',J=',M1
          DO J=1,M1  
           DO I=1,L1
             WRITE(1,*) X(I),R(J),sqrt(U(I,J)**2+V(I,J)**2)
            END DO
          END DO
          CLOSE(1)
         END IF   
           
        IF(LSTOP.OR.(ITER.EQ.LAST)) THEN  
          OPEN(1,FILE='PATHLINE.dat')            !速度分量
          WRITE(1,*) 'TITLE="VEL"'
          WRITE(1,*) 'VARIABLES="X","Y","U","V"'
          WRITE(1,*) 'ZONE I=',L1,',J=',M1
          DO J=1,M1  
           DO I=1,L1
             WRITE(1,*) X(I),R(J),U(I,J),V(I,J)
            END DO
          END DO
          CLOSE(1)
         END IF

        IF(LSTOP.OR.(ITER.EQ.LAST)) THEN 
          OPEN(1,FILE='PRESSURE.dat')             !压力
          WRITE(1,*) 'TITLE="P"'
          WRITE(1,*) 'VARIABLES="X","Y","P"'
          WRITE(1,*) 'ZONE I=',L1,',J=',M1
          DO J=1,M1  
           DO I=1,L1
             WRITE(1,*) X(I),R(J),P(I,J)
            END DO
          END DO
          CLOSE(1)
         END IF

         IF(LSTOP.OR.(ITER.EQ.LAST)) THEN 
          OPEN(1,FILE='BLRB.dat')             !齿位置
          WRITE(1,*) 'TITLE="B"'
          WRITE(1,*) 'VARIABLES="X","Y","BL","RB"'
          WRITE(1,*) 'ZONE K=',NUM,',J=',M1
          DO K=1,NUM  
           DO J=1,M1
             WRITE(1,*) K,J,BL(K,J),RB(K,J)
          END DO
          END DO
          CLOSE(1)
         END IF


        
      RETURN
C----------------------------------------------------------------------
      ENTRY GAMSOR
        
        DO I=1,L1
          DO J=1,M1
            GAM(I,J)=1.862E-5
          END DO
        END DO
        
        IF(NF.EQ.2) THEN
          DO J=3,M2
            DO I=2,L2
              RSWM=FY(J)*WR(I,J)+FYM(J)*WR(I,J-1)
              CON(I,J)=RHO(I,J)*(RSWM**2.)/(RMN(J)**3.)
              AP(I,J)=-2.*GAM(I,J)/(RMN(J)**2.)
            END DO
          END DO
        END IF
        
        IF(NF.EQ.4) THEN
          DO J=2,M2
            DO I=2,L2
              TEMP=2.*GAM(I,J)/(R(J)*YDIF(J))
              CON(I,J)=TEMP*WR(I,J-1)
              AP(I,J)=-TEMP
            END DO
          END DO
        END IF
      RETURN
C----------------------------------------------------------------------
      ENTRY TSTEP
        ITER=0
        LSTOP=.FALSE.
        IF (ITIME.GE.MTIME) LEND=.TRUE.
        TIME=TIME+DT
        ITIME=ITIME+1
      RETURN
C----------------------------------------------------------------------
      ENTRY STEPUP
        DO K=1,NFMAX
          DO J=1,M1
            DO I=1,L1
              FOLD(I,J,K)=F(I,J,K)
            END DO
          END DO
        END DO
      RETURN
C----------------------------------------------------------------------
      ENTRY ISSTOP
        IF(LSOLVE(1)) THEN
          IF (SMAX<RSMAX) LSTOP=.TRUE.
        END IF
      RETURN
C----------------------------------------------------------------------
      END