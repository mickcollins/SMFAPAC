      SUBROUTINE INDEXI(N,IARRIN,INDX)
      implicit real*8 (a-h,o-z)
      DIMENSION IARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          IQ=IARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          IQ=IARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(IARRIN(INDX(J)).LT.IARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(IQ.LT.IARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END
