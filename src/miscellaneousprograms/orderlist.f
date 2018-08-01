      program orderlist
      implicit double precision(a-h,o-z)

      real*8, allocatable :: ch(:),fr(:),nb(:),ab(:),po(:)

      real*8, allocatable :: ti(:)

      character*30, allocatable  :: cch(:),cfr(:),cnb(:),cab(:),cpo(:)

      character*30, allocatable  :: ca(:)

      integer, allocatable  :: indx(:)

      character*80 abset(100)

      character*1 disp

      character*20 ca1

      scalc=3.d0

      open(unit=1,file='IN_PACKAGE',status='old')
      read(1,*)nqmprog
      close(unit=1)

      nchrgs=0
      open(unit=1,file='IN_CHARGES',status='old')
      read(1,*,end=333)
      read(1,*)nchrgs
333   close(unit=1)

c find out if dispersion is to be calculated
      open(unit=1,file='ABSETTINGS',status='old')
      n=1
1     read(1,100,end=2)abset(n)
100   format(a80)
      n=n+1
      go to 1
2     continue
      nabset=n-1
      do i=1,nabset
      if(trim(abset(i)).eq.'Is long range dispersion accounted for?')
     . then
       disp=trim(abset(i+1))
       go to 3
      endif
      enddo
3     continue

      open(unit=1,file='OUT_ELECTRONS',status='old')
      read(1,*)
      read(1,*)nch
      if(nch.gt.0)then
       allocate(ch(nch))
       allocate(cch(nch))
       read(1,*)
       do i=1,nch
        read(1,*)i1,ch(i),cch(i)
        ch(i)=-ch(i)**scalc
       enddo
      endif
      read(1,*)
      read(1,*)nfr
      if(nfr.gt.0)then
       allocate(fr(nfr))
       allocate(cfr(nfr))
       read(1,*)
       do i=1,nfr
        read(1,*)i1,fr(i),cfr(i)
       enddo
      endif
      read(1,*)
      read(1,*)nnb
      if(nnb.gt.0)then
       allocate(nb(nnb))
       allocate(cnb(nnb))
       read(1,*)
       do i=1,nnb
        read(1,*)i1,nb(i),cnb(i)
       enddo
      endif
      read(1,*)
      read(1,*)npo
      if(npo.gt.0)then
       allocate(po(npo))
       allocate(cpo(npo))
       read(1,*)
       do i=1,npo
        read(1,*)i1,po(i),cpo(i)
       enddo
      endif
      read(1,*)
      read(1,*)nab
      if(nab.gt.0)then
       allocate(ab(nab))
       allocate(cab(nab))
       read(1,*)
       do i=1,nab
        read(1,*)i1,ab(i),cab(i)
       enddo
      endif

      close(unit=1)

c out data needed tp parallelize the input deck construction
      open(unit=1,file='JOBSDATA',status='unknown')
      write(1,*)nch
      write(1,*)nfr
      write(1,*)nab
      write(1,*)nnb
      write(1,*)npo
      write(1,*)nqmprog
      write(1,334)disp
334   format(a1)
      close(unit=1)     

      if(nqmprog.eq.2.or.nqmprog.eq.4)then
       if(disp.eq."Y".or.nchrgs.eq.0)then
        ntot=nfr+nnb+npo+nab
       else
        ntot=nfr+nnb+nab
       endif
      else
       ntot=nfr+nnb+nab
      endif

      allocate(ti(ntot))
      allocate(ca(ntot))
      allocate(indx(ntot))


      if(nfr.gt.0)then
       do i=1,nfr
        ti(i)=-fr(i)**scalc
        ca(i)=cfr(i)
       enddo
      endif

      if(nnb.gt.0)then
       do i=1,nnb
        ti(nfr+i)=-nb(i)**scalc
        ca(nfr+i)=cnb(i)
       enddo
      endif

      nsofar=nfr+nnb

      if(nqmprog.eq.2.or.nqmprog.eq.4)then
       if(disp.eq."Y".or.nchrgs.eq.0)then
        if(npo.gt.0)then
         do i=1,npo
          ti(nsofar+i)=-2.5*po(i)**scalc
          ca(nsofar+i)=cpo(i)
         enddo
        endif
        nsofar=nsofar+npo
       endif
      endif

      if(nab.gt.0)then
       do i=1,nab
        ti(nsofar+i)=-ab(i)**scalc
        ca(nsofar+i)=cab(i)
       enddo
      endif

      call indexx(ntot,ti,indx)

      open(unit=1,file='INLIST',status='unknown')
      do i=1,ntot
       write(1,*)ca(indx(i))
       write(8,*)(-ti(indx(i)))**(1.d0/scalc)
      enddo
      close(unit=1)



      nsofar=0

      if(npo.eq.0)go to 4

      if(disp.eq.'Y')then
       do i=1,npo
        ti(i)=-po(i)**scalc
        call filelabel(i,ca1)
        n1=index(ca1,' ')-1
        ca(i)='grp.'//ca1(1:n1)//'.mol'
       enddo
       nsofar=npo
      endif


      if(nqmprog.eq.1)then
       do i=1,npo
        ti(nsofar+i)=-po(i)**scalc
        n1=index(cpo(i),'.com')-1
        ca(nsofar+i)=cpo(i)(1:n1)//'.mol'
       enddo
       nsofar=nsofar+npo
      endif

      if(nqmprog.eq.3)then
       if(disp.eq."Y".or.nchrgs.eq.0)then
        do i=1,npo
         ti(nsofar+i)=-po(i)**scalc
         n1=index(cpo(i),'.com')-1
         ca(nsofar+i)=cpo(i)(1:n1)//'.mol'
        enddo
        nsofar=nsofar+npo
       endif
      endif

      if(nsofar.eq.0)go to 4
       if(nsofar.eq.1)then
        indx(1)=1
       else
        call indexx(nsofar,ti,indx)
       endif

      open(unit=1,file='INLISTDAL',status='unknown')
       if(nsofar.gt.0)then
        do i=1,nsofar
         write(1,*)ca(indx(i))
        enddo
       endif
       close(unit=1)


4     continue

      if(nch.eq.0)go to 5
      if(nch.gt.1)then
       call indexx(nch,ch,indx)
      else
       indx(1)=1
      endif

       open(unit=1,file='INLISTCHG',status='unknown')
       do i=1,nch
        write(1,*)cch(indx(i))
       enddo
       close(unit=1)

5     continue

      end


      SUBROUTINE INDEXX(N,ARRIN,INDX)
      implicit real*8 (a-h,o-z)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
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
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
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

      subroutine filelabel(k,ca)
      implicit double precision(a-h,o-z)

c this subroutines returns a character variable ca
c for the integer k

      character*20 ca
      character*1 I0(0:10),anum(10),ca1,ca2,ca3,ca4,ca5

      I0(10)='0'
      I0(0)='0'
      I0(1)='1'
      I0(2)='2'
      I0(3)='3'
      I0(4)='4'
      I0(5)='5'
      I0(6)='6'
      I0(7)='7'
      I0(8)='8'
      I0(9)='9'

      if(k.le.9)then
       ca=I0(k)
      endif

      if(k.gt.9.and.k.le.99)then
      k1=k/10
      k2=k-k1*10
      if(k2.eq.0)k2=10
      ca1=I0(k1)
      ca2=I0(k2)
      ca=ca1//ca2
      endif

      if(k.gt.99.and.k.le.999)then
      k1=k/100
      k3=k/10 - INT(k1)*10
      k4=k - INT(k1)*100 - INT(k3)*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca=ca1//ca2//ca3
      endif

      if(k.gt.999.and.k.le.9999)then
      k1=k/1000
      k3=k/100 - INT(k1)*10
      k4=k/10 - INT(k1)*100 - INT(k3)*10
      k5=k-INT(k1)*1000-k3*100-k4*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca4=I0(k5)
      ca=ca1//ca2//ca3//ca4
      endif

      if(k.gt.9999.and.k.le.99999)then
      k1=k/10000
      k3=k/1000 - INT(k1)*10
      k4=k/100 - INT(k1)*100 - INT(k3)*10
      k5=k/10-INT(k1)*1000-k3*100-k4*10
      k6=k-INT(k1)*10000-k3*1000-k4*100-k5*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca4=I0(k5)
      ca5=I0(k6)
      ca=ca1//ca2//ca3//ca4//ca5
      endif

      return
      end

