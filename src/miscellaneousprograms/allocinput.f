      program allocinput
      implicit double precision(a-h,o-z)

      integer, allocatable :: m1(:,:),m2(:,:),ndel(:)

      dimension njobs(5)
      character*1 disp
      character*80 code,perlsc

c     read(5,*)nchgs
c     read(5,*)nfrags
c     read(5,*)nab
c     read(5,*)nnb
c     read(5,*)npol
c     read(5,*)nqm
c     read(5,*)disp
c     read(5,*)ncpus

      open(unit=1,file='CODENAME',status='old')
      read(1,80)code
      close(unit=1)
80    format(a80)

      na1=index(code,' ')-1

      do i=1,5
       read(5,*)njobs(i)
      enddo
      read(5,*)nqm
      read(5,*)disp
      read(5,*)ncpus
      allocate(m1(5,ncpus))
      allocate(m2(5,ncpus))
      allocate(ndel(ncpus))

      do j=1,ncpus
      do i=1,5
       m1(i,j)=0
       m2(i,j)=0
      enddo
      enddo

      do i=1,5

       if(njobs(i).gt.0)then
        n1=njobs(i)/ncpus
        n2=njobs(i)-n1*ncpus

        do j=1,ncpus
        ndel(j)=n1
        if(j.le.n2)ndel(j)=ndel(j)+1
        enddo

        m1(i,1)=1
        m2(i,1)=ndel(1)
        do j=2,ncpus
         m1(i,j)=m2(i,j-1)+1
         m2(i,j)=m1(i,j)+ndel(j)-1
         if(m2(i,j).gt.njobs(i))then
          m2(i,j)=njobs(i)
          go to 1
         endif
         if(m1(i,j).gt.njobs(i))then
          m1(i,j)=0
          m2(i,j)=0
          go to 1
         endif
        enddo
1       continue
       endif
      enddo

100   format(a65,6i6,a3)
101   format(a65,8i6,a3)
102   format(a65,2i6,a3)

      if(nqm.eq.1)then
       perlsc=code(1:na1)//'/SMFAgaminputs_MAC.pl '
       na2=index(perlsc,' ')
       
       open(unit=1,file='gaminput',recl=130,status='unknown')
       do j=1,ncpus
c      write(1,100)'SMFAgaminputs_MAC.pl ',(m1(i,j),m2(i,j),i=2,4),'  &'
       write(1,100)perlsc(1:na2),(m1(i,j),m2(i,j),i=2,4),'  &'
       enddo
       write(1,*) "wait"
       close(unit=1)
      endif

      if(nqm.eq.2)then
       perlsc=code(1:na1)//'/SMFAgauinputs_MAC.pl '
       na2=index(perlsc,' ')
       open(unit=1,file='gauinput',recl=130,status='unknown')
       do j=1,ncpus
c      write(1,101)'SMFAgauinputs_MAC.pl ',(m1(i,j),m2(i,j),i=2,5),'  &'
       write(1,101)perlsc(1:na2),(m1(i,j),m2(i,j),i=2,5),'  &'
       enddo
       write(1,*) "wait"
       close(unit=1)
       perlsc=code(1:na1)//'/SMFAmkchinputs_gau.pl '
       na2=index(perlsc,' ')
       open(unit=1,file='gauchinput',recl=130,status='unknown')
       do j=1,ncpus
        if(m1(1,j).gt.0)then
c         write(1,102)'SMFAmkchinputs_gau.pl ', m1(1,j),m2(1,j),'  &'
         write(1,102)perlsc(1:na2),m1(1,j),m2(1,j),'  &'
        endif
       enddo
       write(1,*) "wait"
       close(unit=1)
      endif

      if(nqm.eq.3)then
       perlsc=code(1:na1)//'/SMFAnwcinputs_MAC.pl '
       na2=index(perlsc,' ')
       open(unit=1,file='nwcinput',recl=130,status='unknown')
       do j=1,ncpus
c       write(1,100)'SMFAnwcinputs_MAC.pl ',(m1(i,j),m2(i,j),i=2,4),'  &'
        write(1,100)perlsc(1:na2),(m1(i,j),m2(i,j),i=2,4),'  &'
       enddo
       write(1,*) "wait"
       close(unit=1)
       perlsc=code(1:na1)//'/SMFAmkchinputs_nwc.pl '
       na2=index(perlsc,' ')
       open(unit=1,file='nwcchinput',recl=130,status='unknown')
       do j=1,ncpus
        if(m1(1,j).gt.0)then
c         write(1,102)'SMFAmkchinputs_nwc.pl', m1(1,j),m2(1,j),'  &'
         write(1,102)perlsc(1:na2),m1(1,j),m2(1,j),'  &'
        endif
       enddo
       write(1,*) "wait"
       close(unit=1)
      endif

      if(nqm.eq.4)then
       perlsc=code(1:na1)//'/SMFAqchinputs_MAC.pl '
       na2=index(perlsc,' ')
       open(unit=1,file='qchinput',recl=130,status='unknown')
       do j=1,ncpus
c       write(1,101)'SMFAqchinputs_MAC.pl ',(m1(i,j),m2(i,j),i=2,5),'  &'
        write(1,101)perlsc(1:na2),(m1(i,j),m2(i,j),i=2,5),'  &'
       enddo
       write(1,*) "wait"
       close(unit=1)
       perlsc=code(1:na1)//'/SMFAmkchinputs_qch.pl '
       na2=index(perlsc,' ')
       open(unit=1,file='qchchinput',recl=130,status='unknown')
       do j=1,ncpus
        if(m1(1,j).gt.0)then
c         write(1,102)'SMFAmkchinputs_qch.pl', m1(1,j),m2(1,j),'  &'
         write(1,102)perlsc(1:na2),m1(1,j),m2(1,j),'  &'
        endif
       enddo
       write(1,*) "wait"
       close(unit=1)
      endif


c dalton input

       end
