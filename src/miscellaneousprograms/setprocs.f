      program nprocs
      implicit double precision(a-h,o-z)

      real*8, allocatable :: ch(:),fr(:),nb(:),ab(:),po(:)

      integer, allocatable :: calc(:)

      character*80 abset(100)

      character*1 disp

c  this program simply tries to estimate the number of processors
c  required for single processor jobs
c  using the number of electrons ** scale as an estimate of the
c time for each calculation


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
       read(1,*)
       do i=1,nch
        read(1,*)i1,ch(i)
       enddo
      endif
      read(1,*)
      read(1,*)nfr
      if(nfr.gt.0)then
       allocate(fr(nfr))
       read(1,*)
       do i=1,nfr
        read(1,*)i1,fr(i)
       enddo
      endif
      read(1,*)
      read(1,*)nnb
      if(nnb.gt.0)then
       allocate(nb(nnb))
       read(1,*)
       do i=1,nnb
        read(1,*)i1,nb(i)
       enddo
      endif
      read(1,*)
      read(1,*)npo
      if(npo.gt.0)then
       allocate(po(npo))
       read(1,*)
       do i=1,npo
        read(1,*)i1,po(i)
       enddo
      endif
      read(1,*)
      read(1,*)nab
      if(nab.gt.0)then
       allocate(ab(nab))
       read(1,*)
       do i=1,nab
        read(1,*)i1,ab(i)
       enddo
      endif

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

      allocate(calc(ntot))

      scalc=3.0d0

      if(nfr.gt.0)then
       do i=1,nfr
        calc(i)=idint(fr(i)**scalc)
       enddo
      endif
      nsofar=nfr
      if(nnb.gt.0)then
       do i=1,nnb
        calc(nsofar+i)=idint(nb(i)**scalc)
       enddo
       nsofar=nsofar+nnb
      endif

      if(nqmprog.eq.2.or.nqmprog.eq.4)then
       if(disp.eq."Y".or.nchrgs.eq.0)then
        if(npo.gt.0)then
         do i=1,npo
          calc(nsofar+i)=idint(2.5*po(i)**scalc)
         enddo
        endif
        nsofar=nsofar+npo
       endif
      endif

      if(nab.gt.0)then
       do i=1,nab
        calc(nsofar+i)=idint(ab(i)**scalc)
       enddo
      endif

      n=nsofar+nab

      write(6,*)
      write(6,*)' The total number of ab initio calculations'
      write(6,*)' to be carried out is'
      write(6,*) n
      write(6,*)

c sort njobs jobs into time order, longest time last
      call piksrt(n,calc)
      call alloc(n,calc,nproc)

      write(6,*)
      write(6,*)'           PARALLEL IMPLEMENTATION'
      write(6,*)
      write(6,*)' The bulk of the cputime is employed in performing'
      write(6,*)' electronic structure calculations for the many'
      write(6,*)' molecular fragments, and for fragment properties that'
      write(6,*)' are needed to evaluate the non-bonded interactions.'
      write(6,*)' For these calculations, the recommended number of'
      write(6,*)' cpus for maximum efficiency and and minimum walltime'
      write(6,*)' is'
      write(6,*)
      write(6,*) nproc
      write(6,*)
      write(6,*)' This assumes that each fragment uses only 1 cpu.'
      write(6,*)
c      write(6,*)' If you allocate N cpus to each fragment, then the'
c      write(6,*)' recommended number of cpus for maximum efficiency'
c      write(6,*)' and minimum walltime is N times the number above.'
c      write(6,*)
      write(6,*)' The user sets the number of cpus in main menu item 4'
      write(6,*)
      write(6,*)' The user can specify less than ',nproc,' cpus,'
      write(6,*)' for example, if only a smaller number is available.'
      write(6,*)

      nprocreg=nproc

      nprocch=0
      nprocdal=0
      if(nch.eq.0)go to 200
      if(nch.eq.1)then
       nprocch=1
       go to 200
      endif

      do i=1,nch
       calc(i)=idint(ch(i)**3)
      enddo

      n=nch
      call piksrt(n,calc)
      call alloc(n,calc,nproc)

      if(nproc.gt.nprocreg)nproc=nprocreg

      nprocch=nproc
      if(nch.gt.nprocreg)nprocch=nprocreg
c new code 250618
      if(nprocch.gt.nch)nprocch=nch

c     write(6,*)' There are formally charged groups in the molecule,'
c     write(6,*)' and SMFA must perform calculations to determine the'
c     write(6,*)' charge distribution in all such groups.'
c     write(6,*)' SMFA will use ',nprocch,' cpus for these calculations'
c     write(6,*)' (one cpu for each group), unless this number exceeds'
c     write(6,*)' the limit placed on cpus in main menu item 4.'
c     write(6,*)

200   continue

c now dalton jobs

      if(npo.eq.0)go to 201
      nsofar=0
      if(disp.eq.'Y')then
       do i=1,npo
        calc(i)=idint(po(i)**3)
       enddo
       nsofar=npo
      endif
      if(nqmprog.eq.1.or.nqmprog.eq.3)then
       do i=1,npo
        calc(nsofar+i)=idint(po(i)**3)
       enddo
       nsofar=nsofar+npo
      endif

      if(nsofar.eq.0)go to 201

      if(nsofar.eq.1)then
       nprocdal=1
       go to 201
      endif

      call piksrt(nsofar,calc)
      call alloc(nsofar,calc,nproc)

      if(nproc.gt.nprocreg)nproc=nprocreg

      nprocdal=nproc
      if(nsofar.ge.nprocreg)nprocdal=nprocreg

c     write(6,*)' SMFA needs to perform calculations of dispersion'
c     write(6,*)' coefficients and (possibly) polarizabilities for '
c     write(6,*)' each functional group, using DALTON.'
c     write(6,*)' SMFA will use ',nprocdal,' cpus for these'
c     write(6,*)' calculations (one cpu for each group), unless this'
c     write(6,*)' number exceeds the limit specified in menu item 4.'

201   continue

      open(unit=1,file='NPROCS',status='unknown')
      write(1,*)nprocreg
      write(1,*)nprocch
      write(1,*)nprocdal
      close(unit=1)

      end


      subroutine piksrt(n,arr)
      integer arr(n)
      do 12 j=2,n
        a=arr(j)
        do 11 i=j-1,1,-1
         if(arr(i).le.a) go to 10
         arr(i+1)=arr(i)
11      continue
       i=0
10    arr(i+1)=a
12    continue
      return
      end

      subroutine alloc(n,calc,nproc)
      implicit double precision(a-h,o-z)
      integer calc(n),ifin(n),timelim,thistime

      njobs=n
      do i=1,n
       ifin(i)=0
      enddo

c longest job in proc1, second longest in proc2
      timelim=calc(njobs)
      ifin(njobs)=1
      ifin(njobs-1)=1
      nproc=2
      last=njobs-1
      next=last-1
c      write(7,*)timelim
13    thistime=calc(last)
      do j=next,1,-1
       if(ifin(j).eq.1)go to 10
       if(calc(j)+thistime.lt.timelim)then
        thistime=thistime+calc(j)
        ifin(j)=1
       endif
10    enddo
c     write(7,*)thistime
      nproc=nproc+1
      do k=last-1,1,-1
       if(ifin(k).eq.0)go to 11
      enddo
11    if(k.gt.1)then
       last=k
       next=last-1
      else
       go to 12
      endif
      go to 13
12    continue

      return
      end
