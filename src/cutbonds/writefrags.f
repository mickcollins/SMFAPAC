      subroutine writefrags

      use fractdata
      implicit double precision(a-h,o-z)

      integer, allocatable :: kgroup(:)
      character*20 ca,ca1

      open(unit=3,file='frags.out',status='unknown')
      open(unit=8,file='signs.out',status='unknown')
      open(unit=4,file='OUT_LIST_FRAG',status='unknown')

      write(4,*)' This file contains the coefficients that multiply the'
      write(4,*)' energies of each bonded fragment. The list below'
      write(4,*)' contains the fragment number of coefficient.'

      icap=0
      nffinal=0
      do i=1,nf
      if(nstop(i).le.1.and.isign(i).lt.0)icap=icap+1
      if(nstop(i).le.1)nffinal=nffinal+1
      enddo

      write(3,*)'Total number of fragments is'
      write(3,*)nffinal
      write(3,*)'Total number of overlapping fragments is'
      write(3,*)icap
      write(3,*)'The number of atoms in each fragment is'
      do i=1,nf
      if(nstop(i).le.1.and.isign(i).gt.0)then
      itot=0
      do k=1,numat(i)
       itot=itot+nfam(natstore(i,k))
      enddo
      write(3,*)itot
      write(4,*)isign(i)
      write(8,*)isign(i)
      endif
      enddo
      do i=1,nf
      if(nstop(i).le.1.and.isign(i).lt.0)then
      itot=0
      do k=1,numat(i)
       itot=itot+nfam(natstore(i,k))
      enddo
      write(3,*)itot
      write(4,*)isign(i)
      write(8,*)isign(i)
      endif
      enddo

      write(3,*)'The atoms in each fragment are'

c allocate junk
      allocate(junk(nf,natomall))
      allocate(itotf(nf))

      icc=0
      do i=1,nf
      if(nstop(i).le.1.and.isign(i).gt.0)then
      icc=icc+1
      itot=0
      do k=1,numat(i)
      do j=1,nfam(natstore(i,k))
      itot=itot+1
      junk(icc,itot)=ifam(natstore(i,k),j)
      enddo
      enddo
      itotf(icc)=itot
      write(3,*)(junk(icc,k),k=1,itot)
c     write(3,93)(junk(icc,k),k=1,itot)
      endif
      enddo
      do i=1,nf
      if(nstop(i).le.1.and.isign(i).lt.0)then
      icc=icc+1
      itot=0
      do k=1,numat(i)
      do j=1,nfam(natstore(i,k))
      itot=itot+1
      junk(icc,itot)=ifam(natstore(i,k),j)
      enddo
      enddo
      itotf(icc)=itot
      write(3,*)(junk(icc,k),k=1,itot)
c     write(3,93)(junk(icc,k),k=1,itot)
      endif
      enddo

93    format(20i4)

c get average number of groups
      avgroups=0
      do i=1,nffinal
       avgroups=avgroups+numat(i)
      enddo
       avgroups=avgroups/nffinal
       write(40,*)avgroups,natom

c write the fragments to a file
      open(unit=9,file='OUT_FRAGMENTS',status='unknown')
      write(9,*)nffinal
      do i=1,nf
      if(nstop(i).le.1.and.isign(i).gt.0)then
      write(9,*)numat(i),i
      write(9,*)(natstore(i,k),k=1,numat(i))
      endif
      enddo
      do i=1,nf
      if(nstop(i).le.1.and.isign(i).lt.0)then
      write(9,*)numat(i),i
      write(9,*)(natstore(i,k),k=1,numat(i))
      endif
      enddo
      close(unit=9)
c     close(unit=3)
      close(unit=4)
      close(unit=8)

c create the files which contain the list of charged groups to be
c used as embedded charges for each fragment

      if(Level.lt.1)go to 2113

      open(unit=20,file='OUT_CHARGEDGROUPS',status='old')
      read(20,*)
      read(20,*)kcharge
      if(kcharge.eq.0)go to 2111
      allocate(kgroup(kcharge))
      read(20,*)
      do k=1,kcharge
       read(20,*)kgroup(k)
       read(20,*)
       read(20,*)
       read(20,*)
      enddo

      do n=1,nffinal
       call filelabel(n,ca1)
       n1=index(ca1,' ')-1
       ca='chFRAG.'//ca1(1:n1)
       open(unit=21,file=ca,status='unknown')
       do m=1,kcharge
        do i=1,numat(n)
         if(natstore(n,i).eq.kgroup(m))go to 2112
        enddo
        write(21,*)kgroup(m)
2112   enddo
       close(unit=21)
      enddo

2111  close(unit=20)
2113  continue

      return
      end
