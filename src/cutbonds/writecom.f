      subroutine writecom
      use fractdata
      implicit double precision(a-h,o-z)

      integer, allocatable   :: natf(:),nfat(:),nalloc(:,:)
      real*8, allocatable    :: cf(:,:),cfraction(:,:)
      character*2, allocatable :: fat(:)

      dimension xl(3),xn(3)
      character*1 I0(0:10),anum(10),ca1,ca2,ca3,ca4,ca5
      character*20 ca,caFRAG

      allocate(natf(nffinal))
      allocate(nfat(2*natomall))
      allocate(cf(2*natomall,3))
      allocate(fat(2*natomall))
      allocate(nalloc(2*natomall,2))
      allocate(cfraction(2*natomall,2))

      fact1=0.68d0

      open(unit=20,file='OUT_ALLOC_FRAG',status='unknown')
      write(20,*)' this files contains the list of atom numbers '
      write(20,*)' in the original molecule that'
      write(20,*)' correspond to the atomic coordinates in the '
      write(20,*)' input files for each'
      write(20,*)' fragment whose energy is to be evaluated ab initio'
      write(20,*)' Original atoms have their number followed '
      write(20,*)' by "0    1.    0."'
      write(20,*)' Capping H atoms are listed as two original '
      write(20,*)' atoms followed by '
      write(20,*)' the relative allocations. For example'
      write(20,*)' "23      32      0.7d0   0.3d0"'



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


c link atoms will be added to a fragment when the "bonded" atom
c is missing from the junk array

      write(3,*)'The capping atoms are attached to'

      open(unit=4,file='molecule.com',status='unknown')

c this version accounts for charge
c First find the net charge on the whole molecule

      itchg=0
      numtot=0
      do i=1,natomall
      itchg=itchg+ichg(i)
      numtot=numtot+numa(i)
      enddo
c if numtot is even we make the molecule a singlet,
c if numtot is odd, we make the molecule a doublet
c correct for charge


      numtot=numtot-itchg
      numtot2=numtot/2
      ndg=numtot-2*numtot2
      if(ndg.eq.0)then
       ndegen=1
      else
       ndegen=2
      endif

c If Level  > 0, we will output the number of electrons in the whole
c molecule and in each of the fragments
      open(unit=99,file='ONEONLY',status='old')
      read(99,*)Levelflag
      close(unit=99)
      if(Level.gt.1.or.(Levelflag.eq.2))then
c       open(unit=22,file='OUT_ELECTRONS',status='unknown')
       open(unit=22,file='OUT_ELECTRONS',status='old',action='write',
     .        form='formatted',position="append")

       open(unit=23,file='OUT_ELECTRONS_SUMMARY',status='unknown')
c        write(22,*)'Total numbers of atoms and electrons: ',natomall,numtot
         write(23,*)'Total numbers of atoms and electrons: ',natomall,numtot
         write(22,*)' The number of fragments is'
         write(23,*)' The number of fragments is'
         write(22,*)nffinal
         write(23,*)nffinal
      endif

82    format(a80)
      do m=1,nabitio
      write(4,82)abinitio(m)
      enddo
      write(4,*)
      write(4,*)'Number of atoms is ',natomall
      write(4,*)'The whole molecule'
      write(4,*)
c     write(4,*)'0,1'
      write(4,*)itchg,",",ndegen

      do n=1,natomall
      write(4,100)atoms(n),(c(n,k),k=1,3)
      enddo
100   format(a2,3f13.6)
      write(4,*)

      close(unit=4)

c     write(4,101)'--link1--'
101   format(a9)

c make a file that has the number of atoms in the fragments
      open(unit=15,file='fragnum',status='unknown')
      write(15,8765)ititle
8765  format(1x,a20)
      numav=0
      numav2=0
      numcaps=0

      open(unit=42,file='OUT_NCOORD_FRAG',status='unknown')
      write(42,*)'Number of fragments is ',nffinal

      nemax=0.d0
      nmax=0
      write(22,*)'The numbers of atoms and electrons in each fragment'


c now write out the fragments
c NOTE that subroutine finalcancel put the +ve fragments first
      do k=1,nffinal

c  charge
      itchg=0
      numtot=0

      itot=0
      do i=1,numat(k)
      do j=1,nfam(natstore(k,i))
       itot=itot+1
       junk(itot)=ifam(natstore(k,i),j)
      enddo
      enddo
      itotf(k)=itot


      do m=1,itotf(k)
      n=junk(m)
      itchg=itchg+ichg(n)
      numtot=numtot+numa(n)

c     write(4,100)atoms(n),(c(n,k1),k1=1,3)
c record atom n in fragment k
c record coords of each atom in frag k
      do k1=1,3
      cf(m,k1)=c(n,k1)
      enddo
      fat(m)=atoms(n)
      nfat(m)=numa(n)
      nalloc(m,1)=n
      nalloc(m,2)=0
      cfraction(m,1)=1.d0
      cfraction(m,2)=0.d0
      enddo
      natf(k)=itotf(k)

c find the missing links
c and write out the linking H atoms

c     write(6,*)' Fragment ',k

      do m=1,itotf(k)
      n=junk(m)


      do i=1,nbondso

      if(multstore(i).lt.0)go to 44

      if(mbstore(i).eq.n)then
       match=0
       do j=1,itotf(k)
       if(j.eq.m)go to 20
       if(junk(j).eq.nbstore(i))match=1
20     continue
       enddo
       if(match.eq.0)then
c work out fact1
       sum1=0.23d0+radius(n)
       sum2=radius(n)+radius(nbstore(i))
       fact1=sum1/sum2
        do k1=1,3
        xl(k1)=c(n,k1)+fact1*(c(nbstore(i),k1)-c(n,k1))
        enddo
       write(3,*)k,mbstore(i),nbstore(i),fact1
       write(33,*)' nb ',k,radius(n),radius(nbstore(i))
        natf(k)=natf(k)+1
        numtot=numtot+1
        nalloc(natf(k),1)=n
        nalloc(natf(k),2)=nbstore(i)
        cfraction(natf(k),1)=1.d0-fact1
        cfraction(natf(k),2)=fact1

        do k1=1,3
        cf(natf(k),k1)=xl(k1)
        enddo
        fat(natf(k))='H '
        nfat(natf(k))=1
        numcaps=numcaps+1
       endif
      endif

      if(nbstore(i).eq.n)then
       match=0
       do j=1,itotf(k)
       if(j.eq.m)go to 21
       if(junk(j).eq.mbstore(i))match=1
21     continue
       enddo
       if(match.eq.0)then
       sum1=0.23+radius(n)
       sum2=radius(n)+radius(mbstore(i))
       fact1=sum1/sum2
        do k1=1,3
        xl(k1)=c(n,k1)+fact1*(c(mbstore(i),k1)-c(n,k1))
        enddo
       write(3,*)k,nbstore(i),mbstore(i),fact1
       write(33,*)' mb ',k,radius(n),radius(mbstore(i))
        natf(k)=natf(k)+1
        numtot=numtot+1
        nalloc(natf(k),1)=n
        nalloc(natf(k),2)=mbstore(i)
        cfraction(natf(k),1)=1.d0-fact1
        cfraction(natf(k),2)=fact1
        do k1=1,3
        cf(natf(k),k1)=xl(k1)
        enddo
        fat(natf(k))='H '
        nfat(natf(k))=1
        numcaps=numcaps+1
       endif
      endif
44    continue
c end the loop over bonds
      enddo
c  end the loop over atoms in the fragment
      enddo

c     write(8,*)

c if numtot is even we make the molecule a singlet,
c if numtot is odd, we make the molecule a doublet
c correct for charge
      numtot=numtot-itchg
      numtot2=numtot/2
      ndg=numtot-2*numtot2
      if(ndg.eq.0)then
       ndegen=1
      else
       ndegen=2
      endif


      write(15,*)num1,num2
      numav=numav+num1
      numav2=numav2+num2

c open a separate output file for each fragment
c calculate filename for each fragment
      if(k.le.9)then
       ca='COORD'//I0(k)
       caFRAG='FRAG'//I0(k)//'.com'
      endif
      if(k.gt.9.and.k.le.99)then
      k1=k/10
      k2=k-k1*10
      if(k2.eq.0)k2=10
      ca1=I0(k1)
      ca2=I0(k2)
      ca='COORD'//ca1//ca2
      caFRAG='FRAG'//ca1//ca2//'.com'
      endif
      if(k.gt.99.and.k.le.999)then
      k1=k/100
      k3=k/10 - INT(k1)*10
      k4=k - INT(k1)*100 - INT(k3)*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca='COORD'//ca1//ca2//ca3
      caFRAG='FRAG'//ca1//ca2//ca3//'.com'
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
      ca='COORD'//ca1//ca2//ca3//ca4
      caFRAG='FRAG'//ca1//ca2//ca3//ca4//'.com'
      endif

c added extra 041218
      if(k.gt.9999.and.k.le.99999)then
       k1=k/10000
       k3=k/1000 - int(k1)*10
       k4=k/100 - int(k1)*100 - int(k3)*10
       k5=k/10 - int(k1)*1000 - int(k3)*100 - int(k4)*10
       k6=k - int(k1)*10000 - int(k3)*1000 - int(k4)*100 - int(k5)*10
       ca1=I0(k1)
       ca2=I0(k3)
       ca3=I0(k4)
       ca4=I0(k5)
       ca5=I0(k6)
       ca='COORD'//ca1//ca2//ca3//ca4//ca5
       caFRAG='FRAG'//ca1//ca2//ca3//ca4//ca5//'.com'
      endif


      if(Level.gt.1.or.(Levelflag.eq.2))then
        write(22,*)natf(k),numtot,caFRAG
        if(numtot.gt.nemax)then
         nemax=numtot
         nmax=k
        endif
      endif

      open(unit=44,file=ca,status='unknown')
      write(44,*)itchg,',',ndegen
      do m=1,natf(k)
      if(fat(m).eq.'H ')then
      num1=num1+1
      ncf=ncf+1
      else
      num2=num2+1
      ncf=ncf+1
      endif
      write(44,103)fat(m),(cf(m,k1),k1=1,3)
      enddo
103   format(a2,3f25.16)
      close(unit=44)
      write(42,*)ncf

      write(20,*)' Allocations for ab initio fragemnt ',k
      write(20,*)k,natf(k)
      num1=0
      num2=0
      do m=1,natf(k)
      if(fat(m).eq.'H ')then
      num1=num1+1
      else
      num2=num2+1
      endif

      write(20,601)nalloc(m,1),nalloc(m,2),
     .             cfraction(m,1),cfraction(m,2)
      enddo
601    format(2i6,2F10.5)

c  end the loop over fragments
      enddo

      avnum=dble(float(numav))/dble(float(nffinal))
      avnum2=dble(float(numav2))/dble(float(nffinal))
      write(15,5566)avnum,avnum2,natom,numcaps
5566  format(1x,2f10.2,2i10)
      close(unit=15)
      close(unit=20)

c output geometries for illustrations


      close(unit=42)
      write(23,*)
      write(23,*)' The largest number of electrons in a fragment'
      write(23,*)' is ',nemax,' for fragment number ',nmax
      close(unit=22)
      close(unit=23)

      return
      end
