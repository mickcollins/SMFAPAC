      subroutine families

      use fractdata
      implicit double precision(a-h,o-z)

      integer, allocatable   :: ngp(:),kgroup(:),numchgs(:),natchg(:,:)
      integer, allocatable   :: formal(:)
      integer, allocatable   :: numatoms(:),numelects(:)

      dimension xcap(24,3),xmean(3)
      dimension nattach(24)

c assumes no more than 10 formal charges in a group

      character*1 I0(0:10),ca1,ca2,ca3,ca4
      character*15 ca
      character*15 da
      character*21 ea,fa,pa,panames(1000)
      character*27 ga

      integer, allocatable   ::nchgatoms(:)

      allocate(ngp(natomall))
      allocate(kgroup(natomall))
      allocate(numatoms(natomall))
      allocate(numelects(natomall))
      allocate(numchgs(natomall))
      allocate(natchg(natomall,maxfamily))
      allocate(formal(natomall))

      formal=0

      bohr=1.8897259886d0

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

c read in which atoms are "formally" charged
      open(unit=1,file='IN_CHARGES',status='old')
      read(1,*)
      read(1,*)numcharges
      if(numcharges.eq.0)go to 601
      allocate(nchgatoms(numcharges))
      read(1,*)
      do n=1,numcharges
       read(1,*)nchgatoms(n)
      enddo
601   continue
      close(unit=1)

c  this sub outputs the atoms in each family.
c  these are needed by findHbond and findHandCbond to
c  avoid connecting families multiple times by "-1" bonds

      open(unit=20,file='families.out',status='unknown')

      write(20,*)' The number of groups'
      write(20,*)natom
      write(20,*)' The number of atoms in each group & their identities'
      do n=1,natom
      write(20,*)nfam(n)
      write(20,*)(ifam(n,i),i=1,nfam(n))
      enddo

      close(unit=20)

c evaluate group number
      do n=1,natom
       do i=1,nfam(n)
        ngp(ifam(n,i))=n
       enddo
      enddo

      do n=1,natomall
       write(84,*)n,ngp(n)
      enddo

c if Level = 0 and formal charges exist, ".chg" are created for use
c as background charges

      if(Level.eq.0)then

c check to see if any groups are charged
      match=0
      do n=1,natom
      do i=1,nfam(n)
      if(iabs(ichg(ifam(n,i))).gt.0)then
       match=1
       formal(n)=1
      endif

c extra code to account for atoms in IN_CHARGES
      if(numcharges.gt.0)then
       do j=1,numcharges
        if(ifam(n,i).eq.nchgatoms(j))then
         match=1
         formal(n)=1
        endif
       enddo
      endif

      enddo
      enddo

c allow all groups to have npa charges calculated
c     match=1


      if(match.eq.1)then
       open(unit=90,file='CHARGECOORDS',status='unknown')
      endif

c keep count of the number of groups containing formal charges
      kcharge=0

      do n=1,natom

       nsum=0
       nchg=0
       name1=0
       do i=1,nfam(n)
       nsum=nsum+ichg(ifam(n,i))
       nchg=nchg+iabs(ichg(ifam(n,i)))

c extra code to account for atoms in IN_CHARGES
       if(numcharges.gt.0)then
        do j=1,numcharges
         if(ifam(n,i).eq.nchgatoms(j))nchg=nchg+1
        enddo
       endif

       if(iabs(ichg(ifam(n,i))).gt.0)name1=ifam(n,i)
       enddo

       if(name1.eq.0)name1=ifam(n,1)


c allow all groups to have npa charges calculated
c      nchg=1

c temp skip to create "charge" files for every group"
c       if(iabs(nsum).gt.0)then

         if(nchg.gt.0)then

c as of 060613, use the group number in the name of the file
         name1=n

         kcharge=kcharge+1
         kgroup(kcharge)=n
         ncount=0
         do i=1,nfam(n)
          if(iabs(ichg(ifam(n,i))).gt.0)then
           ncount=ncount+1
           natchg(kcharge,ncount)=ifam(n,i)
          endif
         enddo
         numchgs(kcharge)=ncount

c find caps as these must be added to get the mean group position
       numbercaps=0
       do j=1,3
        xmean(j)=0.d0
       enddo
       do i=1,nfam(n)
        do m=1,nbondso
         if(mbstore(m).eq.ifam(n,i))then
          if(ngp(nbstore(m)).ne.n.and.multstore(m).gt.0)then
            sum1=0.23d0+radius(mbstore(m))
            sum2=radius(mbstore(m))+radius(nbstore(m))
            fact1=sum1/sum2
            numbercaps=numbercaps+1
            do j=1,3
            xmean(j)=xmean(j)+c(mbstore(m),j)
     .      +fact1*(c(nbstore(m),j)-c(mbstore(m),j))
            xcap(numbercaps,j)=c(mbstore(m),j)
     .      +fact1*(c(nbstore(m),j)-c(mbstore(m),j))
            enddo
            nattach(numbercaps)=mbstore(m)
          endif
         endif
         if(nbstore(m).eq.ifam(n,i))then
          if(ngp(mbstore(m)).ne.n.and.multstore(m).gt.0)then
            sum1=0.23d0+radius(nbstore(m))
            sum2=radius(mbstore(m))+radius(nbstore(m))
            fact1=sum1/sum2
            numbercaps=numbercaps+1
            do j=1,3
            xmean(j)=xmean(j)+c(nbstore(m),j)
     .      +fact1*(c(mbstore(m),j)-c(nbstore(m),j))
            xcap(numbercaps,j)=c(nbstore(m),j)
     .      +fact1*(c(mbstore(m),j)-c(nbstore(m),j))
            enddo
            nattach(numbercaps)=nbstore(m)
          endif
         endif
        enddo
       enddo

       do k=1,3
       do i=1,nfam(n)
        xmean(k)=xmean(k)+c(ifam(n,i),k)
       enddo
       enddo
       do k=1,3
        xmean(k)=xmean(k)/dble(float(nfam(n)+numbercaps))
       enddo

      k=name1
      if(k.le.9)then
       ca='charge.'//I0(k)//'.chg'
       da='charge.'//I0(k)//'.cart'
       ea='charge.'//I0(k)//'.coord'
       fa='charge.'//I0(k)//'.coord.cap'
       pa='charge.'//I0(k)//'.grp'
      endif
      if(k.gt.9.and.k.le.99)then
      k1=k/10
      k2=k-k1*10
      ca1=I0(k1)
      ca2=I0(k2)
      ca='charge.'//ca1//ca2//'.chg'
      da='charge.'//ca1//ca2//'.cart'
      ea='charge.'//ca1//ca2//'.coord'
      fa='charge.'//ca1//ca2//'.coord.cap'
       pa='charge.'//ca1//ca2//'.grp'
      endif
      if(k.gt.99.and.k.le.999)then
      k1=k/100
      k3=k/10 - INT(k1)*10
      k4=k - INT(k1)*100 - INT(k3)*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca='charge.'//ca1//ca2//ca3//'.chg'
      da='charge.'//ca1//ca2//ca3//'.cart'
      ea='charge.'//ca1//ca2//ca3//'.coord'
      fa='charge.'//ca1//ca2//ca3//'.coord.cap'
       pa='charge.'//ca1//ca2//ca3//'.grp'
      endif
      if(k.gt.999.and.k.le.9999)then
      k1=k/1000
      k2=k/100 - INT(k1)*10
      k3=k/10 - INT(k2)*10-int(k1)*100
      k4=k - INT(k1)*1000 - INT(k3)*10-int(k2)*100
      ca1=I0(k1)
      ca2=I0(k2)
      ca3=I0(k3)
      ca4=I0(k4)
      ca='charge.'//ca1//ca2//ca3//ca4//'.chg'
      da='charge.'//ca1//ca2//ca3//ca4//'.cart'
      ea='charge.'//ca1//ca2//ca3//ca4//'.coord'
      fa='charge.'//ca1//ca2//ca3//ca4//'.coord.cap'
       pa='charge.'//ca1//ca2//ca3//ca4//'.grp'
      endif

      panames(kcharge)=pa
      nend=index(ea," ")-1
      ga=ea(1:nend)//'_atoms'

c calculate the charge and miltiplicity
      nchargetot=0
      numtot=0
      do i=1,nfam(n)
       numtot=numtot+numa(ifam(n,i))
       nchargetot=nchargetot+ichg(ifam(n,i))
      enddo
      numtot=numtot+numbercaps-nchargetot

      numatoms(kcharge)=nfam(n)+numbercaps
      numelects(kcharge)=numtot

      numtot2=numtot/2
      ndg=numtot-2*numtot2
      if(ndg.eq.0)then
       ndegen=1
      else
       ndegen=2
      endif

c make a blank file for group charges
c which will be filled in other scripts/programs
      open(unit=66,file=pa,status='unknown')
      close(unit=66)

      open(unit=66,file=ca,status='unknown')
      write(66,100)(xmean(k),k=1,3),dble(float(nsum))
      close(unit=66)
      open(unit=67,file=da,status='unknown')
      write(67,*)'  The distributed Cartesian multipoles'
      write(67,*) nfam(n)+numbercaps
      open(unit=68,file=ea,status='unknown')
      open(unit=70,file=ga,status='unknown')
c     write(68,*)nsum," 1"
      write(68,*)nsum,ndegen
      do i=1,nfam(n)
      write(68,666)atoms(ifam(n,i)),(c(ifam(n,i),k1),k1=1,3)
      write(70,*)ifam(n,i)
      enddo
      do i=1,numbercaps
       write(68,666)'H ',(xcap(i,k1),k1=1,3)
       write(70,*)0
      enddo
      close(unit=68)
      close(unit=70)
      open(unit=69,file=fa,status='unknown')
      write(69,*)nfam(n),numbercaps,n
      do i=1,numbercaps
       nsite=0
       do i1=1,nfam(n)
        if(nattach(i).eq.ifam(n,i1))nsite=i1
       enddo
       if(nsite.eq.0)then
        write(6,*)' nsite = 0 in families.f'
        stop
       endif
      write(69,333)(xcap(i,k1),k1=1,3),nsite
      enddo
333   format(3f14.7,I10)
      close(unit=69)

      write(83,*) nfam(n)+numbercaps
      write(83,*)nfam(n)
      do i=1,nfam(n)
       write(67,*)atoms(ifam(n,i))
       do k1=1,3
        write(67,*)c(ifam(n,i),k1)*bohr
       enddo
        write(90,333)(c(ifam(n,i),k1)*bohr,k1=1,3)
      write(83,666)atoms(ifam(n,i)),(c(ifam(n,i),k1),k1=1,3)
666   format(a2,3f13.6,i10)
       write(67,*)'  Charge'
       if(i.eq.1)then
        write(67,*) dble(float(nsum))
       else
        write(67,*) 0.d0
       endif 
       write(67,*)'  Dipole'
       do k1=1,3
        write(67,*) 0.d0
       enddo
       write(67,*)'  Quodrupole'
       do k1=1,3
       do k2=1,k1
        write(67,*) 0.d0
       enddo
       enddo
       write(67,*)'  Octapole'
       do k1=1,3
       do k2=1,k1
       do k3=1,k2
        write(67,*) 0.d0
       enddo
       enddo
       enddo
       write(67,*)'  Hexadecapole'
       do k1=1,3
       do k2=1,k1
       do k3=1,k2
       do k4=1,k3
        write(67,*) 0.d0
       enddo
       enddo
       enddo
       enddo
      enddo
      do i=1,numbercaps
       write(67,*)'H '
       do k1=1,3
        write(67,*)xcap(i,k1)*bohr 
       enddo
        write(90,333)(xcap(i,k1)*bohr,k1=1,3)
       write(83,666)'H ',(xcap(i,k1),k1=1,3)
       write(67,*)'  Charge'
        write(67,*) 0.d0
       write(67,*)'  Dipole'
       do k1=1,3
        write(67,*) 0.d0
       enddo
       write(67,*)'  Quodrupole'
       do k1=1,3
       do k2=1,k1
        write(67,*) 0.d0
       enddo
       enddo
       write(67,*)'  Octapole'
       do k1=1,3
       do k2=1,k1
       do k3=1,k2
        write(67,*) 0.d0
       enddo
       enddo
       enddo
       write(67,*)'  Hexadecapole'
       do k1=1,3
       do k2=1,k1
       do k3=1,k2
       do k4=1,k3
        write(67,*) 0.d0
       enddo
       enddo
       enddo
       enddo
      enddo

      close(unit=67)
      write(83,*)


c end the nchg if
      endif

c close the natom loop
      enddo

      open(unit=30,file='OUT_ELECTRONS',status='unknown')
      write(30,*)' The number of charge calculations is'
      write(30,*)kcharge
      if(kcharge.gt.0)then
       write(30,*)' The number of atoms and electrons for charge calcs'
       do j=1,kcharge
        write(30,*)numatoms(j),numelects(j),panames(j)
       enddo
      endif



c output the identities of the groups which contain formal charges
c and the associated numbers of the charged atoms
      open(unit=66,file='OUT_CHARGEDGROUPS',status='unknown')
      write(66,*)' The number of groups with formal charges is'
      write(66,*)kcharge
      write(66,*)' The group numbers,# of charges, and atom numbers are'
      do j=1,kcharge
       write(66,*)kgroup(j)
       write(66,*)numchgs(j)
       write(66,*)(natchg(j,k),k=1,numchgs(j))
       write(66,*)
      enddo
      close(unit=66)

c output the charged groups for viewing
      open(unit=66,file='VIEWCHARGEDGROUPS',status='unknown')
      write(66,*)natomall
      write(66,*)
      do j=1,natomall
      write(66,666)atoms(j),(c(j,k),k=1,3)
      enddo
      do j=1,kcharge
      if(formal(kgroup(j)).eq.1)then
       write(66,*)nfam(kgroup(j))
       write(66,667)(ifam(kgroup(j),m),m=1,nfam(kgroup(j)))
       do m=1,nfam(kgroup(j))
        write(66,666)atoms(ifam(kgroup(j),m)),(c(ifam(kgroup(j),m),k),k=1,3)
       enddo
      endif
      enddo
667   format(21I6)

c end the Level if
      endif

100   format(4f13.6) 

      deallocate(ngp)
      deallocate(kgroup)
      deallocate(numchgs)
      deallocate(natchg)
      if(numcharges.gt.0)deallocate(nchgatoms)

      return
      end

