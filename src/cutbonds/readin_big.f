      subroutine readin_big

      use fractdata
      implicit double precision(a-h,o-z)

      character*1 marker

      integer, allocatable  :: iproj(:),idouble(:,:),nblank(:),ndoub(:)

      integer, allocatable  :: ifamtemp(:),kbond(:)

      integer, allocatable  :: kgroup(:),jbd(:)

      integer alab2anum

      character*20 ca,ca1

      open(unit=1,file='name.xyz',status='old')

      read(1,*)natomall
      maxdouble=natomall/2 +1

c allocate some arrays
      allocate(atoms(natomall))
      allocate(c(natomall,3))
      allocate(numa(natomall))
      allocate(radius(natomall))
      allocate(idouble(natomall,maxdouble))
      allocate(nblank(natomall))
      allocate(ndoub(natomall))
      allocate(ichg(natomall))
      allocate(iproj(natomall))
      allocate(ifamtemp(natomall))

      allocate(kbond(natomall))

      read(1,90)icomm
90     format(a80)
      do n=1,natomall
c     read(1,91)atoms(n),(c(n,k),k=1,3)
      read(1,*)atoms(n),(c(n,k),k=1,3)
      enddo
91    format(a2,3f13.6)

      close(unit=1)

c use functions to get the atomic number and
c the covalent radius
      do n=1,natomall
       numa(n)=alab2anum(atoms(n))
c check for error
       if(numa(n).eq.0)then
        write(6,*)' An elemental symbol for atom ',n
        write(6,*)' was not recognised'
        stop
       endif
       radius(n)=anum2cov(numa(n))
      enddo

      do n=1,natomall
      nblank(n)=0
      ndoub(n)=0
      iproj(n)=0
      do m=1,maxdouble
      idouble(n,m)=0
      enddo
      enddo
 

      do n=1,natomall
      if(atoms(n).eq.' H')atoms(n)='H'
      enddo

      close(unit=1)

      open(unit=2,file='name.mol',status='old')

      lin=0
1     continue 
      read(2,92,end=2)marker
      lin=lin+1
      go to 1
2     continue
92    format(a1)
      nbonds=lin-4-natomall-1
      close(unit=2)
      open(unit=2,file='name.mol',status='old')
      read(2,303)ititle
303   format(a20)
      write(6,303)ititle
c read in any charges
      do i=1,3
      read(2,90)icomm
      enddo
      do i=1,natomall
      read(2,8001)x1,y1,z1,lunk,ichg(i)
      enddo
8001  format(3f10.4,1x,a2,i4)

c allocate some arrays
      allocate(mb(nbonds))
      allocate(nb(nbonds))
      allocate(mult(nbonds))


c  now read in the bonded atoms
      write(6,*)'bonds'
      do m=1,nbonds
      read(2,*)mb(m),nb(m),mult(m)
      if(mb(m).gt.nb(m))then
      ntem=nb(m)
      nb(m)=mb(m)
      mb(m)=ntem
      endif
      write(6,*)mb(m),nb(m),mult(m)
      enddo
      write(6,*)
      close(unit=2)

      call checkbonds(nbonds,mb,nb,mult)

c allow for Hbonds read in
       nhbonds=0
      open(unit=2,file='IN_HBONDS',status='unknown')
      read(2,90,end=777)icomm
      read(2,*)nhbonds
777   continue

c new allocations
      nnew=nbonds+nhbonds
      allocate(mbstore(nnew))
      allocate(nbstore(nnew))
      allocate(multstore(nnew))

c  make a copy of the original bonding
      do m=1,nbonds
      mbstore(m)=mb(m)
      nbstore(m)=nb(m)
      multstore(m)=mult(m)
      enddo
      if(nhbonds.gt.0)then
       read(2,90)icomm
       do m=nbonds+1,nbonds+nhbonds
c      read(2,*)mb(m),nb(m),mult(m)
       read(2,*)mbstore(m),nbstore(m),multstore(m)
       enddo
       nbonds=nbonds+nhbonds
      endif
      close(unit=2)
125   continue

      nbondso=nbonds

c reset the mb,nb,mult arrays
      deallocate(mb)
      deallocate(nb)
      deallocate(mult)
      allocate(mb(nbonds))
      allocate(nb(nbonds))
      allocate(mult(nbonds))
      do m=1,nbonds
       mb(m)=mbstore(m)
       nb(m)=nbstore(m)
       mult(m)=multstore(m)
      enddo

c look for any atoms that are not bonded to anything
      kbond=0
      do m=1,nbonds
       kbond(mb(m))=1
       kbond(nb(m))=1
      enddo
       

c now we set up the arrays, assuming that multiple bonds
c are never broken, and that XH bonds are never broken

c we make XH one atom and make A=B one atom

111   continue

      do m=1,nbonds
      if(mult(m).le.1)go to 11
      if(mb(m).eq.nb(m))then
      nb(m)=0
      mb(m)=0
      mult(m)=0
      go to 11
      endif
      ndoub(mb(m))=ndoub(mb(m))+1
      if(ndoub(mb(m)).gt.maxdouble)then
      write(6,*)'ndoub > maxdouble for bond ',m,' and atom ',mb(m)
      write(6,*)'ndoub = ',ndoub(mb(m))
      stop
      endif
      idouble(mb(m),ndoub(mb(m)))=nb(m)
      nblank(nb(m))=mb(m)
c patch 060309
      do j=1,natomall
       if(nblank(j).eq.nb(m))then
         nblank(j)=mb(m)
         ndoub(mb(m))=ndoub(mb(m))+1
         ndoub(nb(m))=ndoub(nb(m))-1
         idouble(mb(m),ndoub(mb(m)))=j
       endif
      enddo
c end of patch
c replace nb atom in multiple bond by mb, everywhere
      do k=1,nbonds
      if(k.eq.m)go to 771
      if(mb(k).eq.nb(m))mb(k)=mb(m)
      if(nb(k).eq.nb(m))then
       nb(k)=mb(k)
       mb(k)=mb(m)
      endif
771   continue
      enddo

c set the multiple bond void
      nb(m)=0
      mb(m)=0
      mult(m)=0
11    continue
      enddo
      
      go to 3333

c must check to see if two groups are bonded to each other twice
c this can happen if one set of atoms is conected by double bonds
c to form one group, and the same for another set. Then if different
c memebrs of one set are single bonded to different atoms
c of the other set, it appears that the two sets are bonded twice.
c They must be made one group or adjacent caps will result


      match=0
      do n=1,nbonds-1
      if(mult(n).eq.0)go to 1111
       do k=n+1,nbonds
        if((mb(k).eq.mb(n)).and.(nb(k).eq.nb(n)))then
        mult(n)=2
        mult(k)=0
        mb(k)=0
        nb(k)=0
        match=1
        endif
        if((mb(k).eq.nb(n)).and.(nb(k).eq.mb(n)))then
        mult(n)=2
        mult(k)=0
        mb(k)=0
        nb(k)=0
        match=1
        endif
       enddo
1111   enddo


       if(match.eq.1)go to 111

3333    continue
     
c    temp write
      write(6,*)' bonds after Hbonds and doubles'
      do n=1,nbonds
      write(6,*)mb(n),nb(n),mult(n)
      enddo


c  set up new atom or group numbers, eliminating H atoms 
      natom=0
      do n=1,natomall
c  don't count H atoms
      if(atoms(n).eq.'H')go to 10
c  don't count eliminated atoms
      if(nblank(n).gt.0)go to 10
      natom=natom+1
      ifamtemp(natom)=n
c     ifam(natom,1)=n
      iproj(n)=natom
10    continue
      enddo

      do n=1,natomall
      if(nblank(n).gt.0)then
      iproj(n)=iproj(nblank(n))
      endif
      enddo

c allow for single atoms and H2 molecules as groups
c 14/11/17 disallow this
c      do n=1,nbonds
c       if(mb(n)*nb(n).eq.0)go to 891
c       if(atoms(mb(n)).eq.'H'.and.atoms(nb(n)).eq.'H')then
c        natom=natom+1
c        ifamtemp(natom)=mb(n)
c        iproj(mb(n))=natom
c       endif
c891      enddo

c      do n=1,natomall
c       if(kbond(n).eq.0.and.atoms(n).eq.'H')then
c        natom=natom+1
c        ifamtemp(natom)=n
c        iproj(n)=natom
c       endif
c      enddo

c     write(6,*)' natom = ',natom
c     stop

c natom is now known, so we call allocate some arrays
      allocate(nfam(natom))
      allocate(ifam(natom,maxfamily))
      nfragm=60*natom
c     nsmall2=nsmall/2
      allocate(numat(nfragm))
      allocate(natstore(nfragm,nsmall))
      allocate(ibond(nfragm,nsmall,nsmall))
      allocate(itype(nfragm,nsmall))
      allocate(isign(nfragm))
      allocate(nstop(nfragm))
      allocate(map(nsmall,nsmall))
      allocate(kf(nsmall))

      do n=1,natom
      nfam(n)=0
      do k=1,maxfamily
      ifam(n,k)=0
      enddo
      ifam(n,1)=ifamtemp(n)
      enddo

c  put H atoms in the family of bonded atom
      do n=1,natom
      ic=1
      do m=1,nbonds
      if(mb(m).eq.0)go to 887
      if(mb(m).eq.ifam(n,1).and.atoms(nb(m)).eq.'H')then
      ic=ic+1
      if(ic.gt.maxfamily)then
      write(6,*)' nfam > maxfamily for n = ',n
      stop
      endif
      ifam(n,ic)=nb(m)
      iproj(nb(m))=ifam(n,1)
      mb(m)=0
      nb(m)=0
      mult(m)=0
      endif
      if(nb(m).eq.ifam(n,1).and.atoms(mb(m)).eq.'H')then
      ic=ic+1
      if(ic.gt.maxfamily)then
      write(6,*)' nfam > maxfamily for n = ',n
      stop
      endif
      ifam(n,ic)=mb(m)
      iproj(mb(m))=ifam(n,1)
      mb(m)=0
      nb(m)=0
      mult(m)=0
      endif
887   continue
      enddo
      nfam(n)=ic
      enddo

      do n=1,natomall
      if(ndoub(n).eq.0)go to 888
      do j=1,ndoub(n)
      nfam(iproj(n))=nfam(iproj(n))+1
      if(nfam(iproj(n)).gt.maxfamily)then
      write(6,*)' nfam > maxfamily for double bonds at n = ',n
      stop
      endif

      ifam(iproj(n),nfam(iproj(n)))=idouble(n,j)
889   continue
      enddo
888   continue
      enddo

c temp check
c     write(6,*)' iproj'
c     do n=1,natomall
c      write(6,*)n,iproj(n)
c     enddo

c allow for H atoms or H2 molecules as groups
c       do n=1,natomall
c       if(iproj(n).ne.0.or.atoms(n).ne.'H')go to 890
c       natom=natom+1
c       iproj(n)=natom
c       nfam(natom)=1
c       ifam(natom,1)=n
c       do m=1,nbonds
c        if(mb(m).eq.n)then
c         nfam(natom)=nfam(natom)+1
c         ifam(natom,nfam(natom))=nb(m)
c         iproj(nb(m))=natom
c        endif
c        if(nb(m).eq.n)then
c         nfam(natom)=nfam(natom)+1
c         ifam(natom,nfam(natom))=mb(m)
c         iproj(mb(m))=natom
c        endif
c       enddo
c890   enddo

c     write(6,*)' iproj'
c     do n=1,natomall
c      write(6,*)n,iproj(n)
c     enddo
c     stop

c check for multiples of the same bond
      do m=1,nbonds-1
      do n=m+1,nbonds
      if(mb(n).eq.mb(m).and.nb(n).eq.nb(m))then
      mb(n)=0
      nb(n)=0
      endif
      if(mb(n).eq.nb(m).and.nb(n).eq.mb(m))then
      mb(n)=0
      nb(n)=0
      endif
      enddo
      enddo

c checking
      write(6,*)' bonds with group numbers'
      do m=1,nbonds
      if(mb(m).gt.0.and.iproj(mb(m)).ne.iproj(nb(m)))then
      write(6,*)iproj(mb(m)),iproj(nb(m)),mult(m)
      endif
      enddo
      write(6,*)

c store the Hbond data separately for use in rings (three, four etc)
c     numhbonds=0
c     if(nhbonds.eq.0)go to 778
c     allocate(mh(nhbonds))
c     allocate(nh(nhbonds))
c     do m=1,nbonds
c      if(mult(m).eq.-1)then
c       numhbonds=numhbonds+1
c       mh(numhbonds)=iproj(mb(m))
c       nh(numhbonds)=iproj(nb(m))
c      endif
c     enddo
c778   continue


      ib=0
      do m=1,nbonds
      if(mb(m).eq.0)go to 12
c 221209 remove bonds which bond a composite atom to itself
      if(iproj(mb(m)).eq.iproj(nb(m)))go to 12
      ib=ib+1
      mb(ib)=iproj(mb(m))
      nb(ib)=iproj(nb(m))
      mult(ib)=mult(m)
12    continue
      enddo
      write(6,*)'nbonds = ',ib
      write(6,*)
      write(6,*)' The number of atoms is ',natom
      write(6,*)'The atoms in each group are:'
      do n=1,natom
      write(6,*)nfam(n),(ifam(n,k),k=1,nfam(n))
      enddo
c     stop
c output the families
      call families

      nsum=0
      do n=1,natom
      nsum=nsum+nfam(n)
      enddo
      if(iabs(nsum-natomall).gt.0)then
      write(6,*)
      write(6,*)' The sum of the famililies for each group'
      write(6,*)' does not equal the total number of atoms.'
      write(6,*)' This may be bacause a H atom was not'
      write(6,*)' deemed to be bonded to anything.'
      write(6,*)
      stop
      endif
      write(6,*)


      nbonds=ib

c 271217
       if(Level.gt.0)then
       write(6,*)' nbonds = ',nbonds
       allocate(jbd(natom))
       write(6,*)' nbonds = ',nbonds
       do n=1,natom
        jbd(n)=0
       enddo
       do n=1,nbonds
        jbd(mb(n))=jbd(mb(n))+1
        jbd(nb(n))=jbd(nb(n))+1
       enddo
c      write(6,*) 'jbd'
c      do n=1,natom
c       write(6,*)n,jbd(n)
c      enddo

       allocate(mbn(nbonds))
       allocate(nbn(nbonds))
       allocate(multn(nbonds))
       ic=0
       do n=1,nbonds
        if((jbd(mb(n)).gt.2).and.(jbd(nb(n)).gt.2))then
         ic=ic+1
         mbn(ic)=mb(n)
         nbn(ic)=nb(n)
         multn(ic)=mult(n)
         jbd(mb(n))=jbd(mb(n))-1
         jbd(nb(n))=jbd(nb(n))-1
         mb(n)=0
         nb(n)=0
        endif
       enddo
       nbondsextra=ic

       ic=0
       do n=1,nbonds
        if(mb(n).gt.0)then
         ic=ic+1
         mb(ic)=mb(n)
         nb(ic)=nb(n)
         mult(ic)=mult(n)
        endif
       enddo
      nbonds=ic

       write(6,*)' nbonds = ',nbonds
       write(6,*)' nbondsextra = ',nbondsextra
c      do n=1,natom
c       jbd(n)=0
c      enddo
c      do n=1,nbonds
c       jbd(mb(n))=jbd(mb(n))+1
c       jbd(nb(n))=jbd(nb(n))+1
c      enddo
c      write(6,*) 'jbd'
c      do n=1,natom
c       write(6,*)n,jbd(n)
c      enddo
       endif

c store numbers for groups connected by H bonds
c     nhstore=0
c     do n=1,nbonds
c      if(mult(n).eq.-1)nhstore=nhstore+1
c     enddo

c     allocate(minusbonds(nhstore,2))
c     ic=0
c     do n=1,nbonds
c      if(mult(n).eq.-1)then
c       ic=ic+1
c       minusbonds(ic,1)=mb(n)
c       minusbonds(ic,2)=nb(n)
c      endif
c     enddo

 
      write(6,*)
c  calculate the itype array
      write(6,*)
      write(6,*)
c check to see if all bonds between groups have mult = -1
c if so then "nocapsatall =1" in IN_CAPS
      open(unit=17,file='IN_CAPS',status='unknown')
      write(17,*)'Enter 0 for caps or 1 for no caps'
      nocapsatall=1
      do m=1,nbonds
      write(6,*)mb(m),nb(m),mult(m)
      if(mult(m).gt.0)nocapsatall=0
      enddo
      write(17,*)nocapsatall
      close(unit=17)
c     stop

      allocate(itf(natom,ncomp))
      allocate(ibf(natom,nsmall,ncomp))
      allocate(ngpf(natom,2,ncomp))
      allocate(numgroups(ncomp))

      do n=1,natom
      itf(n,1)=0
      do m=1,nbonds
      if(mb(m).eq.n.or.nb(m).eq.n)itf(n,1)=itf(n,1)+1
      enddo
      itf(n,1)=itf(n,1)-1
      enddo
      write(6,*)' ic and itype'
c  calculate the ibond array
      do n=1,natom
      ic=0
      do m=1,nbonds
      if(mb(m).eq.n)then
      ic=ic+1
      if(ic.gt.nsmall)then
      write(6,*)' ibond overflow'
      stop
      endif
      ibf(n,ic,1)=nb(m)*mult(m)
      endif
      if(nb(m).eq.n)then
      ic=ic+1
      if(ic.gt.nsmall)then
      write(6,*)' ibond overflow'
      stop
      endif
      ibf(n,ic,1)=mb(m)*mult(m)
      endif
      enddo
      write(6,*)n,ic,itf(n,1)+1
      enddo

      write(6,*)'Atom,itype,connected atoms'
      do n=1,natom
      write(6,*)n,itf(n,1),(ibf(n,k,1),k=1,itf(n,1)+1)
      enddo
      write(6,*)

c output the group connectivity
      open(unit=22,file='OUT_GROUPCONNECTIVITY',status='unknown')
      write(22,*)natom
      do n=1,natom
       write(22,*)n,itf(n,1)+1
       write(22,*)(ibf(n,ic,1),ic=1,itf(n,1)+1)
      enddo
      close(unit=22)

      if(Level.eq.0)then
      do n=1,natom
       numat(n)=1
       natstore(n,1)=n
       isign(n)=1
       nstop(n)=1
       itype(n,1)=-1
      enddo
      write(6,*)natom,' groups, one per fragment at Level 0'
      nf=natom
      call writefrags
      call writecom
      call writerawdata

c new code 230417 to output groups that provide
c embedded charges for Level0 fragments
      open(unit=20,file='OUT_CHARGEDGROUPS',status='old')
      read(20,*)
      read(20,*)kcharge
      if(kcharge.le.1)go to 2111
      allocate(kgroup(kcharge))
      read(20,*)
      do k=1,kcharge
       read(20,*)kgroup(k)
       read(20,*)
       read(20,*)
       read(20,*)
      enddo
      do k=1,kcharge
       call filelabel(kgroup(k),ca1)
       n1=index(ca1,' ')-1
       ca='chL0.'//ca1(1:n1)
       open(unit=21,file=ca,status='unknown')
       do m=1,kcharge
        if(m.ne.k)write(21,*)kgroup(m)
       enddo
       close(unit=21)
      enddo
2111  close(unit=20)
c end new code 220417

       stop
      endif

c  writeout the skeleton
c      open(unit=10,file='OUT_SKELETON',status='unknown')
c      write(10,*)natom
c      write(10,*)
c      do n=1,natom
c      m=ifam(n,1)
c      write(10,223)atoms(m),(c(m,k),k=1,3)
c      enddo
c      close(unit=10)
223   format(a2,2x,3f8.3)

c     call neigh5

      deallocate(iproj)
      deallocate(idouble)
      deallocate(nblank)
      deallocate(ndoub)
      deallocate(ifamtemp)

      return
      end


      real*8 function anum2cov(na)


c in Angstrom, from Mick's atomic.dat, other values converted from Chris
c Mick's atomic.dat misses 85-87 incl and from 95 onward

      integer :: na
      real*8 :: covalent_radii(110)

      data covalent_radii/  
     .  0.23, 1.22, 0.68, 0.35, 0.83, 0.68, 0.68, 0.68, 0.64, 1.6,  
     .  0.97, 1.1,  1.35, 1.2,  1.05, 1.02, 0.99, 1.92, 1.33, 0.99, 
     .  1.44, 1.47, 1.33, 0.67, 1.35, 1.34, 1.33, 1.5,  1.52, 1.45, 
     .  1.22, 1.17, 1.21, 1.22, 1.21, 1.98, 1.47, 1.12, 1.78, 1.56, 
     .  1.48, 1.47, 1.35, 1.4,  1.45, 1.5,  1.59, 1.69, 1.63, 1.46, 
     .  1.46, 1.47, 1.4,  2.18, 1.67, 1.34, 1.87, 1.83, 1.82, 1.81, 
     .  1.8,  1.8,  1.99, 1.79, 1.76, 1.75, 1.74, 1.73, 1.72, 1.94, 
     .  1.72, 1.57, 1.43, 1.37, 1.35, 1.37, 1.32, 1.5,  1.5,  1.7,  
     .  1.55, 1.54, 1.54, 1.68, 1.21, 1.50, 1.50, 1.9,  1.88, 1.79, 
     .  1.61, 1.58, 1.55, 1.53, 1.51, 0.99, 1.54, 1.83, 1.50, 1.50, 
     .  1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50/


      anum2cov=covalent_radii(na)

      return

      end function anum2cov

      integer function alab2anum(alab)

c given the atom label this returns the atomic number

      character*2 :: alab,atomic_labels(110)

      integer n

      data atomic_labels/'H','He','Li','Be','B','C','N','O','F','Ne',  
     .'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
     .'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb',
     . 'Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     . 'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',
     . 'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os',
     . 'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     . 'Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',
     . 'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds'/

      do n=1,110
       if(atomic_labels(n).eq.alab)then
        alab2anum=n
        return
       endif
      enddo
      alab2anum=0

      return

      end function alab2anum
