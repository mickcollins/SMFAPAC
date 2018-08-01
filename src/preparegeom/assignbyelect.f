      program assigncharges
      use AtomicData
      implicit double precision(a-h,o-z)
      parameter(maxbond=20)

      character*2,  allocatable  :: lab(:)
      real*8,       allocatable  :: c(:,:),rad(:),elect(:),weight(:)
      integer,      allocatable  :: num(:),nval(:),bonds(:,:),mult(:,:)
      integer,      allocatable  :: nbonds(:),ival(:),ncharge(:)
      integer,      allocatable  :: nshort(:)
      integer,      allocatable  :: ichoice(:),indx(:)
      integer,      allocatable  :: like(:),nused(:)
      integer,      allocatable  :: special(:),natspec(:),ligch(:)
      integer,      allocatable  :: nchspec(:),ivalspec(:)
      integer,      allocatable  :: totbonds(:)
      integer,      allocatable  :: sec(:,:),jmat(:,:),neig(:)
      integer,      allocatable  :: mbex(:),nbex(:)
      integer,      allocatable  :: dmbd(:),dnbd(:)

      real*8,       allocatable  :: shortest(:),dist(:,:)

      dimension neigh(maxbond),ind(maxbond),nchoice(maxbond)
      dimension eneigh(maxbond),dtemp(maxbond),ene(maxbond)
      dimension choice(maxbond),x(3)

      dimension bd(110,110)

      character*1, preserve(2)

c a program to assign charges to atoms in a molecule
c all metals have the charge input by the user, and
c we don't count any bonds to metals for the purpose of 
c ascertaining how many bonds another atom has

      open(unit=10,file='WARNINGS',status='unknown')
c we output any unusual bonding or charges

c open IN_EXTRABONDS, if it exists, and read it
      numextra=0
       allocate(mbex(100))
       allocate(nbex(100))
      open(unit=1,file='IN_EXTRABONDS',status='unknown')
      read(1,*,end=2)
      read(1,*)numextra
      if(numextra.gt.0)then
       read(1,*)
       do n=1,numextra
        read(1,*)mbex(n),nbex(n)
       enddo
      endif
2     continue
      close(unit=1)

c read in any bonds to be made into double bonds
      ndouble=0
      open(unit=1,file='IN_DOUBLE',status='unknown')
      read(1,*,end=4)
      read(1,*)ndouble
      if(ndouble.eq.0)go to 4
      allocate(dmbd(ndouble))
      allocate(dnbd(ndouble))
      read(1,*)
      do n=1,ndouble
       read(1,*)dmbd(n),dnbd(n)
      enddo
4     close(unit=1)

c get the info from IN_PRESERVE_BONDS

c      open(unit=1,file='IN_PRESERVE_BONDS',status='unknown')
c      read(1,*,end=5)
c      do n=1,4
c      read(1,*)preserve(n)
c      enddo
c5     close(unit=1)


c assume standard xyz format for the molecule
      read(5,*)natom
      read(5,*)

      allocate(lab(natom))
      allocate(c(natom,3))
      allocate(rad(natom))
      allocate(elect(natom))
      allocate(weight(natom))
      allocate(num(natom))
      allocate(nbonds(natom))
      allocate(totbonds(natom))
      allocate(bonds(natom,maxbond))
      allocate(mult(natom,maxbond))
      allocate(dist(natom,maxbond))
      allocate(like(maxbond))
      allocate(nval(natom))
      allocate(ival(natom))
      allocate(nused(natom))
      allocate(ncharge(natom))
      allocate(ichoice(natom))
      allocate(indx(natom))
      allocate(shortest(natom))
      allocate(nshort(natom))

      allocate(special(natom))
      allocate(natspec(natom))
      allocate(ligch(natom))
      allocate(nchspec(natom))
      allocate(ivalspec(natom))
      allocate(sec(natom,2))

      allocate(jmat(natom,natom))
      allocate(neig(natom))


      do n=1,natom
c assume standard xyz format
      read(5,*)lab(n),(c(n,k),k=1,3)
c check the element symbols are correct
      call elements(lab(n))
      enddo

c echo the coordinates in standard xyz format
      open(unit=2,file='molecule.xyz',status='unknown')
      write(2,*)natom
      write(2,*)
      do n=1,natom
       write(2,400)lab(n),(c(n,k),k=1,3)
      enddo
400   format(a2,3f10.4)
      close(unit=2)

c assign atomic numbers

      do n=1,natom
       num(n)=aLabToANum(lab(n))
      enddo

c assign covalent radii
      do n=1,natom
       rad(n)=aRad(num(n))
      enddo

c assign electronegativities
      do n=1,natom
       elect(n)=anum2elneg(num(n))
      enddo
 
c initialize the array of bond lenth limits for double bonds
      call bondlimits(bd)


c assign the number of valence electrons
      do n=1,natom
       nval(n)=0
       if(num(n).le.2)nval(n)=num(n)
       if(num(n).gt.2.and.num(n).le.10)nval(n)=num(n)-2
       if(num(n).gt.10.and.num(n).le.18)nval(n)=num(n)-10
       if(num(n).gt.18.and.num(n).le.20)nval(n)=num(n)-18
       if(num(n).gt.30.and.num(n).le.36)nval(n)=num(n)-28
       if(num(n).gt.36.and.num(n).le.38)nval(n)=num(n)-36
       if(num(n).gt.48.and.num(n).le.54)nval(n)=num(n)-46
       if(num(n).gt.54.and.num(n).le.56)nval(n)=num(n)-54
       if(num(n).gt.80.and.num(n).le.86)nval(n)=num(n)-78
      enddo

c the valence of all types of metals is unassigned

c  look for specified charges and valency
      write(66,*)' input specified charges'
      open(unit=1,file='IN_SPECIFIED_CHARGES',status='unknown')
      ist1=0
      read(1,*,end=200)
c an existing IN_SPECIFIED_CHARGES sets nflag=1
      nflag=1
      ist1=1
      nspec=0
      read(1,*,end=200)
      read(1,*)
201   continue
      read(1,*,end=200)n1,n2
      write(66,*)n1,n2
      nspec=nspec+1
      natspec(nspec)=n1
      nchspec(nspec)=n2
      go to 201
200   continue
      if(nflag.eq.1)go to 202
c an empty IN_SPECIFIED_CHARGES sets nflag=0
      nflag=0
      close(unit=1)
      if(ist1.eq.0)nspec=0
202   continue
c echo this input
      if(nspec.gt.0)then
      write(66,*)' The following specific charges were input'
      write(66,*)' Atom number    Charge     '
      do n=1,nspec
      write(66,300)natspec(n),nchspec(n),
     .           '      (',lab(natspec(n)),')' 
300   format(4x,i6,7x,i3,8x,a7,a2,a1)
      enddo
      write(66,*)
      endif

      do n=1,natom
      special(n)=0
      enddo

c temporary? assignment of charges to metals
c we have to decide how a decision is made on
c metal charges
      if(nspec.eq.0.and.nflag.eq.0)then
       do n=1,natom
        if(nval(n).eq.0)then
         nspec=nspec+1
         natspec(nspec)=n
         nchspec(nspec)=2
        endif
       enddo
      endif

      if(nspec.gt.0)then
       do n=1,nspec
        special(natspec(n))=1
        ncharge(natspec(n))=nchspec(n)
c for all "special" atoms, we put the electronegativity to be large
        elect(natspec(n))=99.d0
       enddo
      endif


c check all special cases are included
      do n=1,natom
       if(special(n).eq.0)then
       if(nval(n).eq.0)then
      write(66,*)' Atom ',n,', element ',lab(n),'has no assigned charge'
      write(66,*)' The user must assign a charge for this atom'
      write(66,*)' in the specifying metals option in the input'
c     write(66,*)' in the file IN_SPECIFIED_CHARGES'
c     write(66,*)'  or via the input setup'
      stop
       endif
       endif
      enddo


c  assign the  valency
      do n=1,natom
c       if(special(n).eq.1)go to 13  14/11/17
       if(nval(n).eq.0)go to 13

       if(nval(n).le.4)ival(n)=nval(n)
       if(nval(n).gt.4)ival(n)=8-nval(n)
       if(num(n).eq.1)ival(n)=1
       if(num(n).eq.2)ival(n)=2
13    enddo

      do n=1,natom
       nbonds(n)=0
       do m=1,maxbond
        dist(n,m)=1.d6
        bonds(n,m)=0
       enddo
      enddo

c hydrogens can only be attached to one atom
c so treat hydrogens separately first
c connecting each to its closest atom
c if the closest atom is outside tol range, then add this bond as an
c "EXTRABOND"
c also account for HH bonds
      numberofbonds=0
      do n=1,natom
       if(lab(n).ne.'H')go to 561
       distmax=1.d6
       distH=1.d6
       maxH=0
       mmax=0
       do m=1,natom
        if(m.eq.n)go to 562
        dis=0.d0
        do k=1,3
         dis=dis+(c(n,k)-c(m,k))**2
        enddo
        dis=sqrt(dis)
        if(lab(m).eq.'H'.and.dis.lt.distH)then
         distH=dis
         maxH=m
        endif
        if(lab(m).ne.'H'.and.dis.lt.distmax)then
         distmax=dis
         mmax=m
        endif
562    enddo
       if(maxH.eq.0)go to 565

       tol=rad(n)+rad(maxH)+0.4d0
       if(distH.le.tol)then
       write(10,*)' A very short (unphysical?) H--H distance was found'
       write(10,*)' between atoms ',n,' and ',maxH
c       write(6,*)' HH bond ',n
c       write(6,*)maxH,distH
c       write(6,*)mmax,distmax
c avoid repetition of HH bonds
        if(nbonds(n).gt.0)then
         do k=1,nbonds(n)
          if(bonds(n,k).eq.maxH)go to 565
         enddo
        endif
        nbonds(n)=nbonds(n)+1
        nbonds(maxH)=nbonds(maxH)+1
        bonds(n,nbonds(n))=maxH
        bonds(maxH,nbonds(maxH))=n
        numberofbonds=numberofbonds+1
        dist(n,nbonds(n))=2.d0
        dist(maxH,nbonds(maxH))=2.d0
       endif

565    continue
       if(mmax.eq.0)go to 566

       tol=rad(n)+rad(mmax)+0.4d0
c      if(n.eq.15)write(6,*)tol,distmax,mmax
       if(distmax.le.tol)then
        nbonds(n)=nbonds(n)+1
        nbonds(mmax)=nbonds(mmax)+1
        bonds(n,nbonds(n))=mmax
        bonds(mmax,nbonds(mmax))=n
        numberofbonds=numberofbonds+1
        dist(n,nbonds(n))=2.d0
        dist(mmax,nbonds(mmax))=2.d0
        if(distmax.le.0.6d0*tol)then
         write(10,*)' A very short (unphysical?) atom-atom distance'
         write(10,*)' was found betweem atoms ',n,' and ',mmax,
     .               "  ",lab(n),"--",lab(mmax)

        endif
c add a disjoint H atom to the EXTRABONDS
c bonded to the nearest heavy
       else
       write(66,*)' Disjoint H atom extra bonding allocated'
       write(66,*)n,mmax
       numextra=numextra+1
       mbex(numextra)=mmax
       nbex(numextra)=n
       endif
566    continue
561   enddo

      write(66,*)' Bonds involving hydrogen'
      do n=1,natom
       if(nbonds(n).gt.0)then
       do m=1,nbonds(n)
       write(66,*)n,bonds(n,m)
       enddo
       endif
      enddo

 
c now get the atom-atom distances and assign bonds
      do n=1,natom
      shortest(n)=1.d6
      nshort(n)=0
      enddo
      do n=1,natom-1
       if(lab(n).eq.'H')go to 564
      do m=n+1,natom
       if(lab(m).eq.'H')go to 563
       dis=0.d0
       do k=1,3
        dis=dis+(c(n,k)-c(m,k))**2
       enddo
       dis=sqrt(dis)
       if(dis.lt.shortest(n))then
        shortest(n)=dis
        nshort(n)=m
       endif
       if(dis.lt.shortest(m))then
        shortest(m)=dis
        nshort(m)=n
       endif
c rad(n)+rad(m)+0.4d0 is the maximum allowed length for a single bond
       dis=dis/(rad(n)+rad(m)+0.4d0)
       if(dis.le.1.d0)then
        nbonds(n)=nbonds(n)+1
        nbonds(m)=nbonds(m)+1
        bonds(n,nbonds(n))=m
        bonds(m,nbonds(m))=n
        numberofbonds=numberofbonds+1
        if(dis.le.0.6d0)then
         write(10,*)' A very short (unphysical?) atom-atom distance'
         write(10,*)' was found betweem atoms ',n,' and ',m,
     .               "  ",lab(n),"--",lab(m)
        endif
c bd contains the maximim allowed length for a double bond 
        bdlim=bd(num(n),num(m))
        dist1=dis*(rad(n)+rad(m)+0.4d0)/bdlim
c dist < 1 implies a double bond
        dist(n,nbonds(n))=dist1
        dist(m,nbonds(m))=dist1
       endif
563   continue
      enddo
564   continue
      enddo

c later, if nbonds=0 for some heavy, an extrabond will attach this
c atom to its nearest heavy

c put bond multiplicity to 1
      write(66,*)' single bonds allocated'
      do n=1,natom
       if(nbonds(n).gt.0)then
        do m=1,nbonds(n)
         mult(n,m)=1
      write(66,*)n,bonds(n,m)
        enddo
       endif
      enddo

c get the total number of bonds, 
c not counting bonds to special atoms (metals)
      do n=1,natom
       totbonds(n)=0
       if(nbonds(n).gt.0)then
        do m=1,nbonds(n)
c         if(special(bonds(n,m)).eq.0)then   14/11/17
         if(nval(bonds(n,m)).ne.0)then

          totbonds(n)=totbonds(n)+1
         endif
        enddo
       endif
      enddo

c try to account for pentavalent members of the third group
      do n=1,natom
       if(nval(n).eq.3.and.totbonds(n).gt.3)then
        ival(n)=5
        write(66,*)' Atom number ',n,lab(n)
        write(66,*)' has been set to be pentavalent'
       endif
c try to account for hexavalent group 6 above Oxygen
       if(nval(n).eq.6.and.totbonds(n).gt.2.and.num(n).gt.8)then
        ival(n)=6
        write(66,*)' Atom number ',n,lab(n)
        write(66,*)' has been set to be hexavalent'
       endif
c try to account for heptavalent group 7 above flurine
       if(nval(n).eq.7.and.totbonds(n).gt.1.and.num(n).gt.9)then
        ival(n)=7
        write(66,*)' Atom number ',n,lab(n)
        write(66,*)' has been set to be heptavalent'
       endif
      enddo


c now evaluate a preliminary charge for each atom, except metals
      do n=1,natom
       if(special(n).eq.0)ncharge(n)=totbonds(n)-ival(n)
14    enddo
      write(66,*)
      write(66,*)' preliminary charges for each atom'
      do n=1,natom
       write(66,*)n,ncharge(n)
      enddo

c work out how many ways an unsatisfied atom could make a multiple bond
c atoms with least choice available will be examined first
      do n=1,natom
       ichoice(n)=0
       if(nbonds(n).gt.0)then
        do m=1,nbonds(n)
         m1=bonds(n,m)
         if(ncharge(m1).lt.0.or.special(m1).eq.1)ichoice(n)=ichoice(n)+1
        enddo
       endif
      enddo

c order the atoms according to the charge, choice & electronegativity
      do n=1,natom
       weight(n)=dble(float(ncharge(n)))
     .           +dble(float(ichoice(n)))/4.d0
     .           +elect(n)/10.d0
      enddo

c consider atoms bonded to metals last
      wmax=-1.d6
      do n=1,natom
       if(weight(n).gt.wmax)wmax=weight(n)
      enddo
      do n=1,natom
       if(nbonds(n).gt.0)then
        do m=1,nbonds(n)
         m1=bonds(n,m)
         if(special(m1).eq.1)weight(n)=wmax+1.d0
        enddo
       endif
      enddo

c find the order of ascending weight
      call indexx(natom,weight,indx)

      do n1=1,natom
       n=indx(n1)
111    continue
       if(lab(n).eq.'H') go to 100
       if(special(n).gt.0)go to 100
       if(ncharge(n).ge.0)go to 100
       if(nbonds(n).eq.0)go to 100

       if(nbonds(n).gt.1)then
c order the atoms bonded to n according to charge,choice&electronegavitity
       call reordersubs(n,natom,nbonds,nval,bonds,mult,special,
     .                  ncharge,dist,elect,maxbond)
       endif
c      write(66,*)'atom, bonded atoms'
c      write(66,*)n,(bonds(n,j),j=1,nbonds(n))

c if this atom is bonded to a metal, then leave the charge unchanged
c if all the other substituents are neutral
       neispec=0
       nch=0
        do m=1,nbonds(n)
         m1=bonds(n,m)
         if(special(m1).eq.1)neispec=1
         if(special(m1).eq.0.and.ncharge(m1).ne.0)nch=1
        enddo
       if(neispec.eq.1.and.nch.eq.0)go to 100

c try to form bonds with the "best" neighbour
       do m=1,nbonds(n)
        m1=bonds(n,m)
        if(lab(m1).eq.'H')go to 99
        if(special(m1).gt.0)go to 99
c don't let group 3 make multiple bonds with neutral atoms
        if(nval(n).eq.3.and.ncharge(m1).eq.0)go to 99
        if((ncharge(m1).lt.0).or.
     .    (ncharge(m1).eq.0.and.elect(m1).gt.elect(n)+1.d-2))then
c increase the multiplicity
         mult(n,m)=mult(n,m)+1
         do k=1,nbonds(m1)
         if(bonds(m1,k).eq.n)mult(m1,k)=mult(m1,k)+1
         enddo
         ncharge(n)=ncharge(n)+1
         ncharge(m1)=ncharge(m1)+1
c        write(66,*)' bonded ',n,m1
         if(ncharge(n).lt.0)then
c         write(66,*)' still negative'
          go to 111
         else
c find the neighbour of n and m1 that has least  choice (excepting
c metal neighbours)
          call nextatom(n,m1,natom,nbonds,bonds,ncharge,special,
     .                    elect,dist,next,maxbond)
          if(next.ne.0)then
           n=next
           go to 111
          endif

         endif
        endif
         if(ncharge(n).ge.0)go to 100
c end the neighbour nbonds(n) loop
99     enddo
c end the natom loop
100   enddo
c check for pentavalent nitrogen
      do n=1,natom
       if(lab(n).ne.'N')go to 101
       if(ncharge(n).ne.0)go to 101
       neg=0
       if(nbonds(n).gt.0)then
        do m=1,nbonds(n)
         m1=bonds(n,m)
         if(ncharge(m1).lt.0)neg=neg+1
        enddo
       endif
       if(neg.lt.2)go to 101
       write(66,*)' Atom number ',n,lab(n)
       write(66,*)' has been made pentavalent'
       call reordersubs(n,natom,nbonds,nval,bonds,mult,special,
     .                  ncharge,dist,elect,maxbond)
       neg=0
       if(nbonds(n).gt.0)then
        do m=1,nbonds(n)
         m1=bonds(n,m)
         if(ncharge(m1).lt.0)then
          neg=neg+1
          ncharge(m1)=ncharge(m1)+1
          mult(n,m)=mult(n,m)+1
           do k=1,nbonds(m1)
            if(bonds(m1,k).eq.n)mult(m1,k)=mult(m1,k)+1
           enddo
          if(neg.eq.2)go to 101
         endif
        enddo
       endif
101   enddo

c check for pentavalent Phosphorus (or above)
c making 5 bonds to a group 5, if that removes charges
      do n=1,natom
       if(lab(n).eq.'N')go to 103
       if(nval(n).ne.5)go to 103
       neg=0
       nbtot=0
       if(nbonds(n).gt.0)then
        do m=1,nbonds(n)
         m1=bonds(n,m)
         nbtot=nbtot+mult(n,m)
         if(special(m1).eq.0.and.ncharge(m1).lt.0)then
          neg=neg-ncharge(m1)
         endif
        enddo
       endif
       if(nbtot+neg.ge.5)then
       write(66,*)' Atom number ',n,lab(n)
       write(66,*)' has been made pentavalent'
       call reordersubs(n,natom,nbonds,nval,bonds,mult,special,
     .                  ncharge,dist,elect,maxbond)
        ntot=nbtot
        if(nbonds(n).gt.0)then
         do m=1,nbonds(n)
          m1=bonds(n,m)
          if(special(m1).eq.0.and.ncharge(m1).lt.0)then
           jtot=ntot-ncharge(m1)-5
           if(jtot.le.0)then
            ntot=ntot-ncharge(m1)
            mult(n,m)=mult(n,m)-ncharge(m1)
            do k=1,nbonds(m1)
            if(bonds(m1,k).eq.n)mult(m1,k)=mult(n,m)
            enddo
            ncharge(m1)=0
           endif
          endif
         enddo
        endif
        ncharge(n)=0
       endif
103   enddo  

c this completes the assignment of charges

c we now make some bonds multiples on the basis of bondlength
c for example, up to this point, benzene has 3 single CC bonds
c and 3 double CC bonds. All must be made double.

      write(66,*)' Charges finished'
      write(66,*)' charges before special adjustments are'
      do n=1,natom
       write(66,*)n, ncharge(n)
      enddo


      go to 568
c unless included in IN_SPECIFIEDCHARGES, a negative atom
c of group 4 is made neutral
c and an isolated H atom is made neutral
      do n=1,natom
       if(ncharge(n).eq.0)go to 567
       if(nval(n).eq.4.and.ncharge(n).lt.0)then
        do k=1,nspec
         if(natspec(k).eq.n)go to 567
        enddo
        ncharge(n)=0
        write(10,*)' The following group 4 atom has been made neutral'
        write(10,*)n
        write(10,*)
       endif
567   enddo
568   continue

      do n=1,natom
       if(nval(n).eq.1.and.ncharge(n).lt.0)then
        do k=1,nspec
         if(natspec(k).eq.n)go to 569
        enddo
        ncharge(n)=0
        write(10,*)' The following H atom has been made neutral'
        write(10,*)n
        write(10,*)
       endif        
569   enddo



c      open(unit=1,file='IN_FIXEDCHARGES',status='unknown')
c     nfixed=0
c     read(1,*,end=600)
c     read(1,*)nfixed
c     if(nfixed.gt.0)then
c      read(1,*)
c      do n=1,nfixed
c      read(1,*)n1,i1
c      ncharge(n1)=i1
c      enddo
c     endif

c     write(66,*)' charges after IN_FIXEDCHARGES are'
c     do n=1,natom
c      write(66,*)n, ncharge(n)
c     enddo

600   continue
      
c     stop

c note any unattached atoms
      nobond=0
      do n=1,natom
       if(lab(n).ne.'H'.and.nbonds(n).eq.0)nobond=nobond+1
      enddo
      if(nobond.gt.0)then
       write(66,*)' WARNING'
       write(66,*)' The following atoms have no bonds'
       write(66,*)' Atom    Shortest distance'
       do n=1,natom
        if(lab(n).ne.'H'.and.nbonds(n).eq.0)write(66,*)n,shortest(n),nshort(n)
       enddo

       write(66,*)
       write(66,*)' We attach these atoms to their nearest neighbour'
       do n=1,natom
        if(lab(n).ne.'H'.and.nbonds(n).eq.0)then
         if(numextra.gt.0)then
          do j=1,numextra
           if(mbex(j).eq.n.and.nbex(j).eq.nshort(n))go to 1234
           if(nbex(j).eq.n.and.mbex(j).eq.nshort(n))go to 1234
          enddo
         endif
         numextra=numextra+1
         write(66,*)numextra,n,nshort(n)
         mbex(numextra)=nshort(n)
         nbex(numextra)=n
        endif
1234   continue
       enddo
      endif

c add any extra bonds
      if(numextra.gt.0)then
       write(10,*)' Unusual bonds were added to the structure'
       do n=1,numextra
        nbonds(mbex(n))=nbonds(mbex(n))+1
        bonds(mbex(n),nbonds(mbex(n)))=nbex(n)
        mult(mbex(n),nbonds(mbex(n)))=-1
        nbonds(nbex(n))=nbonds(nbex(n))+1
        bonds(nbex(n),nbonds(nbex(n)))=mbex(n)
        mult(nbex(n),nbonds(nbex(n)))=-1
        write(10,*)mbex(n),nbex(n)
       enddo
      endif
 
c zero arrays
      do n1=1,natom
       neig(n1)=0
      do n2=1,natom
       jmat(n1,n2)=0
      enddo
      enddo

c all bonds to metals are made double
c but not to other atoms which are neutral
      write(66,*)' input metals'
      do k=1,nspec
       write(66,*)natspec(k),nval(natspec(k)),nbonds(natspec(k))
      enddo
      write(66,*)' metal bonds'
      do k=1,nspec
       ligch(k)=0
       na=natspec(k)
       if(nchspec(k).eq.0)go to 888
       if(nbonds(na).gt.0)then
        do m=1,nbonds(na)
         mult(na,m)=2
         m1=bonds(na,m)
         write(66,*)na,m1,mult(na,m)
         do j=1,nbonds(m1)
         j1=bonds(m1,j)
         if(j1.eq.na)mult(m1,j)=2
         enddo
c zero the charge
c 120417: don't do this
c        if(ncharge(m1).eq.-1)ligch(k)=ligch(k)-1
c        if(special(m1).eq.0)then
c         ncharge(m1)=0
c         write(66,*)' zeroed charge for ',m1
c        endif
        enddo
       endif
888   continue
      enddo

      write(66,*)' double bonded to metals'

c all bonds shorter than the allowed maximum are made double
      do n=1,natom
      if(nbonds(n).gt.0)then
       do m=1,nbonds(n)
        m1=bonds(n,m)
        if(dist(n,m).lt.1.d0)then
         mult(n,m)=2
         do k=1,nbonds(m1)
          k1=bonds(m1,k)
          if(k1.eq.n)mult(m1,k)=2
         enddo
c        write(66,*)' short ',n,m1
        endif
       enddo
      endif
      enddo

      write(66,*)' short bonds are double'

c ensure bonds to H atoms are singles
      do n=1,natom
      if(lab(n).eq.'H')then
       if(nbonds(n).gt.0)then
        do m=1,nbonds(n)
         m1=bonds(n,m)
         mult(n,m)=1
         do k=1,nbonds(m1)
          j1=bonds(m1,k)
          if(j1.eq.n)mult(m1,k)=1
         enddo
        enddo
       endif
      endif
      enddo


      write(66,*)' bonds to H are single'


c adjust any other multiplicities as possible
      call makenvalsingle(natom,nval,ival,nbonds,
     . bonds,mult,ncharge,special,maxbond)

      write(66,*)' makenvalsingle completed'

      call adjustsingle(natom,nbonds,bonds,mult)
      write(66,*)' adjustsingle completed'

c c all bonds to charged groups are made double
c       do n=1,natom
c        if(ncharge(n).ne.0)then
c         if(nbonds(n).gt.0)then
c          do m=1,nbonds(n)
c           m1=bonds(n,m)
c           if(lab(m1).ne.'H')then
c            mult(n,m)=2
c            do k=1,nbonds(m1)
c             j1=bonds(m1,k)
c             if(j1.eq.n)mult(m1,k)=2
c            enddo
c           endif
c          enddo
c         endif
c        endif
c       enddo
501   format(2I10,f10.4)

c       write(66,*)' bonds to charged groups made double'
 
c make any additional double bonds that were requested on input
      if(ndouble.gt.0)then
       do n=1,ndouble
        do m=1,nbonds(dmbd(n))
         m1=bonds(dmbd(n),m)
         if(m1.eq.dnbd(n))mult(dmbd(n),m)=2
        enddo
        do m=1,nbonds(dnbd(n))
         m1=bonds(dnbd(n),m)
         if(m1.eq.dmbd(n))mult(dnbd(n),m)=2
        enddo
       enddo
      endif



      write(66,*)' IN_DOUBLE accounted for'

c include input requestes for classes of double bonds
      call preservebonds(natom,nbonds,bonds,mult,ncharge,
     .                         special,nval,lab,maxbond)

      write(66,*)' PRESERVED BONDS done'

c we will output charged atoms and any atoms connected to them
c by multiple bonds + any hydrogens


      do n1=1,natom
       jmat(n1,n1)=1
       nused(n1)=0
      if(nbonds(n1).gt.0)then
       do m=1,nbonds(n1)
        m1=bonds(n1,m)
        if(mult(n1,m).gt.1)then
         jmat(n1,m1)=1
         jmat(m1,n1)=1
        endif
        if(lab(m1).eq.'H')then
         jmat(n1,m1)=1
         jmat(m1,n1)=1
        endif
       enddo
      endif
      enddo

c jmat = 1 for multiple bonds, or hydrogens


c output the metal-containing groups
      open(unit=1,file='OUT_METALGROUPS',status='unknown')
      write(1,*)"The following metal-containing groups were found"
      write(1,*)"Check that the metal charge is correct "
      write(1,*)
      do k=1,nspec
       na=natspec(k)
c find all the atoms multiple bonded to metals (neig =1)
       call specgrp(na,natom,jmat,neig)
c count the number of such atoms
c and note that each jas been assigned to a metal group (nused=1)
       mber=0
       ntot=0
       do n=1,natom
        if(neig(n).eq.1.and.nused(n).eq.0)then
         mber=mber+1
         ntot=ntot+ncharge(n)
         nused(n)=1
         if(nbonds(n).gt.0)then
          do m=1,nbonds(n)
           m1=bonds(n,m)
           if(neig(m1).eq.0)then
            mber=mber+1
            ntot=ntot+ncharge(m1)
           endif
          enddo
         endif
c estimate the total ligand charge
c        if(special(n).eq.0)ntot=ntot+ncharge(n)
        endif
       enddo
       if(mber.eq.0)go to 112
c find the atom number for the first atom in this group
        do m=1,natom
         if(neig(m).eq.1)then
          mfirst=m
          go to 6548
         endif
        enddo
6548    continue
c       write(1,*)' Ligand charge has been estimated as ',ntot+ligch(k)
       write(1,*)"The net charge on this group is ",ntot
       write(1,*)' The first atom in this group is atom number ',mfirst
       write(1,*)"The number of atoms in the group is "
       write(1,*)mber
       write(1,*)' The group coordinates are'
c write out the atoms
       do n=1,natom
        if(neig(n).eq.1)then
         write(1,500)lab(n),(c(n,j),j=1,3)
         if(nbonds(n).gt.0)then
          do m=1,nbonds(n)
           m1=bonds(n,m)
           if(neig(m1).eq.0)then
            dx=0.d0
            do j=1,3
             x(j)=c(m1,j)-c(n,j)
             dx=dx+x(j)**2
            enddo
             dx=sqrt(dx)
            write(1,500)'H ',(c(n,j)+x(j)/dx,j=1,3)
           endif
          enddo
         endif
        endif
       enddo

112   enddo 

      close(unit=1)

c     write(66,*)' after metal groups '
c check to see if charges are connected by a sequence of double bonds
c if this sequence does not involve special (metal) atoms, or
c atoms in metal groups, then if the pair of charges are equal and 
c opposite, we set both to zero.
      call cancel(natom,nbonds,bonds,mult,special,
     .                  ncharge,nused,maxbond)


c now output the known charged groups, other than metals
      open(unit=1,file='OUT_NONMETALCHARGES',status='unknown')
      write(1,*)'The following other charged groups were found'
      write(1,*)
      do n=1,natom
       if(ncharge(n).ne.0.and.nused(n).eq.0)then
        call specgrp(n,natom,jmat,neig)
c in order to add groups adjacent to charges to the charged group,
c we identify these atoms, make the bonds double and
c recalculate the group.
        do n1=1,natom
         if(nbonds(n1).eq.0)go to 222
c find out if n1 is bonded to the charged group
         match=0
         do m=1,nbonds(n1)
          m1=bonds(n1,m)
          if(neig(m1).eq.1)match=1
         enddo
         if(match.eq.0)go to 222
         do m=1,nbonds(n1)
          m1=bonds(n1,m)
          if(neig(m1).eq.1)then
           if(lab(m1).ne.'H'.and.lab(n1).ne.'H')mult(n1,m)=2
           do k=1,nbonds(m1)
            k1=bonds(m1,k)
            if(k1.eq.n1)then
              if(lab(m1).ne.'H'.and.lab(k1).ne.'H')mult(m1,k)=2
            endif
           enddo
           jmat(n1,m1)=1
           jmat(m1,n1)=1
          endif
          if(lab(m1).eq.'H')then
           jmat(n1,m1)=1
           jmat(m1,n1)=1
          endif
         enddo

222     enddo
        call specgrp(n,natom,jmat,neig)

c count the atoms, and mark which they are Unused=1)
        mber=0
        ntot=0
        do m=1,natom
         if(neig(m).eq.1)then
          mber=mber+1
          nused(m)=1
          ntot=ntot+ncharge(m)
          if(nbonds(m).gt.0)then
           do m2=1,nbonds(m)
            m1=bonds(m,m2)
            if(neig(m1).eq.0)mber=mber+1
           enddo
          endif

         endif
        enddo
        if(mber.eq.0)go to 113

c find the atom number for the first atom in this group
        do m=1,natom
         if(neig(m).eq.1)then
          mfirst=m
          go to 6547
         endif
        enddo
6547    continue
        write(1,*)' Total charge for this group has been estimated as ',ntot
        write(1,*)' The first atom in this group is atom number ',mfirst
        write(1,*)'The number of atoms in the group is '
        write(1,*)mber
        write(1,*)'The atomic coordinates are'
        do m=1,natom
         if(neig(m).eq.1)then
          write(1,500)lab(m),(c(m,j),j=1,3)

          if(nbonds(m).gt.0)then
           do m2=1,nbonds(m)
            m1=bonds(m,m2)
            if(neig(m1).eq.0)then
             dx=0.d0
             do j=1,3
              x(j)=c(m1,j)-c(m,j)
              dx=dx+x(j)**2
             enddo
              dx=sqrt(dx)
             write(1,500)'H ',(c(m,j)+x(j)/dx,j=1,3)
            endif
           enddo
          endif

         endif
        enddo
      write(1,*)
       endif
113   enddo

      close(unit=1)
500   format(a2,3F10.4)



c output the bonding
       write(66,*)' The bonding used to allocate formal charges is:'
       do n=1,natom
        if(nbonds(n).gt.0)then
        do m=1,nbonds(n)
         if(bonds(n,m).gt.n)write(66,*)n,bonds(n,m),mult(n,m)
        enddo
        endif
       enddo
       write(66,*)

c output the charges
       numchgs=0
       nchtot=0
       do n=1,natom
       if(iabs(ncharge(n)).gt.0.or.special(n).eq.1)numchgs=numchgs+1
       enddo
       write(66,*)'  The number of charged atoms is'
       write(66,*)numchgs
       write(66,*)' The charged atoms are'
       do n=1,natom
        if(ncharge(n).ne.0.or.special(n).eq.1)then
        write(66,*)n,ncharge(n),lab(n)
        nchtot=nchtot+ncharge(n)
        endif
       enddo
       write(66,*)
       write(6,*)' The total molecular charge = ',nchtot

c if nflag=0, we create an IN_SPECIFIEDCHARGES file
      if(nflag.eq.0)then
       open(unit=3,file='IN_SPECIFIED_CHARGES',status='unknown')
       if(nspec.eq.0)then
        write(3,*)' There are no charges to be specified'
       else
        write(3,*)' The following allocated charges need to be viewed'
        write(3,*)' and altered as necessary'
        write(3,*)' The atom numbers and charges are:'
        do i=1,nspec
         write(3,*)natspec(i),nchspec(i)
        enddo
       endif
      endif

      nchall=0
      if(nchall.eq.0)then

c create the IN_CHARGES file
       open(unit=2,file='IN_CHARGES',status='unknown')
       write(2,*)'  The number of charged atoms is'
       write(2,*)numchgs
       write(2,*)'  The atom numbers and associated charge'
       do n=1,natom
        if(ncharge(n).ne.0.or.special(n).eq.1)then
         write(2,*)n,ncharge(n)
        endif
       enddo
       close(unit=2)

      else

c put charges on amide groups
      ic=0
      if(preserve(1).eq.'F')go to 22
      do n=1,natom

       if(nbonds(n).lt.2)go to 26
       if(lab(n).ne.'C')go to 26
       mark1=0
       mark2=0
       do m=1,nbonds(n)
        m1=bonds(n,m)
        if(lab(m1).eq.'N')mark1=m1
        if(lab(m1).eq.'O'.and.mult(n,m).eq.2)mark2=m1
       enddo
       if(mark1*mark2.eq.0)go to 26
       if(ncharge(mark1).ne.0.or.ncharge(mark2).ne.0)go to 26
       ic=ic+1
26    enddo
22    continue

      numchgstot=ic+numchgs

c create the IN_CHARGES file
       open(unit=2,file='IN_CHARGES',status='unknown')
       write(2,*)'  The number of charged atoms is'
       write(2,*)numchgstot
       write(2,*)'  The atom numbers and associated charge'
       if(numchgstot.eq.0)go to 23
       do n=1,natom
        if(ncharge(n).ne.0.or.special(n).eq.1)then
         write(2,*)n,ncharge(n)
        endif
       enddo
       do n=1,natom
        if(nbonds(n).lt.2)go to 24
        if(lab(n).ne.'C')go to 24
        mark1=0
        mark2=0
        do m=1,nbonds(n)
         m1=bonds(n,m)
         if(lab(m1).eq.'N')mark1=m
         if(lab(m1).eq.'O'.and.mult(n,m).eq.2)mark2=m
        enddo
        if(mark1*mark2.eq.0)go to 24
        if(ncharge(mark1).ne.0.or.ncharge(mark2).ne.0)go to 24
        write(2,*)n,0
24     enddo

23     close(unit=2)


cc put charges on all heteroatoms
c      ic=0
c      do n=1,natom
c       if(ncharge(n).eq.0.and.special(n).eq.0)then
c        if(lab(n).ne.'H'.and.nval(n).ne.4.and.nval(n).ne.8)ic=ic+1
c       endif
c      enddo

c      numchgstot=ic+numchgs

cc create the IN_CHARGES file
c       open(unit=2,file='IN_CHARGES',status='unknown')
c       write(2,*)'  The number of charged atoms is'
c       write(2,*)numchgstot
c       write(2,*)'  The atom numbers and associated charge'
c       do n=1,natom
c        if(ncharge(n).ne.0.or.special(n).eq.1)then
c         write(2,*)n,ncharge(n)
c        endif
c       enddo
c      do n=1,natom
c       if(ncharge(n).eq.0.and.special(n).eq.0)then
c       if(lab(n).ne.'H'.and.nval(n).ne.4.and.nval(n).ne.8)write(2,*)n,0
c       endif
c      enddo

c       close(unit=2)

      endif

      if(nspec.eq.0.or.nflag.eq.1)then

c we create the name.* files
      open(unit=2,file='name.xyz',status='unknown')
      open(unit=3,file='name.cart',status='unknown')
      write(2,*)natom
      write(2,*)
      do n=1,natom
      write(2,306)lab(n),(c(n,k),k=1,3)
      write(3,302)lab(n),(c(n,k),k=1,3)
      enddo
302   format(a2,3f13.6)
306   format(a2,3d25.16)
      close(unit=2)
      close(unit=3)

      open(unit=2,file='name.mol',status='unknown')
      write(2,*)' name.mol generated by Preparegeom'
      write(2,*)
      write(2,*)' Number of atoms and bonds'
      write(2,*)natom,numberofbonds+numextra
      do n=1,natom
      write(2,303)(c(n,k),k=1,3),' ',lab(n),ncharge(n)
      enddo
303   format(3f10.4,a1,a2,i4)
      do n=1,natom-1
       if(nbonds(n).eq.0)go to 305
       do m=1,nbonds(n)
        m1=bonds(n,m)
        if(m1.gt.n)write(2,304)n,m1,mult(n,m)
       enddo
305   enddo
304   format(3I9)

      write(2,*)'M  END'
      close(unit=2)

c close the nflag if
      endif

      end

      subroutine reordersubs(n,natom,nbonds,nval,bonds,mult,special,
     .                      ncharge,dist,elect,maxbond)
      implicit double precision(a-h,o-z)

      dimension ncharge(natom),nbonds(natom),nval(natom)
      integer bonds(natom,maxbond),mult(natom,maxbond)
      integer ind(maxbond),neigh(maxbond)
      integer special(natom)

      dimension dist(natom,maxbond),elect(natom)

      dimension choice(maxbond),eneigh(maxbond),ene(maxbond)
      dimension dtemp(maxbond)


        do m=1,nbonds(n)
         m1=bonds(n,m)
         nb=nbonds(m1)
         tot=0.d0
          do k=1,nb
           k1=bonds(m1,k)
           if(special(k1).eq.1)tot=tot+2.d0
           if(ncharge(k1).lt.0.and.dist(m1,k).lt.1.d0)tot=tot+1.d0
           if(nval(k1).eq.5)then
c calculate total bonds
            j1=0
            do j=1,nbonds(k1)
             j1=j1+mult(k1,j)
            enddo
c ask if this group 5 could become pentavalent
            if(ncharge(m1).lt.0.and.j1.lt.5)tot=tot+1.d0
           endif
          enddo
         choice(m)=tot-1.d0
        enddo
        do m=1,nbonds(n)
         m1=bonds(n,m)
         eneigh(m)=dble(float(ncharge(m1)))
     .            +choice(m)/4.d0+elect(m1)/10.d0
         if(dist(n,m).gt.1.d0)eneigh(m)=eneigh(m)+20.d0
c       if(n.eq.12.or.n.eq.13)then
c        write(66,555)m1,ncharge(m1),dist(n,m),choice(m),elect(m1),
c    .               eneigh(m)
c       endif
555      format(2i10,4f10.4)
        enddo

c if only one neighbour has negative charge then make it the
c lowest weight
        nneg=0
        do m=1,nbonds(n)
         m1=bonds(n,m)
         if(ncharge(m1).lt.0)nneg=nneg+1
        enddo
        if(nneg.eq.1)then
         do m=1,nbonds(n)
          m1=bonds(n,m)
          if(ncharge(m1).lt.0)eneigh(m)=-10.d0
         enddo
        endif

        if(nbonds(n).lt.maxbond)then
         do m=nbonds(n)+1,maxbond
          eneigh(m)=100.d0
         enddo
        endif

        call indexx(maxbond,eneigh,ind)

        do m=1,nbonds(n)
         neigh(m)=bonds(n,ind(m))
         dtemp(m)=dist(n,ind(m))
         ene(m)=eneigh(ind(m))
        enddo
        do m=1,nbonds(n)
         bonds(n,m)=neigh(m)
         dist(n,m)=dtemp(m)
         eneigh(m)=ene(m)
        enddo

      return
      end

      subroutine nextatom(n,m1,natom,nbonds,bonds,ncharge,special,
     .                    elect,dist,next,maxbond)
      implicit double precision(a-h,o-z)

      dimension nbonds(natom),ncharge(natom)

      integer special(natom),bonds(natom,maxbond)

      dimension dist(natom,maxbond),elect(natom)

      dimension nscore(maxbond)

      do i=1,maxbond
       nscore(i)=0
      enddo

      do n1=1,nbonds(n)
       k1=bonds(n,n1)
       if(k1.eq.m1)go to 1
       if(ncharge(k1).ge.0)go to 1
       if(special(k1).eq.1)then
        nscore(n1)=100
        go to 1
       endif
       do j1=1,nbonds(k1)
        k2=bonds(k1,j1)
        if(special(k2).eq.1)then
         nscore(n1)=100
         go to 1
        endif
        if(k2.eq.n)go to 2
        d1=dist(k1,j1)
        if(d1.gt.1.d0)go to 2
        if(ncharge(k2).lt.0.or.
     .   (ncharge(k2).eq.0.and.elect(k2).gt.elect(k1)+1.d-2))then
         nscore(n1)=nscore(n1)+1
        endif
2      enddo
1     enddo

      nposs1=0
      nsc1=100
      do n1=1,nbonds(n)
       if(nscore(n1).gt.0.and.nscore(n1).lt.nsc1)then
        nsc1=nscore(n1)
        nposs1=n1
       endif
      enddo

c other end of the bond
      do i=1,maxbond
       nscore(i)=0
      enddo

      do n1=1,nbonds(m1)
       k1=bonds(m1,n1)
       if(k1.eq.n)go to 3
       if(ncharge(k1).ge.0)go to 3
       if(special(k1).eq.1)then
        nscore(n1)=100
        go to 3
       endif
       do j1=1,nbonds(k1)
        k2=bonds(k1,j1)
        if(k2.eq.m1)go to 4
        d1=dist(k1,j1)
        if(d1.gt.1.d0)go to 4
        if(ncharge(k2).lt.0.or.
     .   (ncharge(k2).eq.0.and.elect(k2).gt.elect(k1)+1.d-2))then
         nscore(n1)=nscore(n1)+1
        endif
4      enddo
3     enddo

      nposs2=0
      nsc2=100
      do n1=1,nbonds(m1)
       if(nscore(n1).gt.0.and.nscore(n1).lt.nsc2)then
        nsc2=nscore(n1)
        nposs2=n1
       endif
      enddo  

c     write(66,*)' nposs1, nsc1, nposs2, nsc2'
c     write(66,*)nposs1, nsc1, nposs2, nsc2

      if(ncs1.eq.100.and.nsc2.eq.100)then
       next=0
       return
      endif

      if(nposs1.eq.0.and.nposs2.eq.0)then
       next=0
       return
      endif

      if(nsc1.eq.100)then
       next=bonds(m1,nposs2)
       return
      endif

      if(nsc2.eq.100)then
       next=bonds(n,nposs1)
       return
      endif

      if(nsc1.le.nsc2)then
       next=bonds(n,nposs1)
      else
       next=bonds(m1,nposs2)
      endif

      return
      end


      subroutine makenvalsingle(natom,nval,ival,nbonds,
     . bonds,mult,ncharge,special,maxbond)
      implicit double precision(a-h,o-z)

      dimension nval(natom),ival(natom),nbonds(natom),ncharge(natom)

      integer bonds(natom,maxbond),special(natom)

      dimension mult(natom,maxbond),normal(natom),jval(natom)


      

c some bonds may have been made multiple on the basis of length.

c However, if the number of bonds to an atom is consistent
c with normal valence, then the bonds are made single


c recompute the normal valence
      do n=1,natom
       if(special(n).eq.1)go to 1
       if(nval(n).le.4)jval(n)=nval(n)
       if(nval(n).gt.4)jval(n)=8-nval(n)
1     enddo

c compare with ival, as some atoms may have been allocated
c extraordinary valence
      do n=1,natom
       normal(n)=0
       if(jval(n).ne.ival(n))normal(n)=1
      enddo

c adjust multiplicity only if both atoms in a bond
c have nbonds = ival
      do n=1,natom
       if(ncharge(n).ne.0)go to 2
       if(special(n).eq.1)go to 2
       if(normal(n).eq.1)go to 2
       if(nbonds(n).eq.ival(n))then
        do m=1,nbonds(n)
         m1=bonds(n,m)
         if(special(m1).eq.1)go to 2
         if(normal(m1).eq.1)go to 2
         if(nbonds(m1).eq.ival(m1))then
          mult(n,m)=1
          do k=1,nbonds(m1)
           k1=bonds(m1,k)
           if(k1.eq.n)mult(m1,k)=1
          enddo
         endif
        enddo
       endif
2     enddo

      return
      end

      subroutine preservebonds(natom,nbonds,bonds,mult,ncharge,
     .                         special,nval,lab,maxbond)
      implicit double precision(a-h,o-z)

      dimension nbonds(natom),ncharge(natom),mult(natom,maxbond),nval(natom)

      integer bonds(natom,maxbond),special(natom)

      character*2 lab(natom)

      character*1, preserve(2)

c get the info from IN_PRESERVE_BONDS

      preserve(1)='T'
      preserve(2)='T'

      open(unit=1,file='IN_PRESERVE_BONDS',status='unknown')
      read(1,*,end=5)
      do n=1,2
      read(1,*)preserve(n)
      enddo
5     close(unit=1)

c where there are more than one halogen bonded to an atom,
c make these double bonds

      if(preserve(2).eq.'T')then

      do n=1,natom
       if(nbonds(n).lt.2)go to 1
       nhal=0
       do m=1,nbonds(n)
        m1=bonds(n,m)
        if(nval(m1).eq.7)nhal=nhal+1
       enddo
       if(nhal.lt.2)go to 1
       do m=1,nbonds(n)
        m1=bonds(n,m)
        if(nval(m1).eq.7)then
         mult(n,m)=2
         do k=1,nbonds(m1)
          k1=bonds(m1,k)
          if(k1.eq.n)mult(m1,k)=2
         enddo
        endif
       enddo
1     enddo

      endif

c make amide bonds double


      do n=1,natom
       if(nbonds(n).lt.2)go to 2
       if(lab(n).ne.'C')go to 2
       mark1=0
       mark2=0
       do m=1,nbonds(n)
        m1=bonds(n,m)
        if(lab(m1).eq.'N')mark1=m
        if(lab(m1).eq.'O'.and.mult(n,m).eq.2)mark2=m
       enddo
       if(mark1*mark2.eq.0)go to 2
       if(preserve(1).eq.'T')then
        do m=1,nbonds(n)
         m1=bonds(n,m)
         if(lab(m1).eq.'N')then
          mult(n,m)=2
          do k=1,nbonds(m1)
           k1=bonds(m1,k)
           if(k1.eq.n)mult(m1,k)=2
          enddo
         endif
        enddo
       endif
       if(preserve(1).eq.'F')then
        do m=1,nbonds(n)
         m1=bonds(n,m)
         if(lab(m1).eq.'N')then
          mult(n,m)=1
          do k=1,nbonds(m1)
           k1=bonds(m1,k)
           if(k1.eq.n)mult(m1,k)=1
          enddo
         endif
        enddo
       endif
2     enddo
      

      return
      end

      subroutine cancel(natom,nbonds,bonds,mult,special,
     .                  ncharge,nused,maxbond)
      implicit double precision(a-h,o-z)

      dimension nbonds(natom),mult(natom,maxbond),ncharge(natom)
      dimension nused(natom),ndone(natom)

      integer special(natom), bonds(natom,maxbond)

      dimension nmat(natom,natom),nat(natom)

c come here after the metal groups are found

      nch=0
      do n=1,natom
       if(nused(n).eq.1)go to 1
       if(ncharge(n).ne.0)nch=nch+1
1     enddo

      if(nch.lt.2)return


      do n=1,natom
      do m=1,natom
       nmat(n,m)=0
      enddo
      enddo

c set up the double bond connectivity
      do n=1,natom
       if(special(n).eq.1)go to 11
       if(nbonds(n).eq.0)go to 11
       do m=1,nbonds(n)
        if(mult(n,m).eq.1)go to 22
        m1=bonds(n,m)
        if(special(m1).eq.1)go to 22
        if(nbonds(m1).eq.1)go to 22
        nmat(n,m1)=1
        nmat(m1,n)=1
22     enddo
11    enddo

      do n=1,natom
       ndone(n)=0
      enddo

      do n=1,natom
       if(nused(n).eq.1)go to 2
       if(nbonds(n).eq.0)go to 2
       if(ncharge(n).ne.0)then
        nat(1)=n
        nf=1
        nfrag=1
       else
        go to 2
       endif
10      do m=1,nbonds(nat(nfrag))
         m1=bonds(nat(nfrag),m)
         if(nmat(nat(nfrag),m1).eq.1)then
          if(nbonds(m1).gt.1.and.ndone(m1).eq.0)then
           nf=nf+1
           nat(nf)=m1
          endif
         endif
        enddo
        ndone(nat(nfrag))=1
        if(nfrag.lt.nf)then
         nfrag=nfrag+1
         go to 10
        endif
        if(nf.lt.2)go to 2

        do i=2,nf
         if(ncharge(nat(1))+ncharge(nat(i)).eq.0)then
          ncharge(nat(1))=0
          ncharge(nat(i))=0
          go to 2
         endif
        enddo
      
2     enddo   

      return
      end
