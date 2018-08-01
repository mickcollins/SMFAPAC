      program convert2dalton
      implicit double precision(a-h,o-z)

      real*8, allocatable :: c(:,:)

      integer, allocatable :: num(:)
      integer miss(1000)

      character*60, allocatable :: blines(:)

      character*2, allocatable :: lab(:),atom(:)

      character*60 line,comm,chcut,dalbasis

      character*1 exclam

      integer alab2anum

      exclam="!"
100   format(a30)

      do n=1,1000
       miss(n)=1
      enddo

c find out how many atoms
      comm="none"
      i=0
      n=1
10    continue
      read(5,100,end=20)line
      if(line(1:1).eq.exclam)then
      comm=line
       miss(n)=0
       n=n+1
       go to 11
      endif
      m=index(line,'.')
      if(m.eq.0)then
       miss(n)=0
       n=n+1
       go to 11
      endif
      i=i+1
      n=n+1
11    go to 10
20    continue
      natom=i
      nlines=n-1

      allocate(lab(natom))
      allocate(atom(natom))
      allocate(c(natom,3))
      allocate(num(natom))
      allocate(blines(natom))

      do n=1,natom
       num(n)=0
      do k=1,3
       c(n,k)=0.d0
      enddo
      enddo

c locate the charge and multiplicity line
      do n=1,nlines-1
      if(miss(n).eq.0.and.miss(n+1).eq.1)nchline=n
      enddo

      rewind(unit=5)

      i=0
      ic=0
      do n=1,nlines
       if(miss(n).eq.1)then
        i=i+1
        read(5,*)lab(i),(c(i,k),k=1,3)
c convert to Bohr
        do k=1,3
         c(i,k)=c(i,k)*1.8897259886d0
        enddo
        if(ic.eq.0)then
         ic=ic+1
         atom(ic)=lab(i)
         num(ic)=num(ic)+1
         go to 1
        endif
        if(ic.gt.0)then
         do m=1,ic
          if(atom(m).eq.lab(i))then
           num(m)=num(m)+1
           go to 1
          endif
         enddo
         ic=ic+1
         atom(ic)=lab(i)
         num(ic)=num(ic)+1
        endif
        else
         if(n.eq.nchline)then
          read(5,*)nch,mult
         else
          read(5,*)
         endif
        endif
1      enddo

      nelements=ic

c is this job for NWChem or GAMESS?
      open(unit=1,file='ABSETTINGS',status='old')
      read(1,*)
      read(1,*)npackage
      close(unit=1)
      if(npackage.eq.3)then
c it must be NWChem
      ic=0
      open(unit=1,file='NWCbasis',status='old')
3      read(1,200,end=2)line
      if(line(1:1).eq."*")then
        blines(1)=line
        nalllines=1
       go to 2
      else
       ic=ic+1
       blines(ic)=line
       go to 3
      endif
2     close(unit=1)
200   format(a60)
 
      if(ic.gt.0)nalllines=ic

      else
c it must be GAMESS
      open(unit=1,file='ABSETTINGS',status='old')
4     read(1,200,end=5)line
      dalbasis=line
      go to 4
5     close(unit=1)
c dalbasis is now the last line of ABSETTINGS
      nend=index(dalbasis," ")
      nend=nend-1
      blines(1)="*   library   "//line(1:nend)
      nalllines=1
c end the package if
      endif

      write(6,501)"ATOMBASIS"
501   format(a9)
      write(6,502)"Polarisation calculation"
502   format(a24)
      write(6,503)"nb polar job"
503   format(a12)
      write(6,500)'Atomtypes=',nelements,' Charge=',nch,' Nosymmetry'
500   format(a10,i2,a8,i2,a11)
      do n=1,nelements
       n1=alab2anum(atom(n))
       a1=float(n1)
       if(nalllines.eq.1)then
        n1=1
       else 
        do k=1,nalllines
         if(blines(k)(1:2).eq.atom(n))n1=k
        enddo
       endif
       nl=len_trim(blines(n1))
       ns=index(blines(n1),"library")
       chcut=blines(n1)(ns+7:nl)
       chcut=adjustl(chcut)

       write(6,300)'Charge=',a1,' Atoms=',num(n),' Basis= ',chcut
       do m=1,natom
        if(lab(m).eq.atom(n))write(6,400)lab(m),(c(m,k),k=1,3)
       enddo
      enddo
300   format(a7,f4.1,a7,i2,a8,a30)
400   format(a2,3f13.7)

      end



      integer function alab2anum(alab)

c given the atom label this returns the atomic number

      character*2 :: alab,atomic_labels(110)

      data atomic_labels/'H','He','Li','Be','B','C','N','O','F','Ne',   
     .'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
     . 'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb',
     . 'Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', 
     . 'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu', 
     . 'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os', 
     . 'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     .   'Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No'
     . ,'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds'/

      integer n

      do n=1,110
        if(atomic_labels(n).eq.alab)then
      alab2anum=n
      return
        endif
      enddo
      alab2anum=0

      return

      end function alab2anum
