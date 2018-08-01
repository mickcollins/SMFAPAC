      subroutine readfrag
      use derivheader

c  read the frags.out file
c which was created by the fragmentation program
      open(unit=3,file='frags.out',status='old')

      read(3,*)
      read(3,*)nfrag
      read(3,*)
      read(3,*)
      read(3,*)
      allocate(nat0(nfrag))
      do m=1,nfrag
c nfrag = the number of fragments (read in readdata)
       read(3,*)nat0(m)
c nat0 = the number of real atoms in the fragment (readdata)
      enddo
      maxatom=0
      do m=1,nfrag
       if(nat0(m).gt.maxatom)maxatom=nat0(m)
      enddo
      maxatom=maxatom+maxcaps

      allocate(nat(nfrag,maxatom))
      allocate(ncap(nfrag))
      allocate(natcap(nfrag,maxcaps,2))
      allocate(fact(nfrag,maxcaps))

      read(3,*)
      do m=1,nfrag
       read(3,*)(nat(m,i),i=1,nat0(m))
c nat = the whole-molecule atom numbers in each fragment (readdata)
      enddo
      read(3,*)
      do m=1,nfrag
       ncap(m)=0
c ncap = the number of caps in each fragment
      enddo
c  read in the identity of the atoms connected to the caps
444   continue
      read(3,*,end=445)k1,m1,n1,fac
c k1 the fragment number
c m1 and n1 are the atoms which define the cap
c fac is the weight for m1 versus n1
c ncap, natcap and fact were allocated in readdata
      ncap(k1)=ncap(k1)+1
      natcap(k1,ncap(k1),1)=m1
      natcap(k1,ncap(k1),2)=n1
      fact(k1,ncap(k1))=fac
      go to 444
445   continue

c total number of atoms in each fragment
      allocate(natomfrag(nfrag))
      do m=1,nfrag
      natomfrag(m)=nat0(m)+ncap(m)
      enddo
      close(unit=3)

      return
      end

