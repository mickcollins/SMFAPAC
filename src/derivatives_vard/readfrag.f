      subroutine readfrag
      use derivheader

c  read the frags.out file
c which was created by the fragmentation program
      open(unit=3,file='frags.out',status='old')

      read(3,*)comlin
      read(3,*)
      read(3,*)comlin
      read(3,*)
      read(3,*)comlin
      do m=1,nfrag
c nfrag = the number of fragments (read in readdata)
      read(3,*)nat0(m)
c nat0 = the number of real atoms in the fragment (readdata)
      enddo
      read(3,*)comlin
      do m=1,nfrag
      read(3,*)(nat(m,i),i=1,nat0(m))
c nat = the whole-molecule atom numbers in each fragment (readdata)
      enddo
      read(3,*)comlin
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

c check the data in frags.out matches that read in from the
c ab initio fragment output, read in readdata
      do m=1,nfrag
      if(iabs(nat0(m)+ncap(m)-natomfrag(m)).gt.0)then
      write(6,*)' The number of atoms in fragment ',m
      write(6,*)' does not match the numbers in frags.out'
      stop
      endif
      enddo

c  find out how many real atoms are in the whole molecule

      natom=0
      do m=1,nfrag
      do n=1,nat0(m)
      if(nat(m,n).gt.natom)natom=nat(m,n)
      enddo
      enddo
      close(unit=3)

      return
      end

