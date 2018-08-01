      subroutine readfrag
      use momentheader

c  read the frags.out file
c which was created by the fragmentation program

      open(unit=3,file='frags.out_Lev1',status='old')

      read(3,*)
      read(3,*)nfrag
      read(3,*)
      read(3,*)
      read(3,*)
      maxatom=0
      do m=1,nfrag
c nfrag = the number of fragments (read in readdata)
       read(3,*)nat0
       if(nat0.gt.maxatom)maxatom=nat0
c nat0 = the number of real atoms in the fragment (readdata)
      enddo
      close(unit=3)
      maxatom=maxatom+maxcaps

       open(unit=1,file='OUT_Lev1_ATOMALLOCATION',status='old')
       read(1,*)
       read(1,*)NL1
       read(1,*)
       read(1,*)

       allocate(nat1(NL1))
       allocate(num1(NL1,maxatom))
       allocate(natnum1(NL1,maxatom,6))
       allocate(w1(NL1,maxatom,6))

       do n=1,NL1
        read(1,*)i1,nat1(n)
        do m=1,nat1(n)
         read(1,*)num1(n,m),
     .           (natnum1(n,m,k),w1(n,m,k),k=1,num1(n,m))
        enddo
       enddo

       close(unit=1)

       do n=1,NL1
       do m=1,nat1(n)
        if(num1(n,m).eq.1)w1(n,m,1)=1.d0
        if(num1(n,m).eq.2)w1(n,m,2)=1.d0-w1(n,m,1)
       enddo
       enddo


      return
      end

