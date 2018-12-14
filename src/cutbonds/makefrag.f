      subroutine makefrag(isgn)

      use fractdata
      implicit double precision(a-h,o-z)

c now set up new fragments

      do k=1,newfrag

      nf=nf+1
      if(nf.gt.nfragm)then
      write(6,*)' nf > nfragm in makefrag'
      stop
      endif
      isign(nf)=isign(nfrag)*isgn
      nstop(nf)=0

      do i=1,nsmall
      itype(nf,i)=0
      enddo
 
      do i=1,nsmall
      do j=1,nsmall
      ibond(i,j)=0
      enddo
      enddo

      numat(nf)=kf(k)
      do i=1,kf(k)
      natstore(nf,i)=natstore(nfrag,map(k,i))
      enddo

      do i=1,kf(k)
      isum=0
      nb1=0
      do j=1,kf(k)
      if(iabs(ilink(map(k,i),map(k,j))).gt.0)then
      nb1=nb1+1
      ibond(i,nb1)=j*ilink(map(k,i),map(k,j))
      endif
      isum=isum+iabs(ilink(map(k,i),map(k,j)))
      enddo
      itype(nf,i)=isum-1
      enddo

      ic=0
      do i=1,numat(nf)
      do j=1,itype(nf,i)+1
       ic=ic+1
       ib1(nf,ic)=ibond(i,j)
      enddo
      enddo

c check to see if this fragment is completely decomposed

      if(kf(k).le.(Level+1))nstop(nf)=1
      it=0
      do i=1,kf(k)
      if(itype(nf,i).gt.it)it=itype(nf,i)
      enddo
      if(kf(k).le.it+Level)nstop(nf)=1


      enddo

 
      return
      end
