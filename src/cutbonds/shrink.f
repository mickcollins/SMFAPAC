      subroutine shrink

      use fractdata
      implicit double precision(a-h,o-z)

      allocate(mats(nsmall))
      allocate(matr(nsmall))
      allocate(itp(nsmall))
      allocate(ibo(nsmall,nsmall))

      do i=1,nsmall
      mats(i)=0
      matr(i)=0
      do j=1,nsmall
      ibond(i,j)=0
      enddo
      enddo

      ic=0
      do i=1,numat(nfrag)
      do j=1,itype(nfrag,i)+1
       ic=ic+1
       ibond(i,j)=ib1(nfrag,ic)
      enddo
      enddo

c record the atoms in the fragment
      nafrag=numat(nfrag)
      ico=0
      do n=1,nsmall
c adjust for negative ibond values
      if(ibond(n,1).ne.0)then
      ico=ico+1
      mats(ico)=n
      matr(n)=ico
      endif
      enddo

      if(iabs(ico-nafrag).gt.0)then
      write(6,*)'ibond and numat disagree at nfrag = ',nfrag
      write(6,*)' This probably means that some groups'
      write(6,*)' are not bonded to any other group'
      write(6,*)' nfrag = ',nfrag
      write(6,*)'numat = ',numat(nfrag)
      write(6,*)(natstore(nfrag,k),k=1,numat(nfrag))
      do n=1,nsmall
      write(6,*)n,ibond(n,1)
      enddo
      stop
      endif

      do i=1,nsmall
      itp(i)=0
      do j=1,nsmall
      ibo(i,j)=0
      ilink(i,j)=0
      enddo
      enddo

c make a copy of the itype array with the new numbers
      do n=1,nafrag
      itp(n)=itype(nfrag,mats(n))
      enddo

c if all bonded atoms are negative, nocapsatall = 1
      nocapsatall=1
      do n=1,nafrag
      do j=1,itp(n)+1
      mm=ibond(mats(n),j)
      ibo(n,j)=matr(iabs(ibond(mats(n),j)))*(mm/iabs(mm))
c maintains the sign
      if(ibo(n,j).gt.0)nocapsatall=0
      enddo
      enddo

c  the link array is not necessary, so we avoid it
c  and set up the ilink array from ibo
      do n=1,nafrag
       do j=1,itp(n)+1
c       ilink(n,iabs(ibo(n,j)))=1
c       ilink(iabs(ibo(n,j)),n)=1
        ilink(n,iabs(ibo(n,j)))=ibo(n,j)/iabs(ibo(n,j))
        ilink(iabs(ibo(n,j)),n)=ibo(n,j)/iabs(ibo(n,j))
       enddo
      enddo

      return
      end
