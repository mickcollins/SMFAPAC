      subroutine findfragelim

      use fractdata
      implicit double precision(a-h,o-z)

      integer, allocatable  :: cmat(:,:),kt(:,:),mt(:,:),ifrag(:),
     .                         ktp(:),kbo(:,:)
      natom2=natom
      allocate(cmat(natom2,natom2))
      allocate(kt(natom2,natom2))
      allocate(mt(natom2,natom2))
      allocate(ifrag(natom2))
      allocate(ktp(natom2))
      allocate(kbo(natom2,32))

c  make up new versions of itp and ibo

      do i=1,nafrag
       ktp(i)=0
      do j=1,nafrag
c      if(ilink(i,j).gt.0)then
       if(ilink(i,j).ne.0)then
        ktp(i)=ktp(i)+1
         if(ktp(i).gt.32)then
          write(6,*)' ktp too big in findfragelim'
          stop
         endif
        kbo(i,ktp(i))=j
       endif
      enddo
      ktp(i)=ktp(i)-1
      enddo 

c  take powers of incidence matrix using
c  [a + a^2 + ... + a^(n-1)]_(i,j) <> 0 iff i,j path connected

      do i=1,nafrag
      do j=1,nafrag
         mt(i,j)=iabs(ilink(i,j))
         cmat(i,j)=iabs(ilink(i,j))
      enddo
      enddo

      do n=1,nafrag-1

         do i=1,nafrag
         do j=i,nafrag
            kt(i,j)=0
          if(ktp(i).lt.0)go to 200
         do m=1,ktp(i)+1
         kt(i,j)=kt(i,j)+mt(kbo(i,m),j)
         enddo
200    continue

         enddo
         enddo


         do i=1,nafrag
         do j=i,nafrag
            mt(i,j)=kt(i,j)
            if(mt(i,j).gt.0)mt(i,j)=1
            cmat(i,j)=cmat(i,j) + mt(i,j)
         enddo
         enddo

         do i=1,nafrag
         do j=i,nafrag
          mt(j,i)=mt(i,j)
          cmat(j,i)=cmat(i,j)
         enddo
         enddo

         do i=1,nafrag-1
         do j=i+1,nafrag
         if(cmat(i,j).eq.0)go to 41
         enddo
         enddo
         go to 42
41       continue

      enddo
42    continue
c now find fragments

      do i=1,nafrag
         ifrag(i)=0
      enddo

c modification 260711 for eliminatrion code
      newfrag=1
      do i=1,nafrag
       if(notthere(i).eq.1)then
        nfirst=i
        ifrag(i)=1
        go to 555
       endif
      enddo
555   continue

      if(nfirst.eq.nafrag)go to 556

      do j=1,nafrag
         if (cmat(nfirst,j).gt.0) ifrag(j)=1
      enddo

      do 10 i=1,nafrag
         if (ifrag(i).gt.0) goto 10
         if(notthere(i).eq.0)go to 10
         newfrag=newfrag+1
         ifrag(i)=newfrag
         do j=1,nafrag
            if (cmat(i,j).gt.0) ifrag(j)=newfrag
         enddo
10    continue

556   continue

c  size of fragments: kf(j) is cardinality of jth fragment 

      do j=1,newfrag
         kf(j)=0
         nn=0
         do n=1,nafrag
            if (ifrag(n).eq.j) then 
               kf(j)=kf(j)+1 
               nn=nn+1
               map(j,nn)=n
            endif
         enddo
      enddo

      deallocate(cmat)
      deallocate(kt)
      deallocate(mt)
      deallocate(ifrag)
      deallocate(ktp)
      deallocate(kbo)


      return
      end
