      program fixnearcharge

c this version of adjust only shifts the charge from the cap to the
c capped atom

      implicit double precision(a-h,o-z)

      real*8, allocatable  :: c(:,:),charge(:),x(:,:),chg(:)
      real*8, allocatable  :: z(:,:),weight(:,:)

      integer, allocatable  :: nattach(:),match(:),ngp(:,:),ncentre(:)
      integer, allocatable  :: numorig(:),norig(:),natominsite(:,:)

      character*2, allocatable   :: lab(:)


      dimension xcent(3)

      character*24 file1, file2,file3,fileout

      character*1 I0(0:10),ca1,ca2,ca3,ca4
      character*15 ca
      character*21 canum
      character*27 ga

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

c read a "charge.XX.coord.npa" file from unit 5

      read(5,*)file1
      read(5,*)natom

c allcocate arrays
      allocate(c(natom,3))
      allocate(charge(natom))
      allocate(match(natom))
      allocate(numorig(natom))
      allocate(lab(natom))

      do n=1,natom
       read(5,*)(c(n,k),k=1,3),charge(n)
      enddo

      nend=index(file1," ")-1
      ga=file1(1:nend)//'_atoms'

      open(unit=1,file=file1,status='old')
      open(unit=2,file=ga,status='old')
      read(1,*)
      do n=1,natom
c     read(1,101)lab(n)
      read(1,666)lab(n)
      read(2,*)numorig(n)
      enddo
      close(unit=1)
      close(unit=2)
101   format(a2)
666   format(a2,3f13.6,i10)

c file1 is a string like "charge.xx.coord". We will need the
c "charge.xx" part

      l1=len(file1)
      ns=index(file1,' ')
      file2=file1(1:ns-7)
      ns1=index(file2,' ')
      fileout=file2(1:ns1-1)//'.chg'
      file3=file1(1:ns-1)//'.cap'
      open(unit=1,file=file3,status='old')

      read(1,*)ncore,ncaps,numfam
c     if(ncaps.eq.0)then
c      write(6,*)' ncaps = 0 in a charge file'
c      stop
c     endif

      if(ncaps.gt.0)then

c allocate memory
      allocate(x(ncaps,3))
      allocate(chg(ncaps))
      allocate(z(ncaps,3))
      allocate(nattach(ncaps))
      allocate(ngp(ncaps,12))
      allocate(ncentre(ncaps))
      allocate(norig(ncaps))
      allocate(natominsite(ncaps,ncaps))
      allocate(weight(ncaps,ncaps))


      do n=1,ncaps
       read(1,*)(x(n,k),k=1,3),nattach(n)
      enddo
      endif
      close(unit=1)

c generate the group file name
      k=numfam
      if(k.le.9)then
       ca='charge.'//I0(k)//'.grp'
      endif
      if(k.gt.9.and.k.le.99)then
      k1=k/10
      k2=k-k1*10
      ca1=I0(k1)
      ca2=I0(k2)
      ca='charge.'//ca1//ca2//'.grp'

      endif
      if(k.gt.99.and.k.le.999)then
      k1=k/100
      k3=k/10 - INT(k1)*10
      k4=k - INT(k1)*100 - INT(k3)*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca='charge.'//ca1//ca2//ca3//'.grp'
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
      ca='charge.'//ca1//ca2//ca3//ca4//'.grp'
      endif
      open(unit=66,file=ca,status='unknown')

      nend=index(ca," ")-1
      canum=ca(1:nend)//'_atoms'
      open(unit=67,file=canum,status='unknown')

      do n=1,natom
       match(n)=0
      enddo
      nsites=0
c allow possible skip to the output
c     go to 888

c determine the number of distinct cap sites (nsites)

      if(ncaps.eq.0)then
       nsites=0
      else
       nsites=1
       ncentre(1)=nattach(1)
       if(ncaps.gt.1)then
       do i=2,ncaps
        do j=1,nsites
         if(ncentre(j).eq.nattach(i))go to 20
        enddo
        nsites=nsites+1
        ncentre(nsites)=nattach(i)
c end the i loop
20      enddo
       endif
      endif


      if(nsites.eq.0)go to 555

c build the sets of H cap atoms near each cap site
      tol=1.8d0
      totalsitecharge=0.d0
      do j=1,nsites

       ngp(j,1)=ncentre(j)
       match(ncentre(j))=1
       m=1

c determine if a H atom is a cap on ncentre(j)

       do n=1,natom
        if(lab(n).ne.'H ')go to 21

        do nn=1,ncaps
         if(nattach(nn).ne.ncentre(j))go to 212
         sum=0.d0
         do k=1,3
          sum=sum+(c(n,k)-x(nn,k))**2
         enddo
         if(sum.lt.1.d-2)then
          m=m+1
          ngp(j,m)=n
          match(n)=1
          go to 21
         endif
212     enddo

21     enddo



c calculate the charge for this site
c and the original atom numbers and weights
       norig(j)=0
       chg(j)=0.d0
       do n=1,m
        chg(j)=chg(j)+charge(ngp(j,n))
        totalsitecharge=totalsitecharge+charge(ngp(j,n))
         if(numorig(ngp(j,n)).gt.0)then
          norig(j)=norig(j)+1
          natominsite(j,norig(j))=numorig(ngp(j,n))
         endif
       enddo

       do k=1,norig(j)
        weight(j,k)=1.d0/dble(float(norig(j)))
       enddo

c calculate the mean position, excluding caps
       do k=1,3
        z(j,k)=0.d0
       enddo
       ic=0
       do n=1,m
c exclude caps

c temp 120213 include caps
c      go to 999

       do i=1,ncaps
        sum=0.d0
        do k=1,3
         sum=sum+(c(ngp(j,n),k)-x(i,k))**2
        enddo
        if(sum.lt.0.01d0)go to 22
       enddo

999   continue
       ic=ic+1
       do k=1,3
        z(j,k)=z(j,k)+c(ngp(j,n),k)
       enddo
c end the n loop
22    enddo
      do k=1,3
       z(j,k)=z(j,k)/dble(float(ic))
      enddo

c end the loop over sites
      enddo

555   continue

      

888   continue

c calculate the total number of atoms not involved in a site
       ntot=0
       do n=1,natom
        if(match(n).eq.0)ntot=ntot+1
       enddo

       open(unit=2,file=fileout,status='unknown')
       do n=1,natom
       if(match(n).eq.0)then
        write(2,100)(c(n,k),k=1,3),charge(n)
        write(66,100)(c(n,k),k=1,3),charge(n)
        write(67,667)1,numorig(n),1.d0
       endif
       enddo
667    format(i4,i10,f10.6)
c      go to 777

      if(nsites.gt.0)then


c output sites
       do n=1,nsites
       write(2,100)(z(n,k),k=1,3),chg(n)
       write(66,100)(z(n,k),k=1,3),chg(n)
       write(67,668)norig(n),(natominsite(n,k),weight(n,k),k=1,norig(n))
       enddo
668    format(i4,6(i6,f10.6))
      endif
777   continue

100    format(4f13.7)
       close(unit=2)
       close(unit=66)
       
       
      end
