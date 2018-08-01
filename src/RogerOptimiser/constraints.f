      subroutine constraints(natom,coord,cvalue,gc,hc,jobtype)
      implicit double precision(a-h,o-z)

      dimension gc(3,natom),hc(3*natom,3*natom),coord(3,natom)

      integer, allocatable  :: mb(:),nb(:),ia(:),ja(:),ka(:)
      integer, allocatable  :: id(:),jd(:),kd(:),ld(:)
      real*8, allocatable   :: r0(:),a0(:),d0(:)

      dimension dr(2,3),d2r(2,3,2,3)
      dimension da(3,3),d2a(3,3,3,3)
      dimension jb(2),kb(2),rb0(2),drb(2)
      dimension dr1(2,3),d2r1(2,3,2,3)
      dimension dcd(4,3),dsd(4,3),d2cd(4,4,3,3),d2sd(4,4,3,3)

c mb and nb are the atom numbers at the ends of constrained bonds
c ia,ja,ka are the atoms in constrained angles ia..ja..ka

      Bohr=1.88972598860
      n3=3*natom

c always evaluate the second derivatives
      jobtype=2

      do n=1,natom
      do k=1,3
      gc(k,n)=0.d0
      enddo
      enddo
      do i=1,n3
      do j=1,n3
      hc(i,j)=0.d0
      enddo
      enddo
      cvalue=0.d0
      nbonds=0
      nangs=0

      pi=2.d0*acos(0.d0)

c read the constraints
c if a bond angle is constrained, so must the two bonds be constrained

      open(unit=1,file='IN_CONSTRAINTS',status='unknown')
      read(1,*,end=100)
      read(1,*)nbonds
      if(nbonds.gt.0)then
       allocate(mb(nbonds))
       allocate(nb(nbonds))
       allocate(r0(nbonds))
      endif
      if(nbonds.gt.0)then
       read(1,*)
       write(6,*)' The following bonds are constrained'
       do n=1,nbonds
        read(1,*)mb(n),nb(n),r0(n)
        write(6,*)mb(n),nb(n),r0(n)
c the coordinates are in Bohr, so convert the bond length to Bohr

        r0(n)=r0(n)*Bohr

       enddo
      endif

      read(1,*,end=101)
      read(1,*,end=101)nangs
      if(nangs.gt.0)then
       allocate(ia(nangs))
       allocate(ja(nangs))
       allocate(ka(nangs))
       allocate(a0(nangs))
      read(1,*)
      write(6,*)' The following angles are constrained'
      do n=1,nangs
       read(1,*)ia(n),ja(n),ka(n),a0(n)
       write(6,*)ia(n),ja(n),ka(n),a0(n)
c convert angles to cosines
       a0(n)=cos(pi*a0(n)/180.d0)
      enddo
      endif
      go to 102
101   nangs=0
102   continue

      ndiheds=0
      read(1,*,end=103)
      read(1,*)ndiheds
      if(ndiheds.gt.0)then
       allocate(id(ndiheds))
       allocate(jd(ndiheds))
       allocate(kd(ndiheds))
       allocate(ld(ndiheds))
       allocate(d0(ndiheds))
       read(1,*)
       write(6,*)' The following dihedrals are constrained'
       do n=1,ndiheds
        read(1,*)id(n),jd(n),kd(n),ld(n),d0(n)
        write(6,*)id(n),jd(n),kd(n),ld(n),d0(n)
        d0(n)=d0(n)*pi/180.d0
       enddo
      endif
103   continue
      

c finished reading
      close(unit=1)

c  fix the energy penalties
      penR=1.d1
      penA=1.d1

c loop over the constrained  bonds, calculating grads and hessians
c and allocating to the correct atoms

      if(nbonds.eq.0)go to 701

      write(6,*)' Bond constraints'
      do n=1,nbonds
201   format(2i10,2f13.6)
      call dbond(coord(:,mb(n)),coord(:,nb(n)),r,dr,d2r,jobtype)
      write(6,201) mb(n),nb(n),r0(n)/Bohr,r/Bohr

      cvalue=cvalue+0.5d0*penR*(r-r0(n))**2
      do k=1,3
       gc(k,mb(n))=gc(k,mb(n))+penR*(r-r0(n))*dr(1,k)
       gc(k,nb(n))=gc(k,nb(n))+penR*(r-r0(n))*dr(2,k)
      enddo

      if(jobtype.eq.2)then

      do n1=1,2
      do k1=1,3
       if(n1.eq.1)then
       j1=3*(mb(n)-1)+k1
       else
       j1=3*(nb(n)-1)+k1
       endif
      do n2=1,2
      do k2=1,3
       if(n2.eq.1)then
       j2=3*(mb(n)-1)+k2
       else
       j2=3*(nb(n)-1)+k2
       endif

       hc(j1,j2)=hc(j1,j2)+penR*dr(n1,k1)*dr(n2,k2)
     .           +penR*(r-r0(n))*d2r(n1,k1,n2,k2)

      enddo
      enddo
      enddo
      enddo

c end jobtype if
      endif

c finished with bonds
      enddo
701   continue

c now the angles
      if(nangs.eq.0)go to 702
      write(6,*)' Angle constraints'
      do n=1,nangs
202   format(3i10,2f13.6)
      call dangle(coord(:,ia(n)),coord(:,ja(n)),coord(:,ka(n)),
     .            a,da,d2a,jobtype)
      write(6,202)ia(n),ja(n),ka(n),a0(n),a

      cvalue=cvalue+0.5d0*penA*(a-a0(n))**2

      do k=1,3
       gc(k,ia(n))=gc(k,ia(n))+penA*(a-a0(n))*da(1,k)
       gc(k,ja(n))=gc(k,ja(n))+penA*(a-a0(n))*da(2,k)
       gc(k,ka(n))=gc(k,ka(n))+penA*(a-a0(n))*da(3,k)
      enddo

      if(jobtype.eq.2)then

      do n1=1,3
      do k1=1,3
       if(n1.eq.1)then
        j1=3*(ia(n)-1)+k1
       else
        if(n1.eq.2)then
         j1=3*(ja(n)-1)+k1
        else
         j1=3*(ka(n)-1)+k1
        endif
       endif
      do n2=1,3
      do k2=1,3
       if(n2.eq.1)then
        j2=3*(ia(n)-1)+k2
       else
        if(n2.eq.2)then
         j2=3*(ja(n)-1)+k2
        else
         j2=3*(ka(n)-1)+k2
        endif
       endif

       hc(j1,j2)=hc(j1,j2)+penA*da(n1,k1)*da(n2,k2)
     .           +penA*(a-a0(n))*d2a(n1,k1,n2,k2)

      enddo
      enddo
      enddo
      enddo

c end jobtype if
      endif

c finished with angles
      enddo
702   continue

      if(ndiheds.eq.0)go to 100
      write(6,*)' Dihedral constraints'
      do n=1,ndiheds
       call ddihed(coord(:,id(n)),coord(:,jd(n)),coord(:,kd(n)),
     .     coord(:,ld(n)),cd,dcd,d2cd,sd,dsd,d2sd,jobtype)

       cvalue=cvalue+0.25d0*penA*(cd-cos(d0(n)))**2
     .              +0.25d0*penA*(sd-sin(d0(n)))**2

      do k=1,3
       gc(k,id(n))=gc(k,id(n))+0.5d0*penA*(cd-cos(d0(n)))*dcd(1,k)
     .                        +0.5d0*penA*(sd-sin(d0(n)))*dsd(1,k)
       gc(k,jd(n))=gc(k,jd(n))+0.5d0*penA*(cd-cos(d0(n)))*dcd(2,k)
     .                        +0.5d0*penA*(sd-sin(d0(n)))*dsd(2,k)
       gc(k,kd(n))=gc(k,kd(n))+0.5d0*penA*(cd-cos(d0(n)))*dcd(3,k)
     .                        +0.5d0*penA*(sd-sin(d0(n)))*dsd(3,k)
       gc(k,ld(n))=gc(k,ld(n))+0.5d0*penA*(cd-cos(d0(n)))*dcd(4,k)
     .                        +0.5d0*penA*(sd-sin(d0(n)))*dsd(4,k)
      enddo

      if(jobtype.eq.2)then

       do n1=1,4
       do k1=1,3
        if(n1.eq.1)j1=3*(id(n)-1)+k1
        if(n1.eq.2)j1=3*(jd(n)-1)+k1
        if(n1.eq.3)j1=3*(kd(n)-1)+k1
        if(n1.eq.4)j1=3*(ld(n)-1)+k1
       do n2=1,4
       do k2=1,3
        if(n2.eq.1)j2=3*(id(n)-1)+k2
        if(n2.eq.2)j2=3*(jd(n)-1)+k2
        if(n2.eq.3)j2=3*(kd(n)-1)+k2
        if(n2.eq.4)j2=3*(ld(n)-1)+k2

        hc(j1,j2)=hc(j1,j2)+0.5d0*penA*dcd(n1,k1)*dcd(n2,k2)
     .                     +0.5d0*penA*dsd(n1,k1)*dsd(n2,k2)
     .                   +0.5d0*pena*(cd-cos(d0(n)))*d2cd(n1,n2,k1,k2)
     .                   +0.5d0*pena*(sd-sin(d0(n)))*d2cd(n1,n2,k1,k2)
       enddo
       enddo
       enddo
       enddo

      endif
      enddo

100   continue
      return
      end


