      subroutine distmin(ndim,n1,n2,c1,c2,rad1,rad2,dmin)
      implicit double precision(a-h,o-z)
      dimension c1(ndim,3),c2(ndim,3),rad1(ndim),rad2(ndim)

c find the minimum distance (rel to VdW radius)
c betweem two groups

      dmin=1.d6
      do j1=1,n1
      do j2=1,n2

      d=0.d0
      do k=1,3
       d=d+(c1(j1,k)-c2(j2,k))**2
      enddo
      d=sqrt(d)/(rad1(j1)+rad2(j2))
      if(d.lt.dmin)dmin=d

      enddo
      enddo

      return
      end
