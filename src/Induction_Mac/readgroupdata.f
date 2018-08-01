      subroutine readgroupdata(ngroups,c0,nat0,alpha,cent)
      implicit double precision(a-h,o-z)

      dimension nat0(ngroups)
      dimension c0(ngroups,60,3),alpha(ngroups,3,3),cent(ngroups,3)

      character*6 ca
      character*16 ca1
      character*2 lab0

c read in the group polarizabilities 
      open(unit=1,file='Polarisabilities.out',status='old')
      do n=1,ngroups
       read(1,*)ii,alpha(n,1,1),alpha(n,2,1),alpha(n,2,2),
     .             alpha(n,3,1),alpha(n,3,2),alpha(n,3,3)
       alpha(n,1,2)=alpha(n,2,1)
       alpha(n,1,3)=alpha(n,3,1)
       alpha(n,2,3)=alpha(n,3,2)
      enddo
      close(unit=1)

c read in the Level 0 (group) coordinates
      do n=1,ngroups
       call filelabel(n,ca)
       ca1='Lev0_COORD'//ca
       open(unit=1,file=ca1,status='old')
       read(1,*)nch0
       m=1
1      continue
       read(1,*,end=2)lab0,(c0(n,m,k),k=1,3)
       m=m+1
       go to 1
2      continue
       nat0(n)=m-1
       close(unit=1)
      enddo
100   format(a2,3f13.6)

c calculate the centre of each group
      do n=1,ngroups
      do k=1,3
       cent(n,k)=0.d0
       do j=1,nat0(n)
        cent(n,k)=cent(n,k)+c0(n,j,k)/dble(float(nat0(n)))
       enddo
      enddo
      enddo

c convert to Bohr
      Bohr=1.88972598860
      do n=1,ngroups
      do k=1,3
       cent(n,k)=cent(n,k)*Bohr
      enddo
      enddo

      return
      end
