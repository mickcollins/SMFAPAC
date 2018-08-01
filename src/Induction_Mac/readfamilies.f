      subroutine readfamilies(ngroups,nfam,ifam)
      implicit double precision(a-h,o-z)

      dimension nfam(ngroups),ifam(ngroups,60)

c need the families
      open(unit=4,file='families.out',status='old')
      read(4,*)
      read(4,*)
      read(4,*)
      do n=1,ngroups
       read(4,*)nfam(n)
c nfam is the number of atoms in the group
       read(4,*)(ifam(n,k),k=1,nfam(n))
c ifam is the real atom number
      enddo
      close(unit=4)
      return
      end
