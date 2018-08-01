      subroutine printdaltons(ngroups,nat0,nch0,lab0,c0)

      use molecule_adt
      use dalton

      implicit none

      integer ngroups,n,i,n1,k

      character*80 daltonname

      character*20 ca

      integer nat0(ngroups),nch0(ngroups)
      real*8 c0(ngroups,60,3),centre(ngroups,3)
      character*2 lab0(ngroups,60)

      type(molecule),allocatable    ::groups_mol(:)


      allocate(groups_mol(ngroups))

      do n=1,ngroups
       call init_mol(groups_mol(n))
       call new_mol(groups_mol(n),nat0(n))
      enddo

c find the group centres
      do n=1,ngroups
      do k=1,3
       centre(n,k)=0.d0
       do i=1,nat0(n)
        centre(n,k)=centre(n,k)+c0(n,i,k)
       enddo
       centre(n,k)=centre(n,k)/dble(float(nat0(n)))
      enddo
      enddo

      do n=1,ngroups
       do i=1,nat0(n)
        groups_mol(n)%atoms(i)%label=lab0(n,i)
        do k=1,3
        groups_mol(n)%atoms(i)%Coords(k)=c0(n,i,k)-centre(n,k)
        enddo
       enddo
      enddo

      do n=1,ngroups
       call filelabel(n,ca)
       n1=index(ca,' ')-1
       daltonname='grp.'//ca(1:n1)//'.mol'    
       call print_dalton(daltonname,groups_mol(n),nch0(n))
      enddo

      return
      end
