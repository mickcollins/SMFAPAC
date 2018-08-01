      subroutine adjustsingle(natom,nbonds,bonds,mult)

      implicit double precision(a-h,o-z)

      dimension nbonds(natom)

      integer bonds(natom,20)

      dimension mult(natom,20)


      open(unit=1,file='IN_SINGLE',status='unknown')
      nsingle=0
      read(1,*,end=2)
      read(1,*)nsingle
      if(nsingle.eq.0)go to 2
      read(1,*)
      do n=1,nsingle
       read(1,*)n1,n2
       if(nbonds(n1).eq.0.or.nbonds(n2).eq.0)then
        write(6,*)' Single, rather than multiple, bonds were requested'
        write(6,*)' for atoms that are not bonded to any atoms.'
        if(nbonds(n1).eq.0)write(6,*)' Inappropriate atom number ',n1
        if(nbonds(n2).eq.0)write(6,*)' Inappropriate atom number ',n2
        write(6,*)' Please adjust the optional bonding input.'
        write(6,*)' Program Preparegeom stopped'
        stop
       endif
       match=0
       do i=1,nbonds(n1)
        if(bonds(n1,i).eq.n2)then
         mult(n1,i)=1
         match=1
        endif
       enddo
       if(match.eq.0)then
        write(6,*)' An incorrect pair of atoms were entered, for which'
        write(6,*)' a single bond was requested, rather than a multiple'
        write(6,*)' bond. However, a multiple bond would not normally'
        write(6,*)' be appropriate for these atoms. The atoms are:'
        write(6,*) n1,n2
        write(6,*)' Please adjust the optional bonding input.'
        write(6,*)' Program Preparegeom stopped'
        stop
       endif

       match=0
       do i=1,nbonds(n2)
        if(bonds(n2,i).eq.n1)then
         mult(n2,i)=1
         match=1
        endif
       enddo
       if(match.eq.0)then
        write(6,*)' An incorrect pair of atoms were entered, for which'
        write(6,*)' a single bond was requested, rather than a multiple'
        write(6,*)' bond. However, a multiple bond would not normally'
        write(6,*)' be appropriate for these atoms. The atoms are:'
        write(6,*) n1,n2
        write(6,*)' Program Preparegeom stopped'
        stop
       endif
      enddo
2     continue
      close(unit=1)
      return
      end
