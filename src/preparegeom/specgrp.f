      subroutine specgrp(na,natom,jmat,neigh)
      implicit double precision(a-h,o-z)

      dimension jmat(natom,natom),neigh(natom)

      neis=0
      neistore=0
      do n=1,natom
       neigh(n)=0
      enddo
      neigh(na)=1

1     continue
      do n=1,natom
       if(neigh(n).ne.0)then

       do m=1,natom
        if(jmat(n,m).ne.0)then
         neigh(m)=1
        endif
       enddo

       endif
      enddo

      neis=0
      do m=1,natom
      if(neigh(m).eq.1)neis=neis+1
      enddo
c     write(6,*)neis

      if(neis.eq.neistore)then
       go to 2
      else
       neistore=neis
       neis=0
       go to 1
      endif

2     continue

      return
      end
