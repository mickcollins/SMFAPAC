      subroutine atomicnumber(natom,lab,numa)
      implicit double precision(a-h,o-z)

      dimension numa(natom),numb(91)
      character*2 lab(natom),sym(91)


      open(unit=7,file='atomic.dat',status='old')

      do n=1,91
      read(7,*)numb(n)
      read(7,21)sym(n)
      read(7,*)am
      read(7,*)rad
      enddo
      close(unit=7)
21    format(a2)

      do n=1,natom
      match=0
      do m=1,91
      if(sym(m).eq.lab(n))then
       match=1
       numa(n)=numb(m)
      endif
      enddo
      if(match.eq.0)then
       write(6,*)' no match in atomicnumber'
       stop
      endif
      enddo

      return
      end
