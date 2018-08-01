      subroutine AtNum(ndim,na,lab,num,sym,numb,nflag)
      implicit double precision(a-h,o-z)
      dimension numb(110)
      character*2 lab(ndim),sym(110)
      dimension num(ndim)

c nflag is an error flag
      nflag=0
c get the atomic number for each atom 
      do n=1,na
       do m=1,110
100     format(1x,a2,1x,a2)
        if(sym(m).eq.lab(n))then
         num(n)=numb(m)
         go to 1
        endif
       enddo
       write(*,*)' subroutine AtNum failed to match'
       nflag=1
1      continue
      enddo
      return
      end
