      subroutine elements(lab)
      implicit none

c converts elemental symbols into the coorect
c upper followed by lower case

      character*2 lab
      character*1 lab1,lab2

      integer L1,L2, capa,capz,lowa,lowz

      capa=iachar('A')
      capz=iachar('Z')

      lowa=iachar('a')
      lowz=iachar('z')

      lab1=lab(1:1)
      lab2=lab(2:2)

      if(ichar(lab1).ge.lowa.and.ichar(lab1).le.lowz)then
       lab1=char(ichar(lab1)-lowa+capa)
      endif

      if(ichar(lab2).ge.capa.and.ichar(lab2).le.capz)then
       lab2=char(ichar(lab2)+lowa-capa)
      endif

      lab=lab1//lab2

      return
      end
