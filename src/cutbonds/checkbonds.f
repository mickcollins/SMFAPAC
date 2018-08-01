      subroutine checkbonds(nbonds,mb,nb,mult)
      implicit double precision(a-h,o-z)

      dimension mb(nbonds), nb(nbonds), mult(nbonds)

      dimension killit(nbonds)

c check to see if some bonds are repeated. This can happen because
c bonds where introduced in IN_EXTRABONDS that already existed on the
c basis of distance. An “extra bond” may have a mult value of -1
c so that it can be broken along a reaction path without caps being 
c introduced in the fragmentation.

      killit=0

      do n=1,nbonds-1
       if(killit(n).eq.1)go to 1
      do m=n+1,nbonds
       if(killit(m).eq.1)go to 2
       if(mb(n).eq.mb(m).and.nb(n).eq.nb(m))killit(n)=1
       if(mb(n).eq.nb(m).and.nb(n).eq.mb(m))killit(n)=1
2     enddo
1     enddo

c remove redundant bonds
      ic=0
      do n=1,nbonds
       if(killit(n).eq.1)go to 3
       ic=ic+1
       mb(ic)=mb(n)
       nb(ic)=nb(n)
       mult(ic)=mult(n)
3     enddo

      if(ic.lt.nbonds)then
      write(6,*)' Removed redundant bond definitions'
      nbonds=ic
      endif

      return
      end
