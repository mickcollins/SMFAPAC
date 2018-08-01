      subroutine match_mg(ndim,mg,numat,natstore,nffinal,
     .                    nbextra,match)
      implicit double precision(a-h,o-z)

      dimension numat(ndim),natstore(ndim,8)

c to check if group mg is already a L1 fragment,
c if not to create an extra L1 fragment for mg
c and return the location of this fragemnt in the L1 list

       if(mg.eq.0)then
        write(6,*)' mg = 0'
        stop
       endif

       match=0
       do j=1,nffinal+nbextra
        if(numat(j).eq.1.and.mg.eq.natstore(j,1))then
         match=j
         go to 114
        endif
       enddo
114    continue
c if no match, we must make an addition to the L1 list
       if(match.eq.0)then
        nbextra=nbextra+1
        numat(nffinal+nbextra)=1
        natstore(nffinal+nbextra,1)=mg
        match=nffinal+nbextra
       endif
       return
       end
