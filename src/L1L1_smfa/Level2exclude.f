      subroutine Level2exclude(nfextra,ndim,k1,k2,numat,natstore,nflag)
      implicit double precision(a-h,o-z)

      dimension numat(ndim), natstore(ndim,8)
      dimension k1store(100),k2store(100)

      nflag=0

c find fragments that contain an atom from fragments k1 and k2

      nk1=0
      nk2=0

      do n1=1,nfextra

       do j1=1,numat(n1)
        do m1=1,numat(k1)
         if(natstore(n1,j1).eq.natstore(k1,m1))then
          nk1=nk1+1
          k1store(nk1)=n1
          go to 10
         endif
        enddo
       enddo
10     continue

       do j1=1,numat(n1)
        do m2=1,numat(k2)
         if(natstore(n1,j1).eq.natstore(k2,m2))then
          nk2=nk2+1
          k2store(nk2)=n1
          go to 20
         endif
        enddo
       enddo
20     continue

      enddo


c see if these two sets  of fragments have a member in common

      do n1=1,nk1
      do n2=1,nk2

        do n3=1,nfextra
         match=0
         do j3=1,numat(n3)
          do m1=1,numat(k1store(n1))
           if(natstore(k1store(n1),m1).eq.natstore(n3,j3))match=1
          enddo
         enddo
         if(match.eq.1)then
         do j3=1,numat(n3)
          do m2=1,numat(k2store(n2))
           if(natstore(k2store(n2),m2).eq.natstore(n3,j3))match=2
          enddo
         enddo
         endif
         if(match.eq.2)go to 1
        enddo
        
      enddo
      enddo

1     continue
      if(match.eq.2)nflag=1

      return
      end
