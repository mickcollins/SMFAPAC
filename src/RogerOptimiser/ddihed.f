      subroutine ddihed(x1,x2,x3,x4,cd,dcd,d2cd,sd,dsd,d2sd,nflag)
      implicit double precision(a-h,o-z)

c checked by finite difference 151217

      dimension x1(3),x2(3),x3(3),x4(3)
      dimension dcd(4,3),dsd(4,3)
      dimension d2cd(4,4,3,3),d2sd(4,4,3,3)

      dimension r12(3),r23(3),r34(3)
      dimension a(3),b(3),c(3)

      dimension dca(3),dcb(3),dsa(3),dsb(3),dsc(3)
      dimension d2a(3,4,4,3,3),d2b(3,4,4,3,3)
      dimension da(3,4,3),db(3,4,3),dc(3,4,3)
      dimension acb(3),acc(3),ccb(3)
      dimension dccbdb(3,3),dccbdc(3,3),daccda(3,3)
      dimension daccdc(3,3),dacbda(3,3),dacbdb(3,3)

      dimension d2caa(3,3),d2cab(3,3),d2cbb(3,3)
      dimension d2saa(3,3),d2sab(3,3),d2sac(3,3)
      dimension d2sbb(3,3),d2sbc(3,3),d2scc(3,3)
      

      do k=1,3
       r12(k)=x2(k)-x1(k)
       r23(k)=x3(k)-x2(k)
       r34(k)=x4(k)-x3(k)
       c(k)=r23(k)
      enddo
      a(1)=r34(2)*r23(3)-r34(3)*r23(2)
      a(2)=r34(3)*r23(1)-r34(1)*r23(3)
      a(3)=r34(1)*r23(2)-r34(2)*r23(1)
      b(1)=r23(2)*r12(3)-r23(3)*r12(2)
      b(2)=r23(3)*r12(1)-r23(1)*r12(3)
      b(3)=r23(1)*r12(2)-r23(2)*r12(1)
      
      da(1,1,1)=0.d0
      da(1,1,2)=0.d0
      da(1,1,3)=0.d0
      da(1,2,1)=0.d0
      da(1,2,2)=r34(3)
      da(1,2,3)=-r34(2)
      da(1,3,1)=0.d0
      da(1,3,2)=-r23(3)-r34(3)
      da(1,3,3)=r34(2)+r23(2)
      da(1,4,1)=0.d0
      da(1,4,2)=r23(3)
      da(1,4,3)=-r23(2)
      
      da(2,1,1)=0.d0
      da(2,1,2)=0.d0
      da(2,1,3)=0.d0
      da(2,2,1)=-r34(3)
      da(2,2,2)=0.d0
      da(2,2,3)=r34(1)
      da(2,3,1)=r34(3)+r23(3)
      da(2,3,2)=0.d0
      da(2,3,3)=-r34(1)-r23(1)
      da(2,4,1)=-r23(3)
      da(2,4,2)=0.d0
      da(2,4,3)=r23(1)

      da(3,1,1)=0.d0
      da(3,1,2)=0.d0
      da(3,1,3)=0.d0
      da(3,2,1)=r34(2)
      da(3,2,2)=-r34(1)
      da(3,2,3)=0.d0
      da(3,3,1)=-r23(2)-r34(2)
      da(3,3,2)=r34(1)+r23(1)
      da(3,3,3)=0.d0
      da(3,4,1)=r23(2)
      da(3,4,2)=-r23(1)
      da(3,4,3)=0.d0

      db(1,1,1)=0.d0
      db(1,1,2)=r23(3)
      db(1,1,3)=-r23(2)
      db(1,2,1)=0.d0
      db(1,2,2)=-r12(3)-r23(3)
      db(1,2,3)=r23(2)+r12(2)
      db(1,3,1)=0.d0
      db(1,3,2)=r12(3)
      db(1,3,3)=-r12(2)
      db(1,4,1)=0.d0
      db(1,4,2)=0.d0
      db(1,4,3)=0.d0

      db(2,1,1)=-r23(3)
      db(2,1,2)=0.d0
      db(2,1,3)=r23(1)
      db(2,2,1)=r23(3)+r12(3)
      db(2,2,2)=0.d0
      db(2,2,3)=-r12(1)-r23(1)
      db(2,3,1)=-r12(3)
      db(2,3,2)=0.d0
      db(2,3,3)=r12(1)
      db(2,4,1)=0.d0
      db(2,4,2)=0.d0
      db(2,4,3)=0.d0

      db(3,1,1)=r23(2)
      db(3,1,2)=-r23(1)
      db(3,1,3)=0.d0
      db(3,2,1)=-r23(2)-r12(2)
      db(3,2,2)=r23(1)+r12(1)
      db(3,2,3)=0.d0
      db(3,3,1)=r12(2)
      db(3,3,2)=-r12(1)
      db(3,3,3)=0.d0
      db(3,4,1)=0.d0
      db(3,4,2)=0.d0
      db(3,4,3)=0.d0

      dc(1,1,1)=0.d0
      dc(1,1,2)=0.d0
      dc(1,1,3)=0.d0
      dc(1,2,1)=-1.d0
      dc(1,2,2)=0.d0
      dc(1,2,3)=0.d0
      dc(1,3,1)=1.d0
      dc(1,3,2)=0.d0
      dc(1,3,3)=0.d0
      dc(1,4,1)=0.d0
      dc(1,4,2)=0.d0
      dc(1,4,3)=0.d0

      dc(2,1,1)=0.d0
      dc(2,1,2)=0.d0
      dc(2,1,3)=0.d0
      dc(2,2,1)=0.d0
      dc(2,2,2)=-1.d0
      dc(2,2,3)=0.d0
      dc(2,3,1)=0.d0
      dc(2,3,2)=1.d0
      dc(2,3,3)=0.d0
      dc(2,4,1)=0.d0
      dc(2,4,2)=0.d0
      dc(2,4,3)=0.d0

      dc(3,1,1)=0.d0
      dc(3,1,2)=0.d0
      dc(3,1,3)=0.d0
      dc(3,2,1)=0.d0
      dc(3,2,2)=0.d0
      dc(3,2,3)=-1.d0
      dc(3,3,1)=0.d0
      dc(3,3,2)=0.d0
      dc(3,3,3)=1.d0
      dc(3,4,1)=0.d0
      dc(3,4,2)=0.d0
      dc(3,4,3)=0.d0

      ccb(1)=c(2)*b(3)-c(3)*b(2)
      ccb(2)=c(3)*b(1)-c(1)*b(3)
      ccb(3)=c(1)*b(2)-c(2)*b(1)

      acc(1)=a(2)*c(3)-a(3)*c(2)
      acc(2)=a(3)*c(1)-a(1)*c(3)
      acc(3)=a(1)*c(2)-a(2)*c(1)

      acb(1)=a(2)*b(3)-a(3)*b(2)
      acb(2)=a(3)*b(1)-a(1)*b(3)
      acb(3)=a(1)*b(2)-a(2)*b(1)

      sum1=0.d0
      sum2=0.d0
      sum3=0.d0
      sum4=0.d0
      sum5=0.d0
      do k=1,3
       sum1=sum1+a(k)*a(k)
       sum2=sum2+b(k)*b(k)
       sum3=sum3+c(k)*c(k)
       sum4=sum4+a(k)*b(k)
       sum5=sum5+c(k)*acb(k)
      enddo
      sum1=sqrt(sum1)
      sum2=sqrt(sum2)
      sum3=sqrt(sum3)
      sumab=sum1*sum2
      cd=sum4/sumab
      sumabc=sumab*sum3
      sd=-sum5/sumabc

      do k=1,3
       dca(k)=b(k)/sumab-cd*a(k)/sum1**2
       dcb(k)=a(k)/sumab-cd*b(k)/sum2**2
       dsa(k)=ccb(k)/sumabc - sd*a(k)/sum1**2
       dsb(k)=acc(k)/sumabc - sd*b(k)/sum2**2
       dsc(k)=-acb(k)/sumabc - sd*c(k)/sum3**2
      enddo

      do n=1,4
      do k=1,3
       dcd(n,k)=dca(1)*da(1,n,k)+dca(2)*da(2,n,k)+dca(3)*da(3,n,k)
     .         +dcb(1)*db(1,n,k)+dcb(2)*db(2,n,k)+dcb(3)*db(3,n,k)
       dsd(n,k)=dsa(1)*da(1,n,k)+dsa(2)*da(2,n,k)+dsa(3)*da(3,n,k)
     .         +dsb(1)*db(1,n,k)+dsb(2)*db(2,n,k)+dsb(3)*db(3,n,k)
     .         +dsc(1)*dc(1,n,k)+dsc(2)*dc(2,n,k)+dsc(3)*dc(3,n,k)
      enddo
      enddo

      if(nflag.lt.2)return


      do j=1,3
      do n1=1,4
      do k1=1,3
      do n2=1,4
      do k2=1,3
       d2a(j,n1,n2,k1,k2)=0.d0
       d2b(j,n1,n2,k1,k2)=0.d0
      enddo
      enddo
      enddo
      enddo
      enddo

      d2a(1,2,3,2,3)=-1.d0
      d2a(1,2,4,2,3)=1.d0
      d2a(1,2,3,3,2)=1.d0
      d2a(1,2,4,3,2)=-1.d0
      d2a(1,3,2,2,3)=1.d0
      d2a(1,3,4,2,3)=-1.d0
      d2a(1,3,2,3,2)=-1.d0
      d2a(1,3,4,3,2)=1.d0
      d2a(1,4,2,2,3)=-1.d0
      d2a(1,4,3,2,3)=1.d0
      d2a(1,4,2,3,2)=1.d0
      d2a(1,4,3,3,2)=-1.d0

      d2a(2,2,3,1,3)=1.d0
      d2a(2,2,4,1,3)=-1.d0
      d2a(2,2,3,3,1)=-1.d0
      d2a(2,2,4,3,1)=1.d0
      d2a(2,3,2,1,3)=-1.d0
      d2a(2,3,4,1,3)=1.d0
      d2a(2,3,2,3,1)=1.d0
      d2a(2,3,4,3,1)=-1.d0
      d2a(2,4,2,1,3)=1.d0
      d2a(2,4,3,1,3)=-1.d0
      d2a(2,4,2,3,1)=-1.d0
      d2a(2,4,3,3,1)=1.d0

      d2a(3,2,3,1,2)=-1.d0
      d2a(3,2,4,1,2)=1.d0
      d2a(3,2,3,2,1)=1.d0
      d2a(3,2,4,2,1)=-1.d0
      d2a(3,3,2,1,2)=1.d0
      d2a(3,3,4,1,2)=-1.d0
      d2a(3,3,2,2,1)=-1.d0
      d2a(3,3,4,2,1)=1.d0 
      d2a(3,4,2,1,2)=-1.d0
      d2a(3,4,3,1,2)=1.d0
      d2a(3,4,2,2,1)=1.d0
      d2a(3,4,3,2,1)=-1.d0


      d2b(1,1,2,2,3)=-1.d0
      d2b(1,1,3,2,3)=1.d0
      d2b(1,1,2,3,2)=1.d0
      d2b(1,1,3,3,2)=-1.d0
      d2b(1,2,1,2,3)=1.d0
      d2b(1,2,3,2,3)=-1.d0
      d2b(1,2,1,3,2)=-1.d0
      d2b(1,2,3,3,2)=1.d0
      d2b(1,3,1,2,3)=-1.d0
      d2b(1,3,2,2,3)=1.d0
      d2b(1,3,1,3,2)=1.d0
      d2b(1,3,2,3,2)=-1.d0

      d2b(2,1,2,1,3)=1.d0
      d2b(2,1,3,1,3)=-1.d0
      d2b(2,1,2,3,1)=-1.d0
      d2b(2,1,3,3,1)=1.d0
      d2b(2,2,1,1,3)=-1.d0
      d2b(2,2,3,1,3)=1.d0
      d2b(2,2,1,3,1)=1.d0
      d2b(2,2,3,3,1)=-1.d0
      d2b(2,3,1,1,3)=1.d0
      d2b(2,3,2,1,3)=-1.d0
      d2b(2,3,1,3,1)=-1.d0
      d2b(2,3,2,3,1)=1.d0


      d2b(3,1,2,1,2)=-1.d0
      d2b(3,1,3,1,2)=1.d0
      d2b(3,1,2,2,1)=1.d0
      d2b(3,1,3,2,1)=-1.d0
      d2b(3,2,1,1,2)=1.d0
      d2b(3,2,3,1,2)=-1.d0
      d2b(3,2,1,2,1)=-1.d0
      d2b(3,2,3,2,1)=1.d0
      d2b(3,3,1,1,2)=-1.d0
      d2b(3,3,2,1,2)=1.d0
      d2b(3,3,1,2,1)=1.d0
      d2b(3,3,2,2,1)=-1.d0

      dccbdb(1,1)=0.d0
      dccbdb(1,2)=-c(3)
      dccbdb(1,3)=c(2)
      dccbdb(2,1)=c(3)
      dccbdb(2,2)=0.d0
      dccbdb(2,3)=-c(1)
      dccbdb(3,1)=-c(2)
      dccbdb(3,2)=c(1)
      dccbdb(3,3)=0.d0

      dccbdc(1,1)=0.d0
      dccbdc(1,2)=b(3)
      dccbdc(1,3)=-b(2)
      dccbdc(2,1)=-b(3)
      dccbdc(2,2)=0.d0
      dccbdc(2,3)=b(1)
      dccbdc(3,1)=b(2)
      dccbdc(3,2)=-b(1)
      dccbdc(3,3)=0.d0

      daccda(1,1)=0.d0
      daccda(1,2)=c(3)
      daccda(1,3)=-c(2)
      daccda(2,1)=-c(3)
      daccda(2,2)=0.d0
      daccda(2,3)=c(1)
      daccda(3,1)=c(2)
      daccda(3,2)=-c(1)
      daccda(3,3)=0.d0

      daccdc(1,1)=0.d0
      daccdc(1,2)=-a(3)
      daccdc(1,3)=a(2)
      daccdc(2,1)=a(3)
      daccdc(2,2)=0.d0
      daccdc(2,3)=-a(1)
      daccdc(3,1)=-a(2)
      daccdc(3,2)=a(1)
      daccdc(3,3)=0.d0

      dacbda(1,1)=0.d0
      dacbda(1,2)=b(3)
      dacbda(1,3)=-b(2)
      dacbda(2,1)=-b(3)
      dacbda(2,2)=0.d0
      dacbda(2,3)=b(1)
      dacbda(3,1)=b(2)
      dacbda(3,2)=-b(1)
      dacbda(3,3)=0.d0
      
      dacbdb(1,1)=0.d0
      dacbdb(1,2)=-a(3)
      dacbdb(1,3)=a(2)
      dacbdb(2,1)=a(3)
      dacbdb(2,2)=0.d0
      dacbdb(2,3)=-a(1)
      dacbdb(3,1)=-a(2)
      dacbdb(3,2)=a(1)
      dacbdb(3,3)=0.d0

      do k1=1,3
      do k2=1,3
       d2caa(k1,k2)=-a(k2)*b(k1)/(sumab*sum1**2)-dca(k2)*a(k1)/sum1**2
     .              + 2.d0*cd*a(k1)*a(k2)/sum1**4
       if(k1.eq.k2)d2caa(k1,k2)=d2caa(k1,k2)-cd/sum1**2
       d2cab(k1,k2)=-b(k1)*b(k2)/(sumab*sum2**2)-dcb(k2)*a(k1)/sum1**2
       if(k1.eq.k2)d2cab(k1,k2)=d2cab(k1,k2)+1.d0/sumab
       d2cbb(k1,k2)=-a(k1)*b(k2)/(sumab*sum2**2)-dcb(k2)*b(k1)/sum2**2
     .              +2.d0*cd*b(k1)*b(k2)/sum2**4
       if(k1.eq.k2)d2cbb(k1,k2)=d2cbb(k1,k2)-cd/sum2**2
      enddo
      enddo

      do k1=1,3
      do k2=1,3
       d2saa(k1,k2)=-a(k2)*ccb(k1)/(sumabc*sum1**2)
     .              -dsa(k2)*a(k1)/sum1**2+2.d0*sd*a(k1)*a(k2)/sum1**4
       if(k1.eq.k2)d2saa(k1,k2)=d2saa(k1,k2)-sd/sum1**2
       d2sab(k1,k2)=dccbdb(k1,k2)/sumabc-ccb(k1)*b(k2)/(sumabc*sum2**2)
     .              -dsb(k2)*a(k1)/sum1**2
       d2sac(k1,k2)=dccbdc(k1,k2)/sumabc-ccb(k1)*c(k2)/(sumabc*sum3**2)
     .              -dsc(k2)*a(k1)/sum1**2
       d2sbb(k1,k2)=-acc(k1)*b(k2)/(sumabc*sum2**2)
     .              -dsb(k2)*b(k1)/sum2**2+2.d0*sd*b(k1)*b(k2)/sum2**4
       if(k1.eq.k2)d2sbb(k1,k2)=d2sbb(k1,k2)-sd/sum2**2
       d2sbc(k1,k2)=daccdc(k1,k2)/sumabc-acc(k1)*c(k2)/(sumabc*sum3**2)
     .              -dsc(k2)*b(k1)/sum2**2
       d2scc(k1,k2)=acb(k1)*c(k2)/(sumabc*sum3**2)
     .              -dsc(k2)*c(k1)/sum3**2+2.d0*sd*c(k1)*c(k2)/sum3**4
       if(k1.eq.k2)d2scc(k1,k2)=d2scc(k1,k2)-sd/sum3**2
      enddo
      enddo

      do n1=1,4
      do k1=1,3
      do n2=1,4
      do k2=1,3
       d2cd(n1,n2,k1,k2)=0.d0
       d2sd(n1,n2,k1,k2)=0.d0
       do j1=1,3
       do j2=1,3
        d2cd(n1,n2,k1,k2)=d2cd(n1,n2,k1,k2)+
     .                    d2caa(j1,j2)*da(j1,n1,k1)*da(j2,n2,k2)+
     .                    d2cab(j1,j2)*da(j1,n1,k1)*db(j2,n2,k2)+
     .                    d2cab(j2,j1)*da(j2,n2,k2)*db(j1,n1,k1)+
     .                    d2cbb(j1,j2)*db(j1,n1,k1)*db(j2,n2,k2)
        d2sd(n1,n2,k1,k2)=d2sd(n1,n2,k1,k2)+
     .                    d2saa(j1,j2)*da(j1,n1,k1)*da(j2,n2,k2)+
     .                    d2sab(j1,j2)*da(j1,n1,k1)*db(j2,n2,k2)+
     .                    d2sab(j2,j1)*da(j2,n2,k2)*db(j1,n1,k1)+
     .                    d2sac(j1,j2)*da(j1,n1,k1)*dc(j2,n2,k2)+
     .                    d2sac(j2,j1)*da(j2,n2,k2)*dc(j1,n1,k1)+
     .                    d2sbb(j1,j2)*db(j1,n1,k1)*db(j2,n2,k2)+
     .                    d2sbc(j1,j2)*db(j1,n1,k1)*dc(j2,n2,k2)+
     .                    d2sbc(j2,j1)*db(j2,n2,k2)*dc(j1,n1,k1)+
     .                    d2scc(j1,j2)*dc(j1,n1,k1)*dc(j2,n2,k2)
       enddo
       enddo
      enddo       
      enddo       
      enddo       
      enddo       


      do n1=1,4
      do k1=1,3
      do n2=1,4
      do k2=1,3
       do j1=1,3
        d2cd(n1,n2,k1,k2)=d2cd(n1,n2,k1,k2)+dca(j1)*d2a(j1,n1,n2,k1,k2)
     .                                     +dcb(j1)*d2b(j1,n1,n2,k1,k2)
        d2sd(n1,n2,k1,k2)=d2sd(n1,n2,k1,k2)+dsa(j1)*d2a(j1,n1,n2,k1,k2)
     .                                     +dsb(j1)*d2b(j1,n1,n2,k1,k2)
       enddo
      enddo
      enddo
      enddo
      enddo


      return
      end 
