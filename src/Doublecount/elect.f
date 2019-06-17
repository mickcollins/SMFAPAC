      subroutine elect(natom1,c1,ch1,d1,qu1,oct1,hex1,
     .                 natom2,c2,ch2,d2,qu2,oct2,hex2,
     .                 elen,del1,del2,secd1,secd2,secd3,jobtype)

      implicit double precision(a-h,o-z)

c calculates the electrostatic interaction of two fragments

c if jobtype=1, part of the cartesian derivatives of the energy is
c evaluated: the part that depends on the "t" tensors

c  output is:
c       the electrostatic, e1
c       del1(n1,k1) - the derivative of el wrt coordinate k1 of atom n1
c       del2(n2,k2) - the derivative of el wrt coordinate k2 of atom n2


      dimension c1(natom1,3),ch1(natom1),d1(natom1,3),
     .          qu1(natom1,3,3),oct1(natom1,3,3,3),
     .          hex1(natom1,3,3,3,3)

      dimension c2(natom2,3),ch2(natom2),d2(natom2,3),
     .          qu2(natom2,3,3),oct2(natom2,3,3,3),
     .          hex2(natom2,3,3,3,3)

      dimension del1(natom1,3),del2(natom2,3)
      dimension gstore1(natom1,3),gstore2(natom2,3)
      dimension del01(3),del02(3),dedR(3,3)

      dimension secd1(natom1,3,natom1,3),secd2(natom1,3,natom2,3),
     .          secd3(natom2,3,natom2,3),saveg(natom1,natom2,3)


      dimension t1(3),t2(3,3),t3(3,3,3),t4(3,3,3,3),
     .          t5(3,3,3,3,3),t6(3,3,3,3,3,3)

      dimension dect0c(3),dect1da(3),dect1db(3),
     .          dect2qa(3),dect2qb(3),dect3oa(3),dect3ob(3),
     .          dect4ha(3),dect4hb(3),dedt2d(3),dedt3qa(3),dedt3qb(3),
     .          dedt4oa(3),dedt4ob(3),deqt4q(3)


      dimension tot(10),tot1(3),tot2(3),dd(3,3)

      eps=0.0005d0


c loop over the atoms
      elen=0.d0
      do n1=1,natom1
      do k=1,3
       del1(n1,k)=0.d0
       gstore1(n1,k)=0.d0
      enddo
      enddo
      do n1=1,natom2
      do k=1,3
       del2(n1,k)=0.d0
       gstore2(n1,k)=0.d0
      enddo
      enddo

      secd1=0.d0
      secd2=0.d0
      secd3=0.d0

      do n1=1,natom1
      do n2=1,natom2

      do i=1,10
       tot(i)=0.d0
      enddo

      call tfactorfast(c1(n1,:),c2(n2,:),t0,t1,t2,t3,t4,t5,t6)
c     call tfactor(c1(n1,:),c2(n2,:),t0,t1,t2,t3,t4,t5,t6)

      call ct1d(ch2(n2),d1(n1,:),t1,ect1da)
      call ct1d(ch1(n1),d2(n2,:),t1,ect1db)
      call ct2q(ch2(n2),qu1(n1,:,:),t2,ect2qa)
      call ct2q(ch1(n1),qu2(n2,:,:),t2,ect2qb)
c next 4 re-added
c     call ct3o(ch2(n2),oct1(n1,:,:,:),t3,ect3oa)
c     call ct3o(ch1(n1),oct2(n2,:,:,:),t3,ect3ob)
c     call ct4h(ch2(n2),hex1(n1,:,:,:,:),t4,ect4ha)
c     call ct4h(ch1(n1),hex2(n2,:,:,:,:),t4,ect4hb)
      call dt2d(d2(n2,:),d1(n1,:),t2,edt2d)
      call dt3q(d2(n2,:),qu1(n1,:,:),t3,edt3qa)
      call dt3q(d1(n1,:),qu2(n2,:,:),t3,edt3qb)
c next 2 re-added
c     call dt4o(d2(n2,:),oct1(n1,:,:,:),t4,edt4oa)
c     call dt4o(d1(n1,:),oct2(n2,:,:,:),t4,edt4ob)
      call qt4q(qu2(n2,:,:),qu1(n1,:,:),t4,eqt4q)

      tot(1)=ch1(n1)*ch2(n2)*t0
      tot(2)=ect1da-ect1db
      tot(3)=(ect2qa+ect2qb)/3.d0
      tot(4)=-edt2d
c     tot(5)=(ect3oa-ect3ob)/15.d0
      tot(6)=(-edt3qa+edt3qb)/3.d0
c     tot(7)=(ect4ha+ect4hb)/105.d0
c     tot(8)=-(edt4oa+edt4ob)/15.d0
      tot(9)=eqt4q/9.d0


c temp patch for testing
c     tot(2)=0.d0
c     tot(3)=0.d0
c     tot(1)=0.d0
c     tot(6)=0.d0
c     tot(9)=0.d0


      do i=1,9
      elen=elen+tot(i)
      enddo

c end the n1,n2 loops
      enddo
      enddo


      if(jobtype.eq.0)go to 1000

      saveg=0.d0
c gradients

      do n1=1,natom1
      do n2=1,natom2
      call tfactorfast(c1(n1,:),c2(n2,:),t0,t1,t2,t3,t4,t5,t6)

      do k=1,3
       dect0c(k)=0.d0
       dect1da(k)=0.d0
       dect1db(k)=0.d0
       dect2qa(k)=0.d0
       dect2qb(k)=0.d0
       dect3oa(k)=0.d0
       dect3ob(k)=0.d0
       dect4ha(k)=0.d0
       dect4hb(k)=0.d0
       dedt2d(k)=0.d0
       dedt3qa(k)=0.d0
       dedt3qb(k)=0.d0
       dedt4oa(k)=0.d0
       dedt4ob(k)=0.d0
       deqt4q(k)=0.d0
      enddo

      call ct0c_d(ch2(n2),ch1(n1),t1,dect0c)
      call ct1d_d(ch2(n2),d1(n1,:),t2,dect1da)
      call ct1d_d(ch1(n1),d2(n2,:),t2,dect1db)
      call ct2q_d(ch2(n2),qu1(n1,:,:),t3,dect2qa)
      call ct2q_d(ch1(n1),qu2(n2,:,:),t3,dect2qb)
      call ct3o_d(ch2(n2),oct1(n1,:,:,:),t4,dect3oa)
      call ct3o_d(ch1(n1),oct2(n2,:,:,:),t4,dect3ob)
      call ct4h_d(ch2(n2),hex1(n1,:,:,:,:),t5,dect4ha)
      call ct4h_d(ch1(n1),hex2(n2,:,:,:,:),t5,dect4hb)
      call dt2d_d(d2(n2,:),d1(n1,:),t3,dedt2d)
      call dt3q_d(d2(n2,:),qu1(n1,:,:),t4,dedt3qa)
      call dt3q_d(d1(n1,:),qu2(n2,:,:),t4,dedt3qb)
      call dt4o_d(d2(n2,:),oct1(n1,:,:,:),t5,dedt4oa)
      call dt4o_d(d1(n1,:),oct2(n2,:,:,:),t5,dedt4ob)
      call qt4q_d(qu2(n2,:,:),qu1(n1,:,:),t5,deqt4q)

c     go to 4444

      do k=1,3
       sume=dect0c(k)+dect1da(k)-dect1db(k)
     .  +(dect2qa(k)+dect2qb(k))/3.d0
     .  +(dect3oa(k)-dect3ob(k))/15.d0
     .  +(dect4ha(k)+dect4hb(k))/105.d0
     .  -dedt2d(k)
     .  +(-dedt3qa(k)+dedt3qb(k))/3.d0
     .  -(dedt4oa(k)+dedt4ob(k))/15.d0
     .  +deqt4q(k)/9.d0
       del1(n1,k)=del1(n1,k)+sume
       del2(n2,k)=del2(n2,k)-sume
       saveg(n1,n2,k)=sume
      enddo

4444  continue
c temp patch for testing
c     do k=1,3
c      sume=dect0c(k)
c      sume=dect1da(k)-dect1db(k)
c      sume=(dect2qa(k)+dect2qb(k))/3.d0
c      sume=-dedt2d(k)
c      del1(n1,k)=del1(n1,k)+sume
c      del2(n2,k)=del2(n2,k)-sume
c      saveg(n1,n2,k)=sume
c     enddo

c end the n1,n2 loops
      enddo
      enddo

      do n1=1,natom1
      do k=1,3
       gstore1(n1,k)=del1(n1,k)
      enddo
      enddo
      do n2=1,natom2
      do k=1,3
       gstore2(n2,k)=del2(n2,k)
      enddo
      enddo

      if(jobtype.eq.1)go to 1000

c here because jobtype=2

      do n1=1,natom1
      do n2=1,natom2

       do j=1,3
       c1(n1,j)=c1(n1,j)+eps

       call tfactorfast(c1(n1,:),c2(n2,:),t0,t1,t2,t3,t4,t5,t6)

      call ct0c_d(ch2(n2),ch1(n1),t1,dect0c)
      call ct1d_d(ch2(n2),d1(n1,:),t2,dect1da)
      call ct1d_d(ch1(n1),d2(n2,:),t2,dect1db)
      call ct2q_d(ch2(n2),qu1(n1,:,:),t3,dect2qa)
      call ct2q_d(ch1(n1),qu2(n2,:,:),t3,dect2qb)
      call ct3o_d(ch2(n2),oct1(n1,:,:,:),t4,dect3oa)
      call ct3o_d(ch1(n1),oct2(n2,:,:,:),t4,dect3ob)
      call ct4h_d(ch2(n2),hex1(n1,:,:,:,:),t5,dect4ha)
      call ct4h_d(ch1(n1),hex2(n2,:,:,:,:),t5,dect4hb)
      call dt2d_d(d2(n2,:),d1(n1,:),t3,dedt2d)
      call dt3q_d(d2(n2,:),qu1(n1,:,:),t4,dedt3qa)
      call dt3q_d(d1(n1,:),qu2(n2,:,:),t4,dedt3qb)
      call dt4o_d(d2(n2,:),oct1(n1,:,:,:),t5,dedt4oa)
      call dt4o_d(d1(n1,:),oct2(n2,:,:,:),t5,dedt4ob)
      call qt4q_d(qu2(n2,:,:),qu1(n1,:,:),t5,deqt4q)

c     go to 4444

      do k=1,3
       tot1(k)=dect0c(k)+dect1da(k)-dect1db(k)
     .  +(dect2qa(k)+dect2qb(k))/3.d0
     .  +(dect3oa(k)-dect3ob(k))/15.d0
     .  +(dect4ha(k)+dect4hb(k))/105.d0
     .  -dedt2d(k)
     .  +(-dedt3qa(k)+dedt3qb(k))/3.d0
     .  -(dedt4oa(k)+dedt4ob(k))/15.d0
     .  +deqt4q(k)/9.d0
      enddo


c temp patch for testing
c      do k=1,3
c       tot1(k)=dect0c(k)
c       tot1(k)=dect1da(k)-dect1db(k)
c       tot1(k)=(dect2qa(k)+dect2qb(k))/3.d0
c       tot1(k)=-dedt2d(k)
c      enddo

       c1(n1,j)=c1(n1,j)-2.d0*eps
       call tfactorfast(c1(n1,:),c2(n2,:),t0,t1,t2,t3,t4,t5,t6)
      call ct0c_d(ch2(n2),ch1(n1),t1,dect0c)
      call ct1d_d(ch2(n2),d1(n1,:),t2,dect1da)
      call ct1d_d(ch1(n1),d2(n2,:),t2,dect1db)
      call ct2q_d(ch2(n2),qu1(n1,:,:),t3,dect2qa)
      call ct2q_d(ch1(n1),qu2(n2,:,:),t3,dect2qb)
      call ct3o_d(ch2(n2),oct1(n1,:,:,:),t4,dect3oa)
      call ct3o_d(ch1(n1),oct2(n2,:,:,:),t4,dect3ob)
      call ct4h_d(ch2(n2),hex1(n1,:,:,:,:),t5,dect4ha)
      call ct4h_d(ch1(n1),hex2(n2,:,:,:,:),t5,dect4hb)
      call dt2d_d(d2(n2,:),d1(n1,:),t3,dedt2d)
      call dt3q_d(d2(n2,:),qu1(n1,:,:),t4,dedt3qa)
      call dt3q_d(d1(n1,:),qu2(n2,:,:),t4,dedt3qb)
      call dt4o_d(d2(n2,:),oct1(n1,:,:,:),t5,dedt4oa)
      call dt4o_d(d1(n1,:),oct2(n2,:,:,:),t5,dedt4ob)
      call qt4q_d(qu2(n2,:,:),qu1(n1,:,:),t5,deqt4q)

       do k=1,3
       tot2(k)=dect0c(k)+dect1da(k)-dect1db(k)
     .  +(dect2qa(k)+dect2qb(k))/3.d0
     .  +(dect3oa(k)-dect3ob(k))/15.d0
     .  +(dect4ha(k)+dect4hb(k))/105.d0
     .  -dedt2d(k)
     .  +(-dedt3qa(k)+dedt3qb(k))/3.d0
     .  -(dedt4oa(k)+dedt4ob(k))/15.d0
     .  +deqt4q(k)/9.d0
       enddo

c temp patch for testing
c      do k=1,3
c       tot2(k)=dect0c(k)
c       tot2(k)=dect1da(k)-dect1db(k)
c       tot2(k)=(dect2qa(k)+dect2qb(k))/3.d0
c       tot2(k)=-dedt2d(k)
c      enddo

c reset c1
       c1(n1,j)=c1(n1,j)+eps

       do k=1,3
       dd(j,k)=(tot1(k)-tot2(k))*0.5d0/eps
c      secd1(n1,j,n1,k)=secd1(n1,j,n1,k)+(tot1(k)-tot2(k))*0.5d0/eps
c      secd2(n1,j,n2,k)=-(tot1(k)-tot2(k))*0.5d0/eps
c      secd3(n2,j,n2,k)=secd3(n2,j,n2,k)+(tot1(k)-tot2(k))*0.5d0/eps
       enddo

c end j loop
       enddo

       do j=1,3
       do k=1,3
        secd1(n1,j,n1,k)=0.5d0*(dd(j,k)+dd(k,j))
        secd2(n1,j,n2,k)=-0.5d0*(dd(j,k)+dd(k,j))
        secd3(n2,j,n2,k)=0.5d0*(dd(j,k)+dd(k,j))
       enddo
       enddo

c  end the n2 loop
      enddo
c end the n1 loop
      enddo

1000   continue
c reset the gradients
      do n1=1,natom1
      do k=1,3
       del1(n1,k)=gstore1(n1,k)
      enddo
      enddo
      do n2=1,natom2
      do k=1,3
       del2(n2,k)=gstore2(n2,k)
      enddo
      enddo

c      write(36,*)natom1
c      do n1=1,natom1
c      do k1=1,3
c      do k2=1,k1
c       write(36,9)secd1(n1,k1,n1,k2),secd1(n1,k1,n1,k2)-secd1(n1,k2,n1,k1),k1,k2
c      enddo
c      enddo
c      enddo
c      write(36,*)natom2
c      do n1=1,natom2
c      do k1=1,3
c      do k2=1,k1
c       write(36,9)secd3(n1,k1,n1,k2),secd3(n1,k1,n1,k2)-secd3(n1,k2,n1,k1),k1,k2
c      enddo
c      enddo
c      enddo
c9     format(2e15.6,2I10)

c     do k=1,3
c      tot(k)=0.d0
c      do n1=1,natom1
c      tot(k)=tot(k)+del1(n1,k)
c      enddo
c      do n2=1,natom2
c      tot(k)=tot(k)+del2(n2,k)
c      enddo
c     enddo

c     write(6,*)'elect ',(tot(k),k=1,3)

c     do n1=1,natom1
c     do n2=1,natom2
c     do j=1,3
c     do k=1,3
c      if(abs(secd1(n1,j,n1,k)).gt.1.d-1)write(6,*)'elect ',n1,j,n1,k,1
c      if(abs(secd2(n1,j,n2,k)).gt.1.d-1)write(6,*)'elect ',n1,j,n2,k,2
c      if(abs(secd3(n2,j,n2,k)).gt.1.d-1)write(6,*)'elect ',n2,j,n2,k,3
c     enddo
c     enddo
c     enddo
c     enddo

c check seconds
c      do n1=1,natom1
c      do k1=1,3
c      do k2=1,3
c       sumk=0.d0
cc      do n2=1,natom1
c       sumk=sumk+secd1(n1,k1,n1,k2)
cc      enddo
c       do n2=1,natom2
c       sumk=sumk+secd2(n1,k1,n2,k2)
c       enddo
c       write(35,*)sumk
c      enddo
c      enddo
c      enddo
c      write(35,*)
c      do n2=1,natom2
c      do k1=1,3
c      do k2=1,3
c       sumk=0.d0
c       do n1=1,natom1
c       sumk=sumk+secd2(n1,k1,n2,k2)
c       enddo
cc      do n1=1,natom2
c       sumk=sumk+secd3(n2,k1,n2,k2)
cc      enddo
c       write(35,*)sumk
c      enddo
c      enddo
c      enddo

      return
      end


