module dalton

use molecule_adt
contains

subroutine print_dalton(file_name,mol,nchg)

implicit none
character*80      ::file_name
type(molecule)    ::mol
real*8            ::anum_real
real*8            ::coords_bohr(3)
integer           ::bin_array(2,mol%natom)
integer           ::acount,anum,freq,current,nchg
integer           ::atype,atom1,m,n1 !counter
character*25      ::fmth,myFmt,fmt1,fmt2,fmt3
character*20      ::ca

call histogram(mol,bin_array,acount)

open(unit=26,file=file_name,status='unknown')
! fmth=get_format_head(acount)
      call filelabel(acount,ca)
      n1=index(ca,' ') - 1
      if(nchg.eq.0.and.n1.eq.1)then
       fmtH='(a10,i1,a11)'
       write(26,fmtH)'Atomtypes=',acount,' Nosymmetry'
      endif
      if(nchg.eq.0.and.n1.eq.2)then
       fmtH='(a10,i2,a11)'
       write(26,fmtH)'Atomtypes=',acount,' Nosymmetry'
      endif
      if(nchg.lt.0.and.n1.eq.1)then
       fmtH='(a10,i1,a8,i2,a11)'
       write(26,fmtH)'Atomtypes=',acount,' Charge=',nchg,' Nosymmetry'
      endif
      if(nchg.lt.0.and.n1.eq.2)then
       fmtH='(a10,i2,a8,i2,a11)'
       write(26,fmtH)'Atomtypes=',acount,' Charge=',nchg,' Nosymmetry'
      endif
      if(nchg.gt.0.and.n1.eq.1)then
       fmtH='(a10,i1,a8,i1,a11)'
       write(26,fmtH)'Atomtypes=',acount,' Charge=',nchg,' Nosymmetry'
      endif
      if(nchg.gt.0.and.n1.eq.2)then
       fmtH='(a10,i2,a8,i1,a11)'
       write(26,fmtH)'Atomtypes=',acount,' Charge=',nchg,' Nosymmetry'
      endif
! write(26,fmtH)'Atomtypes=',acount,' Nosymmetry'
atype_loop: do atype=1,acount
   myFmt=get_format2(bin_array(1,atype),bin_array(2,atype))
   !write(26,*),myFmt
   anum=bin_array(1,atype)
   anum_real=real(anum)
   freq=bin_array(2,atype)
   fmt1=get_format(anum,1)
   fmt2=get_format(anum,10)  !two digit case
   fmt3=get_format(anum,100) ! three digit case
   write(26,myFmt)'Charge=',anum_real,' Atoms=',freq,'Basis=aug-cc-pVDZ'

   current=1 
   write_loop: do m=1,mol%natom
      coords_bohr=mol%atoms(m)%Coords*ang2bohr
      if(alab2anum(mol%atoms(m)%label)==anum)then
         if(current<10)then
            write(26,fmt1)mol%atoms(m)%label,current,coords_bohr
            !write(26,fmt1)mol%atoms(m)%label,current,mol%atoms(m)%Coords(1),&
            !mol%atoms(m)%Coords(2),mol%atoms(m)%Coords(3)
         elseif(current>9.and.current<100)then
            write(26,fmt2)mol%atoms(m)%label,current,mol%atoms(m)%Coords(1),&
            mol%atoms(m)%Coords(2),mol%atoms(m)%Coords(3)
         elseif(current>99.and.current<1000)then
            write(26,fmt3)mol%atoms(m)%label,current,mol%atoms(m)%Coords(1),&
            mol%atoms(m)%Coords(2),mol%atoms(m)%Coords(3)
         endif
         if((current)==freq)cycle atype_loop   
         current=current+1
      endif
   enddo write_loop

enddo atype_loop
close(unit=26)
return
end subroutine print_dalton

function get_format_head(acount) result(fmtH)

character*25      ::fmtH
integer           ::acount,acount_len

acount_len=get_intlength(acount)
write(fmtH,10)acount_len
10 format('(A10,I',I1',A11)')

return
end function get_format_head

function get_format(anum,freq) result(myFmt)

character*25      ::myFmt
integer           ::anum,freq,freq_len

!if(freq<10)then
!   freq_len=1
!elseif((freq>=10).and.(freq<100))then
!   freq_len=2
!elseif((freq>=100).and.(freq<1000))then
!   freq_len=3
!endif
freq_len=get_intlength(freq)
write(myFmt,10)len_Trim(anum2alab(anum)),(len_Trim(anum2alab(anum))+1),freq_len
10 format('(A',I1,',T',I1,',I'I1',3f13.6)')

return
end function get_format


function get_format2(anum,freq) result(myFmt)

character*25      ::myFmt
integer           ::anum,freq,freq_len,anum_len

anum_len=get_intlength(anum)
freq_len=get_intlength(freq)

write(myFmt,10)anum_len+2,freq_len
10 format('(A7,F',I1,'.1,A7,I',I1,',1X,A17)')

return
end function get_format2

function get_intlength(myInteger) result(myLen)

integer        ::myInteger,myLen

myLen=1
if(myInteger<10)then
   myLen=1
elseif((myInteger>=10).and.(myInteger<100))then
   myLen=2
elseif((myInteger>=100).and.(myInteger<1000))then
   myLen=3
endif
return
end function get_intlength

subroutine histogram(mol,bin_array,acount)

implicit none

integer        ::acount ! number of distinct atoms
integer        ::bin_array(:,:)
integer        ::anum,aindex
logical        ::seen(110)
integer        ::n !counter
type(molecule) ::mol


seen=.false.
acount=1
bin_array=0

do n=1,mol%natom 
   anum=alab2anum(mol%atoms(n)%label)
   if(seen(anum) .eqv. .false.)then
      bin_array(1,acount)=anum
      bin_array(2,acount)=1
      seen(anum)=.true.
      acount=acount+1
   else
      aindex=get_index(anum,bin_array)
      write(*,*)aindex, anum
      bin_array(2,aindex)=bin_array(2,aindex)+1
   endif
enddo


acount=acount-1
!write(*,*)'Anum     Freq'
!do n=1,mol%natom
!   write(*,*)bin_array(:,n)
!enddo
!write(*,*)'leaving histo'

return
end subroutine histogram



function get_index(anum,bin_array) result(aindex)

implicit none

integer        ::anum,aindex
integer        ::bin_array(:,:)
integer        ::n,natom

aindex=0  
natom=size(bin_array,2)

do n=1,natom
   if(bin_array(1,n)==anum)then
      aindex=n
      return
   endif
enddo



return
end function get_index

end module dalton
