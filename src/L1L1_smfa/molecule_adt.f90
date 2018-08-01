module molecule_adt

type molecule
   integer :: natom
   type(atom),pointer :: atoms(:)
end type molecule

type atom
   character*2  :: label
   logical     :: real_or_not
   integer     :: orig_atom_no1
   integer     :: orig_atom_no2
   real*8 :: Coords(3)
   real*8 :: Charge
   real*8  :: Dipole(3)
   real*8  :: Quadrupole(6)
   real*8  :: Octapole(10)
   real*8  :: Hexadecapole(15)
end type atom

type two_c_bond
   integer  :: a1
   integer  :: a2
end type two_c_bond

type interaction
   integer        :: frag_no1
   integer        :: frag_no2
   integer        :: groups1(3)
   integer        :: groups2(3)
   type(molecule) :: frag1
   type(molecule) :: frag2
   integer        :: isg_coeff
   real*8         :: least_r
   logical        :: ab_initio
end type interaction

type level 
   integer        ::groups(2)
   integer        ::coeff
end type

real*8,parameter       :: bohr2ang=0.529177249d0
real*8,parameter       :: ang2bohr=1.8897259886d0

interface init
   module procedure init_atom
   module procedure init_mol
end interface

interface operator (-)
   module procedure subtract_atoms
   module procedure subtract_mol
end interface

interface assignment (=)
   module procedure assign_atom
   module procedure assign_mol
  ! module procedure copy_atom
  ! module procedure copy_mol
end interface 

interface copy
   module procedure copy_atom
   module procedure copy_mol
end interface

interface operator (+)
   module procedure add_atoms
   module procedure add_mol 
end interface

interface print
   module procedure print_atom
   module procedure print_mol
end interface

interface operator (*)
   module procedure multiply_atom_real
   module procedure multiply_atom_int
   module procedure multiply_mol
end interface
contains


subroutine init_atom(this)
   type(atom) :: this

   return
end subroutine init_atom

subroutine init_mol(mol)
   type(molecule) :: mol

   mol%natom=0
   nullify(mol%atoms)
   return
end subroutine init_mol

subroutine init_interaction(inter)
   type(interaction) :: inter

   inter%frag_no1=0
   inter%frag_no2=0
   inter%groups1=0
   inter%groups2=0
   call init_mol(inter%frag1)
   call init_mol(inter%frag2)
   inter%isg_coeff=0
   inter%least_r=0.0
   inter%ab_initio=.false.
   return
end subroutine init_interaction

subroutine new_mol(mol,natom)
   type(molecule)    :: mol
   integer           :: natom

   call destroy_mol(mol)
   allocate(mol%atoms(natom))
   mol%natom=natom
   
   return
end subroutine new_mol
subroutine new_mol_inarray(mol,natom)
   type(molecule)    :: mol
   integer           :: natom

   allocate(mol%atoms(natom))
   mol%natom=natom
   
   return
end subroutine new_mol_inarray

subroutine destroy_mol(mol)
   type(molecule)    :: mol

   if(associated(mol%atoms)) deallocate(mol%atoms)
   mol%natom=0

   return
end subroutine destroy_mol

logical function same_array(array1,array2)

implicit none
integer              ::array1(:)
integer              ::array2(:)
integer              ::len1,len2,n

same_array=.false.
len1=size(array1)
len2=size(array2)
if(len1/=len2)return
do n=1,len1
   if(array1(n)/=array2(n))return
enddo
same_array=.true.
return
end function same_array   

logical function sub_array(array1,array2)

implicit none
integer              ::array1(:)
integer              ::array2(:)
integer              ::len1,len2,n
integer              ::el1,el2

write(409,*)'-------------------------'
write(409,'(A10,10I3)')'array1= ',array1
write(409,'(A10,10I3)')'array2= ',array2
sub_array=.false.
len1=size(array1)
len2=size(array2)
if(len2==len1)return
if(len1<len2)then
   !array1 is smaller
   array1_loop: do el1=1,len1
      array2_loop: do el2=1,len2
         if(array1(el1)==array2(el2))then
            cycle array1_loop
         endif
         enddo array2_loop
         return !if we make it here some element in array1 wasn't found in array2
      enddo array1_loop
      sub_array=.true.
      write(409,*)'A1 is contained in A2'
endif

return
end function sub_array

integer function first_zero(array)

integer,intent(in)        ::array(:)
integer                   ::n,length

length=size(array)

do n=1,length
   if(array(n)==0)then
      first_zero=n
      return
   endif
enddo
first_zero=0

return
end function first_zero

logical function same_atom(atom1,atom2)
   type(atom),intent(in) :: atom1,atom2
   real*8,parameter :: tol=1.d-5
   same_atom=.false.

   !write(101,*)
   !write(101,*)atom1%label,atom1%Coords(1),atom1%Coords(2),atom1%Coords(3)
   !write(101,*)atom2%label,atom2%Coords(1),atom2%Coords(2),atom2%Coords(3)

      if(atom1%label.ne.atom2%label) return
         !if(sqrt((atom1%Coords(1)-atom2%Coords(1))**2).gt.tol)return
         !   if(sqrt((atom1%Coords(2)-atom2%Coords(2))**2).gt.tol)return
         !      if(sqrt((atom1%Coords(3)-atom2%Coords(3))**2).gt.tol)return
      if(sqrt(sum((atom1%Coords-atom2%Coords)**2)).gt.tol) return
   same_atom=.true.
   !write(101,*)same_atom
   return
end function same_atom

logical function contains_atom(mol_in, atom1)
   type(molecule),intent(in)  :: mol_in
   type(atom),intent(in) :: atom1
   integer :: n
   real*8,parameter :: tol=1.d-7
   contains_atom=.true.

   do n=1,mol_in%natom
      if(same_atom(mol_in%atoms(n),atom1).eqv..true.) return 
   enddo
   contains_atom=.false.
   
end function contains_atom

logical function same_inter(inter1,inter2)
implicit none
type(interaction),intent(in)     ::inter1,inter2
   
same_inter=.false.
if(inter1%frag_no1/=inter2%frag_no1)return
if(inter1%frag_no2/=inter2%frag_no2)return
same_inter=.true.
return
end function same_inter

logical function exact_same_inter(inter1,inter2)
implicit none
type(interaction),intent(in)     ::inter1,inter2
   
exact_same_inter=.false.
if(inter1%frag_no1/=inter2%frag_no1)return
if(inter1%frag_no2/=inter2%frag_no2)return
if(inter1%ab_initio.neqv.inter2%ab_initio)return
exact_same_inter=.true.
return
end function exact_same_inter

logical function same_interaction(inter1,inter2)
implicit none
type(interaction),intent(in)     ::inter1,inter2

same_interaction=.false.
if(same_array(inter1%groups1,inter2%groups1) .eqv. .false.)return
if(same_array(inter1%groups2,inter2%groups2) .eqv. .false.)return
same_interaction=.true.
return
end function same_interaction

logical function exact_same_interaction(inter1,inter2)
implicit none
type(interaction),intent(in)     ::inter1,inter2

exact_same_interaction=.false.
if(same_array(inter1%groups1,inter2%groups1) .eqv. .false.)return
if(same_array(inter1%groups2,inter2%groups2) .eqv. .false.)return
if(inter1%ab_initio.neqv.inter2%ab_initio)return
exact_same_interaction=.true.
return
end function exact_same_interaction

!subroutine multiply_mol(c,mol)
!   type(molecule),intent(out)  :: mol
!   real,          intent(in)  :: c
!   integer           :: natom
!
!   natom=mol%natom
!   do n=1,natom
!      mol%atoms(n)=c*mol%atoms(n)
!      !call print(mol%atoms(n))
!   !write(*,*) 'Still haven''t died!'
!   enddo
!  return 
!end subroutine multiply_mol

type(molecule) function multiply_mol(c,mol)

   type(molecule),intent(in)  :: mol
   real,          intent(in)  :: c
   integer           :: natom
   
   natom=mol%natom
   call new_mol(multiply_mol,mol%natom)
   !write(*,*) 'Up to here!'
   do n=1,natom
      !multiply_mol%atoms(n)=multiply_atom_real(c,mol%atoms(n))
      multiply_mol%atoms(n)=c*mol%atoms(n)
      !call print(multiply_mol%atoms(n))
   !write(*,*) 'Still haven''t died!'
   enddo
   !write(*,*) 'Aint ded yet'
   return
end function multiply_mol

!subroutine multiply_mol(mol_in,mol_out,coeff)
!   type(molecule),intent(in)  :: mol_in
!   type(molecule),intent(out) :: mol_out
!   real,intent(in)         :: coeff
!   integer           :: natom
!   
!   natom=mol_in%natom
!   call new_mol(mol_out,natom)
!   !mol_out%natom=natom !added in
!   do n=1,natom
!      mol_out%atoms(n)=multiply_atom_real(coeff,mol_in%atoms(n))
!   enddo
!      
!   return
!end subroutine multiply_mol

subroutine copy_mol(mol_in,mol_out)
   type(molecule),intent(in)  :: mol_in
   type(molecule),intent(out) :: mol_out

   integer :: natom,n

   natom=mol_in%natom
   call new_mol(mol_out,natom)
   mol_out%natom=natom
   do n=1,natom
      call copy_atom(mol_in%atoms(n),mol_out%atoms(n))
   enddo

return
end subroutine copy_mol

subroutine copy_atom(atom_in,atom_out)
   type(atom),intent(in)  :: atom_in
   type(atom),intent(out) :: atom_out

   atom_out%label        =atom_in%label
   atom_out%real_or_not  =atom_in%real_or_not
   atom_out%orig_atom_no1=atom_in%orig_atom_no1
   atom_out%orig_atom_no2=atom_in%orig_atom_no2
   atom_out%Coords       =atom_in%Coords
   atom_out%Charge       =atom_in%Charge
   atom_out%Dipole       =atom_in%Dipole
   atom_out%Quadrupole   =atom_in%Quadrupole
   atom_out%Octapole     =atom_in%Octapole
   atom_out%Hexadecapole =atom_in%Hexadecapole
   return
end subroutine copy_atom

type(atom) function multiply_atom_real(c,atom_in)
   type(atom),intent(in):: atom_in
   real,      intent(in):: c

   multiply_atom_real%label        =atom_in%label
   multiply_atom_real%real_or_not  =atom_in%real_or_not
   multiply_atom_real%orig_atom_no1=atom_in%orig_atom_no1
   multiply_atom_real%orig_atom_no2=atom_in%orig_atom_no2
   multiply_atom_real%Coords       =atom_in%Coords
   multiply_atom_real%Charge       =c*atom_in%Charge
   multiply_atom_real%Dipole       =c*atom_in%Dipole
   multiply_atom_real%Quadrupole   =c*atom_in%Quadrupole
   multiply_atom_real%Octapole     =c*atom_in%Octapole
   multiply_atom_real%Hexadecapole =c*atom_in%Hexadecapole   
   return
end function multiply_atom_real

type(atom) function multiply_atom_int(i,atom_in)
   type(atom),intent(in) :: atom_in
   integer,   intent(in) :: i
   multiply_atom_int=multiply_atom_real(real(i),atom_in)
   return
end function multiply_atom_int

subroutine assign_mol(mol_out,mol_in)
   type(molecule),intent(out):: mol_out
   type(molecule),intent(in) :: mol_in

   call copy_mol(mol_in,mol_out)

   return
end subroutine assign_mol

subroutine assign_atom(atom_out,atom_in)
   type(atom),intent(out):: atom_out
   type(atom),intent(in) :: atom_in

   call copy_atom(atom_in,atom_out)
   return
end subroutine assign_atom

type(molecule) function add_mol(mol1,mol2)
!subroutine add_mol(mol1,mol2,mol_out)
   type(molecule),intent(in)  :: mol1,mol2
   type(molecule)             :: mol_temp
   !type(atom)     ::summed_atom

   integer  :: n,m,nu,num_total
   integer  :: num_current

   integer  :: atom_map(2,mol1%natom+mol2%natom)
  
   atom_map=0 
   do n=1,mol1%natom
      atom_map(1,n)=n   
   enddo
   nu=mol1%natom
   
   mol2_loop: do m=1,mol2%natom
      ! search through mol1 for a match
      mol1_loop: do n=1,mol1%natom
        if(same_atom(mol1%atoms(n),mol2%atoms(m)))then
            atom_map(2,n)=m
            cycle mol2_loop
         endif   
      enddo mol1_loop

      ! if we didn't find a match
      nu=nu+1
      atom_map(2,nu)=m
   enddo mol2_loop
 
   !      write(*,*) 'Atom map'
   !do n=1,nu
   !      write(*,*) atom_map(1,n),atom_map(2,n)    
   !enddo
   
    allocate(add_mol%atoms(nu))
    add_mol%natom=nu 
   !call new_mol(mol_temp,nu)
   do n=1,nu
      if((atom_map(1,n)*atom_map(2,n)).ne.0)then
         add_mol%atoms(n)=mol1%atoms(atom_map(1,n))+mol2%atoms(atom_map(2,n))
      elseif(atom_map(1,n).eq.0)then
         add_mol%atoms(n)=mol2%atoms(atom_map(2,n))
      elseif(atom_map(2,n).eq.0)then
         add_mol%atoms(n)=mol1%atoms(atom_map(1,n))
      endif
   enddo 
 
   !call print(mol_temp)
   !call copy(mol_temp,add_mol)
   !add_mol=mol_temp 

   !   
!   
return
end function add_mol  

type(molecule) function complement(mol1,mol2)
   type(molecule),intent(in)  :: mol1,mol2

   integer  :: n,m,nc,nu,ni !nc = num_complement, nu = num_union, ni = num_intersection

   integer  :: atom_map(2,mol1%natom+mol2%natom)

   atom_map=0
   do n=1,mol1%natom
      atom_map(1,n)=n
   enddo
   nu=mol1%natom

   mol2_loop: do m=1,mol2%natom
      ! search through mol1 for a match
      mol1_loop: do n=1,mol1%natom
        if(same_atom(mol1%atoms(n),mol2%atoms(m)))then
            atom_map(2,n)=m
            ni=ni+1
            cycle mol2_loop
         endif
      enddo mol1_loop

      ! if we didn't find a match
      nu=nu+1
      atom_map(2,nu)=m
   enddo mol2_loop   

   nc=nu-ni
  !       write(*,*) 'Atom map'
  !       write(*,*)'nc=', nc,' nu=', nu, ' ni=',ni 
  ! do n=1,nu
  !       write(*,*) atom_map(1,n),atom_map(2,n)    
  ! enddo

   allocate(complement%atoms(nc))
   complement%natom=nc
   !call new_mol(mol_temp,nc)
   m=1
   do n=1,nu
      if((atom_map(1,n)*atom_map(2,n)).eq.0)then
         if(atom_map(1,n).eq.0)then
            complement%atoms(m)=mol2%atoms(atom_map(2,n))
            write(*,*) atom_map(1,n),atom_map(2,n)
            m=m+1    
         elseif(atom_map(2,n).eq.0)then
            complement%atoms(m)=mol1%atoms(atom_map(1,n))
            write(*,*) atom_map(1,n),atom_map(2,n)    
            m=m+1    
         endif
      endif
   enddo

return
end function complement

type(molecule) function remove_h(mol)
   type(molecule),intent(in)     ::mol
   integer                       ::n,m,ns !ns = num_skeleton (ie. non H)

   do m=1,mol%natom
      if(mol%atoms(m)%label/='H')then
         ns=ns+1
      endif
   enddo
   allocate(remove_h%atoms(ns))
   remove_h%natom=ns
   write(*,*)'ns= ', ns

!            write(34,*)'removing h'
   n=1
   do m=1,mol%natom
      if(mol%atoms(m)%label/='H')then
         remove_h%atoms(n)=mol%atoms(m)
!            write(34,*)n,m,remove_h%atoms(n)%orig_atom_no1
         if(remove_h%atoms(n)%orig_atom_no1==0)then
!            write(34,*)'here'
            remove_h%atoms(n)%orig_atom_no1=m
!            write(34,*)n,m
         endif
         n=n+1
      endif
   enddo 

return   
end function remove_h

function orig_atom(heavy_mol,atom1) result(heavy_no)
!takes a heavy_molecule and an atom number and returns the atom number in the heavy molecule
   type(molecule),intent(in)     ::heavy_mol
   integer                       ::atom1
   integer                       ::n
   integer                       ::heavy_no


do n=1,heavy_mol%natom
   if(heavy_mol%atoms(n)%orig_atom_no1==atom1)then
      heavy_no=n
      write(*,*)atom1,n
      return
   endif
enddo
end function orig_atom   

type(molecule) function all_caps(heavy_mol) !TEST
implicit none
type(molecule),intent(in)        ::heavy_mol
integer                          ::atom1,atom2
integer                          ::ncaps,na

ncaps=heavy_mol%natom*(heavy_mol%natom-1)
allocate(all_caps%atoms(ncaps))
all_caps%natom=ncaps

na=1

do atom1=1,heavy_mol%natom-1
 do atom2=atom1+1,heavy_mol%natom
      all_caps%atoms(na)=make_capH(heavy_mol%atoms(atom1),heavy_mol%atoms(atom2)); na=na+1
      all_caps%atoms(na)=make_capH(heavy_mol%atoms(atom2),heavy_mol%atoms(atom1)); na=na+1
   enddo
enddo

return      
end function all_caps

subroutine fix_capstatus(real_mol,mol)
implicit none
type(molecule),intent(in)        ::real_mol
type(molecule)                   ::mol
!type(molecule),intent(inout)     ::mol
integer                          ::atom1,atom2

!write(102,*)mol%natom

mol_loop: do atom1=1,mol%natom
   if(mol%atoms(atom1)%label/='H')then
      mol%atoms(atom1)%real_or_not=.true.
      !loop anyway to get orig_atom_no
      oa_loop: do atom2=1,real_mol%natom
         if(same_atom(mol%atoms(atom1),real_mol%atoms(atom2)).eqv..true.)then
            mol%atoms(atom1)%orig_atom_no1=real_mol%atoms(atom2)%orig_atom_no1
            cycle mol_loop
         endif
      enddo oa_loop
      cycle mol_loop
   else !test if it's a cap or not
      real_mol_loop: do atom2=1,real_mol%natom
         if(same_atom(mol%atoms(atom1),real_mol%atoms(atom2)).eqv..true.)then
            !write(101,*)'Its true'
            mol%atoms(atom1)%real_or_not=.true.
            mol%atoms(atom1)%orig_atom_no1=real_mol%atoms(atom2)%orig_atom_no1
            !write(101,*)mol%atoms(atom1)%real_or_not
            cycle mol_loop
         endif
      enddo real_mol_loop
      mol%atoms(atom1)%real_or_not=.false.
   endif
enddo mol_loop

!   do atom1=1,mol%natom
!      write(102,*)mol%atoms(atom1)%label,mol%atoms(atom1)%Coords,mol%atoms(atom1)%real_or_not
!   enddo

return
end subroutine fix_capstatus


type(molecule) function delete_atom(mol,atom1)

implicit none
type(molecule),intent(in)        ::mol
type(atom),intent(in)            ::atom1
integer                          ::num_atoms
integer                          ::m,alloc

num_atoms=mol%natom-1

allocate(delete_atom%atoms(num_atoms))
delete_atom%natom=num_atoms

alloc=1
do m=1,mol%natom
   if(same_atom(mol%atoms(m),atom1) .eqv. .false.)then
      delete_atom%atoms(alloc)=mol%atoms(m)
      alloc=alloc+1
   endif
enddo

return
end function delete_atom

type(atom) function add_atoms(atom1,atom2)
   type(atom),intent(in)  :: atom1,atom2

   if(same_atom(atom1,atom2).eqv..true.)then !else I want to die with an error
      add_atoms%label         =atom1%label
      add_atoms%real_or_not   =atom1%real_or_not
      add_atoms%orig_atom_no1 =atom1%orig_atom_no1
      add_atoms%orig_atom_no2 =atom1%orig_atom_no2
      add_atoms%Coords        =atom1%Coords
      add_atoms%Charge        =atom1%Charge+atom2%Charge
      add_atoms%Dipole        =atom1%Dipole+atom2%Dipole
      add_atoms%Quadrupole    =atom1%Quadrupole+atom2%Quadrupole
      add_atoms%Octapole      =atom1%Octapole+atom2%Octapole
      add_atoms%Hexadecapole  =atom1%Hexadecapole+atom2%Hexadecapole
   endif
   return
end function add_atoms
      
type(atom) function subtract_atoms(atom1,atom2)!This needs testing
   type(atom),intent(in)  :: atom1,atom2

   if(same_atom(atom1,atom2).eqv..true.)then !else I want to die with an error
      subtract_atoms%label         =atom1%label
      subtract_atoms%real_or_not   =atom1%real_or_not
      subtract_atoms%orig_atom_no1 =atom1%orig_atom_no1
      subtract_atoms%orig_atom_no2 =atom1%orig_atom_no2
      subtract_atoms%Coords        =atom1%Coords
      subtract_atoms%Charge        =atom1%Charge-atom2%Charge
      subtract_atoms%Dipole        =atom1%Dipole-atom2%Dipole
      subtract_atoms%Quadrupole    =atom1%Quadrupole-atom2%Quadrupole
      subtract_atoms%Octapole      =atom1%Octapole-atom2%Octapole
      subtract_atoms%Hexadecapole  =atom1%Hexadecapole-atom2%Hexadecapole
   endif
   return
end function subtract_atoms
      
!subroutine add_atoms(atom1,atom2,atom_out)
!type(atom),intent(in)  :: atom1,atom2
!   type(atom),intent(out) :: atom_out
!
!   if(same_atom(atom1,atom2).eqv..true.)then !else I want to die with an error
!      atom_out%label=atom1%label
!      atom_out%real_or_not=atom1%real_or_not
!      atom_out%orig_atom_no1=atom1%orig_atom_no1
!      atom_out%orig_atom_no2=atom1%orig_atom_no2
!      atom_out%Coords =atom1%Coords
!      atom_out%Charge =atom1%Charge+atom2%Charge
!      atom_out%Dipole =atom1%Dipole+atom2%Dipole
!      atom_out%Quadrupole =atom1%Quadrupole+atom2%Quadrupole
!      atom_out%Octapole =atom1%Octapole+atom2%Octapole
!      atom_out%Hexadecapole =atom1%Hexadecapole+atom2%Hexadecapole
!   endif
!return
!end subroutine add_atoms

type(molecule) function subtract_mol(mol1,mol2)
   type(molecule),intent(in)  :: mol1,mol2
   
   subtract_mol=(-1.0)*mol2+mol1
   return
end function subtract_mol

function centroidm(mol) result(cent)
type(molecule),intent(in)   ::mol
real*8                      ::mass(mol%natom)
real*8                      ::cent(3)
integer                     ::atom1

!Temp vector of masses
mass=0
do atom1=1,mol%natom
   mass(atom1)=anum2amass(alab2anum(mol%atoms(atom1)%label))
enddo

cent(1)=(sum(mass*(mol%atoms%Coords(1))))/sum(mass)
cent(2)=(sum(mass*(mol%atoms%Coords(2))))/sum(mass)
cent(3)=(sum(mass*(mol%atoms%Coords(3))))/sum(mass)

end function centroidm


function centroida(mol) result(cent)
type(molecule),intent(in)   ::mol
integer                     ::anum(mol%natom)
real*8                      ::cent(3)
integer                     ::atom1

!Temp vector of atomic numbers 
anum=0
do atom1=1,mol%natom
anum(atom1)=alab2anum(mol%atoms(atom1)%label)
enddo

cent(1)=(sum(anum*(mol%atoms%Coords(1))))/sum(anum)
cent(2)=(sum(anum*(mol%atoms%Coords(2))))/sum(anum)
cent(3)=(sum(anum*(mol%atoms%Coords(3))))/sum(anum)

end function centroida

function get_mult(mol) result(mult) 
type(molecule),intent(in)   ::mol
integer                     ::anum(mol%natom)
real*8                      ::acharge(mol%natom)
integer                     ::atom1,mult,nelec


!Temp vector of atomic numbers and atomic charges 
anum=0
do atom1=1,mol%natom
anum(atom1)=alab2anum(mol%atoms(atom1)%label)
acharge(atom1)=mol%atoms(atom1)%Charge
enddo
nelec=nint(sum(acharge))-sum(anum)

if(mod(nelec,2)==0)then
   mult=1
else
   mult=2
endif

end function get_mult

function get_charge(mol) result(charge) 
type(molecule),intent(in)   ::mol
real*8                      ::acharge(mol%natom)
integer                     ::atom1
real*8                      ::charge


!Temp vector of atomic charges 
anum=0
do atom1=1,mol%natom
acharge(atom1)=mol%atoms(atom1)%Charge
!write(222,*)mol%atoms(atom1)%label,mol%atoms(atom1)%Charge,mol%atoms(atom1)%orig_atom_no1
write(222,*)mol%atoms(atom1)%Coords(1),mol%atoms(atom1)%Coords(2),mol%atoms(atom1)%Coords(3)
enddo
write(222,*)
charge=sum(acharge)
write(*,*)'Charge is ',nint(charge)
end function get_charge

real*8 function bond_dist(atom1,atom2)
   type(atom),intent(in)  :: atom1,atom2

   if(same_atom(atom1,atom2).eqv..true.)then
      bond_dist = 0.0
   else
      bond_dist = dsqrt(sum((atom1%Coords-atom2%Coords)**2))
   endif
 
return   
end function bond_dist

real*8 function least_dist(mol1,mol2)
   type(molecule),intent(in)     ::mol1,mol2
   integer                       ::m,n
   real*8                        ::tmp_dist

least_dist=bond_dist(mol1%atoms(1),mol2%atoms(1))
   
   do n=1,mol1%natom
      do m=1,mol2%natom
         tmp_dist=bond_dist(mol1%atoms(n),mol2%atoms(m))
         if(tmp_dist.lt.least_dist)least_dist=tmp_dist
      enddo
   enddo
   
return
end function least_dist

real*8 function least_dist_real(mol1,mol2)
   type(molecule),intent(in)     ::mol1,mol2
   integer                       ::m,n
   real*8                        ::tmp_dist

   !Least distance between real atoms (i.e. exlcuding capH)

least_dist_real=bond_dist(mol1%atoms(1),mol2%atoms(1))
   
   do n=1,mol1%natom
      do m=1,mol2%natom
         write(99,*)n,mol1%atoms(n)%label,mol1%atoms(n)%real_or_not,m,mol2%atoms(m)%label,mol2%atoms(m)%real_or_not
         if(mol1%atoms(n)%real_or_not .eqv. .true. .and. mol2%atoms(m)%real_or_not .eqv. .true.)then
            tmp_dist=bond_dist(mol1%atoms(n),mol2%atoms(m))
            if(tmp_dist.lt.least_dist_real)least_dist_real=tmp_dist
         endif
      enddo
   enddo
   
return
end function least_dist_real

   
subroutine read_mol(mol_out,file_no_in)

   integer,intent(in),optional   ::file_no_in
   type(molecule)     :: mol_out

   integer     :: m,natom,file_no

   if(present(file_no_in))then
      file_no=file_no_in
   else
      file_no=5
   endif
   
   read(file_no,*)  ! read the header line
   read(file_no,*)  natom ! maybe read this number and check against num_atom_frag(n)
   call new_mol(mol_out,natom)
   do m=1,natom
      read(file_no,*) mol_out%atoms(m)%label
      read(file_no,*) mol_out%atoms(m)%Coords
      read(file_no,*)  ! read the "Charge" text line 
      read(file_no,*) mol_out%atoms(m)%Charge
      read(file_no,*)  ! read the "Dipole" text line 
      read(file_no,*) mol_out%atoms(m)%Dipole
      read(file_no,*)  ! read the "Quadrupole" text line 
      read(file_no,*) mol_out%atoms(m)%Quadrupole
      read(file_no,*)  ! read the "Octapole" text line 
      read(file_no,*) mol_out%atoms(m)%Octapole
      read(file_no,*)  ! read the "Hexadecapole" text line 
      read(file_no,*) mol_out%atoms(m)%Hexadecapole
   enddo
return
end subroutine read_mol


subroutine read_centrals(mol_out,npoint,file_no_in)

   integer,intent(in),optional   ::file_no_in
   type(molecule)     :: mol_out

   integer     :: m,npoint,file_no
   character*80  :: tmp
   
   if(present(file_no_in))then
      file_no=file_no_in
   else
      file_no=5
   endif

   call new_mol(mol_out,npoint)
   do m=1,npoint
      read(file_no,*) tmp !header ="Fragment n" 
      read(file_no,*) mol_out%atoms(m)%Coords
      read(file_no,*) mol_out%atoms(m)%Charge
      read(file_no,*) mol_out%atoms(m)%Dipole
      read(file_no,*) mol_out%atoms(m)%Quadrupole
      read(file_no,*) mol_out%atoms(m)%Octapole
      read(file_no,*) mol_out%atoms(m)%Hexadecapole
   enddo
return

end subroutine read_centrals

subroutine print_mol(mol,file_no_out)

   integer,intent(in),optional   ::file_no_out
   type(molecule)     :: mol

   integer     :: m,natom,file_no

   if(present(file_no_out))then
      file_no=file_no_out
   else
      file_no=6
   endif
   
   natom=mol%natom
   
   do m=1,natom
      call print_atom(mol%atoms(m),file_no)
   enddo
end subroutine print_mol

subroutine print_atom(a1,file_no_out)

   integer,intent(in),optional   ::file_no_out
   type(atom)     :: a1

   integer     :: file_no

   if(present(file_no_out))then
      file_no=file_no_out
   else
      file_no=6
   endif

   write(file_no,*) a1%label
   !write(file_no,*) a1%real_or_not
   write(file_no,*) a1%Coords
   write(file_no,*)  'Charge' ! read the "Charge" text line 
   write(file_no,*) a1%Charge
   write(file_no,*)  'Dipole'! read the "Dipole" text line 
   write(file_no,*) a1%Dipole
   write(file_no,*)  'Quadrupole'! read the "Quadrupole" text line 
   write(file_no,*) a1%Quadrupole
   write(file_no,*)  'Octopus'! read the "Octapole" text line 
   write(file_no,*) a1%Octapole
   write(file_no,*)  'Hexadecapole' ! read the "Hexadecapole" text line 
   write(file_no,*) a1%Hexadecapole

end subroutine print_atom

function get_coords(mol)

type(molecule),intent(in) :: mol
real*8 :: get_coords(3,mol%natom)

integer :: n

do n=1,mol%natom
   get_coords(:,n)=mol%atoms(n)%Coords
!   write(*,*)mol%atoms(n)%Coords
enddo

end function get_coords

function get_multipoles(atom_in)

type(atom), intent(in) ::atom_in
real*8 :: get_multipoles(35)

integer  ::n

   get_multipoles(1)=atom_in%Charge
   do n=2,4
      get_multipoles(n)=atom_in%Dipole(n-1)
   enddo
   do n=5,10
      get_multipoles(n)=atom_in%Quadrupole(n-4)
   enddo   
   do n=11,20
      get_multipoles(n)=atom_in%Octapole(n-10)
   enddo   
   do n=21,35
      get_multipoles(n)=atom_in%Hexadecapole(n-20)
   enddo   
   
end function get_multipoles


!subroutine conn_matrix(cmatrix,mol,atomic_data)
function conn_matrix(mol) result(cmat)

implicit none

type(molecule),intent(in) :: mol
integer :: cmat(mol%natom,mol%natom)
integer :: anum(mol%natom)
integer :: nhbonds,mb,nb,mult,m,n,n1,n2 !counters
real*8  ::bond,btol

      write(*,*)'In conn_matrix'
      write(6,*)' number of atoms = ',mol%natom
cmat=0

!come up with a vector of atomic numbers
do m=1,mol%natom
       anum(m)=alab2anum(mol%atoms(m)%label)
enddo
!write(*,*)anum

!now do the conn matrix
do n=1,mol%natom
   do m=1,mol%natom
      bond=bond_dist(mol%atoms(n),mol%atoms(m))
      btol=(anum2cov(anum(n))+anum2cov(anum(m)))+0.40d0
      !write(*,*),n,m,bond,btol
      if((bond.lt.btol).and.(n/=m))cmat(n,m)=1 
   enddo
enddo      

! read in the IN_HBONDS file, so that Matt's L1 fragments will account
! for H-bonds. Try this patch on 200513
!      write(*,*)' calling readhbonds'
!      call readhbonds(mol%natom,cmat)
!      write(*,*)' after calling readhbonds'
! subroutine readhbonds is at the end of this file
! end of patch 200513
! avoid subroutine, do it here 050613

      write(6,*)' opening IN_HBONDS'
      open(unit=128,file='IN_HBONDS',status='unknown')

      read(128,*,end=777)
      read(128,*)nhbonds
       if(nhbonds.eq.0)go to 777
      read(128,*)
      do m=1,nhbonds

      read(128,*)mb,nb,mult
      n1=0
      n2=0
      do n=1,mol%natom
       if(mol%atoms(n)%orig_atom_no1.eq.mb)n1=n
       if(mol%atoms(n)%orig_atom_no1.eq.nb)n2=n
      enddo
       
      
      write(6,*)mb,nb,mult,n1,n2
!      if(cmat(mb,nb).lt.2)then
!      cmat(mb,nb)=mult
!      cmat(nb,mb)=mult
!      endif
       if(cmat(n1,n2).lt.2)then
        cmat(n1,n2)=mult
        cmat(n2,n1)=mult
       endif
      enddo
777   continue
      close(unit=128)
! end patch 050613
      write(*,*)' finished conn_matrix'
return
!end subroutine conn_matrix

end function conn_matrix

function make_capH(atom1,atom2) result(capH)

implicit none

type(atom), intent(in)     ::atom1,atom2
type(atom)                 ::capH
real*8                     ::ratio
integer                    ::anum1,anum2

anum1=alab2anum(atom1%label)
anum2=alab2anum(atom2%label)

ratio =(anum2cov(1)+anum2cov(anum1))/(anum2cov(anum1)+anum2cov(anum2))

!write(*,*)'ratio= ',ratio

capH%label='H'
capH%real_or_not=.false.
capH%orig_atom_no1=atom1%orig_atom_no1
capH%orig_atom_no2=atom2%orig_atom_no1
capH%Coords=(atom2%Coords-atom1%Coords)*ratio+atom1%Coords
capH%Charge=0
capH%Dipole=0
capH%Quadrupole=0
capH%Octapole=0
capH%Hexadecapole=0

return
end function make_capH 

function get_bonded(atom_no,cmat) result(nbonded) 

implicit none

integer,dimension(:,:),intent(in)   :: cmat
integer                 :: n,m !counters
integer                 :: atom_no,natom
integer                 :: nbonds
!integer,allocatable     :: nbonded(:)
integer                 :: nbonded(size(cmat,1))

natom=size(cmat,1)
!write(*,*)natom
!allocate(nbonded(natom))
nbonded=0
nbonds=0

do m=1,natom
!   if(cmat(atom_no,m)>=1)then
! patch on 200513 to allow H-bonds
   if(abs(cmat(atom_no,m))>=1)then
      nbonds=nbonds+1
      nbonded(nbonds)=m
!      write(*,*)atom_no,m
   endif
enddo
   
return
end function get_bonded

function get_bondedH(atom_no,mol) result(nbonded)

!returns a vector with the atom numbers of the H atoms attached to the given atom

implicit none

type(molecule),intent(in)           :: mol
integer :: cmat(mol%natom,mol%natom)
integer                 :: n,m !counters
integer                 :: atom_no,natom
integer                 :: nbonds
integer,allocatable     :: nbonded(:)

cmat=conn_matrix(mol)
!write(*,*)natom
allocate(nbonded(mol%natom))
nbonded=0
nbonds=0

do m=1,mol%natom
   if((cmat(atom_no,m)==1).and.mol%atoms(m)%label=='H')then
      nbonds=nbonds+1
      nbonded(nbonds)=m
!      write(*,*)atom_no,m
   endif
enddo

return
end function get_bondedH

function get_conn_num(atom_no,cmat) result(nbonds)
implicit none

integer,dimension(:,:),intent(in)   :: cmat
integer                 :: m !counters
integer                 :: atom_no,natom
integer                 :: nbonds


natom=size(cmat,1)
!natom=nint(sqrt(real(size(cmat))))
!write(*,*)natom

nbonds=0
do m=1,natom
   if(cmat(atom_no,m)>=1)then
      nbonds=nbonds+1
   endif
enddo
   
return
end function get_conn_num


function detect_rings3(cmat) result(cmatd)

implicit none

integer,dimension(:,:),intent(in)   :: cmat
integer,dimension(size(cmat,1),size(cmat,2))        :: cmatd
integer                 :: n,m,l !counters
integer                 :: natom
integer                 :: nbonds
integer,allocatable     :: nbonded(:)

!natom=nint(sqrt(real(size(cmat)))
natom=size(cmat(:,1))
!write(*,*)natom
allocate(nbonded(natom))
cmatd=cmat
   
!Three member rings
!For triangle 1-2-3: 1 is bonded to 2 and 3, Are 2 and 3  bonded?

do n=1,natom
   nbonds=get_conn_num(n,cmat)
   nbonded=get_bonded(n,cmat)
   !write(*,*)n,'(',nbonds,'):', nbonded
   if(nbonds.gt.1)then
!      write(*,*)'n= ',n, 'in if bit' 
!DIR$ NOVECTOR
      do m=1,nbonds-1
!DIR$ NOVECTOR
         do l=m+1,nbonds
!            write(*,*)'Up to here',nbonded(m),nbonded(l)
            if(cmat(nbonded(m),nbonded(l))>=1)then
               cmatd(n,nbonded(m))=2;cmatd(nbonded(m),n)=2
               cmatd(nbonded(l),nbonded(m))=2;cmatd(nbonded(m),nbonded(l))=2
               cmatd(n,nbonded(l))=2;cmatd(nbonded(l),n)=2
               write(*,*)'Triangle!',n,nbonded(m),nbonded(l)
            endif
         enddo
      enddo
   endif
enddo

return
end function detect_rings3

function detect_rings4(cmat) result(cmatd)

implicit none

integer,dimension(:,:),intent(in)   :: cmat
integer,dimension(size(cmat,1),size(cmat,2))        :: cmatd
integer                 :: atom1,atom2,atom3 !counters
integer                 :: n,m !counters
integer                 :: natom
integer                 :: nbonds_atom1,nbonds_atom2,nbonds_atom3
integer                 :: bonded_atom1(size(cmat,1)),bonded_atom2(size(cmat,1)),bonded_atom3(size(cmat,1))

natom=size(cmat(:,1))
cmatd=cmat

!Four member rings
! For square 1-2-3-4: 1 is bonded to 2 and 3. Are 2 and 3 bonded to something common (i.e. 4)

do atom1=1,natom
   nbonds_atom1=get_conn_num(atom1,cmat)
   bonded_atom1=get_bonded(atom1,cmat)
   if(nbonds_atom1>=2)then
      do atom2=1,nbonds_atom1-1
         do atom3=atom2+1,nbonds_atom1
            nbonds_atom2=get_conn_num(bonded_atom1(atom2),cmat)
            bonded_atom2=get_bonded(bonded_atom1(atom2),cmat)
            nbonds_atom3=get_conn_num(bonded_atom1(atom3),cmat)
            bonded_atom3=get_bonded(bonded_atom1(atom3),cmat)
            do n=1,nbonds_atom2
            do m=1,nbonds_atom3
!           write(*,*)'Atom1 =',atom1
!           write(*,*)'Atom2 =',atom2
!           write(*,*)'Bonded2=',bonded_atom2
!           write(*,*)'Atom3 =',atom3
!           write(*,*)'Bonded3=',bonded_atom3
               if((bonded_atom2(n)==bonded_atom3(m)).and.(bonded_atom2(n)/=atom1))then
                  !We found a square!
                  write(*,*)'Found a square!'
                  write(*,*)bonded_atom2(n),bonded_atom3(m),atom1 
                  cmatd(atom1,bonded_atom1(atom2))=2;cmatd(bonded_atom1(atom2),atom1)=2
                  cmatd(atom1,bonded_atom1(atom3))=2;cmatd(bonded_atom1(atom3),atom1)=2
                  cmatd(bonded_atom1(atom2),bonded_atom2(n))=2;cmatd(bonded_atom2(n),bonded_atom1(atom2))=2
                  cmatd(bonded_atom1(atom3),bonded_atom2(n))=2;cmatd(bonded_atom2(n),bonded_atom1(atom3))=2
               endif
            enddo
            enddo
         enddo
      enddo
   endif
enddo

return
end function detect_rings4

function detect_doubles(mol) result(cmatd) 
implicit none
type(molecule),intent(in)           :: mol
integer :: cmat(mol%natom,mol%natom)
integer :: cmatd(mol%natom,mol%natom)
integer :: bonded(mol%natom)
integer :: n,m,l,k !counters
integer :: nbonds 
logical :: cf,co,oco,amide,mkdbl,lexist
integer :: istat
character*10      ::func_group
character*2       ::fgrp

bonded=0
cmat=conn_matrix(mol)
cmatd=cmat

inquire(file='IN_PRESERVE_BONDS',exist=lexist)
if(lexist .eqv. .false.)then
   amide=.false.
   oco=.false.
   co=.true.
   cf=.false.
else
   amide=.false.;cf=.false.;oco=.false.;co=.false.
   open(unit=20,file='IN_PRESERVE_BONDS',status='unknown',iostat=istat)
   do
      read(20,*,iostat=istat)func_group,mkdbl
      if(istat<0) exit
      fgrp=func_group(1:2)
      if(fgrp=='am'.or.fgrp=='Am'.or.fgrp=='AM')then
         amide=mkdbl
      elseif(fgrp=='CF'.or.fgrp=='cf')then
         cf=mkdbl
      elseif(fgrp=='CO'.or.fgrp=='co')then
         co=mkdbl
      elseif(fgrp=='OC'.or.fgrp=='oc')then
         oco=mkdbl
      endif
   enddo
endif


!function get_bonded(atom_no,cmat) result(nbonded) 
!Detect all C-F bonds
if(cf .eqv. .true.)then
   do n=1,mol%natom
      if(mol%atoms(n)%label=='C')then
         bonded=get_bonded(n,cmat)
         do m=1,mol%natom
         if(bonded(m)/=0)then
            if((mol%atoms(bonded(m))%label)=='F')then
               cmatd(n,bonded(m))=2;cmatd(bonded(m),n)=2
            endif
          endif
         enddo
      endif
   enddo
endif

!Detect C=O
mol_loop: do n=1,mol%natom
   if(mol%atoms(n)%label=='C')then
      bonded=get_bonded(n,cmat)
      bonded_loop: do m=1,mol%natom
         if((bonded(m)/=0).and.(mol%atoms(bonded(m))%label)=='O')then
            nbonds=get_conn_num(bonded(m),cmat) !if O is only bonded to the C above, it must be a carbonyl
            !CO_loop: do l=1,mol%natom
               if(nbonds==1)then
                  if(co .eqv. .true.)then
                     cmatd(n,bonded(m))=2;cmatd(bonded(m),n)=2 !sets the C=O bond to double
                     amide_test: if(amide .eqv. .true.)then
                     ON_loop: do k=1,mol%natom !Now test for a N or O attached to the carbon ie NC=O or OC=O
                        if(k==m) cycle
                           !if((bonded(k)/=0).and.((mol%atoms(bonded(k))%label=='O').or.(mol%atoms(bonded(k))%label=='N')))then
                           if((bonded(k)/=0).and.((mol%atoms(bonded(k))%label=='N')))then
                              cmatd(n,bonded(k))=2;cmatd(bonded(k),n)=2
                              exit
                           endif
                     enddo ON_loop
                     endif amide_test
                     oco_test: if(oco .eqv. .true.)then
                     OO_loop: do k=1,mol%natom !Now test for a N or O attached to the carbon ie NC=O or OC=O
                        if(k==m) cycle
                           if((bonded(k)/=0).and.((mol%atoms(bonded(k))%label=='O')))then
                              cmatd(n,bonded(k))=2;cmatd(bonded(k),n)=2
                              exit
                           endif
                     enddo OO_loop
                     endif oco_test
                  endif
               endif
            !enddo CO_loop
         endif
      enddo bonded_loop
   endif
enddo mol_loop

return
end function detect_doubles





function num_breaks(cmat) result(nbreaks)

implicit none
integer,dimension(:,:),intent(in)   ::cmat
integer        :: n,m !counter
integer        :: natom
integer        :: nbreaks

nbreaks=0
natom=size(cmat,1)

do n=1,natom
   do m=n,natom
!      if(cmat(n,m)==1)then
! patch 200513 allow H-bonds
      if(abs(cmat(n,m))==1)then
         nbreaks=nbreaks+1
      endif
   enddo
enddo

return
end function num_breaks

function get_breaks(cmat) result(breaklist)

implicit none
integer,dimension(:,:),intent(in)   ::cmat
integer        :: n,m,nb !counter
integer        :: natom
integer        :: nbreaks
integer,allocatable   ::breaklist(:,:)

nbreaks=num_breaks(cmat)
natom=nint(sqrt(real(size(cmat))))
allocate(breaklist(2,nbreaks))

nb=0
do n=1,natom
   do m=n,natom
!      if(cmat(n,m)==1)then
! patch 200513 allow H-bonds
      if(abs(cmat(n,m))==1)then
         nb=nb+1
         breaklist(1,nb)=n
         breaklist(2,nb)=m
      endif
   enddo
enddo



return
end function get_breaks

subroutine make_groups(mol,groups_array,nag,glist,caps_array,ncapg,num_l0,cmatd)

implicit none

type(molecule),intent(in)           :: mol
integer,optional :: cmatd(:,:)
integer :: cmat(mol%natom,mol%natom)
!integer :: cmatd_local(mol%natom,mol%natom)  !use this one if cmatd is not supplied
integer :: A(mol%natom,mol%natom)
integer              :: nag(:) !num_atoms_group
integer              :: ncapg(:) !num_atoms_group
integer              :: nha    !num_heavy_atoms
integer              :: groups_array(:,:)
integer              :: glist(:) !group list - what group is each atom in?
type(two_c_bond)           :: caps_array(:,:)
integer              :: num_groups
logical              :: assigned(mol%natom)
integer              :: k,l,m,n,r !counters
integer              :: grp,atom1,atom2,row1,grp1,grp2,cap
integer,allocatable  :: GCM_tmp(:,:)
integer              :: num_bonded
integer              :: bonded(mol%natom) 
integer              :: num_l0 !number of actual groups
logical              :: modify !for use with GCM_tmp

num_l0=0
caps_array(:,:)%a1=0
caps_array(:,:)%a2=0
write(*,*)'In make_groups'
do n=1,mol%natom
   write(*,*)mol%atoms(n)%orig_atom_no1
enddo

if(present(cmatd))then
   write(*,*)'cmatd provided'
   cmat=conn_matrix(mol)
   A=cmatd-cmat
else
   write(*,*)'in else block'
   cmat=conn_matrix(mol)
   cmatd=detect_doubles(mol)
   cmatd=detect_rings3(cmatd)
   cmatd=detect_rings4(cmatd)
   A=cmatd-cmat
endif

write(*,*)'up to here'

!cmat=conn_matrix(mol)
!A=cmatd-cmat

write(*,*)'cmat='
do n=1,mol%natom
   write(*,'(30I3)')cmat(n,:)
enddo
write(*,*)
write(*,*)'cmatd='
do n=1,mol%natom
   write(*,'(30I3)')cmatd(n,:)
enddo
write(*,*)
write(*,*)'A='
do n=1,mol%natom
   write(*,'(30I3)')A(n,:)
enddo
write(*,*)

nha=mol%natom
!num_groups=num_breaks(cmatd)+1
num_groups=nha
write(*,*)'num_groups= ',num_groups

!allocate(nag(num_groups))
!allocate(groups_array(num_groups,mol%natom))

groups_array=0
nag=1 !set to 1 and then decrement all values at the end
ncapg=0
glist=0
assigned=.false.

!write(*,*)
!write(*,*)assigned
!write(*,*)
groups_loop: do grp=1,num_groups
   !nag(grp)=1
   !write(*,*)'Line 941 nag= ',nag

!Assign first non-assigned atom to be the first atom in the group
   first_atom_loop: do atom1=1,nha
      write(*,*)'grp= ',grp, ' atom1= ',atom1
      if(nag(grp)==1)then !only executes first time
         if(assigned(atom1) .eqv. .false.)then
            num_l0=num_l0+1 ! we are making a new l0 frag (group)
            groups_array(grp,nag(grp))=mol%atoms(atom1)%orig_atom_no1
            glist(mol%atoms(atom1)%orig_atom_no1)=grp
            nag(grp)=nag(grp)+1
            assigned(atom1)=.true.
            write(*,*)'Assigned atom',atom1,"FAL"
            call make_caps(atom1,grp,mol,cmatd,caps_array,ncapg)
         ! Now assign anything else that should be assigned (i.e. anything double bonded to atom1)
            row_loop: do row1=1,nha 
               if(A(row1,atom1)>=1)then
                  groups_array(grp,nag(grp))=mol%atoms(row1)%orig_atom_no1
                  glist(mol%atoms(row1)%orig_atom_no1)=grp
                  nag(grp)=nag(grp)+1
                  assigned(row1)=.true.
                  write(*,*)'Assigned atom',row1,"ROW"
                  call make_caps(row1,grp,mol,cmatd,caps_array,ncapg)
               endif
            enddo row_loop
         endif
      endif
   enddo first_atom_loop
   ! Looping over anything assigned, to get chains of double bonds
   find_chain_loop: do l=1,nha !originally l=2, changed to 1 to find if one thing is double bonded to two others.
      atom2=orig_atom(mol,groups_array(grp,l))
      !atom2=groups_array(grp,l) !This line will produce the wrong thing if all the hevy atoms aren't first
      if(atom2/=0)then
         row_loop2: do k=1,nha
         !write(*,'(30I3)')A(atom2,:)
         if((A(atom2,k)>=1).and.(assigned(k) .eqv..false.))then
            write(*,*)'Line 983: grp=',grp,' atom1= ',atom1,' atom2= ',atom2,' group= ',groups_array(grp,:) 
            groups_array(grp,nag(grp))=mol%atoms(k)%orig_atom_no1
            glist(mol%atoms(k)%orig_atom_no1)=grp
            assigned(k)=.true.
            write(*,*)'Assigned atom',k,"ROW2"
            nag(grp)=nag(grp)+1
            call make_caps(k,grp,mol,cmatd,caps_array,ncapg)
         endif
      enddo row_loop2
      endif
      enddo find_chain_loop
   !enddo
enddo groups_loop

!Subtract 1 from each num_atom_group as we have incremented it once too many
nag=nag-1

write(42,*)assigned
if(any(assigned .eqv. .false.))then
   write(*,*)'Not all atoms assigned. This is probably due to molecule parts not being connected'
   stop
endif

write(*,*)'finished make_groups'
write(*,*)'groups_array='
do n=1,num_groups
   write(*,'(I3,A,30I3)')n,': ',groups_array(n,1:nag(n))
   !write(*,*)n,': ',groups_array(n,:)
enddo

!Check that two groups aren't stuck to each other at both ends
!Set up a temp GCM
allocate(GCM_tmp(num_groups,num_groups))
call make_GCM_tmp(mol,glist,cmatd,GCM_tmp)

write(42,*)'tmp_GCM'
do n=1,num_groups
   write(42,'(15I3)')GCM_tmp(:,n)
enddo
write(42,*)'end tmp_GCM'

write(42,*)'finished make_groups Hi!'
write(42,*)'groups_array='
do n=1,num_groups
   write(42,'(I3,A,30I3)')n,': ',groups_array(n,1:nag(n))
   !write(*,*)n,': ',groups_array(n,:)
enddo

modify=.true.
!Concatenate any groups that are doubly stuck
test_loop: do while (modify .eqv. .true.)
   do grp1=1,num_groups
      do grp2=grp1,num_groups
         if(GCM_tmp(grp1,grp2)>2)then
            groups_array(grp1,nag(grp1)+1:nag(grp1)+nag(grp2))=groups_array(grp2,1:nag(grp2))
            nag(grp1)=nag(grp1)+nag(grp2)
            caps_array(grp1,ncapg(grp1)+1:ncapg(grp1)+ncapg(grp2))=caps_array(grp2,1:ncapg(grp2))
            ncapg(grp1)=ncapg(grp1)+ncapg(grp2)
            glist(groups_array(grp2,1:nag(grp2)))=grp1
            groups_array(grp2:num_groups-1,:)=groups_array(grp2+1:num_groups,:)
            nag(grp2:num_groups-1)=nag(grp2+1:num_groups)
            caps_array(grp2:num_groups-1,:)=caps_array(grp2+1:num_groups,:)
            ncapg(grp2:num_groups-1)=ncapg(grp2+1:num_groups)
            where(glist>grp2)glist=glist-1
            num_l0=num_l0-1 ! decrement number of groups
            call make_GCM_tmp(mol,glist,cmatd,GCM_tmp) ! update GCM_tmp 
            cycle test_loop
         endif
      enddo
   enddo
   modify=.false. !We've made it through without needing to modify groups
enddo test_loop
write(42,*)'Actual number of groups is',num_l0
!Zero stuff just in case
groups_array(num_l0+1:num_groups,:)=0
nag(num_l0+1:num_groups)=0

write(42,*)'corrected make_groups'
write(42,*)'groups_array='
do n=1,num_groups
   write(42,'(I3,A,I3,A,30I3)')n,': ',nag(n),':',groups_array(n,:)
   !write(*,*)n,': ',groups_array(n,:)
enddo
write(42,*)'glist='
write(42,*)glist

write(42,*)'caps_array='
do n=1,num_groups
   do m=1,ncapg(n)
   write(42,'(I3,A,2I3)')n,': ',caps_array(n,m)%a1,caps_array(n,m)%a2
   !write(*,*)n,': ',groups_array(n,:)
   enddo
enddo

!Fix caps
do grp=1,num_l0
   !do cap=1,ncapg(grp)
   cap=1
   do while (caps_array(grp,cap)%a1>0)
      grp1=glist(caps_array(grp,cap)%a1)
      grp2=glist(caps_array(grp,cap)%a2)
      write(42,*)grp,cap,'|',caps_array(grp,cap)%a1,grp1,'|',caps_array(grp,cap)%a2,grp2
      if(grp1==grp2)then !delete the cap
         caps_array(grp,cap:ncapg(grp)-1)=caps_array(grp,cap+1:ncapg(grp))
         ncapg(grp)=ncapg(grp)-1
         caps_array(grp,ncapg(grp)+1)%a1=0
         cap=cap-1
      endif
      cap=cap+1
   enddo
enddo

write(42,*)'fixed caps_array='
do n=1,num_groups
   do m=1,ncapg(n)
   write(42,'(I3,A,2I3)')n,': ',caps_array(n,m)%a1,caps_array(n,m)%a2
   !write(*,*)n,': ',groups_array(n,:)
   enddo
enddo


write(*,*)'nearly outta make_groups'
deallocate(GCM_tmp)
write(*,*)'outta make_groups'
return
end subroutine make_groups

subroutine make_GCM_tmp(mol,glist,cmatd,GCM_tmp)

type(molecule),intent(in)     ::mol !heavy mol
integer,intent(in)            ::glist(:)
integer,intent(in)            ::cmatd(:,:)

integer                       ::GCM_tmp(:,:)

integer                       ::atom1,atom2,grp1,grp2,nha
integer                       ::num_bonded
integer                       ::bonded(mol%natom)

nha=mol%natom

GCM_tmp=0
do atom1=1,nha
   num_bonded=get_conn_num(atom1,cmatd)
   bonded=get_bonded(atom1,cmatd)
   write(42,*)atom1,':',bonded(1:num_bonded)
   do atom2=1,num_bonded
      grp1=glist(mol%atoms(atom1)%orig_atom_no1)
      grp2=glist(mol%atoms(bonded(atom2))%orig_atom_no1)
      write(42,*)mol%atoms(atom1)%orig_atom_no1,'is in',grp1
      write(42,*)mol%atoms(bonded(atom2))%orig_atom_no1,'is in',grp2
      if(grp1/=grp2)then
         GCM_tmp(grp1,grp2)=GCM_tmp(grp1,grp2)+1
         GCM_tmp(grp2,grp1)=GCM_tmp(grp2,grp1)+1
      endif
   enddo
enddo

return
end subroutine make_GCM_tmp


subroutine group_conn_mat(mol,groups_array,nhg,glist,gcm)

implicit none

type(molecule),intent(in)     ::mol
integer,intent(in)              :: nhg(:) !num_atoms_group. Heavies only
integer,intent(in)              :: groups_array(:,:)
integer,intent(in)              :: glist(:) !group list - what group is each atom in?
integer              :: cmat(mol%natom,mol%natom) 
integer        ::bonded(mol%natom)
integer        ::grp,atom1,atom2,n,m,ngroups

integer        ::gcm(:,:)
!integer        ::gcm_shadow(:,:)
   
!gcm=.false.
gcm=0
!gcm_shadow=0

cmat=conn_matrix(mol)
ngroups=size(groups_array,1)

!write(*,*)'In GCM, cmat='
!do n=1,ngroups
!   write(*,*)cmat(:,n)
!enddo
!
!write(*,*)'In GCM, glist='
!write(*,*)glist


do grp=1,ngroups
   do n=1,nhg(grp)
      atom1=groups_array(grp,n)
      !write(*,*)'atom1=',atom1
      bonded=get_bonded(atom1,cmat)
      do m=1,mol%natom
         atom2=bonded(m)
         if(atom2/=0)then
      !write(*,*)'atom2=',atom2
            if(glist(atom2)/=grp.and.glist(atom2)/=0)then  !if it's connected to atom1 and not in the same group, then it is in a neighbouring grp
               !write(*,*)grp,glist(atom2)
               gcm(grp,glist(atom2))=1
               gcm(glist(atom2),grp)=1
!               gcm_shadow(grp,glist(atom2))=1
!               gcm_shadow(glist(atom2),grp)=1
            endif
         endif
      enddo
   enddo
enddo

end subroutine group_conn_mat

subroutine make_caps(atom_no,group,mol,cmatd,caps_array,ncapg)
implicit none

type(molecule),intent(in)           :: mol
integer :: cmatd(:,:)
integer              :: ncapg(:) !num_caps_group
integer              :: atom_no,group    
type(two_c_bond)           :: caps_array(:,:)
integer              :: num_groups
integer              :: num_atoms
logical              :: assigned(mol%natom)
integer              :: k,l,m,n,r !counters
integer              :: ncap,atom1,atom2,row1

num_groups=(size(caps_array,1))
num_atoms=(size(cmatd,1))


do n=1,num_atoms
   if(cmatd(atom_no,n)==1)then
      caps_array(group,ncapg(group)+1)%a1=mol%atoms(atom_no)%orig_atom_no1
      caps_array(group,ncapg(group)+1)%a2=mol%atoms(n)%orig_atom_no1
      ncapg(group)=ncapg(group)+1
      !write(*,*)mol%atoms(atom_no)%orig_atom_no1,mol%atoms(n)%orig_atom_no1,group,ncapg(group)
   endif
enddo

return
end subroutine make_caps

subroutine cap_charges(mol)
implicit none 

type(molecule)       ::mol
real*8               ::real_charge,cap_charge
integer              ::atom1,num_caps

num_caps=0
real_charge=0

do atom1=1,mol%natom
   if(mol%atoms(atom1)%real_or_not .eqv. .true.)then
      real_charge=real_charge+mol%atoms(atom1)%Charge
   else
      num_caps=num_caps+1
   endif
enddo

cap_charge=-(real_charge/num_caps)

do atom1=1,mol%natom
   if(mol%atoms(atom1)%real_or_not .eqv. .false.)then
      mol%atoms(atom1)%Charge=cap_charge
   endif 
enddo

return
end subroutine cap_charges


subroutine cap_charges2(mol,group,caps_array,num_caps)
implicit none
type(molecule),intent(in)     ::mol
type(molecule)                ::group
integer                       ::num_caps
type(two_c_bond)              :: caps_array(:)
real*8                        ::OA_charges(num_caps) !OrigAtom charges
real*8                        ::capH_charges(num_caps)
real*8                        ::group_charge,cap_charge,OA_sum
integer                       ::n,m   !counters


group_charge=0
capH_charges=0
OA_charges=0


do n=1,num_caps
   OA_charges(n)=mol%atoms(caps_array(n)%a2)%Charge
enddo
write(*,*)'OA q=', OA_charges
OA_sum=sum(OA_charges)

do n=1,group%natom
   if(group%atoms(n)%real_or_not .eqv. .true.)then
      group_charge=group_charge+group%atoms(n)%Charge
   endif
enddo
cap_charge=-group_charge
!write(*,*)cap_charge

capH_charges=(OA_charges/OA_sum)*cap_charge

m=1
write(19,*)capH_charges
do n=1,group%natom
   if(group%atoms(n)%real_or_not .eqv. .false.)then
      group%atoms(n)%Charge=capH_charges(m)
      m=m+1
   endif
enddo

return
end subroutine cap_charges2


subroutine cap_charges_even(mol,group,caps_array,num_caps)
implicit none
type(molecule),intent(in)     ::mol
type(molecule)                ::group
integer                       ::num_caps
type(two_c_bond)              :: caps_array(:)
real*8                        ::OA_charges(num_caps) !OrigAtom charges
real*8                        ::capH_charges(num_caps)
real*8                        ::group_charge,cap_charge,OA_sum
integer                       ::n,m   !counters


group_charge=0
capH_charges=0
OA_charges=0

do n=1,group%natom
   if(group%atoms(n)%real_or_not .eqv. .true.)then
      group_charge=group_charge+group%atoms(n)%Charge
   endif
enddo
cap_charge=-group_charge
!write(*,*)cap_charge

capH_charges=cap_charge/num_caps

m=1
write(19,*)capH_charges
do n=1,group%natom
   if(group%atoms(n)%real_or_not .eqv. .false.)then
      group%atoms(n)%Charge=capH_charges(m)
      m=m+1
   endif
enddo

return
end subroutine cap_charges_even


subroutine addH_groups(mol,groups_array,nag,glist)

implicit none

type(molecule),intent(in)           :: mol
integer :: bondedH(mol%natom)
integer              :: nag(:) !num_atoms_group
integer              :: groups_array(:,:)
integer              :: glist(:)
integer,dimension(size(groups_array,1),size(groups_array,2))  ::heavy_array
integer              :: num_groups
integer              :: k,l,m,n,r !counters
integer              :: grp,Hatom,atom1,atom2,row1

heavy_array=groups_array
num_groups=size(groups_array,1)

group_loop: do grp=1,num_groups
   atom_loop: do atom1=1,mol%natom
      if(heavy_array(grp,atom1)/=0)then
         bondedH=get_bondedH(heavy_array(grp,atom1),mol)
         write(*,*)grp,': atom= ',heavy_array(grp,atom1), ': bondedH=', bondedH
         H_loop: do Hatom=1,mol%natom
         if(bondedH(Hatom)/=0)then
            groups_array(grp,nag(grp)+1)=bondedH(Hatom)
            glist(bondedH(Hatom))=grp
            nag(grp)=nag(grp)+1
         endif
         enddo H_loop
      endif
      enddo atom_loop
   enddo group_loop
         


return
end subroutine addH_groups

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 07/08/09

function count_heavy_atoms(mol_in) result(num_heavy)
implicit none

type(molecule),intent(in)        :: mol_in
integer                          :: num_heavy

num_heavy=0
num_heavy=count(mol_in%atoms%label/='H')


return
end function count_heavy_atoms


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
real*8 function anum2amass(na)
  
 integer na
  
 real*8 atomic_masses(110)
  
 data atomic_masses/ &
    1.0079,   4.0026,   6.9410,   9.0122,  10.8110,  12.0107, &
   14.0067,  15.9994,  18.9984,  20.1797,  22.9898,  24.3050, &
   26.9815,  28.0855,  30.9738,  32.0650,  35.4530,  39.9480, &
   39.0983,  40.0780,  44.9559,  47.8670,  50.9415,  51.9961, &
   54.9380,  55.8450,  58.6934,  58.9332,  63.5460,  65.3800, &
   69.7230,  72.6400,  74.9216,  78.9600,  79.9040,  83.7980, &
   85.4678,  87.6200,  88.9059,  91.2240,  92.9064,  95.9600, &
   98.0000, 101.0700, 102.9055, 106.4200, 107.8682, 112.4110, &
  114.8180, 118.7100, 121.7600, 127.6000, 126.9045, 131.2930, &
  132.9055, 137.3270, 138.9055, 140.1160, 140.9077, 144.2420, &
  145.0000, 150.3600, 151.9640, 157.2500, 158.9254, 162.5000, &
  164.9303, 167.2590, 168.9342, 173.0540, 174.9668, 178.4900, &
  180.9479, 183.8400, 186.2070, 190.2300, 192.2170, 195.0840, &
  196.9666, 200.5900, 204.3833, 207.2000, 208.9804, 210.0000, &
  210.0000, 220.0000, 223.0000, 226.0000, 227.0000, 231.0359, &
  232.0381, 237.0000, 238.0289, 243.0000, 244.0000, 247.0000, &
  247.0000, 251.0000, 252.0000, 257.0000, 258.0000, 259.0000, &
  262.0000, 261.0000, 262.0000, 266.0000, 264.0000, 277.0000, &
  268.0000, 271.0000/

 anum2amass=atomic_masses(na)

 return

 end function anum2amass

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


real*8 function anum2vdw(na)

integer na
real*8 vanderwaals_radii(110)

data vanderwaals_radii/ &
  2.2676727, 2.6456187, 3.4393037, 3.7794547, 3.7794547, 3.2125367, &
  2.9290777, 2.8723857, 2.7778997, 2.9101797, 4.2896807, 3.2692277, &
  3.7794547, 3.9684267, 3.4015087, 3.4015087, 3.3070227, 3.5526877, &
  5.1967497, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.0802557, 2.6456187, 2.6267207, &
  3.5337897, 3.7794547, 3.4959957, 3.5904817, 3.4959957, 3.8172487, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.0802557, 3.2503307, 2.9857687, &
  3.6471737, 4.1007077, 3.7794547, 3.8928377, 3.7416597, 4.0818107, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.2503307, &
  3.1369477, 2.9290777, 3.7038657, 3.8172487, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.5148927, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547/

 anum2vdw=vanderwaals_radii(na)

 return

 end function anum2vdw


!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!real*8 function anum2cov(na)
!
!integer :: na
!real*8 :: covalent_radii(110)
!
!data covalent_radii/  &
!! disagreement for covalent radii of H, some say 0.23 A = 0.43 bohr
!!                                       others   0.37 A = 0.70 bohr
!! 0.6999189 is the original H value, changed to 0.43 to agree with Mick
!! 07/10/08 MAA
!  0.43, 2.8345907, 1.2850147, 0.6614047, 1.5684737, 1.2850147,  &
!  1.2850147, 1.2850147, 1.2094257, 2.8345907, 1.8330357, 2.0787007,  &
!  2.5511317, 2.2676727, 1.9842137, 1.9275217, 1.8708307, 2.8534887,  &
!  2.5133377, 1.8708307, 2.7212077, 2.7778997, 2.5133377, 2.5511317,  &
!  2.5511317, 2.5322347, 2.5133377, 2.8345907, 2.8723857, 2.7401047,  &
!  2.3054677, 2.2109807, 2.2865707, 2.3054677, 2.2865707, 2.8345907,  &
!  2.7778997, 2.1164947, 3.3637147, 2.9479747, 2.7967967, 2.7778997,  &
!  2.5511317, 2.6456187, 2.7401047, 2.8345907, 3.0046667, 3.1936387,  &
!  3.0802557, 2.7590017, 2.7590017, 2.7778997, 2.6456187, 2.8345907,  &
!  3.1558447, 2.5322347, 3.5337897, 3.4582007, 3.4393037, 3.4204067,  &
!  3.4015087, 3.4015087, 3.7605567, 3.3826117, 3.3259197, 3.3070227,  &
!  3.2881257, 3.2692277, 3.2503307, 3.6660707, 3.2503307, 2.9668717,  &
!  2.7023097, 2.5889267, 2.5511317, 2.5889267, 2.4944397, 2.8345907,  &
!  2.8345907, 3.2125367, 2.9290777, 2.9101797, 2.9101797, 3.1747417,  &
!  2.2865707, 2.8345907, 2.8345907, 3.5904817, 3.5526877, 3.3826117,  &
!  3.0424607, 2.9857687, 2.9290777, 2.8912827, 2.8534887, 1.8708307,  &
!  2.9101797, 3.4582007, 2.8345907, 2.8345907, 2.8345907, 2.8345907,  &
!  2.8345907, 2.8345907, 2.8345907, 2.8345907, 2.8345907, 2.8345907,  &
!  2.8345907, 2.8345907/
!
!anum2cov=covalent_radii(na)
!
!return
!
!end function anum2cov
!
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

real*8 function anum2cov(na)

integer :: na
real*8 :: covalent_radii(110)

data covalent_radii/  &

! in Angstrom, from Mick's atomic.dat, other values converted from Chris
! Mick's atomic.dat misses 85-87 incl and from 95 onward

  0.23, 1.22, 0.68, 0.35, 0.83, 0.68, 0.68, 0.68, 0.64, 1.6,  &
  0.97, 1.1,  1.35, 1.2,  1.05, 1.02, 0.99, 1.92, 1.33, 0.99, &
  1.44, 1.47, 1.33, 0.67, 1.35, 1.34, 1.33, 1.5,  1.52, 1.45, &
  1.22, 1.17, 1.21, 1.22, 1.21, 1.98, 1.47, 1.12, 1.78, 1.56, & 
  1.48, 1.47, 1.35, 1.4,  1.45, 1.5,  1.59, 1.69, 1.63, 1.46, &
  1.46, 1.47, 1.4,  2.18, 1.67, 1.34, 1.87, 1.83, 1.82, 1.81, &
  1.8,  1.8,  1.99, 1.79, 1.76, 1.75, 1.74, 1.73, 1.72, 1.94, &
  1.72, 1.57, 1.43, 1.37, 1.35, 1.37, 1.32, 1.5,  1.5,  1.7,  &
  1.55, 1.54, 1.54, 1.68, 1.21, 1.50, 1.50, 1.9,  1.88, 1.79, &
  1.61, 1.58, 1.55, 1.53, 1.51, 0.99, 1.54, 1.83, 1.50, 1.50, &
  1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50/


anum2cov=covalent_radii(na)

return

end function anum2cov


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

integer function alab2anum(alab)

! given the atom label this returns the atomic number

character*2 :: alab,atomic_labels(110)

data atomic_labels/'H','He','Li','Be','B','C','N','O','F','Ne',      &
 'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',  &
 'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb',   &
 'Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',    &
 'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',    &
 'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os',    &
 'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',   &
 'Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',    &
 'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds'/

integer n

do n=1,110
   if(atomic_labels(n).eq.alab)then
      alab2anum=n
      return
   endif
enddo
alab2anum=0

return

end function alab2anum

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

function anum2alab(anum) result(alab)

! given the atomic number this returns the atom label

character*2 :: alab,atomic_labels(110)

data atomic_labels/'H','He','Li','Be','B','C','N','O','F','Ne',      &
 'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',  &
 'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb',   &
 'Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',    &
 'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',    &
 'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os',    &
 'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',   &
 'Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',    &
 'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds'/

integer anum

alab=atomic_labels(anum)

if(anum>110)then
   write(*,*)'Invalid atomic number, only up to 110 supported'
   alab=''
endif

return

end function anum2alab

subroutine checkcharge(ndim,nchgs,xchgs,mol,nchgtot)
      type(molecule),intent(in)   ::mol
      real*8                      ::x(3),a2b
      integer                     ::ndim,k,nchgtot,atom1,j
      real*8                      ::xchgs(ndim,3),xsum

      a2b=1.8897259886d0
      nchgtot=0

      do atom1=1,mol%natom

       do k=1,3
        x(k)=mol%atoms(atom1)%Coords(k)
       enddo

       do j=1,nchgs
       xsum=0.d0
        do k=1,3
         xsum=xsum+(x(k)*a2b-xchgs(j,k))**2
        enddo

       if(xsum.lt.0.1d0)then
        nchgtot=1
        return
       endif

       enddo

      enddo

      return
end subroutine checkcharge

end module molecule_adt
