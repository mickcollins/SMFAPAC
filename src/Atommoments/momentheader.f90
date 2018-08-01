      module momentheader

      implicit double precision(a-h,o-z)

      parameter(maxcaps=20)

      real*8, allocatable   :: c(:,:),ch(:),d(:,:),q(:,:,:)
      real*8, allocatable   :: o(:,:,:,:),h(:,:,:,:,:)

      real*8, allocatable   :: w1(:,:,:)

      integer, allocatable  :: nat1(:),num1(:,:),natnum1(:,:,:)
      integer, allocatable  :: isign(:)

      character*2, allocatable :: lab(:)

      integer nfrag,nL1,natom,nflag,maxatom

      real*8 bohr

      end module momentheader
