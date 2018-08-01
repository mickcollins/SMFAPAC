      module derivheader

      implicit double precision(a-h,o-z)

      parameter(maxcaps=20)

      real*8, allocatable   :: c(:,:,:),force(:,:,:),fc(:,:,:)
      real*8, allocatable   :: amas(:),coord(:,:),energyfrag(:)
      real*8, allocatable   :: fact(:,:)

      integer, allocatable  :: natomfrag(:),nat0(:),nat(:,:)
      integer, allocatable  :: ncap(:),natcap(:,:,:),isign(:)
      integer, allocatable  :: numstore(:,:)

      character*2, allocatable  :: lab(:)

      integer nfrag,natom,nflag

      character*80 comlin

      real*8 bohr

      end module derivheader
