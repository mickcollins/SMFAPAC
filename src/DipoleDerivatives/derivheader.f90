      module derivheader

      implicit double precision(a-h,o-z)

      parameter(maxcaps=20)

      real*8, allocatable   :: dipdr(:,:,:),dipdrfrag(:,:,:)
      real*8, allocatable   :: fact(:,:)

      integer, allocatable  :: natomfrag(:),nat0(:),nat(:,:)
      integer, allocatable  :: ncap(:),natcap(:,:,:),isign(:)

      integer nfrag,natom,nflag,maxatom

      real*8 bohr

      end module derivheader
