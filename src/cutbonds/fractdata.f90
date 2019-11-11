      module fractdata
      implicit double precision(a-h,o-z)
      parameter(nsmall=120)
      parameter(ncomp=12)
      parameter(maxfamily=50)
      save

      character*80 icomm
      character*2 sym(91)
      character*80 inputformat
      character*20 ititle
      dimension numb(91)
      dimension rad(91),subtimes(10)

      character*80, allocatable :: abinitio(:)
      character*2, allocatable  :: atoms(:)
      integer, allocatable      :: numa(:)
      integer, allocatable      :: nstop(:),numat(:),natstore(:,:)
      integer, allocatable      :: isign(:),itype(:,:),ibond(:,:)
      integer, allocatable      :: ib1(:,:)

      integer, allocatable      :: ilink(:,:),map(:,:),kf(:)
      integer, allocatable      :: nfam(:),ifam(:,:),mult(:),ichg(:)
      integer, allocatable      ::  mbstore(:),nbstore(:),multstore(:)
      integer, allocatable     :: junk(:),itotf(:)

      integer, allocatable      :: nb(:),mb(:),newmult(:),matb(:,:)

      integer, allocatable      :: mats(:),matr(:),itp(:),ibo(:,:)

      integer, allocatable      :: ma(:,:),notthere(:),numgroups(:)

      real*8, allocatable       :: coords(:,:),radius(:),c(:,:)

      integer, allocatable      :: ibf(:,:,:),itf(:,:),ngpf(:,:,:)

      integer, allocatable      :: kat(:),mbn(:),nbn(:),multn(:)

      integer newfrag,Level,natom,natomall,nbondso,nbonds,nffinal
      integer itfinal,icycle,nfrag,nafrag,nf,nbig3,nfstore,nhstore
      integer nocapsatall,nelim,nabitio,iterfinal,nfragm,numhbonds
      integer nbondsextra

      real*8 dtol

      end module fractdata
