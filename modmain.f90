module modmain
use modconstants

integer, parameter :: nrmax=3
real(8), parameter :: a=2.46d0*astr
real(8), parameter :: a_1(3)=a*(/1d0,0d0,0d0/)
real(8), parameter :: a_2(3)=a*(/1d0/2d0,dsqrt(3d0)/2d0,0d0/)
real(8), parameter :: b_1(3)=(2*pi/a)*(/1d0,-1d0/dsqrt(3d0),0d0/)
real(8), parameter :: b_2(3)=(2*pi/a)*(/0d0, 2d0/dsqrt(3d0),0d0/)
real(8), parameter :: tau_1_vec(3)=(2d0*a_2-a_1)/3d0
real(8), parameter :: nnv(3,3)=transpose(reshape(a*(/0.d0     ,  1.d0/dsqrt(3.d0)      ,0.d0&
                                                    ,1.d0/2.d0, -1.d0/2.d0/dsqrt(3.d0) ,0.d0&
                                                   ,-1.d0/2.d0, -1.d0/2.d0/dsqrt(3.d0) ,0.d0/),(/3,3/)))
real(8), parameter :: v_pppi_0=-2.7d0*eV
real(8), parameter :: v_ppsgm_0=0.48d0*eV
real(8), parameter :: dlt_0=0.184d0*a
real(8), parameter :: a_0=a/dsqrt(3d0)
real(8), parameter :: d_0=3.35d0*astr

real(8), parameter :: t_1=2.8*eV
! LINEBANDPLOT
integer :: nkline, nhspts
integer, allocatable :: n_p_to_p_grid(:)
real(8), allocatable :: hspts(:,:), kline(:)

! LAPACK
integer :: lwork,lainfo
complex(8), allocatable :: work(:)
real(8), allocatable :: rwork(:)

! PARAMETERS FOR SYSTEM
integer, parameter :: norb=2
integer, parameter :: iorbstart=1
logical :: lineplotmode=.false.
integer, parameter :: npx=125, npy=125, nenrg=1
real(8) :: dt=0.2d0
integer :: itmax=0
integer :: exitmax=0
real(8) :: atnum(norb)

! Final State parameter
integer, parameter :: ngmax=2
integer, parameter :: ngtotax=2*ngmax+1
integer, parameter :: ngtot=ngtotax**2

real(8) :: g_vecs(ngtot,3)

! ARRAYS
real(8) :: orb_pos(norb,3)
real(8) :: pxg(npx),pyg(npy),enrg(nenrg)
integer ::  nppt
real(8) :: efermi=0
real(8), allocatable :: p_vecs(:,:)
complex(8), allocatable :: eigvecs_0(:,:,:)
complex(8), allocatable :: eigvecs_t(:,:,:)
real(8), allocatable :: eigvals(:,:)
complex(8), allocatable :: grd_hk(:,:,:)

!Berry curvature calculation
logical, parameter :: BerryPlot=.false.
end module
