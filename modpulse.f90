module modpulse
use modconstants
implicit none
real(8) :: t_pump_start
real(8) :: t_pump_end
real(8) :: t_error

real(8), parameter :: w_pump=3*eV
real(8), parameter :: e_to_v_pump=eV/astr*137/w_pump ! E-field(eV/astr) to A-field(a.u.) unit conversion factor
real(8), parameter :: afield_vec_pump(3)=(/0.0d0,0.0d0,0.d0/)*e_to_v_pump
real(8), parameter :: phase_pump=0.d0
real(8), parameter :: hd_pump=6*fs ! half-duration

real(8), parameter :: w_probe=100*eV
real(8), parameter :: e_to_v_probe=eV/astr*137/w_probe ! E-field to A-field unit conversion factor
real(8), parameter :: afield_vec_probe(3)=(/0.00000d0,0.00005d0,0.00000d0/)*e_to_v_probe
real(8), parameter :: phase_probe=0.d0
real(8), parameter :: hd_probe=40*fs ! half-duration

real(8) :: delta_t=0*fs ! delta_t=probe pusle center - pump pulse center

real(8) :: tshift

real(8), allocatable :: pump(:,:),probe(:,:) !timeindex,vectorindex
real(8), allocatable :: e_pump(:,:),e_probe(:,:) !timeindex,vectorindex

! Circular polar inputs ::
logical, parameter :: circular_pump=.true.
real(8), parameter :: w_cpump=2*eV
real(8), parameter :: e_to_v_cpump=eV/astr*137/w_cpump ! E-field(eV/astr) to A-field(a.u.) unit conversion factor
real(8), parameter :: afield_intensity_cpump=0.0d0*e_to_v_cpump

logical, parameter :: circular_probe=.true.
real(8), parameter :: w_cprobe=35*eV
real(8), parameter :: e_to_v_cprobe=eV/astr*137/w_cprobe ! E-field(eV/astr) to A-field(a.u.) unit conversion factor
real(8), parameter :: afield_intensity_cprobe=0.00005d0*e_to_v_cprobe
end module
