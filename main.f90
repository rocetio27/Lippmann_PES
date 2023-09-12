program main
  use modmain
  use modmpi
  use modomp
  use modpulse
  use modconstants
  implicit none
  integer :: i1,i2,ie,ipx,ipy,it,itex,i,j,ii,ip,iz
  real(8) :: er, ec, origin_vec(3), width, height
  real(8) :: d_vec(3)
  real(8) :: d
  real(8) :: p_vec(3)
  complex(8) :: pes_tintegrand_vec(3)
  complex(8) :: pes_vec(norb,nenrg,3)
  complex(8) :: pes_rintegral_vecs(nenrg,norb,3)
  complex(8) :: pes_temp
  complex(8) :: pes_rcp_temp
  complex(8) :: pes_lcp_temp
  real(8), allocatable :: pes2(:,:)
  real(8), allocatable :: pes2_rcp(:,:)
  real(8), allocatable :: pes2_lcp(:,:)
  real(8) :: phase(nenrg,norb)
  real(8) :: phase2
  real(8) :: temp
  real(8) :: percent
  complex(8) :: phasefactor
  integer :: nmpi
  integer :: nthd, ompthnum
  complex(8) :: hk_t_1(norb,norb), hk_t_2(norb,norb), hk_t_3(norb,norb)
  complex(8) :: temp_grpvel_mat(norb,norb,3)
  real(8), allocatable :: berryC(:,:)
  real(8) :: vbmax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!MPI AND OMP INITIALIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call mpi_init(ierror)
call mpi_comm_dup(mpi_comm_world,mpicom,ierror)
call mpi_comm_size(mpicom,np_mpi,ierror)
call mpi_comm_rank(mpicom,lp_mpi,ierror)
if (lp_mpi.eq.0) then
  mp_mpi=.true.
else
  mp_mpi=.false.
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!MPI AND OMP INITIALIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

atnum=6d0  
call get_orb_pos(orb_pos)
 
!G VECTOR GENERATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i1=-ngmax,ngmax
do i2=-ngmax,ngmax
  g_vecs((i1+ngmax)*ngtotax+(i2+ngmax+1),:)=i1*b_1+i2*b_2
enddo
enddo
!==========================================================

!K AND E POINT GENERATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------
if(lineplotmode.eqv.(.true.)) then
!----------------------------------------------------------
    nhspts=2 !!
    allocate(hspts(nhspts,3)) !3 means px,py,pz
    allocate(n_p_to_p_grid(nhspts-1))

    hspts(1,:)=0.66666666d0*b_1+0.33333333d0*b_2-(/0.0d0, 0.12d0, 0.0d0/) !!
    hspts(2,:)=0.66666666d0*b_1+0.33333333d0*b_2+(/0.0d0, 0.12d0, 0.0d0/) !!
    n_p_to_p_grid=(/250/) !!
    nkline=sum(n_p_to_p_grid)+1
    nppt=nkline

    allocate(p_vecs(nppt,3))
    allocate(kline(nppt))
    call genkline(hspts,p_vecs,kline)
    !ENERGY----------------------------------------------
    er=1.2d0*eV                       !!energy range
    if(circular_probe.eqv.(.true.)) then
      ec=abs(w_cprobe)-0.4d0*eV              !!energy center
      call gengrid1d(ec-er/2d0,ec+er/2d0,nenrg,enrg)
    elseif(circular_probe.eqv.(.false.)) then
      ec=w_probe                    !!energy center
      call gengrid1d(ec-er/2d0,ec+er/2d0,nenrg,enrg)
    endif
!--------------------------------------------------------
elseif(lineplotmode.eqv.(.false.)) then
!--------------------------------------------------------
    nppt=npx*npy
    allocate(p_vecs(nppt,3))

    origin_vec=0.66666666d0*b_1+0.33333333*b_2 !!
    write(*,*) "(info) origin_vec:", origin_vec
    width=0.22d0                               !!
    height=0.22d0                              !!
    call gengrid1d(origin_vec(1)-width/2d0,origin_vec(1)+width/2d0,npx,pxg)
    call gengrid1d(origin_vec(2)-height/2d0,origin_vec(2)+height/2d0,npy,pyg)
    do ipy=1,npy
    do ipx=1,npx
      p_vecs(npx*(ipy-1)+ipx,:)=(/pxg(ipx),pyg(ipy),0.d0/)
    enddo
    enddo

    !ENERGY----------------------------------------------
    if(circular_probe.eqv.(.true.)) then
      enrg=abs(w_cprobe)-0.4d0*eV
    elseif(circular_probe.eqv.(.false.)) then
      enrg=w_probe-0.4d0*eV
    endif
!----------------------------------------------------------
endif
!==========================================================

!PUMP-PROBE SETTINGS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if((-1d0*hd_pump).le.(delta_t-hd_probe)) then
    tshift=-1d0*hd_pump !coordinate shift to fit the left-side first
  else
    tshift=delta_t-hd_probe !coordinate shift to fit the left-side first
  endif
  itmax=ceiling((max(-tshift-hd_pump,-tshift+hd_pump,delta_t-tshift-hd_probe,delta_t-tshift+hd_probe)&
       -min(-tshift-hd_pump,-tshift+hd_pump,delta_t-tshift-hd_probe,delta_t-tshift+hd_probe))/dt)+1
  t_pump_start=-tshift-hd_pump
  t_pump_end=-tshift+hd_pump
  t_error=ceiling(t_pump_start/dt)*dt-t_pump_start

   exitmax=2*itmax-1
   allocate(pump(exitmax,3),probe(exitmax,3))
   allocate(e_probe(exitmax,3))
   call genafield_pump(pump)
   call genafield_probe(probe)
   call genefield_probe(e_probe)
!==========================================================

!GROUND STATE DIAGONALIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(eigvecs_0(nppt,norb,norb))
    allocate(eigvecs_t(nppt,norb,norb))
    allocate(eigvals(nppt,norb))
    allocate(grd_hk(nppt,norb,norb))
    if (allocated(work)) deallocate(work)
    if (allocated(rwork)) deallocate(rwork)
    lwork = 2*norb-1
    allocate(work(lwork),rwork(3*norb-2))
   !$OMP PARALLEL DEFAULT(SHARED) &
   !$OMP PRIVATE(ip,work, rwork)&
   !$OMP NUM_THREADS(12)
   !$OMP DO
    do ip=1,nppt
      call genhk(p_vecs(ip,:),grd_hk(ip,:,:)) 
      eigvecs_0(ip,:,:)=grd_hk(ip,:,:)
      call zheev('V','U',norb,eigvecs_0(ip,:,:),norb,eigvals(ip,:),work,lwork,rwork,lainfo)
    enddo
   !$OMP END DO
   !$OMP END PARALLEL
    vbmax=maxval(eigvals(:,norb/2))+minval(eigvals(:,norb/2+1)) ! shift to the VBmax
    write(*,*) "Vbmax(meV):", Vbmax/eV*1000
    eigvecs_t=eigvecs_0
!==========================================================

!BERRY CURVATURE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(BerryPlot) then
if(allocated(berryC)) deallocate(berryC)
allocate(berryC(norb,nppt))
berryC(:,:)=0.d0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ip,p_vec,temp_grpvel_mat)&
!$OMP NUM_THREADS(12)
!$OMP DO
  do ip=1,nppt
    if (mod(ip-1,np_mpi).ne.lp_mpi) cycle
    p_vec=p_vecs(ip,:)
    call gengrpvel(p_vec,temp_grpvel_mat)
    call genBerryC(eigvecs_0(ip,:,:),eigvals(ip,:),temp_grpvel_mat,berryC(:,ip))
  enddo
!$OMP END DO
!$OMP END PARALLEL
call mpi_barrier(mpicom,ierror)

if (np_mpi.gt.1) then
  nmpi=norb*nppt
  call mpi_allreduce(mpi_in_place,berryC,nmpi,mpi_double_precision,mpi_sum,mpicom,ierror)
endif

call mpi_barrier(mpicom,ierror)
if (mp_mpi) then
  open(50,file='berryC.txt',form='FORMATTED')
  do i=1,norb
    do ip=1,nppt
       write(50,'(5G18.10)') p_vecs(ip,1), p_vecs(ip,2), berryC(i,ip)
    enddo
  enddo
  close(50)
endif
call mpi_barrier(mpicom,ierror)
write(*,*) "Berry curvature calculation finished"
stop
endif
!==========================================================

!MESSAGES AND PLOTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if(mp_mpi) then
     write(*,27) "(Info) pump start(fs):", t_pump_start/fs
     27 format (a40,g18.10,g18.10)
     write(*,27) "(Info) probe start(fs):", (delta_t-tshift-hd_probe)/fs
     write(*,27) "(Info) pump half duration(fs):",hd_pump/fs
     write(*,27) "(Info) probe half duration(fs):",hd_probe/fs
     write(*,27) "(Info) delay(fs):",delta_t/fs
     write(*,27) "(Info) itmax:", itmax
     write(*,27) "(Info) exitmax:",exitmax
     write(*,27) "(Info) circular_pump:",circular_pump
     write(*,27) "(Info) circular_probe:",circular_probe
     write(*,27) "(Info) lineplotmode:", lineplotmode
     write(*,27) "(Info) nenrg:", nenrg
     if(lineplotmode) then
     write(*,27) "(Info) nkline :",nkline
     else
     write(*,27) "(Info) npx,npy:",npx,npy
     endif
     write(*,27) "(Info) nppt:", nppt
     write(*,27) "(Info) ngmax, ngtot:",ngmax, ngtot
     !-------------------------------------------------------
     open(50,file='pump.txt',form='FORMATTED')
       do itex=1,exitmax
       write(50,'(5G18.10)'),itex, itex*(dt/2.d0), pump(itex,1), pump(itex,2), pump(itex,3)
       enddo
     close(50)

     open(50,file='probe.txt',form='FORMATTED')
       do itex=1,exitmax
       write(50,'(5G18.10)'),itex, itex*(dt/2.d0), probe(itex,1), probe(itex,2), probe(itex,3)
       enddo
     close(50)

     open(50,file='eprobe.txt',form='FORMATTED')
       do itex=1,exitmax
       write(50,'(5G18.10)'),itex, itex*(dt/2.d0), e_probe(itex,1), e_probe(itex,2), e_probe(itex,3)
       enddo
     close(50)
     if(lineplotmode) then
     open(50,file='eigenvalues.txt',form='FORMATTED')
       do i=1,norb
       do ip=1,nppt
           write(50,'(5G18.10)'), kline(ip), eigvals(ip,i)
       enddo
       enddo
     close(50)
     endif
     !-------------------------------------------------------
   endif
!==========================================================

allocate(pes2(nppt,nenrg))
allocate(pes2_rcp(nppt,nenrg))
allocate(pes2_lcp(nppt,nenrg))
pes2(:,:)=0.d0
pes2_rcp(:,:)=0.d0
pes2_lcp(:,:)=0.d0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ip,it,i,ie,p_vec,hk_t_1,hk_t_2,hk_t_3)&
!$OMP PRIVATE(pes_tintegrand_vec,phasefactor,phase)&
!$OMP PRIVATE(j,pes_rintegral_vecs,phase2)&
!$OMP PRIVATE(ompthnum,nthd,percent)&
!$OMP PRIVATE(pes_vec,pes_temp,pes_rcp_temp,pes_lcp_temp)&
!$OMP NUM_THREADS(12)
ompthnum=omp_get_thread_num()
nthd=omp_get_num_threads()
!$OMP DO
do ip=1,nppt
  if (mod(ip-1,np_mpi).ne.lp_mpi) cycle
  it=1
  pes_vec(:,:,:)=0.d0
  phase(:,:)=0.d0
  percent=real(ip)*real(nthd)/real(nppt)*100d0
  if ((lp_mpi.eq.0).and.(ompthnum.eq.0)) write(*,'(a40,g18.5)') "(info) %",percent
  p_vec=p_vecs(ip,:)
  !call genhk(p_vec+pump(2*it-1,:)/sol,hk_t_1)
  do it=1,itmax-1
    !Define H(2*it-1), H(2*it) and H(2*it+1): When it=itmax-1, 2*it+1=2*(itmax-1)+1=2*itmax-1
    p_vec=p_vecs(ip,:) !here p_vec is an inplane vector
    !call genhk(p_vec+pump(2*it,:)/sol,hk_t_2)
    !call genhk(p_vec+pump(2*it+1,:)/sol,hk_t_3)
    do i=iorbstart,norb !!! CAUTION :: COINCIDE THE FIRST BAND INDEX
      do ie=1,nenrg
        if((2*enrg(ie)-p_vec(1)**2-p_vec(2)**2).ge.0.d0) then
          p_vec(3)=sqrt(2*enrg(ie)-p_vec(1)**2-p_vec(2)**2)
          phase(ie,i)=phase(ie,i)&
                      +(((p_vec(1)+pump(2*it-1,1)/sol)**2&
                      +(p_vec(2)+pump(2*it-1,2)/sol)**2&
                      +(p_vec(3)+pump(2*it-1,3)/sol)**2)/2.d0+vbmax-eigvals(ip,i))*dt
          phasefactor=cmplx(dcos(phase(ie,i)),dsin(phase(ie,i)),8)
          do j=1,norb 
          if((it.eq.1).and.(i.eq.iorbstart)) then !!! CAUTION :: COINCIDE THE FIRST BAND INDEX
          !position integral
          call intgrl_crtsn(p_vec,j,pes_rintegral_vecs(ie,j,:))
          endif

            phase2=-p_vec(3)*orb_pos(j,3)
            pes_tintegrand_vec=eigvecs_t(ip,j,i)*cmplx(dcos(phase2),dsin(phase2),8)*phasefactor&
                           *(/pes_rintegral_vecs(ie,j,1)*probe(2*it-1,1),&
                            pes_rintegral_vecs(ie,j,2)*probe(2*it-1,2),&
                            pes_rintegral_vecs(ie,j,3)*probe(2*it-1,3)/)*dt

            pes_vec(i,ie,:)=pes_vec(i,ie,:)+pes_tintegrand_vec
          enddo
        endif
      enddo
      !call rk4n(norb,hk_t_1,hk_t_2,hk_t_3,dt,eigvecs_t(ip,:,i))
    enddo
    !hk_t_1=hk_t_3
  enddo ! t sum
  p_vec=p_vecs(ip,:)
  do ie=1,nenrg
    if((2*enrg(ie)-p_vec(1)**2-p_vec(2)**2).ge.0.d0) then
      do i=iorbstart,norb !!! CUATION :: COINCIDE THE FIRST BAND INDEX
        if (circular_probe) then
        pes_rcp_temp=pes_vec(i,ie,1)-pes_vec(i,ie,2)
        pes_lcp_temp=pes_vec(i,ie,1)+pes_vec(i,ie,2)
        pes2_rcp(ip,ie)=pes2_rcp(ip,ie)+dble(pes_rcp_temp)**2+aimag(pes_rcp_temp)**2
        pes2_lcp(ip,ie)=pes2_lcp(ip,ie)+dble(pes_lcp_temp)**2+aimag(pes_lcp_temp)**2
        else
        pes_temp=pes_vec(i,ie,1)+pes_vec(i,ie,2)+pes_vec(i,ie,3)
        pes2(ip,ie)=pes2(ip,ie)+dble(pes_temp)**2+aimag(pes_temp)**2
        endif
      enddo
    endif
  enddo
enddo
!$OMP END DO
!$OMP END PARALLEL
call mpi_barrier(mpicom,ierror)

if (np_mpi.gt.1) then
  nmpi=nppt*nenrg
  if (circular_probe) then
  call mpi_allreduce(mpi_in_place,pes2_rcp,nmpi,mpi_double_precision,mpi_sum,mpicom,ierror)
  call mpi_allreduce(mpi_in_place,pes2_lcp,nmpi,mpi_double_precision,mpi_sum,mpicom,ierror)
  else
  call mpi_allreduce(mpi_in_place,pes2,nmpi,mpi_double_precision,mpi_sum,mpicom,ierror)
  endif
endif
call mpi_barrier(mpicom,ierror)

if (mp_mpi) then
if (lineplotmode) then
  if (circular_probe) then

    open(50,file='pes_line_rcp.txt',form='FORMATTED')
    do ip=1,nppt
    do ie=1,nenrg
      write(50,'(5G18.10)') kline(ip),enrg(ie),pes2_rcp(ip,ie)
    enddo
    enddo
    close(50)
    open(50,file='pes_line_lcp.txt',form='FORMATTED')
    do ip=1,nppt
    do ie=1,nenrg
      write(50,'(5G18.10)') kline(ip),enrg(ie),pes2_lcp(ip,ie)
    enddo
    enddo
    close(50)

  else

    open(50,file='pes_line.txt',form='FORMATTED')
    do ip=1,nppt
    do ie=1,nenrg
      write(50,'(5G18.10)') kline(ip),enrg(ie),pes2(ip,ie)
    enddo
    enddo
    close(50)

  endif

else
  if (circular_probe) then

    open(50,file='pes_rcp.txt',form='FORMATTED')
    do ip=1,nppt
    do ie=1,nenrg
      write(50,'(5G18.10)') p_vecs(ip,1),p_vecs(ip,2),enrg(ie),pes2_rcp(ip,ie)
    enddo
    enddo
    close(50)
    open(50,file='pes_lcp.txt',form='FORMATTED')
    do ip=1,nppt
    do ie=1,nenrg
      write(50,'(5G18.10)') p_vecs(ip,1),p_vecs(ip,2),enrg(ie),pes2_lcp(ip,ie)
    enddo
    enddo
    close(50)

  else

    open(50,file='pes.txt',form='FORMATTED')
    do ip=1,nppt
    do ie=1,nenrg
      write(50,'(5G18.10)') p_vecs(ip,1),p_vecs(ip,2),enrg(ie),pes2(ip,ie)
    enddo
    enddo
    close(50)

  endif

endif
endif

call mpi_barrier(mpicom,ierror)
if (mp_mpi) write(*,*) "calculation finished"

deallocate(hspts)
deallocate(n_p_to_p_grid)
deallocate(kline)
deallocate(work,rwork)
deallocate(pump,probe,e_probe)
deallocate(p_vecs)
deallocate(eigvecs_0)
deallocate(eigvecs_t)
deallocate(eigvals)
deallocate(grd_hk)
deallocate(pes2)
deallocate(pes2_rcp)
deallocate(pes2_lcp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!MPI AND OMP STOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! terminate MPI execution environment
call mpi_finalize(ierror)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!MPI AND OMP STOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
stop
contains
end program main
