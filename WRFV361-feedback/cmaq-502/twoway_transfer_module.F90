module twoway_transfer_module

  INTEGER :: conc_nvars, conc_nlays

  INTEGER :: conc_jdate = 0, conc_jtime = 0

  real*8 :: start_time(2), end_time(2), total_time(2) = 0.0d0

  REAL, ALLOCATABLE :: conc_data_wrf(:,:,:,:)
  REAL, ALLOCATABLE :: conc_data_cmaq(:,:,:,:)

  contains

! ------------------------------------------------------------------------------------
SUBROUTINE conc_wrf_write (grid)

  USE module_domain           ! WRF module

  USE twoway_data_module
  USE SE_MODULES
  
  IMPLICIT NONE
  
  include 'mpif.h'

  TYPE(domain), INTENT(IN) :: grid

  integer :: v, k

  start_time(1) = mpi_wtime()

  do v = 1, conc_nvars
     do k = 1, conc_nlays
        conc_data_wrf (:,:,k,v) = grid%tracer(sc:ec, k, sr:er, v+1)
     end do
  end do

  call se_wrf_cmaq_comm (twoway_mype, conc_data_wrf, conc_data_cmaq,           &
                         wrf_cmaq_c_send_to, wrf_cmaq_c_recv_from,             &
                         wrf_cmaq_c_send_index_l, wrf_cmaq_c_recv_index_l, 1)

  end_time(1) = mpi_wtime()
  total_time(1) = total_time(1) + (end_time(1) - start_time(1))

END SUBROUTINE conc_wrf_write

! ------------------------------------------------------------------------------------
SUBROUTINE conc_wrf_read (grid)

  USE module_domain                ! WRF module
  USE module_state_description     ! WRF module

  USE twoway_data_module
  USE SE_MODULES
  USE HGRD_DEFN

  IMPLICIT NONE

  include 'mpif.h'

  TYPE(domain), INTENT(IN) :: grid

  LOGICAL, SAVE :: firstime = .TRUE.

  integer :: stat, v, l, c, r, s, d, e

  logical, save :: north_bndy_pe = .false.
  logical, save :: east_bndy_pe  = .false.
  logical, save :: south_bndy_pe = .false.
  logical, save :: west_bndy_pe  = .false.
  character (len = 4), save :: pe_str

  start_time(1) = mpi_wtime()

  if (firstime) then

     if ((nprocs - mype) .le. npcol) then
        north_bndy_pe = .true.
     end if

     if (mod(mype, npcol) .eq. npcol - 1) then
        east_bndy_pe = .true.
     end if
  
     if (mype .lt. npcol) then
        south_bndy_pe = .true.
     end if
  
     if (mod(mype, npcol) .eq. 0) then
        west_bndy_pe = .true.
     end if

     firstime = .false.

  end if

  call se_cmaq_wrf_comm4 (twoway_mype, conc_data_cmaq,                            &
                         conc_data_wrf, cmaq_wrf_c_send_to, cmaq_wrf_c_recv_from, &
                         cmaq_wrf_c_send_index_l, cmaq_wrf_c_recv_index_l, 6)

  if (north_bndy_pe) then
     s = cmaq_c_domain_map(2,2,mype) - sr + 1
     do r = cmaq_c_domain_map(2,2,mype)+1, wrf_c_domain_map(2,2,mype)
        conc_data_wrf(:,r-sr+1,:,:) = conc_data_wrf(:,s,:,:)
     end do
  end if

  if (east_bndy_pe) then
     s = cmaq_c_domain_map(2,1,mype) - sc + 1
     d = wrf_c_domain_map(2,1,mype) - cmaq_c_domain_map(2,1,mype)
     do r = lbound(conc_data_wrf,2), ubound(conc_data_wrf,2)
        do c = s+1, s+d
           conc_data_wrf(c,r,:,:) = conc_data_wrf(s,r,:,:)
        end do
     end do
  end if

  if (south_bndy_pe) then
     do r = 1, delta_y
        conc_data_wrf(:,r,:,:) = conc_data_wrf(:,delta_y+1,:,:)
     end do
  end if

  if (west_bndy_pe) then
     do r = lbound(conc_data_wrf,2), ubound(conc_data_wrf,2)
       do c = 1, delta_x
           conc_data_wrf(c,r,:,:) = conc_data_wrf(delta_x+1,r,:,:)
        end do
     end do
  end if

! print *, ' ==d== pp ', p_ct001, p_ct002, p_ct079, size(grid%tracer, 2), conc_nlays

  do v = 1, conc_nvars
     do l = 1, conc_nlays
        do r = sr, er
           do c = sc, ec
!             grid%tracer(c, l, r, p_ct001) = conc_data_wrf(c-sc+1,r-sr+1,l,1)
              grid%tracer(c, l, r, v+1) = conc_data_wrf(c-sc+1,r-sr+1,l,v)
           end do
        end do
     end do
     grid%tracer(c, conc_nlays+1, r, v+1) = grid%tracer(c, conc_nlays, r, v)
  end do

  end_time(1) = mpi_wtime()
  total_time(1) = total_time(1) + (end_time(1) - start_time(1))

END SUBROUTINE conc_wrf_read
    
! ------------------------------------------------------------------------------------
SUBROUTINE conc_cmaq_write ( cgrid )
 
  USE twoway_data_module

  IMPLICIT NONE

  include 'mpif.h'

  INCLUDE SUBST_FILES_ID    ! file name parameters

  real, intent(in) :: cgrid(:,:,:,:)

  logical, save :: firstime = .true.
  integer       :: stat

  start_time(1) = mpi_wtime()

  IF ( firstime ) THEN

     allocate ( conc_data_cmaq (cmaq_c_ncols, cmaq_c_nrows, conc_nlays, conc_nvars), stat=stat)
     allocate ( conc_data_wrf (wrf_c_ncols, wrf_c_nrows, conc_nlays, conc_nvars), stat=stat)

     firstime = .false.

  ENDIF  ! first time

  conc_data_cmaq = cgrid

! print *, ' ==d== value conc_cmaq_write ', minval(cgrid(:,:,1,1)), maxval(cgrid(:,:,1,1))

  end_time(1) = mpi_wtime()
  total_time(1) = total_time(1) + (end_time(1) - start_time(1))

END SUBROUTINE conc_cmaq_write

! ------------------------------------------------------------------------------------
SUBROUTINE conc_cmaq_read ( cgrid )

  IMPLICIT NONE

  include 'mpif.h'

  real, intent(out) :: cgrid(:,:,:,:)

  start_time(1) = mpi_wtime()

  cgrid = conc_data_cmaq

! print *, ' ==d== value conc_cmaq_read  ', minval(cgrid(:,:,1,1)), maxval(cgrid(:,:,1,1))

  end_time(1) = mpi_wtime()
  total_time(1) = total_time(1) + (end_time(1) - start_time(1))

END SUBROUTINE conc_cmaq_read

end module twoway_transfer_module
