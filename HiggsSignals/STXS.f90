module STXS

! Still to do:
! 1: Read in correlation matrix
! 2: Write chi^2 test


!  use numerics
!  use combinatorics
 use usefulbits_hs
 implicit none

!  integer :: i,j,k
!  double precision,parameter :: pi=3.14159265358979323846264338328D0
!  integer, allocatable :: peakindices_best(:,:)

type STXS_observable
  integer :: id
  character(LEN=100) :: label  ! Reference
  character(LEN=100) :: desc   ! Description
  character(LEN=3)  :: expt    ! Experiment  
  character(LEN=10) :: collider
  character(LEN=10) :: collaboration
  double precision :: lumi,dlumi,energy  
  character(LEN=100) :: assignmentgroup  
  integer :: rate_SM_normalized
  integer :: mhchisq  
  double precision :: massobs, dmassobs ! This one enters the chi^2 for the mass!
  double precision :: mass, dmass ! This one is the mass position for the measurement and the "experimentally allowed assignment range"
  double precision ::  eff_ref_mass ! This is the mass for which the signal efficiencies are given.
  double precision, allocatable :: model_rate_per_Higgs(:,:)
  double precision, allocatable :: inclusive_SM_rate(:)
  integer :: Nc  
  double precision :: model_total_rate
  double precision :: rate, rate_up, rate_low, drate_up, drate_low
  double precision :: SMrate, SMrate_up, SMrate_low, dSMrate_up, dSMrate_low   ! SM rate used/quoted by the experiment
  ! At the moment, interpret STXS observables as "pure" channels (production, or decay rate)
  character(LEN=5), allocatable :: channel_id_str(:)		  ! Channels array as string, dim(Nc)  
!   integer, allocatable :: channel_id(:)
  integer, allocatable :: channel_p_id(:)		  ! Production channels array, dim(Nc)
  integer, allocatable :: channel_d_id(:)		  ! Decay channels array, dim(Nc)    
  double precision, allocatable :: channel_efficiency(:) ! SM signal efficiency of inclusive rates (analysis-specific)
  double precision, allocatable :: relative_efficiency(:,:)  ! Model signal efficiency relative to SM per Higgs 
  double precision :: chisq
!   character(LEN=10),allocatable :: channel_description(:,:)  
  ! TODO: How do we deal with ratio of BRs?  
end type 

type(STXS_observable), allocatable :: STXSlist(:)

type(correlation_info), allocatable :: STXScorrlist(:)
contains
!------------------------------------------------------------------------------------   
subroutine load_STXS(dataset)
!------------------------------------------------------------------------------------   
 use store_pathname_HS
 use usefulbits, only: file_id_common2, file_id_common3, np, Hneut
 use datatables, only : read_in_mass_resolution_and_assignment_group
!  implicit none

 character(LEN=*), intent(in) :: dataset
 character(LEN=100) :: datafile(500) 
  character(LEN=pathname_length+150) :: fullfilename     
 integer, allocatable :: skip(:)
 integer :: i, n, n_datafiles,n_correlations, n_correlations_tmp, ios, k, m, int1, int2
 double precision :: db1
 character(LEN=200) :: comment
 character(LEN=1) :: firstchar
 character(LEN=100) :: line
 integer :: id, posperiod

 call system('basename -a `ls -1 -p '//trim(adjustl(pathname_HS))// &
&	 'Expt_tables/'//trim(adjustl(dataset))//'/*.stxs 2>/dev/null` > STXS_analyses.txt 2>/dev/null')

 open(file_id_common3, file="STXS_analyses.txt",form='formatted')

 print *, "Reading in STXS measurements from analysis-set "//&
  			trim(adjustl(dataset))//":"

 n = 0
 n_datafiles = 0
 do
  n = n+1
  read(file_id_common3,'(A)', iostat=ios) datafile(n)
  if(ios.ne.0) exit
  write(*,'(I4,2X,A)') n, datafile(n)
 enddo
 n_datafiles = n - 1

 close(file_id_common3)
	
 allocate(STXSlist(n_datafiles),skip(n_datafiles))

  do n=1,n_datafiles
   skip(n)=11
     open(file_id_common3, file=trim(adjustl(pathname_HS)) //'Expt_tables/'// &
&	 trim(adjustl(dataset))//'/' // datafile(n))
   do 
   	read(file_id_common3,'(A)') comment
	comment = trim(adjustl(comment))
	write(firstchar,'(A1)') comment   
	if(firstchar.ne.'#') then
	 exit
    else
	 skip(n)=skip(n)+1
    endif  
   enddo 
   backspace(file_id_common3)

   read(file_id_common3,*) STXSlist(n)%id
   read(file_id_common3,'(A)') STXSlist(n)%label
   read(file_id_common3,*) STXSlist(n)%collider,STXSlist(n)%collaboration, &
   & STXSlist(n)%expt
   read(file_id_common3,'(A)') STXSlist(n)%desc
   read(file_id_common3,*) STXSlist(n)%energy, STXSlist(n)%lumi, STXSlist(n)%dlumi
   read(file_id_common3,*) STXSlist(n)%mhchisq, STXSlist(n)%rate_SM_normalized
   if(STXSlist(n)%mhchisq == 1) then
   read(file_id_common3,*) STXSlist(n)%massobs, STXSlist(n)%dmassobs
   else
   read(file_id_common3,*)
   STXSlist(n)%massobs = 0.0D0
   STXSlist(n)%dmassobs = 0.0D0
   endif
!--CHECK FOR ASSIGNMENT GROUP AS SECOND COLUMN:
   read(file_id_common3,*) STXSlist(n)%mass
   read(file_id_common3,'(A)') line
   call read_in_mass_resolution_and_assignment_group(line, STXSlist(n)%dmass,&
&   STXSlist(n)%assignmentgroup)
   read(file_id_common3,*) STXSlist(n)%Nc, STXSlist(n)%eff_ref_mass
   allocate(STXSlist(n)%channel_id_str(STXSlist(n)%Nc))
   allocate(STXSlist(n)%channel_p_id(STXSlist(n)%Nc))
   allocate(STXSlist(n)%channel_d_id(STXSlist(n)%Nc))   
   read(file_id_common3,*) (STXSlist(n)%channel_id_str(i),i=1,STXSlist(n)%Nc)
   do i=1,STXSlist(n)%Nc   
    posperiod = index(STXSlist(n)%channel_id_str(i),'.')
    if(posperiod.eq.0) then
     if(len(trim(adjustl(STXSlist(n)%channel_id_str(i)))).eq.2) then
      read(STXSlist(n)%channel_id_str(i),*) id
      STXSlist(n)%channel_p_id(i) = int((id-modulo(id,10))/dble(10))
      STXSlist(n)%channel_d_id(i) = modulo(id,10)
     else
      write(*,*) " For observable ID = ",STXSlist(n)%id
      stop " Error: Cannot handle channel IDs!"
     endif
    else
     read(STXSlist(n)%channel_id_str(i)(:posperiod-1),*) STXSlist(n)%channel_p_id(i)
     read(STXSlist(n)%channel_id_str(i)(posperiod+1:),*) STXSlist(n)%channel_d_id(i)
    endif
   enddo
!    write(*,*) "Production channels = ",STXSlist(n)%channel_p_id
!    write(*,*) "Decay channels = ",STXSlist(n)%channel_d_id      
!    allocate(STXSlist(n)%channel_id(STXSlist(n)%Nc))
!    read(file_id_common3,*) (STXSlist(n)%channel_id(i),i=1,STXSlist(n)%Nc)
   allocate(STXSlist(n)%channel_efficiency(STXSlist(n)%Nc))
   if(STXSlist(n)%eff_ref_mass.ge.0D0) then
    read(file_id_common3,*) (STXSlist(n)%channel_efficiency(i),i=1,STXSlist(n)%Nc)
   else
    do i=1,STXSlist(n)%Nc
     STXSlist(n)%channel_efficiency(i)=1.0D0
    enddo 
    read(file_id_common3,*)
   endif      
!    read(file_id_common3,*) STXSlist(n)%channel_id
!   read(file_id_common3,*) STXSlist(n)%relative_efficiency
   read(file_id_common3,*) STXSlist(n)%rate_low, STXSlist(n)%rate, STXSlist(n)%rate_up
   if(STXSlist(n)%rate_SM_normalized.eq.1) then
    read(file_id_common3,*) comment
    STXSlist(n)%SMrate_low = 0.0D0
    STXSlist(n)%SMrate = 0.0D0
    STXSlist(n)%SMrate_up = 0.0D0
   else
    read(file_id_common3,*) STXSlist(n)%SMrate_low, STXSlist(n)%SMrate, STXSlist(n)%SMrate_up
   endif 
   STXSlist(n)%drate_low = STXSlist(n)%rate - STXSlist(n)%rate_low
   STXSlist(n)%drate_up = STXSlist(n)%rate_up - STXSlist(n)%rate
   STXSlist(n)%dSMrate_low = STXSlist(n)%SMrate - STXSlist(n)%SMrate_low
   STXSlist(n)%dSMrate_up = STXSlist(n)%SMrate_up - STXSlist(n)%SMrate
   close(file_id_common3)
   
   allocate(STXSlist(n)%relative_efficiency(np(Hneut),STXSlist(n)%Nc))
   
   do k=1, np(Hneut)
    STXSlist(n)%relative_efficiency(k,:)=1.0D0
   enddo 
  enddo

  close(file_id_common3)
!NEW:
   call system('basename -a `ls -1 -p '//trim(adjustl(pathname_HS))// &
&	 'Expt_tables/'//trim(adjustl(dataset))//'/*.stxscorr 2>/dev/null` > STXS_correlations.txt 2>/dev/null')
   call system('rm -rf STXS_ncorrelations.txt')

   open(file_id_common3, file="STXS_correlations.txt",form='formatted')

  print *, "Reading in correlations from the following datafiles in analysis-set "//				&
  			trim(adjustl(Exptdir))//":"

  n = 0
  n_datafiles = 0
  n_correlations = 0
  do
   n = n+1
   read(file_id_common3,'(A)', iostat=ios) datafile(n)
   if(ios.ne.0) exit
   fullfilename=trim(adjustl(pathname_HS))//'Expt_tables/'//trim(adjustl(dataset))//'/'&
&                //trim(datafile(n))
   call system('cat '//trim(adjustl(fullfilename))//' | wc -l > STXS_ncorrelations.txt')
   open(file_id_common2,file="STXS_ncorrelations.txt",form='formatted')
   read(file_id_common2,'(I10)') n_correlations_tmp
   close(file_id_common2)
   write(*,'(2I4,2X,A)') n, n_correlations_tmp, datafile(n)
   n_correlations = n_correlations + n_correlations_tmp
  enddo
  n_datafiles = n - 1
  close(file_id_common3)  
  
  allocate(STXScorrlist(n_correlations))	

  m=0
  do n=1,n_datafiles
   fullfilename=trim(adjustl(pathname_HS))//'Expt_tables/'//trim(adjustl(dataset))//'/'&
&                //trim(datafile(n))  
   open(file_id_common3,file=fullfilename)   
   do
    m= m+1   
    read(file_id_common3,*,iostat=ios) int1, int2, db1
!     write(*,*) m, int1, int2, db1, ios
    if(ios.ne.0) exit
    STXScorrlist(m)%obsID1 = int1
    STXScorrlist(m)%obsID2 = int2
    STXScorrlist(m)%corr = db1
   enddo
   m=m-1
   close(file_id_common3)
  enddo  
  
end subroutine load_STXS
!------------------------------------------------------------------------------------   
subroutine assign_modelefficiencies_to_STXS(obsID, Nc, relative_efficiency)
!------------------------------------------------------------------------------------   
 use usefulbits, only : np, Hneut
 implicit none
 integer, intent(in) :: obsID
 integer, intent(in) :: Nc
 double precision, dimension(np(Hneut),Nc), intent(in) :: relative_efficiency
 
 integer :: i
 logical :: foundid = .False.

 do i=lbound(STXSlist, dim=1),ubound(STXSlist, dim=1)
  if(STXSlist(i)%id.eq.obsID) then
   if(Nc.ne.STXSlist(i)%Nc) then
    stop 'Error: Number of channels does not match!'
   else 
    STXSlist(i)%relative_efficiency = relative_efficiency
    foundid = .True.
   endif 
 endif
 enddo

 
 if(.not.foundid) write(*,*) "WARNING in assign_modelefficiencies_to_STXS: ",&
 & "observable ID ",obsID," not known!"

end subroutine assign_modelefficiencies_to_STXS
!------------------------------------------------------------------------------------   
subroutine get_chisq_from_STXS(chisq_tot, pval)
!------------------------------------------------------------------------------------   
use usefulbits, only : vsmall
use usefulbits_hs, only : Nparam
use numerics, only : invmatrix, matmult, gammp
 implicit none
 
 double precision, intent(out) :: chisq_tot, pval
 integer :: i,j,m,N
 double precision :: cov
 logical :: correlationfound, somecorrelationsmissing
 double precision, allocatable :: covmat(:,:),vmat(:,:),invcovmat(:,:)
 double precision, allocatable :: v(:), v2(:)
 N=size(STXSlist)

 allocate(covmat(N,N),invcovmat(N,N))
 allocate(v(N),v2(N))
 allocate(vmat(N,1))

 somecorrelationsmissing = .False.

 do i=1,N
  do j=1,N

   correlationfound=.False.
   do m=lbound(STXScorrlist,dim=1), ubound(STXScorrlist,dim=1)
    if((STXScorrlist(m)%obsID1.eq.STXSlist(i)%id.and.STXScorrlist(m)%obsID2.eq.STXSlist(j)%id)&
&.or.(STXScorrlist(m)%obsID2.eq.STXSlist(i)%id.and.STXScorrlist(m)%obsID1.eq.STXSlist(j)%id)) then
     covmat(i,j) = STXScorrlist(m)%corr*get_drate(i)*get_drate(j)
     correlationfound=.True.
    endif
   enddo
   if(.not.correlationfound) then
!  Use a unit-matrix for the correlations here.   
!     if(.not.somecorrelationsmissing) then
!      write(*,*) "Warning: Correlation matrix element not found for observable ids: ",STXSlist(i)%id, STXSlist(j)%id
!      write(*,*) "         Suppressing future warnings about missing correlation matrix elements..."    
!     endif 
    covmat(i,j) = 0.0D0
    somecorrelationsmissing = .True.
    if(STXSlist(i)%id.eq.STXSlist(j)%id) then    
     covmat(i,j) = get_drate(i)*get_drate(j)
    endif 
   endif 

  enddo
  v(i) = STXSlist(i)%rate - STXSlist(i)%model_total_rate
  vmat(i,1) = v(i)
 enddo
 
!  if(somecorrelationsmissing) then
!   write(*,*) "Warning: Some correlation matrix elements were not found."
!  endif
 
 call invmatrix(covmat,invcovmat)
 call matmult(invcovmat,vmat,v2,N,1)

 chisq_tot= 0.0D0

 do i=1,N
  STXSlist(i)%chisq = v(i)*v2(i)
  chisq_tot = chisq_tot + STXSlist(i)%chisq
 enddo

 pval = 1.0D0
  if(chisq_tot.gt.vsmall.and.(N-Nparam).gt.0) then
   pval = 1 - gammp(dble(N-Nparam)/2,chisq_tot/2)
  endif  
  deallocate(covmat,invcovmat,v,v2,vmat)



end subroutine get_chisq_from_STXS
!------------------------------------------------------------------------------------   
subroutine get_number_of_STXS_observables(Nobs_rates, Nobs_mh)
 integer, intent(out) ::  Nobs_rates, Nobs_mh
 
 Nobs_rates=size(STXSlist)
 Nobs_mh = 0
 
end subroutine get_number_of_STXS_observables
!------------------------------------------------------------------------------------   
function get_drate(i)
!------------------------------------------------------------------------------------   
implicit none

 integer :: i
 double precision get_drate

 if(STXSlist(i)%model_total_rate.le.STXSlist(i)%rate) then
  get_drate = STXSlist(i)%drate_low
 else
 get_drate = STXSlist(i)%drate_up  
 endif

end function get_drate
!------------------------------------------------------------------------------------   
subroutine calculate_model_predictions_for_STXS()
!------------------------------------------------------------------------------------   
 use usefulbits, only : theo
 use theo_manip, only : HB5_complete_theo
 integer :: i

 call HB5_complete_theo
  
 do i=lbound(STXSlist,dim=1), ubound(STXSlist,dim=1)
  call evaluate_model_for_STXS(STXSlist(i),theo(1))
 enddo

end subroutine calculate_model_predictions_for_STXS
!------------------------------------------------------------------------------------   
subroutine evaluate_model_for_STXS(STXSobs, t)
!------------------------------------------------------------------------------------   
 use usefulbits, only : theo, div, small, np, Hneut, dataset, vsmall
 use usefulbits_HS, only : normalize_rates_to_reference_position, &
&   normalize_rates_to_reference_position_outside_dmtheo, &
&   assignmentrange_STXS
! use_SMrate_at_reference_position_for_STXS,
 use theory_XS_SM_functions
 use theory_BRfunctions
 use theo_manip, only : HB5_complete_theo  
 implicit none
 
 
  type(STXS_observable), intent(inout) :: STXSobs
  type(dataset), intent(in) :: t  
  
  double precision :: norm_rate, SMrate, SMrate_refmass, refmass, BR_SMref
  integer :: i, j, id, p, d

  STXSobs%model_total_rate = 0.0D0
  
  if(.not.allocated(theo))then
   stop 'subroutine HiggsSignals_initialize must be called first'
  endif
  
  if(.not.allocated(STXSobs%model_rate_per_Higgs)) then
   allocate(STXSobs%model_rate_per_Higgs(np(Hneut),STXSobs%Nc))
  endif 
  
  if(.not.allocated(STXSobs%inclusive_SM_rate)) then
!    allocate(STXSobs%inclusive_SM_rate(np(Hneut),STXSobs%Nc))
   allocate(STXSobs%inclusive_SM_rate(STXSobs%Nc))
  endif 

! write(*,*) 'DEBUG HS: id = ', STXSobs%id
! write(*,*) 'DEBUG HS, channel = ',STXSobs%channel_id

  refmass = STXSobs%mass
  do i=1,STXSobs%Nc
!   id = STXSobs%channel_id(i)
!    p = int((id-modulo(id,10))/dble(10))
!    d = modulo(id,10)
   p = STXSobs%channel_p_id(i)
   d = STXSobs%channel_d_id(i)  
  
   do j=1, np(Hneut)
!  write(*,*) 'DEBUG HS, m = ', t%particle(Hneut)%M(j)  
!--Do the production rate for the relevant experiment and cms-energy 
  if(STXSobs%collider.eq.'LHC') then
   if(abs(STXSobs%energy-7.0D0).le.small) then
    if(p.eq.1) then 
     norm_rate=t%lhc7%XS_hj_ratio(j)
     SMrate=t%lhc7%XS_H_SM(j)
     SMrate_refmass=XS_lhc7_gg_H_SM(refmass)+XS_lhc7_bb_H_SM(refmass)
!      STXSobs%channel_description(i,1)='singleH'
    else if(p.eq.2) then
     norm_rate=t%lhc7%XS_vbf_ratio(j)
     SMrate=t%lhc7%XS_vbf_SM(j)
     SMrate_refmass=XS_lhc7_vbf_SM(refmass)
!      STXSobs%channel_description(i,1)='VBF'     
    else if(p.eq.3) then
     norm_rate=t%lhc7%XS_hjW_ratio(j)
     SMrate=t%lhc7%XS_HW_SM(j) 
     SMrate_refmass=XS_lhc7_HW_SM(refmass)
!      STXSobs%channel_description(i,1)='HW'     
    else if(p.eq.4) then
     norm_rate=t%lhc7%XS_hjZ_ratio(j)  
     SMrate=t%lhc7%XS_HZ_SM(j)
     SMrate_refmass=ZH_cpmix_nnlo_ggqqbb(refmass,'LHC7 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
!      STXSobs%channel_description(i,1)='HZ'       
    else if(p.eq.5) then
     norm_rate=t%lhc7%XS_tthj_ratio(j)
     SMrate=t%lhc7%XS_ttH_SM(j)
     SMrate_refmass=XS_lhc7_ttH_SM(refmass)
!      STXSobs%channel_description(i,1)='ttH'
    else if(p.eq.6) then
     norm_rate=t%lhc7%XS_gg_hj_ratio(j)
     SMrate=t%lhc7%XS_gg_H_SM(j)
     SMrate_refmass=XS_lhc7_gg_H_SM(refmass)
!      mutab%channel_description(i,1)='ggH'    
    else if(p.eq.7) then
     norm_rate=t%lhc7%XS_bb_hj_ratio(j)
     SMrate=t%lhc7%XS_bb_H_SM(j)
     SMrate_refmass=XS_lhc7_bb_H_SM(refmass)
!      mutab%channel_description(i,1)='bbH'
    else if(p.eq.8) then
     norm_rate=t%lhc7%XS_thj_tchan_ratio(j)
     SMrate=t%lhc7%XS_tH_tchan_SM(j)
     SMrate_refmass=XS_lhc7_tH_tchan_SM(refmass)
!      mutab%channel_description(i,1)='tH (t-channel)'     
    else if(p.eq.9) then
     norm_rate=t%lhc7%XS_thj_schan_ratio(j)
     SMrate=t%lhc7%XS_tH_schan_SM(j)
     SMrate_refmass=XS_lhc7_tH_schan_SM(refmass)
!      mutab%channel_description(i,1)='tH (s-channel)'
    else if(p.eq.10) then
     norm_rate=t%lhc7%XS_qq_hjZ_ratio(j)  
     SMrate=t%lhc7%XS_qq_HZ_SM(j)
     SMrate_refmass=ZH_cpmix_nnlo_qqbb(refmass,'LHC7 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
!      mutab%channel_description(i,1)='qq-HZ' 
    else if(p.eq.11) then
     norm_rate=t%lhc7%XS_gg_hjZ_ratio(j)  
     SMrate=t%lhc7%XS_gg_HZ_SM(j)
     SMrate_refmass=ZH_cpmix_nnlo_gg(refmass,'LHC7 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
!      mutab%channel_description(i,1)='gg-HZ'      
    else if(p.eq.0) then
     norm_rate=1.0D0
     SMrate=1.0D0
     SMrate_refmass=1.0D0
!      STXSobs%channel_description(i,1)='none'
    else
     write(*,*) "WARNING: Unknown production mode id p=",p," for STXS observable id = ",STXSobs%id 
    endif 
   else if(abs(STXSobs%energy-8.0D0).le.small) then
    if(p.eq.1) then 
     norm_rate=t%lhc8%XS_hj_ratio(j)
     SMrate=t%lhc8%XS_H_SM(j)
     SMrate_refmass=XS_lhc8_gg_H_SM(refmass)+XS_lhc8_bb_H_SM(refmass)
!      STXSobs%channel_description(i,1)='singleH'
    else if(p.eq.2) then
     norm_rate=t%lhc8%XS_vbf_ratio(j)
     SMrate=t%lhc8%XS_vbf_SM(j)
     SMrate_refmass=XS_lhc8_vbf_SM(refmass)     
!      STXSobs%channel_description(i,1)='VBF'     
    else if(p.eq.3) then
     norm_rate=t%lhc8%XS_hjW_ratio(j)
     SMrate=t%lhc8%XS_HW_SM(j) 
     SMrate_refmass=XS_lhc8_HW_SM(refmass)     
!      STXSobs%channel_description(i,1)='HW'     
    else if(p.eq.4) then
     norm_rate=t%lhc8%XS_hjZ_ratio(j)  
     SMrate=t%lhc8%XS_HZ_SM(j)
     SMrate_refmass=ZH_cpmix_nnlo_ggqqbb(refmass,'LHC8 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
!      STXSobs%channel_description(i,1)='HZ'       
    else if(p.eq.5) then
     norm_rate=t%lhc8%XS_tthj_ratio(j)
     SMrate=t%lhc8%XS_ttH_SM(j)
     SMrate_refmass=XS_lhc8_ttH_SM(refmass)     
!      STXSobs%channel_description(i,1)='ttH'
    else if(p.eq.6) then
     norm_rate=t%lhc8%XS_gg_hj_ratio(j)
     SMrate=t%lhc8%XS_gg_H_SM(j)
     SMrate_refmass=XS_lhc8_gg_H_SM(refmass)
!      mutab%channel_description(i,1)='ggH'    
    else if(p.eq.7) then
     norm_rate=t%lhc8%XS_bb_hj_ratio(j)
     SMrate=t%lhc8%XS_bb_H_SM(j)
     SMrate_refmass=XS_lhc8_bb_H_SM(refmass)
!      mutab%channel_description(i,1)='bbH'
    else if(p.eq.8) then
     norm_rate=t%lhc8%XS_thj_tchan_ratio(j)
     SMrate=t%lhc8%XS_tH_tchan_SM(j)
     SMrate_refmass=XS_lhc8_tH_tchan_SM(refmass)
!      mutab%channel_description(i,1)='tH (t-channel)'     
    else if(p.eq.9) then
     norm_rate=t%lhc8%XS_thj_schan_ratio(j)
     SMrate=t%lhc8%XS_tH_schan_SM(j)
     SMrate_refmass=XS_lhc8_tH_schan_SM(refmass)
!      mutab%channel_description(i,1)='tH (s-channel)'
    else if(p.eq.10) then
     norm_rate=t%lhc8%XS_qq_hjZ_ratio(j)  
     SMrate=t%lhc8%XS_qq_HZ_SM(j)
     SMrate_refmass=ZH_cpmix_nnlo_qqbb(refmass,'LHC8 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
!      mutab%channel_description(i,1)='qq-HZ' 
    else if(p.eq.11) then
     norm_rate=t%lhc8%XS_gg_hjZ_ratio(j)  
     SMrate=t%lhc8%XS_gg_HZ_SM(j)
     SMrate_refmass=ZH_cpmix_nnlo_gg(refmass,'LHC8 ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
    else if(p.eq.0) then
     norm_rate=1.0D0
     SMrate=1.0D0
     SMrate_refmass=1.0D0     
!      STXSobs%channel_description(i,1)='none'         
    else
     write(*,*) "WARNING: Unknown production mode id p=",p," for STXS observable id = ",STXSobs%id 
    endif  
   else if(abs(STXSobs%energy-13.0D0).le.small) then
    if(p.eq.1) then 
     norm_rate=t%lhc13%XS_hj_ratio(j)
     SMrate=t%lhc13%XS_H_SM(j)
     SMrate_refmass=XS_lhc13_gg_H_SM(refmass)+XS_lhc13_bb_H_SM(refmass)
!      STXSobs%channel_description(i,1)='singleH'
    else if(p.eq.2) then
     norm_rate=t%lhc13%XS_vbf_ratio(j)
     SMrate=t%lhc13%XS_vbf_SM(j)
     SMrate_refmass=XS_lhc13_vbf_SM(refmass)     
!      STXSobs%channel_description(i,1)='VBF'     
    else if(p.eq.3) then
     norm_rate=t%lhc13%XS_hjW_ratio(j)
     SMrate=t%lhc13%XS_HW_SM(j) 
     SMrate_refmass=XS_lhc13_HW_SM(refmass)     
!      STXSobs%channel_description(i,1)='HW'     
    else if(p.eq.4) then
     norm_rate=t%lhc13%XS_hjZ_ratio(j)  
     SMrate=t%lhc13%XS_HZ_SM(j)
     SMrate_refmass=ZH_cpmix_nnlo_ggqqbb(refmass,'LHC13',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
!      STXSobs%channel_description(i,1)='HZ'       
    else if(p.eq.5) then
     norm_rate=t%lhc13%XS_tthj_ratio(j)
     SMrate=t%lhc13%XS_ttH_SM(j)
     SMrate_refmass=XS_lhc13_ttH_SM(refmass)     
!      STXSobs%channel_description(i,1)='ttH' 
    else if(p.eq.6) then
     norm_rate=t%lhc13%XS_gg_hj_ratio(j)
     SMrate=t%lhc13%XS_gg_H_SM(j)
     SMrate_refmass=XS_lhc13_gg_H_SM(refmass)
!      mutab%channel_description(i,1)='ggH'    
    else if(p.eq.7) then
     norm_rate=t%lhc13%XS_bb_hj_ratio(j)
     SMrate=t%lhc13%XS_bb_H_SM(j)
     SMrate_refmass=XS_lhc13_bb_H_SM(refmass)
!      mutab%channel_description(i,1)='bbH'
    else if(p.eq.8) then
     norm_rate=t%lhc13%XS_thj_tchan_ratio(j)
     SMrate=t%lhc13%XS_tH_tchan_SM(j)
     SMrate_refmass=XS_lhc13_tH_tchan_SM(refmass)
!      mutab%channel_description(i,1)='tH (t-channel)'     
    else if(p.eq.9) then
     norm_rate=t%lhc13%XS_thj_schan_ratio(j)
     SMrate=t%lhc13%XS_tH_schan_SM(j)
     SMrate_refmass=XS_lhc13_tH_schan_SM(refmass)
!      mutab%channel_description(i,1)='tH (s-channel)'
    else if(p.eq.10) then
     norm_rate=t%lhc13%XS_qq_hjZ_ratio(j)  
     SMrate=t%lhc13%XS_qq_HZ_SM(j)
     SMrate_refmass=ZH_cpmix_nnlo_qqbb(refmass,'LHC13',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
!      mutab%channel_description(i,1)='qq-HZ' 
    else if(p.eq.11) then
     norm_rate=t%lhc13%XS_gg_hjZ_ratio(j)  
     SMrate=t%lhc13%XS_gg_HZ_SM(j)
     SMrate_refmass=ZH_cpmix_nnlo_gg(refmass,'LHC13',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
    else if(p.eq.0) then
     norm_rate=1.0D0
     SMrate=1.0D0
     SMrate_refmass=1.0D0     
!      STXSobs%channel_description(i,1)='none'         
    else
     write(*,*) "WARNING: Unknown production mode id p=",p," for STXS observable id = ",STXSobs%id 
    endif 
   endif 
  else if(STXSobs%collider.eq.'TEV') then
    if(p.eq.1) then 
     norm_rate=t%tev%XS_hj_ratio(j)
     SMrate=t%tev%XS_H_SM(j)
     SMrate_refmass=XS_tev_gg_H_SM(refmass)+XS_tev_bb_H_SM(refmass)
!      STXSobs%channel_description(i,1)='singleH'
    else if(p.eq.2) then
     norm_rate=t%tev%XS_vbf_ratio(j)
     SMrate=t%tev%XS_vbf_SM(j)
     SMrate_refmass=XS_tev_vbf_SM(refmass)     
!      STXSobs%channel_description(i,1)='VBF'     
    else if(p.eq.3) then
     norm_rate=t%tev%XS_hjW_ratio(j)
     SMrate=t%tev%XS_HW_SM(j) 
     SMrate_refmass=XS_tev_HW_SM(refmass)
!      STXSobs%channel_description(i,1)='HW'     
    else if(p.eq.4) then
     norm_rate=t%tev%XS_hjZ_ratio(j)  
     SMrate=t%tev%XS_HZ_SM(j)
     SMrate_refmass=ZH_cpmix_nnlo_ggqqbb(refmass,'TEV  ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
!      STXSobs%channel_description(i,1)='HZ'       
    else if(p.eq.5) then
     norm_rate=t%tev%XS_tthj_ratio(j)
     SMrate=t%tev%XS_ttH_SM(j)
     SMrate_refmass=XS_tev_ttH_SM(refmass)
!      STXSobs%channel_description(i,1)='ttH' 
    else if(p.eq.10) then
     norm_rate=t%tev%XS_qq_hjZ_ratio(j)  
     SMrate=t%tev%XS_qq_HZ_SM(j)
     SMrate_refmass=ZH_cpmix_nnlo_qqbb(refmass,'TEV  ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
!      mutab%channel_description(i,1)='qq-HZ' 
    else if(p.eq.11) then
     norm_rate=t%tev%XS_gg_hjZ_ratio(j)  
     SMrate=t%tev%XS_gg_HZ_SM(j)
     SMrate_refmass=ZH_cpmix_nnlo_gg(refmass,'TEV  ',1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,.True.)
    else if(p.eq.0) then
     norm_rate=1.0D0
     SMrate=1.0D0
     SMrate_refmass=1.0D0     
!      STXSobs%channel_description(i,1)='none'         
    else
     write(*,*) "WARNING: Unknown production mode id p=",p," for STXS observable id = ",STXSobs%id 
    endif        
  else if(STXSobs%collider.eq.'ILC') then
!--n.B.: As a first attempt, we use the LHC8 normalized cross sections for ZH, VBF, ttH.
!        In order to do this properly, a separate input for the ILC cross sections
!        has to be provided! It works only for single production mode observables (no
!        correct weighting of channels included!)Then, at least in the effective coupling
!        approximation, there is no difference to a full implementation.
!        The theoretical uncertainty of the ILC production modes will are defined in
!        usefulbits_HS.f90.
    if(p.eq.1.or.p.eq.2) then 
     write(*,*) 'Warning: Unknown ILC production mode (',p,') in table ',STXSobs%id
     norm_rate=0.0D0
     SMrate=1.0D0
     SMrate_refmass=1.0D0     
!      STXSobs%channel_description(i,1)='unknown'
    else if(p.eq.3) then
     norm_rate=t%lhc8%XS_hjW_ratio(j)
     SMrate=t%lhc8%XS_HW_SM(j)
     SMrate_refmass=XS_lhc8_HW_SM(refmass)          
!      STXSobs%channel_description(i,1)='WBF'     
    else if(p.eq.4) then
     norm_rate=t%lhc8%XS_hjZ_ratio(j)  
     SMrate=t%lhc8%XS_HZ_SM(j)
     SMrate_refmass=XS_lhc8_HZ_SM(refmass)     
!      STXSobs%channel_description(i,1)='HZ'       
    else if(p.eq.5) then
     norm_rate=t%lhc8%XS_tthj_ratio(j)
     SMrate=t%lhc8%XS_ttH_SM(j)
     SMrate_refmass=XS_lhc8_ttH_SM(refmass)
!      STXSobs%channel_description(i,1)='ttH' 
    else if(p.eq.0) then
     norm_rate=1.0D0
     SMrate=1.0D0
     SMrate_refmass=1.0D0     
!      STXSobs%channel_description(i,1)='none'
    else
     write(*,*) "WARNING: Unknown production mode id p=",p," for STXS observable id = ",STXSobs%id 
    endif           
   endif
!--Multiply now by the decay rate
  if(d.eq.1) then
   norm_rate=norm_rate*div(t%BR_hjgaga(j),t%BR_Hgaga_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hgaga_SM(j)
   SMrate_refmass = SMrate_refmass*BRSM_Hgaga(refmass)
!    STXSobs%channel_description(i,2)='gammagamma'   
  else if(d.eq.2) then
   norm_rate=norm_rate*div(t%BR_hjWW(j),t%BR_HWW_SM(j),0.0D0,1.0D0)   
   SMrate=SMrate*t%BR_HWW_SM(j)
   SMrate_refmass = SMrate_refmass*BRSM_HWW(refmass)
!    STXSobs%channel_description(i,2)='WW'   
  else if(d.eq.3) then
   norm_rate=norm_rate*div(t%BR_hjZZ(j),t%BR_HZZ_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_HZZ_SM(j)
   SMrate_refmass = SMrate_refmass*BRSM_HZZ(refmass)   
!    STXSobs%channel_description(i,2)='ZZ'
  else if(d.eq.4) then
   norm_rate=norm_rate*div(t%BR_hjtautau(j),t%BR_Htautau_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Htautau_SM(j)
   SMrate_refmass = SMrate_refmass*BRSM_Htautau(refmass)
!    STXSobs%channel_description(i,2)='tautau'
  else if(d.eq.5) then
   norm_rate=norm_rate*div(t%BR_hjbb(j),t%BR_Hbb_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hbb_SM(j)
   SMrate_refmass = SMrate_refmass*BRSM_Hbb(refmass)
!    STXSobs%channel_description(i,2)='bb'  
  else if(d.eq.6) then
   norm_rate=norm_rate*div(t%BR_hjZga(j),t%BR_HZga_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_HZga_SM(j)
   SMrate_refmass = SMrate_refmass*BRSM_HZga(refmass)   
!    STXSobs%channel_description(i,2)='Zgamma'
  else if(d.eq.7) then
   norm_rate=norm_rate*div(t%BR_hjcc(j),t%BR_Hcc_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hcc_SM(j)
   SMrate_refmass = SMrate_refmass*BRSM_Hcc(refmass)
!    STXSobs%channel_description(i,2)='cc'
  else if(d.eq.8) then
   norm_rate=norm_rate*div(t%BR_hjmumu(j),t%BR_Hmumu_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hmumu_SM(j)
   SMrate_refmass = SMrate_refmass*BRSM_Hmumu(refmass)
!    STXSobs%channel_description(i,2)='mumu'
  else if(d.eq.9) then
   norm_rate=norm_rate*div(t%BR_hjgg(j),t%BR_Hgg_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hgg_SM(j)
   SMrate_refmass = SMrate_refmass*BRSM_Hgg(refmass)
!    STXSobs%channel_description(i,2)='gg'
  else if(d.eq.10) then
   norm_rate=norm_rate*div(t%BR_hjss(j),t%BR_Hss_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Hss_SM(j)
   SMrate_refmass = SMrate_refmass*BRSM_Hss(refmass)
!    mutab%channel_description(i,2)='ss'
  else if(d.eq.11) then
   norm_rate=norm_rate*div(t%BR_hjtt(j),t%BR_Htt_SM(j),0.0D0,1.0D0)
   SMrate=SMrate*t%BR_Htt_SM(j)
   SMrate_refmass = SMrate_refmass*BRSM_Htoptop(refmass)
!    mutab%channel_description(i,2)='tt'
  else if(d.eq.0) then
   norm_rate=norm_rate*1.0D0
   SMrate=SMrate*1.0D0
   SMrate_refmass = SMrate_refmass*1.0D0
!    STXSobs%channel_description(i,2)='none'
  endif
  
!-------------------------
! NEW FEATURE (since HB-5.2): Enable to set channelrates directly.
  if(p.ne.0.and.d.ne.0) then
   select case(d)
    case(1)
     BR_SMref = t%BR_Hgaga_SM(j)
!      BR_SMref_mpeak = BRSM_Hgaga(refmass)
    case(2)
     BR_SMref = t%BR_HWW_SM(j)
!      BR_SMref_mpeak = BRSM_HWW(refmass)
    case(3)
     BR_SMref = t%BR_HZZ_SM(j)
!      BR_SMref_mpeak = BRSM_HZZ(refmass)
    case(4)
     BR_SMref = t%BR_Htautau_SM(j)
!      BR_SMref_mpeak = BRSM_Htautau(refmass)
    case(5)
     BR_SMref = t%BR_Hbb_SM(j)
!      BR_SMref_mpeak = BRSM_Hbb(refmass)
    case(6)
     BR_SMref = t%BR_HZga_SM(j)
!      BR_SMref_mpeak = BRSM_HZga(refmass)
    case(7)
     BR_SMref = t%BR_Hcc_SM(j)
!      BR_SMref_mpeak = BRSM_Hcc(refmass)
    case(8)
     BR_SMref = t%BR_Hmumu_SM(j)
!      BR_SMref_mpeak = BRSM_Hmumu(refmass)
    case(9)
     BR_SMref = t%BR_Hgg_SM(j)
!      BR_SMref_mpeak = BRSM_Hgg(refmass)
    case(10)
     BR_SMref = t%BR_Hss_SM(j)
!      BR_SMref_mpeak = BRSM_Hgg(refmass)
    case(11)
     BR_SMref = t%BR_Htt_SM(j)
!      BR_SMref_mpeak = BRSM_Htoptop(refmass)

   end select
   if(STXSobs%collider.eq.'LHC') then
    if(abs(STXSobs%energy-7.0D0).le.small) then
     if(t%lhc7%channelrates(j,p,d).ge.0.0d0) then
      norm_rate=div(t%lhc7%channelrates(j,p,d),BR_SMref,0.0D0,1.0D0)
     endif
    else if(abs(STXSobs%energy-8.0D0).le.small) then  
     if(t%lhc8%channelrates(j,p,d).ge.0.0d0) then
      norm_rate=div(t%lhc8%channelrates(j,p,d),BR_SMref,0.0D0,1.0D0)
     endif
    else if(abs(STXSobs%energy-13.0D0).le.small) then  
     if(t%lhc13%channelrates(j,p,d).ge.0.0d0) then
      norm_rate=div(t%lhc13%channelrates(j,p,d),BR_SMref,0.0D0,1.0D0)
     endif
    endif 
   else if(STXSobs%collider.eq.'TEV') then
     if(t%tev%channelrates(j,p,d).ge.0.0d0) then
      norm_rate=div(t%tev%channelrates(j,p,d),BR_SMref,0.0D0,1.0D0)
     endif
   endif
  endif  
!------------------------- 
 if(abs(t%particle(Hneut)%M(j) - STXSobs%mass).le.(assignmentrange_STXS * &
    sqrt(t%particle(Hneut)%dM(j)**2.0D0+STXSobs%dmass**2.0D0))) then  
!  if(STXSobs%rate_SM_normalized.eq.1) then    
  if(normalize_rates_to_reference_position) then 
   STXSobs%model_rate_per_Higgs(j,i)=norm_rate*SMrate/(SMrate_refmass)

  else
   STXSobs%model_rate_per_Higgs(j,i)=norm_rate  !! OLD WAY
  endif
  
 if(normalize_rates_to_reference_position_outside_dmtheo) then
  if(abs(STXSobs%mass-t%particle(Hneut)%M(j)).ge.t%particle(Hneut)%dM(j)) then
   STXSobs%model_rate_per_Higgs(j,i)=norm_rate*SMrate/(SMrate_refmass)
   endif
  endif  
 
!  else
!   if(use_SMrate_at_reference_position_for_STXS) then
!---
! n.B.: Need to use officially quoted SM prediction here, because HB/HS do not contain
!       SM predictions for exclusive STXS bins (but only inclusive SM rates).
!---
!    STXSobs%model_rate_per_Higgs(j,i)=norm_rate*STXSobs%SMrate   
!   else
!    STXSobs%model_rate=norm_rate*SMrate   
!   endif
!  endif
 else
  STXSobs%model_rate_per_Higgs(j,i) = 0.0D0
!   STXSobs%inclusive_SM_rate(j,i) = 0.0D0
 endif

! Inclusive SM rate must always be evaluated at the mass position of the measurement!
 STXSobs%inclusive_SM_rate(i) = SMrate_refmass * STXSobs%channel_efficiency(i)      


!  write(*,*) "j, i, STXSobs%model_rate_per_Higgs(j,i) = ",j,i, STXSobs%model_rate_per_Higgs(j,i)

! Turn normalized rate into absolute rate (per Higgs per channel) 
 STXSobs%model_rate_per_Higgs(j,i) = STXSobs%model_rate_per_Higgs(j,i) * &
 &                                   STXSobs%relative_efficiency(j,i) *  &
 &                                   STXSobs%inclusive_SM_rate(i)
 !&                                   SMrate * STXSobs%channel_efficiency(i)


!  write(*,*) "j, i, absolute  STXSobs%model_rate_per_Higgs(j,i),  STXSobs%inclusive_SM_rate(i) = ",&
!  & j, i,  STXSobs%model_rate_per_Higgs(j,i), STXSobs%inclusive_SM_rate(j,i)

!---
! Take into account model-dependent signal efficiency (relative to SM).
! These have to be given by the user for each observable using the subroutine
! assign_modelefficiencies_to_STXS:
!--- 
 enddo 
 enddo

! write(*,*) "STXSobs%id = " , STXSobs%id
! write(*,*) " model_rate_per_Higgs = ",STXSobs%model_rate_per_Higgs
! write(*,*) " inclusive SM rate = ",STXSobs%inclusive_SM_rate

if(sum(STXSobs%inclusive_SM_rate).ge.vsmall) then
 STXSobs%model_total_rate = sum(STXSobs%model_rate_per_Higgs)/sum(STXSobs%inclusive_SM_rate) ! Mistake: Don't divide by the sum over SM!!
else 
 STXSobs%model_total_rate = 0.0D0
endif
  
!  write(*,*) "STXSobs%model_total_rate (SM norm)= ",  STXSobs%model_total_rate
 
 if(STXSobs%rate_SM_normalized.eq.0) then
  STXSobs%model_total_rate = STXSobs%model_total_rate * STXSobs%SMrate
!   write(*,*) "STXSobs%model_total_rate (absolute)= ", STXSobs%model_total_rate 
 endif

!  write(*,*) "#--------------- ", STXSobs%id, " ---------------#"
! do i=1,STXSobs%Nc
!  write(*,*) "channel id = ", STXSobs%channel_id(i), " rate = ", &
!  & sum(STXSobs%model_rate_per_Higgs(:,i))/sum(STXSobs%inclusive_SM_rate)*STXSobs%SMrate
! enddo 
!   STXSobs%model_total_rate + STXSobs%relative_efficiency(j) * STXSobs%model_rate_per_Higgs(j) 

 
! write(*,*) "Total rate: ", STXSobs%model_total_rate 

end subroutine evaluate_model_for_STXS
!------------------------------------------------------------------------------------
subroutine print_STXS()
!------------------------------------------------------------------------------------   
 implicit none
  integer :: i
  character(LEN=100) :: formatter
  
  do i=lbound(STXSlist,dim=1), ubound(STXSlist,dim=1)
  write(*,*) "#--------------------------------------------------#"  
  write(*,*) "#-              STXS observable ",i,"             -#"
  write(*,*) "#--------------------------------------------------#"    
  write(*,'(A,I10)') " ID                     = ", STXSlist(i)%id
  write(*,'(A,A)') " Label                  = ", STXSlist(i)%label
  write(*,'(A,A)') " Description            = ", STXSlist(i)%desc  
  write(*,'(A,A)') " Experiment             = ", STXSlist(i)%expt
  write(*,'(A,2F6.2)') " Energy, Luminosity     = ", STXSlist(i)%energy, STXSlist(i)%lumi
  write(*,'(A,F10.5,A,F10.5,A,F10.5)') " Obs Signal rate [pb]   = ",&
  & STXSlist(i)%rate, " + ", STXSlist(i)%drate_up, " - ", STXSlist(i)%drate_low
  write(*,'(A,F10.5,A,F10.5,A,F10.5)') " SM Signal rate [pb]    = ",&
  & STXSlist(i)%SMrate, " + ", STXSlist(i)%dSMrate_up, " - ", STXSlist(i)%dSMrate_low  
  write(*,'(A,F10.5)') " Pred. Signal rate [pb] = ", STXSlist(i)%model_total_rate  
  write(formatter,*) "(A,",STXSlist(i)%Nc,"I10)"
  formatter = trim(adjustl(formatter))
  write(*,'(A)') " Channels               = ", STXSlist(i)%channel_id_str
  write(formatter,*) "(A,",STXSlist(i)%Nc,"F10.5)"
  formatter = trim(adjustl(formatter))
  write(*,formatter) " Channel efficiency     = ", STXSlist(i)%channel_efficiency
 enddo
write(*,*) "#--------------------------------------------------#"

end subroutine print_STXS
!------------------------------------------------------------------------------------   
subroutine print_STXS_to_file
!------------------------------------------------------------------------------------    
 use usefulbits, only : file_id_common3
 use usefulbits_hs, only : StrCompress
 implicit none
 character(LEN=100) :: formatspec
 integer :: i

 formatspec='(I3,7X,I10,1X,F6.2,1X,6F10.6,1X,A3,1X,F6.2,1X,F6.2,1X,A,5X,A)'
 open(file_id_common3,file="STXS_information.txt")
 write(file_id_common3,*) "#HiggsSignals-"//trim(adjustl(HSvers))//						&
&						  " with experimental dataset '"//trim(adjustl(Exptdir))//"'" 
 write(file_id_common3,*) "#Number STXS-ID mass-pos   rate_obs  drate_low  drate_high ",		&
 &				"rate_SM  dSMrate_low  dSMrate_high  collaboration   energy 	luminosity   description   reference"
 write(file_id_common3,*) "#"
 do i=lbound(STXSlist,dim=1),ubound(STXSlist,dim=1)
  write(file_id_common3,formatspec) i ,STXSlist(i)%id,STXSlist(i)%mass,		&
  & STXSlist(i)%rate, STXSlist(i)%drate_low,STXSlist(i)%drate_up,	&
  & STXSlist(i)%SMrate, STXSlist(i)%dSMrate_low,STXSlist(i)%dSMrate_up,	&  
  & STXSlist(i)%collaboration, STXSlist(i)%energy,	&
  & STXSlist(i)%lumi, trim(strcompress(STXSlist(i)%desc)), STXSlist(i)%label
 enddo 
 close(file_id_common3)
end subroutine print_STXS_to_file
!------------------------------------------------------------------------------------  
subroutine clear_STXS()
!------------------------------------------------------------------------------------   
 implicit none
 integer :: i
 
 do i=lbound(STXSlist,dim=1), ubound(STXSlist,dim=1)
  deallocate(STXSlist(i)%model_rate_per_Higgs)
  deallocate(STXSlist(i)%inclusive_SM_rate)
  deallocate(STXSlist(i)%channel_id_str)
  deallocate(STXSlist(i)%channel_p_id)
  deallocate(STXSlist(i)%channel_d_id)
  deallocate(STXSlist(i)%channel_efficiency)
  deallocate(STXSlist(i)%relative_efficiency)
 enddo
 deallocate(STXSlist)
 
 if(allocated(STXScorrlist)) deallocate(STXScorrlist)
 
end subroutine clear_STXS
!------------------------------------------------------------------------------------   
end module STXS
!------------------------------------------------------------------------------------  