MODULE dmrg
    use utilities
    use tensors
    implicit none

    public
!>@brief \b type definition for controlling dmrg program execution
!>@todo most of the variables are dummy ones, some attention is needed
    type dmrg_controls
        !type(time_profiles), allocatable, dimension(:) :: profile       !< store up to 10 time profiles, identified by a user defined strings of at most 512 characters
        type(time_profiles), dimension(12) :: profile       !< store up to 10 time profiles, identified by a user defined strings of at most 512 characters
        integer             ::  time_profiles = -1                      !< how many fields should be allocated for time profiling, if 0 or lower - no time profiles
        real*8              ::  tol = 1D-10          !< tolerance for arpack, lapack and dmrg, @todo remember to pass TOL to arpack and lapack procedures
        integer             ::  max_states = 20     !< maximal number of density matrix eigenvectors stored
        integer             ::  omp_num_threads = 8 !< default number of threads
        integer             ::  lapack_threads = 1  !< number of threads for lapack - poor scaling, only 1 thread
        integer             ::  arpack_threads = 8  !< number of threasd for arpack, little difference between 4 and 8 threads
        integer             ::  error_code          !< pass some error codes from procedures - I have no idea what was it supposed to be
        integer             ::  current_step            !
        integer             ::  spin_selector=0     !< choose states with spin = spin_selector* 1/2
        character(len=512)  ::  error_string        !< pass information about encountered error
        real*8              ::  ev                  !< ground state energy found in current step @todo don't forget to actually save the eigenvalue
        real*8              ::  truncation_err      !< the truncation error of current step
    end type dmrg_controls
    
!>@brief \b type definition for observables, use when observable for current length is only needed explicite
    type observable_8
        type(tensor_2d_8) ::  site                                !< double precision tensor_2d_8 object
        character*512     ::  name                                      !< the name of observable
        logical           ::  update = .false.                          !< should the observalbe be updated? Leave observable update for the last sweep
        logical           ::  calculate = .false.                       !< should the observalbe be updated? Leave observable update for the last sweep
        real*8            ::  expected_val                              !< expected value of an observable in a given state
    end type observable_8
!>@brief \b type definition for observables, members: allocatable vector of tensor_2d_4 type and character*512 string containing name of observable
    type observable_4
        type(tensor_2d_4), allocatable, dimension(:)   ::  site         !< array of single precision tensor_2d_8 objects, preferably containing observable for each site written in current basis
        character*512                                  ::  name         !< the name of observable
    end type observable_4
!>@brief \b type definition with four double precision array-like tensor_2d_8 members: Sz, Sp, Sm and H
!>@todo change the name to spin_ops_8 and add single precision version
    type spin_ops
        type(tensor_2d_8) ::  Sz                                        !< double precision tensor_2d_8 object, preferably containing single site spin \f$S_z\f$ operator
        type(tensor_2d_8) ::  Sp                                        !< double precision tensor_2d_8 object, preferably containing single site spin \f$S_p\f$ operator
        type(tensor_2d_8) ::  Sm                                        !< double precision tensor_2d_8 object, preferably containing single site spin \f$S_m\f$ operator
        type(tensor_2d_8) ::  H                                         !< double precision tensor_2d_8 object, preferably containing single site Hamiltonian \f$H\f$
        integer*4,allocatable, dimension(:) ::  basis_Sz                !< integer vector for storing expected value of Sz in basis states, in havles, that is Sz|up> = basis_Sz/2 |up>
        real*8, allocatable, dimension(:)   ::  guess                   !< guess of the ground state in reduced basis (for single Sz selector)
    end type spin_ops

!>@brief \b subroutine interface for printing spin_ops variable; usage: <em> call print(spin_ops variable) </em>
    interface print
        module procedure printSpinOps
    end interface    
!>@brief \b subroutine interface for deallocation of matrices within spin_ops objects; usage: <em> call deallocate(spin_ops variable)</em>
    interface deallocate
        module procedure deallocate_spin_ops
    end interface
!>@brief \b operator .is. for checking whether two spin_ops variables have the same address; 
    interface operator( .is. )
        module procedure isTheSame_spin_ops
    end interface
!>@brief \b assignment = overloading for fast tensor_2d_4 and tensor_2d_8 initialization; usage: <em> variable='Sz'; variable='Sp'; variable='Sm'; variable='H'</em>
    interface assignment(=)
        module procedure initTensor2d_8, initTensor2d_4
    end interface
!>@brief \b subroutine interface for performing single DMRG step
    interface singleStepDMRG
        module procedure singleStepDMRG_8
    end interface
!>@brief \b subroutine interface for performing all observable updates (expandig by 2 and rotating and truncating by matrix O)    
!>@todo change expanding by 2 to expanding by an integer defined somewhere, perhaps in another parameter (single site tensor object?)
    interface updateAllObservables
        module procedure updateAllObservables_8
    end interface
!>@brief \b subroutine interface for performing single observable update (expandig by 2 and rotating and truncating by matrix O)    
!>@todo change expanding by 2 to expanding by an integer defined somewhere, perhaps in another parameter (single site tensor object?)
    interface updateSingleObservable
        module procedure updateSingleObservable_8
    end interface 
    contains
!>@brief Check whether addresses of two objects are the same.
    function isTheSame_spin_ops(ops1,ops2) result(same)
        implicit none
        type(spin_ops), intent(in)      ::  ops1                        !< @param - object for checking
        type(spin_ops), intent(in)      ::  ops2                        !< @param - object for checking
        logical                         ::  same                        !< @return If objects have the same address in memory: .true., otherwise: .false.
        same=.false.
        if (loc(ops1) == loc(ops2)) same=.true.
    end function isTheSame_spin_ops

!>@brief Deallocate allocatable parts of spin_ops fields: Sp, Sm, Sz and H. Can be run on objects with unallocated internal fields.
    subroutine deallocate_spin_ops(ops)
        implicit none
        type(spin_ops), intent(inout)   ::  ops                         !< @param - object of type spin_ops, containing Sp, Sm, Sz and H fields of type tensor_2d_8
        call deallocate_tensor_2d_8(ops%Sp)
        call deallocate_tensor_2d_8(ops%Sm)
        call deallocate_tensor_2d_8(ops%Sz)
        call deallocate_tensor_2d_8(ops%H)
    end subroutine deallocate_spin_ops 

!>@brief \b subroutine interface for performing single DMRG step
!>@todo add current_site variable to controls, then use it for assignment to observables
    function singleStepDMRG_8(l_b,r_b,l_ss,r_ss,arpack,lapack,observable,controls) result(next)
        use   ::  lapack_simple
        use   ::  arpack_simple
        implicit none
        type(spin_ops), intent(in)          ::  l_b                     !< @param left block operators
        type(spin_ops), intent(in)          ::  r_b                     !< @param right block operators (identical with left block operator in IDMRG algorithm)
        type(spin_ops), intent(in)          ::  l_ss                    !< @param single site operator for left block
        type(spin_ops), intent(in)          ::  r_ss                    !< @param single site operator for right block
        type(DS_UPD), intent(inout)         ::  arpack                  !< @param variable containing all data for arpack run
        type(DSYEVX_L), intent(inout)       ::  lapack                  !< @param variable containing all data for lapack run
        type(dmrg_controls), intent(inout)  ::  controls                !< @param dmrg_controls defining parameters for a single step of the procedure
        type(spin_ops)                      ::  l_e                     !< left enlarged set of operators 
        type(spin_ops)                      ::  r_e                     !< right enlarged set of operators 
        type(observable_8), intent(inout)   ::  observable              !< @param observables, vector of observables for each site or 
        type(spin_ops)                      ::  next                    !< @result updated (rotated and truncated) operators \f$S_z\f$, \f$S_p\f$, \f$S_m\f$ and \f$H\f$
        type(tensor_2d_8)                   ::  O, lapmattemp, Sz_sup, Hsup  
        integer                             ::  prob_dim, l_arp_in, r_arp_in, l_arp_out, r_arp_out, a, b, cnt
        type :: subspace_indices
            sequence
            integer, dimension(:), allocatable  ::  list        !spin list in halves
            integer, dimension(:), allocatable  ::  l           !system index
            integer, dimension(:), allocatable  ::  r           !
        end type subspace_indices
        real*8, allocatable, dimension(:,:,:)   ::  psi1
        type(subspace_indices), dimension(-20:20)   ::  spin_basis
        integer, dimension(-20:20)                  ::  spin_basis_size, spin_basis_counter
        integer                             ::  total_spin, spin_selector
        real*8 , allocatable, dimension(:)  ::  eigvec_temp,eigvec_out
        spin_basis_size = 0
        spin_selector = controls % spin_selector
        l_e = enlargeBlock2(l_b,l_ss, controls) 
        if ( ( r_b .is. l_b ) .and. (l_ss .is. r_ss) ) then
            r_e = l_e
        else
            r_e = enlargeBlock2(r_b,r_ss, controls)
        endif
        
        call timeProfilingSimple(controls % profile)
        cnt=0
        do a=1, l_e % H % X
            do b=1, r_e % H % X
           !     write (*,'(f5.1f16.12)', advance = 'no') l_e % basis_Sz(a)+r_e % basis_Sz(b),&
            !    &arpack%EIGVEC((a-1)*l_e % H % X + b,1)
                total_spin = l_e % basis_Sz(a) + r_e % basis_Sz(b)
                if ( abs(total_spin) < 21 ) then
                    spin_basis_size(total_spin) = spin_basis_size(total_spin) + 1
                    cnt = cnt + 1 
                endif
            enddo
        enddo
!        print *,'number of coupled Sz basis states: ',spin_basis_size(0), ' out of total' , l_e % H % X * r_e % H % X , 'states'
        do a=-20,20
            allocate(spin_basis(a) % list(spin_basis_size(a)))
            allocate(spin_basis(a) % l(spin_basis_size(a)))
            allocate(spin_basis(a) % r(spin_basis_size(a)))
        enddo
        spin_basis_counter = 0
        do a=1, l_e % H % X
            do b=1, r_e % H % X
                total_spin = l_e % basis_Sz(a) + r_e % basis_Sz(b)
                if ( abs(total_spin) < 21 ) then
                    spin_basis_counter(total_spin) = spin_basis_counter(total_spin) +1
                    spin_basis(total_spin) % list(spin_basis_counter(total_spin)) = (a-1)*r_e % H % X + b
                    spin_basis(total_spin) % l(spin_basis_counter(total_spin)) = a  !not neccasary, can be computed from above
                    spin_basis(total_spin) % r(spin_basis_counter(total_spin)) = b  !not necesasry, can be computer
                endif
            enddo
        enddo
        do a=-5,5
!           print *,'tot spin ', a/2. , 'states' , spin_basis_size(a)
!            print *,'states: ', spin_basis(a) % list
        enddo
        call timeProfilingSimple(controls % profile,'restrictedBasisPreparation')
        
    !    prob_dim = l_e%H%X*r_e%H%X
        prob_dim = spin_basis_size(spin_selector)
        allocate(eigvec_temp( l_e%H%X*r_e%H%X ))
        allocate(eigvec_out( l_e%H%X*r_e%H%X ))

        call dsaupd_reinit(arpack, prob_dim, 1, .true., 20)
        do
            call arpack_iteration(arpack, controls)
            l_arp_in  = arpack%ipntr(1)
            l_arp_out = arpack%ipntr(2)

            call timeProfilingSimple(controls % profile)                 
            eigvec_temp = 0D0
            do a = 1, size(spin_basis(spin_selector)%list)
                eigvec_temp(spin_basis(spin_selector)%list(a)) = arpack%workd( l_arp_in + a  -1 )
            enddo
            call timeProfilingSimple(controls % profile, 'transformingArpackBasis')

            call tensorProdMatMul3(l_e,r_e,eigvec_temp,eigvec_out, controls)

            call timeProfilingSimple(controls % profile)                    !profiling tensor product matrix vector multiplication - start
            do a = 1, size(spin_basis(spin_selector)%list)
                arpack%workd( l_arp_out + a - 1 ) = eigvec_out(spin_basis(spin_selector)%list(a))
            enddo
            call timeProfilingSimple(controls % profile, 'transformingArpackBasis')

!            call tensorProdMatMul2(spin_basis,l_e,r_e, arpack%workd(l_arp_in:r_arp_in), arpack%workd(l_arp_out:r_arp_out), &
 !               &spin_selector) 

            if ( dsaupd_finished(arpack) .eqv. .true. ) exit
        enddo
!        print *,arpack%IPARAM
        call dseupd_reinit(arpack)                                      ! initialize and allocate variables for final eigenvector and eigenvalue evaluation
        call arpack_evaluate(arpack,controls)                           ! evaluate eigenvalue and eigenvector

        call reinitDSYEVX(lapack, l_e%H%x, controls % max_states,'last')
        call reallocate(lapmattemp,r_e%H%X,l_e%H%x)                         !initialize the size of lapmattemp temporary array object
        controls % ev = arpack%EIGVAL(1)    
        call timeProfilingSimple(controls % profile)
        eigvec_temp = 0D0
        do a = 1, size(spin_basis(spin_selector)%list)
            eigvec_temp(spin_basis(spin_selector)%list(a)) = arpack%EIGVEC(a,1)
        enddo

        lapmattemp % X = r_e % H % x
        lapmattemp % Y = l_e % H % x
        lapmattemp % M = reshape(eigvec_temp, [ lapmattemp % X , lapmattemp % Y])
        call timeProfilingSimple(controls % profile, 'copyingEigenvector')

        call timeProfilingSimple(controls % profile)
        lapack%A = matmul(transpose(lapmattemp%M), (lapmattemp%M))
        call timeProfilingSimple(controls % profile, 'matMulDensMatrCreate')

 !       call deallocate(lapmattemp)
        
        do a = 1, l_e % H % x                                           ! add random noise on diagonal
            lapack%A (a,a) = lapack%A (a,a) + (rand()-5D-1)*1D-8
        enddo
        call lapack_evaluate(lapack, controls)
        
        call reallocate(O, l_e%H%x, l_e%H%y, l_e%H%x, controls % max_states)           !initialize the size of O array object, let the Y dimension be not larger than max_states
        controls % truncation_err=1D0
        call timeProfilingSimple(controls % profile)
            do a=1,lapack%M
                O%M(:,a) = lapack%Z(:,lapack%M-a+1)
                controls % truncation_err = controls % truncation_err-lapack%W(lapack%M-a+1)
            enddo
        call timeProfilingSimple(controls % profile, 'sortingEigenvalues')

        call timeProfilingSimple(controls % profile)
            call updateAllOperators(next, l_e, O)
            next%basis_Sz = updateVector(l_e % basis_Sz,O)  !trzeba jeszcze obrócić i obciąć!
     !       lapmattemp % M = matmul(transpose(O % M),lapmattemp % M )
        call timeProfilingSimple(controls % profile, 'updatingOperators')
  !      allocate(psi1(lapmattemp % X , lapmattemp % Y/4,2))
  !      psi1 = reshape (lapmattemp % M, [ lapmattemp % X , lapmattemp % Y/4,2], order = [1,3,2])
        call timeProfilingSimple(controls % profile)
        if (observable % update .and. .not. allocated(observable % site % M )) observable % site = l_ss %Sz !single site should already be provided
        if (observable % update) then
            observable%site = (observable%site .expand.  l_ss%Sz%X) + (observable%site%X .expand. l_ss%Sz)
            if (observable%calculate) then
                observable%expected_val = vecTensorIdentityVecDotSymm(observable,eigvec_temp(:)) &
                    & + vecIdentityTensorVecDotSymm(observable,eigvec_temp(:))
            endif
            call updateSingleObservable(observable, O)
        endif
        call timeProfilingSimple(controls % profile, 'observableCalculations')

        call deallocate(arpack)
        call deallocate(lapack)
    end function singleStepDMRG_8
    
    subroutine lapack_evaluate(lapack, controls)
        use   ::  lapack_simple    
        type(DSYEVX_L), intent(inout)       ::  lapack                  !< @param variable containing all data for lapack run
        type(dmrg_controls), intent(inout)  ::  controls     
        call timeProfilingSimple(controls % profile)
        call omp_set_num_threads(controls % lapack_threads)                 !set optimal number of threads
        call runDSYEVX(lapack)                                              !run lapack
        call omp_set_num_threads(controls % omp_num_threads)                !set default number of threads
        call timeProfilingSimple(controls % profile, 'lapackRun')
    end subroutine lapack_evaluate
    
    subroutine arpack_iteration(arpack, controls)
        use   ::  arpack_simple    
        type(DS_UPD), intent(inout)         ::  arpack                  !< @param variable containing all data for arpack run
        type(dmrg_controls), intent(inout)  ::  controls
        call timeProfilingSimple(controls % profile)                    !profiling arpack - start
        call omp_set_num_threads(controls % arpack_threads)             !set optimal number of threads
        call dsaupd_run(arpack)                                         !run arpack
        call omp_set_num_threads(controls % omp_num_threads)            !set default number of threads
        call timeProfilingSimple(controls % profile, 'arpackRun')       !profiling arpack - stop
    end subroutine arpack_iteration
    
    subroutine arpack_evaluate(arpack, controls)
        use   ::  arpack_simple
        type(DS_UPD), intent(inout)         ::  arpack                  !< @param variable containing all data for arpack run
        type(dmrg_controls), intent(inout)  ::  controls
        call timeProfilingSimple(controls % profile)                        !profiling calculating eigenvalues and eigenvectors procedure - start
        call dseupd_run(arpack)
        call timeProfilingSimple(controls % profile, 'arpackFindingSolution')!profiling calculating eigenvalues and eigenvectors procedure - stop
    end subroutine arpack_evaluate
!>@brief Subroutine for updating (expanding, rotating and truncating) all observables with `update` field set to TRUE
    subroutine updateAllObservables_8(observables, O)
        type(observable_8), dimension(:), allocatable, intent(inout) :: observables  !< @param all observables, each containing vector of operators for all sites
        type(tensor_2d_8), intent(in) :: O                              !< @param rotation and truncation tensor_2d_8 object
        integer             ::  a
        if ( .not. allocated(observables) ) call die('observables not allocated!')
        do a=1, size(observables)
            if (observables(a) % update) call updateSingleObservable(observables(a),O)
        enddo
    end subroutine updateAllObservables_8
!>@brief Soubroutine for updating a single observable, the value of `update` field is ignored here
!>@todo its only a stub of a subroutine
    subroutine updateSingleObservable_8(observable, O)
        type(observable_8), intent(inout)  ::  observable               !< @param single observable containing vector with operators for all sites
        type(tensor_2d_8), intent(in) :: O                              !< @param rotation and truncation tensor_2d_8 object
        integer             ::  a
        if (allocated(observable % site % M )) then
             observable % site = rotateTruncate(observable % site, O)    
        endif
    end subroutine updateSingleObservable_8
    
    subroutine tensorProdMatMul2(spin_basis, l_b,r_b, vec_in, vec_out, spin_selector)
        use tensors
        implicit none
        type :: subspace_indices
            sequence
            integer, dimension(:), allocatable  ::  list        !spin list in halves
            integer, dimension(:), allocatable  ::  l           !system index
            integer, dimension(:), allocatable  ::  r           !
        end type subspace_indices
        type(subspace_indices), dimension(-20:20)   ::  spin_basis
        type(spin_ops), intent(in)          :: l_b, r_b
        real*8, dimension(:), intent(in)    :: vec_in
        real*8, dimension(:), intent(inout) :: vec_out
        integer, intent(in)                 :: spin_selector
        integer     ::  a,b,c,d,e,f
        real*8      ::  ZERO=1D-20
        vec_out = 0D0
        !$OMP PARALLEL DO NUM_THREADS(8) PRIVATE(a,b,e,f)
    do c=1, size(spin_basis(0)%list)
            a = spin_basis(0) % l(c)
            e = spin_basis(0) % r(c)
        do d=1, size(spin_basis(0)%list)
            b = spin_basis(0) % l(d)
            f = spin_basis(0) % r(d)
            !something like:
            if (a == b) vec_out(c) = vec_out(c) + r_b%H%M(f,e)*vec_in(d)
            if (e == f) vec_out(c) = vec_out(c) + l_b%H%M(b,a)*vec_in(d)
            
            vec_out(c) = vec_out(c) + ( &
                         &   + 5D-1*l_b%Sp%M(b,a) * r_b%Sm%M(f,e) &
                         &   + 5D-1*l_b%Sm%M(b,a) * r_b%Sp%M(f,e) &
                         &   + l_b%Sz%M(b,a) * r_b%Sz%M(f,e) &
                         &   ) * vec_in(d)
        enddo
    enddo
        !$OMP END PARALLEL DO

    end subroutine tensorProdMatMul2


    subroutine tensorProdMatMul3(l_b,r_b, vec_in, vec_out, controls)
        use tensors
        implicit none
        type(spin_ops), intent(in)          :: l_b, r_b
        real*8, dimension(:), intent(in)    :: vec_in
        real*8, dimension(:), intent(inout) :: vec_out
        type(dmrg_controls), intent(inout)  ::  controls                !< @param dmrg_controls defining parameters for a single step of the procedure
        integer     ::  a,b, l_in, r_in, l_out, r_out
        real*8      ::  ZERO=1D-20
        call timeProfilingSimple(controls % profile)                    !profiling tensor product matrix vector multiplication - start
        !$OMP PARALLEL DO NUM_THREADS(8) PRIVATE(l_out, r_out, l_in, r_in)
        do a=1, l_b%H%X
            l_out = 1+(a-1)*r_b%H%y
            r_out = a*r_b%H%Y
            vec_out(l_out:r_out) = matmul(r_b%H%M,vec_in(l_out:r_out))
            do b=1, l_b%H%Y
                l_in = 1+(b-1)*r_b%H%y
                r_in = b*r_b%H%Y
                if (abs(l_b%H%M(b,a))  > ZERO) vec_out(l_out:r_out) = vec_out(l_out:r_out) &
                                                    &+ l_b%H%M(b,a)*vec_in(l_in:r_in)
                if (abs(l_b%Sp%M(b,a)) > ZERO) vec_out(l_out:r_out) = vec_out(l_out:r_out) &
                        &    + 5D-1*l_b%Sp%M(b,a)*matmul(transpose(r_b%Sm%M),vec_in(l_in:r_in))
                if (abs(l_b%Sm%M(b,a)) > ZERO) vec_out(l_out:r_out) = vec_out(l_out:r_out) &
                        &    + 5D-1*l_b%Sm%M(b,a)*matmul(transpose(r_b%Sp%M),vec_in(l_in:r_in))
                if (abs(l_b%Sz%M(b,a)) > ZERO) vec_out(l_out:r_out) = vec_out(l_out:r_out) &
                        &    +      l_b%Sz%M(b,a)*matmul(transpose(r_b%Sz%M),vec_in(l_in:r_in))
            enddo
        enddo
        !$OMP END PARALLEL DO
    call timeProfilingSimple(controls % profile, 'matMulArpack')    !profiling tensor product matrix vector multiplication - stop
    end subroutine tensorProdMatMul3

    subroutine tensorProdMatMul(l_b,r_b, vec_in, vec_out, spin_selector)
        use tensors
        implicit none
        type(spin_ops), intent(in)          :: l_b, r_b
        real*8, dimension(:), intent(in)    :: vec_in
        real*8, dimension(:), intent(inout) :: vec_out
        integer, intent(in)                 :: spin_selector
        integer     ::  a,b, l_in, r_in, l_out, r_out
        real*8      ::  ZERO=1D-20
        !$OMP PARALLEL DO NUM_THREADS(8) PRIVATE(l_out, r_out, l_in, r_in)
        do a=1, l_b%H%X
            l_out = 1+(a-1)*r_b%H%y
            r_out = a*r_b%H%Y
            vec_out(l_out:r_out) = matmul(r_b%H%M,vec_in(l_out:r_out))
            do b=1, l_b%H%Y
                l_in = 1+(b-1)*r_b%H%y
                r_in = b*r_b%H%Y
                if (abs(l_b%H%M(b,a))  > ZERO) vec_out(l_out:r_out) = vec_out(l_out:r_out) &
                                                    &+ l_b%H%M(b,a)*vec_in(l_in:r_in)
                if (abs(l_b%Sp%M(b,a)) > ZERO) vec_out(l_out:r_out) = vec_out(l_out:r_out) &
                        &    + 5D-1*l_b%Sp%M(b,a)*matmul(transpose(r_b%Sm%M),vec_in(l_in:r_in))
                if (abs(l_b%Sm%M(b,a)) > ZERO) vec_out(l_out:r_out) = vec_out(l_out:r_out) &
                        &    + 5D-1*l_b%Sm%M(b,a)*matmul(transpose(r_b%Sp%M),vec_in(l_in:r_in))
                if (abs(l_b%Sz%M(b,a)) > ZERO) vec_out(l_out:r_out) = vec_out(l_out:r_out) &
                        &    +      l_b%Sz%M(b,a)*matmul(transpose(r_b%Sz%M),vec_in(l_in:r_in))
            enddo
        enddo
        !$OMP END PARALLEL DO
    end subroutine tensorProdMatMul


    function vecTensorIdentityVecDotSymm(observable, vec_in) result(dot)
        implicit none
        type(observable_8), intent(inout)       :: observable
        real*8, dimension(:), intent(in)        :: vec_in
        integer     ::  a,b, l_in, r_in, l_out, r_out
        real*8      ::  ZERO=1D-20, dot
        dot = 0D0
        !$OMP PARALLEL DO NUM_THREADS(8) PRIVATE(l_out, r_out, l_in, r_in) REDUCTION(+ : dot)
        do a=1, observable % site % X
            l_out = 1+(a-1)*observable % site % y
            r_out = a*observable % site % Y
            do b=1, observable % site % Y
                l_in = 1+(b-1)*observable % site % y
                r_in = b*observable % site % Y
                if (abs(observable % site % M(b,a)) > ZERO) dot = dot + &
                    &observable % site % M(b,a)*sum(vec_in(l_out:r_out)*vec_in(l_in:r_in)) !it can be improved, write vec_in as a matrix, do matmul, multiply with observable and do summation
            enddo
        enddo
        !$OMP END PARALLEL DO
    end function vecTensorIdentityVecDotSymm 
    
    function vecIdentityTensorVecDotSymm(observable, vec_in) result(dot)
        implicit none
        type(observable_8), intent(inout)       :: observable
        real*8, dimension(:), intent(in)        :: vec_in
        integer     ::  a, l_in, r_in
        real*8      ::  dot
        dot = 0D0
        !$OMP PARALLEL DO NUM_THREADS(8) PRIVATE(l_in, r_in) REDUCTION(+ : dot)
        do a=1, observable%site%X
            l_in = 1+(a-1)*observable%site%y
            r_in = a*observable%site%Y
            dot = dot + sum(vec_in(l_in:r_in)*matmul(transpose(observable%site%M),vec_in(l_in:r_in)))
        enddo
        !$OMP END PARALLEL DO
    end function vecIdentityTensorVecDotSymm 
    
!> @brief 
    subroutine updateAllOperators(truncated, enlarged, O)
        implicit none
        type(spin_ops), intent(in)          :: enlarged
        type(spin_ops), intent(inout)       :: truncated
        type(tensor_2d_8), intent(in) :: O
        
        call updateOperator(truncated%H,  enlarged%H,  O)
        call updateOperator(truncated%Sp, enlarged%Sp, O)
        call updateOperator(truncated%Sm, enlarged%Sm, O)
        call updateOperator(truncated%Sz, enlarged%Sz, O)
        
    end subroutine updateAllOperators

    function updateVector(vect_in, rot) result(vect_out)
        implicit none
        type(tensor_2d_8), intent(in)       ::  rot
        integer, dimension(:), intent(in)    ::  vect_in
        integer, dimension(:), allocatable   ::  vect_out
        integer :: a
        real*8  :: er, er2
        allocate(vect_out(rot%y))
        
        
        vect_out = nint(matmul(transpose(rot%M**2), 1D0*vect_in))
        
        do a =1, rot%Y
            er = sum(rot%M(:,a)**2*1D0*vect_in)
            er2 = abs(1D0*nint(er) - er)
            if ( abs(er2) > 1D-9) then
                        print *,'admixture from different spin states: ', abs(er2)

            endif
        enddo
    end function updateVector
    
    subroutine updateOperator(oper_out, oper_in, rot)
        implicit none
        type(tensor_2d_8), intent(in)    ::  oper_in, rot
        type(tensor_2d_8), intent(inout) ::  oper_out
        if (allocated(oper_out%M)) deallocate(oper_out%M)
        oper_out%X=rot%y
        oper_out%Y=rot%y
        allocate(oper_out%M(oper_out%x,oper_out%y))
        oper_out%M = matmul(transpose(rot%M), matmul(oper_in%M,rot%M))
    end subroutine updateOperator

    function rotateTruncate(oper_in, rot) result(oper_out)
        implicit none
        type(tensor_2d_8), intent(in)    ::  oper_in, rot
        type(tensor_2d_8)                ::  oper_out
        if (allocated(oper_out%M)) deallocate(oper_out%M)
        oper_out%X=rot%y
        oper_out%Y=rot%y
        allocate(oper_out%M(oper_out%x,oper_out%y))
        oper_out%M = matmul(transpose(rot%M), matmul(oper_in%M,rot%M))
    end function rotateTruncate
  
    subroutine reallocate(realloc, x, y, max_x, max_y)
        implicit none
        type(tensor_2d_8), intent(inout)  ::  realloc
        integer, intent(in)                     ::  x,y
        integer, optional, intent(in)           ::  max_x, max_y
        if ( x < 0 .or. y < 0 ) call die('reallocated array needs to have positive dimensions x and y!')
        if (allocated(realloc%M)) deallocate(realloc%M)
        realloc%X=x
        realloc%Y=y
        if (present(max_x)) then
            if ( max_x < x ) realloc%X=max_x
        endif
        if (present(max_y)) then
            if ( max_y < y ) realloc%Y=max_y
        endif
        allocate(realloc%M(realloc%X,realloc%Y))
    end subroutine reallocate  

!>brief przenieść do dmrg.f90
    subroutine enlargeBlock(cur,enlarged,ss)    
        implicit none
        type(spin_ops), intent(in)      ::  cur                             ! current list of operators
        type(spin_ops), intent(in)      ::  ss                              ! single site operators
        type(spin_ops), intent(inout)   ::  enlarged                        ! enlarged block operators
        
   !     call tensor_identity_product(cur%H, ss%H%X, enlarged%H)             ! He = H_b (x) 1_(n)
   !     call tensor_tensor_product(cur%Sp, ss%Sm, enlarged%H, '+', 0.5D0)   ! He = He + 0.5d0 * S_pb (x) S_m
   !     call tensor_tensor_product(cur%Sm, ss%Sp, enlarged%H, '+', 0.5D0)   ! He = He + 0.5d0 * S_mb (x) S_p
  !      call tensor_tensor_product(cur%Sz, ss%Sz, enlarged%H, '+', 1.0D0)   ! He = He + 1.0d0 * S_zb (x) S_z
 !       call identity_tensor_product(cur%H%X, ss%Sp, enlarged%Sp)           ! 1_(b) (x) Sp = Srpe
 !       call identity_tensor_product(cur%H%X, ss%Sm, enlarged%Sm)           ! 1_(b) (x) Sm = Srme
 !       call identity_tensor_product(cur%H%X, ss%Sz, enlarged%Sz)           ! 1_(b) (x) Sz = Srze
        
        enlarged%H = (cur%H .expand. ss%H%X) + .5*cur%Sp*ss%Sm + .5*cur%Sm*ss%Sp + cur%Sz*ss%Sz
        enlarged%Sp = cur%H%X .expand. ss%Sp
        enlarged%Sm = cur%H%X .expand. ss%Sm
        enlarged%Sz = cur%H%X .expand. ss%Sz
        
        !old H_superblock hamiltonian generation procedure
            !call tensor_identity_product(He, He%x, Hsup)                    ! Hsup = He (x) 1_(e)
            !call identity_tensor_product(He%x, He, Hsup, '+', 1.0D0)        ! Hsup = Hsup + 1. * 1_n (x) He
            !call tensor_tensor_product(Srpe, Srme, Hsup, '+', 0.5D0)        ! Hsup = Hsup + 0.5 * Srpe (x) Srme
            !call tensor_tensor_product(Srme, Srpe, Hsup, '+', 0.5D0)        ! Hsup = Hsup + 0.5 * Srme (x) Srpe
            !call tensor_tensor_product(Srze, Srze, Hsup, '+', 1.0D0)        ! Hsup = Hsup + 1.0 * Srze (x) Srze
    end subroutine enlargeBlock
!>@brief function for enlarging all spi operators and Sz of basis states
!>@todo there is an assumption regarding the size of the single site system, it is neccessary to remove it at some point, also move the basisi_Sz enlargement to another function
    function enlargeBlock2(cur,ss, controls) result(enlarged)
        implicit none
        type(spin_ops), intent(in)      ::  cur                             ! current list of operators
        type(spin_ops), intent(in)      ::  ss                              ! single site operators
        type(dmrg_controls), intent(inout)  ::  controls                !< @param dmrg_controls defining parameters for a single step of the procedure
        type(spin_ops)                  ::  enlarged                        ! enlarged block operators
        integer ::  a
        call timeProfilingSimple(controls % profile)
        enlarged%H = (cur%H .expand. ss%H%X) + .5*cur%Sp*ss%Sm + .5*cur%Sm*ss%Sp + cur%Sz*ss%Sz
        enlarged%Sp = cur%H%X .expand. ss%Sp
        enlarged%Sm = cur%H%X .expand. ss%Sm
        enlarged%Sz = cur%H%X .expand. ss%Sz
        allocate(enlarged%basis_Sz(cur%Sz%x * ss%Sz%x))
        do a = 1, cur%Sz%x
            enlarged%basis_Sz((a-1)*ss%Sz%x+1:a*ss%Sz%x) = ss%basis_Sz + cur%basis_Sz(a)
        enddo
        call timeProfilingSimple(controls % profile, 'enlargingBlocks')
    end function enlargeBlock2

    subroutine initTensor2d_8(m_in, what)
        type(tensor_2d_8), intent(out)  ::      m_in
        character(len=*), intent(in)    ::      what

        if ( allocated(m_in%M) )       deallocate(m_in%M)
        m_in%X=2
        m_in%Y=2
        allocate( m_in%M(m_in%X, m_in%Y) )
   
        if (trim(what) == 'Sz') then
            m_in%M(:,1)=[0.5D0, 0.0D0]
            m_in%M(:,2)=[0.0D0,-0.5D0]
        else if (trim(what) == 'Sp') then
            m_in%M(:,1)=[0.0D0, 1.0D0]
            m_in%M(:,2)=[0.0D0, 0.0D0]
        else if (trim(what) == 'Sm') then
            m_in%M(:,1)=[0.0D0, 0.0D0]
            m_in%M(:,2)=[1.0D0, 0.0D0]
        else if (trim(what) == 'Zero') then
            m_in%M=0D0
        endif
    end subroutine initTensor2d_8
  
    subroutine initTensor2d_4(m_in, what)
        type(tensor_2d_4), intent(out)  ::      m_in
        character(len=*), intent(in)    ::      what

        if ( allocated(m_in%M) )       deallocate(m_in%M)
        m_in%X=2
        m_in%Y=2
        allocate( m_in%M(m_in%X, m_in%Y) )
   
        if (trim(what) == 'Sz') then
            m_in%M(:,1)=[0.5E0, 0.0E0]
            m_in%M(:,2)=[0.0E0,-0.5E0]
        else if (trim(what) == 'Sp') then
            m_in%M(:,1)=[0.0E0, 1.0E0]
            m_in%M(:,2)=[0.0E0, 0.0E0]
        else if (trim(what) == 'Sm') then
            m_in%M(:,1)=[0.0E0, 0.0E0]
            m_in%M(:,2)=[1.0E0, 0.0E0]
        else if (trim(what) == 'Zero') then
            m_in%M=0E0
        endif
    end subroutine initTensor2d_4

    subroutine printSpinOps(ops)
        type(spin_ops), intent(in)     ::      ops
        print *,'Sz:'
        call print(ops%Sz)
        print *,'Sp:'
        call print(ops%Sp)
        print *,'Sm:'
        call print(ops%Sm)
        print *,'H:'
        call print(ops%H)
        print *
    end subroutine printSpinOps     

END MODULE dmrg
