MODULE tensors
  use :: utilities
  implicit none
  public
  

  type tensor_2d_8
    real*8, allocatable, dimension(:,:) ::  M
    integer                             ::  X
    integer                             ::  Y
  end type tensor_2d_8
    
   type tensor_2d_4
    real*4, allocatable, dimension(:,:) ::  M
    integer                             ::  X
    integer                             ::  Y
  end type tensor_2d_4
 
  type spin_ops
    type(tensor_2d_8) ::  Sz, Sp, Sm, H
  end type spin_ops
   
   
    interface assignment(=)
        module procedure initTensor2d_8, initTensor2d_4
    end interface
    
    interface print
        module procedure print_array, printSpinOps
    end interface
    
    interface deallocate
        module procedure deallocate_tensor_2d_4,deallocate_tensor_2d_8, deallocate_spin_ops
    end interface
    
    interface checkTensor
        module procedure checkTensor_tensor_2d_4, checkTensor_tensor_2d_8
    end interface checkTensor
    
    interface checkConformity
        module procedure checkConformity_tensor_2d_4, checkConformity_tensor_2d_8
    end interface checkConformity
    
    interface operator( .is. )
        module procedure isTheSame_tensor_2d_8,isTheSame_tensor_2d_4,isTheSame_spin_ops
    end interface
    
    interface operator( .expand. )
        module procedure identityTensorProduct_tensor_2d_4, identityTensorProduct_tensor_2d_8,&
            & tensorIdentityProduct_tensor_2d_4, tensorIdentityProduct_tensor_2d_8
    end interface 
    
    interface operator (*)
        module procedure tensorProduct_tensor_2d_4, tensorProduct_tensor_2d_8, &
                        &tensorMultiply_tensor_2d_8_r8, tensorMultiply_tensor_2d_4_r4,&
                        &tensorMultiply_tensor_2d_8_r4, tensorMultiply_tensor_2d_4_r8,&
                        &tensorMultiply_tensor_2d_8_i8, tensorMultiply_tensor_2d_4_i4,&
                        &tensorMultiply_tensor_2d_8_i4, tensorMultiply_tensor_2d_4_i8,&
                &multiplyTensor_r4_tensor_2d_8,multiplyTensor_r8_tensor_2d_8,&
                &multiplyTensor_r4_tensor_2d_4,multiplyTensor_r8_tensor_2d_4,&
                &multiplyTensor_i4_tensor_2d_8,multiplyTensor_i8_tensor_2d_8,&
                &multiplyTensor_i4_tensor_2d_4,multiplyTensor_i8_tensor_2d_4
    end interface
    
    interface operator (+)
        module procedure tensorAdd_tensor_2d_8, tensorAdd_tensor_2d_4
    end interface
    
    interface operator (-)
        module procedure tensorSubstract_tensor_2d_8, tensorSubstract_tensor_2d_4
    end interface
    
  contains
  
!przenieść do dmrg.f90
    subroutine enlargeBlock(cur,enlarged,ss)    
        implicit none
        type(spin_ops), intent(in)      ::  cur                             ! current list of operators
        type(spin_ops), intent(in)      ::  ss                              ! single site operators
        type(spin_ops), intent(inout)   ::  enlarged                        ! enlarged block operators
        
        call tensor_identity_product(cur%H, ss%H%X, enlarged%H)             ! He = H_b (x) 1_(n)
        call tensor_tensor_product(cur%Sp, ss%Sm, enlarged%H, '+', 0.5D0)   ! He = He + 0.0d0 * S_pb (x) S_m
        call tensor_tensor_product(cur%Sm, ss%Sp, enlarged%H, '+', 0.5D0)   ! He = He + 0.0d0 * S_mb (x) S_p
        call tensor_tensor_product(cur%Sz, ss%Sz, enlarged%H, '+', 1.0D0)   ! He = He + 1.0d0 * S_zb (x) S_z
        call identity_tensor_product(cur%H%X, ss%Sp, enlarged%Sp)           ! 1_(b) (x) Sp = Srpe
        call identity_tensor_product(cur%H%X, ss%Sm, enlarged%Sm)           ! 1_(b) (x) Sm = Srme
        call identity_tensor_product(cur%H%X, ss%Sz, enlarged%Sz)           ! 1_(b) (x) Sz = Srze
        
        !old H_superblock hamiltonian generation procedure
            !call tensor_identity_product(He, He%x, Hsup)                    ! Hsup = He (x) 1_(e)
            !call identity_tensor_product(He%x, He, Hsup, '+', 1.0D0)        ! Hsup = Hsup + 1. * 1_n (x) He
            !call tensor_tensor_product(Srpe, Srme, Hsup, '+', 0.5D0)        ! Hsup = Hsup + 0.5 * Srpe (x) Srme
            !call tensor_tensor_product(Srme, Srpe, Hsup, '+', 0.5D0)        ! Hsup = Hsup + 0.5 * Srme (x) Srpe
            !call tensor_tensor_product(Srze, Srze, Hsup, '+', 1.0D0)        ! Hsup = Hsup + 1.0 * Srze (x) Srze
    end subroutine enlargeBlock

!przenieść do dmrg.f90

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

!przenieść do dmrg.f90

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
  
    subroutine print_array(m_in)
        type(tensor_2d_8), intent(in)     ::      m_in
        integer                                 ::      a
        do a=1, m_in%Y
            print "(*(F6.2))",m_in%M(:,a)
        enddo
        print *
    end subroutine print_array
    subroutine printSpinOps(ops)
        type(spin_ops), intent(in)     ::      ops
        print *,'Sz:'
        call print_array(ops%Sz)
        print *,'Sp:'
        call print_array(ops%Sp)
        print *,'Sm:'
        call print_array(ops%Sm)
        print *,'H:'
        call print_array(ops%H)
        print *
    end subroutine printSpinOps
        
   
    subroutine identity_tensor_product(n, m_in, m_out, operation, mul)
        integer, intent(in)                     ::      n
        type(tensor_2d_8), intent(in)           ::      m_in
        type(tensor_2d_8), intent(inout)        ::      m_out
        integer                                 ::      a, x, y
        character(len=*), optional, intent(in)  ::      operation
        real*8, optional, intent(in)            ::      mul             !multiplicative constant for tensor multiplication out = 1_n * m_in * mul
        real*8                                  ::      mul_const
        call checkTensor(m_in)
        if ( .not. present(mul) ) then
            mul_const = 1D0
        else
            mul_const = mul
        endif
        if ( present(operation) .eqv. .true. ) then
            if (trim(operation) == '+') then
                if ( m_in%X*n /= m_out%X .or. m_in%Y*n /= m_out%Y ) &
                                                call die('Sizes of `m_in * identity matrix` and `m_out` do not conform!')
                if ( .not. allocated(m_out%M) ) call die('m_out needs to be allocated for `+` operation!')                                                
                do a=1, n
                    x = m_in%X
                    y = m_in%Y
                    m_out%M(    1+(a-1)*x : a*x   ,   1+(a-1)*y : a*y   ) = &
                        m_out%M(    1+(a-1)*x : a*x   ,   1+(a-1)*y : a*y   ) + mul_const*m_in%M
                enddo
            else
                call die('No such operation in `identity_tensor_product`')
            endif
        else
            if ( allocated(m_out%M) )       deallocate(m_out%M)
            m_out%X=m_in%X*n
            m_out%Y=m_in%Y*n
            allocate( m_out%M(m_out%X, m_out%Y) )
            m_out%M = 0D0
            do a=1, n
                x = m_in%X
                y = m_in%Y
                m_out%M(    1+(a-1)*x : a*x   ,   1+(a-1)*y : a*y   ) = m_in%M
            enddo
        endif
    end subroutine identity_tensor_product

    subroutine tensor_identity_product(m_in, n, m_out, operation, mul)
        implicit none
        integer, intent(in)                     ::      n
        type(tensor_2d_8), intent(in)           ::      m_in
        type(tensor_2d_8), intent(inout)        ::      m_out
        character(len=*), optional, intent(in)  ::      operation
        real*8, optional, intent(in)            ::      mul             !multiplicative constant for tensor multiplication out = 1_n * m_in * mul
        real*8                                  ::      mul_const
        integer                                 ::      a, b, c
        call checkTensor(m_in)
        if ( .not. present(mul) ) then
            mul_const = 1D0
        else
            mul_const = mul
        endif
                
        if ( present(operation) .eqv. .true. ) then
            if (trim(operation) == '+') then
                if ( m_in%X*n /= m_out%X .or. m_in%Y*n /= m_out%Y ) &
                                                call die('Sizes of `m_in * identity matrix` and `m_out` do not conform!')
                if ( .not. allocated(m_out%M) ) call die('m_out needs to be allocated for `+` operation!')                                                
                do a=1, m_in%X
                    do b=1, m_in%Y
                        do c=1,n
                            m_out%M(1+(a-1)*n+c-1,1+(b-1)*n+c-1) = m_out%M(1+(a-1)*n+c-1,1+(b-1)*n+c-1) + m_in%M(a,b)*mul_const
                        enddo
                    enddo
                enddo
            else
                call die('No such operation in `tensor_identity_product`')
            endif
        else
            if ( allocated(m_out%M) )       deallocate(m_out%M)
            m_out%X=m_in%X*n
            m_out%Y=m_in%Y*n
            allocate( m_out%M(m_out%X, m_out%Y) )
            m_out%M = 0D0
            do a=1, m_in%X
                do b=1, m_in%Y
                    do c=1,n
                        m_out%M(1+(a-1)*n+c-1,1+(b-1)*n+c-1) = m_in%M(a,b)
                    enddo
                enddo
            enddo
        endif   
    end subroutine tensor_identity_product

    subroutine tensor_tensor_product(m_in1, m_in2, m_out, operation, mul)
        type(tensor_2d_8), intent(in)     ::      m_in1, m_in2
        type(tensor_2d_8), intent(inout)  ::      m_out
        character(len=*), optional, intent(in)  ::      operation
        real*8, optional, intent(in)            ::      mul             !multiplicative constant for tensor multiplication out = 1_n * m_in * mul
        real*8                                  ::      mul_const        
        integer                                 ::      a, b, c, d
        call checkTensor(m_in1)
        call checkTensor(m_in2)
        if ( .not. present(mul) ) then
            mul_const = 1D0
        else
            mul_const = mul
        endif
                
        if ( present(operation) .eqv. .true. ) then
            if (trim(operation) == '+') then
                if ( m_in1%X*m_in2%X /= m_out%X .or. m_in1%Y*m_in2%Y /= m_out%Y ) &
                                                call die('Sizes of `m_in1 * m_in2` and `m_out` do not conform!')
                if ( .not. allocated(m_out%M) ) call die('m_out needs to be allocated for `+` operation!')                                                
                do a=1, m_in1%X
                    do b=1, m_in1%Y
                    if (abs(m_in1%M(a,b)) <1D-15) cycle
                        do c=1,m_in2%X
                            do d=1,m_in2%Y
                                m_out%M(1+(a-1)*m_in2%X+c-1,1+(b-1)*m_in2%Y+d-1) = &
                                    m_out%M(1+(a-1)*m_in2%X+c-1,1+(b-1)*m_in2%Y+d-1)  + mul_const * m_in1%M(a,b)*m_in2%M(c,d)
                            enddo
                        enddo
                    enddo
                enddo
            else
                call die('No such operation in `tensor_tensor_product`')
            endif
        else
            if ( allocated(m_out%M) )       deallocate(m_out%M)
            m_out%X=m_in1%X*m_in2%X
            m_out%Y=m_in1%Y*m_in2%Y
            allocate( m_out%M(m_out%X, m_out%Y) )
            m_out%M = 0D0
            do a=1, m_in1%X
                do b=1, m_in1%Y
                if (abs(m_in1%M(a,b)) <1D-15) cycle
                    do c=1,m_in2%X
                        do d=1,m_in2%Y
                            m_out%M(1+(a-1)*m_in2%X+c-1,1+(b-1)*m_in2%Y+d-1) = m_in1%M(a,b)*m_in2%M(c,d)
                        enddo
                    enddo
                enddo
            enddo
        endif
        
    end subroutine tensor_tensor_product
    
    function tensorProduct_tensor_2d_8(m_in1, m_in2) result (m_out)
        implicit none
        type(tensor_2d_8), intent(in)           ::      m_in1, m_in2
        type(tensor_2d_8)                       ::      m_out
        integer                                 ::      a, b, c, d
        print *,'prod'
        call checkTensor(m_in1)
        call checkTensor(m_in2)
        m_out%X=m_in1%X*m_in2%X
        m_out%Y=m_in1%Y*m_in2%Y
        allocate( m_out%M(m_out%X, m_out%Y) )
        m_out%M = 0D0
        do a=1, m_in1%X
            do b=1, m_in1%Y
            if (abs(m_in1%M(a,b)) <1D-20) cycle
                do c=1,m_in2%X
                    do d=1,m_in2%Y
                        m_out%M(1+(a-1)*m_in2%X+c-1,1+(b-1)*m_in2%Y+d-1) = m_in1%M(a,b)*m_in2%M(c,d)
                    enddo
                enddo
            enddo
        enddo
        
    end function tensorProduct_tensor_2d_8

    function tensorProduct_tensor_2d_4(m_in1, m_in2) result (m_out)
        implicit none
        type(tensor_2d_4), intent(in)           ::      m_in1, m_in2
        type(tensor_2d_4)                       ::      m_out
        integer                                 ::      a, b, c, d
        call checkTensor(m_in1)
        call checkTensor(m_in2)
        m_out%X=m_in1%X*m_in2%X
        m_out%Y=m_in1%Y*m_in2%Y
        allocate( m_out%M(m_out%X, m_out%Y) )
        m_out%M = 0D0
        do a=1, m_in1%X
            do b=1, m_in1%Y
            if (abs(m_in1%M(a,b)) <1D-15) cycle
                do c=1,m_in2%X
                    do d=1,m_in2%Y
                        m_out%M(1+(a-1)*m_in2%X+c-1,1+(b-1)*m_in2%Y+d-1) = m_in1%M(a,b)*m_in2%M(c,d)
                    enddo
                enddo
            enddo
        enddo
    end function tensorProduct_tensor_2d_4

    function identityTensorProduct_tensor_2d_8(n, m_in) result (m_out)
        implicit none
        integer, intent(in)               ::      n
        type(tensor_2d_8), intent(in)     ::      m_in
        type(tensor_2d_8)                 ::      m_out
        integer                           ::      a, x, y
        call checkTensor(m_in)
        m_out%X=m_in%X*n
        m_out%Y=m_in%Y*n
        allocate( m_out%M(m_out%X, m_out%Y) )
        m_out%M = 0D0
        do a=1, n
            x = m_in%X
            y = m_in%Y
            m_out%M(    1+(a-1)*x : a*x   ,   1+(a-1)*y : a*y   ) = m_in%M
        enddo
    end function identityTensorProduct_tensor_2d_8

    function identityTensorProduct_tensor_2d_4(n, m_in) result (m_out)
        implicit none
        integer, intent(in)               ::      n
        type(tensor_2d_4), intent(in)     ::      m_in
        type(tensor_2d_4)                 ::      m_out
        integer                           ::      a, x, y
        call checkTensor(m_in)
        m_out%X=m_in%X*n
        m_out%Y=m_in%Y*n
        allocate( m_out%M(m_out%X, m_out%Y) )
        m_out%M = 0D0
        do a=1, n
            x = m_in%X
            y = m_in%Y
            m_out%M(    1+(a-1)*x : a*x   ,   1+(a-1)*y : a*y   ) = m_in%M
        enddo
    end function identityTensorProduct_tensor_2d_4

    function tensorIdentityProduct_tensor_2d_8(m_in, n) result(m_out)
        implicit none
        integer, intent(in)               ::      n
        type(tensor_2d_8), intent(in)     ::      m_in
        type(tensor_2d_8)                 ::      m_out
        integer                           ::      a, b, c
        call checkTensor(m_in)
        m_out%X=m_in%X*n
        m_out%Y=m_in%Y*n
        allocate( m_out%M(m_out%X, m_out%Y) )
        m_out%M = 0D0
        do a=1, m_in%X
            do b=1, m_in%Y
                do c=1,n
                    m_out%M(1+(a-1)*n+c-1,1+(b-1)*n+c-1) = m_in%M(a,b)
                enddo
            enddo
        enddo
    end function tensorIdentityProduct_tensor_2d_8

    function tensorIdentityProduct_tensor_2d_4(m_in, n) result(m_out)
        implicit none
        integer, intent(in)               ::      n
        type(tensor_2d_4), intent(in)     ::      m_in
        type(tensor_2d_4)                 ::      m_out
        integer                           ::      a, b, c
        call checkTensor(m_in)
        m_out%X=m_in%X*n
        m_out%Y=m_in%Y*n
        allocate( m_out%M(m_out%X, m_out%Y) )
        m_out%M = 0D0
        do a=1, m_in%X
            do b=1, m_in%Y
                do c=1,n
                    m_out%M(1+(a-1)*n+c-1,1+(b-1)*n+c-1) = m_in%M(a,b)
                enddo
            enddo
        enddo
    end function tensorIdentityProduct_tensor_2d_4


    
    function tensorAdd_tensor_2d_8(m_in1, m_in2) result(m_out)
        implicit none
        type(tensor_2d_8), intent(in)     ::      m_in1, m_in2
        type(tensor_2d_8)                 ::      m_out
        call checkTensor(m_in1)
        call checkTensor(m_in2)
        call checkConformity(m_in1,m_in2)
        m_out%X=m_in1%X
        m_out%Y=m_in1%Y
        allocate( m_out%M(m_out%X, m_out%Y) )
        m_out%M = m_in1%M + m_in2%M
    end function tensorAdd_tensor_2d_8    

    function tensorAdd_tensor_2d_4(m_in1, m_in2) result(m_out)
        implicit none
        type(tensor_2d_4), intent(in)     ::      m_in1, m_in2
        type(tensor_2d_4)                 ::      m_out
        call checkTensor(m_in1)
        call checkTensor(m_in2)
        call checkConformity(m_in1,m_in2)
        m_out%X=m_in1%X
        m_out%Y=m_in1%Y
        allocate( m_out%M(m_out%X, m_out%Y) )
        m_out%M = m_in1%M + m_in2%M
    end function tensorAdd_tensor_2d_4  

    function tensorSubstract_tensor_2d_8(m_in1, m_in2) result(m_out)
        implicit none
        type(tensor_2d_8), intent(in)     ::      m_in1, m_in2
        type(tensor_2d_8)                 ::      m_out
        call checkTensor(m_in1)
        call checkTensor(m_in2)
        call checkConformity(m_in1,m_in2)
        m_out%X=m_in1%X
        m_out%Y=m_in1%Y
        allocate( m_out%M(m_out%X, m_out%Y) )
        m_out%M = m_in1%M - m_in2%M
    end function tensorSubstract_tensor_2d_8    

    function tensorSubstract_tensor_2d_4(m_in1, m_in2) result(m_out)
        implicit none
        type(tensor_2d_4), intent(in)     ::      m_in1, m_in2
        type(tensor_2d_4)                 ::      m_out
        call checkTensor(m_in1)
        call checkTensor(m_in2)
        call checkConformity(m_in1,m_in2)
        m_out%X=m_in1%X
        m_out%Y=m_in1%Y
        allocate( m_out%M(m_out%X, m_out%Y) )
        m_out%M = m_in1%M - m_in2%M
    end function tensorSubstract_tensor_2d_4  
    
    function tensorMultiply_tensor_2d_8_r8(m_in, mul) result(m_out)
        implicit none
        real*8, intent(in)                ::      mul
        type(tensor_2d_8), intent(in)     ::      m_in
        type(tensor_2d_8)                 ::      m_out
        call checkTensor(m_in)
        m_out%X=m_in%X
        m_out%Y=m_in%Y
        allocate( m_out%M(m_out%X, m_out%Y) )
        m_out%M = mul * m_in%M
    end function tensorMultiply_tensor_2d_8_r8
        
    function tensorMultiply_tensor_2d_4_r4(m_in, mul) result(m_out)
        implicit none
        real*4, intent(in)                ::      mul
        type(tensor_2d_4), intent(in)     ::      m_in
        type(tensor_2d_4)                 ::      m_out
        call checkTensor(m_in)
        m_out%X=m_in%X
        m_out%Y=m_in%Y
        allocate( m_out%M(m_out%X, m_out%Y) )
        m_out%M = mul * m_in%M
    end function tensorMultiply_tensor_2d_4_r4

    function tensorMultiply_tensor_2d_8_r4(m_in, mul) result(m_out)
        implicit none
        real*4, intent(in)                ::      mul
        type(tensor_2d_8), intent(in)     ::      m_in
        type(tensor_2d_8)                 ::      m_out
        call checkTensor(m_in)
        m_out%X=m_in%X
        m_out%Y=m_in%Y
        allocate( m_out%M(m_out%X, m_out%Y) )
        m_out%M = mul * m_in%M
    end function tensorMultiply_tensor_2d_8_r4

    function tensorMultiply_tensor_2d_4_r8(m_in, mul) result(m_out)
        implicit none
        real*8, intent(in)                ::      mul
        type(tensor_2d_4), intent(in)     ::      m_in
        type(tensor_2d_4)                 ::      m_out
        call checkTensor(m_in)
        m_out%X=m_in%X
        m_out%Y=m_in%Y
        allocate( m_out%M(m_out%X, m_out%Y) )
        m_out%M = real(mul,4) * m_in%M
    end function tensorMultiply_tensor_2d_4_r8    

    function tensorMultiply_tensor_2d_8_i8(m_in, mul) result(m_out)
        implicit none
        integer*8, intent(in)             ::      mul
        type(tensor_2d_8), intent(in)     ::      m_in
        type(tensor_2d_8)                 ::      m_out
        m_out = tensorMultiply_tensor_2d_8_r8(m_in, real(mul,8))
    end function tensorMultiply_tensor_2d_8_i8
        
    function tensorMultiply_tensor_2d_8_i4(m_in, mul) result(m_out)
        implicit none
        integer*4, intent(in)             ::      mul
        type(tensor_2d_8), intent(in)     ::      m_in
        type(tensor_2d_8)                 ::      m_out
        m_out = tensorMultiply_tensor_2d_8_r4(m_in, real(mul,4))
    end function tensorMultiply_tensor_2d_8_i4

    function tensorMultiply_tensor_2d_4_i8(m_in, mul) result(m_out)
        implicit none
        integer*8, intent(in)             ::      mul
        type(tensor_2d_4), intent(in)     ::      m_in
        type(tensor_2d_4)                 ::      m_out
        m_out = tensorMultiply_tensor_2d_4_r8(m_in, real(mul,8))
    end function tensorMultiply_tensor_2d_4_i8

    function tensorMultiply_tensor_2d_4_i4(m_in, mul) result(m_out)
        implicit none
        integer*4, intent(in)             ::      mul
        type(tensor_2d_4), intent(in)     ::      m_in
        type(tensor_2d_4)                 ::      m_out
        m_out = tensorMultiply_tensor_2d_4_r4(m_in, real(mul,4))
    end function tensorMultiply_tensor_2d_4_i4  

    function multiplyTensor_r4_tensor_2d_8(mul, m_in) result(m_out)
        implicit none
        real*4, intent(in)                ::      mul
        type(tensor_2d_8), intent(in)     ::      m_in
        type(tensor_2d_8)                 ::      m_out
        m_out = tensorMultiply_tensor_2d_8_r4(m_in, mul)
    end function multiplyTensor_r4_tensor_2d_8

    function multiplyTensor_r8_tensor_2d_8(mul, m_in) result(m_out)
        implicit none
        real*8, intent(in)                ::      mul
        type(tensor_2d_8), intent(in)     ::      m_in
        type(tensor_2d_8)                 ::      m_out
        m_out = tensorMultiply_tensor_2d_8_r8(m_in, mul)
    end function multiplyTensor_r8_tensor_2d_8

    function multiplyTensor_r4_tensor_2d_4(mul, m_in) result(m_out)
        implicit none
        real*4, intent(in)                ::      mul
        type(tensor_2d_4), intent(in)     ::      m_in
        type(tensor_2d_4)                 ::      m_out
        m_out = tensorMultiply_tensor_2d_4_r4(m_in, mul)
    end function multiplyTensor_r4_tensor_2d_4

    function multiplyTensor_r8_tensor_2d_4(mul, m_in) result(m_out)
        implicit none
        real*8, intent(in)                ::      mul
        type(tensor_2d_4), intent(in)     ::      m_in
        type(tensor_2d_4)                 ::      m_out
        m_out = tensorMultiply_tensor_2d_4_r8(m_in, mul)
    end function multiplyTensor_r8_tensor_2d_4

    function multiplyTensor_i4_tensor_2d_8(mul, m_in) result(m_out)
        implicit none
        integer*4, intent(in)             ::      mul
        type(tensor_2d_8), intent(in)     ::      m_in
        type(tensor_2d_8)                 ::      m_out
        m_out = tensorMultiply_tensor_2d_8_r4(m_in, real(mul,4))
    end function multiplyTensor_i4_tensor_2d_8

    function multiplyTensor_i8_tensor_2d_8(mul, m_in) result(m_out)
        implicit none
        integer*8, intent(in)             ::      mul
        type(tensor_2d_8), intent(in)     ::      m_in
        type(tensor_2d_8)                 ::      m_out
        m_out = tensorMultiply_tensor_2d_8_r8(m_in, real(mul,8))
    end function multiplyTensor_i8_tensor_2d_8

    function multiplyTensor_i4_tensor_2d_4(mul, m_in) result(m_out)
        implicit none
        integer*4, intent(in)             ::      mul
        type(tensor_2d_4), intent(in)     ::      m_in
        type(tensor_2d_4)                 ::      m_out
        m_out = tensorMultiply_tensor_2d_4_r4(m_in, real(mul,4))
    end function multiplyTensor_i4_tensor_2d_4

    function multiplyTensor_i8_tensor_2d_4(mul, m_in) result(m_out)
        implicit none
        integer*8, intent(in)             ::      mul
        type(tensor_2d_4), intent(in)     ::      m_in
        type(tensor_2d_4)                 ::      m_out
        m_out = tensorMultiply_tensor_2d_4_r8(m_in, real(mul,8))
    end function multiplyTensor_i8_tensor_2d_4

    function isTheSame_tensor_2d_4(tensor1,tensor2) result(same)
        implicit none
        type(tensor_2d_4), intent(in)   ::  tensor1, tensor2
        logical                         ::  same
        same=.false.
        if (loc(tensor1) == loc(tensor2)) same=.true.
    end function isTheSame_tensor_2d_4
    
    function isTheSame_tensor_2d_8(tensor1,tensor2) result(same)
        implicit none
        type(tensor_2d_8), intent(in)   ::  tensor1, tensor2
        logical                         ::  same
        same=.false.
        if (loc(tensor1) == loc(tensor2)) same=.true.
    end function isTheSame_tensor_2d_8

    function isTheSame_spin_ops(ops1,ops2) result(same)
        implicit none
        type(spin_ops), intent(in)      ::  ops1, ops2
        logical                         ::  same
        same=.false.
        if (loc(ops1) == loc(ops2)) same=.true.
    end function isTheSame_spin_ops
    
    subroutine deallocate_tensor_2d_4(tensor1)
        implicit none
        type(tensor_2d_4), intent(inout)   ::  tensor1    
        if (allocated(tensor1%M)) deallocate(tensor1%M)
    end subroutine deallocate_tensor_2d_4

    subroutine deallocate_tensor_2d_8(tensor1)
        implicit none
        type(tensor_2d_8), intent(inout)   ::  tensor1    
        if (allocated(tensor1%M)) deallocate(tensor1%M)
    end subroutine deallocate_tensor_2d_8 

    subroutine deallocate_spin_ops(ops)
        implicit none
        type(spin_ops), intent(inout)   ::  ops    
        call deallocate_tensor_2d_8(ops%Sp)
        call deallocate_tensor_2d_8(ops%Sm)
        call deallocate_tensor_2d_8(ops%Sz)
        call deallocate_tensor_2d_8(ops%H)
    end subroutine deallocate_spin_ops 

    subroutine checkConformity_tensor_2d_8(t1,t2)!trzeba poprzerabiać to na funkcje, bo debugowanie będzie bardzo utrudnione
        implicit none
        type(tensor_2d_8)   ::  t1,t2
        if ( t1%X /= t2%X .or. t1%Y /= t2%Y) call die('Sizes do not conform!')
    end subroutine checkConformity_tensor_2d_8
    
    subroutine checkConformity_tensor_2d_4(t1,t2)
        implicit none
        type(tensor_2d_4)   ::  t1,t2
        if ( t1%X /= t2%X .or. t1%Y /= t2%Y) call die('Sizes do not conform!')
    end subroutine checkConformity_tensor_2d_4
                                                        
    subroutine checkTensor_tensor_2d_8(t)
        implicit none
        type(tensor_2d_8)   ::  t
        if ( .not. allocated(t%M) )        call die('Structure tensor_2d_8 not allocated')
        if ( size(t % M) /= t % X * t % Y )call die('Structure dimensions do not conform with factual size')
    end subroutine checkTensor_tensor_2d_8

    subroutine checkTensor_tensor_2d_4(t)
        implicit none
        type(tensor_2d_4)   ::  t
        if ( .not. allocated(t%M) )        call die('Structure tensor_2d_4 not allocated')
        if ( size(t % M) /= t % X * t % Y )call die('Structure dimensions do not conform with factual size')
    end subroutine checkTensor_tensor_2d_4

        
END MODULE tensors
