program su2
    use omp_lib
    implicit none
    ! define the precision
    integer, parameter :: r8=selected_real_kind(8,15)
    integer, parameter :: i4=selected_int_kind(8)

    ! define the parameters
    real(kind=r8), parameter :: pi=acos(-1.0_r8)
    !integer(kind=i4), parameter :: nr=10 ! number of sites in space
    !integer(kind=i4), parameter :: nt=10 !  number of sites in time

    ! define the matrices
    real(kind=r8), allocatable, dimension(:,:,:,:,:,:) :: ulink, ugauge,ur,z
    real(kind=r8), allocatable, dimension(:,:) :: w, wr, wv, we, wo

    ! defining our variables
    real(kind=r8) :: beta, tol
    integer(kind=i4) :: nr, nt, rate,itimes1,itimes2
    integer(kind=i4) :: i, j, nmc, nterm, ncorr, iprint, a, b
    character(len=1) :: function
    !call cpu_time(start)
    call system_clock(count_rate=rate)
    call system_clock(itimes1)
    open(unit=1,file='su2-in.dat')
    read(1,*) nr, nt            ! lattice size space time
    read(1,*) beta              ! beta
    read(1,*) nmc, nterm        ! number of monte carlo sweeps number of thermalization sweeps
    read(1,*) a, b              ! maximum size of the wilson loops we want
    read(1,*) iprint            ! if 1 we print the configurations on the file suN-lattice
    read(1,*) ncorr             ! number of discarted configurations
    read(1,*) tol               ! gauge fix tolerance
    read(1,*) function
    close(1)

    allocate(ulink(nr,nr,nr,nt,4,4),ur(nr,nr,nr,nt,4,4))
    allocate(ugauge(nr,nr,nr,nt,4,4),z(nr,nr,nr,nt,4,4))
    allocate(w(a,b),wv(a,b),wr(a,b),we(a,b),wo(a,b))
    call init(ulink,beta,2)
    
    select case(function)
    case('a') ! termilize the lattice and find the correlation time
        ! termilize the lattice
        do i=1,nterm
            call hbstep(ulink,beta)
        enddo
        do i=1,nmc
            call hbstep(ulink,beta)
            call measurewilson(ulink,1,1,w)
            write(101,*) i, w(1,1)
        enddo
    case('b') ! termilize the lattice and generate uncorrelated configurations
        do i=1,nterm
            call hbstep(ulink,beta)
        enddo
        do i=1,nmc
            ! save the actual configuration
            if(iprint.eq.1)call save_lattice(ulink,nr,i)

            ! compute the desired observables
            call maximal_center_gauge(ulink,ugauge,tol)
            call centerprojection(ugauge,z)
            call centerremotion(ugauge,ur)
            call measurewilson(ugauge,a,b,w)
            call measurewilson(z,a,b,wv)
            call measurewilson(ur,a,b,wr)
            call pvortexdensity(z,a,b)
            call vortexsperatedwilson(ugauge,z,a,b,we,wo)
            call vortexcomponets(w,we,wo,a,b)
            call vortexlimitedwilson(ugauge,z,a,b)
            !call wilsonvortex(ugauge,ur,z,w,wr,wv,we,wo,a,b)
            
            do j=1,a
                write(100+j,*) i, w(j,:)
                write(200+j,*) i, wv(j,:)
                write(300+j,*) i, wr(j,:)
            enddo
            ! get the next uncorrelated configurations
            do j=1,ncorr ! number of correlated configurations
                ! discart ncorr configurations
                call hbstep(ulink,beta)
            enddo
        enddo
    case('c')   ! read the configurations and make measures
        do i=1,nmc-1
            ! save the actual configuration
            call read_lattice(ulink,nr,i)
            call maximal_center_gauge(ulink,ugauge,tol)
            call centerprojection(ugauge,z)
            call centerremotion(ugauge,ur)
            call wilsonvortex(ugauge,ur,z,w,wr,wv,we,wo,a,b)
        enddo
    case default
        stop 'no valid case was selected'
    endselect

    deallocate(ulink,ugauge,z,wr,w,wv)
    !call cpu_time(end)
    call system_clock(itimes2)
    print*,"cpu time [seconds]: ", (itimes2-itimes1)/float(rate)!(end-start)/60
    write(9000,*) "cpu time [seconds]: ", (itimes2-itimes1)/float(rate)
    write(9000,*) "cpu time [minutes]: ", (itimes2-itimes1)/float(60*rate)
    write(9000,*) "cpu time [hours]: ", (itimes2-itimes1)/float(3600*rate)
    contains
    function ident()
        real(kind=r8) :: ident(4)
        ident=(/0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8/)
    endfunction ident

    function crossprod(a,b)
        real(kind=r8), dimension(4) :: a, b
        real(kind=r8), dimension(3) ::  c, crossprod
        c(1)=a(2)*b(3)-a(3)*b(2)
        c(2)=a(3)*b(1)-a(1)*b(3)
        c(3)=a(1)*b(2)-a(2)*b(1)
        crossprod=c
    endfunction crossprod

    function dot3vec(a,b)
        real(kind=r8), dimension(4) :: a, b
        real(kind=r8):: dot3vec
        dot3vec=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
    endfunction dot3vec

    function linkmult(a,b)
        real(kind=r8), dimension(4) :: a, b, c, linkmult
        real(kind=r8) :: ab(3)
        ab=crossprod(a,b)
        c(4)=a(4)*b(4)-dot3vec(a,b)

        c(1)=b(4)*a(1)+a(4)*b(1)-ab(1)
        c(2)=b(4)*a(2)+a(4)*b(2)-ab(2)
        c(3)=b(4)*a(3)+a(4)*b(3)-ab(3)
        linkmult=c
    endfunction linkmult

    function trlink(u)
        real(kind=r8) :: trlink, u(4)
        trlink=2.0_r8*u(4)
    endfunction trlink

    ! function that computes the determinant of su(2) matrices
    function detlink(u)
        real(kind=r8) :: detlink, u(4)
        detlink=u(1)*u(1)+u(2)*u(2)+u(3)*u(3)+u(4)*u(4)
    endfunction detlink

    ! function that returns the dagger of the link
    function dlink(u)
        real(kind=r8), dimension(4) :: dlink, u, v
        v(1)=-u(1)
        v(2)=-u(2)
        v(3)=-u(3)
        v(4)=u(4)
        dlink=v
    endfunction dlink

    ! subroutine that returns the neighbors of a given link
    ! ONLY FOR HYPERCUBIC LATTICES
    subroutine neighborhood(nzz,npz,nmz,nzp,nzm,npp,npm,nmp,mi,ni)
        integer(kind=i4), dimension(4) :: nzz,npz,nmz,nzp,nzm,npp,npm,nmp
        integer(kind=i4) :: mi, ni
        ! start the neighboors vectors
        npz=nzz
        nmz=nzz
        nzp=nzz
        nzm=nzz
        npp=nzz
        npm=nzz
        nmp=nzz
        ! in mi direction
        ! plus mi
        npz(mi)=npz(mi)+1
        npp(mi)=npp(mi)+1
        npm(mi)=npm(mi)+1
        if(npz(mi).gt.nr)then
            npz(mi)=1
            npp(mi)=1
            npm(mi)=1
        endif
        ! minus mi
        nmz(mi)=nmz(mi)-1
        nmp(mi)=nmp(mi)-1
        if(nmz(mi).lt.1)then
            nmz(mi)=nr
            nmp(mi)=nr
        endif

        ! in ni direction
        ! plus ni
        nzp(ni)=nzp(ni)+1
        npp(ni)=npp(ni)+1
        nmp(ni)=nmp(ni)+1
        if(nzp(ni).gt.nr)then
            nzp(ni)=1
            npp(ni)=1
            nmp(ni)=1
        endif

        ! minus ni
        nzm(ni)=nzm(ni)-1
        npm(ni)=npm(ni)-1
        if(nzm(ni).lt.1)then
            nzm(ni)=nr
            npm(ni)=nr
        endif
    endsubroutine neighborhood

    ! function that computes the plaquette
    function plaquette(u)
        real(kind=r8) :: u(nr,nr,nr,nt,4,4), sum, w(4), plaquette
        integer(kind=i4) :: e1, e2, e3, e4, mi, ni
        integer(kind=i4), dimension(4) :: nzz,npz,nmz,nzp,nzm,npp,npm,nmp
        sum=0.0_r8
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nt
            nzz=(/e1,e2,e3,e4/)
            do mi=1,3
            do ni=mi+1,4
                call neighborhood(nzz,npz,nmz,nzp,nzm,npp,npm,nmp,mi,ni)
                w=ident()
                !w=linkmult(w,dlink(u(nzz(1),nzz(2),nzz(3),nzz(4),ni,:)))
                !w=linkmult(w,dlink(u(nzp(1),nzp(2),nzp(3),nzp(4),mi,:)))
                !w=linkmult(w,u(npz(1),npz(2),npz(3),npz(4),ni,:))
                !w=linkmult(w,u(nzz(1),nzz(2),nzz(3),nzz(4),mi,:))
                w=linkmult(w,u(nzz(1),nzz(2),nzz(3),nzz(4),mi,:))
                w=linkmult(w,u(npz(1),npz(2),npz(3),npz(4),ni,:))
                w=linkmult(w,dlink(u(nzp(1),nzp(2),nzp(3),nzp(4),mi,:)))
                w=linkmult(w,dlink(u(nzz(1),nzz(2),nzz(3),nzz(4),ni,:)))
                sum=sum+w(4)
            enddo
            enddo
        enddo
        enddo
        enddo
        enddo
        plaquette=sum/(6.0_r8*nr*nr*nr*nt)
    endfunction plaquette

    ! funtion that compute the sum of the stamples
    function sumstamples(u,e1,e2,e3,e4,mi)
        real(kind=r8) :: u(nr,nr,nr,nt,4,4)
        real(kind=r8), dimension(4) :: sumstamples, temp, sum
        integer(kind=i4) :: e1, e2, e3, e4, mi, ni
        integer(kind=i4), dimension(4) :: nzz, npz, nmz, nzp, nzm, npp, npm, nmp

        nzz=(/e1,e2,e3,e4/)
        sum=0.0_r8
        do ni=1,4
            if (ni.ne.mi) then
                call neighborhood(nzz,npz,nmz,nzp,nzm,npp,npm,nmp,mi,ni)
                temp=ident()
                temp=linkmult(temp,u(npz(1),npz(2),npz(3),npz(4),ni,:))
                temp=linkmult(temp,dlink(u(nzp(1),nzp(2),nzp(3),nzp(4),mi,:)))
                temp=linkmult(temp,dlink(u(nzz(1),nzz(2),nzz(3),nzz(4),ni,:)))
                
                sum=sum+temp
                temp=ident()
                temp=linkmult(temp,dlink(u(npm(1),npm(2),npm(3),npm(4),ni,:)))
                temp=linkmult(temp,dlink(u(nzm(1),nzm(2),nzm(3),nzm(4),mi,:)))
                temp=linkmult(temp,u(nzm(1),nzm(2),nzm(3),nzm(4),ni,:))
                sum=sum+temp
            endif
        enddo
        sumstamples=sum
    endfunction sumstamples

    ! function that computes the heatbath matrix
    function hbmatrix(a,beta)
        real(kind=r8) :: hbmatrix(4), x(4), r(3)
        real(kind=r8) :: lambda2, rr, a, beta
        integer(kind=i4) :: i, j

        do j=1,200
            do i=1,3
                call random_number(r(i))
                r(i)=1.0_r8-r(i)
            enddo
            lambda2=-0.5_r8*(dlog(r(1))+(dcos(2.0_r8*pi*r(2)))**2.0_r8*dlog(r(3)))/(a*beta)
            call random_number(rr)
            rr=1.0_r8-rr
            if(rr*rr.lt.(1.0_r8-lambda2))exit
        enddo
        x(4)=1.0_r8-2.0_r8*lambda2

        do j=1,2000
            do i=1,3
                call random_number(r(i))
                r(i)=1.0_r8-2.0_r8*r(i)
            enddo
            if((r(1)*r(1)+r(2)*r(2)+r(3)*r(3)).lt.1.0_r8)exit
        enddo

        do i=1,3
            x(i)=dsqrt((1.0_r8-x(4)*x(4)))*r(i)/dsqrt((r(1)*r(1)+r(2)*r(2)+r(3)*r(3)))
        enddo
        hbmatrix=x
    endfunction hbmatrix

    ! function that makes the heatbath steap
    subroutine hbstep(u,beta)
        real(kind=r8), intent(inout) :: u(nr,nr,nr,nt,4,4)
        real(kind=r8) ::  sums(4), a, beta
        integer(kind=i4) :: e1, e2, e3, e4, mi

        !=============================================================================!
        
        !$omp parallel sections private(e1,e2,e3,e4,mi,a,sums)
        !$omp section
        mi=1
        do e1=2,nr,2
        do e2=1,nr
        do e3=1,nr
        do e4=1,nt
            !!$omp ordered
                sums=sumstamples(u,e1,e2,e3,e4,mi)
                a=dsqrt(detlink(sums))
                u(e1,e2,e3,e4,mi,:)=linkmult(hbmatrix(a,beta),dlink(sums))/a
                u(e1,e2,e3,e4,mi,:)=u(e1,e2,e3,e4,mi,:)/dsqrt(detlink(u(e1,e2,e3,e4,mi,:)))
            !!$omp end ordered
        enddo
        enddo
        enddo
        enddo
        !!$omp end parallel do

        !!$omp parallel do ordered
        do e1=1,nr-1,2
        do e2=1,nr
        do e3=1,nr
        do e4=1,nt
         !   !$omp ordered
                sums=sumstamples(u,e1,e2,e3,e4,mi)
                a=dsqrt(detlink(sums))
                u(e1,e2,e3,e4,mi,:)=linkmult(hbmatrix(a,beta),dlink(sums))/a
                u(e1,e2,e3,e4,mi,:)=u(e1,e2,e3,e4,mi,:)/dsqrt(detlink(u(e1,e2,e3,e4,mi,:)))
          !  !$omp end ordered
        enddo
        enddo
        enddo
        enddo
        !!$omp end parallel do
        !=============================================================================!
        !$omp section
        mi=2
        !!$omp parallel do ordered
        do e1=1,nr
        do e2=2,nr,2
        do e3=1,nr
        do e4=1,nt
            !!$omp ordered
                sums=sumstamples(u,e1,e2,e3,e4,mi)
                a=dsqrt(detlink(sums))
                u(e1,e2,e3,e4,mi,:)=linkmult(hbmatrix(a,beta),dlink(sums))/a
                u(e1,e2,e3,e4,mi,:)=u(e1,e2,e3,e4,mi,:)/dsqrt(detlink(u(e1,e2,e3,e4,mi,:)))
            !!$omp end ordered
        enddo
        enddo
        enddo
        enddo
        !!$omp end parallel do

        do e1=1,nr
        do e2=1,nr-1,2
        do e3=1,nr
        do e4=1,nt
                sums=sumstamples(u,e1,e2,e3,e4,mi)
                a=dsqrt(detlink(sums))
                u(e1,e2,e3,e4,mi,:)=linkmult(hbmatrix(a,beta),dlink(sums))/a
                u(e1,e2,e3,e4,mi,:)=u(e1,e2,e3,e4,mi,:)/dsqrt(detlink(u(e1,e2,e3,e4,mi,:)))
        enddo
        enddo
        enddo
        enddo
        !=============================================================================!
        !$omp section
        mi=3
        !!$omp parallel do ordered
        do e1=1,nr
        do e2=1,nr
        do e3=2,nr,2
        do e4=1,nt
         !   !$omp ordered
                sums=sumstamples(u,e1,e2,e3,e4,mi)
                a=dsqrt(detlink(sums))
                u(e1,e2,e3,e4,mi,:)=linkmult(hbmatrix(a,beta),dlink(sums))/a
                u(e1,e2,e3,e4,mi,:)=u(e1,e2,e3,e4,mi,:)/dsqrt(detlink(u(e1,e2,e3,e4,mi,:)))
          !  !$omp end ordered
        enddo
        enddo
        enddo
        enddo
        !!$omp end parallel do

        !!$omp parallel do ordered
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr-1,2
        do e4=1,nt
         !   !$omp ordered
                sums=sumstamples(u,e1,e2,e3,e4,mi)
                a=dsqrt(detlink(sums))
                u(e1,e2,e3,e4,mi,:)=linkmult(hbmatrix(a,beta),dlink(sums))/a
                u(e1,e2,e3,e4,mi,:)=u(e1,e2,e3,e4,mi,:)/dsqrt(detlink(u(e1,e2,e3,e4,mi,:)))
          !  !$omp end ordered
        enddo
        enddo
        enddo
        enddo
        !!$omp end parallel do
        !=============================================================================!
        !$omp section
        mi=4
        !!$omp parallel do ordered
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=2,nt,2
            !!$omp ordered
                sums=sumstamples(u,e1,e2,e3,e4,mi)
                a=dsqrt(detlink(sums))
                u(e1,e2,e3,e4,mi,:)=linkmult(hbmatrix(a,beta),dlink(sums))/a
                u(e1,e2,e3,e4,mi,:)=u(e1,e2,e3,e4,mi,:)/dsqrt(detlink(u(e1,e2,e3,e4,mi,:)))
            !!$omp end ordered
        enddo
        enddo
        enddo
        enddo
        !!$omp end parallel do

        !!$omp parallel do ordered
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nt-1,2
         !   !$omp ordered
                sums=sumstamples(u,e1,e2,e3,e4,mi)
                a=dsqrt(detlink(sums))
                u(e1,e2,e3,e4,mi,:)=linkmult(hbmatrix(a,beta),dlink(sums))/a
                u(e1,e2,e3,e4,mi,:)=u(e1,e2,e3,e4,mi,:)/dsqrt(detlink(u(e1,e2,e3,e4,mi,:)))
          !  !$omp end ordered
        enddo
        enddo
        enddo
        enddo
        !$omp end parallel sections
        !=============================================================================!
    endsubroutine hbstep

    subroutine init(u,beta,type)
        real(kind=r8), intent(inout) :: u(nr,nr,nr,nt,4,4)
        real(kind=r8) :: beta
        integer(kind=i4) :: e1, e2, e3, e4, mi, type
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nt
            do mi=1,4
                if(type==1)then         ! cold init
                    u(e1,e2,e3,e4,mi,1:3)=0.0_r8
                    u(e1,e2,e3,e4,mi,4)=1.0_r8
                elseif(type==2)then         ! hot init
                    u(e1,e2,e3,e4,mi,:)=hbmatrix(1.0_r8,beta)
                !elseif(type==3)then    ! continue from a previous configuration
                endif
            enddo
        enddo
        enddo
        enddo
        enddo
    endsubroutine init

    ! function that computes a square wilson loop on a site n
    function wilsonplanar(u,n,mi,ni,a,b)
        real(kind=r8) :: u(nr,nr,nr,nr,4,4), w(4), t(4), sum, wilsonplanar
        integer(kind=i4) :: mi, ni, n(4), a, b, i
        sum=0.0_r8
        
        w=ident()
        do i=1,a
            t=u(n(1),n(2),n(3),n(4),mi,:)
            w=linkmult(w,t)
            n(mi)=n(mi)+1
            if(n(mi).gt.nr) n(mi)=1
        enddo
                
        do i=1,b
            t=u(n(1),n(2),n(3),n(4),ni,:)
            w=linkmult(w,t)
            n(ni)=n(ni)+1
            if(n(ni).gt.nr) n(ni)=1
        enddo
        do i=1,a
            n(mi)=n(mi)-1
            if(n(mi).lt.1) n(mi)=nr
            t=u(n(1),n(2),n(3),n(4),mi,:)
            w=linkmult(w,dlink(t))
        enddo
        do i=1,b
            n(ni)=n(ni)-1
            if(n(ni).lt.1) n(ni)=nr
            t=u(n(1),n(2),n(3),n(4),ni,:)
            w=linkmult(w,dlink(t))
        enddo
        wilsonplanar=w(4)
    endfunction wilsonplanar

    ! function that makes the measure of the wilson loop in the lattice
    subroutine measurewilson(u,a,b,w)
        real(kind=r8), intent(in) :: u(nr,nr,nr,nt,4,4)
        real(kind=r8), intent(out) :: w(a,b)
        real(kind=r8) :: sum
        integer(kind=i4) :: a,b, i, j, n(4), e1, e2, e3, e4, mu, nu
        ! we compute the loops w(a,b) and save the result
        do i=1,a
        do j=1,b
            sum=0.0_r8
            !$omp parallel do ordered 
            ! compute the mean wilson loop
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nt
                n=(/e1,e2,e3,e4/)
                do mu=1,3
                do nu=mu+1,4
                    !$omp ordered
                    sum=sum+wilsonplanar(u,n,mu,nu,i,j)
                    ! the lattice is simmetric, so we
                    ! must consider the case a->b and b->a
                    ! now we interchange a and b
                    if(i.ne.j)then 
                    sum=sum+wilsonplanar(u,n,mu,nu,j,i)                    
                    endif
                    ! the loop for this site is calculated!
                    !$omp end ordered
                enddo
                enddo
            enddo
            enddo
            enddo
            enddo
            !$omp end parallel do
            if(i.eq.j)then
                w(i,j)=sum/(6.0_r8*nt*nr**3)    ! quadratic loops
            else
                w(i,j)=sum/(12_r8*nt*nr**3) ! rectangular loops
            endif
        enddo
        enddo
    endsubroutine measurewilson

    ! subroutine that does the center projection
    subroutine centerprojection(u,z)
        real(kind=r8) :: u(nr,nr,nr,nt,4,4), z(nr,nr,nr,nt,4,4)
        integer(kind=i4) :: e1,e2,e3,e4,mi!, z(nr,nr,nr,nt,4)
        !$omp parallel do
        do mi=1,4
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nr
            
                if(u(e1,e2,e3,e4,mi,4).gt.0.0_r8)then
                    z(e1,e2,e3,e4,mi,:)=(/0.0_r8,0.0_r8,0.0_r8,1.0_r8/)
                else
                    z(e1,e2,e3,e4,mi,:)=(/0.0_r8,0.0_r8,0.0_r8,-1.0_r8/)
                endif
            enddo
        enddo
        enddo
        enddo
        enddo
        !!$omp parallel enddo
    endsubroutine centerprojection

    ! subroutine that remove the vortices
    subroutine centerremotion(u,ur)
        real(kind=r8) :: u(nr,nr,nr,nt,4,4),ur(nr,nr,nr,nt,4,4)
        integer(kind=i4) :: e1,e2,e3,e4,mi
        !$omp parallel do
        do mi=1,4
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nr
            !do mi=1,4
                if(u(e1,e2,e3,e4,mi,4).gt.0.0_r8) ur(e1,e2,e3,e4,mi,:)=u(e1,e2,e3,e4,mi,:)
                if(u(e1,e2,e3,e4,mi,4).lt.0.0_r8) ur(e1,e2,e3,e4,mi,:)=-u(e1,e2,e3,e4,mi,:)
            enddo
        enddo
        enddo
        enddo
        enddo
        !$omp end parallel do
    endsubroutine centerremotion

    ! function that says if we have a vortex in a given plaquette
    function vortexinplaq(z,n,mi,ni)
        real(kind=r8) ::  vortexinplaq,z(nr,nr,nr,nt,4,4)
        integer(kind=i4) :: mi, ni
        integer(kind=i4), dimension(4) :: n,npz,nmz,nzp,nzm,npp,npm,nmp

        call neighborhood(n,npz,nmz,nzp,nzm,npp,npm,nmp,mi,ni)

        vortexinplaq=z(n(1),n(2),n(3),n(4),mi,4)*z(npz(1),npz(2),npz(3),npz(4),ni,4)*&
        &z(nzp(1),nzp(2),nzp(3),nzp(4),mi,4)*z(n(1),n(2),n(3),n(4),ni,4)
    endfunction vortexinplaq

    ! function that returns the number of vortices in a flat surface in the lattice
    function numberofvortices(z,n,mi,ni,lmi,lni)
        real(kind=r8) :: z(nr,nr,nr,nt,4,4), nv,  numberofvortices
        integer(kind=i4) ::  mi, ni, lmi, lni, ini, imi
        integer(kind=i4) :: n(4), m(4)
        nv=0
        m=n
        ! we need to check plaquete per plaquette of the lattice
        do imi=1,lmi
        ! move in the ni direction
        do ini=1,lni
            if(vortexinplaq(z,m,mi,ni).lt.0.0_r8) nv=nv+1.0_r8
            m(ni)=m(ni)+1
            if(m(ni).gt.nr) m(ni)=1
        enddo
        ! move in the mi direction for all the ni's
        m(mi)=m(mi)+1
        if(m(mi).gt.nr) m(mi)=1
        m(ni)=n(ni)
        enddo

        numberofvortices=nv
    endfunction numberofvortices

    ! subroutine that fixes the direct maximal center gauge
    subroutine maximal_center_gauge(u,ug,tol)
        real(kind=r8), dimension(nr,nr,nr,nr,4,4) :: u, ug
        real(kind=r8), dimension(nr,nr,nr,nr,4) :: g
        real(kind=r8) :: d(4,8), bigd(4,4), x(4,1), norm, tracenew, traceold, ngs, omega, tol
        integer(kind=i4) :: e1, e2, e3, e4, mi, n(4)
        integer(kind=i4) :: i, j, l, igs
        ! initialize the gauge fixed configuration
        ug=u
        traceold=0.0_r8
        g=1.0_r8
        omega=1.70_r8
        
        do igs=1,10000
            ! visit each lattice site
            !!$omp parallel do ordered
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
                !!$omp ordered
                ! we must pick the values of d(i,l)
                do l=1,4
                    n=(/e1,e2,e3,e4/)
                    d(1,l)=-ug(n(1),n(2),n(3),n(4),l,1)
                    d(2,l)=-ug(n(1),n(2),n(3),n(4),l,2)
                    d(3,l)=-ug(n(1),n(2),n(3),n(4),l,3)
                    d(4,l)=ug(n(1),n(2),n(3),n(4),l,4)
                    n(l)=n(l)-1
                    if(n(l).lt.1) n(l)=nr
                    d(1,l+4)=ug(n(1),n(2),n(3),n(4),l,1)
                    d(2,l+4)=ug(n(1),n(2),n(3),n(4),l,2)
                    d(3,l+4)=ug(n(1),n(2),n(3),n(4),l,3)
                    d(4,l+4)=ug(n(1),n(2),n(3),n(4),l,4)
                enddo

                ! now we can find the matrix bigd, tha must be diagonalized
                do i=1,4
                do j=1,4
                    bigd(i,j)=0.0_r8
                    do l=1,8
                    bigd(i,j)=bigd(i,j)+d(i,l)*d(j,l)
                    enddo
                enddo
                enddo

                ! now we must find the eigenvector associated with the bigest eigenvalue of this matrix
                ! we will use the power method
                x=1.0_r8
                do i=1,50
                    x=matmul(bigd,x)
                    ! we must find the bigger component of x
                    norm=x(1,1)
                    do j=2,4
                        if(x(j,1).gt.norm) norm=x(j,1)
                    enddo
                    x=x/norm
                enddo

                ! we found the gauge transformation
                g(e1,e2,e3,e4,1)=x(1,1)
                g(e1,e2,e3,e4,2)=x(2,1)
                g(e1,e2,e3,e4,3)=x(3,1)
                g(e1,e2,e3,e4,4)=x(4,1)

                ! overrelaxation step
                g(e1,e2,e3,e4,:)=g(e1,e2,e3,e4,:)-ident()
                g(e1,e2,e3,e4,:)=g(e1,e2,e3,e4,:)*omega
                g(e1,e2,e3,e4,:)=g(e1,e2,e3,e4,:)+ident()
                norm=dsqrt(detlink(g(e1,e2,e3,e4,:)))
 
                g(e1,e2,e3,e4,1)=g(e1,e2,e3,e4,1)/norm
                g(e1,e2,e3,e4,2)=g(e1,e2,e3,e4,2)/norm
                g(e1,e2,e3,e4,3)=g(e1,e2,e3,e4,3)/norm
                g(e1,e2,e3,e4,4)=g(e1,e2,e3,e4,4)/norm
                
                do mi=1,4
                n=(/e1,e2,e3,e4/)
                n(mi)=n(mi)-1
                if(n(mi).lt.1) n(mi)=nr
                ug(e1,e2,e3,e4,mi,:)=linkmult(g(e1,e2,e3,e4,:),ug(e1,e2,e3,e4,mi,:))
                ug(n(1),n(2),n(3),n(4),mi,:)=linkmult(ug(n(1),n(2),n(3),n(4),mi,:),dlink(g(e1,e2,e3,e4,:)))
                enddo
                !!$omp end ordered
            enddo
            enddo
            enddo
            enddo
            !!$omp end parallel do
            ! we gauge transform the configuration
            tracenew=0.0_r8
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
                do mi=1,4
                    !n=(/e1,e2,e3,e4/)
                    !n(mi)=n(mi)+1
                    !if(n(mi).gt.nr) n(mi)=1
                    !ug(e1,e2,e3,e4,mi,:)=linkmult(g(e1,e2,e3,e4,:),&
                    !&linkmult(ug(e1,e2,e3,e4,mi,:),dlink(g(n(1),n(2),n(3),n(4),:))))

                    ! compute the trace of the transformed configuration
                    tracenew=tracenew+ug(e1,e2,e3,e4,mi,4)**2.0_r8/4.0_r8
                enddo
            enddo
            enddo
            enddo
            enddo
            if(abs(tracenew-traceold).lt.tol)then
                print*,'Maximal Center Gauge fixed, tolerance=',abs(tracenew-traceold)
                exit
            endif
            traceold=tracenew
        enddo
        if(i==10000) print*,'MAG NOT fixed, tolerance=',abs(tracenew-traceold)
    endsubroutine maximal_center_gauge
    
    subroutine pvortexdensity(z,a,b)
        real(kind=r8) :: z(nr,nr,nr,nt,4,4), probeven, probodd, nv
        integer(kind=i4) :: i, j, a, b, lamba, area
        integer(kind=i4) :: e1, e2, e3, e4, mi, ni
        lamba=nr*nr*nr*nt
        do i=1,a
        do j=1,b
            probeven=0.0_r8
            probodd=0.0_r8
            area=i*j
            !$omp parallel do ordered shared(probeven,probodd)
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nt
                do mi=1,3
                do ni=mi+1,4
                    !$omp ordered
                    ! count the number of vortices
                    nv=numberofvortices(z,(/e1,e2,e3,e4/),mi,ni,i,j)
                    ! nv is odd
                    if(mod(nv,2.0_r8).gt.0.0_r8)then
                        probodd=probodd+1.0_r8
                    else
                        probeven=probeven+1.0_r8
                    endif
                    !$omp end ordered
                enddo
                enddo
            enddo
            enddo
            enddo
            enddo
            !$omp end parallel do
            ! we have the probability for this area, get out with it
            probeven=probeven/(6.0_r8*lamba)
            probodd=probodd/(6.0_r8*lamba)
            write(400+area,*) area, probeven, probodd, probeven+probodd
        enddo
        enddo
    endsubroutine pvortexdensity

    subroutine vortexsperatedwilson(u,z,a,b,we,wo)
        real(kind=r8), intent(in) :: u(nr,nr,nr,nt,4,4),z(nr,nr,nr,nt,4,4)
        real(kind=r8), intent(out) :: we(a,b), wo(a,b)
        real(kind=r8) :: sumo, sume
        integer(kind=i4) :: a,b, i, j, n(4), ne, no
        integer(kind=i4) :: e1, e2, e3, e4, mu, nu
        ! we compute the loops w(a,b) and save the result
        do i=1,a
        do j=1,b
            sumo=0.0_r8
            sume=0.0_r8
            ne=0
            no=0
            ! compute the mean wilson loop
            !$omp parallel do
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nt
                n=(/e1,e2,e3,e4/)
                do mu=1,3
                do nu=mu+1,4
                    if(mod(numberofvortices(z,n,mu,nu,i,j),2.0_r8).gt.0.5_r8)then
                        sumo=sumo+wilsonplanar(u,n,mu,nu,i,j)
                        no=no+1
                    else
                        sume=sume+wilsonplanar(u,n,mu,nu,i,j)
                        ne=ne+1
                    endif
                    ! the lattice is simmetric, so we
                    ! must consider the case a->b and b->a
                    ! now we interchange a and b
                    if(i.ne.j)then 
                        if(mod(numberofvortices(z,n,mu,nu,j,i),2.0_r8).gt.0.5_r8)then
                            sumo=sumo+wilsonplanar(u,n,mu,nu,j,i)
                            no=no+1
                        else
                            sume=sume+wilsonplanar(u,n,mu,nu,j,i)
                            ne=ne+1
                        endif
                    endif
                    ! the loop for this site is calculated!
                enddo
                enddo
            enddo
            enddo
            enddo
            enddo
            !$omp end parallel do
            we(i,j)=sume/ne
            wo(i,j)=sumo/no
        enddo
        enddo
    endsubroutine vortexsperatedwilson

    subroutine vortexcomponets(wf,we,wo,a,b)
        real(kind=r8), dimension(a,b) :: wf, wo, we
        integer(kind=i4) :: a, b, i, j, area
        do i=1,a
        do j=1,b
            area=i*j
            write(500+area,*) area, wf(i,j), we(i,j), wo(i,j)
        enddo
        enddo
    endsubroutine vortexcomponets
    
    subroutine vortexlimitedwilson(u,z,a,b)
        real(kind=r8), dimension(nr,nr,nr,nt,4,4) :: u, z
        real(kind=r8) :: w0, w1, w2
        integer(kind=i4) :: ia, ib, a, b, n0, n1, n2, nv
        integer(kind=i4) :: e1, e2, e3, e4, mi, ni
        do ia=1,a
        do ib=1,b
            w0=0.0_r8
            w1=0.0_r8
            w2=0.0_r8
            n0=0
            n1=0
            n2=0
            ! we sweep the lattice
            !$omp parallel do
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
                do mi=1,3
                do ni=mi+1,3
                    ! we select the data
                    nv=int(numberofvortices(z,(/e1,e2,e3,e4/),mi,ni,ia,ib))
                    if(nv.eq.0)then
                        w0=w0+wilsonplanar(u,(/e1,e2,e3,e4/),mi,ni,ia,ib)
                        n0=n0+1
                    elseif(nv.eq.1)then
                        w1=w1+wilsonplanar(u,(/e1,e2,e3,e4/),mi,ni,ia,ib)
                        n1=n1+1
                    elseif(nv.eq.2)then
                        w2=w2+wilsonplanar(u,(/e1,e2,e3,e4/),mi,ni,ia,ib)
                        n2=n2+1
                    endif
                enddo
                enddo
            enddo
            enddo
            enddo
            enddo
            !$omp end parallel do
            write(600+ia*ib,*) ia*ib, w0, w1, w2
        enddo
        enddo
    endsubroutine vortexlimitedwilson

    ! subroutine that measures all kind of wilson loops
    subroutine wilsonvortex(u,ur,z,w,wr,wv,we,wo,a,b)
        real(kind=r8), dimension(nr,nr,nr,nt,4,4), intent(in) :: u, ur, z
        real(kind=r8), dimension(a,b), intent(out) :: w, wr, wv, we, wo
        real(kind=r8) :: sum, sumr, sumv, sume, sumo, w0, w1, w2
        real(kind=r8) :: probeven, probodd, nv
        integer(kind=i4) :: a, b, i, j, ne, no, n(4)
        integer(kind=i4) :: e1, e2, e3, e4, mi, ni
        ! we compute the loops w(a,b) and save the result
        
        
        do i=1,a
        do j=1,b
            ! initialize the sums
            sum=0.0_r8
            sumr=0.0_r8
            sumv=0.0_r8
            sume=0.0_r8
            sumo=0.0_r8
            w0=0.0_r8
            w1=0.0_r8
            w2=0.0_r8
            probeven=0.0_r8
            probodd=0.0_r8
            no=0
            ne=0
            ! sweep the lattice
            ! compute the mean all the kinds of wilson loops
            !$OMP parallel do ordered
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nt
                n=(/e1,e2,e3,e4/)
                do mi=1,3
                do ni=mi+1,4
                    !$omp ordered
                    ! =============================================================== !
                    ! compute the normal wilson loops
                    ! full  wlson

                    sum=sum+wilsonplanar(u,n,mi,ni,i,j)
                    ! vortex projected wilson
                    sumv=sumv+wilsonplanar(z,n,mi,ni,i,j)
                    ! vortex removed wilson
                    sumr=sumr+wilsonplanar(ur,n,mi,ni,i,j)

                    ! the lattice is simmetric, so we
                    ! must consider the case a->b and b->a
                    ! now we interchange a and b
                    if(i.ne.j)then 
                        ! full  wlson
                        sum=sum+wilsonplanar(u,n,mi,ni,i,j)
                        ! vortex projected wilson
                        sumv=sumv+wilsonplanar(z,n,mi,ni,i,j)
                        ! vortex removed wilson
                        sumr=sumr+wilsonplanar(ur,n,mi,ni,i,j)
                    endif
                    !$omp end ordered  
                    ! =============================================================== !
                    ! =============================================================== !
                    ! compute the wilson loops with even and odd vortex numbers
                    if(mod(numberofvortices(z,n,mi,ni,i,j),2.0_r8).gt.0.5_r8)then
                        sumo=sumo+wilsonplanar(u,n,mi,ni,i,j)
                        no=no+1
                    else
                        sume=sume+wilsonplanar(u,n,mi,ni,i,j)
                        ne=ne+1
                    endif
                    ! the lattice is simmetric, so we
                    ! must consider the case a->b and b->a
                    ! now we interchange a and b
                    if(i.ne.j)then 
                        if(mod(numberofvortices(z,n,mi,ni,j,i),2.0_r8).gt.0.5_r8)then
                            sumo=sumo+wilsonplanar(u,n,mi,ni,j,i)
                            no=no+1
                        else
                            sume=sume+wilsonplanar(u,n,mi,ni,j,i)
                            ne=ne+1
                        endif
                    endif
                    ! =============================================================== !
                    ! =============================================================== !
                    ! compute the vortex limited wilson
                    ! we select the data
                    nv=int(numberofvortices(z,(/e1,e2,e3,e4/),mi,ni,i,j))
                    if(nv.lt.0.50_r8)then
                        w0=w0+wilsonplanar(u,(/e1,e2,e3,e4/),mi,ni,i,j)
                    elseif((nv.gt.0.5_r8).and.(nv.lt.1.5_r8))then
                        w1=w1+wilsonplanar(u,(/e1,e2,e3,e4/),mi,ni,i,j)
                    elseif((nv.lt.2.5_r8).and.(nv.gt.1.5_r8))then
                        w2=w2+wilsonplanar(u,(/e1,e2,e3,e4/),mi,ni,i,j)
                    endif
                    ! =============================================================== !
                    ! =============================================================== !
                    ! count the fraction of vortex pierced loops
                    ! count the number of vortices
                    
                    nv=numberofvortices(z,(/e1,e2,e3,e4/),mi,ni,i,j)
                    ! nv is odd
                    if(mod(nv,2.0_r8).gt.0.0_r8)then
                        probodd=probodd+1.0_r8
                    else
                        probeven=probeven+1.0_r8
                    endif
                enddo
                enddo
            enddo
            enddo
            enddo
            enddo
            !$omp end parallel do
            if(i.eq.j)then
                w(i,j)=sum/(6.0_r8*nt*nr**3)    ! quadratic loops
                wv(i,j)=sumv/(6.0_r8*nt*nr**3)    ! quadratic loops
                wr(i,j)=sumr/(6.0_r8*nt*nr**3)    ! quadratic loops
            else
                w(i,j)=sum/(12_r8*nt*nr**3) ! rectangular loops
                wv(i,j)=sumv/(12_r8*nt*nr**3) ! rectangular loops
                wr(i,j)=sumr/(12_r8*nt*nr**3) ! rectangular loops
            endif
            we(i,j)=sume/ne
            wo(i,j)=sumo/no
            probeven=probeven/(6.0_r8*nt*nr**3)
            probodd=probodd/(6.0_r8*nt*nr**3)
            write(400+i*j,*) i*j, probeven, probodd, probeven+probodd
            write(600+i*j,*) i*j, w0, w1, w2
        enddo
        enddo
        do i=1,a
        write(100+i,*) i, w(i,:)
        write(200+i,*) i, wv(i,:)
        write(300+i,*) i, wr(i,:)
        do j=1,b
            write(500+i*j,*) i*j, w(i,j), we(i,j), wo(i,j)
        enddo
        enddo
    endsubroutine wilsonvortex

    subroutine save_lattice(u,nr,iconf)
        real(kind=r8) :: u(nr,nr,nr,nr,4,4)
        integer(kind=i4) :: nr, iconf, e1, e2, e3, e4, mi, seed(1)
        character(len=60) :: conf, str, b

        ! create the file of the form
        ! conf_nre4_sweep=i.dat
        write(b,'(f4.2)') beta
        conf='configurations/beta='
        conf=trim(conf)//trim(b)
        conf=trim(conf)//'_lattice='
        if(nr.lt.10) write(str,'(i1)') nr
        if(nr.ge.10) write(str,'(i2)') nr

        conf=trim(conf)//trim(str)
        conf=trim(conf)//'e4_sweep='

        if(iconf.lt.10) write(str,'(i1)') iconf
        if((iconf.ge.10).and.(iconf.lt.100)) write(str,'(i2)') iconf
        if((iconf.ge.100).and.(iconf.lt.1000)) write(str,'(i3)') iconf
        if((iconf.ge.1000).and.(iconf.lt.10000)) write(str,'(i4)') iconf

        conf=trim(conf)//trim(str)
        conf=trim(conf)//'.dat'
        print*,conf
        open(unit=10,file=conf)

        ! write the configuration
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nr
            do mi=1,4
                write(10,*) e1, e2, e3, e4, mi, u(e1,e2,e3,e4,mi,:)
            enddo
        enddo
        enddo
        enddo
        enddo
        close(10)
    endsubroutine save_lattice

    subroutine read_lattice(u,nr,iconf)
        real(kind=r8), intent(inout) :: u(nr,nr,nr,nr,4,4)
        real(kind=r8) :: uaux(4)
        integer(kind=i4) :: nr, iconf, e1, e2, e3, e4, mi, i, seed(1)
        character(len=60) :: conf, str, b

        ! create the file of the form
        ! conf_nre4_sweep=i.dat
        write(b,'(f4.2)') beta
        conf='configurations/beta='
        conf=trim(conf)//trim(b)
        conf=trim(conf)//'_lattice='
        if(nr.lt.10) write(str,'(i1)') nr
        if(nr.ge.10) write(str,'(i2)') nr

        conf=trim(conf)//trim(str)
        conf=trim(conf)//'e4_sweep='

        if(iconf.lt.10) write(str,'(i1)') iconf
        if((iconf.ge.10).and.(iconf.lt.100)) write(str,'(i2)') iconf
        if((iconf.ge.100).and.(iconf.lt.1000)) write(str,'(i3)') iconf
        if((iconf.ge.1000).and.(iconf.lt.10000)) write(str,'(i4)') iconf

        conf=trim(conf)//trim(str)
        conf=trim(conf)//'.dat'
        print*,conf
        open(unit=10,file=conf)

        ! write the configuration
        do i=1,nr**3*nt*4
            read(10,*) e1,e2,e3,e4,mi,uaux(:)
            u(e1,e2,e3,e4,mi,:)=uaux
        enddo
        close(10)
    endsubroutine read_lattice

    ! subroutine that makes one step of smearing
    subroutine smearing(u,usmeared)
        real(kind=r8), dimension(nr,nr,nr,nr,4,4) :: u, usmeared
        real(kind=r8) :: sums(4), temp(4)
        integer(kind=i4) :: x, y, z, t, mi, ni
        integer(kind=i4), dimension(4) :: nzz, npz, nmz, nzp, nzm, npp, npm, nmp
        usmeared=0.0_r8
        ! run over the lattice
        do x=1,nr
        do y=1,nr
        do z=1,nr
        do t=1,nr
            do mi=1,4
                ! compute the sum of the stamples
                nzz=(/x,y,z,t/)
                sums=0.0_r8
                do ni=1,4
                    if (ni.ne.mi) then
                        call neighborhood(nzz,npz,nmz,nzp,nzm,npp,npm,nmp,mi,ni)
                        temp=ident()
                        temp=linkmult(temp,u(nzz(1),nzz(2),nzz(3),nzz(4),ni,:))
                        temp=linkmult(temp,u(nzp(1),nzp(2),nzp(3),nzp(4),mi,:))
                        temp=linkmult(temp,dlink(u(npz(1),npz(2),npz(3),npz(4),ni,:)))
                        
                        sums=sums+temp
                        temp=ident()
                        temp=linkmult(temp,dlink(u(nzm(1),nzm(2),nzm(3),nzm(4),ni,:)))
                        temp=linkmult(temp,u(nzm(1),nzm(2),nzm(3),nzm(4),mi,:))
                        temp=linkmult(temp,u(npm(1),npm(2),npm(3),npm(4),ni,:))
                        sums=sums+temp
                    endif
                enddo
                usmeared(x,y,z,t,mi,:)=ulink(x,y,z,t,mi,:)+sums
                usmeared(x,y,z,t,mi,:)=usmeared(x,y,z,t,mi,:)/dsqrt(detlink(usmeared(x,y,z,t,mi,:)))
            enddo
        enddo
        enddo
        enddo
        enddo
    endsubroutine smearing

    ! subroutine that fixes the landau gauge
    subroutine landau_gauge(u,ug,tol)
        real(kind=r8), dimension(nr,nr,nr,nr,4,4) :: u, ug
        real(kind=r8) :: epsilon, omega, divmi, g(nr,nr,nr,nr,4)
        real(kind=r8):: tol, r(4), h(4)
        integer(kind=i4) :: e1, e2, e3, e4, mi, int, a, n(4), np(4), nm(4)
        g=1.0_r8
        omega=1.70_r8
        ! loop over the interactions
        do int=1,10000
            ! find the matrix of the transformation
            ! loop over the lattice
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
                ! compute the h matrix
                n=(/e1,e2,e3,e4/)
                h=0.0_r8
                do mi=1,4
                    np=n
                    nm=n
                    np(mi)=np(mi)+1
                    nm(mi)=nm(mi)-1
                    if(np(mi).gt.nr) np(mi)=1
                    if(nm(mi).lt.1) nm(mi)=nr

                    h=h+linkmult(u(n(1),n(2),n(3),n(4),mi,:),dlink(g(np(1),np(2),np(3),np(4),:)))
                    h=h+linkmult(dlink(u(nm(1),nm(2),nm(3),nm(4),mi,:)),dlink(g(nm(1),nm(2),nm(3),nm(4),:)))
                enddo

                ! compute the matrix of the overrelaxation step
                !r=linkmult(dlink(h),dlink(g(n(1),n(2),n(3),n(4),:)))/dsqrt(detlink(h))
                !r=r+(1.0_r8-omega)*ident()
                !r=r/dsqrt(detlink(r))
                !g(e1,e2,e3,e4,:)=linkmult(r,g(e1,e2,e3,e4,:))
                !g(e1,e2,e3,e4,:)=g(e1,e2,e3,e4,:)/dsqrt(detlink(g(e1,e2,e3,e4,:)))
                ! los alamos
                g(e1,e2,e3,e4,:)=dlink(h)/dsqrt(detlink(h))

            enddo
            enddo
            enddo
            enddo

            ! make the transformation
            epsilon=0.0_r8
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
                n=(/e1,e2,e3,e4/)
                do mi=1,4
                    nm=n
                    nm(mi)=nm(mi)+1
                    if(nm(mi).gt.nr) nm(mi)=1

                    ug(e1,e2,e3,e4,mi,:)=linkmult(g(e1,e2,e3,e4,:),&
                    linkmult(u(e1,e2,e3,e4,mi,:),dlink(g(nm(1),nm(2),nm(3),nm(4),:))))

                    !ug(e1,e2,e3,e4,mi,:)=ug(e1,e2,e3,e4,mi,:)/detlink(ug(e1,e2,e3,e4,mi,:))
                enddo
            enddo
            enddo
            enddo
            enddo

            epsilon=0.0_r8
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
                n=(/e1,e2,e3,e4/)
                do a=1,3
                    divmi=0.0_r8
                    do mi=1,4
                        nm=n
                        nm(mi)=nm(mi)-1
                        if(nm(mi).lt.1) nm(mi)=nr
                        divmi=divmi+ug(e1,e2,e3,e4,mi,a)-ug(nm(1),nm(2),nm(3),nm(4),mi,a)
                    enddo
                    epsilon=epsilon+divmi**2.0_r8
                enddo
            enddo
            enddo
            enddo
            enddo
            epsilon=epsilon/(nr**4.0_r8)
            if(epsilon.lt.tol) exit
        enddo
        print*,'gauge fixed, tolerance=',epsilon
    endsubroutine landau_gauge
end program su2