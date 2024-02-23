program su2
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
    real(kind=r8) :: beta, start, end
    integer(kind=i4) :: nr, nt
    integer(kind=i4) :: i, j, nmc, nterm, ncorr, iprint, a, b

    call cpu_time(start)

    open(unit=1,file='su2-in.dat')
    read(1,*) nr, nt            ! lattice size space time
    read(1,*) beta              ! beta
    read(1,*) nmc, nterm        ! number of monte carlo sweeps number of thermalization sweeps
    read(1,*) a, b              ! maximum size of the wilson loops we want
    read(1,*) iprint            ! if 1 we print the configurations on the file suN-lattice
    read(1,*) ncorr             ! number of correlation
    close(1)

    allocate(ulink(nr,nr,nr,nt,4,4),ur(nr,nr,nr,nt,4,4))
    allocate(ugauge(nr,nr,nr,nt,4,4),z(nr,nr,nr,nt,4,4))
    allocate(w(a,b),wv(a,b),wr(a,b),we(a,b),wo(a,b))
    call init(ulink,beta,1)
    
    ! termilize the lattice
    do i=1,nterm
        call hbstep(ulink,beta)
    enddo
    
    do i=1,nmc
        ! number of correlated configurations
        do j=1,ncorr
            ! discart ncorr configurations
            call hbstep(ulink,beta)
        enddo
        call measurewilson(ulink,a,b,w)

        ! project the lattice and measure the links
        call maximal_center_gauge(ulink,ugauge,0.0001_r8)
        call centerprojection(ugauge,z)
        call centerremotion(ugauge,ur)
        call wilsonvortex(ugauge,ur,z,w,wr,wv,we,wo,a,b)
    enddo

    deallocate(ulink,ugauge,z,wr,w,wv)
    call cpu_time(end)
    print*,"cpu time [seconds]: ", (end-start)/60
    write(9000,*) "cpu time [minutes]: ", (end-start)/60
    write(9000,*) "cpu time [hours]: ", (end-start)/3600
    
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
    ! ONLY FOR SQUARE LATTICES
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
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nt
            do mi=1,4
                sums=sumstamples(u,e1,e2,e3,e4,mi)
                a=dsqrt(detlink(sums))
                u(e1,e2,e3,e4,mi,:)=linkmult(hbmatrix(a,beta),dlink(sums))/a
                u(e1,e2,e3,e4,mi,:)=u(e1,e2,e3,e4,mi,:)/dsqrt(detlink(u(e1,e2,e3,e4,mi,:)))
            enddo
        enddo
        enddo
        enddo
        enddo
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
        real(kind=r8) :: u(nr,nr,nr,nr,4,4), w(4), sum, wilsonplanar
        integer(kind=i4) :: mi, ni, n(4), a, b, i
        sum=0.0_r8
        
        w=ident()
        do i=1,a
            w=linkmult(w,u(n(1),n(2),n(3),n(4),mi,:))
            n(mi)=n(mi)+1
            if(n(mi).gt.nr) n(mi)=1
        enddo
                
        do i=1,b
            w=linkmult(w,u(n(1),n(2),n(3),n(4),ni,:))
            n(ni)=n(ni)+1
            if(n(ni).gt.nr) n(ni)=1
        enddo
        do i=1,a
            n(mi)=n(mi)-1
            if(n(mi).lt.1) n(mi)=nr
            w=linkmult(w,dlink(u(n(1),n(2),n(3),n(4),mi,:)))
        enddo
        do i=1,b
            n(ni)=n(ni)-1
            if(n(ni).lt.1) n(ni)=nr
            w=linkmult(w,dlink(u(n(1),n(2),n(3),n(4),ni,:)))
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
            ! compute the mean wilson loop
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nt
                n=(/e1,e2,e3,e4/)
                do mu=1,3
                do nu=mu+1,4
                    sum=sum+wilsonplanar(u,n,mu,nu,i,j)
                    ! the lattice is simmetric, so we
                    ! must consider the case a->b and b->a
                    ! now we interchange a and b
                    if(i.ne.j)then 
                    sum=sum+wilsonplanar(u,n,mu,nu,j,i)                    
                    endif
                    ! the loop for this site is calculated!
                enddo
                enddo
            enddo
            enddo
            enddo
            enddo
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
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nr
            do mi=1,4
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
    endsubroutine centerprojection

    ! subroutine that remove the vortices
    subroutine centerremotion(u,ur)
        real(kind=r8) :: u(nr,nr,nr,nt,4,4),ur(nr,nr,nr,nt,4,4)
        integer(kind=i4) :: e1,e2,e3,e4,mi
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nr
            do mi=1,4
                if(u(e1,e2,e3,e4,mi,4).gt.0.0_r8) ur(e1,e2,e3,e4,mi,:)=u(e1,e2,e3,e4,mi,:)
                if(u(e1,e2,e3,e4,mi,4).lt.0.0_r8) ur(e1,e2,e3,e4,mi,:)=-u(e1,e2,e3,e4,mi,:)
            enddo
        enddo
        enddo
        enddo
        enddo
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
    subroutine maximal_center_gauge(u,ug,ngs)
        real(kind=r8), dimension(nr,nr,nr,nr,4,4) :: u, ug
        real(kind=r8), dimension(nr,nr,nr,nr,4) :: g
        real(kind=r8) :: d(4,8), bigd(4,4), x(4,1), norm, tracenew, traceold, ngs, omega
        integer(kind=i4) :: e1, e2, e3, e4, mi, n(4)
        integer(kind=i4) :: i, j, l, igs
        ! initialize the gauge fixed configuration
        ug=u
        traceold=0.0_r8
        g=1.0_r8
        omega=1.70_r8
        do igs=1,5000
            ! visit each lattice site
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nr
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
            enddo
            enddo
            enddo
            enddo

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
            if(abs(tracenew-traceold).lt.ngs) exit
            traceold=tracenew
        enddo
        !print*,abs(tracenew-traceold)
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
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nt
                do mi=1,3
                do ni=mi+1,4
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
            do e1=1,nr
            do e2=1,nr
            do e3=1,nr
            do e4=1,nt
                n=(/e1,e2,e3,e4/)
                do mi=1,3
                do ni=mi+1,4
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

        ! get the data out
        do i=1,a
         write(100+i,*) i, w(i,:)
        write(200+i,*) i, wv(i,:)
        write(300+i,*) i, wr(i,:)
        do j=1,b
            write(500+i*j,*) i*j, w(i,j), we(i,j), wo(i,j)
        enddo
        enddo
    endsubroutine wilsonvortex
end program su2