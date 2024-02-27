module su2lattice
implicit none
    ! define the precision
    integer, parameter :: r8=selected_real_kind(8,15)
    integer, parameter :: i4=selected_int_kind(8)

    ! define the parameters
    real(kind=r8), parameter :: pi=acos(-1.0_r8)
    !integer(kind=i4), parameter :: nr=10 ! number of sites in space
    !integer(kind=i4), parameter :: nt=10 !  number of sites in time

    ! define the matrices
    real(kind=r8), allocatable, dimension(:,:,:,:,:,:) :: ulink, ugauge, ur, z
    real(kind=r8), allocatable, dimension(:,:) :: w, wr, wv
    
    ! defining our variables
    real(kind=r8) :: beta, start, end
    integer(kind=i4) :: nr, nt, e1, e2, e3, e4, mi, ni
    !integer(kind=i4) :: a, b
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
        real(kind=r8), intent(in), dimension(4) :: a, b
        real(kind=r8), dimension(4) :: linkmult
        real(kind=r8) :: ab(3), c(4)
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
        real(kind=r8), intent(in) :: u(4)
        real(kind=r8) :: detlink
        detlink=u(1)*u(1)+u(2)*u(2)+u(3)*u(3)+u(4)*u(4)
    endfunction detlink

    ! function that returns the dagger of the link
    function dlink(u)
        real(kind=r8), dimension(4) :: dlink
        real(kind=r8), intent(in), dimension(4) :: u!, v
        !v(1)=-u(1)
        !v(2)=-u(2)
        !v(3)=-u(3)
        !v(4)=u(4)
        dlink=(/-u(1),-u(2),-u(3),u(4)/)
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
endmodule su2lattice