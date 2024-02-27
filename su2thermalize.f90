module su2thermalize
    use su2lattice
    implicit none
    contains
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
        do e4=1,nr
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

endmodule su2thermalize