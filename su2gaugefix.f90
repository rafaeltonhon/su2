module su2gaugefix
    use su2lattice
    implicit none
    contains
    ! subroutine that makes a random gauge transform
    subroutine random_gauge(u,ugauge)
        real(kind=r8), dimension(nr,nr,nr,nr,4,4) :: u, ugauge
        real(kind=r8), dimension(nr,nr,nr,nr,4) :: omega
        real(kind=r8) :: r(4)
        integer(kind=i4) :: e1, e2, e3, e4, mi, n(4)
        ! create the random matrix
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nr
            call random_number(r)
            omega(e1,e2,e3,e4,:)=r/dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3)+r(4)*r(4))
        enddo
        enddo
        enddo
        enddo

        ! do the gauge transformation
        do e1=1,nr
        do e2=1,nr
        do e3=1,nr
        do e4=1,nr
            do mi=1,4
                n=(/e1,e2,e3,e4/)
                n(mi)=n(mi)+1
                if(n(mi).gr.nr) n(mi)=1
                ugauge(e1,e2,e3,e4,mi,:)=linkmult(omega(e1,e2,e3,e4,:),&
                linkmult(u(e1,e2,e3,e4,mi,:),dlink(omega(n(1),n(2),n(3),n(4),:))))
            enddo
        enddo
        enddo
        enddo
        enddo
    endsubroutine random_gauge

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
module su2gaugefix