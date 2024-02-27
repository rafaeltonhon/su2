module su2measures
    use su2lattice
    implicit none
    contains
    function wilsonplanar(u,n,mi,ni,a,b)
        real(kind=r8) :: u(nr,nr,nr,nr,4,4), w(4), sum, wilsonplanar
        integer(kind=i4) :: mi, ni, n(4), a, b, i
        sum=0.0_r8
        w=ident()
        do i=1,a
            w=linkmult(w,u(n(1),n(2),n(3),n(4),mi,:))
            !if((a.eq.nr).and.(b.eq.nr)) print*,u(n(1),n(2),n(3),n(4),mi,:)
            !if((a.eq.nr).and.(b.eq.nr)) print*,n
            !if((a.eq.nr).and.(b.eq.nr)) print*,w
            n(mi)=n(mi)+1
            if(n(mi).gt.nr) n(mi)=1
        enddo
                
        do i=1,b
            w=linkmult(w,u(n(1),n(2),n(3),n(4),ni,:))
            !if((a.eq.nr).and.(b.eq.nr)) print*,u(n(1),n(2),n(3),n(4),ni,:)
            !if((a.eq.nr).and.(b.eq.nr)) print*,n
            !if((a.eq.nr).and.(b.eq.nr)) print*,w
            n(ni)=n(ni)+1
            if(n(ni).gt.nr) n(ni)=1
        enddo
        do i=1,a
            n(mi)=n(mi)-1
            if(n(mi).lt.1) n(mi)=nr
            w=linkmult(w,dlink(u(n(1),n(2),n(3),n(4),mi,:)))
            !if((a.eq.nr).and.(b.eq.nr)) print*,u(n(1),n(2),n(3),n(4),mi,:)
            !if((a.eq.nr).and.(b.eq.nr)) print*,n
            !if((a.eq.nr).and.(b.eq.nr)) print*,w
        enddo
        do i=1,b
            n(ni)=n(ni)-1
            if(n(ni).lt.1) n(ni)=nr
            w=linkmult(w,dlink(u(n(1),n(2),n(3),n(4),ni,:)))
            !if((a.eq.nr).and.(b.eq.nr)) print*,u(n(1),n(2),n(3),n(4),ni,:)
            !if((a.eq.nr).and.(b.eq.nr)) print*,n
            !if((a.eq.nr).and.(b.eq.nr)) print*,w
        enddo
        !print*,''
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
    
    
endmodule su2measures