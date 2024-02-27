program su2main
    use su2lattice
    use su2thermalize
    use su2measures
    use su2center
    implicit none
    integer(kind=i4) :: i, j, nmc, nterm, ncorr, iprint, a, b

    call cpu_time(start)

    open(unit=1,file='su2-in.dat')
    read(1,*) nr, nt            ! lattice size space time
    read(1,*) beta              ! beta
    read(1,*) nmc, nterm        ! number of monte carlo sweeps number of thermalization sweeps
    read(1,*) a, b              ! maximum size of the wilson loops we want
    read(1,*) iprint            ! if 1 we print the configurations on the file suN-lattice
    close(1)

    allocate(ulink(nr,nr,nr,nt,4,4),ur(nr,nr,nr,nt,4,4))
    allocate(ugauge(nr,nr,nr,nt,4,4),z(nr,nr,nr,nt,4,4))
    allocate(w(a,b),wv(a,b),wr(a,b),we(a,b),wo(a,b))
    call init(ulink,beta,1)
    
    ! termilize the lattice
    do i=1,nterm
        call hbstep(ulink,beta)
    enddo
    ncorr=1
    do i=1,nmc
        ! number of correlated configurations
            call hbstep(ulink,beta)
        call measurewilson(ulink,a,b,w)
        ! project the lattice and measure the links
        !call maximal_center_gauge(ulink,ugauge,0.001_r8)
        !call centerprojection(ugauge,z)

        ! remove the vortices
        !call centerremotion(ur,ugauge)
        
        ! make the measures
        !call wilsonvortex(ulink,ur,z,w,wr,wv,we,wo,a,b)
    enddo

    deallocate(ulink,ugauge,z,wr,w,wv)
    call cpu_time(end)
    print*,"cpu time [seconds]: ", (end-start)/60
    write(9000,*) "cpu time [minutes]: ", (end-start)/60
    write(9000,*) "cpu time [hours]: ", (end-start)/3600
endprogram su2main