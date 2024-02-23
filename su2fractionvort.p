set parametric
set xr[0:30]
set tr[0:30]
set grid lw 3 lc 'gray' lt 2
set key box lw 2 lc 'black'
set yl ' p-vortex fraction'
set xl ' Loop area'
set xtics 0,5,30
set mxtics
set title 'Fraction of loops pierced by p-vortives'
set terminal epslatex color colortext standalone size 7,6 font ',20'
set output 'pvortexfraction.tex'
plot 'pvortex-probability.dat' u 1:2 lw 5 lc 'red' ps 2 pt 12 t ' Even',\
     '' u 1:3 lw 5 lc 'blue' ps 2 pt 6 t ' Odd',\
     t,0.5 lw 15 lt '--' lc 'black' t ''
set output

set xr[10:]
set pointintervalbox 3
set yl ' $W(C)$'
set title ' $W(C)$ for full, even and odd number of vortices configurations'
set output 'wilsoncomp.tex'
plot 'wilson-vortexcomp.dat' u 1:2 lw 5 lc 'black' ps 2 pt 12 t ' Full',\
     '' u 1:3 lw 5 lc 'blue' ps 1 pt 23 t ' Even',\
     '' u 1:4 lw 5 lc 'red' ps 1 pt 34 t ' Odd',\
     t,0 lw 5 lc 'black' t ''
set output

set xr[5:30]
set yr[-2:2]
set tr[0:30]
set title ' Ratio between the vortex limited Wilson loops'
set yl '$W_n(C)/W_0(C)$'
set key left
set output 'vortexlim.tex
plot t,1 lw 10 lt '--' lc 'black' t '', t,-1 lw 10 lt '--' lc 'black' t '', t,0 lw 6 lt '-.' lc 'black' t '',\
     'voterxlimied-wilson.dat' u 1:($3/$2) lw 2 lc 'red' pt 12 ps 2 t ' $W_1(C)/W_0(C)$',\
     'voterxlimied-wilson.dat' u 1:($4/$2) lw 2 lc 'blue' pt 34 ps 1.5 t ' $W_2(C)/W_0(C)$'
set output
