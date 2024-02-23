# plot and fit for the potentials
set xl '$r/a$'
set yl '$a\,V(r/a)-a\,V_0$'
set key left
set key Left
set key box lw 3 lc 'black'
set grid lw 3 lc 'gray'
set dummy x, y
set parametric
set print 'su2pot-fit.dat'
set output 'potentials.tex'
set terminal epslatex color colortext standalone size 7,6 font ",20"
set title "$q\\bar{q}$ potential for $\\beta=2.30$ on a $16^4$ lattice" 
#set style data lines

# fiting the data
v(t)=a+b*t+c/t
vp(t)=d+e*t
vr(t)=f+g/t

fit v(t) 'potqq.dat' via a, b, c
V0=a
print 'The full potential'
print "V0=",a
print 'a^2 sigma=',b
print 'c=',c
print ''

fit vp(t) 'potqq-proj.dat' via d, e
print 'The center-projected potential'
print 'a^2 sigma=',e
print 'c=',d
print ''

fit vr(t) 'potqq-rem.dat' via f, g
print 'The center-removed potential'
print 'a^2 sigma=',f
print 'c=',g

# fit the sommer parameter
V0=0

# now we make the plots
set tr[0:8]
set yr[-0.5:]
set xr [0:8]
plot 'potqq.dat' u 1:($2-a):3 w yerr lw 5 lc 'red' ps 2 pt 6 title ' Full',\
     t,v(t)-a w l lw 5 lc 'red' lt '.-.' t '',\
     'potqq-rem.dat' u 1:($2-a):3 w yerr lw 3 lc 'black' ps 2 pt 34 title ' Removed',\
     t, vr(t)-a w l lw 5 lc 'black' lt '.-.' t '',\
     'potqq-proj.dat' u 1:($2-d):3 w yerr lw 5 lc 'blue' ps 2 pt 12 title ' Projected',\
     t,vp(t)-d w l lw 8 lc 'blue' lt '.-.' t ''
    
set output

# now we plot the creutz ratios
unset xr
unset yr
set xr [0:5]
set output 'creutz.tex'
set title 'Creutz ratios for full, projected and removed configurations'
set yl ' $\chi(l,l)$'
set xl ' $l$'
plot 'fullcreutz.dat' lw 4 lc 'black' pt 6 ps 3 t ' Full',\
     'projcreutz.dat' lw 4 lc 'red' pt 12 ps 3 t ' Projected',\
     'remcreutz.dat' lw 4 lc 'blue' pt 34 ps 2 t ' Removed'
set output