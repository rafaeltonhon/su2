import numpy as np
import matplotlib.pyplot as plt
import os
#====================================================================#
# FUNCTIONS FOR THE FITS
# function that fits to a curve like a*x1+b*x2+c
def fit_3parameters(x,y,sigmay,x1,x2):
    sumx1=sumx2=sumy=sumsigma2=sumx1x2=sumx1y=sumx2y=sumx12=sumx22=0
    n=len(x)
    if(sigmay==-1):
        sigmay=list()
        for i in range(0,n+1):
            sigmay.append(1)

            
    for i in range(0,n):
        sumx1+=x1(x[i])/(sigmay[i])**2
        sumx2+=x2(x[i])/(sigmay[i])**2
        sumy+=y[i]/(sigmay[i])**2
        
        sumx1x2+=x1(x[i])*x2(x[i])/(sigmay[i])**2
        sumx1y+=x1(x[i])*y[i]/(sigmay[i])**2
        sumx2y+=x2(x[i])*y[i]/(sigmay[i])**2

        sumx12+=(x1(x[i]))**2/(sigmay[i])**2
        sumx22+=(x2(x[i]))**2/(sigmay[i])**2
        sumsigma2+=1/(sigmay[i])**2

    summatrix=[[sumx12,sumx1x2,sumx1],[sumx1x2,sumx22,sumx2],[sumx1,sumx2,sumsigma2]]
    coeffmatrix=[[sumx1y],[sumx2y],[sumy]]

    summatrix=np.array(summatrix)
    coeffmatrix=np.array(coeffmatrix)

    inverse=np.linalg.inv(summatrix)
    coeffmatrix=np.matmul(inverse,coeffmatrix)
    return coeffmatrix[0][0],coeffmatrix[1][0],coeffmatrix[2][0]


# function that fits to a 2 parameter curve a*x1+b
def fit_2parameters(x,y,sigmay,x1):
    sumx1=sumy=sumx1y=sumx12=sumsigma2=0
    n=len(x)
    if(sigmay==-1):
        sigmay=list()
        for i in range(0,n):
            sigmay.append(1)
    
    for i in range(0,n):
        sumx1+=x1(x[i])/(sigmay[i])**2
        sumy+=y[i]/(sigmay[i])**2
        sumx1y+=x1(x[i])*y[i]/(sigmay[i])**2
        sumx12+=(x1(x[i]))**2/(sigmay[i])**2
        sumsigma2+=1/(sigmay[i])**2

    summatrix=[[sumx12,sumx1],[sumx1,sumsigma2]]
    coeffmatrix=[[sumx1y],[sumy]]

    summatrix=np.array(summatrix)
    coeffmatrix=np.array(coeffmatrix)

    summatrix=np.linalg.inv(summatrix)
    coeffmatrix=np.matmul(summatrix,coeffmatrix)

    return coeffmatrix[0][0],coeffmatrix[1][0]


# x to the minus 1 power f(x)=1/x
def xminus1(x):
    return 1/x


# x to the fist power f(x)=x
def xfirst(x):
    return x


#=============================================================+++++++#
#   STATISTICS FUNCTIONS
#=============================================================+++++++#
def mean(data):
    n=len(data)
    databar=0
    for i in range(0,n):
        databar+=data[i]
    return databar/n


# varianve of the random variable
def varianceX(data):
    n=len(data)
    databar=mean(data)
    sigma=0
    # make the mean value
    for i in range(0,n):
        sigma+=(data[i]-databar)**2
    return sigma/(n)


# variance of the mean estimator
def varianceXhat(data):
    n=len(data)
    databar=mean(data)
    sigma=0
    # make the mean value
    for i in range(0,n):
        sigma+=(data[i]-databar)**2
    return sigma/(n*(n-1))


# function that given a correlation number, make the bining of a set of data
def bining(filename):
    file=open(filename,'r')
    data=list()
    sigmak=list()
    lenk=list()
    # read the data from the file
    for line in file.readlines():
        l=[float(x) for x in line.split()]
        data.append(l[1])
    file.close()

    # number of correlated points
    nc=len(data)

    # now we construct the blocks
    for k in range(2,nc):
        blockeddata=[]
        # construct the blocks
        block=0
        nk=0
        for i in range(0,nc-nc%k):
            block+=data[i]
            nk+=1
            if(nk==k): # we have a block
                blockeddata.append(block/k)
                block=0
                nk=0
        # we have the data on blocks, that the standart deviation
        if(len(blockeddata)!=1): 
            sigmak.append(varianceXhat(blockeddata))
            lenk.append(k)
    
    # now we try to fit a a curve a/k+b to this set
    # the value kcorr for with the data is good, is our nof correlation
    file=open('bining.dat','w')
    for ns in range(0,len(lenk)-1):
        # fit the data
        x=lenk[ns::]
        y=sigmak[ns::]
        a,b=fit_2parameters(x,y,-1,xminus1)

        # compute the variance of the data with respect to the fit
        sum=0
        for i in range(0,len(x)//2): # we consider only half of the points to not consider blocks to large
            sum+=(y[i]-a/x[i]-b)**2
        file.write(' '.join([str(x) for x in [x[0],sum,]]))
        file.write('\n')
        plt.scatter(ns,sum)
        #plt.plot(x,a/np.array(x)+b)
        
    file.close()
    #plt.scatter(lenk,sigmak)
    plt.show()


# function that given a set of uncorrelated date, compute the error by the jackkinfe
def jackkinfe(data):
    thetahat=mean(data)
    thetan=list()
    n=len(data)
    sthetatilde=0
    # construct the sets if the n-th element removed
    for i in range(0,n):
        l=[]
        l=data[::]
        del l[i]
        thetan.append(mean(l))

    for i in range(1,n):
        sthetatilde+=(thetan[i]-thetahat)**2
    return (n-1)*sthetatilde/n


# jackknife for a plot data
def jackkinfe_ploteddata(datax,datay,np,x1,x2):
    thetan1=list()
    thetan2=list()
    n=len(datax)
    sthetatilde1=0
    sthetatilde2=0

    if np==2: # 2 parameter fit
        thetahat1,thetahat2=fit_2parameters(datax,datay,-1,x1)
        
    if np==3: # 2 parameter fit
        thetahat1,thetahat2, thetahat3=fit_3parameters(datax,datay,-1,x1,x2)
        thetan3=list()
        sthetatilde3=0

    # construct the sets if the n-th element removed
    for i in range(0,n):
        x=[]
        y=[]
        x=datax[::]
        y=datay[::]
        del x[i]
        del y[i]
        
        if np==2:
            a,b=fit_2parameters(x,y,-1,x1)
            thetan1.append(a)
            thetan2.append(b)
        if np==3:
            a,b,c=fit_3parameters(x,y,-1,x1,x2)
            thetan1.append(a)
            thetan2.append(b)
            thetan3.append(c)

    # compute the variations
    for i in range(1,n):
        sthetatilde1+=(thetan1[i]-thetahat1)**2
        sthetatilde2+=(thetan2[i]-thetahat2)**2
        if(np==3): sthetatilde3+=(thetan3[i]-thetahat3)**2

    sthetatilde1*=(n-1)/n
    sthetatilde2*=(n-1)/n
    if(np==3): sthetatilde3*=(n-1)/n

    if(np==2): return sthetatilde1, sthetatilde2
    if(np==3): return sthetatilde1, sthetatilde2, sthetatilde3


#=============================================================+++++++#
# function that, given a number of correlation, put the data of a file in blocks
# the output is a file with name unc'name of the file'
# with the blocked, uncorrelated data
def block_data(filein,ncorr):
    data=list()
    file=open(filein,'r')
    # read the data from a file
    n=0
    for line in file.readlines():
        l=[float(x) for x in line.split()]

        # we will put the data of the same column in the same array
        if (n==0):
            for i in range(1,len(l)):
                # we take the element from that column and append on data
                col=[]
                col.append(l[i])
                data.append(col)
        else:
            # we append the data in the respect column of data array
            for i in range(0,len(l)-1):
                data[i].append(l[i+1])
        n+=1
    file.close()
    
    # we have read the data, now we make blocks with if
    databloked=list()
    for i in range(0,len(data)): # we block each data column
        block=0
        blockcol=[]
        nblock=0
        for j in range(0,n-n%ncorr): # we block this column
            block+=data[i][j]
            nblock+=1
            if(nblock==ncorr):
                blockcol.append(block/ncorr)
                block=0
                nblock=0
        databloked.append(blockcol)
    data=[]
    # now we make the mean and the standar deviation of each new set of uncorrelated points
    for i in range(0,len(databloked)):
        data.append([mean(databloked[i]),np.sqrt(jackkinfe(databloked[i]))])
    
    
    # we get out with the uncorrelated means
    file1=open("unc"+filein,'w')
    file2=open("mean"+filein,'w')
    for i in range(0,len(databloked)):
        file1.write(' '.join([str(x) for x in databloked[i]]))
        file1.write('\n')

    for i in range(0,len(data)):
        file2.write(' '.join([str(x) for x in data[i]]))

        file2.write('\n')
    file1.close()
    file2.close()


# function that block the data for the wilson loops
def blockwilson(a,p,ncorr):
    # make the mean value of the wilson loops
    for i in range(1,a+1):
        # open the file
        if(p=='full'): file=f"fort.{100+i}" # full loops
        if(p=='proj'): file=f"fort.{200+i}" # vortex projected loops  
        if(p=='rem'): file=f"fort.{300+i}" # vortex remov loops
        if(p=='even'): file=f"fort.{500+i}" # vortex remov loops
        if(p=='odd'): file=f"fort.{600+i}" # vortex remov loops
        block_data(file,ncorr)


# function that computes the potential
def pot(a,b,p,nt):    
    w=list()
    sw=list()
    v=list()
    sv=list()
    t=list()
    r=list()
    for i in range(0,b):
        t.append(i+1)

    if(p=='full'): fileout=open(f"potqq.dat",'w') # full loops
    if(p=='proj'): fileout=open(f"potqq-proj.dat",'w') # vortex projected loops  
    if(p=='rem'): fileout=open(f"potqq-rem.dat",'w') # vortex remov loops

    # read data from file meanfort.10+a
    for i in range(1,a+1):
        w=[]
        sw=[]
        r.append(i)

        # open the file
        if(p=='full'): file=open(f"meanfort.{100+i}",'r') # full loops
        if(p=='proj'): file=open(f"meanfort.{200+i}",'r') # vortex projected loops  
        if(p=='rem'): file=open(f"meanfort.{300+i}",'r') # vortex remov loops

        # read the expetation values for the wilson loop of size rxt, with t=1,...b]]
        for line in file.readlines():
            l=[]
            l=[float(x) for x in line.split()]
            w.append(np.log(abs(l[0])))
            sw.append(abs(l[1])/l[0])
        file.close()
        # having the values we must make a linear regression to comput the potential
        # fit w=aexp(-v t)=:log(w)=log(a)-v*t
        av, bv =fit_2parameters(t[nt::],w[nt::],sw[nt::],xfirst)
        da, db=jackkinfe_ploteddata(t[nt::],w[nt::],2,xfirst,xfirst)
        v.append(-av)
        sv.append(np.sqrt(da))
    
        # get out with the data
        fileout.write(' '.join([str(x) for x in [i,-av,da]]))
        fileout.write('\n')
    plt.errorbar(r,v,yerr=sv,label=p,fmt='D')

    fileout.close()
    if p=='full':
        # fit the curve v(x)=sigma*r+k/r+v0 to the potential
        sigma, k, v0=fit_3parameters(r,v,sv,xfirst,xminus1)
        dsigma, dk, dv0=jackkinfe_ploteddata(r,v,3,xfirst,xminus1)
        return sigma, dsigma, k, dk, v0, dv0
    if p=='proj':
        # fit the curve v(x)=sigma*r+k/r+v0 to the potential
        sigma,  v0=fit_2parameters(r,v,sv,xfirst)
        dsigma, dv0=jackkinfe_ploteddata(r,v,2,xfirst,xminus1)
        return sigma, dsigma, v0, dv0
    if p=='rem':
        # fit the curve v(x)=sigma*r+k/r+v0 to the potential
        sigma,  v0=fit_2parameters(r,v,sv,xminus1)
        dsigma, dv0=jackkinfe_ploteddata(r,v,2,xminus1,-1)
        return sigma, dsigma, v0, dv0


# function that computes the creutz rations
def creutz(a,p):
    chi=list()
    dchi=list()
    for i in range(1,a+1):
        w=[]
        sw=[]
        # open the file
        if(p=='full'): file=open(f"meanfort.{100+i}",'r') # full loops
        if(p=='proj'): file=open(f"meanfort.{200+i}",'r') # vortex projected loops  
        if(p=='rem'): file=open(f"meanfort.{300+i}",'r') # vortex remov loops

        # read the data
        wi=[]
        swi=[]
        for line in file.readlines():
            l=[float(x) for x in line.split()]
            wi.append(l[0])
            swi.append(l[1])
        w.append(wi)
        sw.append(swi)
        file.close()
        
    # once we have the data, we compute the rations
    for i in range(0,a):
        if i==0:
            chi.append(-np.log(2*w[i][i]))
        else:
            print(i,w)
            print(w[i][i])
            print(w[i-1][i-1])
            print(w[i-1][i])
            print(w[i][i-1])
            chi.append(-np.log((w[i][i]*w[i-1][i-1])/(w[i-1][i]*w[i][i-1])))

    # get out with the data
    fileout=open(p+'creutz.dat','w')
    for i in range(0,a):
        fileout.write(' '.join([str(x) for x in [i+1,chi[i],'\n']]))
    fileout.close()


def wilsoncomponens(a,b):
    even=list()
    odd=list()
    full=list()
    area=list()
    for i in range(1,a+1):
        for j in range(1,b+1):
            file=open(f'fort.{500+i*j}','r')
            eveni=[]
            oddi=[]
            fulli=[]
            for line in file.readlines():
                l=[float(x) for x in line.split()]
                fulli.append(l[1])
                eveni.append(l[2])
                oddi.append(l[3])
            file.close()
            full.append(np.mean(fulli))
            even.append(np.mean(eveni))
            odd.append(np.mean(oddi))
            area.append(i*j)

    file=open('wilson-vortexcomp.dat','w')
    for i in range(0,len(area)):
        file.write(' '.join([str(x) for x in [area[i],full[i],even[i],odd[i],'\n']]))
    file.close()

# function that sepate
# function that deals with the projected-vortices (p-vortices) density
def pvorticesprob(a,b):
    even=list()
    odd=list()
    area=list()
    for i in range(1,a+1):
        for j in range(1,b+1):
            file=open(f'fort.{400+i*j}','r')
            eveni=[]
            oddi=[]
            for line in file.readlines():
                l=[float(x) for x in line.split()]
                eveni.append(l[1])
                oddi.append(l[2])
            file.close()
            even.append(np.mean(eveni))
            odd.append(np.mean(oddi))
            area.append(i*j)

    file=open('pvortex-probability.dat','w')
    for i in range(0,len(area)):
        file.write(' '.join([str(x) for x in [area[i],even[i],odd[i],'\n']]))
    file.close()


def vortexlimitedwilson(a,b):
    w0=list()
    w1=list()
    w2=list()
    area=list()
    for i in range(1,a+1):
        for j in range(1,b+1):
            file=open(f'fort.{600+i*j}','r')
            w0i=[]
            w1i=[]
            w2i=[]
            for line in file.readlines():
                l=[float(x) for x in line.split()]
                w0i.append(l[1])
                w1i.append(l[2])
                w2i.append(l[3])
            file.close()
            w0.append(np.mean(w0i))
            w1.append(np.mean(w1i))
            w2.append(np.mean(w2i))
            area.append(i*j)

    file=open('voterxlimied-wilson.dat','w')
    for i in range(0,len(area)):
        file.write(' '.join([str(x) for x in [area[i],w0[i],w1[i],w2[i],'\n']]))
    file.close()
#=============================================================================================================#
#=============================================================================================================#
#=============================================================================================================#
#               MAIN PROGRAM
#=============================================================================================================#
#=============================================================================================================#
#=============================================================================================================#
nc=2
nt=1
a=5
b=5

pvorticesprob(5,5)
wilsoncomponens(5,5)
vortexlimitedwilson(5,5)
blockwilson(a,'full',nc)
sigma, dsigma, k,dk,v0,dv0=pot(a,b,'full',nt)
print(f'Full configurations string tension \nsigma a2 = {sigma}+-{dsigma}\n')

blockwilson(a,'proj',nc)
sigma, dsigma,v0,dv0=pot(a,b,'proj',nt)
print(f'Vortex projected configurations string tension \nsigma a2 = {sigma}+-{dsigma}\n')

blockwilson(a,'rem',nc)
k,dk,v0,dv0=pot(a,b,'rem',nt)
#print(f'Vortex removed configurations string tension \nsigma a2 = {sigma}+-{dsigma}')

plt.legend()
#plt.show()

os.system('gnuplot su2pot.p')
os.system('gnuplot su2fraction.p')
os.system('pdflatex potentials.tex')
"""
creutz(a,'full')
creutz(a,'proj')
creutz(a,'rem')"""
