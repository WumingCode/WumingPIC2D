#load modules
from h5py import *
from numpy import *
from matplotlib.pyplot import *
from matplotlib.colors import *
from matplotlib.cm import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#parameters
dir    ='./'                                    #directory for output
param  = File(dir+'init_param.h5' ,'r')         #read initial condition file
nx     = param.attrs['nx']                      #number of x grird
ny     = param.attrs['ny']                      #number of y grird
np     = param.attrs['np']                      #maximum number of particles in y column
c      = param.attrs['c']                       #speed of light
m      = param.attrs['m']                       #mass
me     = m[1]                                   #electron mass
mi     = m[0]                                   #ion mass
q      = param.attrs['q']                       #charge
e      = q[1]                                   #electron charge
vte    = param.attrs['vte']                     #electron thermal velocity in the upstream frame
vti    = param.attrs['vti']                     #ion thermal velocity in the upstream frame
wpe    = param.attrs['fpe']                     #electron plasma frequency
wpi    = param.attrs['fpi']                     #ion plasma frequency
wge    = param.attrs['fge']                     #electron gyro frequency
wgi    = param.attrs['fgi']                     #ion gyro frequency
va     = param.attrs['va']                      #Alfven velocity
rgi    = param.attrs['rgi']                     #ion gyro radius (Vth,i/wgi)
dx     = param.attrs['delx']                    #grid size
dt     = param.attrs['delt']                    #time interval
beta   = param.attrs['beta']                    #plasma beta
rtemp  = param.attrs['rtemp']                   #Te/Ti
n0     = param.attrs['n0']                      #particle number/cell in the upstream
ls     = param.attrs['ls']                      #electron skin depth in proper frame
b0     = abs(wge*me*c/e)                        #initial magnetic field
gam0   = (ls*wpe/c)**2                          #upstream Lorentz factor
u0     = sqrt(gam0**2-1)                        #upstream four velocity


#read moment data
it   = 1000
dat  = File(dir+'%07d_mom.h5' % it,'r')
bx   = dat['bx']
by   = dat['by']
bz   = dat['bz']
ex   = dat['ex']
ey   = dat['ey']
ez   = dat['ez']
den  = dat['den']
vx   = dat['vx']
vy   = dat['vy']
vz   = dat['vz']
txx  = dat['txx']
tyy  = dat['tyy']
tzz  = dat['tzz']


#plot moment
lsize = 12
tsize = 12
pad   = 0.1
xmin  = 0
xmax  = nx*dx/ls
ymin  = 0
ymax  = ny*dx/ls

figure(figsize=(8,12))
subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.2)

subplot(911)
title(r'$\omega _{pe}t = %5d$' % (wpe/sqrt(gam0)*it*dt) ,fontsize = lsize,y=1.05)
imshow(den[1,:,:]/n0,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=0,vmax=5,aspect='auto')
xlim(xmin,xmax)
ylim(ymin,ymax)
ax = gca()
ax.set_xticklabels([])
yticks(fontsize=tsize)
ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
minorticks_on()
cax = make_axes_locatable(ax).append_axes('right',size='1%', pad=pad)
cbar = colorbar(cax=cax)
cbar.set_ticks(linspace(0,5,6))
cbar.ax.tick_params(labelsize=tsize)
cbar.set_label(r'$N_e/N_1$',fontsize=lsize)

subplot(912)
imshow(den[0,:,:]/n0,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=0, vmax=5,aspect='auto')
xlim(xmin,xmax)
ylim(ymin,ymax)
ax = gca()
ax.set_xticklabels([])
yticks(fontsize=tsize)
ylabel(r'$y /(c/\omega_{pe})$',fontsize=lsize)
minorticks_on()
cax = make_axes_locatable(ax).append_axes('right',size='1%',pad=pad)
cbar = colorbar(cax=cax)
cbar.set_ticks(linspace(0,5,6))
cbar.ax.tick_params(labelsize=tsize)
cbar.set_label(r'$N_i/N_1$',fontsize=lsize)

subplot(913)
imshow(ex/b0,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=-2.5, vmax=2.5,cmap=seismic,aspect='auto')
xlim(xmin,xmax)
ylim(ymin,ymax)
ax = gca()
ax.set_xticklabels([])
yticks(fontsize=tsize)
ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
minorticks_on()
cax = make_axes_locatable(ax).append_axes('right',size='1%', pad=pad)
cbar = colorbar(cax=cax)
cbar.set_ticks(linspace(-2,2,5))
cbar.ax.tick_params(labelsize=tsize)
cbar.set_label(r'$E_x/B_1$',fontsize=lsize)

subplot(914)
imshow(ey/b0,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=-3.5, vmax=1.5,cmap=seismic,aspect='auto')
xlim(xmin,xmax)
ylim(ymin,ymax)
ax = gca()
ax.set_xticklabels([])
yticks(fontsize=tsize)
ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
minorticks_on()
cax = make_axes_locatable(ax).append_axes('right',size='1%', pad=pad)
cbar = colorbar(cax=cax)
cbar.set_ticks(linspace(-3,1,5))
cbar.ax.tick_params(labelsize=tsize)
cbar.set_label(r'$E_y/B_1$',fontsize=lsize)

subplot(915)
imshow(bz/b0,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=-1.5, vmax=3.5,cmap=seismic,aspect='auto')
xlim(xmin,xmax)
ylim(ymin,ymax)
ax = gca()
ax.set_xticklabels([])
yticks(fontsize=tsize)
ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
minorticks_on()
cax = make_axes_locatable(ax).append_axes('right',size='1%', pad=pad)
cbar = colorbar(cax=cax)
cbar.set_ticks(linspace(-1,3,5))
cbar.ax.tick_params(labelsize=tsize)
cbar.set_label(r'$B_z/B_1$',fontsize=lsize)

subplot(916)
imshow(vx[1,:,:]/c,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=-1,vmax=1,cmap=seismic,aspect='auto')
xlim(xmin,xmax)
ylim(ymin,ymax)
yticks(fontsize=tsize)
ylabel(r'$y /(c/\omega_{pe})$',fontsize=lsize)
minorticks_on()
ax = gca()
ax.set_xticklabels([])
cax = make_axes_locatable(ax).append_axes('right',size='1%',pad=pad)
cbar = colorbar(cax=cax)
cbar.ax.tick_params(labelsize=tsize)
cbar.set_label(r'$\beta_{xe}$',fontsize = lsize)

subplot(917)
imshow(vy[1,:,:]/c,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=-1,vmax=1,cmap=seismic,aspect='auto')
xlim(xmin,xmax)
ylim(ymin,ymax)
ax = gca()
ax.set_xticklabels([])
yticks(fontsize=tsize)
ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
minorticks_on()
cax = make_axes_locatable(ax).append_axes('right',size='1%',pad=pad)
cbar = colorbar(cax=cax)
cbar.ax.tick_params(labelsize=tsize)
cbar.set_label(r'$\beta_{ye}$',fontsize=lsize)

subplot(918)
imshow(txx[1,:,:]/gam0/me/c**2,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=0,vmax=1,cmap=hot,aspect='auto')
xlim(xmin,xmax)
ylim(ymin,ymax)
ax = gca()
ax.set_xticklabels([])
yticks(fontsize=tsize)
ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
minorticks_on()
cax = make_axes_locatable(ax).append_axes('right',size='1%', pad=pad)
cbar = colorbar(cax=cax)
cbar.ax.tick_params(labelsize=tsize)
cbar.set_label(r'$k_BT^{xx}_e/\gamma_1m_ec^2$',fontsize=lsize)

subplot(919)
imshow(tyy[1,:,:]/gam0/me/c**2,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=0,vmax=1,cmap=hot,aspect='auto')
xlim(xmin,xmax)
ylim(ymin,ymax)
xticks(fontsize=tsize)
yticks(fontsize=tsize)
xlabel(r'$x /(c/\omega_{pe})$',fontsize=lsize)
ylabel(r'$y /(c/\omega_{pe})$',fontsize=lsize)
minorticks_on()
ax = gca()
cax = make_axes_locatable(ax).append_axes('right',size='1%', pad=pad)
cbar = colorbar(cax=cax)
cbar.ax.tick_params(labelsize=tsize)
cbar.set_label(r'$k_BT^{yy}_e/\gamma_1m_ec^2$',fontsize = lsize)


#read particle data
it   = 1000
dat  = File(dir+'%07d.h5' % it,'r')
up   = dat['up']
xpe  = up[1,:,:,0]
ype  = up[1,:,:,1]
upxe = up[1,:,:,2]
upye = up[1,:,:,3]
xpi  = up[0,:,:,0]
ypi  = up[0,:,:,1]
upxi = up[0,:,:,2]
upyi = up[0,:,:,3]

tmp  = where(xpe>0)
xpe  = xpe[tmp]
ype  = ype[tmp]
upxe = upxe[tmp]
upye = upye[tmp]
tmp  = where(xpi>0)
xpi  = xpi[tmp]
ypi  = ypi[tmp]
upxi = upxi[tmp]
upyi = upyi[tmp]

#plot phase space density integrated over y axis
lsize = 16
tsize = 16
pad   = 0.1
binv  = [1000,100]
xmin  = 0
xmax  = nx*dx/ls

figure(figsize=(10,8))
subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.2)

subplot(411)
title(r'$\omega _{pe}t = %5d$' % (wpe/sqrt(gam0)*it*dt),fontsize=lsize,y=1.05)
hist2d(xpe/ls,upxe/u0,bins=binv,norm=LogNorm(),density=True)
xlim(xmin,xmax)
ylim(-5,5)
ax = gca()
ax.set_xticklabels([])
yticks(linspace(-4,4,5),fontsize=tsize)
ylabel(r'$u_{xe}/u_1$',fontsize=lsize)
minorticks_on()
grid()
cax = make_axes_locatable(ax).append_axes('right',size='1%',pad=pad)
cbar = colorbar(cax=cax)
cbar.ax.tick_params(labelsize=tsize)
cbar.ax.minorticks_off()
cbar.set_label(r'$f_e$',fontsize=lsize)

subplot(412)
hist2d(xpe/ls,upye/u0,bins=binv,norm=LogNorm(),density=True)
xlim(xmin,xmax)
ylim(-5,5)
ax = gca()
ax.set_xticklabels([])
yticks(linspace(-4,4,5),fontsize=tsize)
ylabel(r'$u_{ye}/u_1$',fontsize=lsize)
minorticks_on()
grid()
cax = make_axes_locatable(ax).append_axes('right',size='1%',pad=pad)
cbar = colorbar(cax=cax)
cbar.ax.tick_params(labelsize=tsize)
cbar.ax.minorticks_off()
cbar.set_label(r'$f_e$',fontsize=lsize)

subplot(413)
hist2d(xpi/ls,upxi/u0,bins=binv,norm=LogNorm(),density=True)
xlim(xmin,xmax)
ylim(-5,5)
ax = gca()
ax.set_xticklabels([])
yticks(linspace(-4,4,5),fontsize=tsize)
ylabel(r'$u_{xi}/u_1$',fontsize=lsize)
minorticks_on()
grid()
ax = gca()
cax = make_axes_locatable(ax).append_axes('right',size='1%',pad=pad)
cbar = colorbar(cax=cax)
cbar.ax.tick_params(labelsize=tsize)
cbar.ax.minorticks_off()
cbar.set_label(r'$f_i$',fontsize=lsize)

subplot(414)
hist2d(xpi/ls,upyi/u0,bins=binv,norm=LogNorm(),density=True)
xlim(xmin,xmax)
ylim(-5,5)
xticks(fontsize=tsize)
yticks(linspace(-4,4,5),fontsize=tsize,)
xlabel(r'$x /(c/\omega_{pe})$',fontsize=lsize)
ylabel(r'$u_{yi}/u_1$',fontsize=lsize)
minorticks_on()
grid()
ax = gca()
cax = make_axes_locatable(ax).append_axes('right',size='1%',pad=pad)
cbar = colorbar(cax=cax)
cbar.ax.tick_params(labelsize=tsize)
cbar.set_label(r'$f_i$',fontsize=lsize)
cbar.ax.minorticks_off()


#read orbit data from its to ite
its   = 300
ite   = 1000
intvl = 5
dit   = range(its,ite+intvl,intvl)
dats  = File(dir+'%07d_orb.h5' % its,'r')
ids   = dats['ide'][4700]  #choose tracer particle's id
ndim  = dats['upe'][()].shape[1]
nsize = int((ite-its)/intvl)+1
tpe   = []
for it in dit:
    dat = File(dir+'%07d_orb.h5' % it,'r')
    tmp = dat['upe'][()]
    tpe = append(tpe,tmp[where(dat['ide'][()]==ids)])

tpe  = tpe.reshape(nsize,ndim)

#plot electron orbit
lsize = 16
tsize = 16

figure()
plot(tpe[:,0]/ls,tpe[:,1]/ls,'-k',lw=1, alpha=0.4)
scatter(tpe[:,0]/ls,tpe[:,1]/ls,marker='.',s=10,c=dit*wpe/sqrt(gam0))
xticks(fontsize=tsize)
yticks(fontsize=tsize)
xlabel(r'$x /(c/\omega_{pe})$',fontsize=lsize)
ylabel(r'$y /(c/\omega_{pe})$',fontsize=lsize)
minorticks_on()
grid()
ax = gca()
ax.set_aspect('equal', adjustable='box')
cax = inset_axes(ax,width='40%',height='5%',loc=1,borderpad=2)
cbar = colorbar(cax=cax,orientation='horizontal')
cbar.ax.tick_params(labelsize=tsize)
cbar.set_label(r'$\omega_{pe}t$',fontsize=lsize)

show()
