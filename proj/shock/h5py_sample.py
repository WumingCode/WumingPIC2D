# For batch execution on a headless machine, run as
# $ MPLBACKEND=agg python h5py_sample.py

#load modules
import numpy as np
import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def test_param(fn):
    #parameters
    print_param = {
        'nx' : 'number of grid in x',
        'ny' : 'number of grid in y',
        'np' : 'maximum number of particles in y',
        'c' : 'speed of light',
        'r' : 'mass',
        'q' : 'charge',
        'fpi' : 'ion plasma frequency',
        'fpe' : 'electron plasma frequency',
        'fgi' : 'ion gyro frequency',
        'fge' : 'electron gyro frequency',
        'vti' : 'ion thermal velocity',
        'vte' : 'electron thermal velocity',
        'vai' : 'ion Alfven velocity',
        'vae' : 'electron Alfven velocity',
        'delx' : 'grid size',
        'delt' : 'time step',
        'n0' : 'number of particle / cell',
        'ls' : 'electron skin depth',
        'b0' : 'initial magnetic field strength',
        'gam0' : 'upstram Lorentz factor',
    }

    # read all parameters
    with h5py.File(fn, 'r') as f:
        param = dict(f.attrs)

    # some additional parameters
    c     = param['c']
    qi    = param['q'][0]
    qe    = param['q'][1]
    mi    = param['r'][0]
    me    = param['r'][1]
    wpe   = param['fpe']
    wpi   = param['fpi']
    wge   = param['fge']
    wgi   = param['fgi']
    param['b0']   = abs(wpe*me*c/qe)
    param['gam0'] = (param['ls']*wpe/c)**2
    param['u0']   = np.sqrt(param['gam0']**2-1)

    print('*** print parameters ***')
    for key, desc in print_param.items():
        print('- {:40s} : {:}'.format(desc, param[key]))

    return param


def test_moment(fn, it, param, batch=True):
    # read moment data
    with h5py.File(fn, 'r') as dat:
        uf   = dat['uf']
        bx   = uf[...,0]
        by   = uf[...,1]
        bz   = uf[...,2]
        ex   = uf[...,3]
        ey   = uf[...,4]
        ez   = uf[...,5]
        den  = dat['den'][...]
        vx   = dat['vel'][:,0,...]
        vy   = dat['vel'][:,1,...]
        vz   = dat['vel'][:,2,...]
        txx  = dat['temp'][:,0,...]
        tyy  = dat['temp'][:,1,...]
        tzz  = dat['temp'][:,2,...]

    # plot moment
    nx    = param['nx']
    ny    = param['ny']
    n0    = param['n0']
    dt    = param['delt']
    dx    = param['delx']
    ls    = param['ls']
    c     = param['c']
    mi    = param['r'][0]
    me    = param['r'][1]
    wpe   = param['fpe']
    wpi   = param['fpi']
    wge   = param['fge']
    wgi   = param['fgi']
    b0    = param['b0']
    gam0  = param['gam0']
    u0    = param['u0']
    lsize = 12
    tsize = 12
    pad   = 0.1
    xmin  = 0
    xmax  = nx*dx/ls
    ymin  = 0
    ymax  = ny*dx/ls

    plt.figure(figsize=(8,12))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.2)
    cmap1 = plt.cm.seismic
    cmap2 = plt.cm.hot

    plt.subplot(911)
    plt.title(r'$\omega _{pe}t = %5d$' % (wpe/np.sqrt(gam0)*it*dt) ,fontsize = lsize,y=1.05)
    plt.imshow(den[1,:,:]/n0,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=0,vmax=5,aspect='auto')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    ax.set_xticklabels([])
    plt.yticks(fontsize=tsize)
    plt.ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='1%', pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(0,5,6))
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$N_e/N_1$',fontsize=lsize)

    plt.subplot(912)
    plt.imshow(den[0,:,:]/n0,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=0, vmax=5,aspect='auto')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    ax.set_xticklabels([])
    plt.yticks(fontsize=tsize)
    plt.ylabel(r'$y /(c/\omega_{pe})$',fontsize=lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='1%',pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(0,5,6))
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$N_i/N_1$',fontsize=lsize)

    plt.subplot(913)
    plt.imshow(ex/b0,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=-2.5, vmax=2.5,cmap=cmap1,aspect='auto')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    ax.set_xticklabels([])
    plt.yticks(fontsize=tsize)
    plt.ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='1%', pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(-2,2,5))
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$E_x/B_1$',fontsize=lsize)

    plt.subplot(914)
    plt.imshow(ey/b0,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=-3.5, vmax=1.5,cmap=cmap1,aspect='auto')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    ax.set_xticklabels([])
    plt.yticks(fontsize=tsize)
    plt.ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='1%', pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(-3,1,5))
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$E_y/B_1$',fontsize=lsize)

    plt.subplot(915)
    plt.imshow(bz/b0,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=-1.5, vmax=3.5,cmap=cmap1,aspect='auto')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    ax.set_xticklabels([])
    plt.yticks(fontsize=tsize)
    plt.ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='1%', pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.set_ticks(np.linspace(-1,3,5))
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$B_z/B_1$',fontsize=lsize)

    plt.subplot(916)
    plt.imshow(vx[1,:,:]/c,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=-1,vmax=1,cmap=cmap1,aspect='auto')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.yticks(fontsize=tsize)
    plt.ylabel(r'$y /(c/\omega_{pe})$',fontsize=lsize)
    plt.minorticks_on()
    ax = plt.gca()
    ax.set_xticklabels([])
    cax = make_axes_locatable(ax).append_axes('right',size='1%',pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$\beta_{xe}$',fontsize = lsize)

    plt.subplot(917)
    plt.imshow(vy[1,:,:]/c,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=-1,vmax=1,cmap=cmap1,aspect='auto')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    ax.set_xticklabels([])
    plt.yticks(fontsize=tsize)
    plt.ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='1%',pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$\beta_{ye}$',fontsize=lsize)

    plt.subplot(918)
    plt.imshow(txx[1,:,:]/gam0/me/c**2,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=0,vmax=1,cmap=cmap2,aspect='auto')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax = plt.gca()
    ax.set_xticklabels([])
    plt.yticks(fontsize=tsize)
    plt.ylabel(r'$y /(c/\omega_{pe})$',fontsize= lsize)
    plt.minorticks_on()
    cax = make_axes_locatable(ax).append_axes('right',size='1%', pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$k_BT^{xx}_e/\gamma_1m_ec^2$',fontsize=lsize)

    plt.subplot(919)
    plt.imshow(tyy[1,:,:]/gam0/me/c**2,extent=[0,nx*dx/ls,ny*dx/ls,0],vmin=0,vmax=1,cmap=cmap2,aspect='auto')
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.xticks(fontsize=tsize)
    plt.yticks(fontsize=tsize)
    plt.xlabel(r'$x /(c/\omega_{pe})$',fontsize=lsize)
    plt.ylabel(r'$y /(c/\omega_{pe})$',fontsize=lsize)
    plt.minorticks_on()
    ax = plt.gca()
    cax = make_axes_locatable(ax).append_axes('right',size='1%', pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$k_BT^{yy}_e/\gamma_1m_ec^2$',fontsize = lsize)

    if batch:
        plt.savefig('moment.png')


def test_particle(fn, it, param, batch=True):
    # read particle data
    with h5py.File(fn, 'r') as dat:
        up   = dat['up']
        xpe  = up[1,:,:,0]
        ype  = up[1,:,:,1]
        upxe = up[1,:,:,2]
        upye = up[1,:,:,3]
        xpi  = up[0,:,:,0]
        ypi  = up[0,:,:,1]
        upxi = up[0,:,:,2]
        upyi = up[0,:,:,3]

    tmp  = np.where(xpe>0)
    xpe  = xpe[tmp]
    ype  = ype[tmp]
    upxe = upxe[tmp]
    upye = upye[tmp]
    tmp  = np.where(xpi>0)
    xpi  = xpi[tmp]
    ypi  = ypi[tmp]
    upxi = upxi[tmp]
    upyi = upyi[tmp]

    # plot phase space density integrated over y axis
    nx    = param['nx']
    ny    = param['ny']
    n0    = param['n0']
    dt    = param['delt']
    dx    = param['delx']
    ls    = param['ls']
    c     = param['c']
    mi    = param['r'][0]
    me    = param['r'][1]
    wpe   = param['fpe']
    wpi   = param['fpi']
    wge   = param['fge']
    wgi   = param['fgi']
    b0    = param['b0']
    gam0  = param['gam0']
    u0    = param['u0']
    lsize = 16
    tsize = 16
    pad   = 0.1
    binv  = [1000,100]
    xmin  = 0
    xmax  = nx*dx/ls
    norm  = mpl.colors.LogNorm

    plt.figure(figsize=(10,8))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, hspace=0.2)

    plt.subplot(411)
    plt.title(r'$\omega _{pe}t = %5d$' % (wpe/np.sqrt(gam0)*it*dt),fontsize=lsize,y=1.05)
    plt.hist2d(xpe/ls,upxe/u0,bins=binv,norm=norm(),density=True)
    plt.xlim(xmin,xmax)
    plt.ylim(-5,5)
    ax = plt.gca()
    ax.set_xticklabels([])
    plt.yticks(np.linspace(-4,4,5),fontsize=tsize)
    plt.ylabel(r'$u_{xe}/u_1$',fontsize=lsize)
    plt.minorticks_on()
    plt.grid()
    cax = make_axes_locatable(ax).append_axes('right',size='1%',pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.ax.tick_params(labelsize=tsize)
    cbar.ax.minorticks_off()
    cbar.set_label(r'$f_e$',fontsize=lsize)

    plt.subplot(412)
    plt.hist2d(xpe/ls,upye/u0,bins=binv,norm=norm(),density=True)
    plt.xlim(xmin,xmax)
    plt.ylim(-5,5)
    ax = plt.gca()
    ax.set_xticklabels([])
    plt.yticks(np.linspace(-4,4,5),fontsize=tsize)
    plt.ylabel(r'$u_{ye}/u_1$',fontsize=lsize)
    plt.minorticks_on()
    plt.grid()
    cax = make_axes_locatable(ax).append_axes('right',size='1%',pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.ax.tick_params(labelsize=tsize)
    cbar.ax.minorticks_off()
    cbar.set_label(r'$f_e$',fontsize=lsize)

    plt.subplot(413)
    plt.hist2d(xpi/ls,upxi/u0,bins=binv,norm=norm(),density=True)
    plt.xlim(xmin,xmax)
    plt.ylim(-5,5)
    ax = plt.gca()
    ax.set_xticklabels([])
    plt.yticks(np.linspace(-4,4,5),fontsize=tsize)
    plt.ylabel(r'$u_{xi}/u_1$',fontsize=lsize)
    plt.minorticks_on()
    plt.grid()
    ax = plt.gca()
    cax = make_axes_locatable(ax).append_axes('right',size='1%',pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.ax.tick_params(labelsize=tsize)
    cbar.ax.minorticks_off()
    cbar.set_label(r'$f_i$',fontsize=lsize)

    plt.subplot(414)
    plt.hist2d(xpi/ls,upyi/u0,bins=binv,norm=norm(),density=True)
    plt.xlim(xmin,xmax)
    plt.ylim(-5,5)
    plt.xticks(fontsize=tsize)
    plt.yticks(np.linspace(-4,4,5),fontsize=tsize,)
    plt.xlabel(r'$x /(c/\omega_{pe})$',fontsize=lsize)
    plt.ylabel(r'$u_{yi}/u_1$',fontsize=lsize)
    plt.minorticks_on()
    plt.grid()
    ax = plt.gca()
    cax = make_axes_locatable(ax).append_axes('right',size='1%',pad=pad)
    cbar = plt.colorbar(cax=cax)
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$f_i$',fontsize=lsize)
    cbar.ax.minorticks_off()

    if batch:
        plt.savefig('particle.png')


def test_orbit(fns, its, param, batch=True):
    ptcl_id   = 4800
    ptcl_name = 'up02'
    with h5py.File(fns[0], 'r') as f:
        upe  = f[ptcl_name]
        ndim = upe.shape[1]
        ide  = upe[ptcl_id,ndim-1]

    tpe = np.zeros((len(fns), ndim))
    for i, f in enumerate(fns):
        with h5py.File(f, 'r') as f:
            ptcl = f[ptcl_name][()]
            tpe[i,:] = ptcl[np.where(ptcl[:,ndim-1]==ide)]

    # plot electron orbit
    dt    = param['delt']
    ls    = param['ls']
    wpe   = param['fpe']
    gam0  = param['gam0']
    lsize = 16
    tsize = 16

    plt.figure()
    plt.plot(tpe[:,0]/ls,tpe[:,1]/ls,'-k',lw=1, alpha=0.4)
    plt.scatter(tpe[:,0]/ls,tpe[:,1]/ls,marker='.',s=10,c=dt*its*wpe/np.sqrt(gam0))
    plt.xticks(fontsize=tsize)
    plt.yticks(fontsize=tsize)
    plt.xlabel(r'$x /(c/\omega_{pe})$',fontsize=lsize)
    plt.ylabel(r'$y /(c/\omega_{pe})$',fontsize=lsize)
    plt.minorticks_on()
    plt.grid()
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    cax = inset_axes(ax,width='40%',height='5%',loc=1,borderpad=2)
    cbar = plt.colorbar(cax=cax,orientation='horizontal')
    cbar.ax.tick_params(labelsize=tsize)
    cbar.set_label(r'$\omega_{pe}t$',fontsize=lsize)

    if batch:
        plt.savefig('orbit.png')


if __name__ == '__main__':
    dir  = './'
    it1  = 300
    it2  = 1000
    itv  = 5
    its  = np.arange(it1, it2+itv, itv)
    fn0  = dir+'init_param.h5'
    fn1  = dir+'{:07d}_mom.h5'.format(it2)
    fn2  = dir+'{:07d}.h5'.format(it2)
    fns  = [dir+'{:07d}_orb.h5'.format(it) for it in its]

    param = test_param(fn0)
    test_moment(fn1, it2, param)
    test_particle(fn2, it2, param)
    test_orbit(fns, its, param)

    if mpl.get_backend() != 'agg':
        plt.show()
