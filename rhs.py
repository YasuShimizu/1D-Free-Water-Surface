import numpy as np

def un_cal(un,u,nx,dx,cfx,hn,g,dt):
    for i in np.arange(0,nx+1):
        dhdx=(hn[i+1]-hn[i])/dx
        un[i]=u[i]+(cfx[i]-g*dhdx)*dt
    return un

def qu_cal(qu,un,hs_up,nx):
    for i in np.arange(0,nx+1):
        qu[i]=un[i]*hs_up[i]
    return qu