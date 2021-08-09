import numpy as np

def eta_init(eta,eta0,eta_up,eta_up0,nx,dx, \
    slope,xl,xb1,xb2,xb3,dbed):
    zb0=xl*slope
    for i in np.arange(0,nx+2):
        xx=dx*float(i)
        eta_up[i]=zb0-xx*slope
        eta_up0[i]=eta_up[i]
#        print(i,nx,eta_up[i])
        if xx>xb1 and xx<xb2:
            ss=xx-xb1
            deta=dbed*ss/(xb2-xb1)
            eta_up[i]=eta_up[i]+deta
        elif xx>=xb2 and xx<xb3:
            ss=xb3-xx
            deta=dbed*ss/(xb3-xb2)
            eta_up[i]=eta_up[i]+deta

    for i in np.arange(1,nx+2):
        eta[i]=(eta_up[i]+eta_up[i-1])*.5
        eta0[i]=(eta_up0[i]+eta_up0[i-1])*.5
    eta[0]=2.*eta[1]-eta[2]
    eta0[0]=2.*eta0[1]-eta0[2]
    return eta,eta0,eta_up,eta_up0

def eta_init_2(eta,eta0,eta_up,eta_up0,nx,dx, \
    xl,x_slope,slope1,slope2):
    zb0=x_slope*slope1+(xl-x_slope)*slope2
    zb1=zb0-x_slope*slope1
    for i in np.arange(0,nx+2):
        xx=dx*float(i)
        if xx <= x_slope:
            eta_up[i]=zb0-xx*slope1
        else:
            eta_up[i]=zb1-(xx-x_slope)*slope2    
        eta_up0[i]=eta_up[i]

    for i in np.arange(1,nx+2):
        eta[i]=(eta_up[i]+eta_up[i-1])*.5
        eta0[i]=(eta_up0[i]+eta_up0[i-1])*.5
    eta[0]=2.*eta[1]-eta[2]
    eta0[0]=2.*eta0[1]-eta0[2]
    return eta,eta0,eta_up,eta_up0

def h_init(eta,eta0,eta_up,eta_up0,h,hs,h_up,hs_up, \
        hs_upstm,hs_dwstm,nx,dx,xl):
    xhalf=xl*.95
    for i in np.arange(0,nx+1):
        xlen=i*dx
        if xlen<xhalf:
            hs_up[i]=hs_upstm
        else:
            hs_up[i]=hs_dwstm

#        hs_up[i]=hs_upstm+(hs_dwstm-hs_upstm)*xlen/xl
        h_up[i]=eta_up0[i]+hs_up[i]

    for i in np.arange(1,nx+1):
        hs[i]=(hs_up[i]+hs_up[i-1])*.5
        h[i]=eta0[i]+hs[i]

    for i in np.arange(0,nx+1):
        hs_up[i]=h_up[i]-eta_up[i]

    for i in np.arange(1,nx+1):
        hs[i]=h[i]-eta[i]

        hs[0]=hs_upstm; h[0]=eta[0]+hs_upstm
        hs[nx+1]=hs_dwstm; h[nx+1]=eta[nx+1]+hs_dwstm

    return h,hs,h_up,hs_up

def h_init_2(eta,eta0,eta_up,eta_up0,h,hs,h_up,hs_up, \
        hs_upstm,hs_dwstm,nx,dx,xl,x_slope):

    for i in np.arange(0,nx+1):
        xlen=i*dx
        if xlen<x_slope:
            hs_up[i]=hs_upstm
        else:
            hs_up[i]=hs_dwstm

#        hs_up[i]=hs_upstm+(hs_dwstm-hs_upstm)*xlen/xl
        h_up[i]=eta_up0[i]+hs_up[i]

    for i in np.arange(1,nx+1):
        hs[i]=(hs_up[i]+hs_up[i-1])*.5
        h[i]=eta0[i]+hs[i]

    for i in np.arange(0,nx+1):
        hs_up[i]=h_up[i]-eta_up[i]

    for i in np.arange(1,nx+1):
        hs[i]=h[i]-eta[i]

        hs[0]=hs_upstm; h[0]=eta[0]+hs_upstm
        hs[nx+1]=hs_dwstm; h[nx+1]=eta[nx+1]+hs_dwstm

    return h,hs,h_up,hs_up

def u_init(g,qp,u,hs_up,fr,nx):
    for i in np.arange(0,nx+1):
        u[i]=qp/hs_up[i]
        fr[i]=u[i]/np.sqrt(g*hs_up[i])
#        print(i,hs_up[i],u[i],fr[i])

    return u,fr

def x_cell_init(x_cell,x,dx,nx):
    for i in np.arange(1,nx+1):
        x_cell[i]=(x[i]+x[i-1])*.5
    x_cell[0]=x_cell[1]-dx
    x_cell[nx+1]=x_cell[nx]+dx

    return x_cell

def h0_cal(eta,eta_up,nx,dx,qp,snm,h0_up):
    for i in np.arange(1,nx):
        slope=(eta[i]-eta[i+1])/dx
        if slope<=0. :
            h0=0.
        else:
            h0=(qp*snm/np.sqrt(slope))**(.6)
        h0_up[i]=eta_up[i]+h0
        if i==1:
            h0_up[0]=eta_up[0]+h0
        elif i==nx-1:
            h0_up[nx]=eta_up[nx]+h0

    return h0_up
