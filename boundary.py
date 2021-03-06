import numpy as np
# j_upstm, j_dwstm .... 1=Wall, 2=Free, 3=Fixed Value(Flowing)

def h_bound(h,hs,eta,hs_upstm,hs_dwstm,nx,j_upstm,j_dwstm):
    if j_upstm==2:
        hs[0]=hs[1]

    if j_dwstm==3:
        hs[nx+1]=hs_dwstm
    else:
        hs[nx+1]=hs[nx]

    h[0]=eta[0]+hs[0]
    h[nx+1]=eta[nx+1]+hs[nx+1]
    return h,hs

def u_bound(u,hs_up,qp,nx,j_upstm,j_dwstm,u_upstm,u_dwstm):
    if j_upstm==1:
        u[0]=0.
    else:
        u[0]=qp/hs_up[0]

    if j_dwstm==1:
        u[nx]=0.
    return u

def hs_up_cal(hs,hs_up,nx,hs_upstm,hs_dwstm,j_upstm,j_dwstm):
    for i in np.arange(0,nx+1):
        hs_up[i]=(hs[i]+hs[i+1])*.5

#    if j_upstm==3: 
#        hs_up[0]=hs_upstm
#    if j_dwstm==3:
#        hs_up[nx]=hs_dwstm

    return hs_up

def h_up_cal(hs_up,eta_up,h_up,nx):
    for i in np.arange(0,nx+1):
        h_up[i]=eta_up[i]+hs_up[i]
    return h_up

def gbound_u(gux,nx):
    gux[0]=0.
    gux[nx]=0.
    return gux
