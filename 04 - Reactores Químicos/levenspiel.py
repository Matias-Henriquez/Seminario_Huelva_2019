import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
import scipy.integrate as integrate

Kinetics = dict()
Kinetics['1er Orden'] = lambda x, k, cA0: 1/(k*cA0*(1-x))
Kinetics['2º Orden'] = lambda x, k, cA0: 1/(k*cA0**2*(1-x)**2)
Kinetics['Langmuir-Hinshelwood'] = lambda x, k, cA0: 1/(k*cA0*(1-x)/(0.5+1.2*k*cA0*(1-x))**2)


def plot(kinetics='Langmuir-Hinshelwood', k=0.5, X1=0.4, X2=0.8, X3=0.9, Reactor1='CSTR', Reactor2='CSTR', Reactor3='CSTR'):

    x=np.linspace(0,0.99,100)
    fig, ax = plt.subplots()
    r=Kinetics[kinetics](x,1.0,k)
    plt.plot(x,r)
    plt.xlabel("$X_A$")
    plt.ylabel("$F_{A0}/r_A$")
    plt.ylim(bottom=0, top=20)
    plt.xlim(0,1)
    
    if (Reactor1=='CSTR'):
        y1=Kinetics[kinetics](X1,1.0,k)
        V1=X1*y1
        rect = Rectangle((0, 0), X1, y1, linewidth=1,edgecolor='y',facecolor='y', alpha=0.5)
        ax.add_patch(rect)
        print("REACTOR1 CSTR Volume=%.3g L"%V1)
    else:
        V1=integrate.quad(lambda x: Kinetics[kinetics](x,1.0,k), 0, X1)[0]
        ix=np.linspace(0,X1,100)
        iy=Kinetics[kinetics](ix,1.0,k)
        verts = [(0, 0), *zip(ix, iy), (X1, 0)]
        poly = Polygon(verts, facecolor='y', edgecolor='y', alpha=0.5)
        ax.add_patch(poly)
        print("REACTOR1 PFR Volume=%.3g L"%V1)

    if (Reactor2=='CSTR'):
        y2=Kinetics[kinetics](X2,1.0,k)
        V2=(X2-X1)*y2
        rect = Rectangle((X1, 0), X2-X1, y2, linewidth=1,edgecolor='c',facecolor='c', alpha=0.5)
        ax.add_patch(rect)
        print("REACTOR2 CSTR Volume=%.3g L"%V2)
    else:
        V2=integrate.quad(lambda x: Kinetics[kinetics](x,1.0,k), X1, X2)[0]
        ix=np.linspace(X1,X2,100)
        iy=Kinetics[kinetics](ix,1.0,k)
        verts = [(X1, 0), *zip(ix, iy), (X2, 0)]
        poly = Polygon(verts, facecolor='c', edgecolor='c', alpha=0.5)
        ax.add_patch(poly)
        print("REACTOR2 PFR Volume=%.3g L"%V2)
      
    if (Reactor3=='CSTR'):
        y3=Kinetics[kinetics](X3,1.0,k)
        V3=(X3-X2)*y3
        rect = Rectangle((X2, 0), X3-X2, y3, linewidth=1,edgecolor='r',facecolor='r', alpha=0.5)
        ax.add_patch(rect)
        print("REACTOR3 CSTR Volume=%.3g L"%V3)
    else:
        V3=integrate.quad(lambda x: Kinetics[kinetics](x,1.0,k), X2, X3)[0]
        ix=np.linspace(X2,X3,100)
        iy=Kinetics[kinetics](ix,1.0,k)
        verts = [(X2, 0), *zip(ix, iy), (X3, 0)]
        poly = Polygon(verts, facecolor='r', edgecolor='r', alpha=0.5)
        ax.add_patch(poly)
        print("REACTOR3 PFR Volume=%.3g L"%V3)
        
        
    print("TOTAL Volume=%.3g L"%(V1+V2+V3))
        
    plt.show()
