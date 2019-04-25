"""
===================================
Shaded & power normalized rendering
===================================

The Mandelbrot set rendering can be improved by using a normalized recount
associated with a power normalized colormap (gamma=0.3). Rendering can be
further enhanced thanks to shading.

The `maxiter` gives the precision of the computation. `maxiter=200` should
take a few seconds on most modern laptops.

At this point I have customized and beaten on this program so much to get it to do things that the code mash-up can be attributed to myselff as well, now.
-Devin Sweet
"""
import numpy as np
import os
import time
import matplotlib
from matplotlib import colors
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d as mpd
import math

def mandelbrot_set(xmin, xmax, ymin, ymax, zmin, zmax, xn, yn, zn, maxiter, horizon=1.15):
    X = np.linspace(xmin, xmax, xn) #, dtype=np.float32)
    Y = np.linspace(ymin, ymax, yn) #, dtype=np.float32)
    B = np.linspace(zmin, zmax, zn) #, dtype=np.float32)
    C = X + Y[:, None] * B[:, None]*1j  #the "j" is used in Python to denote a complex number.
    N = np.ndarray(C.shape, dtype=int)
    Z = np.ndarray(C.shape, np.complex64)
    for n in range(maxiter):
        I = np.less(abs(Z), horizon)
        N[I] = n
        Z[I] = Z[I]**2 + C[I]
    N[N == maxiter-1] = 0
    return Z, N

def fourleaf(Th):
    ra = math.cos(2*Th)
    return ra

if __name__ == '__main__':
    plt.xscale('linear')
    plt.yscale('linear')
    xmin, xmax, xn = -0.725, -1/1.3, 5000
    ymin, ymax, yn = -1/2.25, -0.1, 5000
    zmin, zmax, zn = -1, 1, 5000
    maxiter = 1000
    horizon = (math.e/2.0) ** (math.pi/2)
    try:
        log_horizon = np.log(np.log(horizon))/np.log(2)
        Z, N = mandelbrot_set(xmin, xmax, ymin, ymax, zmin, zmax, xn, yn, zn, maxiter, horizon)
        with np.errstate(invalid='ignore'):
            M = np.nan_to_num(N + 1 - np.log(np.log(abs(Z)))/np.log(3.1415) + log_horizon)
        Ra = fourleaf(math.e)
        #Rx = M*Ra
        dpi = 800
        width = 6*xn/yn
        height = 6*yn/xn
        #depth = 6
        fig = plt.figure(dpi=dpi)
        fig.set_size_inches(width,height) #, forward=True)
        #ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], frameon=True, aspect=1.5)
        light = colors.LightSource(azdeg=92, altdeg=22)
        M = light.shade(M, cmap=plt.cm.inferno, vert_exag=0.0, norm=colors.PowerNorm(1.1))
        plt.imshow(M, extent=[xmin, xmax, ymin, ymax], aspect='auto', interpolation="Bessel")
        #fig = plt.figure(dpi=dpi) #this line is correct and essential.
        #ax = plt.axes(projection='3d') #this line is essential for a 3d graph
        #ax.scatter(M[0],M[1],M[2], cmap='seismic')
        #ax.set_xlim(-5, 5); ax.set_ylim(-5, 5); ax.set_zlim(-5, 5);
        #ax.set_xlabel('x')
        #ax.set_ylabel('y')
        #ax.set_zlabel('z');
        #ax.set_xticks([-1,-0.5,0,0.5,1])
        #ax.set_yticks([-1,-0.5,0,0.5,1])
        #ax.set_zticks([-1,-0.5,0,0.5,1])
        text=str("figure coords: " + str(xmin) + " xmin, " + str(xmax) + " xmax, " + str(ymin) + " ymin, " + str(ymax) + " ymax, " + str(zmin) + " zmin, " + str(zmax) + " zmax.\n" +
                 "xn = " + str(xn) + ", yn = " + str(yn) + ", zn = " + str(zn) + ", horizon = " + str(horizon) + ", n.iter = " + str(maxiter) + ", dpi = " + str(dpi) ) 
        ax = fig.add_axes = plt.subplots(1)
        #ax = fig.add_subplot(111, projection='polar')
        fig.subplots_adjust(bottom=0.3)
        fig.text(0.1,0.20,text, color="red", fontsize=2, alpha=0.75)
        fig.savefig('mbot.png')
        #plt.subplot_tool(targetfig=fig)
        print("Done.")
        #plt.draw()
        plt.ion()
    except:
        os.exit()
