import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import animation
import time

# some interesting places in the set
# http://www.nahee.com/Derbyshire/manguide.html

N = 500
nIts = 250
nZooms = 50
x0=0
y0=-1

movie = np.zeros([N,N,nZooms])

x = np.linspace(-2,1,N)
y = np.linspace(-1,1,N)
X,Y = np.meshgrid(x,y)
c = X + 1j*Y
z = 0*c
for i in range(nIts):
    z = z**2 + c

mask = np.abs(z) < 1
z[z>1]=0
z[np.isnan(z)]=0
movie[:,:,0] = mask


# plotting stuff
for j in range(1,nZooms):
    h=1./(2**j)
    print("Plot number ", j)
    x = x0+h*np.linspace(-1,1,N)
    y = y0+h*np.linspace(-1,1,N)
    X,Y = np.meshgrid(x,y)
    c = X + 1j*Y
    z = 0*c
    for i in range(nIts):
        z = z**2 + c

    mask = np.abs(z) < 1
    z[z>1]=0
    z[np.isnan(z)]=0
    movie[:,:,j] = mask

fig = plt.figure()

for j in range(nZooms):
    name = "mandelbrotimage%d.png" % j
    plt.imshow(movie[:,:,j], cmap = 'RdBu')
    plt.gray()
    plt.axis('equal')
    plt.axis('off')
    #plt.show()
    fig.savefig(name)
    time.sleep(.001)
