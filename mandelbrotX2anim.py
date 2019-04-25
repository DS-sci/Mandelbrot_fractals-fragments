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
#THe reason you have to use numpy is to define matrix arrays of X and Y  coordinates.  It's an Nx2 matrix normally, where N is the number of iterations.

def mandelbrot_set(xmin, xmax, ymin, ymax, zmin, zmax, xn, yn, zn, maxiter, horizon=1.15):
    X = np.linspace(xmin, xmax, xn)#, dtype=np.float32)
    Y = np.linspace(ymin, ymax, yn)#, dtype=np.float32)
    B = np.linspace(zmin, zmax, zn)#, dtype=np.float32)
    C = X[:, None]*1j + Y[:, None]*1j + B[:, None]*1j  #the "j" is used in Python to denote a complex number.
    N = np.ndarray.zeros(C.shape, dtype=int)
    Z = np.ndarray.zeros(C.shape, np.complex64)
    for n in range(maxiter):
        I = np.less(abs(Z), horizon)
        N[I] = n
        Z[I] = Z[I]**3 + C[I]
    N[N == maxiter-1] = 0
    return Z, N


if __name__ == '__main__':
    import time
    import matplotlib
    from matplotlib import colors
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d as mpd
    from matplotlib.tri import triangulation as tria
    plt.xscale('log')
    plt.yscale('linear')
    #the xn,yn,zn values are the number of points calculated between max and min.  This is essentially the resolution of the image.
    xmin, xmax, xn = -2.15, -1.75, 2000
    ymin, ymax, yn = 0, +2.75, 2000
    zmin, zmax, zn = 0, +2.75, 2000
    maxiter = 100
    horizon = 2.0 ** 3.15
    try:
        log_horizon = np.log(np.log(horizon))/np.log(2)
        Z, N = mandelbrot_set(xmin, xmax, ymin, ymax, zmin, zmax, xn, yn, zn, maxiter, horizon)

        # Normalized recount as explained in:
        # https://linas.org/art-gallery/escape/smooth.html
        # https://www.ibm.com/developerworks/community/blogs/jfp/entry/My_Christmas_Gift

        # This line will generate warnings for null values but it is faster to
        # process them afterwards using the nan_to_num
        with np.errstate(invalid='ignore'):
            M = np.nan_to_num(N + 1 -
                              np.log(np.log(abs(Z)))/np.log(3.1415) +
                              log_horizon)

        dpi = 320
        width = 6
        height = 6*yn/xn
        depth = 6
        fig = plt.figure(dpi=dpi) #this line is correct and essential.
        ax = plt.axes(projection='3d') #this line is essential for a 3d graph
        #tri = tria.Triangulation(np.ravel(M[0]), np.ravel(M[1]), np.ravel(M[2])) #This produces an image made up of triangles.  It is perhaps not the best solution.
        tri = tria.Triangulation(M[0], M[1], M[2]) #second variation.  I don't actually know what np.ravel does.
        ax.plot_trisurf(M[0],M[1],M[2], triangles=tri.triangles, cmap='viridis', edgecolor='none');
        ax.set_xlim(-5, 5); ax.set_ylim(-5, 5); ax.set_zlim(-5, 5);
        #pd = mplot3d.art3d.Poly3DCollection(M[0],M[1],M[2]) #this line never worked
        #ax.scatter(M[0],M[1],M[2], cmap='seismic') #this line worked but didn't look right
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z');
        #fig.add_axes([0.0, 0.0, 0.0, 1.0, 1.0, 1.0], frameon=False, aspect=1) #this line no longer works because of teh 3d plot

        # Shaded rendering
        #light = colors.LightSource(azdeg=180, altdeg=15) #this line no longe works because of the 3d plot
        #M = light.shade(data=M[::], cmap=plt.cm.seismic, blend_mode="overlay", vmin=None, vmax=None, vert_exag=0, dx=0.001, dy=0.001, fraction=1.15)
        #plt.imshow(M, extent=[xmin, xmax, ymin, ymax, zmin, zmax], interpolation="gaussian")
        ax.set_xticks([-1,-0.5,0,0.5,1])
        ax.set_yticks([-1,-0.5,0,0.5,1])
        ax.set_zticks([-1,-0.5,0,0.5,1])
        fig.savefig('mbot.tiff')
        
        #fig = pd.do_3d_projection()

        # Some advertisement for matplotlib
        #year = time.strftime("%Y")
        #major, minor, micro = matplotlib.__version__.split('.', 2)
        #text = ("The Mandelbrot fractalset\n" "Rendered with matplotlib %s.%s, %s - http://matplotlib.org"
        #% (major, minor, year))
        #ax.text(xmin+.025, ymin+.025, text,color="white", fontsize=12, alpha=0.5)

    except:
        #fig.savefig('mbotEXCEPTION.tiff')
        #print("an exception occurred",exception)
        
    #fig.savefig('mbot.tiff')
    #fig2.savefig('mbottt.tiff')
    '''
    def animate(i):
        x = np.linspace(0, 2, 1000)
        y = np.sin(2 * np.pi * (x - 0.01 * i))
        line.set_data(x, y)
        return line,
    try:
        # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=20, interval=20, blit=True)

        # save the animation as an mp4.  This requires ffmpeg or mencoder to be
        # installed.  The extra_args ensure that the x264 codec is used, so that
        # the video can be embedded in html5.  You may need to adjust this for
        # your system: for more information, see
        # http://matplotlib.sourceforge.net/api/animation_api.html
        anim.save('basic_animation.mp4', fps=10, extra_args=['-vcodec', 'libx264'])
    except:
        print("an exception ocurred in animate().")
        '''
    plt.subplot_tool(targetfig=fig)
    print("Done.") #plt.show()
    
