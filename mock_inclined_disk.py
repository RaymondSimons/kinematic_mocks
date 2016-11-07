import numpy as np

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm


def atan_rc(R, vflat = 150., rv = 1.):
	return vflat*2/pi*arctan(R/rv)
#	return vflat*R/rv



incl = 50.*pi/180.
PA = 0

N = 100
X, Y = meshgrid(linspace(0,N,N), linspace(0,N,N))
X0 = X[len(X)/2, len(X)/2]
Y0 = Y[len(Y)/2, len(Y)/2]
R = sqrt((X-X0)**2.+(Y-Y0)**2.)#/sin(incl)


costheta = ((-(X-X0)*sin(PA))+(Y-Y0)*cos(PA))/R

V = atan_rc(R, vflat = 200., rv = 0.1*max(X.ravel()))*sin(incl)*costheta

#V+=np.random.normal(0,40,shape(V))


#imshow(V, interpolation = 'nearest', vmin = -100, vmax = 100)

fig = plt.figure(1)
clf()
ax = fig.gca(projection='3d')
ax2 = ax.zaxis.get_axes()
#good = [0.4*X.shape[0]:0.6*X.shape[1],:,:]
width = 0.5
minx = int((0.5-width)*X.shape[0])
maxx = int((0.5+width)*X.shape[0])

#ax2.add_image(imshow(V, alpha = 0.1))




ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Velocity')
ax.set_zlim(-200,200)
#ax.contourf(X[:,minx:maxx], Y[:, minx:maxx], V[:, minx:maxx], zdir = 'x', offset=ax.get_xlim()[0], cmap = 'Greys')
ax.contourf(X[:,minx:maxx], Y[:, minx:maxx], V[:, minx:maxx], zdir = 'z', offset=ax.get_zlim()[0], cmap = 'jet', interpolation = 'nearest')

ax.plot_surface(X[:,minx:maxx], Y[:, minx:maxx], V[:, minx:maxx], alpha=0.2)

ax.view_init(elev = -1, azim = -1)
