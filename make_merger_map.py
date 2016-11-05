import numpy as np
import pyfits
import astropy
from astropy.convolution import convolve_fft, Gaussian1DKernel, Gaussian2DKernel
import scipy
from scipy.optimize import curve_fit

def gauss(x, *p):
    A, mu, sigma= p
    return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))



plt.ioff()


plt.close('all')



N = 50
q = 40

orig_flux_map = zeros((N,N))
orig_vel_map  = zeros((N,N))
orig_disp_map = zeros((N,N))

new_vel_map  = zeros((N,N,2))
new_disp_map = zeros((N,N,2))


f1 = 100
f2 = 100

orig_flux_map[N/2,int(0.38*N)] = f1
orig_flux_map[N/2,int(0.62*N)] = f2

m = 0.68*f2/(sqrt(pi)*N/20.)**2./10.
print m

kern = Gaussian2DKernel(N/20.)
orig_flux_map = convolve_fft(orig_flux_map, kern)
orig_flux_map += np.random.normal(0, m/4., shape(orig_flux_map))

orig_vel_map[:, N/2:][orig_flux_map[:, N/2:] > m] = np.random.normal(13,0.5,shape(orig_vel_map[:, N/2:][orig_flux_map[:, N/2:] > m]))
orig_vel_map[:, :N/2][orig_flux_map[:, :N/2] > m] = np.random.normal(17,0.5,shape(orig_vel_map[:, :N/2][orig_flux_map[:, :N/2] > m]))

orig_disp_map[:, N/2:][orig_flux_map[:, N/2:] > m] = np.random.normal(1.5,0.3,shape(orig_disp_map[:, N/2:][orig_flux_map[:, N/2:] > m]))
orig_disp_map[:, :N/2][orig_flux_map[:, :N/2] > m] = np.random.normal(1.5,0.3,shape(orig_disp_map[:, :N/2][orig_flux_map[:, :N/2] > m]))


orig_cube = zeros((N,N,q))
new_cube  = zeros((N,N,q))
for i in arange(N):
	for j in arange(N):
		orig_cube[i,j,orig_vel_map[i,j].astype('int')] = orig_flux_map[i,j]
		krn_1d = Gaussian1DKernel(orig_disp_map[i,j])
		orig_cube[i,j,:] = convolve_fft(orig_cube[i,j,:], krn_1d)

for i in arange(q):
	krn_2d = Gaussian2DKernel(N/8.)
	new_cube[:,:,i] = convolve_fft(orig_cube[:,:,i], krn_2d)
	new_cube[:,:,i] += np.random.normal(0, m/100., shape(orig_flux_map))


new_flux_map = sum(new_cube, axis = 2)

for i in arange(N):
	for j in arange(N):
		try:
			p_a, v_a = curve_fit(gauss, arange(q), new_cube[i,j], p0 = [max(new_cube[i,j]), 15, 2])
			if sqrt(v_a[1,1]) < 0.5:
				new_vel_map[i,j,0] = p_a[1]
				new_disp_map[i,j,0] = p_a[2]
				new_vel_map[i,j,1] = sqrt(v_a[1,1])
				new_disp_map[i,j,1] = sqrt(v_a[2,2])

		except:
			pass



fig = figure()

ax1 = fig.add_subplot(341)
ax2 = fig.add_subplot(342)
ax3 = fig.add_subplot(343)
ax7 = fig.add_subplot(344)

ax4 = fig.add_subplot(345)
ax5 = fig.add_subplot(346)
ax6 = fig.add_subplot(347)
ax8 = fig.add_subplot(348)

ax9 = fig.add_subplot(3,4,10)
ax10 = fig.add_subplot(3,4,11)


ax1.imshow(orig_flux_map, interpolation = 'nearest', cmap = 'Greys_r')
ax2.imshow(orig_vel_map, interpolation = 'nearest', vmin = 12, vmax = 18)
ax3.imshow(orig_disp_map, interpolation = 'nearest', cmap = 'viridis', vmin = 0.5, vmax = 2.5)
ax4.imshow(new_flux_map, interpolation = 'nearest', cmap = 'Greys_r')
ax5.imshow(new_vel_map[:,:,0], interpolation = 'nearest',  vmin = 12, vmax = 18)
ax6.imshow(new_disp_map[:,:,0], interpolation = 'nearest', cmap = 'viridis', vmin = 0.5, vmax = 2.5)

flip_new_vel_map = fliplr(new_vel_map[:,:,0])
flip_new_disp_map = fliplr(new_disp_map[:,:,0])

ax9.imshow(new_vel_map[:,:,1], interpolation = 'nearest',  cmap = 'Greys_r', vmin = 0, vmax = 0.5)
ax10.imshow(new_disp_map[:,:,1], interpolation = 'nearest',  cmap = 'Greys_r',  vmin = 0., vmax = 0.5)

#ax9.imshow(flip_new_vel_map, interpolation = 'nearest',  vmin = 12, vmax = 18)
#ax10.imshow(flip_new_disp_map, interpolation = 'nearest', cmap = 'viridis', vmin = 0.5, vmax = 2.5)




ax7.plot(new_vel_map[N/2,:,0],'o')
ax7.plot(orig_vel_map[N/2,:],'k.')

ax8.plot(new_disp_map[N/2,:,0],'o')
ax8.plot(orig_disp_map[N/2,:],'k.')

ax7.set_ylim(10,20)
ax8.set_ylim(0,3.5)


ax1.set_axis_off()
ax2.set_axis_off()
ax3.set_axis_off()
ax4.set_axis_off()
ax5.set_axis_off()
ax6.set_axis_off()

ax9.set_axis_off()
ax10.set_axis_off()


ax7.set_xticklabels([''])
ax8.set_xticklabels([''])
ax7.set_yticklabels([''])
ax8.set_yticklabels([''])


#plt.subplots_adjust(hspace = 0, wspace = 0)
#ax8.axis('equal')

fig.tight_layout()
plt.subplots_adjust(hspace = 0, wspace = 0)
plt.savefig('test.png', dpi = 200)













