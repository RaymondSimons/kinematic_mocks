import pyfits
plt.close('all')
plt.ioff()


fig = figure(figsize = (10,6))

n = 20

for i in arange(n):
	ax = fig.add_subplot(n,1,i+1)
	ax.set_xticklabels([''])
	ax.set_yticklabels([''])
	ax.imshow(np.random.random((100,2000)), interpolation = 'nearest', cmap = 'Greys_r')

fig.subplots_adjust(hspace = 0, right = 1, left = 0.0, top = 0.99, bottom = 0.01)
savefig('../figures/jwst.png', dpi = 100)
plt.close('all')




