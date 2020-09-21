
# This was the test code for the noise realisations with QU covariances.
# n = 10000
# x, y = np.random.multivariate_normal([0,0], [[1, 0.9], [0.9, 1]], n).T
# plt.plot(x,y,'.')
# plt.savefig('noisetest2.png')
# plt.clf()
# exit()

# A = np.asarray([[1.0, 0.9], [0.9, 1.0]])
# B = np.asarray([A]*n)
# C = np.linalg.cholesky(B)
# vals = np.random.normal(0,1, size = (n,2))
# x,y = np.einsum('nij,nj->ni', C, vals).T
# plt.plot(x,y,'.')
# plt.savefig('noisetest.png')
# exit()

# C = precalc_C(np.asarray([1.0]*n), np.asarray([1.0]*n), np.asarray([0.9]*n))
# x,y = noiserealisation_QU(C,np.asarray([1.0]*n),np.asarray([1.0]*n))
# plt.plot(x,y,'.')
# plt.savefig('noisetest3.png')
# exit()