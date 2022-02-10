import numpy as np
import matplotlib.pyplot as plt
from smoothnoisemap import *

rescalefactor = 10.0

# This was the test code for the noise realisations with QU covariances.
n = 100000
x, y = np.random.multivariate_normal([0,0], [[1, 0.9], [0.9, 1]], n).T
x *= rescalefactor
y *= rescalefactor
print(np.std(x))
plt.plot(x,y,'.')
plt.savefig('smoothnoise_test1_'+str(rescalefactor)+'.png')
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

C = precalc_C(np.asarray([rescalefactor*1.0]*n), np.asarray([rescalefactor*1.0]*n), np.asarray([rescalefactor*0.9]*n))
x_new,y_new = noiserealisation_QU(C)
plt.plot(x_new,y_new,'.')
print(np.std(x_new))
print(np.std(y_new))
plt.savefig('smoothnoise_test2_'+str(rescalefactor)+'.png')

C = precalc_C(np.asarray([rescalefactor**2*1.0]*n), np.asarray([rescalefactor**2*1.0]*n), np.asarray([rescalefactor**2*0.9]*n))
x_new,y_new = noiserealisation_QU(C)
plt.plot(x_new,y_new,'.')
print(np.std(x_new))
print(np.std(y_new))
plt.savefig('smoothnoise_test3_'+str(rescalefactor)+'.png')

C = precalc_C(np.asarray([1.0]*n), np.asarray([1.0]*n), np.asarray([0.9]*n))
x,y = noiserealisation_QU(C)
x *= rescalefactor
y *= rescalefactor
plt.plot(x,y,'.')
plt.savefig('smoothnoise_test4_'+str(rescalefactor)+'.png')

exit()