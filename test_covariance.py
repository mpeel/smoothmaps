#python3 test.py
# [[9.66517857e+00 2.23214286e-02]
#   [1.11607143e-02 1.43325893e+01]]
# [[0.33314732 0.02779018]
#   [0.02779018 0.33314732]]
# 9.666666666666666
# 14.333333333333334
# 0.027777777777777776

import numpy as np

# Check the inverse weights
a = np.array([[1., 2.], [3., 4.]])
ainv = np.linalg.inv(a)
print(ainv)
print(a[0,1])
det = a[0,0]*a[1,1] - a[0,1]*a[1,0]
print(det)
ainv2 =  [ [a[0,0]/det, -a[0,1]/det], [-a[1,0]/det, a[1,1]/det] ] 
print(ainv2)
exit()

# Main test, comparing mine with Alberto's
Q1 = 11.0
Q2 = 9.0
U1 = 15.0
U2 = 14.0
Q1var = 1.0
U1var = 1.0
QU1var = 0.05
Q2var = 0.5
U2var = 0.5
QU2var = 0.05

W1 = np.linalg.inv(np.asarray([[Q1var,QU1var],[QU1var,U1var]]).T)
D1 = np.asarray([[Q1,0.0],[0.0,U1]]).T
W2 = np.linalg.inv(np.asarray([[Q2var,QU2var],[QU2var,U2var]]).T)
D2 = np.asarray([[Q2,0.0],[0.0,U2]]).T

combine = D1 @ W1 + D2 @ W2
weight = W1 + W2

inv_W = np.linalg.inv(weight)
combine = combine @ inv_W

print(combine)
print(inv_W)

w_Q_corr = 1.0/Q1var
w_Q_uncorr = 1.0/Q2var
w_U_corr = 1.0/U1var
w_U_uncorr = 1.0/U2var


Q = ((Q1/Q1var) + (Q2/Q2var)) / ((1/Q1var) + (1/Q2var))
U = ((U1/U1var) + (U2/U2var)) / ((1/U1var) + (1/U2var))
print(Q)
print(U)

QU = 1/(w_Q_corr+w_Q_uncorr) * 1/(w_U_corr+w_U_uncorr) * (w_Q_corr * w_U_corr * QU1var + w_Q_uncorr * w_U_uncorr * QU2var  )
print(QU)
