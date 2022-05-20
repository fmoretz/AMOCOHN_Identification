import numpy as np
import matplotlib.pyplot as plt

from Identification import*
from dataimport import*

plt.figure(1)
plt.title('Kinetic Parameters X1')
plt.plot(Xr1, Yr1,'o',label = 'Original Data')
plt.plot(Xr1, mdl1.intercept_ + mdl1.coef_[0]*X11 + mdl1.coef_[1]*X12, label = 'fitting')
plt.legend()

plt.figure(2)
plt.title('Kinetic Parameters X2')
plt.plot(func2(Xr2,beta2[0],beta2[1],beta2[2],beta2[3]), 'o',label = 'Estimated')
plt.plot(Yr2,'*', label = 'Real data')
plt.legend()

plt.figure(3)
plt.title('Hydrolyis')
plt.plot(X_hydr, Y_hydr, 'o', label='original data')
plt.plot(X_hydr, mdl_hyd.intercept_ + mdl_hyd.coef_*X_hydr, 'r', label='fitted line')
plt.legend()

plt.figure(4)
plt.title('L/G transfer')
plt.plot(X3r, Y3r, 'o', label = 'data')
plt.plot(X3r, mdl3.intercept_ + mdl3.coef_*X3r, 'r')
plt.legend()

plt.figure(5)
plt.title('Regression for k1')
plt.plot(X4r,Y4r,'o', label = 'data')
plt.plot(X4r, mdl4.intercept_ + mdl4.coef_*X4r, 'r', label = 'Fitting')

plt.figure(6)
plt.title('Regression for k6')
plt.plot(X5r,Y5r,'o', label = 'data')
plt.plot(X5r, mdl5.intercept_ + mdl5.coef_*X5r, 'r', label = 'Fitting')

plt.figure(7)
plt.suptitle('Regression for k6/k3 and k2/k1')

plt.subplot(1,2,1)
plt.plot(X61, Y6r,'o', label = 'data')
plt.plot(X61, mdl6.intercept_ + mdl6.coef_[0]*X61 +mdl6.coef_[1]*X62, 'r', label = 'Fitting')
plt.ylabel('k6/k3*(D*(S2,in -S2))+ k6/k3*k2/k1*(D*(S1,in -S1)+K_hyd*XT)')
plt.xlabel('(*D(S2,in -S2))')

plt.subplot(1,2,2)

plt.plot(X62, Y6r,'o', label = 'data')
plt.plot(X62, mdl6.intercept_ + mdl6.coef_[0]*X61 +mdl6.coef_[1]*X62, 'r', label = 'Fitting')
plt.xlabel('(D*(S1,in -S1)+K_hyd*XT)')


plt.figure(8)
plt.suptitle('Regression for k4/k1 and k5/k6')

plt.subplot(1,2,1)
plt.plot(X71,Y7r,'o', label = 'data')
plt.plot(X71, mdl7.intercept_ + mdl7.coef_[0]*X71 + mdl7.coef_[1]*X72, 'r', label = 'Fitting')
plt.ylabel('k4/k1*(D*(S1,in -S1)+K_hyd*XT) + k5/k6*q_M')
plt.xlabel('(D*(S1,in -S1)+K_hyd*XT)')

plt.subplot(1,2,2)

plt.plot(X72,Y7r,'o', label = 'data')
plt.plot(X72, mdl7.intercept_ + mdl7.coef_[0]*X71 + mdl7.coef_[1]*X72, 'r', label = 'Fitting')
plt.xlabel('qM')

print(mdl6.score(X6r,Y6r))
print(mdl7.score(X7r,Y7r))
plt.show()
