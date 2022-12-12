#-------------------------------------------------
import numpy as np

class RidgeReg:

    def __init__(self, lambda_ = 1.0):
        self.lambda_ = lambda_
        self.coef_ = None
        self.intercept_ = None

    def fit(self, X, y):
        #Add a column of all ones to the first row of the matrix of explanatory variables to include the calculation of the intercept
        #X = np.insert(X, 0, 1, axis=1)
        #Create Identity Matrix
        i = np.eye(X.shape[1])
        #formula for weight
        temp = np.linalg.inv(X.T @ X + self.lambda_ * i) @ X.T @ y
        #value of the regression coefficient
        #self.coef_ = temp[1:]
        self.coef_ = temp[0:]
        #intercept value
        #self.intercept_ = temp[0]

    def predict(self, X):
        #Returns predicted values from a ridge regression model
        #return (X @ self.coef_ + self.intercept_)
        return (X @ self.coef_)

#--------------------------------------------------------------------
import pandas as pd

dfy = pd.read_csv('case.csv')
pp0 = dfy['PP0'].values
y = dfy['gg'].values / (2.0*1.00784) * 1000.0 # g/g to mmol/g for H2

dfx = pd.read_csv('kernel.csv')
x = dfx.values

X = np.zeros((len(pp0),len(x[0])-1))
SX = np.zeros(len(pp0))
RSX = np.zeros(len(pp0))
xm = []
i = 0
for p in pp0:
	flag = 0
	for xx in x:
		if xx[0]>p and flag==0:
			flag = 1
			x2 = xx
			xm = (x2-x1)/(xx[0]-xx0)*(p-xx0)+x1
			X[i] = xm[1:]
			i = i + 1
			break
		x1 = xx
		xx0 = xx[0]
#print(X)

#from sklearn.preprocessing import StandardScaler
#X = StandardScaler().fit_transform(X)

from sklearn.linear_model import Ridge

#------------------
ndata = 50
Ly = np.zeros(ndata)
Lx = np.zeros(ndata)
Lam = np.zeros(ndata)
for i in range(ndata):
	Lam[i] = i*1+4
	#
	# self-made
	#linear = RidgeReg(lambda_ = Lam[i])
	#
	# sklearn
	linear = Ridge(alpha=Lam[i])
	
	linear.fit(X, y)
	#print(linear.coef_)
	Ly[i] = np.log10( np.linalg.norm(linear.coef_, ord=1) )
	Lx[i] = np.log10( np.linalg.norm((X @ linear.coef_ - y), ord=1) )
	print(Lam[i], Lx[i], Ly[i])

Lam_x0 = 200
#--------------------------------------------------------------------
import matplotlib.pyplot as plt
plt.plot(Lx,Ly)
plt.show()

# self-made
#linear = RidgeReg(lambda_ = Lam_x0)
#
# sklearn
linear = Ridge(alpha = Lam_x0)

linear.fit(X, y)
#print(linear.coef_)

wdata = len(x[0])-1
a = np.zeros(wdata)
i = 0
for ww in dfx[:0]:
	if not i==0:
		a[i-1] = float(ww)
	i = i + 1
#
cw = np.zeros(wdata)
cw[0] = 0.0
for i in range(1,wdata):
	cw[i] = cw[i-1] + linear.coef_[i]

fig = plt.figure()

ax1 = fig.subplots()
ax2 = ax1.twinx()

ax1.plot(a,linear.coef_, color="blue", label="dVdw")
ax2.plot(a, cw, c="r", label="cum")
ax1.set_xlabel("Pore width, w [nm]", fontsize=10)
ax1.set_ylabel("Distribution, dV/dw [cm3/g]", fontsize=10)
ax2.set_ylabel("Volume, cm3/g", fontsize=10)
plt.hlines(y=0, xmin=a[0], xmax=a[wdata-1])
plt.title("Title")
plt.legend()
#plt.xscale('log')

plt.show()