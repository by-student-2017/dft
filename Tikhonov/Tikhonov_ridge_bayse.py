#-------------------------------------------------
# pip install scikit-optimize
import numpy as np
from skopt import gp_minimize

class RidgeReg:

    def __init__(self, lambda_ = 1.0):
        self.lambda_ = lambda_
        self.coef_ = None
        self.intercept_ = None

    def fit(self, X, y, L):
        #Add a column of all ones to the first row of the matrix of explanatory variables to include the calculation of the intercept
        #X = np.insert(X, 0, 1, axis=1)
        #Create Identity Matrix
        #i = np.eye(X.shape[1])
        #formula for weight
        #temp = np.linalg.inv(X.T @ X + self.lambda_ * i) @ X.T @ y
        temp = np.linalg.inv(X.T @ X + self.lambda_ * L.T @ L) @ X.T @ y
        #value of the regression coefficient
        #self.coef_ = temp[1:]
        self.coef_ = temp[0:]
        #intercept value
        #self.intercept_ = temp[0]

    def predict(self, X):
        #Returns predicted values from a ridge regression model
        #return (X @ self.coef_ + self.intercept_)
        return (X @ self.coef_)

def Sigma1(lam,U,S,y):
	s1 = 0.0
	for idx in range((len(S))):
		s1 = s1 + S[idx]**2*(S[idx]**2-2.0*lam) / (S[idx]**2+lam)**4 * np.dot(U,y)**2
	return s1

#--------------------------------------------------------------------
import pandas as pd

dfy = pd.read_csv('case.csv')
pp0 = dfy['PP0'].values
if dfy.columns[1] == "gg":
	print("convert g/g to mmol/g for H2")
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
U0, S0, V0 = np.linalg.svd(X, full_matrices=True)
idx = S0.argsort()[::-1]
S = S0[idx]
U = U0[idx]
V = V0[idx]
print("2nd order square matrix with left singular value vectors arranged as columns, U")
print(U)
print("Quadratic square matrix filled with 0 except for singular values, Sigma")
print(S)
print(np.diag(S))
print("quadratic square matrix with right singular value vectors arranged as columns, V")
print(V)
#
#print( Sigma1(150,U,S,y) )
#
#
#Check Singular Value Decomposition (U x Sigma x V, Matrix A is calculated)
#print("Eigenvalue decomposition Check A = U x Sigma x V")
#print(np.dot(np.dot(U,np.diag(S)),V))
#------------------
#--------------------------------------------------------------------

ndata = 250
Ly = np.zeros(ndata)
Lx = np.zeros(ndata)
Lam = np.zeros(ndata)
kappa = np.zeros(ndata)

iLy = 0.0
iLx = 0.0
Lam_x0 = 0.0

#L = np.eye(X.shape[1])
#print(np.diag(L))
L1 = np.ones(X.shape[1])
L = np.diag(L1)

def objective(param):

	#L1 = np.ones(X.shape[1])
	L1 = param
	L = np.diag(L1)
	
	kappa_max = 0.0
	chk = 0
	iLam_x0 = 0
	for i in range(ndata):
		Lam[i] = i*2+100
		#
		# self-made
		linear = RidgeReg(lambda_ = Lam[i])
		#
		# sklearn
		#linear = Ridge(alpha=Lam[i])
		
		linear.fit(X, y, L)
		#print(linear.coef_)
		xlamn1 = np.linalg.norm(linear.coef_, ord=1)
		rlamn1 = np.linalg.norm((y - X @ linear.coef_), ord=1)
		Ly[i] = np.log10(xlamn1)
		Lx[i] = np.log10(rlamn1)
		#
		chs = S[idx]/(S[idx]**2+Lam[i])*np.dot(U,y)*V.T
		if np.max(chs) > 1:
			print("  This system is unstable. You must increase lambda. lambda=", Lam[i])
		#
		#print(S)
		#print(S**2)
		sig1 = np.sum( (S**2 * (S**2-2.0*Lam[i]))/((S**2+Lam[i])**4) * np.dot(U,y)**2 )
		sig2 = np.sum( (S**2)/((S**2+Lam[i])**4) * np.dot(U,y)**2 )
		sig3 = np.sum( (S**2)/((S**2+Lam[i])**3) * np.dot(U,y)**2 )
		xlamn2 = np.linalg.norm(linear.coef_, ord=1)**2
		rlamn2 = np.linalg.norm((y - X @ linear.coef_), ord=1)**2
		kappa[i] = (1.0/(rlamn2 + Lam[i]**2*xlamn2)**(3/2))*abs( rlamn2*xlamn2*(sig1+3.0*Lam[i]*sig2)/(sig3**2) - Lam[i]*(Lam[i]*rlamn2 + xlamn2) )
		#
		print(Lam[i], Lx[i], Ly[i], kappa[i])
		#
		if np.min(linear.coef_)>0 and kappa[i]>kappa_max and chk<=1:
			iLam_x0 = i
			iLx = Lx[i]
			iLy = Ly[i]
			kappa_max = kappa[i]
			chk = 1
		#if np.min(linear.coef_)>0 and kappa[i]<=kappa_max and chk==1:
		#	chk = 2

	Lam_x0 = Lam[iLam_x0]
	
	# self-made
	linear = RidgeReg(lambda_ = Lam_x0)

	linear.fit(X, y, L)
	#print(linear.coef_)

	ret = np.sum(linear.coef_)
	return ret
    
#------------------
space = []
for i in range(len(L1)):
	if i<=17:
		space.append( (0.8, 1.2) )
	else:
		space.append( (1.0, 3.5) )

# scikit-optimize 0.8.1
res = gp_minimize(objective, 
    space,
    acq_func="EI",      # the acquisition function
    n_calls=200,         # the number of evaluations of f
    n_random_starts=5,  # the number of random initialization points
    noise=0.1**2,       # the noise level (optional)
    random_state=1234)  # the random seed

print(res.fun) 
print(res.x)   

L = np.diag(res.x)
#------------------

iLam_x0 = 0
kappa_max = 0
chk = 0
for i in range(ndata):
	Lam[i] = i*2+100
	#
	# self-made
	linear = RidgeReg(lambda_ = Lam[i])
	#
	# sklearn
	#linear = Ridge(alpha=Lam[i])
	
	linear.fit(X, y, L)
	#print(linear.coef_)
	xlamn1 = np.linalg.norm(linear.coef_, ord=1)
	rlamn1 = np.linalg.norm((y - X @ linear.coef_), ord=1)
	Ly[i] = np.log10(xlamn1)
	Lx[i] = np.log10(rlamn1)
	#
	chs = S[idx]/(S[idx]**2+Lam[i])*np.dot(U,y)*V.T
	if np.max(chs) > 1:
		print("  This system is unstable. You must increase lambda. lambda=", Lam[i])
	#
	#print(S)
	#print(S**2)
	sig1 = np.sum( (S**2 * (S**2-2.0*Lam[i]))/((S**2+Lam[i])**4) * np.dot(U,y)**2 )
	sig2 = np.sum( (S**2)/((S**2+Lam[i])**4) * np.dot(U,y)**2 )
	sig3 = np.sum( (S**2)/((S**2+Lam[i])**3) * np.dot(U,y)**2 )
	xlamn2 = np.linalg.norm(linear.coef_, ord=1)**2
	rlamn2 = np.linalg.norm((y - X @ linear.coef_), ord=1)**2
	kappa[i] = (1.0/(rlamn2 + Lam[i]**2*xlamn2)**(3/2))*abs( rlamn2*xlamn2*(sig1+3.0*Lam[i]*sig2)/(sig3**2) - Lam[i]*(Lam[i]*rlamn2 + xlamn2) )
	#
	print(Lam[i], Lx[i], Ly[i], kappa[i])
	#
	if np.min(linear.coef_)>0 and kappa[i]>kappa_max and chk<=1:
		iLam_x0 = i
		iLx = Lx[i]
		iLy = Ly[i]
		kappa_max = kappa[i]
		chk = 1
	if np.min(linear.coef_)>0 and kappa[i]<=kappa_max and chk==1:
		chk = 2

Lam_x0 = Lam[iLam_x0]
#--------------------------------------------------------------------
import matplotlib.pyplot as plt
#------------------
#plt.plot(Lx, Ly, "--", color="blue", label="Series")
#plt.plot(Lx, kappa, "--", color="green", label="curvature")
#plt.plot(iLx, iLy, "o", c="r", label="Fit(lambda={0})".format(Lam_x0))
#plt.legend()
#plt.xlabel("log10||Ax-y||", fontsize=10)
#plt.ylabel("log10||x||", fontsize=10)
#plt.show()
#
fig = plt.figure()
ax1 = fig.subplots()
ax2 = ax1.twinx()
ax1.plot(Lx, Ly, "--", color="blue", label="Series")
ax1.plot(iLx, iLy, "o", c="r", label="Fit(lambda={0})".format(Lam_x0))
ax2.plot(Lx, kappa, "--", color="green", label="curvature")
ax1.set_xlabel("log10||Ax-y||", fontsize=10)
ax1.set_ylabel("log10||x||", fontsize=10)
ax2.set_ylabel("curvature", fontsize=10)
plt.title("Title")
plt.legend()
#plt.xscale('log')
plt.show()
#------------------
# self-made
linear = RidgeReg(lambda_ = Lam_x0)
#
# sklearn
#linear = Ridge(alpha = Lam_x0)

linear.fit(X, y, L)
#print(linear.coef_)
#------------------
wdata = len(x[0])-1
dvdw = np.zeros(wdata)
a = np.zeros(wdata)
i = 0
for ww in dfx[:0]:
	if i>=1:
		a[i-1] = float(ww)
		dvdw[i-1] = linear.coef_[i-1] / a[i-1]
	i = i + 1
#
cw = np.zeros(wdata)
cw[0] = 0.0
for i in range(1,wdata):
	cw[i] = cw[i-1] + linear.coef_[i]
#------------------
plt.plot(pp0,y, color="blue", label="Exp.")
plt.plot(pp0,(X @ linear.coef_), c="r", label="Fit(lambda={0})".format(Lam_x0))
plt.legend()
plt.xlabel("Relative pressure, P/P0", fontsize=10)
plt.ylabel("Adsorption, [mmol/g]", fontsize=10)
plt.xscale('log')
plt.show()
#------------------
fig = plt.figure()
ax1 = fig.subplots()
ax2 = ax1.twinx()
ax1.plot(a,linear.coef_, color="blue", label="dVdw")
#ax1.plot(a,dvdw, color="blue", label="dVdw")
ax2.plot(a, cw, c="r", label="cum")
ax1.set_xlabel("Pore width, w [nm]", fontsize=10)
#ax1.set_ylabel("Distribution, dV/dw [cm3/g]", fontsize=10)
ax1.set_ylabel("Distribution, dV/dw [cm3/g/nm]", fontsize=10)
ax2.set_ylabel("Cumulative Volume, cm3/g", fontsize=10)
plt.hlines(y=0, xmin=a[0], xmax=a[wdata-1])
plt.title("Title")
plt.legend()
#plt.xscale('log')
plt.show()
#--------------------------------------------------------------------