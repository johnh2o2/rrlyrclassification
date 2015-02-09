import numpy as np
from sklearn.svm import SVC
from sklearn.lda import LDA
from sklearn.qda import QDA
from sklearn.metrics import roc_curve, auc, accuracy_score
from math import *
import matplotlib.pyplot as plt

model = SVC(kernel='rbf', probability=True)
#model = QDA()
sigs = np.array([ [ 0.1, 0.2 ], [ 0.3, 0.1 ] ])
mus = np.array([ [ 0.5, 0.5 ], [ -0.5, -0.5] ])

def rand(a, b, num=None):
	return (b - a)*np.random.sample(num) + a

def pfunc(x, y):
	P = 0.
	X = np.array([x,y])
	for i in range(len(sigs)):
		A = 0.
		for a, u, s in zip(X, mus[i], sigs[i]):
			A += (a - u)**2/(2*s**2)
		P += np.exp(-A)
	return min( [ P, 1. ])

def MC(num, prob, bounds = [ -1., -1., 1., 1. ]):
	Rx = lambda : rand( bounds[0], bounds[2] )
	Ry = lambda : rand( bounds[1], bounds[3] )
	R1 = lambda : rand( 0, 1. )
	points = []

	while len(points) < num:
		x, y, p = Rx(), Ry(), R1()
		#print x, y, p, prob(x, y)
		if p < prob(x, y): points.append([ x, y ])

	return points

positive_points = MC(200, pfunc)
negative_points = MC(200, lambda x, y : 1 - pfunc(x, y))

all_points = []
all_points.extend(positive_points)
all_points.extend(negative_points)

def lbl(i):
	if i < len(positive_points): return 1
	else: return 0
labels = [ lbl(i) for i in range(len(all_points)) ]
data = zip(labels, all_points)
np.random.shuffle(data)

labels = [ l for l, d in data ]
all_points = [ d for l, d in data ]

xdist = np.linspace(-1, 1)
ydist = np.linspace(-1, 1)
X, Y = np.meshgrid(xdist, ydist)
print "fitting model"
model.fit(all_points, labels)
Z = np.zeros((len(xdist), len(ydist)))
for i, x in enumerate(xdist):
	for j, y in enumerate(ydist):
		#print model.predict_proba([ x, y ])[0][0]
		Z[i][j] = model.predict_proba([ x, y ])[0][1]

pred = [ model.predict_proba(ap)[0][1] for ap in all_points ]
fpr,tpr,_ = roc_curve( labels,  pred )
		
#print len(fpr), len(tpr)
roc_auc = auc(fpr, tpr)
print roc_auc

px = [ p[0] for p in positive_points ]
py = [ p[1] for p in positive_points ]

nx = [ n[0] for n in negative_points ]
ny = [ n[1] for n in negative_points ]
f = plt.figure()
ax = f.add_subplot(111)
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
cf = ax.contourf(Y, X, Z, cmap=plt.get_cmap('OrRd'), alpha=0.5)
plt.colorbar(cf)
ax.scatter(px, py, facecolor='g', alpha = 0.8)
ax.scatter(nx, ny, facecolor='k', alpha = 0.8)
plt.show()


