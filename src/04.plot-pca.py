import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt

f=open("derived-genetic/grm-3011/pca.eigenvec")
l=f.readlines()
f.close()

X=[]
Y=[]
for i in range(len(l)):
    X.append(float(l[i].split(' ')[2]))
    Y.append(float(l[i].split(' ')[3]))
fig=plt.figure()
plt.plot(X,Y,'ro',alpha=0.1)
#plt.show()
plt.savefig('derived-genetic/grm-3011/pca.png')

