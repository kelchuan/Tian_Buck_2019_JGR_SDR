import numpy as np

data_file = open("data.txt", "r")
data = data_file.readlines()

count_line = len(data)
print(count_line)
#for line in data_file:
#    count_line = count_line + 1
    #print(float(line.split()[0]))
    
x = np.zeros(count_line)
y = np.zeros(count_line)
z = np.zeros(count_line)

for i in range(count_line):
    #print(i)
    x[i] = float(data[i].split()[1])
    y[i] = float(data[i].split()[0])
    z[i] = float(data[i].split()[2])
num_profiles = 20
rows = int(count_line/num_profiles)
#X = np.reshape(x, (rows, num_profiles))
X = np.reshape(x, (num_profiles, rows))
X = np.transpose(X)
Y = np.reshape(y, (num_profiles, rows))
Y = np.transpose(Y)
Z = np.reshape(z, (num_profiles, rows))
Z = np.transpose(Z)
Z = -Z


import matplotlib.pyplot as plt
fig = plt.figure()
pcolor_test = plt.pcolormesh(X, Y, Z, cmap='jet')
fig.colorbar(pcolor_test, shrink=0.5, aspect=5)
#plt.imshow(Z, interpolation = 'bilinear')
plt.show()



from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
#from matplotlib.colors import LightSource

fig = plt.figure()
ax = fig.gca(projection='3d')

# create light source object.
#ls = LightSource(azdeg=160, altdeg=30)
#rgb = ls.shade(Z, plt.cm.jet)



#surf = ax.plot_surface(X, Y, Z, cmap=cm.jet, facecolors=rgb)
surf = ax.plot_surface(X, Y, Z, cmap=cm.jet)
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
