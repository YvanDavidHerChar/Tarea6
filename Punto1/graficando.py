import numpy as np
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages

datos = sys.argv[1]

data = np.loadtxt(datos)

fig1 = plt.figure()
plt.plot(data[:,1], data[:,2])

fig2 = plt.figure()
ax =fig2.add_subplot(1,2,1, projection = '3d')
plt.plot(data[:,1],data[:,2],data[:,3])

pp = PdfPages("DaddysHome.pdf")
pp.savefig(fig1)
pp.savefig(fig2)
pp.close()
