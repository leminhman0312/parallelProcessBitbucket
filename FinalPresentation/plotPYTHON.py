import numpy as np 
import matplotlib.pyplot as plt 
import sod 

#LOAD NUMERICAL DATA
data = np.loadtxt('result.txt')
x  = data[:,0]
rho= data[:,1]
u= data[:,2]
p= data[:,3]
mach = data[:,4]
figwidth = 15
figheight = 10

#LOAD ANALYTICAL DATA
gamma = 1.4
npts = 1000
pleft = 1./(gamma-1)
pright = 0.125/(gamma-1.)
positions, regions, values = sod.solve(left_state=(pleft, 1, 0),
        right_state=(pright, 0.125, 0.),
        geometry=(-10., 10., 0.0), t=2.00, gamma=gamma, npts=npts)
p_ana = values['p']
rho_ana = values['rho']
u_ana = values['u']
speedsound_ana = np.sqrt((gamma*p_ana)/(rho_ana))
mach_ana = u_ana/speedsound_ana  

#PLOTTING 
plt.subplots(figsize=(figwidth,figheight))
plt.subplot(4,1,1)
plt.title("DENSITY",fontsize=20)
plt.plot(x,rho,'-xb',label='Numerical')
plt.plot(x,rho_ana,'-r',label='Analytical')
plt.xticks([], [])
plt.yticks(fontsize=14) 
plt.legend()
plt.subplot(4,1,2)
plt.title("VELOCITY",fontsize=20)
plt.plot(x,u,'-xb',label='Numerical')
plt.plot(x,u_ana,'-r',label='Analytical')
plt.xticks([], [])
plt.yticks(fontsize=14)
plt.legend()
plt.subplot(4,1,3)
plt.title("PRESSURE",fontsize=20)
plt.plot(x,p,'-xb',label='Numerical')
plt.plot(x,p_ana,'-r',label='Analytical')
plt.xticks([],[])
plt.yticks(fontsize=14)
plt.legend()
plt.subplot(4,1,4)
plt.title("MACH",fontsize=20)
plt.plot(x,mach,'-xb',label='Numerical')
plt.plot(x,mach_ana,'-r',label='Analytical')
plt.yticks(fontsize=14)
plt.xticks(fontsize=14) 
plt.legend()
plt.savefig("shock.png")
plt.show()


