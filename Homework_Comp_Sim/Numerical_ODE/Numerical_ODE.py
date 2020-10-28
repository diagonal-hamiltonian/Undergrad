#!/usr/bin/env python
# coding: utf-8

# In[162]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.stats import linregress
import scipy.optimize as optimization
import itertools
from scipy.integrate import odeint


# In[163]:


def f(x,t):
    return (1 + t)*x + 1 - 3*t + t**2


# In[164]:


def Euler(dt,x0,i):
    x = x0 + f(x0 , i*dt)*dt
    return x

def im_Euler(dt,x0,i):
    x = x0 +( f(x0, (i)*dt ) + f(x0 +f(x0 ,(i)*dt )*dt, (i+1)*dt)) * dt/2
    return x


def rkk(dt,x0,i):

    k1 = dt * f( x0, i*dt) 
    k2 = dt* f( x0 + 0.5 * k1, i*dt + 0.5 *dt) 
    k3 =  dt* f( x0 + 0.5 * k2, i*dt + 0.5 * dt) 
    k4 = dt * f( x0 + k3, i*dt + dt) 
    x = x0 + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
    return x


# In[165]:


dt=0.04

x0=0.0655
tf=5.

t=np.linspace(0,5.,int(tf/dt))
sim_eu = np.zeros(int(tf/dt))
im_eu = np.zeros(int(tf/dt))
rk = np.zeros(int(tf/dt))


for i in range(0,int(tf/dt)):
    sim_eu[i]=Euler(dt,x0,i)
    x0=sim_eu[i]

x0=0.0655
for i in range(0,int(tf/dt)):
    im_eu[i]=im_Euler(dt,x0,i)
    x0=im_eu[i]

    
x0=0.0655
for i in range(0,int(tf/dt)):
    rk[i]=rkk(dt,x0,i)
    x0=rk[i]
    


# In[166]:


dt2=0.02

x0=0.0655
tf=5.

t2=np.linspace(0,5.,int(tf/dt2))
sim_eu2 = np.zeros(int(tf/dt2))
im_eu2 = np.zeros(int(tf/dt2))
rk2 = np.zeros(int(tf/dt2))


for i in range(0,int(tf/dt2)):
    sim_eu2[i]=Euler(dt2,x0,i)
    x0=sim_eu2[i]

x0=0.0655
for i in range(0,int(tf/dt2)):
    im_eu2[i]=im_Euler(dt2,x0,i)
    x0=im_eu2[i]

    
x0=0.0655
for i in range(0,int(tf/dt2)):
    rk2[i]=rkk(dt2,x0,i)
    x0=rk2[i]


# In[167]:


tt, xx = np.arange(0 ,5 ,0.2 ), np.arange(-3. ,3. ,0.24 )
T, X = np.meshgrid(tt, xx)
slope = f(X ,T )

tslopes=np.cos(np.arctan(slope))
xslopes=np.sin(np.arctan(slope))


fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111)

ax.set_xlabel('time')
ax.set_ylabel('x')

ax.quiver(T ,X ,tslopes ,xslopes)
ax.plot(t,sim_eu,label='Simple Euler')
ax.plot(t,im_eu,label='Improved Euler')
ax.plot(t,rk, label='Runge-Kutta')


ax.set_ylim(-3,3)
ax.set_xlim(0,5)
ax.legend()
plt.show()


# In[ ]:





# In[168]:


fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111)
ax.streamplot(T ,X ,tslopes ,xslopes, linewidth=0.5 ,color='k')
ax.plot(t,sim_eu,label='Simple Euler')
ax.plot(t,im_eu,label='Improved Euler')
ax.plot(t,rk, label='Runge-Kutta')
ax.set_xlabel('time')
ax.set_ylabel('x')
ax.set_ylim(-3,3)
ax.set_xlim(0,5)
ax.legend()


# In[169]:


fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111)
x0=0.0655
ax.plot(t,sim_eu,label='Simple Euler, $\Delta t=0.04$')
ax.plot(t2,sim_eu2,label='Simple Euler,  $\Delta t=0.02$')

ax.plot(t,im_eu,label='Improved Euler, $\Delta t=0.04$')
ax.plot(t2,im_eu2,label='Improved Euler,  $\Delta t=0.02$')

ax.plot(t,rk, label='Runge-Kutta, $\Delta t=0.04$')
ax.plot(t2,rk2, label='Runge-Kutta, $\Delta t=0.02$')
ax.plot(0,x0,'mo',label='x=0.0655', markersize=10)

ax.set_ylim(-3,3)
ax.grid()
ax.legend()
plt.show()

print('Aaron Miller, 17322856')
