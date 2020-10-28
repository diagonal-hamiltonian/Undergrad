#!/usr/bin/env python
# coding: utf-8

# In[10]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.stats import linregress
import scipy.optimize as optimization


# In[11]:


class Least:
    '''estimates pi by comparing the number of random points inside a unit circle and unit square'''
    
    def __init__(self,time,decays):
         
        self.time = time     
        self.decays = decays
        

    def print_info(self):
        '''Output information about current estimate'''
        print('Number of decays', (len(self.decays)))
        
    def create_lin_lin_fig(self):
        
        fig, ax = plt.subplots(figsize=(7,7))
        ax.plot(self.time, self.decays, 'bo', ms='4', markerfacecolor='none',label='data points')
        ax.set_title('Decays, against time', fontsize=14)
        ax.set_xlabel('t', fontsize=14)
        ax.set_ylabel('Decays per second', fontsize=14)

        return fig, ax


    def create_lin_log_fig(self):

        fig, ax = plt.subplots(figsize=(7,7))
        ax.semilogy(self.time, self.decays, 'bo', ms='4', markerfacecolor='none',label='data points')
        ax.set_xlabel('t', fontsize=14)
        ax.set_ylabel('Ln(Decays per second)', fontsize=14)

        return fig, ax
    
    def parameters(self):
        N=len(self.decays)
        logdecays=np.log(self.decays)
        s1=sum(self.time)
        s2=sum((self.time)**2)
        s3=sum((logdecays))
        s4=sum((logdecays*self.time))
        
        a = (-s2 * s3 + s1 * s4) / (s1**2 - N * s2)
        b = (s1 * s3 - N * s4) / (s1**2 - N * s2)
        
        return a , b


# In[14]:
time = np.array([5,15,25,35,45,55,65,75,85,95,105,115])
decays = np.array([32,17,21,7,8,6,5,3,4,1,5,1])


fig, ax = plt.subplots(figsize=(7,7))
#ax.plot(time, np.exp(fit) , 'b--o', ms='4', markerfacecolor='none',label='Least square fit') #used exp() to invert log
ax.plot(time, decays , 'go' , ms='6', label= 'Data points')
#ax.plot(time, func(time,par1,par2),color='r',label='Scipy fit') #used exp() to invert log


ax.set_xlabel('t', fontsize=14)
ax.set_ylabel('Decays per second', fontsize=14)
ax.set_xlim(0,120)
ax.set_ylim(0,35)
plt.legend()
plt.show()


# In[13]:


data = Least(time,decays)

data.print_info()
fig, ax = data.create_lin_log_fig()
a,b = data.parameters()
print(np.exp(a), -b)


plt.legend()
plt.show()


# In[ ]:





# In[37]:


fig, ax = plt.subplots(figsize=(7,7))
ldec=np.log(decays)
fit=a + b*time 

ax.plot(time, fit , 'b--o', ms='4', markerfacecolor='none',label='Least square fit') #used exp() to invert log
ax.plot(time, ldec , 'gd' , label= 'Data points')

ax.set_xlabel('t', fontsize=14)
ax.set_ylabel('Ln(Decays per second)', fontsize=14)
plt.legend()
plt.show()


# In[15]:



print('slope=',b, a)


# In[16]:



fig, ax = plt.subplots(figsize=(7,7))
ldec=np.log(decays)
fit=a + b*time 

ax.plot(time, np.exp(fit) , 'b--o', ms='4', markerfacecolor='none',label='Least square fit') #used exp() to invert log
ax.plot(time, np.exp(ldec) , 'gd' , label= 'Data points')

ax.set_title('Linear fit ', fontsize=14)
ax.set_xlabel('t', fontsize=14)
ax.set_ylabel('Decays per second', fontsize=14)
plt.show()


# In[37]:



xdata = time
ydata = decays

x0    = np.array([0.1,0.1])

def func(x, a1, b1 ):
    return np.exp(a1+b1*x)

fitparameter=optimization.curve_fit(func, xdata, ydata, x0)[0]

par1=fitparameter[0]
par2=fitparameter[1]

print(par1,par2)


# In[18]:


fig, ax = plt.subplots(figsize=(7,7))
ax.plot(time, np.exp(fit) , 'b--o', ms='6',label='Least square fit') #used exp() to invert log
ax.plot(time, decays , 'go' , ms='6', label= 'Data points')
ax.plot(time, func(time,par1,par2),'r--o', ms='6',label='Scipy fit') #used exp() to invert log


ax.set_xlabel('t', fontsize=14)
ax.set_ylabel('Decays per second', fontsize=14)
ax.set_xlim(0,120)
ax.set_ylim(0,35)
plt.legend()
plt.show()


# In[40]:


print((par1-a)/a, (par2-b)/b)


# In[1]:


print('Aaron Miller\n17322856')


# In[ ]:




