#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 29 17:43:23 2019

@author: macuser
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import plotly.graph_objects as go
from scipy import stats

def function(x):
    return np.exp(x)+x*np.exp(x)

def analytic(x):
    return 2*np.exp(x)+x*np.exp(x)


def abs_error_num1(x,h):
    return abs(D1(x,h)-analytic(x))

def rel_error_num1(x,h):
    return abs(abs_error_num1(x,h)/analytic(x))


def abs_error_num2(x,h):
    return abs(D2(x,h)-analytic(x))

def rel_error_num2(x,h):
    return abs(abs_error_num2(x,h)/analytic(x))


def abs_error_num3(x,h):
    return abs(D3(x,h)-analytic(x))

def rel_error_num3(x,h):
    return abs(abs_error_num3(x,h)/analytic(x))


def abs_error_num4(x,h):
    return abs(D4(x,h)-analytic(x))

def rel_error_num4(x,h):
    return abs(abs_error_num4(x,h)/analytic(x))


def D1(x,h):
    return (function(x+h)-function(x-h))/(2*h)

def D2(x,h):
    return (4*D1(x,h)-D1(x,2*h))/3

def D3(x,h):
    return (16*D2(x,h)-D2(x,2*h))/15

def D4(x,h):
    return (64*D3(x,h)-D3(x,2*h))/63


#plot of f

fig1 = plt.figure()

ax1 = fig1.add_subplot(111)                                    #111 = 1 row, 1 column, 1st axes
x1 = np.linspace(0,2,20)
line_fn, = ax1.plot(x1, function(x1),'g-o')  
line_fn.set_markersize(3)
line_fn.set_label('$f(x)$')
line_ddfn, = ax1.plot(x1, analytic(x1),'b-o')  
line_ddfn.set_markersize(3)
line_ddfn.set_label('$f^{\prime\prime}(x)$')
ax1.set_xlabel('x', fontsize=16)
ax1.set_ylabel('y', fontsize=16)
ax1.set_xticks(np.linspace(0, 2.0, 5))
ax1.legend()
plt.show()


#table of values 

N=6
h = np.array([0.05,0.1,0.2,0.4])
d1,d2,d3,d4,=np.zeros(len(h)),np.zeros(len(h)),np.zeros(len(h)),np.zeros(len(h))
d1=D1(2,h)
d2=D2(2,h)
d3=D3(2,h)
d4=D4(2,h)


fig = go.Figure(data=[go.Table(header=dict(values=['h', 'D_1', 'D_2', 'D_3', 'D_4']),
                 cells=dict(values=[h,d1, d2, d3, d4]))
                     ])
fig.write_html('richardson_table.html', auto_open=True)

#highest accuracy of h D4finder

p=2
k=0.4
while rel_error_num4(p,k)>=rel_error_num4(p,k-0.001):
    k=k-0.001
print ('\nMax relative error =',rel_error_num4(p,k),11,'\nh =',k,'\n2nd derivative approximation x=2 is ',D4(p,k))

h=np.linspace(0.001,0.039,10000)
plt.plot(h,rel_error_num4(2,h),'ko', markersize=0.9)
plt.plot(0.028999999999999693,rel_error_num4(2,0.028999999999999693),'ro', markersize=6,label='Derived point')
plt.xlabel('h')
plt.ylabel('Relative error')
plt.legend()
plt.show()
#graph of scaling relationships

x=2.0
N=1000
start=0.028
stop=0.4
h=np.linspace(start,stop,N)
plt.loglog(h,abs_error_num1(2,h),'rd', markersize=0.9,label='D1')
plt.loglog(h,abs_error_num2(2,h),'ko', markersize=0.9,label='D2')
plt.loglog(h,abs_error_num3(2,h),'gd', markersize=0.9,label='D3')
plt.loglog(h,abs_error_num4(2,h),'yd', markersize=0.9,label='D4')
plt.ylabel('log(absolute error)')
plt.xlabel('log(h)')
plt.legend()
plt.show()

N=100000
start=0.000001
stop=0.028
h=np.linspace(start,stop,N)
plt.loglog(h,abs_error_num4(2,h),'yd', markersize=0.9,label='D4')
plt.loglog(h,abs_error_num3(2,h),'gd', markersize=0.9,label='D3')
plt.loglog(h,abs_error_num2(2,h),'ko', markersize=0.9,label='D2')
plt.loglog(h,abs_error_num1(2,h),'rd', markersize=0.9,label='D1')
plt.ylabel('log(absolute error)')
plt.xlabel('log(h)')
plt.legend()
plt.show()

h=np.linspace(0.028,0.4,1000)
H=np.log(h)
y = np.log(abs_error_num1(2,h))
slope, intercept, r_value, p_value, std_err = stats.linregress(H,y)

print('D1 slope = ',round(slope,4))

h=np.linspace(0.028,0.4,1000)
H=np.log(h)
y = np.log(abs_error_num2(2,h))
slope, intercept, r_value, p_value, std_err = stats.linregress(H,y)

print('D2 slope = ',round(slope,4))

h=np.linspace(0.028,0.4,1000)
H=np.log(h)
y = np.log(abs_error_num3(2,h))
slope, intercept, r_value, p_value, std_err = stats.linregress(H,y)

print('D3 slope = ',round(slope,4))

h=np.linspace(0.028,0.4,1000)
H=np.log(h)
y = np.log(abs_error_num4(2,h))
slope, intercept, r_value, p_value, std_err = stats.linregress(H,y)

print('D4 slope = ',round(slope,4))

h=np.linspace(0.000001,0.00028,1000)
H=np.log(h)
y = np.log(abs_error_num4(2,h))
slope, intercept, r_value, p_value, std_err = stats.linregress(H,y)

print('Rounding error slope = ',round(slope,4))

x=2.0
N=1000
start=0.028
stop=0.4

plt.show()
print('\n\n\n\nAaron Miller, 17322865')



