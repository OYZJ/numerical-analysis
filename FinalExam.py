
# coding: utf-8

# author:Zhangjian Ouyang
# date:12/13/2017

# In[1]:


import matplotlib
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.mlab as mlab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt


# # problem1

# In[2]:


def FD(t):
    a=0
    b=10
    deltay=0.05
    deltat=0.005
    v=1
    h=990
    den=1000
    cp=4180
    d=0.07
    def f(x):
        return 3+10*np.e**(-1/(2*x))
    N=int((b-a)/deltay)
    n=int(t/deltat)
    R=np.zeros((N+1,n+1))
    for i in range(1,N+1):
        R[i,0] = f(i*deltay)
    R[0,:] = 3
    for j in range(1,n+1):
        for i in range(1,N+1):
            R[i,j] = (R[i,j-1]/deltat+v*R[i-1,j]/deltay+4*h*R[i,0]/(den*cp*d))/(1/deltat+v/deltay+4*h/(den*cp*d))
                                                                                       
    y = np.linspace(a,b,N+1)                                                                                       
    plt.plot(y,R[:,n],'green')
    xmajorLocator = MultipleLocator(1)
    patch = mpatches.Patch(color='blue',label='t={0}'.format(t))
    plt.legend(handles=[patch])
    plt.grid(True)
    plt.show()


# In[3]:


for i in range(0,11):
    FD(i)


# In[ ]:





# # problem2

# In[12]:


delt=0.0005
delx=0.01

def FTCS(t):    
    N = int(1/delx)
    n = int(t/delt)
    I = np.eye(N-1)
    B = np.zeros(N-1)
    def f(x):
        return 1
    def ua(t):
        return 0
    def ub(t):
        return 1
    B[0] = -ua(t)
    B[N-2] = -ub(t)
    A = np.zeros((N-1,N-1))
    A[0,0] = A[N-2,N-2] = 2
    A[0,1] = A[N-2,N-3] = -1
    for i in range(1,N-2):
        A[i,i] = 2
        A[i,i-1] = A[i,i+1] = -1
    lam = 0.05*delt/delx**2
    w = np.zeros((N-1,n+1))
    x = np.zeros(N+1)
    x[0] = 0.001
    for i in range(1,N+1):
        x[i] = x[i-1] + delx
    for i in range(0,N-1):
        w[i,0] = f(x[i+1])
    for j in range(1,n+1):
        t_temp = (j-1)*delt
        B[0] = -ua(t_temp)
        B[N-2] = -ub(t_temp)
        w[:,j] = (I-lam*A).dot(w[:,j-1])-lam*B
    w1 = np.zeros(N+1)
    w1[0] = ua(t)
    w1[N] = ub(t)
    for i in range(1,N):
        w1[i] = w[i-1,n]
    return w1


# In[13]:


N = int(1/delx)
n = int(1/delt)
R = np.zeros((N+1,n+1))
for i in range(0,n+1):
    t= delt*i
    R[:,i] = FTCS(t)
    
x= np.linspace(0,1,N+1)
t= np.linspace(0,1,n+1)
plt.contourf(t, x, R, 10, alpha = 0.75, cmap = plt.cm.hot)
C = plt.contour(t, x, R, 10, colors = 'blue', linewidth = 1)
plt.clabel(C, inline = True, fontsize = 10)
plt.show()


# In[19]:


def Sl(x,y):
    z=y
    for i in range(1,100):
        z=z+2*np.e**(-0.05*i**2*np.pi**2*x)*np.sin(i*np.pi*y)/(i*np.pi)
    return z


# In[22]:


def FTCS(t):
    def f(x):
        return 1
    def ua(t):
        return 0
    def ub(t):
        return 1
    ht=0.0005
    hx=0.01
    N = int(1/hx)
    n = int(t/ht)
    I = np.eye(N-1)
    B = np.zeros(N-1)
    B[0] = -ua(t)
    B[N-2] = -ub(t)
    A = np.zeros((N-1,N-1))
    A[0,0] = A[N-2,N-2] = 2
    A[0,1] = A[N-2,N-3] = -1
    for i in range(1,N-2):
        A[i,i] = 2
        A[i,i-1] = A[i,i+1] = -1
    lam = 0.05*ht/hx**2
    w = np.zeros((N-1,n+1))
    x = np.zeros(N+1)
    x[0] = 0.001
    for i in range(1,N+1):
        x[i] = x[i-1] + hx
    for i in range(0,N-1):
        w[i,0] = f(x[i+1])
    for j in range(1,n+1):
        t_temp = (j-1)*ht
        B[0] = -ua(t_temp)
        B[N-2] = -ub(t_temp)
        w[:,j] = (I-lam*A).dot(w[:,j-1])-lam*B
   
    w1 = np.zeros(N+1)
    w2 = np.zeros(N+1)
    w1[0] = w2[0] = ua(t)
    w1[N] = w2[N] = ub(t)
    for i in range(1,N):
        w1[i] = w[i-1,n]
        w2[i] = Sl(t,hx*i)
    
    plt.plot(x,w1,'b')
    patch1 = mpatches.Patch(color='blue',label='Approximate Solution x={0}'.format(t))
    plt.legend(handles=[patch1])
    plt.xlabel('y')
    plt.ylabel('u')
    plt.grid(True)
    plt.show()
                                  
    plt.plot(x,w2,'r')
    patch2 = mpatches.Patch(color='red',label='Exact Solution x={0}'.format(t))
    plt.legend(handles=[patch2])
    plt.xlabel('y')
    plt.ylabel('u')
    plt.grid(True)
    plt.show()
                                  
                                  
FTCS(0.1)


# In[ ]:




