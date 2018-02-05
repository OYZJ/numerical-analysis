
# coding: utf-8

# # problem1 

# In[100]:


import matplotlib.pyplot as plt
import numpy as np
import sympy
import matplotlib.patches as mpatches


# In[106]:


def Lm(x,n,X):
    m=len(x)
    a=1
    b=1
    for i in range(0,m):
        if i==n:
            a=a
            b=b
        else:
            ua=a*(X-x[i])
            b=b*(x[n]-x[i])
    return a/b


# In[108]:


x=np.array([-3.0,0.0,np.e,np.pi])
X=sympy.Symbol('X')
l30=Lm(x,0,X)
l31=Lm(x,1,X)
l32=Lm(x,2,X)
l33=Lm(x,3,X)
print 'L3,0(X)=',l30
print 'L3,1(X)=',l31
print 'L3,2(X)=',l32
print 'L3,3(X)=',l33


# In[105]:


x=np.array([-3.0,0.0,np.e,np.pi])


# In[103]:


def L30(x):
    return (x*(x-np.e)*(x-pi))/(-3*(-3-np.e)*(-3-pi))
def L31(x):
    return ((x+3)*(x-np.e)*(x-pi))/(3*(-np.e)*(-pi))
def L32(x):
    return ((x+3)*x*(x-pi))/((3+np.e)*(np.e)*(-pi+np.e))
def L33(x):
    return ((x+3)*x*(x-np.e))/((pi+3)*pi*(pi-np.e))


# In[104]:


pi = np.pi
t = np.linspace(-3.0,pi,101)
r0 = np.zeros(101)
r1 = np.zeros(101)
r2 = np.zeros(101)
r3 = np.zeros(101)
r0 = L30(t)
r1 = L31(t)
r2 = L32(t)
r3 = L33(t)
plt.plot(t,r0,'r',t,r1,'b',t,r2,'g',t,r3,'y')
red_patch = mpatches.Patch(color='red',label='L30')
blue_patch = mpatches.Patch(color='blue',label='L31')
green_patch = mpatches.Patch(color='green',label='L32')
yellow_patch = mpatches.Patch(color='yellow',label='L33')
plt.legend(handles=[red_patch,blue_patch,green_patch,yellow_patch])
plt.xlabel('x')
plt.ylabel('L(x)')
plt.show()


# In[ ]:





# # problem2

# In[12]:


import matplotlib.pyplot as plt
import numpy as np
import sympy


# In[13]:


x=np.array([1.0,2.0,3.0])
y=np.array([np.log(1.0),np.log(2),np.log(3)])


# In[14]:


def L20(t):
    return (t-x[1])*(t-x[2])/((x[0]-x[1])*(x[0]-x[2]))
def L21(t):
    return (t-x[0])*(t-x[2])/((x[1]-x[0])*(x[1]-x[2]))
def L22(t):
    return (t-x[0])*(t-x[1])/((x[2]-x[0])*(x[2]-x[1]))
def P2(t):
    return y[0]*L20(t)+y[1]*L21(t)+y[2]*L22(t)


# In[15]:


X=sympy.Symbol('x')
P=P2(X)
print(P)


# In[22]:


a=np.arange(1.0,3.0,0.1)
p=P2(a)
plt.plot(a,np.log(a), label='y=ln(x)')
plt.plot(a,p,label='Lagrange Interpolation')
plt.legend()
plt.show()


# In[23]:


f1=P2(1.5)
error1=np.abs(f1-np.log(1.5))
print error1


# In[24]:


# so the error of ln(1.5) is 0.0229


# In[25]:


f2=P2(2.4)
error2=np.abs(f2-np.log(2.4))
print error2


# In[26]:


# so the error of ln(2.4) is 0.0144


# In[27]:


def d4f(t):
    return -6/t**4
def error4(t,ksi):
    return np.abs(d4f(ksi)*(t-x[0])*(t-x[1])*(t-x[2])/(4*3*2*1))


# In[28]:


emin=error4(1.5,3)
print 'The minimum error is',emin
emax=error4(1.5,1)
print 'The maximum error is',emax
print 'The error bound is [',emin,',',emax,']'


# In[29]:


# so the error is in the bound.


# In[ ]:





# # problem3

# In[34]:


import matplotlib.pyplot as plt
import numpy as np


# In[55]:


def Neville(x,y,x0):
    n=len(x)
    p=np.zeros((n,n))
    p[:,0]=y
    for i in range(n-1,0,-1):
        for j in range(0,i):
            p[n-i+j,n-i]=((x0-x[j])*p[n-i+j,n-i-1]-(x0-x[n-i+j])*p[n-i+j-1,n-i-1])/(x[n-i+j]-x[j])
    val=p[n-1,n-1]
    return p,val


# In[56]:


x=np.array([1.0,4,16])
y=np.array([np.sqrt(1.0),np.sqrt(4),np.sqrt(16)])
p,val=Neville(x,y,9.)
print 'The matrix of P is:'
print p
print 'The value of the interpolating polynomial at x=9 that uses all data points is',val


# In[ ]:





# # problem4

# In[57]:


import matplotlib.pyplot as plt
import numpy as np


# In[58]:


def Neville(x,y,x0):
    n=len(x)
    p=np.zeros((n,n))
    p[:,0]=y
    for i in range(n-1,0,-1):
        for j in range(0,i):
            p[n-i+j,n-i]=((x0-x[j])*p[n-i+j,n-i-1]-(x0-x[n-i+j])*p[n-i+j-1,n-i-1])/(x[n-i+j]-x[j])
    val=p[n-1,n-1]
    return p,val


# In[59]:


x=np.array([0.005,0.010,0.020,0.050,0.100,0.200,0.500,1.000,2.000])
y=np.array([0.924,0.896,0.859,0.794,0.732,0.656,0.536,0.430,0.316])
p1,val1=Neville(x,y,0.032)
print 'The mean activity coefficient for a molality of 0.032 is',val1
p2,val2=Neville(x,y,1.682)
print 'The mean activity coefficient for a molality of 1.682 is',val2


# In[60]:


# It is obvious that the mean activity coefficient for a molality of 1.682 is inappropriate.
# It is out of the domain of Neville's algorithm.


# In[ ]:





# # problem5

# In[64]:


import matplotlib.pyplot as plt
import numpy as np
import sympy


# In[67]:


x=np.array([1.0,2,3])
y=np.array([np.cos(1),np.cos(2),np.cos(3)])


# In[68]:


def Newton(x,y,x0):
    n=len(x)
    f=np.zeros((n,n))
    f[:,0]=y
    for i in range(n-1,0,-1):
        for j in range(0,i):
            f[n-i+j,n-i]=(f[n-i+j,n-i-1]-f[n-i+j-1,n-i-1])/(x[n-i+j]-x[j])
    P=0
    for i in range(0,n):
        a=f[i,i]
        for j in range(0,i):
            a=a*(x0-x[j])
        P=P+a
    return P


# In[69]:


X=sympy.Symbol('X')
P=Newton(x,y,X)
print P


# In[ ]:





# # problem6

# In[70]:


import matplotlib.pyplot as plt
import numpy as np
import sympy


# In[71]:


x=np.array([10.,25.,50.,75.,100.])
y=np.array([488.55,485.48,480.36,475.23,470.11])


# In[73]:


def Newton(x,y,x0):
    n=len(x)
    f=np.zeros((n,n))
    f[:,0]=y
    for i in range(n-1,0,-1):
        for j in range(0,i):
            f[n-i+j,n-i]=(f[n-i+j,n-i-1]-f[n-i+j-1,n-i-1])/(x[n-i+j]-x[j])
    P=0
    for i in range(0,n):
        a=f[i,i]
        for j in range(0,i):
            a=a*(x0-x[j])
        P=P+a
    return P


# In[75]:


X=sympy.Symbol('X')
P=Newton(x,y,X)
print 'The interpolating polynomial is:'
print sympy.expand(P)


# In[88]:


t=np.arange(5,101,5)
ST=np.zeros(len(t))
for i in range(0,20):
    ST[i]=P=Newton(x,y,t[i])
print 'The temperature is:'
print t
print 'Surface Tension(dyn/cm):'
print ST


# In[89]:


plt.plot(t,ST)
plt.plot(x,y,"*")
plt.grid(True)
plt.show()


# In[ ]:





# # problem7

# In[90]:


import matplotlib.pyplot as plt
import numpy as np


# In[91]:


x=np.array([100.,200.,300.,400.,500.,600.])
y=np.array([9.4,18.4,26.2,33.3,39.7,45.7])


# In[93]:


def Newton(x,y,x0):
    n=len(x)
    f=np.zeros((n,n))
    f[:,0]=y
    for i in range(n-1,0,-1):
        for j in range(0,i):
            f[n-i+j,n-i]=(f[n-i+j,n-i-1]-f[n-i+j-1,n-i-1])/(x[n-i+j]-x[j])
    P=0
    for i in range(0,n):
        a=f[i,i]
        for j in range(0,i):
            a=a*(x0-x[j])
        P=P+a
    return P


# In[95]:


P1=Newton(x,y,240)
P2=Newton(x,y,485)
print 'The thermal conductivity of air when T=240K is:',P1
print 'The thermal conductivity of air when T=485K is:',P2


# In[ ]:




