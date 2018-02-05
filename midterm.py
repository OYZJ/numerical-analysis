
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt  
import numpy as np  


# In[2]:


## problem 1


# In[3]:


def f(x):
    return 2*x**4-12*x**2+16*x-6


# In[4]:


def df(x):
    return 8*x**3-24*x+16


# In[5]:


def g(x,f,df):
    return x-f(x)/df(x)


# In[45]:


Nmax=11
xn=np.zeros(Nmax)
xn[0]=0


# In[46]:


for i in range(1,Nmax):
    xn[i]=g(xn[i-1],f,df)


# In[47]:


print (xn)


# In[48]:


iternum = np.arange(0,11,1)
plt.plot(iternum,xn,"x" ) 
plt.grid(True) 
plt.show()


# In[49]:


er=np.zeros(Nmax)
for i in range(0,Nmax-1):
    er[i]=np.abs(xn[i+1]-xn[i])


# In[50]:


print(er)


# In[51]:


plt.plot(iternum,er,"X") 
plt.grid(True) 
plt.show() 


# In[52]:


ratio = np.zeros(Nmax)
for i in range(0,Nmax-1):
    ratio[i]=er[i+1]/er[i]
print(ratio)


# # So the order of convergence is linear

# In[54]:


Nmax=10
m=0.0
n=3.0
an=np.zeros(Nmax)
an[0]=1.5
for i in range(0,Nmax-1):
    if f(m)*f(an[i])<0:
        m=m
        n=an[i]
        an[i+1]=(m+n)/2
    else:
        m=an[i]
        n=n
        an[i+1]=(m+n)/2
    i=i+1
print an


# In[24]:


er1=np.zeros(Nmax)
for i in range(0,Nmax-1):
    er1[i]=np.abs(an[i+1]-an[i])
print(er1)


# In[25]:


ratio = np.zeros(Nmax)
for i in range(0,Nmax-1):
    ratio[i]=er1[i+1]/er1[i]
print(ratio)


# # So the order of convergence is linear. The same as the Newton Method

# In[28]:


Nmax=11
xn=np.zeros(Nmax)
xn[0]=-4


# In[29]:


for i in range(1,Nmax):
    xn[i]=g(xn[i-1],f,df)


# In[30]:


print (xn)


# In[31]:


iternum = np.arange(0,11,1)
plt.plot(iternum,xn,"x" ) 
plt.grid(True) 
plt.show()


# In[34]:


er=np.zeros(Nmax)
for i in range(0,Nmax-1):
    er[i]=np.abs(xn[i+1]-xn[i])


# In[36]:


ratio = np.zeros(Nmax)
for i in range(0,Nmax-1):
    if er[i]!=0:
        ratio[i]=er[i+1]/er[i]
    else:
        break
print(ratio)


# # So the order of convergence is quadratic, from the ratio we can learn the order of convergence.
# 

# In[ ]:





# In[38]:


## problem 2


# In[40]:


import numpy as np
import matplotlib.pyplot as plt


# In[63]:


A = np.array([[1.0,0.91,0.82,0.70,0.69,0.60,0.70,0.50],
              [0.91,1.0,0.85,0.80,0.77,0.60,0.69,0.60],
              [0.82,0.85,1.00,0.90,0.83,0.77,0.78,0.67],
              [0.70,0.80,0.90,1.00,0.97,0.85,0.87,0.79],
              [0.69,0.77,0.83,0.97,1.00,0.92,0.95,0.80],
              [0.60,0.60,0.77,0.85,0.92,1.00,0.97,0.92],
              [0.70,0.69,0.78,0.87,0.95,0.97,1.00,0.94],
              [0.50,0.60,0.67,0.79,0.80,0.92,0.94,1.00]])
x = [1.0,1,1,1,1,1,1,1]
n = 10


# In[64]:


def PowerMethod(A,x,n):
    m=0
    y=np.ones(len(x))
    lamda=np.ones(n)
    conv=np.ones(n-1)
    while m<n:
        y=A.dot(x)
        for i in range(0,len(y)):
            if np.abs(y[i])==max(np.abs(y)):
                lamda[m]=y[i]
                p=i
                break
        x=y/lamda[m]
        m=m+1
    for i in range(0,n-1):
        conv[i]=lamda[i+1]-lamda[i]
    return x,lamda,conv


# In[65]:


eigenvector,eigenvalue, conv = PowerMethod(A1,x1,n)
print eigenvector,eigenvalue, conv


# # So the largest eigenvalue is 6.56, principal component is eigenvector.

# In[66]:


m = np.linalg.norm(eigenvector, ord=2)
print m


# In[67]:


NormE = eigenvector/m
print NormE


# In[68]:


PV = 6.56/8
print PV


# # So the the	 percentage of variation accounted for by the principal component is 0.82
