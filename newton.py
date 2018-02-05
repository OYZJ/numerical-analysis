
# coding: utf-8

# In[4]:


import matplotlib.pyplot as plt  
import numpy as np  


# ## Example 1: Find the root of $f(x)=x^2-1$

# In[5]:


def f(x):
    return x**2-1


# In[6]:


def df(x):
    return 2*x


# In[7]:


t = np.arange(0.0, 3.0, 0.1)  
s = f(t)  
plt.plot(t, s)  
plt.grid(True) 
plt.show()   


# In[8]:


def g(x,f,df):
    return x-f(x)/df(x)


# Initial guess is $x_0=2$

# In[9]:


x0=2
x1=g(x0,f,df)
print (x1)


# In[10]:


x2=g(x1,f,df)
print (x2)


# you can do one by one .... or even better 

# In[11]:


Nmax=6
xn=np.zeros(Nmax)
xn[0]=2
print (xn)


# In[12]:


for i in range(1,Nmax):
    xn[i]=g(xn[i-1],f,df)


# In[13]:


print (xn)


# In[14]:


t = np.arange(0.0, 3.0, 0.1)  
s = f(t)  
fn= f(xn)
plt.plot(t, s)  
plt.plot(xn, fn,"X") 
plt.grid(True) 
plt.show() 


# For plot tutorials Go to https://plot.ly/python/ipython-notebook-tutorial/ 

# To find the order of convergence define the error 

# In[15]:


er=np.zeros(Nmax)
for i in range(0,Nmax-1):
    er[i]=np.abs(xn[i+1]-xn[i])


# In[16]:


print(er)
print(xn)


# In[17]:


erp=np.zeros(Nmax)
for i in range(0,Nmax-1):
    erp[i]=er[i+1]


# In[18]:


print(erp)


# In[19]:


plt.plot(er, erp,"X") 
plt.grid(True) 
plt.show() 


# In[20]:


plt.plot(np.log(er), np.log(erp),"X") 
plt.grid(True) 
plt.show() 


# In[21]:


ratio = np.zeros(Nmax)
for i in range(0,Nmax-1):
    ratio[i]=er[i+1]/er[i]**2


# In[22]:


print (ratio)


# ## Example 2: Find the root of  $f(x)=x^2−2x+1$
#  

# In[23]:


def f2(x):
    return x**2-2*x+1


# In[24]:


def df2(x):
    return 2*x-2


# In[25]:


t = np.arange(0.0, 3.0, 0.1)  
s = f2(t)  
plt.plot(t, s)  
plt.grid(True) 
plt.show()   


#  $f(x^*=1)=0$... BUT $\frac{\partial f}{\partial x}(x^*)=0$ too!! Do you expect quadratic convergence order?
#  

# In[26]:


def g2(x,f2,df2):
    return x-f2(x)/df2(x)


# In[27]:


for i in range(1,Nmax):
    xn[i]=g2(xn[i-1],f2,df2)


# In[28]:


print(xn)


# In[29]:


er=np.zeros(Nmax)
for i in range(0,Nmax-1):
    er[i]=np.abs(xn[i+1]-xn[i])


# In[30]:


print(er)


# In[31]:


erp=np.zeros(Nmax)
for i in range(0,Nmax-1):
    erp[i]=er[i+1]


# In[32]:


plt.plot(er, erp,"X") 
plt.grid(True) 
plt.show() 


# You can see that the convergence  is slower than in Example 1. 

# In[33]:


ratio = np.zeros(Nmax)
for i in range(0,Nmax-1):
    ratio[i]=er[i+1]/er[i]


# In[34]:


print(ratio)


# the convergence order is 1 (linear)

# In[ ]:




