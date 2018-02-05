
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt  
import numpy as np  



# ## Example 1: Find the root of $f(x)=x^2-4sin(x)$

# In[3]:


def f(x):
    return x**2-4*np.sin(x)


# In[5]:


t = np.arange(-1.0, 4.0, 0.2)  
s = f(t)  
plt.plot(t, s)  
plt.grid(True) 
plt.show()   


# For secant method we need two initial guesses is $x_0=1$ and $x_1=3$. We will perform 9 iterations, $N_{max}=9$.

# In[33]:


Nmax=12
xn=np.zeros(Nmax)
xn[0]=1
xn[1]=3
print (xn)


# In[34]:


for i in range(1,Nmax-1):
    xn[i+1]=xn[i]-f(xn[i])*(xn[i]-xn[i-1])/(f(xn[i])-f(xn[i-1]))


# In[35]:


print (xn)


# In[36]:


print(f(xn))


# In[12]:


t = np.arange(-1.0, 4.0, 0.2)  
s = f(t)  
fn= f(xn)
plt.plot(t, s)  
plt.plot(xn, fn,"X") 
plt.grid(True) 
plt.show() 


# To find the order of convergence define the error 

# In[37]:


er=np.zeros(Nmax)
for i in range(0,Nmax-1):
    er[i]=np.abs(xn[i+1]-xn[i])


# In[38]:


print(er)


# In[39]:


erp=np.zeros(Nmax)
for i in range(0,Nmax-1):
    erp[i]=er[i+1]


# In[18]:


plt.plot(er, erp,"X") 
plt.grid(True) 
plt.show() 


# In[41]:


plt.plot(np.log(er), np.log(erp),"X") 
plt.plot(np.log(er), np.log(er))
plt.plot(np.log(er), np.log(er*er))
plt.grid(True) 
plt.show() 

