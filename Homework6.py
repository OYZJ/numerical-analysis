
# coding: utf-8

# # problem1

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# In[6]:


def f(y,t):
    return -5*y


# 1.

# the ODE is stable

# 2.

# Eluer's method is unstable with h = 0.5 for 1+h = -1.5<-1

# 3.

# In[7]:


def Euler(f,a,b,h,f0):
    n = (b-a)/h
    N = int(n)
    w = np.zeros(N+1)
    t = np.zeros(N+1)
    w[0] = f0
    t[0] = a
    print "if h = {0:f}".format(h)
    print "t\t\tw"
    for i in range (1, N+1):
        w[i] = w[i-1]+h*f(w[i-1],t[i-1])
        t[i] = t[i-1]+h
        print "{0:f}\t{1}".format(t[i],w[i])


# In[8]:


Euler(f, 0, 5, 0.5, 1)
Euler(f, 0, 5, 0.1, 1)


# 4.

# the backward Euler method is stable for 1-h = 3.5>1

# 5.

# In[9]:


def backward_Euler(f,a,b,h,f0):
    n = (b-a)/h
    N = int(n)
    w = np.zeros(N+1)
    t = np.zeros(N+1)
    w[0] = f0
    t[0] = a
    print "if h = {0:f}".format(h)
    print "t\t\tw"
    for i in range (1, N+1):
        w[i] = w[i-1]/(1+5*h)
        t[i] = t[i-1]+h
        print "{0:f}\t{1}".format(t[i],w[i])


# In[10]:


backward_Euler(f, 0, 5, 0.5, 1)
backward_Euler(f, 0, 5, 0.1, 1)


# In[ ]:





# # problem2

# In[11]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# 1.

# the analytical solution is about y(t) = 20.5e^(-2.8t)-19.5e^(-3.2t) 

# In[13]:


def RK(f1,f2,a,b,h,f10,f20):
    n = (b-a)/h
    N = int(n)
    global w11
    w11 = np.zeros(N+1)
    w21 = np.zeros(N+1)
    t = np.zeros(N+1)
    w11[0] = f10
    w21[0] = f20
    t[0] = a
    print "t\t\tw"
    for i in range(1,N+1):
        k11=h*f1(w11[i-1],w21[i-1],t[i-1])
        k12=h*f1(w11[i-1]+(k11/2),w21[i-1]+(k11/2),t[i-1]+(h/2))
        k13=h*f1(w11[i-1]+(k12/2),w21[i-1]+(k12/2),t[i-1]+(h/2))
        k14=h*f1(w11[i-1]+k13,w21[i-1]+k13,t[i-1]+h)
        w11[i]=w11[i-1]+(k11+2*k12+2*k13+k14)/6
        k21=h*f2(w11[i-1],w21[i-1],t[i-1])
        k22=h*f2(w11[i-1]+(k11/2),w21[i-1]+(k11/2),t[i-1]+(h/2))
        k23=h*f2(w11[i-1]+(k12/2),w21[i-1]+(k12/2),t[i-1]+(h/2))
        k24=h*f2(w11[i-1]+k23,w21[i-1]+k23,t[i-1]+h)
        w21[i]=w21[i-1]+(k21+2*k22+2*k23+k24)/6
        t[i]=t[i-1]+h
        print "{0:f}\t{1}".format(t[i],w11[i])


# In[14]:


def f1(x1,x2,t):
    return x2
def f2(x1,x2,t):
    return -6*x2-8.96*x1


# In[15]:


RK(f1, f2, 0, 8, 0.08, 1, 5)


# 3.

# In[37]:


def Euler(f1,f2,a,b,h,f10,f20):
    n = (b-a)/h
    N = int(n)
    global w12
    w12 = np.zeros(N+1)
    w22 = np.zeros(N+1)
    t = np.zeros(N+1)
    w12[0] = f10
    w22[0] = f20
    t[0] = a
    print "t\t\tw"
    for i in range(1,N+1):
        w12[i] = w12[i-1]+h*f1(w12[i-1],w22[i-1],t[i-1])
        w22[i] = w22[i-1]+h*f2(w12[i-1],w22[i-1],t[i-1])
        t[i] = t[i-1]+h
        print "{0:f}\t{1}".format(t[i],w12[i])


# In[38]:


def f1(x1,x2,t):
    return x2
def f2(x1,x2,t):
    return -6*x2-8.96*x1
Euler(f1,f2,0,8,0.08,1,5)


# 4.

# In[42]:


def Taylor(f1,f2,df1,df2,a,b,h,f10,f20):
    n = (b-a)/h
    N = int(n)
    global w13
    w13 = np.zeros(N+1)
    w23 = np.zeros(N+1)
    t = np.zeros(N+1)
    w13[0] = f10
    w23[0] = f20
    t[0] = a
    print "t\t\tw"
    for i in range(1,N+1):
        w13[i]=w13[i-1]+h*f1(w13[i-1],w23[i-1],t[i-1])+(h**2/2)*df1(w13[i-1],w23[i-1],t[i-1])
        w23[i]=w23[i-1]+h*f2(w13[i-1],w23[i-1],t[i-1])+(h**2/2)*df2(w13[i-1],w23[i-1],t[i-1])
        t[i] = t[i-1]+h
        print "{0:f}\t{1}".format(t[i],w13[i])


# In[44]:


def df1(x1,x2,t):
    return f2(x1,x2,t)
def df2(x1,x2,t):
    return -6*f2(x1,x2,t)-8.96*x2

Taylor(f1,f2,df1,df2,0,8,0.08,1,5)


# 5.

# In[45]:


t = np.arange(0,8.08,0.08)
def x(t):
    return 20.5*np.e**(-2.8*t)-19.5*np.e**(-3.2*t)

w10 = np.zeros(101)
n = np.zeros(101)
for i in range (0,100):
    w10[i] = x(n[i-1])
    n[i] = n[i-1]+0.08


# In[47]:


plt.plot(t,w10,'y')
plt.plot(t,w11,'r')
plt.plot(t,w12,'b')
plt.plot(t,w13,'g')
yellow_patch = mpatches.Patch(color='yellow',label='Analytical result')
red_patch = mpatches.Patch(color='red',label='Runge_Kutta method')
blue_patch = mpatches.Patch(color='blue',label='Euler method')
green_patch = mpatches.Patch(color='green',label='Taylor series method')
plt.legend(handles=[yellow_patch,red_patch,blue_patch,green_patch])
plt.grid(True)
plt.show()


# In[ ]:





# # problem3

# In[50]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# In[53]:


def RK(f1,f2,f3,a,b,h,f10,f20,f30):
    n = (b-a)/h
    N = int(n)
    w1 = np.zeros(N+1)
    w2 = np.zeros(N+1)
    w3 = np.zeros(N+1)
    t = np.zeros(N+1)
    w1[0] = f10
    w2[0] = f20
    w3[0] = f30
    t[0] = a
    for i in range(1,N+1):
        k11=h*f1(w1[i-1],w2[i-1],w3[i-1],t[i-1])
        k12=h*f1(w1[i-1]+(k11/2),w2[i-1]+(k11/2),w3[i-1]+(k11/2),t[i-1]+(h/2))
        k13=h*f1(w1[i-1]+(k12/2),w2[i-1]+(k12/2),w3[i-1]+(k12/2),t[i-1]+(h/2))
        k14=h*f1(w1[i-1]+k13,w2[i-1]+k13,w3[i-1]+k13,t[i-1]+h)
        w1[i]=w1[i-1]+(k11+2*k12+2*k13+k14)/6
        k21=h*f2(w1[i-1],w2[i-1],w3[i-1],t[i-1])
        k22=h*f2(w1[i-1]+(k21/2),w2[i-1]+(k21/2),w3[i-1]+(k21/2),t[i-1]+(h/2))
        k23=h*f2(w1[i-1]+(k22/2),w2[i-1]+(k22/2),w3[i-1]+(k22/2),t[i-1]+(h/2))
        k24=h*f2(w1[i-1]+k23,w2[i-1]+k23,w3[i-1]+k23,t[i-1]+h)
        w2[i]=w2[i-1]+(k21+2*k22+2*k23+k24)/6
        k31=h*f3(w1[i-1],w2[i-1],w3[i-1],t[i-1])
        k32=h*f3(w1[i-1]+(k31/2),w2[i-1]+(k31/2),w3[i-1]+(k31/2),t[i-1]+(h/2))
        k33=h*f3(w1[i-1]+(k32/2),w2[i-1]+(k32/2),w3[i-1]+(k32/2),t[i-1]+(h/2))
        k34=h*f3(w1[i-1]+k33,w2[i-1]+k33,w3[i-1]+k33,t[i-1]+h)
        w3[i]=w3[i-1]+(k31+2*k32+2*k33+k34)/6
        t[i]=t[i-1]+h
    print w1
    print w2
    print w3


# In[54]:


def f1(x1,x2,x3,t):
    return -0.013*x1-1000*x1*x3
def f2(x1,x2,x3,t):
    return -2500*x2*x3
def f3(x1,x2,x3,t):
    return -0.013*x1-1000*x1*x3-2500*x2*x3
RK(f1,f2,f3,0,5,0.0001,1,1,0)


# 2.

# In[63]:


def Euler(f1, f2, f3, a, b, h, f10, f20, f30):
    n = (b-a)/h
    N = int(n)
    w1 = np.zeros(N+1)
    w2 = np.zeros(N+1)
    w3 = np.zeros(N+1)
    t = np.zeros(N+1)
    w1[0] = f10
    w2[0] = f20
    w3[0] = f30
    t[0] = a
    for i in range(1,N+1):
        w1[i] = w1[i-1]+h*f1(w1[i-1],w2[i-1],w3[i-1],t[i-1])
        w2[i] = w2[i-1]+h*f1(w1[i-1],w2[i-1],w3[i-1],t[i-1])
        w3[i] = w3[i-1]+h*f1(w1[i-1],w2[i-1],w3[i-1],t[i-1])
        t[i] = t[i-1]+h
    print 'w1:\n',w1
    print 'w2:\n',w2
    print 'w3:\n',w3


# In[64]:


Euler(f1,f2,f3,0,5,0.0001,1,1,0)


# In[ ]:





# # problem4

# 1.

# h<0.01

# In[67]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# In[68]:


def f(x,t):
    return -200*x+200*np.sin(t)+np.cos(t)
def df(x,t):
    return -200*f(x,t)+200*np.cos(t)-np.sin(t)


# In[69]:


def Taylor(f,df,a,b,h,f10):
    n = (b-a)/h
    N = int(n)
    w = np.zeros(N+1)
    t = np.zeros(N+1)
    w[0] = f10
    t[0] = a
    for i in range(1,N+1):
        w[i]=w[i-1]+h*f(w[i-1],t[i-1])+(h**2/2)*df(w[i-1],t[i-1])
        t[i] = t[i-1]+h
    print 'w:\n',w


# In[70]:


Taylor(f,df,0,5,0.08,1)


# In[71]:


Taylor(f,df,0,5,0.12,1)


# 2.

# -2<h<0
# h<0.0513

# In[79]:


def f1(x1,x2,t):
    return 9*x1+24*x2+5*np.cos(t)-np.sin(t)/3
def f2(x1,x2,t):
    return -24*x1-51*x2-9*np.cos(t)+np.sin(t)/3
def df1(x1,x2,t):
    return 9*f1(x1,x2,t)+24*f2(x1,x2,t)-5*np.sin(t)-np.cos(t)/3
def df2(x1,x2,t):
    return -24*f1(x1,x2,t)-51*f2(x1,x2,t)+9*np.sin(t)+np.cos(t)/3


# In[80]:


def Taylor0(f1,f2,df1,df2,a,b,h,f10,f20):
    n = (b-a)/h
    N = int(n)
    w1 = np.zeros(N+1)
    w2 = np.zeros(N+1)
    t = np.zeros(N+1)
    w1[0] = f10
    w2[0] = f20
    t[0] = a
    for i in range(1,N+1):
        w1[i]=w1[i-1]+h*f1(w1[i-1],w2[i-1],t[i-1])+(h**2/2)*df1(w1[i-1],w2[i-1],t[i-1])
        w2[i]=w2[i-1]+h*f2(w1[i-1],w2[i-1],t[i-1])+(h**2/2)*df2(w1[i-1],w2[i-1],t[i-1])
        t[i] = t[i-1]+h
    print 'w1:\n',w1
    print 'w2:\n',w2


# In[81]:


Taylor0(f1,f2,df1,df2,0,5,0.04,4/3,2/3)


# In[82]:


Taylor0(f1,f2,df1,df2,0,5,0.06,4/3,2/3)


# 3.

# -2<h<0  
# h<0.0513

# In[84]:


def f1(x1,x2,t):
    return -20*x1-19*x2
def f2(x1,x2,t):
    return -19*x1-20*x2
def df1(x1,x2,t):
    return -20*f1(x1,x2,t)-19*f2(x1,x2,t)
def df2(x1,x2,t):
    return -19*f1(x1,x2,t)-20*f1(x1,x2,t)


# In[86]:


Taylor0(f1,f2,df1,df2,0,5,0.04,2,0)


# In[87]:


Taylor0(f1,f2,df1,df2,0,5,0.06,2,0)


# In[ ]:




