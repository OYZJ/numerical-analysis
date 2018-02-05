
# coding: utf-8

# # problem1 

# In[1]:


import numpy as np


# In[5]:


def trapezoidal(f,a,b):
    I = (b-a)*(f(a)+f(b))/2
    print "the evaluate solution of trapezoidal rule is {0:f}".format(I)
def simpson(f,a,b):
    I = (b-a)*(f(a)+4*f((a+b)/2)+f(b))/6
    print "the evaluate solution of simpson rule is {0:f}".format(I)
def midpoint(f,a,b):
    I = (b-a)*f((a+b)/2)
    print "the evaluate solution of midpoint rule is {0:f}".format(I)    


# In[6]:


def f(x):
    return 5*x**2+3


# In[7]:


a = 0
b = 5
trapezoidal(f,a,b)
simpson(f,a,b)
midpoint(f,a,b)


# In[ ]:


# So the analytical result is 223.33
# Only Simpson’s rule gives a exact result. Because Simpson’s rule’s degree of precision is 3, trapezoidal rule and midpoint rule’s degree of precision is 1. So in these 3 methods, only Simpson’s rule can integrate xˆ2 exactly.


# In[10]:


def f(x):
    return 2*x+1
a = 1
b = 4
trapezoidal(f,a,b)
simpson(f,a,b)
midpoint(f,a,b)


# #So the analytical result is 18.
# #All 3 rules give the exact result. Because Simpson’s rule’s degree of precision is 3, trapezoidal rule and midpoint rule’s degree of precision is 1. So all these 3 methods can integrate x exactly.

# In[12]:


def f(x):
    return np.e**x


# In[13]:


a = 1
b = 2
trapezoidal(f,a,b)
simpson(f,a,b)
midpoint(f,a,b)


# #The analytical result is 4.67077427
# #All 3 rules can not give the exact result. Because Simpson’s rule’s degree of precision is 3, trapezoidal rule and midpoint rule’s degree of precision is 1. If using Taylor Series to expand f(x)= eˆx, we can find the largest degree of polynominal goes infinite, so all 3 rule can not work here.

# In[ ]:





# # problem2

# In[22]:


def gauss_quadrature(f,a,b):
    def g(x):
        return f((b-a)*x/2+(a+b)/2)
    I = ((b-a)/2)*(g(-np.sqrt(1/3))+g(np.sqrt(1/3)))
    print "the evaluate of solution2-point Gauss quadrature rule is {0:f}".format(I)


# In[23]:


def f(x):
    return 3*x**4+5*x**2


# In[24]:


a = 1
b = 6
gauss_quadrature(f,a,b)


#  So the exact solution is not expected because the 2-point Gauss quadrature rule’s degree of precision is 3

# In[25]:


def g(x):
    return np.e**(x-5)


# In[26]:


a = 1
b = 4
gauss_quadrature(g,a,b)


#  So the exact solution is not expected because the 2-point Gauss quadrature rule’s degree of precision is 3

# In[28]:


def general_gauss_quadrature(f,a,b):
    I = (b-a)*(f((a+b)/2-np.sqrt(1/3)*(b-a)/2)+f((a+b)/2+np.sqrt(1/3)*(b-a)/2))/2
    print "the evaluate solution of 2-point Gauss quadrature rule is {0:f}".format(I)
def f(x):
    return 3*x**4+5*x**2
a = 1
b = 6
general_gauss_quadrature(f,a,b)


# In[29]:


def f(x):
    return np.e**(x-5)
a = 1
b = 4
general_gauss_quadrature(f,a,b)


# So the evaluate solution of 2-point Gauss quadrature rule is 0.344518

# In[ ]:





# # problem3

# ∫1 dt = 3w = 2
# ∫2t dt =2w ∗ (x1 + x2 + x3)=0
# ∫4t^2 − 2dt =w(4(x_1^2+x_2^2+x_3ˆ2)-6)=-4/3
# ∫8t^3 − 12tdt = w(8(x_1^3+x_2^3+x_3ˆ3)-12*(x_1+x_2+x_3))= 0
# 
# w = 2/3, x1 = 0, x2 = 1/p2, x3 = −1/p2

# In[ ]:





# # problem4

# In[33]:


def gauss_quadrature_2d(f,x11,y11,x12,y12,x21,y21,x22,y22):
    def x_uv(u,v):
        return x11*((1-u)*(1-v)/4)+x12*((1-u)*(1+v)/4)+x21*((1+u)*(1-v)/4)+x22*((1+u)*(
    def y_uv(u,v):
        return y11*((1-u)*(1-v)/4)+y12*((1-u)*(1+v)/4)+y21*((1+u)*(1-v)/4)+y22*((1+u)*(
    def dx_du(x11,x12,x21,x22,v):
        return -x11*(1-v)/4-x12*(1+v)/4+x21*(1-v)/4+x22*(1+v)/4
    def dx_dv(x11,x12,x21,x22,u):
        return -x11*(1-u)/4+x12*(1-u)/4-x21*(1+u)/4+x22*(1+u)/4
    def dy_du(y11,y12,y21,y22,v):
        return -y11*(1-v)/4-y12*(1+v)/4+y21*(1-v)/4+y22*(1+v)/4
    def dy_dv(y11,y12,y21,y22,u):
        return -y11*(1-u)/4+y12*(1-u)/4-y21*(1+u)/4+y22*(1+u)/4
    def J(x11,y11,x12,y12,x21,y21,x22,y22,u,v):
        return dx_du(x11,x12,x21,x22,v)*dy_dv(y11,y12,y21,y22,u)-dx_dv(x11,x12,x21,x22,
    n = [1/np.sqrt(3),-1/np.sqrt(3)]
    I = 0
    for i in range(0,2):
        for j in range(0,2):
            u = n[i]
            v = n[j]
            x = x_uv(u,v)
            y = y_uv(u,v)
            I = I + f(x,y)*J(x11,y11,x12,y12,x21,y21,x22,y22,u,v)
    print("the integrate result is {0:f}".format(I))


# In[35]:


def f(x,y):
    return -(x+3)**2+y**2
x11=0
y11=0
x12=0
y12=6
x21=2
y21=0
x22=2
y22=6


# gauss_quadrature_2d(f,x11,y11,x12,y12,x21,y21,x22,y22)
# the integrate result is -52.000000

# In[ ]:


def f(x,y):
    return np.e**(-(x**2+y**2))
x11=0
y11=0
x12=0
y12=6
x21=2
y21=0
x22=2
y22=6


# gauss_quadrature_2d(f,x11,y11,x12,y12,x21,y21,x22,y22)
# the integrate result is 0.552654

# In[ ]:





# # problem5

# In[ ]:


def gauss_quadrature_2d(f,x1,y1,x2,y2,x3,y3,x4,y4):
    def x_uv(u,v):
    return x1*((1-u)*(1-v)/4)+x2*((1-u)*(1+v)/4)+x3*((1+u)*(1-v)/4)+x4*((1+u)*(1+v)/
    def y_uv(u,v):
    return y1*((1-u)*(1-v)/4)+y2*((1-u)*(1+v)/4)+y3*((1+u)*(1-v)/4)+y4*((1+u)*(1+v)/
    def dxdu(x1,x2,x3,x4,v):
    return -x1*(1-v)/4-x2*(1+v)/4+x3*(1-v)/4+x4*(1+v)/4
    def dxdv(x1,x2,x3,x4,u):
    return -x1*(1-u)/4+x2*(1-u)/4-x3*(1+u)/4+x4*(1+u)/4
    def dydu(y1,y2,y3,y4,v):
    return -y1*(1-v)/4-y2*(1+v)/4+y3*(1-v)/4+y4*(1+v)/4
    def dydv(y1,y2,y3,y4,u):
    return -y1*(1-u)/4+y2*(1-u)/4-y3*(1+u)/4+y4*(1+u)/4
    def J(x1,y1,x2,y2,x3,y3,x4,y4,u,v):
    return dxdu(x1,x2,x3,x4,v)*dydv(y1,y2,y3,y4,u)-dxdv(x1,x2,x3,x4,u)*dydu(y1,y2,y3
    n = [1/np.sqrt(3),-1/np.sqrt(3)]
    Q = 0
    k=0
    for i in range(0,2):
    for j in range(0,2):
    u = n[i]
    v = n[j]
    x = x_uv(u,v)
    y = y_uv(u,v)
    Q = Q + f(x,y)*J(x1,y1,x2,y2,x3,y3,x4,y4,u,v)
    print "the integrate result of 2-point Gauss Quadrature is {0:f}".format(I)


# In[ ]:


def f(x,y):
    return x*y
x1=0
y1=0
x2=0
y2=4
1x3=2
y3=0
x4=4
y4=4
gauss_quadrature_2d(f,x1,y1,x2,y2,x3,y3,x4,y4)


# Using the 2-point Gauss Quadrature rule: 45.333333333333336

# The Jacobian: [[(v+3)/2,(u+1)/2],[0,2]

# In[ ]:





# # problem6

# In[ ]:


def gauss_quadrature_2d(f,x11,y11,x12,y12,x21,y21,x22,y22):
    def x_uv(u,v):
        return x11*((1-u)*(1-v)/4)+x12*((1-u)*(1+v)/4)+x21*((1+u)*(1-v)/4)+x22*((1+u)*(
    def y_uv(u,v):
        return y11*((1-u)*(1-v)/4)+y12*((1-u)*(1+v)/4)+y21*((1+u)*(1-v)/4)+y22*((1+u)*(
    def dx_du(x11,x12,x21,x22,v):
        return -x11*(1-v)/4-x12*(1+v)/4+x21*(1-v)/4+x22*(1+v)/4
    def dx_dv(x11,x12,x21,x22,u):
        return -x11*(1-u)/4+x12*(1-u)/4-x21*(1+u)/4+x22*(1+u)/4
    def dy_du(y11,y12,y21,y22,v):
        return -y11*(1-v)/4-y12*(1+v)/4+y21*(1-v)/4+y22*(1+v)/4
    def dy_dv(y11,y12,y21,y22,u):
        return -y11*(1-u)/4+y12*(1-u)/4-y21*(1+u)/4+y22*(1+u)/4
    def J(x11,y11,x12,y12,x21,y21,x22,y22,u,v):
        return dx_du(x11,x12,x21,x22,v)*dy_dv(y11,y12,y21,y22,u)-dx_dv(x11,x12,x21,x22,
    n = [1/np.sqrt(3),-1/np.sqrt(3)]
    I = 0
    k = 0
    for i in range(0,2):
        for j in range(0,2):
        u = n[i]
        v = n[j]
        x = x_uv(u,v)
        y = y_uv(u,v)
        print(J(x11,y11,x12,y12,x21,y21,x22,y22,u,v))
        I = I + f(x,y)*J(x11,y11,x12,y12,x21,y21,x22,y22,u,v)
    print "the integrate result of 2-point Gauss Quadrature is {0:f}".format(I)


# In[38]:


def f(x,y):
    return x*y
x11=0
y11=0
x12=0
y12=4
x21=2
y21=0
x22=4
y22=4
gauss_quadrature_2d(f,x11,y11,x12,y12,x21,y21,x22,y22)


# the integrate result of 2-point Gauss Quadrature is 125.33333334

# In[43]:


def f(x,y):
    return np.e**(x*y)
x11=-3
y11=-2
x12=-1
y12=2
x21=1
y21=-2
x22=3
y22=2
gauss_quadrature_2d(f,x11,y11,x12,y12,x21,y21,x22,y22)


# the integrate result of 2-point Gauss Quadrature is 81.3336142045678

# In[ ]:





# In[ ]:




