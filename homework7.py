
# coding: utf-8

# # problem1

# In[33]:


import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


# In[22]:


def FDM(p,q,r,a,b,a0,b0,h):
    n = int((b-a)/h)
    A = np.zeros((n+1,n+1))
    w = np.zeros(n+1)
    B = np.zeros(n+1)
    x = np.zeros(n+1)
    A[0,0] = A[n,n] = 1
    A[0,1] = A[n,n-1] = 0
    B[0] = a0
    B[n] = b0
    x[0] = a
    for i in range(1,n):
        x[i] = x[i-1]+h
        B[i] = -h**2*r(x[i])
        A[i,i-1] = -1-h*p(x[i])/2
        A[i,i] = 2+h**2*q(x[i])
        A[i,i+1] = -1+h*p(x[i])/2
    print A,B
    gaussElim_pp(A,B,x)


# In[23]:


def gaussElim_pp(h,k,x):
    n = len(k)
    r = np.arange(n)
    p = np.zeros(n)
    for j in range(0, n - 1):
        max=0
        for c in range(j,n):
            if np.abs(h[r[c],j])>max:
                max=np.abs(h[r[c],j])
        for m in range(j, n):
            if np.abs(h[r[m], j]) == max:
                u = r[j]
                r[j] = r[m]
                r[m] = u
                break
        for i in range(j + 1, n):
            lam = h[r[i], j] / h[r[j], j]
            h[r[i], j:n] = h[r[i], j:n] - lam * h[r[j],j:n]
            k[r[i]] = k[r[i]] - lam * k[r[j]]
    for j in range(n-1, -1, -1):
        k[r[j]] = (k[r[j]] - h[r[j], j+1:n].dot(k[r[j+1:n]])) / h[r[j], j]
        p[j] = k[r[j]]
    plt.plot(x[0:n-1],p[0:n-1])
    plt.grid(True)
    plt.show()


# In[24]:


def p(x):
    return -(x+1)
def q(x):
    return 2
def r(x):
    return (1-x**2)*np.e**(-x)
FDM(p,q,r,0,1,-1,0,0.1)


# In[ ]:





# # problem2

# In[26]:


def p(x):
    return -3
def q(x):
    return 0
def r(x):
    return x**2+np.sin(x)
FDM(p,q,r,-5,13.2,10,23,0.1)


# In[ ]:





# # problem3

# In[28]:


def p(x):
    return -2/x
def q(x):
    return 0
def r(x):
    return -1
FDM(p,q,r,1,2,0,-0.5,0.1)


# In[ ]:





# # problem4

# In[36]:


def possion(f,g,a,b,c,d,h):
    N=int((b-a)/h)
    M=int((d-c)/h)
    A=np.zeros(((N-1)*(M-1),(N-1)*(M-1)))
    A[0,0]=4
    A[N-2,N-2]=4
    A[0,1]=-1
    A[N-2,N-3]=-1
    A[0,N-1]=-1
    A[N-2,2*N-3]=-1
    x=np.zeros(N+1)
    y=np.zeros(M+1)
    for i in range(0,N+1):
        x[i]=a+i*h
    for i in range(0,M+1):
        y[i]=c+i*h
    for i in range(1,N-2):
        A[i,i-1]=-1
        A[i,i+1]=-1
        A[i,i]=4
        A[i,i+N-1]=-1
    A[(N-1)*(M-2),(N-1)*(M-2)]=4
    A[(N-1)*(M-1)-1,(N-1)*(M-1)-1]=4
    A[(N-1)*(M-2),(N-1)*(M-2)+1]=-1
    A[(N-1)*(M-1)-1,(N-1)*(M-1)-2]=-1
    A[(N-1)*(M-2),(N-1)*(M-3)]=-1
    A[(N-1)*(M-1)-1,(N-1)*(M-2)-1]=-1
    for i in range((N-1)*(M-2)+1,(N-1)*(M-1)-1):
        A[i,i-1]=-1
        A[i,i+1]=-1
        A[i,i]=4
        A[i,i-N+1]=-1
    for i in range(1,M-2):
        A[(N-1)*i,(N-1)*i]=4
        A[(N-1)*(i+1)-1,(N-1)*(i+1)-1]=4
        A[(N-1)*i,(N-1)*i+1]=-1
        A[(N-1)*(i+1)-1,(N-1)*(i+1)-2]=-1
        for j in range(0,N-1):
            A[(N-1)*i+j,(N-1)*(i-1)+j]=-1
            A[(N-1)*i+j,(N-1)*(i+1)+j]=-1
        for j in range(0,N-3):
            A[(N-1)*i+1+j,(N-1)*i+1+j]=4
            A[(N-1)*i+1+j,(N-1)*i+j]=-1
            A[(N-1)*i+1+j,(N-1)*i+2+j]=-1
    B=np.zeros((N-1)*(M-1))
    for j in range(1,M):
        for i in range(1,N):
            B[(j-1)*(N-1)+i-1]=-h**2*f(x[i],y[j])
            if j==1:
                B[(j-1)*(N-1)+i-1]=B[(j-1)*(N-1)+i-1]+g(x[i],y[0])
            if j==M-1:
                B[(j-1)*(N-1)+i-1]=B[(j-1)*(N-1)+i-1]+g(x[i],y[M])
            if i==1:
                B[(j-1)*(N-1)+i-1]=B[(j-1)*(N-1)+i-1]+g(x[0],y[j])
            if i==N-1:
                B[(j-1)*(N-1)+i-1]=B[(j-1)*(N-1)+i-1]+g(x[N],y[j])
        w=np.linalg.solve(A,B)
        return w,x,y


# In[37]:


def f(x,y):
    return x+y
def g(x,y):
    if x == 0:
        return 0
    if x == 1:
        return 1
    if y == 0:
        return 0
    if y == 1:
        return 1
w,x,y=possion(f,g,0.,1.,0.,1.,0.1)
print(w)


# In[38]:


W=np.zeros((11,11))
for i in range(0,9):
    W[i+1,1:9]=w[9*i:9*i+8]
X, Y = np.meshgrid(x, y)
plt.figure()
CS = plt.contour(X, Y, W)
plt.clabel(CS, inline=1, fontsize=10)
plt.show()


# In[ ]:





# # problem5

# In[39]:


def h(f,ua,ub,a,b,t,hx,ht,D):
    N=int((b-a)/hx)
    n=int(t/ht)
    I=np.eye(N-1)
    B=np.zeros(N-1)
    B[0]=-ua(t)
    B[N-2]=-ub(t)
    A=np.zeros((N-1,N-1))
    A[0,0]=2
    A[N-2,N-2]=2
    A[0,1]=-1
    A[N-2,N-3]=-1
    for i in range(1,N-2):
        A[i,i]=2
        A[i,i-1]=-1
        A[i,i+1]=-1
    l=D*ht/hx**2
    w=np.zeros((N-1,n+1))
    x=np.zeros(N+1)
    x[0]=a
    for i in range(1,N+1):
        x[i]=x[i-1]+hx
    for i in range(0,N-1):
        w[i,0]=f(x[i+1])
    for i in range(1,n+1):
        tt=(i-1)*ht
        B[0]=-ua(tt)
        B[N-2]=-ub(tt)
        w[:,i]=(I-l*A).dot(w[:,i-1])-l*B
    W=np.zeros(N+1)
    W[0]=ua(t)
    W[N]=ub(t)
    W[1:N-1]=w[0:N-2,n]
    return W,x


# In[40]:


def f(x):
    return 0
def ua(t):
    if t<=3:
        return 5*t/3
    if t>3:
        return 5*np.e**(-(t-3)/5)
def ub(t):
    if t<=3:
        return 5*t/3
    if t>3:
        return 5*np.e**(-(t-3)/5)
W,x=h(f,ua,ub,0.,1100.,3.,5.,0.02,509.76)
print(W)
plt.plot(x,W,'*')
plt.grid(True)
plt.show()


# In[41]:


W1,x=h(f,ua,ub,0.,1100.,10.,5.,0.02,509.76)
W2,x=h(f,ua,ub,0.,1100.,15.,5.,0.02,509.76)
W3,x=h(f,ua,ub,0.,1100.,20.,5.,0.02,509.76)
plt.plot(x,W1,'r')
plt.plot(x,W2,'y')
plt.plot(x,W3,'b')
plt.grid(True)
plt.show()


# In[ ]:




