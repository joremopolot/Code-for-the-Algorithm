#!/usr/bin/env python
# coding: utf-8

# In[45]:


from sympy import *
import numpy as np
from sympy import Matrix,ones
from sympy.ntheory import discrete_log
import sympy as sp

#1:maximum upper bound for the number integers encrypted and decrypted
m=6
n=sp.binomial(2*m-1,m-1)
print("n=",n)

#2:Preliminarymatrix PM
from sympy.matrices import Matrix, ones
m=6
PM=ones(m)
for r in range(1,m):
    for c in range(m-2,-1,-1):
        PM[r,c]=PM[r-1,c]+PM[r,c+1]
pprint(PM)

#3:Secret matrix M
from sympy.matrices import Matrix
from sympy.ntheory import randprime
m=6
primes=[]
while len(primes) < m**2:
    P=randprime(100,1000)
    if P not in primes:
        primes.append(P)
M=Matrix(m,m,primes)
pprint(M)

#4: Find the largest prime in the matrix and determine its number of digits
largest_prime = np.max(M)
d = len(str(largest_prime))
print("largest_prime=",largest_prime)
print("d=",d)

#5: Generate a prime with a specific number of digits

p = sp.randprime(10**((m-1)*d), 10**(((m-1)*d)+1))
d_prime =len(str(p))
print("p=",p)
print("d_prime=",d_prime)

#6:primitive root of p
import sympy as sp

g=sp.ntheory.residue_ntheory.primitive_root(p)
print("g=",g)


# In[46]:


#7:public key matrix PK.
from sympy.ntheory import discrete_log
PK=M.copy()
for i in range(m):
    for j in range(m):
        PK[i,j] = discrete_log(p, M[i,j] ,g)
PK


# In[47]:


# 8. The method of encryption
import numpy as np
def Ciphertext(N, PM, PK):
    N=180
    S = []
    Total = 0
    PK_value = []
    Loc_Addends = []
    if N == 0:
        Loc_Addends.append(0)
        S.append(0)
        return Loc_Addends, S
    r = m-1
    c = 0
    while N > 0:
        e = PM[r,c]
        if e <= N:
            S.append(e)
            Loc_Addends.append([r,c])
            PK_value.append(PK[r,c])
            Total += PK[r,c]

            N = N-e
            if N != 0:
                c = c + 1
                if c < m:
                    e = PM[r,c]
                else:
                    return 0
        else:
            r = r - 1
            if r >= 0:
                e = PM[r,c] 
            else:
                return 0
#         Total = sum(list(PK_value))#.sum()
    print("C=",Total)
    return Total


# In[48]:


C = Ciphertext(180, PM, PK)


# In[65]:



def decryption(T, M, PM):
    T = pow(g,int(C),p)
    Plaintext =0
    r=m-1
    c=0
    while T > 1:
        a=M[r,c]
        if T%M[r,c]==0:
            Plaintext +=PM[r,c]
            T=T//M[r,c]
            if T !=1:
                c+=1
                if c < m:
                    a=M[r,c]
                else:
                    return 0
        else:
            r-=1
            if r>=0:
                a=M[r,c]
            else:
                return 0
    print("T=",T)       
    print("N=",Plaintext) 
                    


# In[66]:


decryption(Exp, M, PM)


# In[ ]:




