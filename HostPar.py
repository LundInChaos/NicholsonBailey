# A script to iterate populations of hosts and parasites using the Nicholson-
# -Bailey model. It Includes three populations (genotypes) of hosts (H[1,*], H[2,*], H[3,*])
# and three populations of parasites (P[1,*], P[2,*], P[3,*]) where 
# P[1,*] only infects H[1,*] etc. The mixing facors then determine 
# how much genetic mixing there is in host and parasite populations, normally
# 1 (complete mixing) for hosts and somewhere in between 0 and 1 for parasites.
# The functions used are in the file "HostPar_funcs.py" and come from an artice
# by Flatt, J. theor. Biol. (2001) 212, 345}354 (doi:10.1006/jtbi.2001.2380).

import numpy as np
import matplotlib.pyplot as plt
from HostPar_funcs import NB_Mix, NB_H, NB_P      # import user-made functions

K = 10             # Host carrying capacity
a = 0.45           # Parasite search efficiency
c = 1              # Parasite fecundity

# Mixing factors:
SFH = 1               # assume a sexual host, like most mammals etc.
SFP = 0.5             # parasites can have variable levels of mixing


# Lambda is the base host growth factor. Here we let it vary between 1 and 50
# like in the artice by Flatt et al.

# Lambdas:
L1 = 1
L2 = 50
Ls = 1/2
L = np.arange(L1,L2,Ls)

N = 1000                   # Total number of steps
H = np.zeros( (3, N) )     # Hosts
P = np.zeros( (3, N) )     # Parasites
Htot = np.zeros(N)         # Total hosts
Ptot = np.zeros(N)         # Total parasites

# Initial values:
H[0,0] = 4
H[1,0] = 3
H[2,0] = 3
P[0,0] = 0.3
P[1,0] = 0.36
P[2,0] = 0.34
Htot[0] = H[0,0] + H[1,0] + H[2,0]
Ptot[0] = P[0,0] + P[1,0] + P[2,0]

# The following vectors are used to store the last population densities of each
# population. They are used to plot the bifurcation diagrams:
# End values:
EN = 10
E = np.zeros( (EN,len(L)) )
F = np.zeros( (EN,len(L)) )

# The mixing variable is calculated with the population densities:
# Initial mixing variable:
Hp = NB_Mix(H[0,0], H[1,0], Htot[0])
Pq = NB_Mix(P[0,0], P[1,0], Ptot[0])


# Print conditions:
print("K = ",K)
print("a = ",a)
print("c = ",c)
print("H0 = ",Htot[0])
print("P0 = ",Ptot[0])
print("p = ",Hp)
print("q = ",Pq)
print("Host mixing factor = ",SFH)
print("Para mixing factor = ",SFP)

# Calculate for each lambda:
for i in range(0,len(L)):
    Lam = L[i]
    # For each time step:
    for ii in range(1,N-1):
        # For all genotypes:
        for h in range(0,2):
            # Base population densities:
            H[h,ii] = NB_H(H[h,ii-1],Htot[ii-1],P[h,ii-1],Lam,K,a)
            P[h,ii] = NB_P(H[h,ii-1],P[h,ii-1],a,c)
        # Total pops:
        Htot[ii] = H[0,ii] + H[1,ii] + H[2,ii]
        Ptot[ii] = P[0,ii] + P[1,ii] + P[2,ii]
        # Mixing factors:
        Hp = NB_Mix(H[0,ii],H[1,ii],Htot[ii])
        Pq = NB_Mix(P[0,ii],P[1,ii],Ptot[ii])
        # New, mixed pops:
        H[0,ii] = (1-SFH)*H[0,ii] + SFH* Hp**2       *Htot[ii]
        H[1,ii] = (1-SFH)*H[1,ii] + SFH* 2*Hp*(1-Hp) *Htot[ii]
        H[2,ii] = (1-SFH)*H[2,ii] + SFH* (1-Hp)**2   *Htot[ii]
        P[0,ii] = (1-SFP)*P[0,ii] + SFP* Pq**2       *Ptot[ii]
        P[1,ii] = (1-SFP)*P[1,ii] + SFP* 2*Pq*(1-Pq) *Ptot[ii]
        P[2,ii] = (1-SFP)*P[2,ii] + SFP* (1-Pq)**2   *Ptot[ii]
        # New total pops:
        Htot[ii] = H[0,ii] + H[1,ii] + H[2,ii]
        Ptot[ii] = P[0,ii] + P[1,ii] + P[2,ii]
    # Store last values:
    for e in range(0,EN):
        E[e,i] = Htot[N-e-1]
        F[e,i] = Ptot[N-e-1]

# time for plotting (if you want to se the evolution of a certain pop):
plott = np.arange(N)

# Plotts:
plt.figure(1)
for e in range(1,EN):
    plt.scatter(L, E[e,:],c='k')
    plt.title("Hosts")
    plt.xlabel("lambda")
    plt.ylabel("Density")
plt.figure(2)
for e in range(1,EN):
    plt.scatter(L, F[e,:],c='r')
    plt.title("Parasites")
    plt.xlabel("lambda")
    plt.ylabel("Density")
