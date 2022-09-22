import numpy as np
from scipy.integrate import *
import matplotlib.pyplot as plt
import math
from math import exp
import random
import copy

KM = 100 #Carrying capacity
d = 0.05 #Background death rate
r = 0.25 #Intrinsic growth rate
k = 1 #Evolvability
B = 0 #Aggressiveness coefficient
sk = 12.5 #Range of resources (affects functional carrying capacity, stabilizing selection)
sa=100 #Species niche width (affects competition coefficient, divergent/convergent evolution)
gamma = 0 #Strategy that maximizes carrying capacity
num = 5 #Number of morphs

time = 5000 #Time for simulation
start = 1500 #Start of environmental perturbation
end = 2500 #End of environmental perturbation

IC = [] #Initial Conditions

def gen_morphs(n):
    for i in range(n):
        var1 = random.uniform(1,10)
        var2 = random.uniform(0,1)
        IC.append(var1)
        IC.append(var2)
        
gen_morphs(num)

def carry(v):
    return KM * math.exp(-((v - gamma) ** 2) / (2 * sk))

def evoLV(X, t):       

    strats = X[1::2]
    pops = X[::2]
    
    temp = []
    
    if t>start and t<end:
        sa=100
    else:
        sa=2

    def comp(n,strats):
        a=0
        i=0
        while i < num:
            a += (1 + math.exp(-(n - strats[i] + B) ** 2 / (2 * sa)) - math.exp(-(B ** 2 / (2 * sa))))*pops[i]
            i+=1
        return a
    
    def pop_update():
        i=0
        while i <num:
            G = r/carry(strats[i]) * (carry(strats[i]) - comp(strats[i],strats)) - d*k
            a = pops[i]*G
            pert = copy.deepcopy(strats)
            pert[i] = strats[i]+.0001 
            Gmod = r/carry(pert[i]) * (carry(pert[i]) - comp(pert[i],pert)) - d*k
            b = k*(Gmod-G)/.0001 
            temp.append(a)
            temp.append(b)
            i+=1
            
    pop_update()
    
    final=np.array(temp)

    return final

intxv = np.array(IC)
time_sp = np.linspace(0,time,time*10-1)
pop = odeint(evoLV, intxv, time_sp,hmax=1)
plt.figure()
plt.subplot(211)
plt.title('Expansion and Niche Collapse')
i=0
while i<num*2:
    if i%2==0:
        if pop[:,i][-1]<1:
            plt.plot(time_sp,pop[:,i],lw=3,color='k')
        else:
            plt.plot(time_sp,pop[:,i],lw=3)
    i+=1
plt.grid(True)
plt.ylabel('Pop Size, x')
plt.subplot(212)
j=0
while j<num*2:
    if j%2==1:
        if pop[:,j-1][-1]<1:
            plt.plot(time_sp,pop[:,j],lw=3,color='k')
        else:
            plt.plot(time_sp,pop[:,j],lw=3)
    j+=1
plt.grid(True)
plt.ylabel('Indv Strategy, v')
plt.xlabel('Time')
plt.show()


def AL(time_G):
    
    popn = pop[time_G-1][::2]
    strats = pop[time_G-1][1::2]
    ranup = max(max(strats)+.1,.5)
    ranlo = min(min(strats)-.1,-.5)

    if time_G>start*10 and time_G<end*10:
        sa=100
    else:
        sa=2

    def comp(n):
        a=0
        i=0
        while i <num:
            a += (1 + math.exp(-(n - strats[i] + B) ** 2 / (2 * sa)) - math.exp(-(B ** 2 / (2 * sa))))*popn[i]
            i+=1
        return a


    def G(n):
        return r/carry(n) * (carry(n) - comp(n)) - d*k

    tem = []

    for n in np.arange(ranlo,ranup,.1):
        tem.append(G(n))


    plt.plot(np.arange(ranlo,ranup,.1),tem,lw=3,color='darkgrey')

    i=0
    while i<num*2:
        if i%2!=0:
            plt.plot(pop[time_G][i],G(pop[time_G][i]),marker='o',lw=20,markersize=10)
        i+=1

    plt.title('Expansion and Niche Collapse: Time '+str(int(time_G/10)))
    plt.xlabel('Evolutionary Strategy: v')
    plt.ylabel('Fitness: G')
    plt.show()
    
for i in np.arange(0,time*10,500):
    AL(i)