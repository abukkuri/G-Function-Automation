import numpy as np
from scipy.integrate import *
import matplotlib.pyplot as plt
import math
from math import exp
import random
import copy

KM = 100 #Carrying capacity
d = 0 #Background death rate
r = 0.25 #Intrinsic growth rate: prey
r2 = 0.25 #Intrinsic growth rate: predator
k = .4 #Evolvability
B = 0 #Aggressiveness coefficient
bm = 0.2 #Maximal probability of prey capture, important
c = 0.25 #Conversion efficiency
sk = 12.5 #Range of resources
sb = 10 #10 #Range of predation, important
sa = 100 #Species niche width (affects competition coefficient, divergent/convergent evolution)
gamma = 0 #Strategy that maximizes carrying capacity


IC_prey = [] #Initial conditions: prey
IC_pred = [] #Initial conditions: predator
num_prey = 3 #Number of prey morphs
num_pred = 2 #Number of predator morphs
total_num = num_pred+num_prey #Total number of morphs

time = 5000 #Time for simulation
start = 1500 #Start of environmental perturbation
end = 2500 #End of environmental perturbation

def gen_prey(n):
    for i in range(n):
        var1 = random.uniform(2,10)
        var2 = random.uniform(-.1,.1)
        IC_prey.append(var1)
        IC_prey.append(var2)
        
def gen_pred(n):
    for i in range(n):
        var1 = random.uniform(2,5)
        var2 = random.uniform(-.1,.1)
        IC_pred.append(var1)
        IC_pred.append(var2)
        
gen_prey(num_prey)
gen_pred(num_pred)

IC = IC_prey+IC_pred

def carry(v):
    return KM * math.exp(-((v - gamma) ** 2) / (2 * sk))

def evoLV(X, t):       

    preys = X[:num_prey*2]
    preds = X[num_prey*2:]
    
    prey_strats = preys[1::2]
    prey_pops = preys[::2]
    
    pred_strats = preds[1::2]
    pred_pops = preds[::2]
    
    temp_prey = []
    temp_pred = []
    
    if t>start and t<end:
        sa=100
    else:
        sa=2

    def comp(n,strats):
        a=0
        i=0
        while i <num_prey:
            a += (1 + math.exp(-(n - strats[i] + B) ** 2 / (2 * sa)) - math.exp(-(B ** 2 / (2 * sa))))*prey_pops[i]
            i+=1
        return a
    
    def die(n,pred_strats):
        b=0
        i=0
        while i < num_pred:
            b += bm*(math.exp(-(n - pred_strats[i]) ** 2 / sb))*pred_pops[i]
            i+=1
        return b
    
    def nom(n,prey_strats):
        b=0
        i=0
        while i < num_prey:
            b += c*bm*(math.exp(-(n - prey_strats[i]) ** 2 / sb))*prey_pops[i]
            i+=1
        return b
    
    def prey_update():
        i=0
        while i <num_prey:
            G = r/carry(prey_strats[i]) * (carry(prey_strats[i]) - comp(prey_strats[i],prey_strats)) - d*k - die(prey_strats[i],pred_strats)
            a = prey_pops[i]*G
            pert = copy.deepcopy(prey_strats)
            pert[i] = prey_strats[i]+.0001 
            Gmod = r/carry(pert[i]) * (carry(pert[i]) - comp(pert[i],pert)) - d*k - die(pert[i],pred_strats)
            b = k*(Gmod-G)/.0001 
            temp_prey.append(a)
            temp_prey.append(b)
            i+=1
            
    def pred_update():
        i=0
        while i <num_pred:
            G = r2 * (1-sum(pred_pops)/(nom(pred_strats[i],prey_strats))) - d*k
            a = pred_pops[i]*G
            pert = copy.deepcopy(pred_strats)
            pert[i] = pred_strats[i]+0.0001 
            Gmod = r2 * (1-sum(pred_pops)/(nom(pert[i],prey_strats))) - d*k
            b = k*(Gmod-G)/.0001 
            temp_pred.append(a)
            temp_pred.append(b)
            i+=1
                
    prey_update()
    pred_update()
    
    temp = temp_prey+temp_pred
    
    final=np.array(temp)

    return final

intxv = np.array(IC)
time_sp = np.linspace(0,time,time*10-1)
pop = odeint(evoLV, intxv, time_sp,hmax=1)
plt.figure()
plt.subplot(211)
plt.title('Adaptive Radiation')
i=0
while i<total_num*2:
    if i%2==0:
        if i<num_prey*2:
            plt.plot(time_sp,pop[:,i],lw=3)
        else:
            plt.plot(time_sp,pop[:,i],lw=3,color='k')
        '''if pop[:,i][-1]<1:
            plt.plot(time_sp,pop[:,i],lw=3,color='k')
        else:
            plt.plot(time_sp,pop[:,i],lw=3)'''
    i+=1
plt.grid(True)
plt.ylabel('Pop Size, x')
plt.subplot(212)
j=0
while j<total_num*2:
    if j%2==1:
        if j<num_prey*2:
            plt.plot(time_sp,pop[:,j],lw=3)
        else:
            plt.plot(time_sp,pop[:,j],lw=3,color='k')
        '''if pop[:,j-1][-1]<1:
            plt.plot(time_sp,pop[:,j],lw=3,color='k')
        else:
            plt.plot(time_sp,pop[:,j],lw=3)'''
    j+=1
plt.grid(True)
plt.ylabel('Indv Strategy, v')
plt.xlabel('Time')
plt.show()


def AL(time_G):
    preys = pop[time_G-1][:num_prey*2]
    preds = pop[time_G-1][num_prey*2:]
    
    print(preys)
    print(preds)
    
    prey_strats = preys[1::2]
    prey_pops = preys[::2]
    
    pred_strats = preds[1::2]
    pred_pops = preds[::2]
        
    ranup = max(max(max(prey_strats),max(pred_strats))+.1,.5)
    ranlo = min(min(min(prey_strats),min(pred_strats))-.1,-.5)

    if time_G>start*10 and time_G<end*10:
        sa=100
    else:
        sa=2
    
    def comp(n):
        a=0
        i=0
        while i < num_prey:
            a += (1 + math.exp(-(n - prey_strats[i] + B) ** 2 / (2 * sa)) - math.exp(-(B ** 2 / (2 * sa))))*prey_pops[i]
            i+=1
        return a
    
    def die(n):
        b=0
        i=0
        while i < num_pred:
            b += bm*(math.exp(-(n - pred_strats[i]) ** 2 / sb))*pred_pops[i]
            i+=1
        return b
    
    def nom(n):
        b=0
        i=0
        while i < num_prey:
            b += c*bm*(math.exp(-(n - prey_strats[i]) ** 2 / sb))*prey_pops[i]
            i+=1
        return b


    def Gprey(n):
        return r/carry(n) * (carry(n) - comp(n)) - d*k - die(n)
    
    def Gpred(n):
        return r2 * (1-sum(pred_pops)/(nom(n))) - d*k
    
    tem1 = []
    tem2 = []

    for n in np.arange(ranlo,ranup,.1):
        tem1.append(Gprey(n))
        tem2.append(Gpred(n))


    plt.plot(np.arange(ranlo,ranup,.1),tem1,lw=3,color='darkgrey',label='Prey')
    plt.plot(np.arange(ranlo,ranup,.1),tem2,lw=3,color='c',label='Predator')

    for i in range(total_num*2):
        if i%2!=0:
            ind = int((i-1)/2)
            if i<num_prey*2:
                if pop[:,ind][time_G]<.1:
                    #plt.plot(pop[time_G][i],Gprey(pop[time_G][i]),marker='o',lw=20,markersize=10)
                    plt.plot(pop[time_G][i],Gprey(pop[time_G][i]),marker='o',color='r',lw=20,markersize=10)
                else:
                    plt.plot(pop[time_G][i],Gprey(pop[time_G][i]),marker='o',color='k',lw=20,markersize=10)
            else:
                if pop[:,ind][time_G]<.1:
                    #plt.plot(pop[time_G][i],Gpred(pop[time_G][i]),marker='o',lw=20,markersize=10)
                    plt.plot(pop[time_G][i],Gpred(pop[time_G][i]),marker='*',color='r',lw=20,markersize=10)
                else:
                    plt.plot(pop[time_G][i],Gpred(pop[time_G][i]),marker='*',color='k',lw=20,markersize=10)
                
        
    tmin = min(min(tem1),min(tem2))
    tmax = max(max(tem1),max(tem2))
        
    plt.title('Expansion and Niche Collapse: Time '+str(int(time_G/10)))
    plt.xlabel('Evolutionary Strategy: v')
    plt.legend()
    plt.ylim(bottom=tmin)#max(tmin,-.05))
    plt.ylim(top=tmax)
    #plt.ylim(top=min(tmax,.2))
    #plt.ylim(-.05,.125)
    plt.ylabel('Fitness: G')
    plt.show()
    #plt.savefig(str(int(time_G/10)))
    
for i in np.arange(0,time*10,500):
    AL(i)