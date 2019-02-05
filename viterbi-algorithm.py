# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 12:58:27 2018

@author: Khali Mohsin 260631318
"""

#Fasta import 
from Bio import SeqIO

fasta_sequences = list(SeqIO.parse(open("hw3_proteins.fa.txt"),'fasta'))

#=======================================[ SET Up]=========================================#
#Transmission_Prob
T = {
     'Hydrophobic':{'Hydrophobic': math.log(0.8), 'Hydrophilic': math.log(0.04), 'Mixed': math.log(0.16)},
     'Hydrophilic':{'Hydrophobic': math.log(0.0375), 'Hydrophilic': math.log(0.875), 'Mixed': math.log(0.0875)},
     'Mixed':{'Hydrophobic': math.log((1/7) * 0.50), 'Hydrophilic': math.log((1/7) * 0.50), 'Mixed': math.log(6/7)}
     }


H_AA = ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
P_AA = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P']

#States
S = ('Hydrophobic','Hydrophilic','Mixed')

#Start_Prob
start = {'Hydrophobic': math.log(1.0 / 3), 'Hydrophilic': math.log(1.0 / 3), 'Mixed': math.log(1.0 / 3)}


#Emission_Prob
E = {'Hydrophobic': {}, 'Hydrophilic': {}, 'Mixed': {}}
E['Hydrophobic'] = E['Hydrophobic'].fromkeys(H_AA, math.log(0.9 * 1.0 / 8))
E['Hydrophobic'].update(E['Hydrophobic'].fromkeys(P_AA, math.log(0.1 * 1.0 / 12)))
E['Hydrophilic'] = E['Hydrophilic'].fromkeys(H_AA, math.log(0.2 * 1.0 / 8))
E['Hydrophilic'].update(E['Hydrophilic'].fromkeys(P_AA, math.log(0.8 * 1.0 / 12)))
E['Mixed'] = E['Mixed'].fromkeys(H_AA, math.log(1.0 / 20))
E['Mixed'].update(E['Mixed'].fromkeys(P_AA, math.log(1.0 / 20)))


#=======================================[ Globals]=========================================#
#Max Value of Segment size
maxva = 0;

L_mixed_frac=0;
L_mixed_frac_seq="";

#sequence with max segment
maxva_seq="";

# Total number of H P and Ms
total_H=0;
total_P=0;
total_M=0;


# The distribution of H P M 
h_dist=[]
p_dist=[]
m_dist=[]

for i in range (0,5000):
    h_dist.append(0)
    p_dist.append(0)
    m_dist.append(0)


#Frequency dictionary 
h_freq = {'A': 0, 'V': 0, 'I': 0, 'L': 0, 'M': 0, 'F': 0, 'Y': 0, 'W': 0, 'R': 0, 'H': 0, 'K': 0, 'D': 0,
                    'E': 0, 'S': 0, 'T': 0, 'N': 0, 'Q': 0, 'S': 0, 'C': 0, 'G': 0, 'P': 0}
p_freq = {'A': 0, 'V': 0, 'I': 0, 'L': 0, 'M': 0, 'F': 0, 'Y': 0, 'W': 0, 'R': 0, 'H': 0, 'K': 0, 'D': 0,
                    'E': 0, 'S': 0, 'T': 0, 'N': 0, 'Q': 0, 'S': 0, 'C': 0, 'G': 0, 'P': 0}
m_freq = {'A': 0, 'V': 0, 'I': 0, 'L': 0, 'M': 0, 'F': 0, 'Y': 0, 'W': 0, 'R': 0, 'H': 0, 'K': 0, 'D': 0, 'E': 0,
              'S': 0, 'T': 0, 'N': 0, 'Q': 0, 'S': 0, 'C': 0, 'G': 0, 'P': 0}

# VITERBI (O,S,\Pi ,Y,A,B):X} :
def viterbi(obs, states, S, T, E):
    lenseq= len(obs)
    lenstate= len(states)
   
    #print(lenseq)
 
    viterbiMatrix = [[0.0 for x in range(lenseq)] for y in range(lenstate)]
    backtrack = [[-1 for x in range(lenseq)] for y in range(lenstate)]
       
    ind=0
    
    #Base Case when only one 
    #T1 [i,1] <--Start[i]*Emission[i,y1] where y1 is first observation
    #T2 [i,1] == 0
    
    for state in states:
        viterbiMatrix[ind][0] = S[state] + E[state][obs[0]]
        ind=ind+1

    #for each observation 2 to T
    #for each state 1 to k 
    
    #Note to self use + rather than * for math logs 
    for column in range(1, lenseq):
        for row in range(lenstate):
           # print(states[0],states[row])
            #print(row)
            
            #first max [T1[0],i-1] + Transition[k to j + E[ from state j observing obs i]
            max_Val = viterbiMatrix[0][column - 1] + T[states[0]][states[row]] + E[states[row]][obs[column]]
            path = 0 #traceback path
            
            #max value (of the statement )
            for i in range(1, lenstate):
                current_value = viterbiMatrix[i][column - 1] + T[states[i]][states[row]] + E[states[row]][obs[column]]
                if (max_Val < current_value):
                    max_Val = current_value
                    path = i 
                    
            # Inputing the max_value calculated with the path for trackback
            viterbiMatrix[row][column] = max_Val
            backtrack[row][column] = path

    Tr = []
    max_i = 0
    current_Max = 0        
        
    
    #zt < --- argmax T1[k,T]  get most probable state
    for i in range(0, lenstate):
        if (current_Max < viterbiMatrix[i][lenseq - 1]):
            current_Max = viterbiMatrix[i][lenseq - 1]
            max_i = i    
            
    # traceback i ← T,T-1,...,2
    Tr.append(states[max_i])
    new_i= max_i
    for column in range(lenseq-2, 0, -1):
        new_i = backtrack[new_i][column]
        Tr.append(states[new_i])
        
    Tr.reverse()

    #Run Question Functions 
    countHPM(Tr,obs)
    largestfraction(Tr,obs)
    lengthDistributionCal(Tr,obs)
    frequencyCalculation(Tr,obs)
 
    return Tr

#=======================================[Run]=========================================#
def example(seq):
    return viterbi(seq,S,start,T,E)


#A=example(fasta_sequences[1].seq._data) #testing singles

a=[] #collects all the results together

#for all the fasta sequences run viterbi
for i in range(0, len(fasta_sequences)):
    a.append(example(fasta_sequences[i].seq._data))
    
    



#=======================================[Question Functions]=========================================#
    
#Counts the maximum number of hydrophobics
def countHPM(V,s):
    global maxva;
    global maxva_seq
    H=0
    M=0
    P=0
    largesthydrophobics=0
    for state in V:
        if (state == "Hydrophobic"):
            H=H+1
            if (largesthydrophobics<=H):
                largesthydrophobics=H
        if(state != "Hydrophobic"):
            H=0;
    
    if(maxva<largesthydrophobics):
        maxva=largesthydrophobics
        maxva_seq = s
    print(largesthydrophobics,P,M)
 
#Question 1 part c 2 counts the maximum fraction of mixed
def largestfraction(V,s):
    global maxva;
    global L_mixed_frac
    global L_mixed_frac_seq
    global total_M;
    global total_P;
    global total_H;
    
    H=0
    M=0
    P=0
    mixedfraction = 0;
    for state in V:
        if (state == "Hydrophobic"):
            H=H+1
        if(state == "Hydrophilic"):
            P=P+1;
        if(state == "Mixed"):
            M=M+1
        
    total_M = total_M+M;
    total_P=  total_P+P;
    total_H=  total_H+H;
    
    mixedfraction = M/(M+H+P)
    
    if(mixedfraction>L_mixed_frac):
        L_mixed_frac = mixedfraction
        L_mixed_frac_seq = s
        
        

# Question 1 part C 3 Calculates the lengeth distributions of all the H,P M 
def lengthDistributionCal(V,s):
    global h_dist;
    global p_dist;
    global m_dist;
    
    H=0
    M=0
    P=0
    
    for i in range (0,len(V)):
        if (V[i] == "Hydrophobic"):
            H=H+1
            if( i< len(V)-1 and len(V)!=0):
                if(V[i+1] !="Hydrophobic"):
                    h_dist[H]=h_dist[H]+1;
                    H=0;
    
    for i in range (0,len(V)):   
            if (V[i] == "Hydrophilic"):
                M=0;
                H=0;
                P=P+1
            if( i< len(V)-1 and len(V)!=0):
                if(V[i+1] !="Hydrophilic"):
                    p_dist[P]=p_dist[P]+1;
                    P=0;

    for i in range (0,len(V)): 
        if (V[i] == "Mixed"):
            H=0;
            P=0;
            M=M+1
            if( i< len(V)-1 and len(V)!=0):
                if( V[i+1] !="Mixed"):
                    m_dist[M]=m_dist[M]+1;
                    M=0
          
# Question 1 Part C4 Calcluates the frequnny of each H,P,M Amino Acids as it occurs 
def frequencyCalculation(V,s):
    global h_freq
    global m_freq
    global p_freq
    
    
    #for each hydrophobic for each Amino acid +1 
    for i in range(0,len(V)):
        if (V[i] == "Hydrophobic"):
            h_freq[s[i]]=h_freq[s[i]]+1
        if (V[i] == "Hydrophilic"):
            p_freq[s[i]]=p_freq[s[i]]+1
        if (V[i] == "Mixed"):
            m_freq[s[i]]=m_freq[s[i]]+1
            


#testings Used spyder to get data and get individual points 
print("Hydrophobic Frequnecy" + str(h_freq))
print("Hydrophilic Frequnecy" + str(p_freq))
print("Mixed Frequnecy " + str(m_freq))

print(h_dist[1][0])

print(fasta_sequences[403].id)