"""
Created on Sat Jul  6 19:33:29 2019

@author: chonpca
"""


import os
import math
import time

import random
import itertools
#import Bio
import numpy as np
import pandas as pd


def documentmove():
    file_move=[]
    f_list = os.listdir(os.getcwd())
    for f in f_list:
        if os.path.splitext(f)[1]  != '.py' and os.path.splitext(f)[1]  != '.txt':
            file_move.append(f)
    for j in file_move:
        os.remove(j)



def NUPACK_documentgeneration(i,j,type):
    if type ==1:        
        f = open('multi_' + i + '_' + j + '.in','w')
        f.write('2')
        f.write('\n')                                          

        f.write(i)
        f.write('\n')
        f.write(j)
        f.write('\n')
        f.write('1 2')
        f.close()
        
    elif type ==2:
        f = open('single_' + i + '.in','w')
        f.write(i)
        f.close()
        


        
def NUPACK_mfe(i,j,type):
    #global code
    if type == 1:
        code= "mfe -multi -material rna -degenerate  " + 'multi_' + i + '_' +j 
        os.system(code)
        time.sleep(0.05)
    elif type == 2:
        code= "mfe -material rna  -degenerate  " + 'single_' + i
        os.system(code)
        time.sleep(0.05)
        
    #print(code)




def NUPACK_Gread(i,j,type):
    #print(type)
    #print(i)
    
    G=0
    
    if type ==1:
        
        path = 'multi_' + i + '_' + j + '.mfe'
        f = open(path,'r')
        counter = 0
        for j in f:
            
            counter = counter + 1
            if counter == 15:
                #print(j)
                G = j.replace('\n','')
        f.close()
        
    elif type ==2:
        
        path = 'single_' + i + '.mfe'
        f = open(path,'r')
        counter = 0
        for j in f:
            counter = counter + 1
            if counter == 15:
                #print(j)
                G = j.replace('\n','')
        f.close()
        
    G = float(G)
    
    return G




def spacing_energy(i):        #spacing
    if len(i) == 5:
        spacing_energy = 0
    elif len(i) >5:
        spacing_energy = 0.048*(len(i)-5)**2+0.24*(len(i)-5)
    else:
        spacing_energy = 12.2 / ((1+math.e**(2.5*(len(i)-5+2))))**3
        
    return spacing_energy





def mrna_energy(i):    #s2
    NUPACK_documentgeneration(i,'',2)
    NUPACK_mfe(i,'',2)
    j = NUPACK_Gread(i,'',2)
    return j


def mrna_rna_energy(i,ASD='ACCTCCTTA'):
    #s1 
    #print(i)
    NUPACK_documentgeneration(i,ASD,1)   
    #print(i)
    NUPACK_mfe(i,ASD,1)
    #print(i)
    j = NUPACK_Gread(i,ASD,1)
    return j


def standby_energy(i,j,k): # nstart_35ton3_4 + - s2
    a = mrna_energy(i)
    b = mrna_energy(j)
    c = mrna_energy(k)
    return a+b-c


def k_re(i):#d
    re = 0.0075
    if 0<= i <= 25:
        re = 0.0072
    elif i == -4:
        re = 0.022
    elif -25 <= i <= -10:
        re = 0.0072 + 0.0004*(i+10)
        
    return re


def noncoupling_energy(i):   #noncoupling
    try:
        a = mrna_energy(i)
    except:
        a = 0
    
    return a


def coupling_energy(i): #coupling
    return mrna_energy(i)


def F_coupling(): 
    F_coupling=1/(1+0.81*4600)
    return F_coupling
    

def R_re(i):
    re=10*i*4600
    return re


def dGm(sequence,start,end):
        
    for j in range(-1,-1-len(sequence),-3):
        if j > (end-1-len(sequence)):
            #print(i[j-3+1:j+1])
            if sequence[j-3+1:j+1] == 'ATG':
                s1=sequence[0:j + len(sequence) -3+1]
                #s1s1.append(s1)
                #print(j)
                break
                                        #J : g fanxiang
                                          
    
                 
    for k in range(0,len(sequence),3):
        if k > 10:
            if sequence[k:k+3] == 'TAA' or sequence[k:k+3] == 'TGA' or sequence[k:k+3] == 'TAG':
                coupling=sequence[0:k+13]
                noncoupling=sequence[k+13:j+len(sequence)-2+35+1]
                break
                                        
        
    dd = j+len(sequence)-2 - k -3  # from zero
    
    #dd.append(d)
    
    
    
    
    
    s2 = sequence[0:j+len(sequence)-2 +35 +1]
    sd = sequence[start-1:end]
    standby=sequence[start-1-4:start-1]
    nstart_35ton3_4 = sequence[0:start-1-4]
    n3tonstart__35 = sequence[start-1:j+len(sequence)-2+35+1]
    
    
    
    spacing =sequence[end-1+1:j+len(sequence)-2]    
    
    #spacingspacing.append(spacing)
    
    #print(s1)
    

    '''
    couplingcoupling.append(coupling)
    noncouplingnoncoupling.append(noncoupling)
    s1s1.append(s1)
    s2s2.append(s2)
    standbystandby.append(standby)
    nstart_35ton3_4nstart_35ton3_4.append(nstart_35ton3_4)
    n3tonstart__35n3tonstart__35.append(n3tonstart__35)
    '''
    
    
    A =  k_re(dd)
    #print(A)
    B = R_re(A)
    #print(B)
    C = mrna_rna_energy(s1)
    #print(C)
    D = spacing_energy(spacing)
    #print(D)
    E = standby_energy(nstart_35ton3_4, n3tonstart__35, s2)
    #print(E)
    F = noncoupling_energy(noncoupling)
    #print(F)
    G = coupling_energy(coupling)
    #print(G)
    H = F_coupling()
    #print(H) 
    #print(s1)
    #print(mrna_energy(s1))
    dG2 = math.e**(-0.45*(C + D -1.194 + E - F - G*H))
    #print(dG2)
    dG = B + dG2
    #print(dG)
    
    #energy_total.append(dG)
    
    
    #counter = counter + 1
    
    #print(counter/len(sequence))    
    #documentmove()
    

    return dG



def RCBS(sequence):
    
    #sequence = 'ttagagtacttctggtgccagagagataattttcatgaacttctcactacgaagctcacgagttaccggcccaaaaatacgcgtaccgataggctgctcgctgttgttgttcagaagaacacaagcattaccatcgaagcgaatgacagaaccgtccgggcgacgaacacccttcttggtgcgcaccactaccgccttcagcacatcaccttttttgaccttaccacgcggaattgcttctttgatggtgatcttgatgatgtcgcctacgcctgcgtagcgacggtgcgagccacccagaaccttgatacacattacgcgacgtgcaccggagttgtcggcgacgttcagcatagtctgttcttggatcat'
    total_coden = len(sequence)/3
    
    a = set()
    for i in range(0,len(sequence),3):
        coden = sequence[i:i+3]
        a.add(coden)
    b = list(a)
    
    #f(xyz)
    frequence = []
    for i in b:
        counter = 0
        for j in range(0,len(sequence),3):
            if i == sequence[j:j+3]:
                counter = counter + 1
        frequence.append(counter)
    
    
    #dxyz
    dxyzs = []
    
    for i in b:
        #f123xyz
        counter = 0
        f123xyz=[]
        for j in i:
            f123xyz_counter = 0
            for k in range(0,len(sequence),3):
                if j == sequence[k:k+3][counter]:
                    f123xyz_counter = f123xyz_counter + 1            
            f123xyz.append(f123xyz_counter)    
            counter = counter + 1
        dxyz= (((frequence[b.index(i)]/total_coden)-((f123xyz[0]/total_coden)*(f123xyz[1]/total_coden)*(f123xyz[2]/total_coden)))/(((f123xyz[0]/total_coden)*(f123xyz[1]/total_coden)*(f123xyz[2]/total_coden))))
        dxyzs.append(dxyz)
    
    RCBS = 1    
    for i in range(0,len(sequence),3):    
        RCBS = RCBS*(1+dxyzs[b.index(sequence[i:i+3])])

    return (RCBS**(1/total_coden)-1)


def RNAStructure_documentgeneration(i):        
    f = open(i + '.fasta','w')
    f.write('>'+i)
    f.write('\n')
    f.write(i)
    f.close()
    f = open(i + '.ct','w')
    f.close()


    
def RNAStructure_mfe(i):

    code= "./Fold " + i  + '.fasta' + ' ' + i  + '.ct'
    os.system(code)



def TUN_H_SD_structure(k,start_pos,end_pos):  #k STA_length
    
    path = k + '.ct'               
    f = open(path,'r')
    counter = 0
    connect = 0
    zero = 0
    
    for i in f:
        
        if counter !=0:
            
            i = i.split('\n')[0].split()
            link_pos = int(i[-2])


            if counter < start_pos:
                if link_pos > end_pos:
                    break
                
            elif start_pos <= counter <= end_pos:
                if link_pos < end_pos:
                    if link_pos == 0:
                        zero += 1
                    connect += 1
                    
                else:
                    break
                
            elif counter > end_pos:
                if link_pos < start_pos:
                    break
        counter = counter + 1


        
        if counter > len(k):
            break
        
    if (connect== 9) and not(zero == 9):
        #print(k)
        #print(zero)
        #kk.append(k)
        #start.append(start_pos)
        #end.append(end_pos)
        #cover.append(9-zero)
        return 0
    else:
        return 1
        
    f.close()



def TUN_M_SD_structure(k,start_pos,end_pos,ATG_pos):
    
    path = k + '.ct'               
    f = open(path,'r')
    counter = 0
    connect = 0
    zero = 0
    
    for i in f:
        if counter != 0:
            
            i = i.split('\n')[0].split()
            link_pos = int(i[-2])


            if counter < start_pos:
                if link_pos > ATG_pos:
                    break
                
            elif start_pos <= counter <= end_pos:
                if link_pos < ATG_pos:
                    if link_pos == 0:
                        zero += 1
                    connect += 1
                    
                else:
                    break
                
            elif  end_pos< counter < ATG_pos:
                if link_pos > ATG_pos:
                    break
        
            elif counter >= ATG_pos:
                if link_pos < ATG_pos:
                    break
        
        counter = counter + 1


        
        if counter > len(k):
            break
        
    if (connect== 9) and not(zero == 9):
        #print(k)
        #print(zero)
        #kk.append(k)
        #start.append(start_pos)
        #end.append(end_pos)
        #cover.append(9-zero)
        return 0
    else:
        return 1
        
    f.close()



def TUN_L_SD_structure(k,start_pos,end_pos):
    
    path = k + '.ct'               
    f = open(path,'r')
    counter = 0
    connect = 0
    zero = 0
    
    for i in f:
        if counter != 0:
            
            i = i.split('\n')[0].split()
            link_pos = int(i[-2])


            if start_pos <= counter <= end_pos:
                if link_pos < end_pos and link_pos > end_pos:
                    if link_pos == 0:
                        zero += 1
                connect += 1
        counter = counter + 1
        
        if counter > len(k):
            break
        
    if (connect== 9) and not(zero == 9):
        #print(k)
        #print(zero)
        #kk.append(k)
        #start.append(start_pos)
        #end.append(end_pos)
        #cover.append(9-zero)
        return 0
    else:
        return 1
        
    f.close()



    
def Vienna_documentgeneration(i,j):
    f = open(i+'.fasta','w')
    f.write('>'+i)
    f.write('\n')
    f.write(j)
    f.close()



def Vienna_mfe(i):
    #time.sleep(1)
    code= "RNAfold " + '-i ' + i + '.fasta' + ' -o'
    
    os.system(code)
    #global a
    time.sleep(2)
    f = open(i+'.fold','r')
    aa = f.readlines()
    f.close()
    #print(a)
    #print(a[0])
    #time.sleep(1)
    
    
    ff = open(i+'.fold','w')
    ff.write('>')
    #f.write('\n')
    ff.write(aa[1])
    #f.write('\n')
    ff.write(aa[2])
    #f.write('\n')
    ff.close()
    #time.sleep(1)
    


    
def Vienna_RNAStructure_dot2ct(i):
    
    code= "./dot2ct " + i  + '.fold' + ' ' + i  + '.ct'  
    #print(cide)
    os.system(code)
    time.sleep(2)


    
def Vienna_Gread(i):
    
    f = open(i+'.ct','r')
    a = f.readlines()
    f.close()
    G = str(a[0].split()[1].replace('(','').replace(')',''))       
    f.close()
    return G


def RNA_Structure_Vienna_Structure(k,i,j,sty):
    
    f = open(k +'.ct','r')
    counter = 0
    
    #global structure
    
    if sty == 1:

        structure = np.zeros((len(i),len(i)))        
        
        for i in f:
            if counter != 0 :
                    i = i.split('\n')[0].split()
                    link_pos = int(i[-2])
                    if link_pos !=0:                        
                        if link_pos != 0:                         
                            structure[counter-1,link_pos-1] = 1
                            structure[link_pos-1,counter-1] = 1             
            counter = counter + 1

    elif sty == 2:
        
        ori = j.shape[0] 
        
        structure = np.zeros((ori,ori))
        
        for i in f:
            if counter !=0:
                i = i.split('\n')[0].split()
                link_pos = int(i[-2])

                if link_pos != 0:
                    
                    if link_pos > ori:
                        structure[:,counter-1] = 10**10
                        structure[counter-1,:] = 10**10
                    
                    else:
                        if 1 in j[counter-1,:]:   
                            
                            link_pos_ori = j[counter-1,:].tolist().index(1)
                            
                            if link_pos == link_pos_ori:
                                structure[counter-1,link_pos-1] = 1
                                structure[link_pos-1,counter-1] = 1
                            
                            elif link_pos != link_pos_ori:
                                structure[counter-1,link_pos-1] = 10**10
                                structure[link_pos-1,counter-1] = 10**10
                        
                        else:
                            structure[counter-1,link_pos-1] = 10**10
                            structure[link_pos-1,counter-1] = 10**10
                    
                else:
                    if 1 in j[counter-1,:]:
                        link_pos_ori = j[counter-1,:].tolist().index(1)
                        structure[counter-1,link_pos_ori] = 10**10
                        structure[link_pos_ori,counter-1] = 10**10
            
            
            counter = counter + 1            
            
            if counter > j.shape[0]:
                break
            
    return structure


def Distance(i,j):
    return ((i-j)**2).sum()  



def STA_H_SD_structure(Ribo,k,start_pos,end_pos):  #k:STA_length Ribo:Ribo sequence
    
    path = str(k) + '.ct'               
    f = open(path,'r')
    counter = 0
    connect = 0
    zero = 0
    
    for i in f:
        if counter > len(Ribo) + k:
            
            i = i.split('\n')[0].split()
            link_pos = int(i[-2])


            if counter < start_pos:
                if link_pos > end_pos:
                    break
                
            elif start_pos <= counter <= end_pos:
                if link_pos < end_pos:
                    if link_pos == 0:
                        zero += 1
                    connect += 1
                    
                else:
                    break
                
            elif counter > end_pos:
                if link_pos < start_pos:
                    break
        counter = counter + 1


        
        if counter > len(Ribo) + k + 38:
            break
        
    if (connect== 9) and not(zero == 9):
        #print(k)
        #print(zero)
        #kk.append(k)
        #start.append(start_pos)
        #end.append(end_pos)
        #cover.append(9-zero)
        return 0
    else:
        return float("inf")
        
    f.close()



def STA_M_SD_structure(Ribo,k,start_pos,end_pos,ATG_pos):
    
    path = str(k) + '.ct'               
    f = open(path,'r')
    counter = 0
    connect = 0
    zero = 0
    
    for i in f:
        if counter > len(Ribo) + k:
            
            i = i.split('\n')[0].split()
            link_pos = int(i[-2])


            if counter < start_pos:
                if link_pos > end_pos:
                    break
                
            elif start_pos <= counter <= end_pos:
                if link_pos < ATG_pos:
                    if link_pos == 0:
                        zero += 1
                    connect += 1
                    
                else:
                    break
                
            elif counter > end_pos:
                if link_pos < ATG_pos:
                    break
        counter = counter + 1


        
        if counter > len(Ribo) + k + 38:
            break
        
    if (connect== 9) and not(zero == 9):
        #print(k)
        #print(zero)
        #kk.append(k)
        #start.append(start_pos)
        #end.append(end_pos)
        #cover.append(9-zero)
        return 0
    else:
        return float("inf")
        
    f.close()



def STA_L_SD_structure(Ribo,k,start_pos,end_pos):
    
    path = str(k) + '.ct'               
    f = open(path,'r')
    counter = 0
    connect = 0
    zero = 0
    
    for i in f:
        if counter > len(Ribo) + k:
            
            i = i.split('\n')[0].split()
            link_pos = int(i[-2])


            if start_pos <= counter <= end_pos:
                if link_pos < end_pos:
                    if link_pos == 0:
                        zero += 1
                connect += 1
        counter = counter + 1
        
        if counter > len(Ribo) + k + 38:
            break
        
    if (connect== 9) and not(zero == 9):
        #print(k)
        #print(zero)
        #kk.append(k)
        #start.append(start_pos)
        #end.append(end_pos)
        #cover.append(9-zero)
        return 0
    else:
        return float("inf")
        
    f.close()


    
def STA_Length(Ribo,ori_GOI,poly_GOI):    #poly_GOI:
    
    try:    
        Vienna_documentgeneration('Ribo',Ribo)
        Vienna_mfe('Ribo')
        Vienna_RNAStructure_dot2ct('Ribo')
        Ribo_mat = RNA_Structure_Vienna_Structure('Ribo',Ribo,'',1)
        
        Vienna_documentgeneration('RGOI',Ribo+ori_GOI)
        Vienna_mfe('RGOI')
        Vienna_RNAStructure_dot2ct('RGOI')
        Ribo_GOI_mat = RNA_Structure_Vienna_Structure('RGOI',Ribo+ori_GOI,Ribo_mat,2)
    
    except:
        Vienna_documentgeneration('Ribo',Ribo)
        Vienna_mfe('Ribo')
        Vienna_RNAStructure_dot2ct('Ribo')
        Ribo_mat = RNA_Structure_Vienna_Structure('Ribo',Ribo,'',1)
        
        Vienna_documentgeneration('RGOI',Ribo+ori_GOI)
        Vienna_mfe('RGOI')
        Vienna_RNAStructure_dot2ct('RGOI')
        Ribo_GOI_mat = RNA_Structure_Vienna_Structure('RGOI',Ribo+ori_GOI,Ribo_mat,2)
        
        
        
        

        
    
    try:
        
        lengths=[]
        scores=[]
        kkk=[]
        
        if Distance(Ribo_GOI_mat,Ribo_mat) != 0:
            
            for i in range(6,len(ori_GOI),3):
                
                Ribo_STA = Ribo + ori_GOI[0:i] + poly_GOI
                
                Vienna_documentgeneration(str(i),Ribo_STA)
                Vienna_mfe(str(i))
                Vienna_RNAStructure_dot2ct(str(i))
                Ribo_STA_mat = RNA_Structure_Vienna_Structure(str(i),Ribo_STA,Ribo_mat,2)
                #global Ribo_mat
                #global Ribo_GOI_mat
                #global Ribo_STA
                score = (Distance(Ribo_STA_mat,Ribo_mat))/(Distance(Ribo_GOI_mat,Ribo_mat))
                score = score + STA_L_SD_structure(Ribo,i,len(Ribo)+i+26,len(Ribo)+i+34)
                
                
                lengths.append(i)
                scores.append(score)
                kkk.append(Distance(Ribo_STA_mat,Ribo_mat))
                
                print(i)
                print(score)
                #print(Distance(Ribo_STA_mat,Ribo_GOI_mat))
                
                #print(Distance(Ribo_STA_mat,Ribo_mat))
                
                
        elif Distance(Ribo_mat,Ribo_GOI_mat) == 0:
            
            
            for i in range(6,len(ori_GOI),3):
                
                Ribo_STA = Ribo + ori_GOI[3:i] + poly_GOI                
                Vienna_documentgeneration(str(i),Ribo_STA)
                Vienna_mfe(str(i))
                Vienna_RNAStructure_dot2ct(str(i))
                Ribo_STA_mat = RNA_Structure_Vienna_Structure(str(i),Ribo_STA,Ribo_mat,2)            
                score = (Distance(Ribo_STA_mat,Ribo_mat))
                score = score + STA_L_SD_structure(Ribo,i,len(Ribo)+i+26,len(Ribo)+i+34)
                
                lengths.append(i)
                scores.append(score)
                
                
                print(i)
                print(score)
        
        
                
        
        return [lengths,scores,kkk]

    
    except:
        
        lengths=[]
        scores=[]
        kkk=[]
        
        
        
        
        if Distance(Ribo_GOI_mat,Ribo_mat) != 0:
            
            for i in range(6,len(ori_GOI),3):
                
                Ribo_STA = Ribo + ori_GOI[0:i] + poly_GOI
                Vienna_documentgeneration(str(i),Ribo_STA)
                Vienna_mfe(str(i))
                Vienna_RNAStructure_dot2ct(str(i))
                Ribo_STA_mat = RNA_Structure_Vienna_Structure(str(i),Ribo_STA,Ribo_mat,2)
                #global Ribo_mat
                #global Ribo_GOI_mat
                #global Ribo_STA
                score = ((Distance(Ribo_STA_mat,Ribo_mat))/(Distance(Ribo_GOI_mat,Ribo_mat)))
                score = score + STA_L_SD_structure(Ribo,i,len(Ribo)+i+26,len(Ribo)+i+34)
                lengths.append(i)
                scores.append(score)
                #kkk.append(Distance(Ribo_STA_mat,Ribo_mat))
                print(i)
                print(score)
                #print(Distance(Ribo_STA_mat,Ribo_GOI_mat))
                #print(Distance(Ribo_STA_mat,Ribo_mat))
                
                
        elif Distance(Ribo_mat,Ribo_GOI_mat) == 0:
            
            
            for i in range(6,len(ori_GOI),3):
                
                Ribo_STA = Ribo + ori_GOI[3:i] + poly_GOI
                
                Vienna_documentgeneration(str(i),Ribo_STA)
                Vienna_mfe(str(i))
                Vienna_RNAStructure_dot2ct(str(i))
                Ribo_STA_mat = RNA_Structure_Vienna_Structure(str(i),Ribo_STA,Ribo_mat,2)            
                score = (Distance(Ribo_STA_mat,Ribo_mat))
                score = score + STA_L_SD_structure(Ribo,i,len(Ribo)+i+26,len(Ribo)+i+34)
                lengths.append(i)
                scores.append(score)
                
                print(i)
                print(score)
        
        
                


    return [lengths,scores]



def Vienna_pos_lock_documentgeneration(i,j,k):    #k: dot-point
    f = open(i+'.fasta','w')
    f.write('>'+i)
    f.write('\n')
    f.write(j)
    f.write('\n')
    f.write(k)
    f.close()
    




def Vienna_pos_lock_mfe(i):
    code= "RNAfold " + '-C ' + i + '.fasta' + ' -o'
    os.system(code)    




def DPH_Structure(i,start_pos):
    
    loop = [0]       
       
    f = open(i+'.ct','r')
    
    counter = 0
    ss = 0
    
    for i in f:
        if counter != 0:
            i = i.split('\n')[0].split()
            #print(i)
            link_pos = int(i[-2])   
            if counter < start_pos and link_pos < start_pos:                
                if (link_pos != 0) and (counter >loop[ss]):
                    ss = ss + 2
                    loop.append(counter)
                    loop.append(link_pos)
        counter = counter + 1
    f.close()
    print(loop)
    
    return loop 
    
    
    
def RiboSwitch(sequence,ligard_domain,start_pos,end_pos):    # sequence: Ribo + ATG + 66bp ; ligard_domain: 0:?
    
    
    ribosome_footprint = 29
    
    #1
    Ginit = mrna_energy(sequence)
    print(Ginit)
    
    #2
    Gligard_structure = ''
    Gligard_structure  = Gligard_structure + ligard_domain
    
    for i in range(len(sequence)-len(ligard_domain)):
        Gligard_structure += '.'
    
    Vienna_pos_lock_documentgeneration('Gligard',sequence,Gligard_structure)
    Vienna_pos_lock_mfe('Gligard')
    Vienna_RNAStructure_dot2ct('Gligard')
    Gligard = Vienna_Gread('Gligard')                   
    print(Gligard)
    
    #3
    
    Gmrna_rna_structure = ''
    for i in range(start_pos-1):
        Gmrna_rna_structure += '.'
    if len(sequence) < ribosome_footprint + len(Gmrna_rna_structure):        
        for i in range(len(sequence)-len(Gmrna_rna_structure)):
            Gmrna_rna_structure += 'x'
    else:
        for i in range(ribosome_footprint):
            Gmrna_rna_structure += 'x'
        for i in range(len(sequence)-len(Gmrna_rna_structure)):
            Gmrna_rna_structure += '.'

    Vienna_pos_lock_documentgeneration('Gmrna_rna',sequence,Gmrna_rna_structure)
    Vienna_pos_lock_mfe('Gmrna_rna')
    Vienna_RNAStructure_dot2ct('Gmrna_rna')
    #Gmrna_rna = float(Vienna_Gread('Gmrna_rna'))
    #print(Gmrna_rna)        
    
    
    #Gmrna_rna += mrna_rna_energy(sequence[start_pos-1:end_pos])
    Gmrna_rna = mrna_rna_energy(sequence)
    print(Gmrna_rna)
    
    
    
    #4
    Gmrna_rna_ligard_structure = ''
    Gmrna_rna_ligard_structure = Gmrna_rna_ligard_structure + ligard_domain
    for i in range(start_pos-len(ligard_domain)-1):
        Gmrna_rna_ligard_structure += '.'
    if len(sequence) < ribosome_footprint + len(Gmrna_rna_ligard_structure):        
        for i in range(len(sequence)-len(Gmrna_rna_ligard_structure)):
            Gmrna_rna_ligard_structure += 'x'
    else:
        for i in range(ribosome_footprint):
            Gmrna_rna_ligard_structure += 'x'
        for i in range(len(sequence)-len(Gmrna_rna_ligard_structure)):
            Gmrna_rna_ligard_structure += '.'
                
    Vienna_pos_lock_documentgeneration('Gmrna_rna_ligard',sequence,Gmrna_rna_ligard_structure)
    Vienna_pos_lock_mfe('Gmrna_rna_ligard')
    Vienna_RNAStructure_dot2ct('Gmrna_rna_ligard')
    #Gmrna_rna_ligard = float(Vienna_Gread('Gmrna_rna_ligard'))
    #Gmrna_rna_ligard += mrna_rna_energy(sequence[start_pos-1:end_pos])
    Gmrna_rna_ligard = mrna_rna_energy(sequence)
    print(Gmrna_rna_ligard)
    
    #5
    
    Gstart = -1.194
    
    #6 
    
    for j in range(-1,-1-len(sequence),-3):
        if j > (end_pos-1-len(sequence)):
            #print(i[j-3+1:j+1])
            if sequence[j-3+1:j+1] == 'ATG':
                #s1=sequence[0:j + len(sequence) -3+1]
                #s1s1.append(s1)
                #print(j)
                break        
                                           #J : g fanxiang
    
    j = j + len(sequence) -2 +1            #from 1
    
    spacing = j - end_pos -1
    #print(spacing)
    spacing =sequence[end_pos-1+1:end_pos-1+1+spacing]  
    print(spacing)
    Gspacing = spacing_energy(spacing)
    
    #print('fuck',Gspacing)
    print(Gspacing)
    
    #7
    
    ribo_loop = DPH_Structure('Gmrna_rna',start_pos)
    
    D = start_pos - ribo_loop[-1] -1
    
    P = ribo_loop[1] -1
    
    H = []
    
    for i in range(1,len(ribo_loop),2):
        H.append(0.5*(ribo_loop[i+1] - ribo_loop[i] +1))
    
    H = np.array(H)
        
    AS = 15 + P + D - H.mean()    
    
    Gdis = 0.038*(AS)**2 - 1.629*(AS) + 17.359
    
    Gunfolding = 0
    
    Gsliding  = 0.2*H.sum()
    
    Gstandby = Gdis + Gunfolding + Gsliding
    print(D,P,H)
    print(Gstandby)
    
    #8
    
    ribo_loop_ligard = DPH_Structure('Gmrna_rna_ligard',start_pos)
    
    D_ligard = start_pos - ribo_loop_ligard[-1] -1
    
    P_ligard = ribo_loop_ligard[1] -1
    
    H_ligard = []
    
    for i in range(1,len(ribo_loop_ligard),2):
        H_ligard.append(0.5*(ribo_loop_ligard[i+1] - ribo_loop_ligard[i] +1))
    
    H_ligard = np.array(H_ligard)
        
    AS_ligard = 15 + P_ligard + D_ligard - H_ligard.mean()    
    
    Gdis_ligard = 0.038*(AS_ligard)**2 - 1.629*(AS_ligard) + 17.359
    
    Gunfolding_ligard = 0
    
    Gsliding_ligard  = 0.2*H_ligard.sum()
        
    G_ligard_standby = Gdis_ligard + Gunfolding_ligard + Gsliding_ligard
    print(D_ligard,P_ligard,H_ligard)
    print(G_ligard_standby)
    
    E_on = float(G_ligard_standby) + float(Gspacing) + float(Gstart) + float(Gmrna_rna_ligard) - float(Gligard)
    
    E_off = float(Gstandby) + float(Gspacing) + float(Gstart) + float(Gmrna_rna) - float(Ginit)
    
    print(E_on,E_off)
    
    R_on = math.e**(-0.45*E_on)
    
    R_off = math.e**(-0.45*E_off)
    
    
    return [R_on,R_off]




'''
def PSO(function,sty,particle_number,max_iter,pra_spc,hig,low):
    global X
    global V
    w = 0.8    
    c1 = 2     
    c2 = 2     
    r1= 0.6  
    r2=0.3  
    pN = particle_number               
    dim = pra_spc                
    max_iter = max_iter      
    X = np.zeros((pN,dim))         
    V = np.zeros((pN,dim))  
    pbest = np.zeros((pN,dim))    
    gbest = np.zeros((1,dim))  
    p_fit = np.zeros(pN)                

    if sty == 'min':
        
        fit = float('inf')          
        
        #init
        
        for i in range(pN):  
            for j in range(dim):  
                X[i][j] = random.uniform(low[j],hig[j])  
                V[i][j] = random.uniform(0,1)  
            
            pbest[i] = X[i].copy()  
            p_fit[i] = function(X[i])
            
            
        for i in range(pN):
            if(p_fit[i] < fit):
                
                fit = p_fit[i]  
                gbest = X[i].copy()        
            
            
        #fitness = []
        
        for t in range(max_iter):  
            for i in range(pN):
                
                if(function(X[i]) < p_fit[i]):        
                    p_fit[i] = function(X[i])  
                    pbest[i] = X[i].copy()  
                
            for i in range(pN):                
                if(function(X[i]) < fit):    
                    fit = function(X[i])
                    gbest = X[i].copy()
                    print(fit,gbest)
                      
            for i in range(pN):  
                V[i] = w*V[i] + c1*r1*(pbest[i] - X[i]) + c2*r2*(gbest - X[i])  
                X[i] = X[i] + V[i]
                
                for j in range(dim):
                    if X[i][j] >= hig[j] or X[i][j] < low[j]:
                        X[i][j] = random.uniform(low[j],hig[j])
            
            #fitness.append(fit)  
                             
        #return [fit,gbest,fitness]

    elif sty == 'max':
        
        fit = -float('inf')          

        #init

        for i in range(pN):  
            for j in range(dim):  
                X[i][j] = random.uniform(low[j],hig[j])  
                V[i][j] = random.uniform(low[j],hig[j])  
                pbest[i] = X[i].copy()  
                p_fit[i] = function(X[i])  

        for i in range(pN):
            if(p_fit[i] > fit):  
                fit = p_fit[i]  
                gbest = X[i].copy()        
                print(fit)
                print(gbest)

        fitness = []  
        #counter = 1
        for t in range(max_iter):  
            for i in range(pN):

                if(function(X[i]) > p_fit[i]):        
                    p_fit[i] = function(X[i])  
                    pbest[i] = X[i].copy()

            for i in range(pN):                
                if(function(X[i]) > fit):    
                    fit = function(X[i])
                    gbest = X[i].copy()
                    #print(fit)
                    #print(gbest)
                    #print(fun1(gbest)==fit)
            

            for i in range(pN):  
                V[i] = w*V[i] + c1*r1*(pbest[i] - X[i]) + c2*r2*(gbest - X[i])  
                X[i] = X[i] + V[i]
                #print(X[i])

                for j in range(dim):
                    if X[i][j] > hig[j] or X[i][j] < low[j]:
                        X[i][j] = random.uniform(low[j],hig[j])
                #print(X[i])
            fitness.append(fit)  

            #counter = counter + 1


    return [fit,gbest]
'''
'''
def Coden_Rich(species):    # coden rich
    
    Coden_Rich_table = []
'''    



def con(sequence):
    con_mat=[]
    for i in range(len(sequence)):
        if sequence[i] == 'A':
            con_mat.append(random.uniform(0,0.25))
        elif sequence[i] =='G':
            con_mat.append(random.uniform(0.25,0.5))
        elif sequence[i] == 'C':
            con_mat.append(random.uniform(0.5,0.75))
        elif sequence[i] == 'T':
            con_mat.append(random.uniform(0.75,1))
        
    return con_mat
        
        
        
def dis(sequence):
    dis_mat = ''
    for i in range(len(sequence)):
        if 0 <= sequence[i] < 0.25:
            dis_mat += 'A'
        elif 0.25 <= sequence[i] < 0.5:
            dis_mat += 'G'
        elif 0.5 <= sequence[i] < 0.75:
            dis_mat += 'C'
        elif 0.75 <= sequence[i] < 1:
            dis_mat += 'T'
            
    return dis_mat



def sequence_length(d):
    length=[]
    if d >= 0:
        length.append((32-d) + 3*(1-(((32-d)/3) - math.floor(((32-d)/3)))))
        length.append(d)
        length.append(33)
    elif d < 0:
        length.append(35 + 3*(1-(((35-d)/3) - math.floor(((35-d)/3)))))
        length.append(-d-6)
        length.append(33 + 3*(1-(((-d)/3) - math.floor(((-d)/3)))))
    
    return np.array(length)



def par_space():
    
    par_num = 0
    
    lengths=[]
    
    for i in range(-25,26):
        lengths.append(sequence_length(i).sum())
        
    par_num = par_num + int(np.array(lengths).max())
    
    
    par_num = par_num + 4 + 6    
    
    return par_num
    
    #d,sd_start_pos,sd_end_pos，m3_library


    
''' 
def get_ASD_oriSD(species):    #species:
    
    a = pd.read_excel('ASD_library.xlsx',header=None)    
    b = np.array(a)
    c = b.tolist()
    
    
    if spe in c:
        
        
        
    elif:
'''         
    
'''
def sd_generation(ASD='ACCTCCTTA',SD='TAAGGAGGT'):
    
    #nu=['A','G','C','T']
    rbs_all = []
    rbs_select = []
    rbs_index = itertools.product('AGCT',repeat=9)
     
    for i in rbs_index:
        rbs=''
        for j in i:
            rbs = rbs + j
        rbs_all.append(rbs)
        
    #energy_all = []
    energy_select = []
    stand = mrna_rna_energy(ASD,SD)
    for i in rbs_all:
        tmp = mrna_rna_energy(i,ASD)
        if  abs(tmp - stand) <= 0.5:
            rbs_select.append(i)
            energy_selece.append(tmp)
        
    print('SD_generation is OK')    
    return [rbs_select,energy_select][0]    
    
'''   




                
def m3_library_generation(Coden_Rich_table,GOI):
    #global m3_library
    
    m3_library=[]
    
    
    
    #GOI = GOI[3:34]
    
    Met=['ATG']
    Arg=['AGG','CGT','CGG']
    Lys=['AAA']
    Gly=['GGA','GGC']
    Glu=['GAA','GAG']
    Leu=['TTA','TTG','CTA','CTG']
    Phe=['TTC']
    Thr=['ACC','ACA']
    Val=['GTA']
    
    
    A=[['AGG','CGT','CGG'],
       ['AAA'],
       ['GGA','GGC'],
       ['GAA','GAG'],
       ['GAA','GAG'],
       ['TTA','TTG','CTA','CTG'],
       ['TTC'],
       ['ACC','ACA'],
       ['GGA','GGC'],
       ['GTA'],
       ['GTA']]
    
    
    check=[]
    for i in range(3):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    for m in range(4):
                        for n in range(2):
                            for o in range(2):
                                check.append([i,j,k,l,m,n,o])
        
    AUGplus35=[]
    
    for i in check:
        seq=''
        seq += Arg[i[0]]
        seq += 'AAA'
        seq += Gly[i[1]]
        seq += Glu[i[2]] 
        seq += Glu[i[3]]
        seq += Leu[i[4]]
        seq += 'TTC'
        seq += Thr[i[5]]
        seq += Gly[i[6]]
        seq += 'GTAGTA'
        AUGplus35.append(seq)
        
        
    m3_library=AUGplus35
    
    print('m3_library is generation OK')
    return m3_library




def Tuner_goal(x,G_goal=3000,species='',GOI=''):
    
    score = 0
    #print(x)
    d = int(math.floor(x[0]))
    start_pos = int(math.floor(x[1]))
    end_pos = int(math.floor(x[2]))
    m3_number = int(math.floor(x[3]))
    if m3_number == len(m3_library):
        m3_number -= 1
    print(m3_number)
    #print(d)
    #print(start_pos)
    #print(end_pos)
    #print(x)
    
    #Coden_Richs = Coden_Rich(species)

    
    if not( -25<= d <= -10 or d == -4 or 0 <= d <= 25):
        score += 100**150
    else:
        #if d >=0:
        #print('a')
        lengths = sequence_length(d)
        sequence_mat = x[4:4+int(lengths.sum())+6 -33]
        sequence = dis(sequence_mat) + m3_library[m3_number]
        print(sequence)
        
        print(d)
        
        if d < 0:
            

            #print('b')
            #print(sequence[int(lengths[0]):int(lengths[0])+3])
            #print(int(lengths[0]))
            if not(sequence[int(lengths[0]):int(lengths[0])+3] == 'ATG'):
                score += 100**140        
            else:
                #print('c')
                if not(sequence[int(lengths[0])+3+int(lengths[1]):int(lengths[0])+3+int(lengths[1])+3] == 'TAA'):
                    score += 100**130
                else:
                    
                    if not(end_pos < lengths[0] + 1):
                        if end_pos-lengths[0] !=0:
                            
                            score += abs(end_pos-lengths[0])*100**120
                        else:
                            score += 100**119
                    else:
                        #print('d')
                        if not(start_pos - end_pos < 0):
                            if start_pos - end_pos !=0:                        
                                score += (start_pos-end_pos)*100**110
                            else:
                                score += 100**109                  
                        else:
                            
                            #print('e')
                            Start_first_pos = []        # from 1
                            Stop_first_pos = []         # from 1
                            
                            if sequence[-3::] == 'ATG':                    
                                Start_first_pos.append(len(sequence)-2) 
                                
                            for j in range(-1,-1-len(sequence)-1,-3):
                                if sequence[j-3+1:j+1] == 'ATG':
                                    Start_first_pos.append(j+len(sequence)-2+1)
                                    
                                                                #J : g fanxiang
                                                         
                            for k in range(0,len(sequence),3):
                                #if k > 10:
                                if sequence[k:k+3] == ('TAA' or 'TGA' or 'TAG'):
                                        #coupling=sequence[0:k+13]
                                        #noncoupling=sequence[k+13:j+len(sequence)-2+35+1]
                                        #break
                                        Stop_first_pos.append(k+1)
                                
                            #d = j+len(sequence)-2 - k -3  # from zero\
                            
                            
                            if not(len(Start_first_pos) ==1):
                                score += 100**100
                            
                            else:
                                #print('f')
                                #print(Start_first_pos[0])
                                #print(sequence[Start_first_pos[0]:Start_first_pos[0]+3])
                                if not(Start_first_pos[0] == int(lengths[0]) + 1 and len(Stop_first_pos) == 1):
                                    score += 100**90
                                else:
                                    #print('g')
                                    print(start_pos,end_pos)
                                    #if end_pos
                                    if not(Stop_first_pos[0] == int(lengths[0]) + 3 + int(lengths[1]) + 1 and start_pos - end_pos == -8):
                                        score += abs(start_pos-end_pos + 8)*100**80
                                    else:
                                                                                
                                        Start_first_pos_reverse = []        # from 1
                                        Stop_first_pos_reverse = []         # from 1
                                        
                                        if sequence[-3::] == ('TAA' or 'TGA' or 'TAG'):                    
                                            Stop_first_pos_reverse.append(len(sequence)-2) 
                                            
                                        for j in range(-1,-1-len(sequence)-1,-3):
                                            if sequence[j-3+1:j+1] == ('TAA' or 'TGA' or 'TAG'):
                                                Stop_first_pos_reverse.append(j+len(sequence)-2+1)
                                            if j+len(sequence)-2+1 >=  Start_first_pos[0]:
                                                break
                                                                            #J : g fanxiang
                                                                     
                                        for k in range(0,len(sequence),3):
                                            #if k > 10:
                                            if sequence[k:k+3] == 'ATG':
                                                    #coupling=sequence[0:k+13]
                                                    #noncoupling=sequence[k+13:j+len(sequence)-2+35+1]
                                                    #break
                                                    Start_first_pos_reverse.append(k+1)                                        
                                            if k+1 >= Stop_first_pos[0]:
                                                break
                                        #print(Start_first_pos_reverse)
                                        #print(Stop_first_pos_reverse)
                                            
                                        if not(len(Start_first_pos_reverse) == 0 and len(Stop_first_pos_reverse) ==0):
                                            score += abs(len(Start_first_pos_reverse)) + abs(len(Stop_first_pos_reverse))*100**70
                                        else:
                                            #print('h')
                                            if not(len(sequence[start_pos-1:end_pos]) == 9):
                                                score += abs(len(sequence[start_pos-1:end_pos])-9)*100**60
                                            else:
                                                #print('i')
                                                ASD='ACCTCCTTA'
                                                SD='TAAGGAGGT'
                                                stand = mrna_rna_energy(SD,ASD)
                                                print(sequence[start_pos-1:end_pos])
                                                tmp = mrna_rna_energy(sequence[start_pos-1:end_pos],ASD)
                                                if not(abs(stand-tmp) <=1):
                                                    score += abs(abs(stand-tmp)-1)*100**50
                                                    
                                                
                                                #if not(sequence[start_pos-1,end_pos] in SD_library):
                                                #    return float('inf')
                                                else:
    
                                                    print('j')                                        
                                                    RNAStructure_documentgeneration(sequence)
                                                    RNAStructure_mfe(sequence)
                                                    b = TUN_M_SD_structure(sequence,start_pos,end_pos,Start_first_pos[0])
                                                    
                                                    if b == 1:
                                                        score += 100**40
                                                    else:
                                                        print('k')
                                                        a = dGm(sequence,start_pos,end_pos)
                                                        ori_150bp = 'ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGTACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACT'
                                                        sfgfp_mins_36bp = 'CCTATTCTGGTGGAACTGGATGGTGATGTCAACGGTCATAAGTTTTCCGTGCGTGGCGAGGGTGAAGGTGACGCAACTAATGGTAAACTGACGCTGAAGTTCATCTGTACTACTGGTAAACTGCCGGTACCTTGGCCGACTCTGGTAACGACGCTGACTTATGGTGTTCAGTGCTTTGCTCGTTATCCGGACCATATGAAGCAGCATGACTTCTTCAAGTCCGCCATGCCGGAAGGCTATGTGCAGGAACGCACGATTTCCTTTAAGGATGACGGCACGTACAAAACGCGTGCGGAAGTGAAATTTGAAGGCGATACCCTGGTAAACCGCATTGAGCTGAAAGGCATTGACTTTAAAGAAGACGGCAATATCCTGGGCCATAAGCTGGAATACAATTTTAACAGCCACAATGTTTACATCACCGCCGATAAACAAAAAAATGGCATTAAAGCGAATTTTAAAATTCGCCACAACGTGGAGGATGGCAGCGTGCAGCTGGCTGATCACTACCAGCAAAACACTCCAATCGGTGATGGTCCTGTTCTGCTGCCAGACAATCACTATCTGAGCACGCAAAGCGTTCTGTCTAAAGATCCGAACGAGAAACGCGATCATATGGTTCTGCTGGAGTTCGTAACCGCAGCGGGCATCACGCATGGTATGGATGAACTGTACAAATGA'
                                                        a = a*RCBS(ori_150bp[3::]+sequence[:Stop_first_pos[0]-1])
                                                        a = a*RCBS(sequence[Start_first_pos[0]+3-1::]+sfgfp_mins_36bp[:-3])
                                                                                                                
                                                        print(a)
                                                        score = 1/a
        
        elif d>=0:
            #a

            #print('b')
            #print(sequence[int(lengths[0])])
            #print(sequence[int(lengths[0]):int(lengths[0])+3])
            if not(sequence[int(lengths[0]):int(lengths[0])+3] == 'TAA'):
                score += 100**140        
            else:
                #print('c')
                if not(sequence[int(lengths[0])+3+int(lengths[1]):int(lengths[0])+3+int(lengths[1])+3] == 'ATG'):
                    score += 100**130
                else:
                    
                    if not(end_pos < lengths[0] +3 + lengths[1] + 1):
                        if end_pos -lengths[0] -3 -lengths[1] !=0:
                            
                            score += (end_pos-lengths[0])*100**120
                        else:
                            score += 100**119
                    else:
                        #print('d')
                        if not(start_pos - end_pos <0):
                            if start_pos-end_pos !=0:
                                
                                score += (start_pos - end_pos)*100**110
                            else:
                                score += 100**109
                        else:
                            
                            #print('e')
                            Start_first_pos = []        # from 1
                            Stop_first_pos = []         # from 1
                            
                            if sequence[-3::] == 'ATG':                    
                                Start_first_pos.append(len(sequence)-2) 
                                
                            for j in range(-1,-1-len(sequence)-1,-3):
                                if sequence[j-3+1:j+1] == 'ATG':
                                    Start_first_pos.append(j+len(sequence)-2 +1)   #from 1
                                    
                                                                #J : g fanxiang 
                                                         
                            for k in range(0,len(sequence),3):
                                #if k > 10:
                                if sequence[k:k+3] == ('TAA' or 'TGA' or 'TAG'):
                                        #coupling=sequence[0:k+13]
                                        #noncoupling=sequence[k+13:j+len(sequence)-2+35+1]
                                        #break
                                        Stop_first_pos.append(k+1)   #from 1
                                
                            #d = j+len(sequence)-2 - k -3  # from zero\
                            
                            
                            if not(len(Stop_first_pos) ==1):
                                score += 100**100
                            
                            else:
                                #print('f')
                                #print(Stop_first_pos[0])
                                #print(sequence[Stop_first_pos[0]:Stop_first_pos[0]+3])
                                if not(Stop_first_pos[0] == int(lengths[0]) +1  and len(Start_first_pos) == 1):
                                    score += 100**90
                                else:
                                    #print('g')
                                    print(start_pos,end_pos)
                                    #if end_pos
                                    if not(Start_first_pos[0] == int(lengths[0]) + 3 + int(lengths[1]) +1 and start_pos - end_pos == -8):
                                        score += abs(start_pos-end_pos + 8)*100**80
                                    else:
                                        
                                        
                                        
                                        Start_first_pos_reverse = []        # from 1
                                        Stop_first_pos_reverse = []         # from 1
                                        
                                        if sequence[-3::] == ('TAA' or 'TGA' or 'TAG'):                    
                                            Stop_first_pos_reverse.append(len(sequence)-2) 
                                            
                                        for j in range(-1,-1-len(sequence)-1,-3):
                                            if sequence[j-3+1:j+1] == ('TAA' or 'TGA' or 'TAG'):
                                                Stop_first_pos_reverse.append(j+len(sequence)-2+1)
                                            if j+len(sequence)-2+1 >=  Start_first_pos[0]:
                                                break
                                                                            #J : g fanxiang
                                                                     
                                        for k in range(0,len(sequence),3):
                                            #if k > 10:
                                            if sequence[k:k+3] == 'ATG':
                                                    #coupling=sequence[0:k+13]
                                                    #nTATCTTGAAAGCTTTTGCTAAGGGAGAGGGAGGGAAGTATGCGGAAAGGAGAGGAGTTGTTCACCGGCGTAGTAoncoupling=sequence[k+13:j+len(sequence)-2+35+1]
                                                    #break
                                                    Start_first_pos_reverse.append(k+1)                                        
                                            if k+1 >= Stop_first_pos[0]:
                                                break
                                        #print(Start_first_pos_reverse)
                                        #print(Stop_first_pos_reverse)
                                            
                                        if not(len(Start_first_pos_reverse) == 0 and len(Stop_first_pos_reverse) ==0):
                                            score += abs(len(Start_first_pos_reverse)) + abs(len(Stop_first_pos_reverse))*100**70
                                        else:
                                            #print('h')
                                            if not(len(sequence[start_pos-1:end_pos]) == 9):
                                                score += abs(len(sequence[start_pos-1:end_pos])-9)*100**60
                                            else:
                                                #print('i')
                                                ASD='ACCTCCTTA'
                                                SD='TAAGGAGGT'
                                                stand = mrna_rna_energy(SD,ASD)
                                                print(sequence[start_pos-1:end_pos])
                                                tmp = mrna_rna_energy(sequence[start_pos-1:end_pos],ASD)
                                                if not(abs(stand-tmp) <=1):
                                                    score += abs(abs(stand-tmp)-1)*100**50
                                                    
                                                
                                                #if not(sequence[start_pos-1,end_pos] in SD_library):
                                                #    return float('inf')
                                                else:
    
                                                    print('j')                                        
                                                    RNAStructure_documentgeneration(sequence)
                                                    RNAStructure_mfe(sequence)
                                                    b = TUN_M_SD_structure(sequence,start_pos,end_pos,Start_first_pos[0])
                                                    
                                                    if b == 1:
                                                        score += 100**40
                                                    else:
                                                        print('k')
                                                        a = dGm(sequence,start_pos,end_pos)
                                                        ori_150bp = 'ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGTACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACT'
                                                        sfgfp_mins_36bp = 'CCTATTCTGGTGGAACTGGATGGTGATGTCAACGGTCATAAGTTTTCCGTGCGTGGCGAGGGTGAAGGTGACGCAACTAATGGTAAACTGACGCTGAAGTTCATCTGTACTACTGGTAAACTGCCGGTACCTTGGCCGACTCTGGTAACGACGCTGACTTATGGTGTTCAGTGCTTTGCTCGTTATCCGGACCATATGAAGCAGCATGACTTCTTCAAGTCCGCCATGCCGGAAGGCTATGTGCAGGAACGCACGATTTCCTTTAAGGATGACGGCACGTACAAAACGCGTGCGGAAGTGAAATTTGAAGGCGATACCCTGGTAAACCGCATTGAGCTGAAAGGCATTGACTTTAAAGAAGACGGCAATATCCTGGGCCATAAGCTGGAATACAATTTTAACAGCCACAATGTTTACATCACCGCCGATAAACAAAAAAATGGCATTAAAGCGAATTTTAAAATTCGCCACAACGTGGAGGATGGCAGCGTGCAGCTGGCTGATCACTACCAGCAAAACACTCCAATCGGTGATGGTCCTGTTCTGCTGCCAGACAATCACTATCTGAGCACGCAAAGCGTTCTGTCTAAAGATCCGAACGAGAAACGCGATCATATGGTTCTGCTGGAGTTCGTAACCGCAGCGGGCATCACGCATGGTATGGATGAACTGTACAAATGA'
                                                        a = a*RCBS(ori_150bp[3::]+sequence[:Stop_first_pos[0]-1])
                                                        a = a*RCBS(sequence[Start_first_pos[0]+3-1::]+sfgfp_mins_36bp[:-3])
                                                        print(a)
                                                        #score += (a-G_goal)**2
                                                        score = 1/a
    
    return score




#Tuner Design            

#from time import time
from sopt.util.functions import *
from sopt.util.pso_config import *
from sopt.PSO.PSO import PSO
from sopt.util.constraints import *
from sopt.GA.GA import GA
from sopt.util.functions import *
from sopt.util.ga_config import *
from sopt.util.constraints import *

def Tuner_design():
    global m3_library
    m3_library = m3_library_generation('','')
    global a
    global b
    a=[]
    #b=[]
    low=[-25,1,1,0]
    high=[25,par_space()-4-33,par_space()-4-33,len(m3_library)]        

    for i in range(par_space()-4-33):
        low.append(0)
        high.append(1)

    class TestGA:
        def __init__(self):
            self.func = Tuner_goal
            self.func_type = 'min'
            self.variables_num = len(high)
            self.lower_bound = np.array(low)
            self.upper_bound = np.array(high)
            self.cross_rate = 0.8
            self.mutation_rate = 0.05
            self.generations = 50
            self.population_size = 1000
            self.binary_code_length = 10
            self.cross_rate_exp = 1.0005
            self.mutation_rate_exp = 1.0005
            self.code_type = code_type.binary
            self.cross_code = False
            self.select_method = select_method.keep_best
            self.rank_select_probs = None
            self.tournament_num = 2
            self.cross_method = cross_method.uniform
            self.arithmetic_cross_alpha = 0.1
            self.arithmetic_cross_exp = 1
            self.mutation_method = mutation_method.uniform
            self.none_uniform_mutation_rate = 1
            #self.complex_constraints = [constraints1,constraints2,constraints3]
            self.complex_constraints = None
            self.complex_constraints_method = complex_constraints_method.penalty
            self.complex_constraints_C = 1e6
            self.M = 1e8
            self.GA = GA(**self.__dict__)
    
        def test(self):
            start_time = time.time()
            self.GA.run()
            print("GA costs %.4f seconds!" % (time.time()-start_time))
            #self.GA.save_plot()
            self.GA.show_result()
            a.append(self.GA.__dict__)
            #b.append(self.GA.global_best_point)
    
    
    
    if __name__ == '__main__':
        TestGA().test()



    print(a[0]['global_best_point'].tolist())
    #print(b)


'''                    
def main():
    
    m3_library = m3_library_generation('','')
    
    #SD_library = sd_generation()
    Tuners=[]
    #for i in range(1,6):
    #    for j in range(5):
    low=[-25,1,1]
    high=[25,par_space()-3,par_space()-3]        

    for i in range(par_space()-3):
        low.append(0)
        high.append(1)



    Tuners.append(PSO(Tuner_goal,'min',1000,1000,par_space(),high,low))
    
    return Tuners
'''





def asRNA_NUPACK_documentgeneration(i,type):
    
    if type ==1:        
        f = open('TBR_Pair_' + i + '.in','w')
        f.write('2')
        f.write('\n')
        f.write(i)
        f.write('\n')
        f.write(TIR)
        f.write('\n')
        f.write('1 2')
        f.close()
        
    elif type ==2:
        f = open('TBR_Only_' + i + '.in','w')
        f.write(i)
        f.close()
        
    elif type ==3:
        f = open('asRNA_' + i + '.in','w')
        f.write(i)
        f.close() 
        
def asRNA_NUPACK_mfe(i,type):
    #global code
    if type == 1:        
        code= "mfe -multi -material rna -degenerate " + 'TBR_Pair_' + i
        os.system(code)
        
    elif type == 2:
        code= "mfe -material rna -degenerate " + 'TBR_Only_' + i 
        os.system(code)
        
    elif type == 3:
        code= "mfe -material rna -degenerate " + 'asRNA_' + i
        os.system(code)
    
    #print(code)

def asRNA_NUPACK_G(i,type):
    G=0
    
    if type ==1:
        path = 'TBR_Pair_' + i + '.mfe'
        f = open(path,'r')
        counter = 0
        for i in f:
            counter = counter + 1
            if counter == 15:
                G = i
        f.close()
        
    elif type ==2:
        path = 'TBR_Only_' + i + '.mfe'
        f = open(path,'r')
        counter = 0
        for i in f:
            counter = counter + 1
            if counter == 15:
                G = i
        f.close()
        
    elif type ==3:
        path = 'asRNA_' + i + '.mfe'
        f = open(path,'r')
        counter = 0
        for i in f:
            counter = counter + 1
            if counter == 15:
                G = i
        f.close()
    
    return G




def asRNA_NUPACK_mfe_todot_toct(i):   # file name
    
    path = 'asRNA_' + i + '.mfe'
    f = open(path,'r')
    counter = 0
    for i in f:
        counter = counter + 1
        if counter == 16:
            dot = i
        elif counter == 4:
            sequence = i.split(':')[-1]
    f.close()
    
    path = 'asRNA_' + i + '.fold'
    f = open(path,'w')
    f.write('>')
    f.write(sequence)
    #f.write('\n')
    f.write(sequence)
    f.write(dot)
    f.close()
    
    path = 'asRNA_' + i + '.ct'
    f = open(path,'w')
    f.close()
    path = 'asRNA_' + i 
    Vienna_RNAStructure_dot2ct(path)





def asRNA_hfq_Structure(k,i,j,sty):   #k:sequence;i:sequence;j:ori_hfq_matrix 
    
    f = open('asRNA' + k +'.ct','r')
    counter = 0
    
    #global structure
    
    if sty == 1:

        structure = np.zeros((len(i),len(i)))        
        
        for i in f:
            if counter != 0 :
                    i = i.split('\n')[0].split()
                    link_pos = int(i[-2])
                    if link_pos !=0:                        
                        if link_pos != 0:                         
                            structure[counter-1,link_pos-1] = 1
                            structure[link_pos-1,counter-1] = 1                
            counter = counter + 1

    elif sty == 2:
        
        ori = j.shape[0]
        sequence = i
        structure = np.zeros((ori,ori))
        
        for i in f:
            if counter > len(sequence) - ori:
                i = i.split('\n')[0].split()
                counter_mat = counter - len(sequence)
                link_pos = int(i[-2])-len(sequence)

                if link_pos != 0:                    #judge_link is exit
                     
                    if link_pos > ori:               #外接
                        structure[:,counter_mat-1] = 10**10
                        structure[counter_mat-1,:] = 10**10
                    
                    else:                            #内接
                        
                        if 1 in j[counter_mat-1,:]:      #org_link is exit
                            
                            link_pos_ori = j[counter_mat-1,:].tolist().index(1)
                            
                            if link_pos == link_pos_ori:            #内接正确
                                structure[counter_mat-1,link_pos-1] = 1
                                structure[link_pos-1,counter_mat-1] = 1
                            
                            elif link_pos != link_pos_ori:         #内接错误之原来有连接
                                structure[counter_mat-1,link_pos-1] = 10**10
                                structure[link_pos-1,counter_mat-1] = 10**10
                        
                        else:                                      #内接错误之原来无连接
                            structure[counter_mat-1,link_pos-1] = 10**10
                            structure[link_pos-1,counter_mat-1] = 10**10
                    
                else:                                        # judge_link is not exit
                    if 1 in j[counter_mat-1,:]:
                        link_pos_ori = j[counter_mat-1,:].tolist().index(1)
                        structure[counter_mat-1,link_pos_ori] = 10**10
                        structure[link_pos_ori,counter_mat-1] = 10**10
            
            
            counter = counter + 1            
            
            #if counter > j.shape[0]:
            #    break
            
    return structure
 
    
    
def asRNA_deactivate_goal(x):
    #print(x)
    TIR = 'TTTATAAAGAGAAGACTATGAATGCCATGATCATGAGTAAA'
    TBR = 'TTTACTCATGATCATGGCATTCATAGTCTTCTCTTTATAAA'
    hfq='CGTCCCGCAAGGATGCGGGTCTGTTTACCCCTATTTCAACCGGCCGCCTCGCGGCCGGTTTTTTTTT'
    aaa = x
    counter1 = 0
    
    for i in range(len(TIR)):
        if 0 <= aaa[i] < 0.5:
            aaa[i] = 0
        elif 0.5<= aaa[i] <= 1:
            aaa[i] = 1
            counter1 = counter1 + 1
    
    print(counter1/len(TIR))
    
    if counter1/len(TIR) < 0.85:
        #print('a')
        return abs((counter1/len(TIR) - 0.85))*100*100**100
    else:
        print('a')
        pra1 = '['
        for i in range(len(TIR)):
           pra1 = pra1 + ('aaa[' + str(i) +'],' ) 
        pra1 = pra1[:-1]
        pra1 = pra1 + ']'        
        k1 = eval(pra1)

        print(k1)
        #data1 = FFFF[0,:,ii].tolist()
        
        new1 = ''

        for i in range(len(k1)):
            new1 = new1 + TBR[i]*int(k1[i])
        
        
        print(new1)
        
        
        asRNA_NUPACK_documentgeneration(new1,1)
        asRNA_NUPACK_mfe(new1,1)
        a = asRNA_NUPACK_G(new1,1)
        #documentmove(i)
        asRNA_NUPACK_documentgeneration(new1,2)
        asRNA_NUPACK_mfe(new1,2)
        b = asRNA_NUPACK_G(new1,2)
        #documentmove(i)
        if not((float(a)-float(b)) < (float(-40))):
            #print(i+':work')        
            #k_o = structure('')        
            #k = structure(i)
            #print(((np.array(k)-np.array(k_o))**2).sum())
            #documentmove(i)
            return ((float(a)-float(b)) - (float(-40)))*100*50
        
        else:
            print(0.3848 - 0.0068*(float(a)-float(b)) - 0.0125*(1-(counter1/len(TIR))) + 0.123)
            return 1/(0.3848 - 0.0068*(float(a)-float(b)) - 0.0125*(1-(counter1/len(TIR))) + 0.123)
        


def asRNA_deactivate_design():
    global TIR
    global TBR
    TIR = 'TTTATAAAGAGAAGACTATGAATGCCATGATCATGAGTAAA'
    TBR = 'TTTACTCATGATCATGGCATTCATAGTCTTCTCTTTATAAA'
    #def main_Tuner_design():
    #global m3_library
    #m3_library = m3_library_generation('','')
    global a
    a = []
    low = []
    high = []
    
    #low=[-25,1,1,0]
    #high=[25,par_space()-4-33,par_space()-4-33,len(m3_library)]        

    for i in range(len(TBR)):
        low.append(0)
        high.append(1)
        
        
    #for i in range(10):
        

    class TestGA:
        def __init__(self):
            self.func = asRNA_deactivate_goal
            self.func_type = 'min'
            self.variables_num = len(high)
            self.lower_bound = np.array(low)
            self.upper_bound = np.array(high)
            self.cross_rate = 0.8
            self.mutation_rate = 0.05
            self.generations = 500
            self.population_size = 100
            self.binary_code_length = 20
            self.cross_rate_exp = 1
            self.mutation_rate_exp = 1
            self.code_type = code_type.binary
            self.cross_code = False
            self.select_method = select_method.keep_best
            self.rank_select_probs = None
            self.tournament_num = 2
            self.cross_method = cross_method.uniform
            self.arithmetic_cross_alpha = 0.1
            self.arithmetic_cross_exp = 1
            self.mutation_method = mutation_method.uniform
            self.none_uniform_mutation_rate = 1
            #self.complex_constraints = [constraints1,constraints2,constraints3]
            self.complex_constraints = None
            self.complex_constraints_method = complex_constraints_method.penalty
            self.complex_constraints_C = 1e6
            self.M = 1e8
            self.GA = GA(**self.__dict__)
    
        def test(self):
            start_time = time.time()
            self.GA.run()
            print("GA costs %.4f seconds!" % (time.time()-start_time))
            #self.GA.save_plot()
            self.GA.show_result()
            a.append(self.GA.__dict__)
    
    
    
    if __name__ == '__main__':
        TestGA().test()

def asRNA_activate_goal(x):
    #print(x)
    TIR = 'AATGATCTCCTTTTTAAGTGTAAGGGCCCAAGTTCACTTAAAAAGGAGATCAACTAATG'
    #TBR = 'TTTACTCATGATCATGGCATTCATAGTCTTCTCTTTATAAA'
    hfq='CGTCCCGCAAGGATGCGGGTCTGTTTACCCCTATTTCAACCGGCCGCCTCGCGGCCGGTTTTTTTTT'
    
    sequence = dis(sequence_mat) + hfq
    
    print(sequence)
    


        
    asRNA_NUPACK_documentgeneration(sequence,1)
    asRNA_NUPACK_mfe(sequence,1)
    a = asRNA_NUPACK_G(sequence,1)
    #documentmove(i)
    
    
    
    asRNA_NUPACK_documentgeneration(new1,2)
    asRNA_NUPACK_mfe(new1,2)
    c = asRNA_NUPACK_G(new1,2)
        

        


def asRNA_activate_design():
    global TIR
    global TBR
    TIR = 'TTTATAAAGAGAAGACTATGAATGCCATGATCATGAGTAAA'
    #TBR = 'TTTACTCATGATCATGGCATTCATAGTCTTCTCTTTATAAA'
    #def main_Tuner_design():
    #global m3_library
    #m3_library = m3_library_generation('','')
    global a
    a = []
    low = []
    high = []
    
    #low=[-25,1,1,0]
    #high=[25,par_space()-4-33,par_space()-4-33,len(m3_library)]        

    for i in range(34):
        low.append(0)
        high.append(1)
        
        
    #for i in range(10):
        

    class TestGA:
        def __init__(self):
            self.func = asRNA_activate_goal
            self.func_type = 'min'
            self.variables_num = len(high)
            self.lower_bound = np.array(low)
            self.upper_bound = np.array(high)
            self.cross_rate = 0.8
            self.mutation_rate = 0.05
            self.generations = 500
            self.population_size = 100
            self.binary_code_length = 20
            self.cross_rate_exp = 1
            self.mutation_rate_exp = 1
            self.code_type = code_type.binary
            self.cross_code = False
            self.select_method = select_method.keep_best
            self.rank_select_probs = None
            self.tournament_num = 2
            self.cross_method = cross_method.uniform
            self.arithmetic_cross_alpha = 0.1
            self.arithmetic_cross_exp = 1
            self.mutation_method = mutation_method.uniform
            self.none_uniform_mutation_rate = 1
            #self.complex_constraints = [constraints1,constraints2,constraints3]
            self.complex_constraints = None
            self.complex_constraints_method = complex_constraints_method.penalty
            self.complex_constraints_C = 1e6
            self.M = 1e8
            self.GA = GA(**self.__dict__)
    
        def test(self):
            start_time = time.time()
            self.GA.run()
            print("GA costs %.4f seconds!" % (time.time()-start_time))
            #self.GA.save_plot()
            self.GA.show_result()
            a.append(self.GA.__dict__)
    
    
    
    if __name__ == '__main__':
        TestGA().test()

