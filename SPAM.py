from Core_Definition import *
import numpy as np
import os

class SPAM_Operator:
    def __init__(self,size):
        self.size=size
        self.matrix=[0]*(2**size)
        for i in range(2**size):
            self.matrix[i]=[0]*(2**size)
        for i in range(2**size):
            self.matrix[i][i]=1
     
    def Import(self,file_name):
        read_file=open(file_name,'r')
        for i in range(2**self.size):
            temp=read_file.readline().split(',')
            for j in range(2**self.size):
                self.matrix[i][j]=float(temp[j])
        read_file.close()
        self.matrix=np.transpose(self.matrix)
        
    def Apply_SPAM(self,state):
        temp=Quantum_State(state.size)
        temp.population_only=True
        temp.population=np.matmul(self.matrix,state.population)
        return temp
    
    def Correct_SPAM(self,state):
        temp=Quantum_State(state.size)
        temp.population_only=True
        temp.population=np.linalg.lstsq(self.matrix,state.population)[0]
        for i in range(len(temp.population)):
            if temp.population[i]<0:
                temp.population[i]=0
        temp.population=np.divide(temp.population,np.sum(temp.population))
        temp.error_lower=np.matmul(self.matrix,state.error_lower)
        temp.error_upper=np.matmul(self.matrix,state.error_upper)
        return temp
   
        
