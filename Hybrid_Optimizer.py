#Daiwei Zhu Aug-9-2018
#A general framework for hybrid optimization using gradient descent (for now, at least)
__author__="Daiwei Zhu"


import numpy as np
import os
from Core_Definition import *


class Hybrid_Quantum_Circuit(Quantum_Circuit):
    def __init__(self,size,name="current"):
        self.size=size
        self.depth=0
        self.gates=[]
        self.probe_distance=0.01
        self.step_size=0.001
        self.iteration_number=0
        self.name=name
        if os.path.isdir(self.name+"_cache")==False:
            os.makedirs(self.name+"_cache")
        if os.path.isdir(self.name+"_readout")==False:
            os.makedirs(self.name+"_readout")
       
     
    def Log_Parameters(self):
        file_name=self.name+"_"+"log.txt"
        write_file=open(file_name,'a')
        for key in param_table:
            write_file.write(str(self.iteration_number)+':parameter:'+key+':'+str(param_table[key])+"\n")
        write_file.close()
        
    def Log_Functions(self,function,value):
        file_name=self.name+"_"+"log.txt"
        write_file=open(file_name,'a')
        write_file.write(str(self.iteration_number)+':function:'+function+':'+str(value)+"\n")
        write_file.close()
        
        
    def Load_Parameters(self,num_iterations):
        file_name=self.name+"_"+"log.txt"
        read_file=open(file_name,"r")
        while True:
            temp=read_file.readline().split(":")
            if len(temp)==0:
                break
            if len(temp[0])==0:
                break
            current_index=int(temp[0])
            if current_index==num_iterations:
                self.iteration_number=num_iterations
                if temp[1]=="parameter":
                    param_table[temp[2]]=float(temp[3])
                   
        read_file.close()
        assert (self.iteration_number==num_iterations)
        
         
        
    def Add_Gate(self,quantum_gate):
        self.depth+=1
        self.gates.append(quantum_gate)
        
    
    def Generate_Gradient_Probe(self):
        file_name=self.name+"_cache/circuit_"+str("{:03d}").format(self.iteration_number)+".txt"  
        write_file=open(file_name,'w')
        write_file.write(self.GatesLab_Sequence()+'\n')
        for key in param_table:
            param_table[key]+=self.probe_distance
            Regulate_Params()
            write_file.write(self.GatesLab_Sequence()+'\n')
            param_table[key]-=self.probe_distance
            Regulate_Params()
        write_file.close()
        return self.name+"_cache/circuit_"+str("{:03d}").format(self.iteration_number)+".txt"
           
    
    
    def Process_Gradient_Result(self,cost_function):
        prefix=self.name+"_readout/"
        file_list=sorted(os.listdir(prefix),key=lambda x: int(x[x.find("_Line_")+6:x.find(".txt")]))
        temp_state=Quantum_State(self.size)
        assert (int(file_list[0][file_list[0].find("_Line_")+6:file_list[0].find(".txt")])==1)
        temp_state.Import(prefix+file_list[0])
        energy_now=temp_state.Mean_Value(cost_function)
        i=1
        for key in param_table:
            temp_state.Import(prefix+file_list[i])
            energy_temp=temp_state.Mean_Value(cost_function)
            param_table[key]-=(energy_temp-energy_now)*self.stepsize/self.probe_distance
            i+=1
        self.iteration_number+=1
        Regulate_Params()
        Self.Log_Parameters()

