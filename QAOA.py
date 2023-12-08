# Daiwei Zhu Aug-9-2018
# an add-on to the hybrid optimizer framework for QAOA problem

__author__="Daiwei Zhu"

import numpy as np
from Hybrid_Optimizer import *
from Core_Definition import *

def Point_Correlation(state,i,j):
    temp_state=np.multiply(np.add(state,-0.5),-2)
    return temp_state[i]*temp_state[j]

def Point_Correlation_Matrix(state):
    temp=[[0]*state.size for i in range(state.size)]
    for i in range(state.size):
        for j in range(state.size):
            for k in range(len(state.population)):
                temp[i][j]+=state.population[k]*Point_Correlation(Binary_Decoding(k,state.size),i,j)
    return temp
            

class QAOA_Circuit(Hybrid_Quantum_Circuit):

    def Add_Layer(self,name,value=0):
        Define_Parameter(name,value)
        if name.find("a")>-1:
            for i in range(self.size-1):
                self.Add_Gate(Quantum_Gate("XA",i,i+1,angle=name))
            self.Add_Gate(Quantum_Gate("XA",0,self.size-1,angle=name))
        elif name.find("b")>-1:
            for i in range(self.size):
                self.Add_Gate(Quantum_Gate("AZ",i,angle=name))
   
    def Generate_Gradient_Probe(self):
        file_name=self.name+"_cache/circuit_"+str("{:03d}").format(self.iteration_number)+".txt"  
        write_file=open(file_name,'w')
        write_file.write(self.GatesLab_Sequence()+'\n')
        write_file.write(self.GatesLab_Sequence_XMeasurement()+'\n')
        for key in param_table:
            param_table[key]+=self.probe_distance
            Regulate_Params()
            write_file.write(self.GatesLab_Sequence()+'\n')
            write_file.write(self.GatesLab_Sequence_XMeasurement()+'\n')
            param_table[key]-=self.probe_distance
            Regulate_Params()
        write_file.close()
        return self.name+"_cache/circuit_"+str("{:03d}").format(self.iteration_number)+".txt"  
           
    
    
    def Process_Gradient_Result(self):
        prefix=self.name+"_readout/"
        file_list=sorted(os.listdir(prefix),key=lambda x: int(x[x.find("_Line_")+6:x.find(".txt")]))
        temp_state=Quantum_State(self.size)
        temp_state_x=Quantum_State(self.size)
        assert (int(file_list[0][file_list[0].find("_Line_")+6:file_list[0].find(".txt")])==1)
        temp_state.Import(prefix+file_list[0])
        temp_state_x.Import(prefix+file_list[1])
        energy_now=temp_state.Mean_Value(Energy_Z)+temp_state_x.Mean_Value(Energy_ZZ)
        self.Log_Functions("Energy",energy_now)
        print(energy_now)
        i=1;
        for key in param_table:
            temp_state.Import(prefix+file_list[i*2])
            temp_state_x.Import(prefix+file_list[i*2+1])
            energy_temp=temp_state.Mean_Value(Energy_Z)+temp_state_x.Mean_Value(Energy_ZZ)
            print(energy_temp)
            self.Log_Functions("Probe_"+key,energy_temp)
            param_table[key]-=(energy_temp-energy_now)*self.step_size/self.probe_distance
            i+=1
        self.iteration_number+=1
        Regulate_Params()
        self.Log_Parameters()


    
     
    