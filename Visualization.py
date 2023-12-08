#from Hybrid_Optimizer_OPTaaS import *   
from Core_Definition import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
import numpy as np

#fig=plt.figure(1,figsize=(15,10))
#fig.suptitle("test for multi bar chart")
#a=plt.subplot(1,2,1)
#b=plt.subplot(1,2,2)


#value1=tempz.population
#xindex1=range(len(tempz.population))
#a.bar(xindex1,value1,alpha=0.5)
#a.set_title("zz correlation")
#a.set_xlabel("|i-j|")


#value2=tempx.population
#xindex2=range(len(tempx.population))
#b.bar(xindex1,value1,alpha=0.5)
#b.bar(xindex1,value2,alpha=0.5)
#b.set_title("zz correlation")
#b.set_xlabel("|i-j|")

plt.show()

def Display_Optimization(log,**kwargs):
    params_history=dict()
    functions_history=dict()
    iterations=0
    read_file=open(log,"r")
    while True:
        temp=read_file.readline().split(":")
        if len(temp[0])==0:
            break
        if len(temp)==0:
            break
        current_index=int(temp[0])
        if temp[1]=="parameter":
            if temp[2] in params_history:
                if len(params_history[temp[2]])>current_index:
                    params_history[temp[2]][current_index]=float(temp[3])
                else:
                    params_history[temp[2]].append(float(temp[3]))
            else:
                params_history[temp[2]]=[float(temp[3])]
                
        elif temp[1]=="function":
            if temp[2] in functions_history:
                if len(functions_history[temp[2]])>current_index:
                    functions_history[temp[2]][current_index]=float(temp[3])
                else:
                    functions_history[temp[2]].append(float(temp[3]))
            else:
                functions_history[temp[2]]=[float(temp[3])]
                    
    for key in functions_history:
        length=len(functions_history[key])

    index=range(length)
    
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20) 
          
    fig=plt.figure(1,figsize=(20,12))
    fig.suptitle("optimization visualization")
    a=plt.subplot(2,1,1)
    b=plt.subplot(2,1,2)
    
    print(functions_history)
    
    print(params_history)

    a.set_ylabel("cost function",fontsize=20)
    a.set_xlabel("iterations",fontsize=20)
    a.xaxis.set_major_locator(MaxNLocator(integer=True))
    for i in functions_history:
        if i=="Energy":
            a.plot(index,functions_history[i][:length],linestyle='-',linewidth=5,markersize=10, marker='d',color="k",label=i)  
        else:
            a.plot(index,functions_history[i][:length],linestyle='-',linewidth=5,markersize=10, marker='o',label=i)                                   
    a.legend(fontsize=20)
    
    b.set_ylabel("paramters",fontsize=20)
    b.set_xlabel("iterations",fontsize=20)
    b.xaxis.set_major_locator(MaxNLocator(integer=True))
    for i in params_history:
        b.plot(index,params_history[i][:length],linestyle='-',linewidth=5,markersize=10, marker='o',label=i)
    b.legend(fontsize=20)         
    plt.show()
      
    if "savefile" in kwargs:
        fig.savefig(kwargs["savefile"]+".svg",format='svg')
        
def Display_States_Population(**kwargs):
    plt.rc('xtick', labelsize=20) 
    plt.rc('ytick', labelsize=20) 
    if "plot_size" in kwargs:
        fig=plt.figure(1,figsize=(kwargs["plot_size"][0],kwargs["plot_size"][1]))
    else:
        fig=plt.figure(1,figsize=(10,10))
        
    fig.suptitle("")
    a=plt.subplot(1,1,1)
    if "title" in kwargs:
        a.set_title(kwargs["title"],fontsize=20)
    a.set_xlabel("states",fontsize=20)
    a.set_ylabel("probability",fontsize=20)
    #a.set_ylim([0,0.25])
    tot=len(kwargs["states"])
    width=0.8/tot
    shift=0
    for i in range(tot):
        temp=kwargs["states"][i]
        value=temp.population
        xindex=np.arange(len(value))+shift
        a.bar(xindex,value,width,alpha=0.9,label=kwargs["label"][i])
        shift+=width
    a.legend(fontsize=20)
    a.set_xticks(np.arange(len(value))+(shift-width)/2)
    a.set_xticklabels(np.arange(len(value)))
    plt.show()
    
    if "savefile" in kwargs:
        fig.savefig(kwargs["savefile"]+".svg",format='svg')
    
def Display_States_Population_3D(x_title,y_title,**kwargs):
    fig = plt.figure(figsize=(15,15))
    ax = plt.subplot(1,1,1, projection='3d')
    
    x=[]
    y=[]
    z=[]
    top=[]
    x_counter=0
    y_counter=0
    for key in kwargs:
        for i in kwargs[key].population:
            x.append(int(x_counter))
            y.append(int(y_counter))
            z.append(0)
            top.append(i)
            x_counter+=1
        x_counter=0
        y_counter+=1
     
    dx=dy=0.5
    x=np.array(x)
    y=np.array(y)
    z=np.array(z)
    
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    #ax.w_xaxis.set_ticklabels(y_mark)
    #ax.w_yaxis.set_ticklabels(y_mark)
    ax.bar3d(y, x, z, dx, dy, top,  shade=True)
    ax.set_xlabel(x_title)
    ax.set_ylabel(y_title)
    plt.show()
    
    if "savefile" in kwargs:
        fig.savefig(kwargs["savefile"]+".svg",format='svg')
        
def Display_List(**kwargs):
    fig=plt.figure(1,figsize=(10,10))
    a=plt.subplot(1,1,1)
    
    for key in kwargs:
        if key.find("list")>-1:
            a.scatter(kwargs[key][0],kwargs[key][1],label=key[5:])
           
    a.legend(fontsize=20)
    
    if "title" in kwargs:
        a.set_title(kwargs["title"])
    if "xlabel" in kwargs:
        a.set_xlabel(kwargs["xlabel"])
    if "ylabel" in kwargs:
        a.set_ylabel(kwargs["ylabel"])

    plt.show()
    if "savefile" in kwargs:
        fig.savefig(kwargs["savefile"]+".svg",format="svg")