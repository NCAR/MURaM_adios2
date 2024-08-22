import numpy as np
import re

#eos_var=['eosT','eosP', 'eosne','eosrhoi','eosamb','Qtot','tau','Jtot','Stot','QxCor']
eos_var=['eosT', 'eosP', 'Qtot',  'tau','QxCor']
diag_var=['1', '2', '3', '4', '5',  '6',  '7',  '8',  'Qres',  'Qvis', 'Qamb']
tau_lev=['1.0','0.1','0.01']
def parsed_param_dict(filename):
    if not filename:
        filename = open('parameters.dat', 'r')
    Lines = filename.readlines()


    parameters={}
    count=0
    for line in Lines:
        key1=""
        #print("Line{}: {}".format(count, line.partition('=')))
        res=[]
        if not line.startswith("#"):
            count=0
            sindex=line.find('|')
            line=line[0:sindex]
            for i in line.split(' '):
                i=i.strip('=\n\t; ')
                if i.isnumeric():
                    #print('num i=',i)
                    res.append(int(i))
                elif len(i)!=0:
                    if i.find('.')!=-1 and (i[0].isnumeric() or i.find('e')!=-1):
                       #print('float i=',i)
                       res.append(i.strip('\n\t;'))
                    elif i.find('.')!=-1 and i[0]=="-" and i[1].isnumeric() :
                       res.append(i.strip('\n\t;'))
                    #elif i[0].isnumeric() and len(i)==2:
                       #print('int i=',i,len(i))
                    #   res.append(int(i))
                    else:
                       if count !=0 :
                           res.append(i)
                       #print('string=',i)
                if count==0:
                    if i!="" and i.find('\t')==-1 and i.find('\n')==-1:
                        #print('count 0 ',i)
                        key1=i
                count += 1
      
        if len(key1)!=0: 
            parameters[key1]=res 

    parameters['eos_var']=eos_var
    parameters['diag_var']=diag_var
    #for k,i in parameters.items():
    #  print(k,i)
    return parameters


def get_output_varname(eos_var,eos_output):
    eos_output_list=[]
    for c,i in enumerate(eos_output):
        if i:
            eos_output_list.append(eos_var[c])
    #print(eos_output_list)
    return eos_output_list
 

def format_level(lev):
    lev_new = []
    #lev = tau_lev
    for i in lev:
        if float(i) >= 0.001:
            tmp_str =  "%.3f" % float(i)
        else:
            tmp_str =  "%.6f" % float(i)
        lev_new.append(tmp_str)
    return lev_new

def parsed_backup_dict(filename):
    if not filename:
        filename = open('backup.dat', 'r')
    Lines = filename.readlines()

    keys=[]
    count=0
    for line in Lines:
        for i in line.split(' '):
            keys.append(i)
    return keys[0] 

par=parsed_param_dict("")
backup=parsed_backup_dict("")
#print(par['tau_lev'])
#get_output_varname(par['eos_var'],par['eos_output'])
#get_output_varname(par['diag_var'],par['diag_output'])
format_yz_lev=[]
yz_lev=par['yz_lev']
for i in yz_lev:
    format_yz_lev.append( "{:04d}".format(i))
print(format_yz_lev)
