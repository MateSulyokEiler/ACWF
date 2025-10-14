#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 11:10:59 2022

@author: raxis
"""

import os
import math
import numpy as np
import json
import re
from .basic import *

#to import read_save, write_save, rama, rama_types, amy_struct_list, amy_struct, amy_interface, amy_chain, atlag, atlag_class, distribution
# amy_listafull, amy_lista, amy_lista_prot, amy_lista_protunique


def read_save(name):
    
    if name.split('.')[-1] == 'json':
        output = read_save_json(name)
    else:
        output = read_save_txt(name)

    return output


def write_save(name, amy_lista):
    
    if name.split('.')[-1] == 'json':
        write_save_json(name, amy_lista)
    else:
        write_save_txt(name, amy_lista)

    return 


def read_save_txt(name):
    
    with open(name, 'r+') as f:
        amy_lista = amy_struct_list()
        
        savefileraw = f.read().splitlines()
        savefile = []
        for line in savefileraw:
            if line.split()[0]=='-':
                savefile[-1]+=line[1:]
            else:
                savefile.append(line)


        for line in savefile:

            data = line.split()
            
            if data[0] == '#':
                continue

            if data[0] == 'i':
                interface = amy_interface()

                n = 1
                while n < len(data):
                    setattr(interface, data[n], data[n+1])
                    n +=2
                
                amy_lista.lista[-1].add_interface(interface)

            elif data[0] == 'c':
                chain = amy_chain()
                n = 1
                while n < len(data):
                    setattr(chain, data[n], data[n+1])
                    n +=2
                    
                amy_lista.lista[-1].add_chain(chain)

            elif data[0] == 'phi':
                if len(data) == 1:
                    continue

                phi = data[1:]
                amy_lista.lista[-1].chain[-1].phi=phi

            elif data[0] == 'psi':
                if len(data) == 1:
                    continue

                psi = data[1:]
                amy_lista.lista[-1].chain[-1].psi=psi

            elif data[0] == 'e':

                er = data[1:]
                amy_lista.lista[-1].error = er

            elif data[0] == 'A':
                amy = amy_struct()
                n = 1
                while n < len(data):

                    setattr(amy, data[n], data[n+1])
                    n +=2
                
                amy_lista.add(amy)

    return amy_lista


def write_save_txt(name, amy_lista):

    with open(name, 'w+') as f:

        for amy in amy_lista.lista:

            line = 'A '
            x = 0
            for attr in vars(amy).keys():
                
                if attr in ['interface', 'chain', 'error']:
                    continue
                if x%6 == 5:
                    line += '\n- '
                x += 1
                line += attr + ' ' + str(getattr(amy, attr))+' '
            
            line += '\n'
            f.write(line)
            
            if hasattr(amy, 'error'):
                line = 'e'
                if amy.error == []:
                    line += ' -'
                else:   
                    for er in amy.error:
                        line += ' ' + er
        
                line += '\n'
                f.write(line)

            for interface in amy.interface:
                line='i '
                x = 0
                for attr in vars(interface).keys():
                    if x%6 == 5:
                        line += '\n- '
                    x += 1
                    line += attr + ' ' + str(getattr(interface, attr))+' '
                
                line += '\n'
                f.write(line)

            for chain in amy.chain:
                line = 'c '
                x = 0
                for attr in vars(chain).keys():
                    if attr in ['phi', 'psi']:
                        continue
                    if x%6 == 5:
                        line += '\n- '
                    x +=1
                    line += attr + ' ' + str(getattr(chain, attr))+' '
                
                line += '\n'
                f.write(line)

                line = 'phi '
                for a in chain.phi:
                    line += str(a)+' '
                f.write(line+'\n')

                line = 'psi '
                for a in chain.psi:
                    line += str(a)+' '
                f.write(line+'\n')

    return


def read_save_json(name):
    
    with open(name, "r") as f: 
        json_data = json.load(f)
    
    amy_lista = amy_struct_list()
    
    for amy_data in json_data['lista']:
        amy = amy_struct()
        
        for var in amy_data:
            
            if var == 'interface':
                for i_data in amy_data[var]:
                    interface = amy_interface()
                    for i_var in i_data:
                        setattr(interface, i_var, i_data[i_var])
                    amy.add_interface(interface)
                
            elif var == 'chain':
                for c_data in amy_data[var]:
                    chain = amy_chain()
                    for c_var in c_data:
                        
                        if c_var in ['phi', 'psi', 'seq', 'resi', 'ramach']:
                            setattr(chain, c_var, c_data[c_var][:])
                        else:
                            setattr(chain, c_var, c_data[c_var])
                    amy.add_chain(chain)
                
            elif var == 'error':
                amy.error = amy_data[var][:]

            else:
                setattr(amy, var, amy_data[var])
                
        amy_lista.add(amy)

    return amy_lista


def write_save_json(name, amy_lista):

    with open(name, 'w+') as f:
        json_txt = json.dumps(amy_lista, default=lambda obj: obj.__dict__, indent=4)
        #json_txt = json_txt.replace('\n                        ','')
        #json_txt = json_txt.replace('\n                    ]',']')
        json_txt = re.sub(r'\n                        ', '', json_txt)
        json_txt = re.sub(r'\n                    ]', ']', json_txt)
        f.write(json_txt)

    return




def rama(phi, psi, rama_type='gen'):
    #Molprobity/pyrama
    #ramachandran general aminoacids

    
    rama_map=rama_types[rama_type]    
    d = rama_map[math.floor(phi)][math.floor(psi)]
    
    #B = beta favoured
    #A = alpha favoured
    #L = left-handed
    #a = allowed
    #o = outlier
        
    if d > 0.02:
        if phi > 0 and phi < 120:
            res = 'L'
        else:
            if psi<49 and psi>-120:
                res ='A'
            else:
                res = 'B'    
    elif d > 0.0005:
        if phi > 0 and phi < 120:
            res = 'l'
        else:
            if psi<49 and psi>-120:
                res ='a'
            else:
                res = 'b' 
    else:
        res = 'o'
        
    return res




class amy_struct_list:

    def __init__(self):
        self.lista = []


    def add(self, struct):
        self.lista.append(struct)
        
        
    def add_blank(self):
        self.lista.append(amy_struct())
        self[-1].code = '-'
        self[-1].source = 'personal'

    def select(self, attr, select):
        #creates a selection 
        #accepts a list for select

        output = amy_struct_list()

        for amy in self.lista:
            if getattr(amy, attr) in select:
                output.add(amy)

        return output


    def change(self, attr, select, new):
        #changes a value if select is true(could do multiple)

        for amy in self.lista:
            if getattr(amy, attr) in select:
                setattr(amy, attr, new)

        return


    def data(self, attr):
        #expose a single attr for all amyloid structure

        output = []

        for amy in self.lista:
            output.append(getattr(amy, attr))

        return output


    def data_i(self, attr):
        #expose a single attr for all amyloid structure interface

        output = []

        for amy in self.lista:
            for i in amy.interface:
                try:
                    output.append(float(getattr(i, attr)))
                except Exception:
                    print(amy.code + ' error during data_i')
        #print(len(output))
        return output


    def sort(self, attr):
        #sorts itself according to attr

        self.lista = sorted(self.lista, key=lambda atom: getattr(atom, attr))


    def prnt(self):

        for n, amy in enumerate(self.lista):
            print(n+1)
            amy.prnt()
            print()


    def __getitem__(self, key):
        
        if isinstance(key, str):
            for amy in self.lista:
                if amy.code==key:
                    return amy
            raise KeyError(key)
        
        elif isinstance(key, int):
            return self.lista[key]
        
        raise TypeError(key)


    def l(self, number=6):

        selection=self.select('seq_len', str(number))

        return selection


    def c(self, classes='12345678'):

        selection=self.select('clas', classes)

        return selection



#amyloid structure
class amy_struct:

    def __init__(self):

        self.interface = []
        self.chain = []

    def  __getitem__(self, key):
        try:
            return getattr(self, key)
        except AttributeError:
            return '-'
    
    
    def add_interface(self, interface):

        self.interface.append(interface)


    def add_chain(self, chain):

        self.chain.append(chain)


    def data(self, attr, num=False):
        #expose a single attr for all interface

        output = []

        for interf in self.interface:
            try:
                if num:
                    output.append(float(getattr(interf, attr)))
                else:
                    output.append(getattr(interf, attr))
            except:
                if num:
                    output.append(0)
                else:
                    output.append('0')
        return output


    def select(self, attr, select):
        #creates a selection within a model
        #accepts a list for select

        output = amy_struct()

        for amy in self.lista:
            if getattr(amy, attr) in select:
                output.add(amy)

        return output


    def get_chain(self, attr, select):
        #returns the first matching chain
        
        output = next((x for x in self.chain if getattr(x, attr) == select), None)
        
        return output


    def prnt(self):
        line=''

        line += self.code+' '
        line += self.seq+' '
        line += self.clas+' '
        line += self.lat_type+' '
        line += self.seq_len+' '
        
        print(line)



class amy_interface:

    def __init__(self):
        pass


    def  __getitem__(self, key):
        try:
            return getattr(self, key)
        except AttributeError:
            return '-'


    def prnt(self):

        line='i '
        line += self.sc+' '
        line += self.area+' '
        line += self.areaimol

        print(line)



class amy_chain:

    def __init__(self):
        #necessary

        self.phi = []
        self.psi = []


    def  __getitem__(self, key):
        try:
            return getattr(self, key)
        except AttributeError:
            return '-'


    def prnt(self):

        line='chain '
        line += self.letter+' '
        line += self.phi+' '
        line += self.psi+' '

        print(line)



def atlag(attr, amy_lista):

    data = []
    for a in amy_lista.data(attr):

        try:
            data.append(float(a))
        except:
            pass

    atlag =  sum(data)/len(data)
    szor = statistics.stdev(data)

    return atlag, szor



def atlag_class(attr, amy_lista):

    for clas in ['1', '2', '4', '5', '6', '7', '8', 'L', 'O']:

        data = []
        for a in amy_lista.select('clas', clas).data(attr):
            try:
                data.append(float(a))
            except:
                pass

        try:
            atlag =  sum(data)/len(data)
        except:
            print(clas, '-')

        try:
            szor = statistics.stdev(data)
            print(clas, round(atlag, 4), round(szor, 4))
        except:
            szor = '-'
            print(clas, round(atlag, 4), szor)

    return



def distribution(attr, amy_lista, n):

    data = []
    for a in amy_lista.data(attr):
        try:
            data.append(float(a))
        except:
            pass

    minim = min(data)
    step = (max(data)-min(data))/n

    dist = {}
    for a in range(n+1):
        dist[a]=0


    for a in data:
        dist[int((a-minim)//step)]+=1

    for a in dist.keys():
        print(round(minim+a*step, 4), dist[a])

    return


