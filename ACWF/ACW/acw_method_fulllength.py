#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 11:10:59 2022

@author: raxis
"""

import traceback

import os
import math
import numpy as np
import platform
import subprocess
import shutil
from itertools import combinations
from copy import deepcopy
#import scipy


import pymol
from pymol import cmd
from pymol import stored
from pymol.wizard import Wizard
from pymol.cgo import *

from .basic import *
from .pdbdata import *
from .pdbcoord import *
from .acw_method_common import *



def straigthen_pca_check(name, name2, sheet_count, count=1):
    """Checks if the current selection is aligned to Y."""

    """
    #sheet 0 middle CA
    sheet = pdb_reader('Sheet_0.pdb').splice_same()

    seq_len = len(sheet.splice()[0].residues())

    atoms = []
    for atom in sheet.select('resi', str(int(seq_len/2.0))).select('name', ['CA']).lista:
        atoms.append(atom.xyz())
    """

    #both sheet centrums
    
    if len(pdb_reader('Sheet_0.pdb').splice(sort=False)) == sheet_count:
        sheet_0 = pdb_reader('Sheet_0.pdb').splice()[int(sheet_count/2.0)-1:int(sheet_count/2.0)+2]
    else:
        sheet_0 = pdb_reader('Sheet_0.pdb').splice_equal(sheet_count)[int(sheet_count/2.0)-1:int(sheet_count/2.0)+2]
    if len(pdb_reader('Sheet_%s.pdb'% str(count)).splice(sort=False)) == sheet_count:
        sheet_1 = pdb_reader('Sheet_%s.pdb'% str(count)).splice()[int(sheet_count/2.0)-1:int(sheet_count/2.0)+2]
    else:
        sheet_1 = pdb_reader('Sheet_%s.pdb'% str(count)).splice_equal(sheet_count)[int(sheet_count/2.0)-1:int(sheet_count/2.0)+2]
    
    atoms = []
    for n in range(len(sheet_0)):
        atoms.append(get_midpoint_coord(sheet_0[n].center(),sheet_1[n].center()))

    #calc for rotate axis to Y
    data , eigenvector_subset = PCA(atoms, 1)

    vector = [x[0] for x in eigenvector_subset]

    v = np.array(vector)
    v = v / np.sqrt(np.sum(v ** 2)) #norm
    angle = - math.degrees(math.acos(np.dot((0, 1, 0), v))) #angle from dot
    v_cross = list(np.cross((0, 1, 0), v)) #axis from cross

    #limit to skip
    if angle > -0.5 or angle < -179.5:
        return None, None

    cmd.rotate(v_cross, angle=angle, selection='all', camera=0, object=None, origin=[0, 0, 0])
    cmd.save('Sheet_0.pdb', name)
    cmd.save('Sheet_%s.pdb'% str(count), name2)
    cmd.rotate(v_cross, angle=-angle, selection='all', camera=0, object=None, origin=[0, 0, 0])
    
    return v_cross, angle


def save_coords(sele):
    """gets coordinates of named selection from pymol"""
    
    data = {'x_coord': {}, 'y_coord': {},'z_coord': {}}

    cmd.iterate_state(-1, sele, 'x_coord[(model,chain,resi,name, index)] = x', space=data)
    cmd.iterate_state(-1, sele, 'y_coord[(model,chain,resi,name, index)] = y', space=data)
    cmd.iterate_state(-1, sele, 'z_coord[(model,chain,resi,name, index)] = z', space=data)
    
    return data


def reset_coords(sele, data):
    """resets coordinates of named selection from pymol"""

    cmd.alter_state(-1, sele, 'x=x_coord[model,chain,resi,name, index]', space=data)
    cmd.alter_state(-1, sele, 'y=y_coord[model,chain,resi,name, index]', space=data)
    cmd.alter_state(-1, sele, 'z=z_coord[model,chain,resi,name, index]', space=data)

    return 


def color_methods(self, method):
    """Recolors the sheets with various methods."""

    cmd.select('sele', 'all')
    if method == 'even-odd__':

        cmd.color('atomic', 'all')
        cmd.color('green', 'elem C')
        for x in self.cur_pdb().chain[0].resi:
            if int(x) % 2 == 0:
                cmd.color('orange', 'sele and elem C and resi %s'%str(x))

    elif method == 'ramachandran__':
        
        cmd.util.cbag('all')
        #colorful
        types = {
            'B':'green', 
            'b':'green', 
            'A':'yellow', 
            'L':'red', 
            'l':'red', 
            'a':'yellow', 
            'o':'grey', 
            'X':'grey',
            '-':'grey'}
        

        #stereo colors
        """
        types = {
            'B':'green', 
            'b':'green', 
            'A':'orange', 
            'L':'orange', 
            'l':'orange', 
            'a':'orange', 
            'o':'orange', 
            'X':'orange',
            '-':'orange'}
        """
        """
        #only for continous chains!!
        
        with open('ramachandran.txt', 'r+') as f:
            raw_data=f.read().splitlines()
            if raw_data == 'error':
                print('error')
                raise IndexError

            data_rama = []
            for a, n in enumerate(raw_data):
                x=n.split()
                data_rama.append(x[-1])
        for x, n in enumerate(self.resvrange[1:-1]):
            cmd.color(types[data_rama[x]], 'all and resi %s'%(n))
        """

        for chain in self.cur_pdb().chain:
            for x, n in enumerate(chain.resi):
                if types[chain.ramach[x]] == 'green':
                    continue
                cmd.color(types[chain.ramach[x]], 'chain %s and resi %s'%(chain.letter, n))
            
        cmd.color('atomic', 'not elem C')

    elif method == 'recolor-sc' or method == 'recolor-area':
        #old method
        #recolor
        #can change between full and half hexa
        #colors all monomer!
        with open('acwi-out4.txt', 'r+') as f:
            raw_data=f.read().splitlines()
            if raw_data == 'error':
                print('error')
                raise IndexError

            data_sc = {}
            data_area = {}

            for n in raw_data:
                x=n.split()
                if not x[0].isdigit():
                    continue
                if len(x)>=5:
                    data_sc[x[0]]=x[2]
                    data_area[x[0]]=x[4]
                else:
                    data_sc[x[0]]=x[2]
                    data_area[x[0]]=x[3]

        cmd.color('atomic', 'all')
        cmd.color('gray', 'elem C')
        cmd.select('sele', 'none')
        chains = pdb_reader("monomer_0.pdb").splice()
        for chain in chains:
            x, y, z = chain.lista[0].xyz()
            cmd.select('sele', 'bychain sele or ( X>%s and X<%s and Y>%s and Y<%s and Z>%s and Z<%s)'%(x-0.005, x+0.005, y-0.005, y+0.005, z-0.005, z+0.005), 0, 1, 0, 1)

        cmd.alter('sele', 'b = -1')

        if method == 'recolor-sc':
            for x, n in enumerate(self.resvrange):
                try:
                    cmd.alter('sele and resi %s'%(n), 'b=%s'%(data_sc[str(n)]))
                except:
                    continue

            cmd.spectrum('b', 'rainbow_rev', 'sele', 0, 1, 0)
            #cmd.spectrum('b', 'rainbow_rev', 'sele')
            cmd.color('gray', 'sele and b = -1')

        elif method == 'recolor-area':
            for x, n in enumerate(self.resvrange):
                try:
                    cmd.alter('sele and resi %s'%(n), 'b=%s'%(data_area[str(n)]))
                except:
                    continue


            cmd.spectrum('b', 'rainbow_rev', 'sele', 0, 300, 0)
            #cmd.spectrum('b', 'rainbow_rev', 'sele')
            cmd.color('gray', 'sele and b = -1')


    elif method == 'sheets':
        #reads in from monomer_X
        
        cmd.color('atomic', 'all')
        try:
            for n in range(int(self.pdb_list.lista[self.pdb_count].protofilament)):
                

                cmd.select('sele_%s'%n, 'none')
                
                chains = pdb_reader("monomer_%s.pdb"%n).splice()
                for i in range(len(chains)):
                    x, y, z = chains[i].lista[0].xyz()
                    cmd.select('sele_%s'%n, 'bychain sele_%s or ( X>%s and X<%s and Y>%s and Y<%s and Z>%s and Z<%s )'%(n, x-0.005, x+0.005, y-0.005, y+0.005, z-0.005, z+0.005), 0, 1, 0, 1)
                    cmd.select('sele_%s'%n, 'sele_%s and polymer.protein'%n) 
                cmd.select('sele_%s'%n, 'sele_%s and polymer'%n)
                cmd.color(self.colors[n], 'sele_%s and elem C'%n)
        except:
            traceback.print_exc()
            pass

    elif method == 'protofilaments__':
        #using existing sele_X
        cmd.color('atomic', 'all')
        try:
            for n in range(int(self.pdb_list.lista[self.pdb_count].protofilament)):
                
                cmd.color(self.colors[n], 'sele_%s and elem C'%n)
                
        except:
            traceback.print_exc()
            pass

    elif method[:7] == 'recolor':
        #recolor-metric-number
        
        cmd.color('atomic', 'all')
        cmd.color('gray', 'elem C')
        cmd.select('sele', 'none')
        cmd.alter('all', 'b = -1')
        method = method.split('-')
        data = []
        
        inter = {'inter':'', 'intra':'i'}[self.color_inter]
        
        for pf in range(int(self.cur_pdb().protofilament)):
            
            with open('acwi-out%s%s_%s.txt'%(method[2], inter, pf), 'r+') as f:
                raw_data=f.read().splitlines()
                if raw_data == 'error':
                    print('error')
                    raise IndexError

                data_pf = {}
                
                for n in raw_data:
                    x=n.split()
                    
                    if not x[0].isdigit():
                        continue
                    
                    if method[1] == 'sc':
                        data_pf[x[0]]=x[2]
                    elif method[1] == 'area':
                        data_pf[x[0]]=x[4]
                    elif method[1] == 'sdi':
                        data_pf[x[0]]=x[5]
                
                for x, n in enumerate(self.cur_pdb().get_chain('pf', pf).resi):
                    if data_pf[str(n)] == '-':
                        continue
                    try:
                        cmd.alter('sele_%s and resi %s'%(pf, n), 'b=%s'%(data_pf[str(n)]))
                    except:
                        continue

            data.append(data_pf)
        
        if method[1] == 'sc':

            if self.color_range == 'absolute':
                print('Sc, absolute:')
                print('min: 0, max: 1')
                cmd.spectrum('b', 'rainbow', 'all and b > 0', 0, 1, 0)
            else:
                print('Sc, relative:')
                cmd.spectrum('b', 'rainbow', 'all and b > 0')
                 
        elif method[1] == 'area':

            if self.color_range == 'absolute':
                print('Area, absolute:')
                print('min: 50, max: 150')
                cmd.spectrum('b', 'rainbow', 'all and b > 0', 50, 150, 0)
            else:
                print('Area, relative:')
                cmd.spectrum('b', 'rainbow', 'all and b > 0')

        elif method[1] == 'sdi':

            if self.color_range == 'absolute':
                print('SDi, absolute:')
                print('min: 1, max: 2')
                cmd.spectrum('b', 'rainbow', 'all and b > 0', 1, 2, 0)
            else:
                print('SDi, relative:')
                cmd.spectrum('b', 'rainbow', 'all and b > 0')
                
        cmd.color('gray', 'all and b < 0')
        
        if self.color_range == 'relative':
            all_data = []
            for pf in data:
                all_data += list(map(float, list(filter('-'.__ne__, pf.values()))))
            print('min: %s, max:%s'%(round(min(all_data),2), round(max(all_data),2)))
        

    else:
        print('error')
    cmd.select('none')
    return


def atomic_area_prot(count, name, which_monomer, chain_count):
    """Calculates atomic area values, but no contribution."""

    midA = pdb_reader('areaimol_%s_sheet_%s_atomic.pdb'%(name, 0))
    midB = pdb_reader('areaimol_%s_sheet_%s_atomic.pdb'%(name, count))
    midAB = pdb_reader('areaimol_%s_sheet_0_%s_atomic.pdb'%(name, count)).sort_coord()

    #print('calc')
    #calculta area in place of B factor
    for atom in midB.lista:
        midA.add(atom)
    midA.sort_coord()

    try:
        for n, atom in enumerate(midAB.lista):
            if is_same(atom, midA.lista[n]):
                atom.b = str(round(float(midA.lista[n].b)-float(atom.b), 1))
            else:
                print('error')
                raise KeyError

        midAB.sort_float('ID')

    except KeyError:

        #old and slow
        midA = pdb_reader('areaimol_%s_sheet_%s_atomic.pdb'%(name, 0))
        midB = pdb_reader('areaimol_%s_sheet_%s_atomic.pdb'%(name, count))
        midAB = pdb_reader('areaimol_%s_sheet_0_%s_atomic.pdb'%(name, count))

        for atomA in midA.lista:
            for atomB in midAB.lista:
                if is_same(atomA, atomB):
                    atomB.b = str(round(float(atomA.b)-float(atomB.b), 1))
                    break
        for atomA in midB.lista:
            for atomB in midAB.lista:
                if is_same(atomA, atomB):
                    atomB.b = str(round(float(atomA.b)-float(atomB.b), 1))
                    break


    """
    #print('coloring')
    #coloring
    cmd.select('sele', 'none')
    cmd.select('atom', 'none')
    cmd.alter('all', 'b=0')
    area = 0
    for atom in midAB.lista:
        if atom.b == '0.0':
            continue
        x, y, z = atom.xyz()
        cmd.select('atom', 'X>%s and X<%s and Y>%s and Y<%s and Z>%s and Z<%s'%(x-0.005, x+0.005, y-0.005, y+0.005, z-0.005, z+0.005), 0, 1, 0, 1)
        cmd.select('sele', 'bychain sele or atom', 0, 1, 0, 1)
        cmd.alter('atom', 'b=%s'%(atom.b))
        area += float(atom.b)
    cmd.spectrum('b', 'rainbow_rev', 'sele')
    cmd.color('gray', 'b = 0')
    """

    #print('calc2')
    contacts = []
    b_sum = 0
    b_dict = {}
    for res in midAB.residues():
        ba = 0
        for atom in res.backbone(False).lista:
            if atom.b == '0.0':
                continue
            ba += float(atom.b)
            #b_sum += float(atom.b)
        """
        try:
            b_dict[lettercodes[res.lista[0].resn]+res.lista[0].resi] += ba
        except KeyError:
            b_dict[lettercodes[res.lista[0].resn]+res.lista[0].resi] = ba
        """
        bb = 0
        for atom in res.backbone(True).lista:
            if atom.b == '0.0':
                continue
            bb += float(atom.b)
            #b_sum += float(atom.b)
        """
        try:
            b_dict['B'+res.lista[0].resi] += bb
        except KeyError:
            b_dict['B'+res.lista[0].resi] = bb
        """
        if ba+bb > 0.6: #limit was originally 0
            #find and list contacts: residue number, monomer number, name, area aa, area bb
            contacts.append([res.lista[0].resi, which_monomer[res.lista[0].chain], lettercodes[res.lista[0].resn], round(ba, 3), round(bb, 3)])

    #contribution calc
    #ommited for speed
    contribution = ''
    """
    for aa in b_dict:
        if b_dict[aa] == 0:
            continue
        contribution += aa + '_' + str(round(b_dict[aa]/b_sum*100, 1))+'_'
    contribution += 'total_' + str(round(b_sum/2/chain_count, 1))
    """

    return contribution, contacts


def new_entry_analyser(self):
    
    cmd.remove('not polymer.protein')
    for n in range(int(self.cur_pdb().protofilament)):
        chains = pdb_reader('monomer_%s.pdb'%n).splice()
        for chain in chains:
            
            chain_data = amy_chain()
            chain_data.pf = n
            chain_data.letter = chain.lista[0].chain
            x,y,z = chains[0].lista[0].xyz()
            cmd.select('sele', 'none')
            cmd.select('sele', 'bychain sele or ( X>%s and X<%s and Y>%s and Y<%s and Z>%s and Z<%s)'%(x-0.005,x+0.005,y-0.005,y+0.005,z-0.005,z+0.005),0,1,0,1)
            chain_data.resi = get_data('(bychain sele) and name O ','resv')
            chain_data.seq = get_data('(bychain sele) and name O ','oneletter')
            
            self.cur_pdb().add_chain(chain_data)

    
def old_new_entry_analyser(resvrange):
    """For new polipeptide protein structures, dssp and sequence."""
    #where it is used? acwp?
    
    cmd.dss()

    with open('dss.txt', 'w+') as g:
        for x, n in enumerate(resvrange):
            line =  str(n)+' '+get_data('(bychain sele) and name O and resi %s'%(n), 'oneletter')[0]+' '
            dss=get_data('(bychain sele) and name O and resi %s'%(n), 'ss')[0]
            if dss == '':
                line += '-\n'
            else:
                line += dss+'\n'

            g.write(line)

    #sequence
    with open('seq.txt', 'w+') as g:
        line = ''
        for x, n in enumerate(resvrange):
            line += get_data('(bychain sele) and name O and resi %s'%(n), 'oneletter')[0]

        g.write(line)


def section_maker(cc_d, out, start , first = True):
    #algorithm find the connected section in a graph of interface interactions
    
    try:
        current = cc_d[start]
        cc_d.pop(start)
    except KeyError:
        return cc_d, out
    out[-1].append(start)
    
    for c in current:
        cc_d, out = section_maker(cc_d, out, c, False)
    
    if first and len(cc_d) > 0 :
        out.append([])
        cc_d, out = section_maker(cc_d, out, list(cc_d)[0])
      
    return cc_d, out



def find_connected_sections():



    lista = [] #all nodes
    lista_main = [] #all nodes in main chain
    listapar = [] #interactions
    listapar2 = [] #mainchain backbone
    listapar3 = [] #close interactions in backbone
    main_contacts = []


    
    try:
        #raise KeyError()
        with open('contacts.txt', 'r+') as g:
            raw_data=g.read().splitlines()
            for line in raw_data:
                lista_main.append(line.split()[0]+'_0')
                main_contacts.append([])
                for n in line.split()[1:]:
                    lista.append(n)
                    main_contacts[-1].append(n)
    except Exception as error:
        print(error)
        return

    #beolvasas
    for x, aa in enumerate(lista_main):
        for contact in main_contacts[x]:
            if contact.split('_')[1] != '0':
                if int(aa.split('_')[0])- int(contact.split('_')[0]) > 0:
                    listapar.append(contact + ' ' + aa)
                else:
                    listapar.append(aa + ' ' + contact) 
            elif abs(int(aa.split('_')[0])- int(contact.split('_')[0])) > 2:
                if int(aa.split('_')[0])- int(contact.split('_')[0]) > 0:
                    listapar.append(contact + ' ' + aa)
                else:
                    listapar.append(aa + ' ' + contact)
            else:
                if int(aa.split('_')[0])- int(contact.split('_')[0]) > 0:
                    listapar3.append(contact + ' ' + aa)
                else:
                    listapar3.append(aa + ' ' + contact)
    
    
    for x, aa in enumerate(lista_main[0:-1]):
        if int(aa.split('_')[0])+1  == int(lista_main[x+1].split('_')[0]):
            listapar2.append([aa,lista_main[x+1]])
        
        
    lista = list(dict.fromkeys(lista))
    
    #listapar = list(dict.fromkeys(listapar))
    
    
    listapar = list(map(lambda x: x.split() , listapar))
    listapar.sort(key = lambda x: int(x[0].split('_')[0]))
    
    #find double connections 
    cc = [] #close contacts
    for n, par in enumerate(listapar[:-1]):
        if par in listapar[n+1:]:
            cc.append(par)
    
    #create dictionary and get rid of monomer marker
    cc_d = {} #close contact dictionary
    for c in cc:
        try:
            cc_d[c[0].split("_")[0]].append(c[1].split("_")[0])
        except:
            cc_d[c[0].split("_")[0]] = [c[1].split("_")[0]]
        try:
            cc_d[c[1].split("_")[0]].append(c[0].split("_")[0])
        except:
            cc_d[c[1].split("_")[0]] = [c[0].split("_")[0]]
    
        
    #break down into connected sections with recursion
    if len(cc_d):
        cc_d, out = section_maker(cc_d,[[]],list(cc_d)[0])
        out = list(map(lambda x: sorted(x, key = int),out))
    else:
        with open('connected_section.txt', 'w+') as f:
            f.write("")
            return
    
    #print(out)
    #complete with missing aa and delete small sections
    sections = []
    for section in out:
        if len(section)<4:
            continue
        
        sections.append([])
        for x, aa in enumerate(section):
            sections[-1].append(aa)
            try:
                if str(int(aa)+2) == section[x+1]:
                    sections[-1].append(str(int(aa)+1))
                elif str(int(aa)+3) == section[x+1]:
                    sections[-1].append(str(int(aa)+1))
                    sections[-1].append(str(int(aa)+2))
            except IndexError:
                pass
    
    sections_final=[]
    #cut in half
    for section in sections:
        continuous = True
        for x, aa in enumerate(section[:-1]):
            #if there is a gap
            if str(int(aa)+1) != section[x+1]:
                sections_final.append(deepcopy(section))
                sections_final[-1].insert(x+1,'-')
                continuous = False
        if continuous:
            length = len(section)
            #delete the middle 2/3
            if length % 2 == 0:
                sections_final.append(deepcopy(section))
                sections_final[-1].pop(int(length/2-1))
                sections_final[-1].pop(int(length/2-1))
                sections_final[-1].insert(int(length/2-1),'-')
            else:
                sections_final.append(deepcopy(section))
                sections_final[-1].pop(int(length/2-1))
                sections_final[-1].pop(int(length/2-1))
                sections_final[-1].pop(int(length/2-1))
                sections_final[-1].insert(int(length/2-1),'-')
    
    #print(sections_final)
    
    with open('connected_section.txt', 'w+') as f:
        for section in sections_final:
            line = ""
            for n in section:
                line += n+'_'
            f.write(line[:-1]+'\n')






#notes for calc methods

#areaimol_input for acw, acwperf, f1,f2,f3,f4,f5
#areaimol_calc_1 for acw, acwperf, f1,f2
#areaimol_calc_2 for f3,f4,f5
#areaimol_out_1 for acw, acwperf
#areaimol_out_2 for f1, f2
#areaimol_out_3 for f3, f4, f5
#sc_input for all
#sc_calc for all
#sc_out for all
#sc_calc_out subprocess style sc
#areaimol_calc_3 subprocess style area calc_1

#A-D
#B-E
#C-F

#A-G
#B-H
#C-I
#D-J
#E-K
#F-L


#A-K
#B-L
#C-M
#D-N
#E-O
#F-P
#G-Q
#H-R
#I-S
#J-T







def areaimol_calc_1(count, name, path, atomic=False):
    """Write and execute Areaimol calculation via CCP4 trough bash."""
    #overwritten by areaimol_calc_3, old version in development


    if atomic:
        output = 'OUTPUT\n'
    else:
        output = ''

    if platform.system() == 'Windows':
        p = subprocess.Popen(['cmd.exe'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
        p.stdin.write('%s/ccp4.setup.bat\n'%path)
        result = p.communicate('%s/bin/areaimol.exe XYZIN %s_slicedsheet_%s.pdb > areaimol_%s_sheet_%s.txt\n%sEND\n'%(path, name, count, name, count, output))
    elif platform.system() == 'Linux':
        p = subprocess.Popen(['/bin/sh'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        result = p.communicate('%s/bin/areaimol XYZIN %s_slicedsheet_%s.pdb > areaimol_%s_sheet_%s.txt\n%sEND\n'%(path, name, count, name, count, output))
    if atomic:
        shutil.move('XYZOUT', 'areaimol_%s_sheet_%s_atomic.pdb'%(name, count))

    if count == 0:
        pass
    else:
        if platform.system() == 'Windows':
            p = subprocess.Popen(['cmd.exe'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
            p.stdin.write('%s/ccp4.setup.bat\n'%path)
            result = p.communicate('%s/bin/areaimol.exe XYZIN %s_slicedsheet_0_%s.pdb > areaimol_%s_sheet_0_%s.txt\n%sEND\n'%(path, name, count, name, count, output))
        elif platform.system() == 'Linux':
            p = subprocess.Popen(['/bin/sh'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
            result = p.communicate('%s/bin/areaimol XYZIN %s_slicedsheet_0_%s.pdb > areaimol_%s_sheet_0_%s.txt\n%sEND\n'%(path, name, count, name, count, output))
        if atomic:
            shutil.move('XYZOUT', 'areaimol_%s_sheet_0_%s_atomic.pdb'%(name, count))


def areaimol_calc_2(count, name, path, atomic=False):
    # write and execute Areaimol calculation via CCP4

    if atomic:
        output = 'OUTPUT\n'
    else:
        output = ''

    if platform.system() == 'Windows':
        p = subprocess.Popen(['cmd.exe'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
        p.stdin.write('%s/ccp4.setup.bat\n'%path)
        result = p.communicate('%s/bin/areaimol.exe XYZIN Sheet_%s.pdb > areaimol_%s_sheet_%s.txt\n%sEND\n'%(path, count, name, count, output))
    elif platform.system() == 'Linux':
        p = subprocess.Popen(['/bin/sh'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        result = p.communicate('%s/bin/areaimol XYZIN Sheet_%s.pdb > areaimol_%s_sheet_%s.txt\n%sEND\n'%(path, count, name, count, output))
    if atomic:
        shutil.move('XYZOUT', 'areaimol_%s_sheet_%s_atomic.pdb'%(name, count))

    if count == 0:
        pass
    else:
        if platform.system() == 'Windows':
            p = subprocess.Popen(['cmd.exe'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
            p.stdin.write('%s/ccp4.setup.bat\n'%path)
            result = p.communicate('%s/bin/areaimol.exe XYZIN Sheet_0-%s.pdb > areaimol_%s_sheet_0_%s.txt\n%sEND\n'%(path, count, name, count, output))
        elif platform.system() == 'Linux':
            p = subprocess.Popen(['/bin/sh'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
            result = p.communicate('%s/bin/areaimol XYZIN Sheet_0-%s.pdb > areaimol_%s_sheet_0_%s.txt\n%sEND\n'%(path, count, name, count, output))
        if atomic:
            shutil.move('XYZOUT', 'areaimol_%s_sheet_0_%s_atomic.pdb'%(name, count))


def areaimol_out_2(count, name, sheet_number):
    # read the output file to get areaimol and area
    #middle chain from each sheet

    area1 = -1
    area1X = -1
    areaX = -1
    carea1 = -1
    carea1X = -1
    careaX = -1
    print('Total area of chain %s'%ch_n(int(sheet_number/2+1)-1))
    with open("areaimol_%s_sheet_%s.txt"%(name, count), 'r+') as f:
        for line in f.read().splitlines():
            if line[0:21] == 'Total area of chain %s'%ch_n(int(sheet_number/2+1)-1):
                areaX = float(line.split()[-1])
            if line[0:29] == 'Total contact area of chain %s'%ch_n(int(sheet_number/2+1)-1):
                careaX = float(line.split()[-1])

    with open("areaimol_%s_sheet_0_%s.txt"%(name, count), 'r+') as f:
        for line in f.read().splitlines():
            if line[0:21] == 'Total area of chain %s'%ch_n(int(sheet_number/2+1)-1):
                area1X = float(line.split()[-1])
            if line[0:29] == 'Total contact area of chain %s'%ch_n(int(sheet_number/2+1)-1):
                carea1X = float(line.split()[-1])
            if line[0:21] == 'Total area of chain %s'%ch_n(int(sheet_number/2+1)+sheet_number-1):
                area1X = area1X+float(line.split()[-1])
            if line[0:29] == 'Total contact area of chain %s'%ch_n(int(sheet_number/2+1)+sheet_number-1):
                carea1X = carea1X+float(line.split()[-1])

    with open("areaimol_%s_sheet_0.txt"%(name), 'r+') as f:
        for line in f.read().splitlines():
            if line[0:21] == 'Total area of chain %s'%ch_n(int(sheet_number/2+1)-1):
                area1 = float(line.split()[-1])
            if line[0:29] == 'Total contact area of chain %s'%ch_n(int(sheet_number/2+1)-1):
                carea1 = float(line.split()[-1])

    return (area1, area1X, areaX, carea1, carea1X, careaX)


def areaimol_out_3(count, name, sheet_number):
    # read the output file to get areaimol and area
    # total area

    area1 = -1
    area1X = -1
    areaX = -1
    carea1 = -1
    carea1X = -1
    careaX = -1

    with open("areaimol_%s_sheet_%s.txt"%(name, count), 'r+') as f:
        for line in f.read().splitlines():
            if line[0:21] == 'Total overall area   ':
                areaX = float(line.split()[-1])
            if line[0:29] == 'Total overall contact area   ':
                careaX = float(line.split()[-1])

    with open("areaimol_%s_sheet_0_%s.txt"%(name, count), 'r+') as f:
        for line in f.read().splitlines():
            if line[0:21] == 'Total overall area   ':
                area1X = float(line.split()[-1])
            if line[0:29] == 'Total overall contact area   ':
                carea1X = float(line.split()[-1])

    with open("areaimol_%s_sheet_0.txt"%(name), 'r+') as f:
        for line in f.read().splitlines():
            if line[0:21] == 'Total overall area   ':
                area1 = float(line.split()[-1])
            if line[0:29] == 'Total overall contact area   ':
                carea1 = float(line.split()[-1])

    return (area1, area1X, areaX, carea1, carea1X, careaX)


def areaimol_out_4(count, name, which_monomer):
    #areaimol out from atomic area

    midA = pdb_reader('areaimol_%s_sheet_%s_atomic.pdb'%(name, 0))
    midB = pdb_reader('areaimol_%s_sheet_%s_atomic.pdb'%(name, count))
    midAB_all = pdb_reader('areaimol_%s_sheet_0_%s_atomic.pdb'%(name, count))

    mid_obj = pdb_reader('Sheet_midobj.pdb')
    mid_obj_dict = dict()
    for atom in mid_obj.lista:
        mid_obj_dict[atom.xyz_key()] = atom
    
    midAB_dict = dict()
    for atom in midAB_all.lista:
        if atom.xyz_key() in mid_obj_dict:
            midAB_dict[atom.xyz_key()] = atom

    res_dictA = {}
    for atom in midA.lista:
        if atom.xyz_key() in midAB_dict:
            midAB_dict[atom.xyz_key()].b = atom.b - midAB_dict[atom.xyz_key()].b 
            try:
                res_dictA[(which_monomer[atom.chain], atom.resi, atom.resn)] += midAB_dict[atom.xyz_key()].b
            except KeyError:
                res_dictA[(which_monomer[atom.chain], atom.resi, atom.resn)] = midAB_dict[atom.xyz_key()].b
            
    res_dictB = {}
    for atom in midB.lista:
        if atom.xyz_key() in midAB_dict:
            midAB_dict[atom.xyz_key()].b = atom.b - midAB_dict[atom.xyz_key()].b
            try:
                res_dictB[(which_monomer[atom.chain], atom.resi, atom.resn)] += midAB_dict[atom.xyz_key()].b
            except KeyError:
                res_dictB[(which_monomer[atom.chain], atom.resi, atom.resn)] = midAB_dict[atom.xyz_key()].b

    area = 0
    for atom in midAB_dict.values():
        area += float(atom.b)
    """
    midAB = mAtom()
    for atom in midAB_dict.values():
        midAB.add(atom)
    select_model(midAB, 'midAB')
    """
    
    fingerprint = [[], []]
    for res in res_dictA:
        if res_dictA[res] > 10 or (res_dictA[res] > 5 and res[2] == 'GLY'):
            fingerprint[0].append(res)
    for res in res_dictB:
        if res_dictB[res] > 10 or (res_dictB[res] > 5 and res[2] == 'GLY'):
            fingerprint[1].append(res)
    
    return area/2.0, fingerprint



def areaimol_out_4_atomic(count, name, which_monomer, chain_count):
    #areaimol out from atomic area
    
    midA = pdb_reader('areaimol_%s_sheet_%s_atomic.pdb'%(name, 0))
    midB = pdb_reader('areaimol_%s_sheet_%s_atomic.pdb'%(name, count))
    midAB_all = pdb_reader('areaimol_%s_sheet_0_%s_atomic.pdb'%(name, count))

    mid_obj = pdb_reader('Sheet_midobj.pdb')
    mid_obj_dict = dict()
    for atom in mid_obj.lista:
        mid_obj_dict[atom.xyz_key()] = atom
    
    midAB_dict = dict()
    for atom in midAB_all.lista:
        if atom.xyz_key() in mid_obj_dict:
            midAB_dict[atom.xyz_key()] = atom

    for atom in midA.lista:
        if atom.xyz_key() in midAB_dict:
            midAB_dict[atom.xyz_key()].b = atom.b - midAB_dict[atom.xyz_key()].b 
    for atom in midB.lista:
        if atom.xyz_key() in midAB_dict:
            midAB_dict[atom.xyz_key()].b = atom.b - midAB_dict[atom.xyz_key()].b
    
    midAB = mAtom()
    for atom in midAB_dict.values():
        midAB.add(atom)
    midAB.sort_float('ID')

    contacts = []
    b_sum = 0
    b_dict = {}
    for res in midAB.residues():
        ba = 0
        for atom in res.backbone(False).lista:
            if atom.b == '0.0':
                continue
            ba += float(atom.b)
            #b_sum += float(atom.b)
        """
        try:
            b_dict[lettercodes[res.lista[0].resn]+res.lista[0].resi] += ba
        except KeyError:
            b_dict[lettercodes[res.lista[0].resn]+res.lista[0].resi] = ba
        """
        bb = 0
        for atom in res.backbone(True).lista:
            if atom.b == '0.0':
                continue
            bb += float(atom.b)
            #b_sum += float(atom.b)
        """
        try:
            b_dict['B'+res.lista[0].resi] += bb
        except KeyError:
            b_dict['B'+res.lista[0].resi] = bb
        """
        if ba+bb > 0.6: #limit was originally 0
            #find and list contacts: residue number, monomer number, name, area aa, area bb
            contacts.append([res.lista[0].resi, which_monomer[res.lista[0].chain], lettercodes[res.lista[0].resn], round(ba, 3), round(bb, 3)])

    #contribution calc
    #ommited for speed
    contribution = ''
    """
    for aa in b_dict:
        if b_dict[aa] == 0:
            continue
        contribution += aa + '_' + str(round(b_dict[aa]/b_sum*100, 1))+'_'
    contribution += 'total_' + str(round(b_sum/2/chain_count, 1))
    """

    return contribution, contacts


def perf_calc_chains(count, sheet_count, to_generate):
    
    modelA = pdb_reader('Sheet_0.pdb')
    modelB = pdb_reader('Sheet_%s.pdb'% str(count))
    
    if len(modelA.splice(sort=False)) == sheet_count:
        model1 = modelA.splice()
    else:
        model1 = modelA.splice_equal(sheet_count)
    if len(modelB.splice(sort=False)) == sheet_count:
        model2 = modelB.splice()
    else:
        model2 = modelB.splice_equal(sheet_count)
    
    start_n = int(to_generate//2*sheet_count/to_generate)
    end_n = int(to_generate//2*sheet_count/to_generate + sheet_count/to_generate)
    
    output = [] 
    for n in range(start_n, end_n):
        
        select_model(model1[n-1: n+2], 'calc_0_c', output='Sheet_0.pdb')
        select_model(model2[n-1: n+2], 'calc_1_c', output='Sheet_%s.pdb'% str(count))
        #output: length, heigth, length_coord
        output.append(perf_calc_fast(count, 3, S='calc_0_c', B='calc_1_c'))

    return output
    

def perf_calc_fast(count, sheet_count, S=None, B=None):
    """ """

    #import acwperf
    #cleverly find points of possible water positions
    #adds a circle of points to the surface points, then with union finds the output
    
    #checks the orientation
    if S is not None:
        v_cross, angle = straigthen_pca_check(S, B, sheet_count, count=count)
    else:
        angle = None

    name='tmp'

    modelA = pdb_reader('Sheet_0.pdb')
    modelB = pdb_reader('Sheet_%s.pdb'% str(count))

    if len(modelA.splice(sort=False)) == sheet_count:
        model1 = modelA.splice()[int(sheet_count/2.0)]
        modelA_spliced = modelA.splice()
    else:
        model1 = modelA.splice_equal(sheet_count)[int(sheet_count/2.0)]
        modelA_spliced = modelA.splice_equal(sheet_count)
    if len(modelB.splice(sort=False)) == sheet_count:
        model2 = modelB.splice()[int(sheet_count/2.0)]
    else:
        model2 = modelB.splice_equal(sheet_count)[int(sheet_count/2.0)]

    axis = find_axis_PCA_spliced(modelA_spliced)
    project = {'x':0,'y':1,'z':2}[axis]
    avg_heigth = modelA.average(axis)
    heigth = get_length_coord(modelA_spliced[0].center(), modelA_spliced[1].center())

    flat1 = flatten(model1, axis)
    flat2 = flatten(model2, axis)

    radius = 3.2 # 1.4+1.8=3.2
    limit = 0.01 #rounding limit was set to 0.05   

    grid_1d = np.arange(-1*radius - limit, radius + limit*2, limit)
    X, Y = np.meshgrid(grid_1d, grid_1d)
    grid_points = np.column_stack((X.ravel(), Y.ravel()))
    grid_distance = np.sum(np.square(grid_points), axis=1)
    mask_d = grid_distance > (radius-limit)**2
    mask_D = grid_distance < (radius+limit)**2
    mask = np.all([mask_d, mask_D], axis=0)

    circle = grid_points[mask]
    
    mask_radius = (radius - limit)**2
    #for each surface generate the outer surface points
    surface1 = np.empty((0,2))
    for atom1 in flat1.lista:
        x1 = round(atom1.x * 100) / float(100)
        y1 = round(atom1.y * 100) / float(100)
        z1 = round(atom1.z * 100) / float(100)
        c_atom1 = np.delete(np.array([x1, y1, z1]), project)
        closest = []
        for atom2 in flat1.lista: #find close atoms
            d = ((x1-atom2.x)**2 + (y1-atom2.y)**2 + (z1-atom2.z)**2 )** 0.5
            if d < 2*radius+limit:
                closest.append(atom2.xyz(axis))
                closest[-1].append(d)
        
        closest.sort(key=lambda x:x[2])
        closest = np.array([x[0:2] for x in closest])
        shifted_circle = circle + c_atom1 #put the cicrle around the atom
        
        for point in closest:
            #print(shifted_circle.shape)
            #distance = np.linalg.norm(shifted_circle - point, axis=1)
            #distance = np.sqrt(np.sum((shifted_circle - point) ** 2, axis=0))
            a_min_b = shifted_circle - point
            distance = np.einsum("ij,ij->i", a_min_b, a_min_b)
            mask = distance > mask_radius
            shifted_circle = shifted_circle[mask]
            if shifted_circle.size == 0:
                break
        else:
            surface1 = np.append(surface1, shifted_circle, axis=0)

    surface2 = np.empty((0,2))
    for atom1 in flat2.lista:
        x1 = round(atom1.x * 100) / float(100)
        y1 = round(atom1.y * 100) / float(100)
        z1 = round(atom1.z * 100) / float(100)
        c_atom1 = np.delete(np.array([x1, y1, z1]), project)
        closest = []
        for atom2 in flat2.lista: #find close atoms
            d = ((x1-atom2.x)**2 + (y1-atom2.y)**2 + (z1-atom2.z)**2 )** 0.5
            if d < 2*radius+limit:
                closest.append(atom2.xyz(axis))
                closest[-1].append(d)
        
        closest.sort(key=lambda x:x[2])
        closest = np.array([x[0:2] for x in closest])
        shifted_circle = circle + c_atom1 #put the cicrle around the atom
        
        for point in closest:
            a_min_b = shifted_circle - point
            distance = np.einsum("ij,ij->i", a_min_b, a_min_b)
            mask = distance > mask_radius
            shifted_circle = shifted_circle[mask]
            if shifted_circle.size == 0:
                break
        else:
            surface2 = np.append(surface2, shifted_circle, axis=0)
    
    tangents = []
    surface1 = set([tuple(x) for x in np.round(surface1,2)])
    surface2 = set([tuple(x) for x in np.round(surface2,2)])
    for point in surface1: #find common points
        if point in surface2:
            tangents.append(point)
    tangents = list(set(tangents))

    """
    #with savefile
    dist = []
    with open("perf_coord_%s_fast.txt" %str(count), 'w+') as f:
        for n, coord1 in enumerate(tangents):
            atom1 = sAtom(name, coord1[0], coord1[1], coord1[2])
            f.write(str(atom1.x)+' '+str(atom1.y)+' '+str(atom1.z)+'\n')
            for coord2 in tangents[n+1:]:
                atom2 = sAtom(name, coord2[0], coord2[1], coord2[2])
                dist.append(get_length(atom1, atom2, 'lengths'))

        if len(tangents) == 0:
            length = 0
        else:
            length = max(dist)-2.8  #testing 2.8
        f.write(str(round(length, 3)))
    """

    #with output
    dist = []
    dist_c = []
    for n, coord1 in enumerate(tangents):
        coord1 = list(coord1)
        coord1.insert(project, 0)
        atom1 = sAtom(name, coord1[0], coord1[1], coord1[2])
        #(str(atom1.x)+' '+str(atom1.y)+' '+str(atom1.z)+'\n')
        for coord2 in tangents[n+1:]:
            coord2 = list(coord2)
            coord2.insert(project, 0)
            atom2 = sAtom(name, coord2[0], coord2[1], coord2[2])
            dist.append(get_length(atom1, atom2))
            dist_c.append([atom1, atom2])

    if len(tangents) == 0:
        length = 0
        length_coord = ""
    else:
        #length = max(dist)-2.8  #two times water 1.4*2

        atom1, atom2 = dist_c[dist.index(max(dist))]
        #cmd.pseudoatom('water2', pos=atom2.xyz())
        for atom in (atom1, atom2):
            flat1_c = coord_reader(vector_segment(atom, get_nearest(atom, flat1, True)[2], 1.4))
            flat2_c = coord_reader(vector_segment(atom, get_nearest(atom, flat2, True)[2], 1.4))
            contact_p = coord_reader(vector_segment(atom, get_midpoint(flat1_c, flat2_c), 1.4))
            for xyz in ('x', 'y', 'z'):
                setattr(atom, xyz, getattr(contact_p, xyz))
        
        #cmd.pseudoatom('contact_p', pos=contact_p.xyz())
        #cmd.pseudoatom('flat1_c', pos=flat1_c.xyz())
        #cmd.pseudoatom('flat2_c', pos=flat2_c.xyz())
        length = get_length(atom1, atom2)
        
        #save coordinates
        setattr(atom1, axis, avg_heigth)
        setattr(atom2, axis, avg_heigth)
        cmd.pseudoatom('tmp1', pos=atom1.xyz())
        cmd.pseudoatom('tmp2', pos=atom2.xyz())
        if angle is not None:
            cmd.rotate(v_cross, angle=-angle, selection='tmp1 or tmp2', camera=0, object=None, origin=[0, 0, 0])

        atom1 = coord_reader(get_data('tmp1','[x, y, z]')[0], name)
        atom2 = coord_reader(get_data('tmp2','[x, y, z]')[0], name)
        length_coord = '%s_%s_%s_%s_%s_%s'%(round(atom1.x, 3), round(atom1.y, 3), round(atom1.z, 3),
                                            round(atom2.x, 3), round(atom2.y, 3), round(atom2.z, 3))

        #drawing a line
        cmd.distance('lengths','tmp1', 'tmp2')
        cmd.delete('tmp1')
        cmd.delete('tmp2')
        cmd.hide('label', 'lengths')

    return length, heigth, length_coord


def cgoCircle(x, y, z, r=8.0, cr=1.0, cg=0.4, cb=0.8, w=2.0):
    """
    Create a CGO circle

    PARAMS
        x, y, z
            X, Y and Z coordinates of the origin

        r
            Radius of the circle

        cr, cg, cb
            Color triplet, [r,g,b] where r,g,b are all [0.0,1.0].

        w
            Line width of the circle

    RETURNS
        the CGO object (it also loads it into PyMOL, too).

    """
    x = float(x)
    y = float(y)
    z = float(z)
    r = abs(float(r))
    cr = abs(float(cr))
    cg = abs(float(cg))
    cb = abs(float(cb))
    w = float(w)

    obj = [ BEGIN, LINES, COLOR, cr, cg, cb ]
    for i in range(180):
        obj.append( VERTEX )
        obj.append(r*math.cos(i) + x )
        obj.append(r*math.sin(i) + y )
        obj.append(z)
        obj.append( VERTEX )
        obj.append(r*math.cos(i+0.1) + x )
        obj.append(r*math.sin(i+0.1) + y )
        obj.append(z)
    obj.append(END)

    cName = cmd.get_unused_name("circle_")
    cmd.load_cgo( obj, cName )
    cmd.set("cgo_line_width", w, cName )
    return obj



def lstsq_cylinder_calc(self, show=False):
    """Least square fit a cyclynder around all CA with various radii"""
    
    print('started')
    centers = []
    first_gues = [[],[]]
    
    for monomer in range(int(self.cur_pdb().protofilament)):
        
        cmd.save('Sheet_0.pdb','sele_%s and polymer'%monomer)  
        p_filament = pdb_reader("Sheet_0.pdb").splice()
        first_gues[0].append(p_filament[0].center())
        first_gues[1].append(p_filament[-1].center())
        
        for resi in self.cur_pdb().get_chain('pf', monomer).resi:
            """
            #show selected atoms 
            test = mAtom()
            for n in range(len(p_filament)):
                test.add(p_filament[n].select('resi',[str(resi)]).select('name', ['CA']).lista[0])
            select_model(test, 'test')
            """
            
            points = []
            for n in range(len(p_filament)):
                xyz = p_filament[n].select('resi',[str(resi)]).select('name', ['CA']).lista[0].xyz()
                points.append(xyz)
            
            centers.append(points)
            
    centers = np.array(centers)
    centers = np.transpose(centers, (1, 2, 0))
        
    #x0 = [100., 100., 0., 100., 100., 100.] + [2.]*centers.shape[2]
    start = [sum(col) / len(col) for col in zip(*first_gues[0])]
    end = [sum(col) / len(col) for col in zip(*first_gues[1])]
    
    x0 = start + end + [10.]*centers.shape[2]
    #x0 = p_filament[0].center() + p_filament[-1].center() + [15.]*centers.shape[2]
    bounds = ([-1000., -1000., -1000., -1000., -1000., -1000.] + [0.05]*centers.shape[2], [1000., 1000., 1000., 1000., 1000., 1000.] + [500.]*centers.shape[2])
    
    p = scipy.optimize.least_squares(lstsq_cylinder_func, x0=x0, bounds=bounds, args=(centers,))
    
    start = p.x[0:3]
    end = p.x[3:6]
    rmax = np.max(p.x[6:])
    
    if show == 'sample':
        for x in range(6, len(p.x), 10):
            #r = np.max(p.x[6:])
            r = p.x[x]
            #print(p.x)
            cmd.load_cgo( [CYLINDER] + list(start) + list(end) + [r, 1,1,0, 1,1,0], "cylinder%s"%x )
            #cgoCircle(x, y, z, r=8.0, cr=1.0, cg=0.4, cb=0.8, w=2.0) #test circle_donut
    elif show == 'max': 
        cmd.load_cgo( [CYLINDER] + list(start) + list(end) + [rmax, 1,1,0, 1,1,0], "cylinder%s"%x )
    
    return p.x.tolist(), p.cost, p.status


def lstsq_cylinder_func(p, xyz):
    # Extract the start and end points of the common axis
    start_point = np.array(p[0:3])  # Shape (3,)
    end_point = np.array(p[3:6])    # Shape (3,)
    radii = np.array(p[6:6 + xyz.shape[2]])  # Array of radii for each cylinder (shape N,)

    # Calculate the axis vector and its norm (length)
    axis_vector = end_point - start_point
    axis_length = np.linalg.norm(axis_vector)

    # Reshape start_point to (1, 3, 1) so it can broadcast with xyz
    start_point = start_point.reshape(1, 3, 1)
    end_point = end_point.reshape(1, 3, 1)

    # Calculate the cross product between each point vector and the axis vector
    # This gives the perpendicular vector from each point to the axis
    cross_product = np.cross(xyz - start_point, xyz - end_point, axis=1)  # Shape (P, 3, N)

    # Calculate the norm (magnitude) of the cross product to get the area of the parallelogram
    cross_norm = np.linalg.norm(cross_product, axis=1)  

    # Divide by the length of the axis vector to get the perpendicular distance
    distances = cross_norm / axis_length  

    # Square the distances and subtract the squared radii for each cylinder
    residuals = (np.power(distances, 2) - np.power(radii, 2)) * distances
    #original: residuals = np.power(distances, 2) - np.power(radii, 2) 
    #residuals = residuals * ((radii - 3) / np.sqrt(1 + np.power(radii-3, 2)) + 2)

    return np.ravel(residuals)


def lstsq_cylinder_calc_axis(self, show=False):
    """Least square fit a cyclynder around all CA with various radii, fixed z axis"""
    
    print('started')
    centers = []
    rs = mAtom()
    first_gues = [[],[]]
    
    for monomer in range(int(self.cur_pdb().protofilament)):
        cmd.save('Sheet_0.pdb','sele_%s and polymer'%monomer)  
        p_filament = pdb_reader("Sheet_0.pdb").splice()
        first_gues[0].append([float(self.cur_pdb().map_x), float(self.cur_pdb().map_y), min(p_filament[0].data('z'))])
        first_gues[1].append([float(self.cur_pdb().map_x), float(self.cur_pdb().map_y), max(p_filament[-1].data('z'))])
        
        rs.addmodel(p_filament[0])

        for resi in self.cur_pdb().get_chain('pf', monomer).resi:
            """
            #show selected atoms 
            test = mAtom()
            for n in range(len(p_filament)):
                test.add(p_filament[n].select('resi',[str(resi)]).select('name', ['CA']).lista[0])
            select_model(test, 'test')
            """
            
            points = []
            for n in range(len(p_filament)):
                xyz = p_filament[n].select('resi',[str(resi)]).select('name', ['CA']).lista[0].xyz()
                points.append(xyz)
            
            centers.append(points)
            
    centers = np.array(centers)
    centers = np.transpose(centers, (1, 2, 0))
        
    #x0 = [100., 100., 0., 100., 100., 100.] + [2.]*centers.shape[2]
    start = [sum(col) / len(col) for col in zip(*first_gues[0])]
    end = [sum(col) / len(col) for col in zip(*first_gues[1])]
    
    rs = rs.select('name', ['CA']).coord_list()
    rs = np.reshape(rs, (-1, 3))
    x0 = start[0:2] + list(get_dist_3d_line(np.array(start), np.array(end), rs)) #[10.]*centers.shape[2]

    #x0 = p_filament[0].center() + p_filament[-1].center() + [15.]*centers.shape[2]
    bounds = ([-1000., -1000.] + [0.0005]*centers.shape[2], [1000., 1000.] + [500.]*centers.shape[2])

    p = scipy.optimize.least_squares(lstsq_cylinder_func_axis, x0=x0, bounds=bounds, args=(centers,))
    
    p.x = np.array(list(p.x[0:2]) + [0.] + list(p.x[0:2]) + [100.0] + list(p.x[2:])) 
    
    start = p.x[0:3]
    end = p.x[3:6]
    rmax = np.max(p.x[6:])
    
    if show == 'sample':
        for x in range(6, len(p.x), 10):
            #r = np.max(p.x[6:])
            r = p.x[x]
            #print(p.x)
            cmd.load_cgo( [CYLINDER] + list(start) + list(end) + [r, 1,1,0, 1,1,0], "cylinder%s"%x )
            #cgoCircle(x, y, z, r=8.0, cr=1.0, cg=0.4, cb=0.8, w=2.0) #test circle_donut
    elif show == 'max': 
        cmd.load_cgo( [CYLINDER] + list(start) + list(end) + [rmax, 1,1,0, 1,1,0], "cylinder%s"%x )
    
    return p.x.tolist(), p.cost, p.status


def lstsq_cylinder_func_axis(p, xyz):
    # Extract the start and end points of the common axis
    start_point = np.append(p[0:2], 0.)  # Shape (3,)
    end_point = np.append(p[0:2], 100.)    # Shape (3,)
    radii = np.array(p[2:2 + xyz.shape[2]])  # Array of radii for each cylinder (shape N,)

    # Calculate the axis vector and its norm (length)
    axis_vector = end_point - start_point
    axis_length = np.linalg.norm(axis_vector)

    # Reshape start_point to (1, 3, 1) so it can broadcast with xyz
    start_point = start_point.reshape(1, 3, 1)
    end_point = end_point.reshape(1, 3, 1)

    # Calculate the cross product between each point vector and the axis vector
    # This gives the perpendicular vector from each point to the axis
    cross_product = np.cross(xyz - start_point, xyz - end_point, axis=1)  # Shape (P, 3, N)

    # Calculate the norm (magnitude) of the cross product to get the area of the parallelogram
    cross_norm = np.linalg.norm(cross_product, axis=1)  

    # Divide by the length of the axis vector to get the perpendicular distance
    distances = cross_norm / axis_length  

    # Square the distances and subtract the squared radii for each cylinder
    residuals = (np.power(distances, 2) - np.power(radii, 2)) * distances
    #original: residuals = np.power(distances, 2) - np.power(radii, 2) 
    #residuals = residuals * ((radii - 3) / np.sqrt(1 + np.power(radii-3, 2)) + 2)

    return np.ravel(residuals)
