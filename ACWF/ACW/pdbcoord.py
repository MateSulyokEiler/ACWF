# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 11:10:59 2022

@author: raxis
"""

import pymol
from pymol import cmd
from pymol import stored
from pymol.wizard import Wizard
from itertools import combinations
import numpy as np
import math

from .basic import *


def pdb_reader(textin):
    #sheet from pymol reder

    with open(textin, 'r+') as f:
        model=mAtom()
        textinlist = f.read().splitlines()
        for line in textinlist:
            if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':

                name = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                b = float(line[60:66].strip())
                elem = line[76:78].strip()
                resn = line[17:20].strip()
                resi = line[22:26].strip()
                chain = line[21:22].strip()
                ID = line[6:11].strip()
                alt = line[16:17].strip()

                atom = sAtom(name, x, y, z, b, elem, resn, resi, chain, ID, alt)

                model.add(atom)

    return model


def coord_reader(coord,name='-'):
    #coord to atom
    output = sAtom(name,float(coord[0]),float(coord[1]),float(coord[2]))

    return output


def PCA(X , num_components):

    X_meaned = X - np.mean(X , axis = 0)

    cov_mat = np.cov(X_meaned , rowvar = False)

    eigen_values , eigen_vectors = np.linalg.eigh(cov_mat)

    sorted_index = np.argsort(eigen_values)[::-1]
    sorted_eigenvalue = eigen_values[sorted_index]
    sorted_eigenvectors = eigen_vectors[:,sorted_index]

    eigenvector_subset = sorted_eigenvectors[:,0:num_components]

    X_reduced = np.dot(eigenvector_subset.transpose() , X_meaned.transpose() ).transpose()

    return X_reduced, eigenvector_subset


def find_axis(model):
    #attempts to find main axis

    chain = model.lista[0].chain
    resi = model.lista[0].resi

    """
    main_axis = ['x','y','z']
    for axis in ['x','y','z']:

        coords = model.select('chain',chain).select('resi', ['1']).select('name',['CA']).data(axis)
        print(coords)
        #cheks if the coord is same in the next sheet
        for coord2 in coords[1:]:
            print(abs(coords[0]-coord2))
            if abs(coords[0]-coord2) < 0.05:
                main_axis.remove(axis)
                break

    if len(main_axis)==1:
        main_axis=main_axis[0]
        #print(main_axis)

    else:
        print('No main axis found!')

    """

    """
    #alternative axis search with tolerance
    diff = []
    for axis in ['x','y','z']:

        coords = model.select('chain',chain).select('resi', resi).select('name',['CA']).data(axis)
        #print(coords)
        diff.append(max(coords)-min(coords))
    #print(diff)
    main_axis=['x','y','z'][diff.index(max(diff))]

    """
    """
    #tolerance for antiparallel high symmetry
    #heigth needs to be around 4.8 A
    chain_count = len(model.splice(sort=False)) -1
    diff = []
    for axis in ['x','y','z']:

        coords = model.select('chain',chain).select('resi', resi).select('name',['CA']).data(axis)
        #print(coords)
        diff.append(abs((max(coords)-min(coords))/chain_count-4.8))

    main_axis=['x','y','z'][diff.index(min(diff))]
    """

    #from straigthen method:
    #In the sheet find vectors that are parallel
    #A = model.select('chain',chain).select('resi', resi).select('name', ['C']).lista[0].xyz()
    A = model.select('resi', resi).select('name', ['C']).lista[0].xyz()
    #print(len(model.select('resi', resi).select('name', ['C']).lista))
    vectors = []
    for atom in model.select('resi', resi).select('name', ['C']).lista[1:]:
        vec=[A[0] - atom.x, A[1] - atom.y, A[2] - atom.z]
        vectors.append(vec)
    #print(vectors)
    for comb in combinations(vectors,2):
        #print(get_lc(np.cross(comb[0],comb[1]),[0,0,0]))
        if get_lc(np.cross(comb[0],comb[1]), [0, 0, 0]) < 0.01:
            vector = comb[0][:]
            vector_p = comb[0][:]
            break
    else:
        main_axis = find_axis_PCA(model)
        return main_axis
        #print('No main axis found!')
        #raise IndexError
    #print(vector)
    vector_abs = [abs(i) for i in vector]
    main_axis = ['x', 'y', 'z'][vector_abs.index(max(vector_abs))]
    return main_axis


def find_axis_PCA(model):

    sheet = model.splice(sort=False)

    atoms = []
    for chain in sheet:
        atoms.append(chain.center())

    data , eigenvector_subset = PCA(atoms,1)

    vector = [x[0] for x in eigenvector_subset]
    vector_abs = [abs(i) for i in vector]
    main_axis = ['x', 'y', 'z'][vector_abs.index(max(vector_abs))]

    return main_axis


def find_axis_PCA_spliced(sheet):

    atoms = [sheet[0].center(), sheet[int(len(sheet)/2)].center(), sheet[-1].center()]

    data , eigenvector_subset = PCA(atoms,1)

    vector = [x[0] for x in eigenvector_subset]
    vector_abs = [abs(i) for i in vector]
    main_axis = ['x', 'y', 'z'][vector_abs.index(max(vector_abs))]

    return main_axis


def get_length(a1,a2,show=''):
    #distance between atoms

    d = ( (a1.x-a2.x)**2 + (a1.y-a2.y)**2 + (a1.z-a2.z)**2 )** 0.5

    if show:
        show_distance(a1,a2,show)

    return d


def get_length_coord(a1,a2,show=''):
    #distance between atoms (coordinates)

    d = ( (a1[0]-a2[0])**2 + (a1[1]-a2[1])**2 + (a1[2]-a2[2])**2 )** 0.5

    if show:
        show_distance(sAtom('atom1',a1[0],a1[1],a1[2]),sAtom('atom2',a2[0],a2[1],a2[2]),show)

    return d


def get_lc(a1,a2):
    #distance between 2D coordinates

    d = ( (a1[0]-a2[0])**2 + (a1[1]-a2[1])**2 )** 0.5

    return d


def get_mlength(m1,m2):
    #list of distances of two models pairwise

    d = []

    for n in range(m1.lista):
        d.append(get_length(m1.lista[n], m2.lista[n]))

    return d


def get_midpoint(a1,a2,name='',show=''):
    #get the midpoint between to atoms

    model = mAtom()
    model.add(a1)
    model.add(a2)
    coord = model.center()

    mp = sAtom(name,coord[0], coord[1], coord[2])

    if show:
        show_distance(a1,a2,show)

    return mp


def get_midpoint_coord(a1,a2):
    #get the midpoint between to coordinates

    mp = []
    for n in range(len(a1)):
        mp.append((a1[n]+a2[n])/2)

    return mp


def get_line(a1,a2,axis):
    #calculates the line through two atoms that are projected along the main axis
    # a*x + b*y + c = 0

    x1, y1 = a1.xyz(axis)
    x2, y2 = a2.xyz(axis)

    #          a      b         c
    line = [ y1-y2, x2-x1 , x1*y2-x2*y1 ]

    return line


def vector_segment(a1,a2,d):
    #along a vector(defined by two points) gives a point that is d distance away

    if isinstance(a1,sAtom):
        a1 = a1.xyz()
    if isinstance(a2,sAtom):
        a2 = a2.xyz()

    l = ( (a1[0]-a2[0])**2 + (a1[1]-a2[1])**2 + (a1[2]-a2[2])**2 )** 0.5

    coord = [a1[0] + (a2[0]-a1[0])*d/l, a1[1] + (a2[1]-a1[1])*d/l, a1[2] + (a2[2]-a1[2])*d/l]

    return coord


def get_distline(a1,line,axis):
    #calculates the distance from a projected line to a projected point

    x1, y1 = a1.xyz(axis)
    a,b,c = line

    d = abs(a*x1 + b*y1 + c)/(a*a + b*b)**0.5

    return d


def get_mdistline(m1, line, axis):
    #list of distances from a projected line to a projected point

    d = []

    for n,a1 in enumerate(m1.lista):
        d.append([get_distline(a1, line, axis), n, a1])

    d.sort(key=lambda close: close[0])

    return d


def get_perpendicular_bisector(a1, a2):
    #2D coordinated
    mid = get_midpoint_coord(a1, a2)
    slope = -1 / ((a2[1] - a1[1]) / (a2[0] - a1[0]))
    intercept = mid[1] - slope * mid[0]
    return slope, intercept


def get_intersection(a1, a2):
    #slope, intercept
    x = (a2[1] - a1[1]) / (a1[0] - a2[0])
    y = a1[0] * x + a1[1]
    return x, y


def triangulate(a1,a2,a3,axis):
    #input atoms, a1-2 is the line, a3 is the external point

    x1, y1 = a1.xyz(axis)
    x2, y2 = a2.xyz(axis)
    x3, y3 = a3.xyz(axis)

    v = (x2-x1,y2-y1)

    d = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)

    return v, d


def add_axis(point, axis, name):
    if axis == 'z':
        a1 = sAtom(name,point[0],point[1],0)

    elif axis == 'y':
        a1 = sAtom(name,point[0],0,point[1])

    else:
        a1 = sAtom(name,0,point[0],point[1])

    return a1


def flatten(model,axis):
    #flatten along axis

    output=mAtom()

    for atom in model.lista:
        output.add(atom.axis(axis))

    return output


def get_closestpoint(a1, line, axis, name):
    # calculates the closest point on the line 2D

    [x1, y1] = a1.xyz(axis)
    [a,b,c] = line

    cp = [(b*(b*x1-a*y1)-a*c)/(a*a+b*b), (a*(-b*x1+a*y1)-b*c)/(a*a+b*b)]

    cpa = add_axis(cp, axis, name)

    return cpa


def get_dist_3d_line(p, q, rs):
    # distance from a line, numpy input
    line_vector = q - p  # Vector along the line
    line_vector_norm = np.dot(line_vector, line_vector)  # Squared norm of the line vector
    t = np.dot(rs - p, line_vector) / line_vector_norm # Projection scalar for each point in rs

    closest_points = p + np.outer(t, line_vector)
    distances = np.linalg.norm(closest_points - rs, axis=1)

    return distances


def get_closest_point_3d(p, q, r):
    #closest point on p,q line from r, numpy input
    
    line_vector = q - p  # Vector along the line
    line_vector_norm = np.dot(line_vector, line_vector)  # Squared norm of the line vector
    t = np.dot(r - p, line_vector) / line_vector_norm # Projection scalar for each point in r

    closest_points = p + np.outer(t, line_vector)

    return closest_points


def get_nearest(atom1, model, fast = False):
    #list from closest to furthest
    #distance, Nth atom, atom

    nearest = []
    x,y,z = atom1.xyz()

    for n, atom2 in enumerate(model.lista):
        nearest.append([( (x-atom2.x)**2 + (y-atom2.y)**2 + (z-atom2.z)**2 )** 0.5,n,atom2])

    if fast:
        return min(nearest, key=lambda close: close[0])
    else:
        nearest.sort(key=lambda close: close[0])
        return nearest


def get_mnearest(model1, model2, fast = False):

    nearest = []
    for n1, atom1 in enumerate(model1.lista):
        x,y,z = atom1.xyz()
        for n2, atom2 in enumerate(model2.lista):
            nearest.append([( (x-atom2.x)**2 + (y-atom2.y)**2 + (z-atom2.z)**2 )** 0.5,n1,atom1,n2,atom2])

    if fast:
        return min(nearest, key=lambda close: close[0])
    else:
        nearest.sort(key=lambda close: close[0])
        return nearest


def get_dihedral(*args):
    
    a0, a1, a2, a3 = map(np.array,args)
    
    
    b0 = -1.0*(a1 - a0) #vectors between points
    b1 = a2 - a1
    b2 = a3 - a2

    b1 =  b1 / np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.arctan2(y, x)


def is_same(atom1,atom2):
    #compares the coordinates of two atoms

    for coord in ['x','y','z']:
        if getattr(atom1, coord) != getattr(atom2, coord):
            return False

    return True


def get_zshift(model1,model2):
    #calculates Z shift from backbone atoms

    axis = find_axis(model1)

    sheet1 = model1.splice()
    sheet2 = model2.splice()
    sheet1 = sorted([sheet1[i].select('name',['CA']).average(axis) for i in range(len(sheet1))])
    sheet2 = sorted([sheet2[i].select('name',['CA']).average(axis) for i in range(len(sheet2))])

    middle = sheet1[int(len(sheet1)/2)]
    #calculates the Z shift for the middle chain to the closest
    sheet2 = sorted([abs(sheet2[i]-middle) for i in range(len(sheet2))])

    return (sheet2[0])


def get_CAdistance(model1, model2, show):
    #avaerage distance of the nearest CA

    d = 0
    connections = []
    for atom1 in model1.lista:
        connections.append(get_nearest(atom1,model2)[0])
        if show:
            show_distance(atom1,connections[-1][2])
    for nearest_atom in connections:
        d += nearest_atom[0]
    d = d/len(connections)

    return d, connections


def get_sheetdistance(modelA, modelB):
    #calculates CA distances between the middle elements of two sheets

    axis = find_axis(modelA)
    model1 = modelA.splice()[int(len(modelA.splice())/2)]
    model2 = sorted(modelB.splice(), key=lambda model: abs(model.backbone().average(axis) - model1.backbone().average(axis)))[0]

    d1A2A = get_CAdistance( model1.select('name',['CA']), model2.select('name',['CA']), True)[0]
    d2A1A = get_CAdistance( model2.select('name',['CA']), model1.select('name',['CA']), True)[0]
    d = (d1A2A+d2A1A)/2

    d1A2 = get_CAdistance( model1.select('name',['CA']), model2, False)[0]
    d2A1 = get_CAdistance( model2.select('name',['CA']), model1, False)[0]
    d_overlap = d1A2A - d1A2 + d2A1A - d2A1 - d

    d_lista = [d1A2A, d1A2, d2A1A, d2A1]

    return d, d_overlap, d_lista


def show_distance(atom1, atom2, name='dist'):

    cmd.pseudoatom('tmp1', pos=atom1.xyz())
    cmd.pseudoatom('tmp2', pos=atom2.xyz())
    cmd.distance(name, 'tmp1', 'tmp2')
    cmd.delete('tmp1')
    cmd.delete('tmp2')

    return


def merge_spliced(models):
    
    output = mAtom()
    
    for model in models:
        output.addmodel(model)
    
    return output


#------------------------------------

#multi atoms
class mAtom:

    def __init__(self):
        self.lista = []

    def add(self,atom):

        self.lista.append(atom)

    def addmodel(self, model):

        for atom in model.lista:
            self.lista.append(atom)

    def select(self, attr, select):
        #creates a selection within a model
        #accepts a list for select

        output = mAtom()

        for atom in self.lista:
            if getattr(atom, attr) in select:
                output.add(atom)

        return output

    def backbone(self, backbone=True):

        output = mAtom()

        for atom in self.lista:
            if getattr(atom, 'name') in ['C','O','N','CA','OXT', 'OX']: #OX becouse areaimol output is bugged
                if backbone:
                    output.add(atom)
            else:
                if not backbone:
                    output.add(atom)

        return output

    def sidechain(self, sidechain=True):
        #CA inclusive for glycine
        output = mAtom()

        for atom in self.lista:
            if getattr(atom, 'name') not in ['C','O','N','OXT', 'OX']: #OX becouse areaimol output is bugged
                if sidechain:
                    output.add(atom)
            else:
                if not sidechain:
                    output.add(atom)

        return output

    def data(self, attr):
        #expose a single attr for all atoms

        output = []

        for atom in self.lista:
            output.append(getattr(atom, attr))

        return output

    def sort(self, attr):
        #sorts itself according to attr

        self.lista = sorted(self.lista, key=lambda atom: getattr(atom, attr))
        return self

    def sort_float(self, attr):
        #sorts itself according to attr converted to float

        self.lista = sorted(self.lista, key=lambda atom: float(getattr(atom, attr)))
        return self

    def sort_coord(self):
        #sorts itself according to all coordinates

        #self.lista = sorted(self.lista, key=lambda atom: 100000001000000010000000 + float(atom.x)*10000000000000000000 + float(atom.y)*100000000000 + float(atom.z)*1000)
        self.lista = sorted(self.lista, key=lambda atom: (float(atom.x), float(atom.y), float(atom.z)))
        return self


    def prnt(self):

        for n, atom in enumerate(self.lista):
            print(n)
            atom.prnt()
            print()

    def average(self, attr):
        #average of all coord or b

        output = 0

        for atom in self.lista:
            output += getattr(atom, attr)

        output = output / len(self.lista)

        return output

    
    def median(self, attr):
        #median of all coord or b

        output = []

        for atom in self.lista:
            output.append(getattr(atom, attr))

        output = output[int(len(output)/2)]

        return output


    def splice(self, sort=True):
        #makes a list of chain-models of a sheet, only if ordered
        #detects the restart of numbering
        #try:
        output = []
        model = mAtom()
        res = int(self.lista[0].resi)

        for n, atom in enumerate(self.lista):
            if int(atom.resi)<res:
                output.append(model)
                model = mAtom()

            model.lista.append(atom)
            res = int(atom.resi)

        output.append(model)

        if sort:
            axis = find_axis(self)
            output = sorted(output, key=lambda i : float(i.median(axis)) )

        return output
        #except:
        #    return []

    def splice_mid(self, sort=True):
        #splice than return the middle sheets as one model

        output = self.splice(sort)

        output2 = mAtom()
        for chain in output[int(len(output)/3):int(2*len(output)/3)]:
            output2.addmodel(chain)

        return output2


    def splice_same(self, sort=True):
        #splice than return every nth sheets as one model
        #in crystal data those should be the same

        output = self.splice(sort)

        output2 = mAtom()
        for chain in output[::int(len(output)/3)]:
            output2.addmodel(chain)

        return output2


    def splice_equal(self, parts, sort=True):
        #splice into equal parts, only if each chain has the same length

        output = []
        length = int(len(self.lista)/parts) #chain
        
        if len(self.lista) % parts:
            raise ValueError
        
        for n in range(parts):
            chain = mAtom()
            for atom in self.lista[n*length:n*length+length]:
                chain.add(atom)
            output.append(chain)

        if sort:
            axis = find_axis_PCA_spliced(output)
            output = sorted(output, key=lambda i : float(i.median(axis)) )

        return output


    def residues(self):
        #makes a list of residue-grouped models, only if ordered
        #detects the jump of numbering

        output = []
        model = mAtom()
        res = int(self.lista[0].resi)
        chain = self.lista[0].chain
        
        for n, atom in enumerate(self.lista):

            if int(atom.resi)!=res or atom.chain!=chain:
                output.append(model)
                model = mAtom()

            model.add(atom)

            res = int(atom.resi)
            chain = atom.chain

        output.append(model)

        return output


    def center(self):
        #returns the center of mass

        center_point = [0, 0, 0]

        for atom in self.lista:
            center_point=[center_point[i]+atom.xyz()[i] for i in range(3)]

        center_point=[center_point[i]/len(self.lista) for i in range(3)]

        return center_point


    def shift(self, vector):

        for atom in self.lista:
            atom.vector(vector)

        return self


    def copy(self):

        model = mAtom()
        model.addmodel(self)
        return model
    
    def coord_list(self):
        
        coords = []
        for atom in self.lista:
            coords.append(atom.xyz())
        
        return coords
    


#single atom
class sAtom:

    def __init__(self, name, x, y, z, b=0, elem='', resn='', resi='', chain='', ID='', alt=''):
        #necessery
        self.name = name #N, CA, OXT ect
        self.x = x
        self.y = y
        self.z = z
        #optional
        self.b = b
        self.elem = elem
        self.resn = resn #name amino acid 3 char
        self.resi = resi #number
        self.chain = chain
        self.ID = ID #pdb index
        self.alt = alt

    def xyz(self, axis=''):
        output=[]
        output.append(self.x)
        output.append(self.y)
        output.append(self.z)

        if axis:
            project = {'x':0,'y':1,'z':2}[axis]
            output.pop(project)

        return output


    def xyz_key(self):

        return tuple(round(c, 2) for c in self.xyz())


    def axis(self, axis):
        #projection along axis

        if axis == 'z':
            self.z = 0

        elif axis == 'y':
            self.y = 0

        elif axis == 'x':
            self.x = 0

        return self


    def vector(self, vector):
        #translate by vector

        self.x += vector[0]
        self.y += vector[1]
        self.z += vector[2]

        return self


    def prnt(self):

        for attr in vars(self):
            print(getattr(self, attr)),
            print(' '),
