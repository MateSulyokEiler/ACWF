# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 18:15:56 2021
Last edited on Mon Dec 04 12:14:30 2023
Final pre-release version
@author: mate
"""

import traceback
try:
    from importlib import reload
except:
    pass


import os
import sys
import glob
import shutil
import threading
import time
import inspect
import platform
import json
import numpy as np
from pymol import cmd
from pymol.wizard import Wizard
import math
import pymol
from pymol.cgo import *
import random

# adding install folder to the system path
acw_database_path = os.path.dirname(inspect.getfile(inspect.currentframe()))
if acw_database_path == '':
    acw_database_path = os.getcwd()
sys.path.insert(0, acw_database_path)

from ACW import *
from ACW.acw_method_fulllength import *


def suspend_decorator(func):
    def inner(self, *args, **kwargs):

        cmd.set('suspend_updates', 1, quiet=1)
        try:
            returned_value = func(self, *args, **kwargs)
        finally:
            cmd.set('suspend_updates', 0, quiet=1)
        return returned_value

    return inner


def busy_decorator(func):
    def inner(self, *args, **kwargs):

        if self.busy:
            print('Wait for completion')
            return None
        try:
            self.busy = 1
            returned_value = func(self, *args, **kwargs)
        finally:
            self.busy = 0
        return returned_value

    return inner


class AcWizardFibrilMultiSheet(Wizard):
    """Pymol wizard for the evaluation and display of amyloid fibril structures.

    The wizard offers two mode:
    Evaluation mode to analyze a new amyloid fibril structure.
    Database mode to display PDB structures.
    """
    def __init__(self, location):
        Wizard.__init__(self)

        self.main_path = location
        try:

            cmd.cd(self.main_path)
        except Exception:
            print('Installation folder is missing')
            raise OSError

        self.advanced = False
        if self.advanced:
            self.pdb_list = read_save('amy_flength_saveout_35.json')

            #self.pdb_list = read_save('amy_flength_saveout_10_bad2.json')
            #self.pdb_list_old = read_save('amy_flength_saveout_25.json')
            #self.pdb_list = read_save('final_pdb_search10.json')
            print(len(self.pdb_list.lista))
        else:
            try:
                self.pdb_list = read_save('acw_personal_database.json')
                print('Personal ACW database is now loaded')
            except Exception:
                try:
                    self.pdb_list = read_save('acwf_database.json')
                    print('Original ACW database is now loaded')
                except Exception:
                    print('Database files are missing or corrupted')
                    self.pdb_list = amy_struct_list()

        self.path = None
        config_reader(self)
        if self.path:
            if not glob.glob('%s/bin/areaimol*'%self.path):
                self.path = None
                print('Areaimol not found at CCP4 path, the evaulation mode is not available')
            elif not glob.glob('%s/bin/sc*'%self.path):
                self.path = None
                print('Sc not found at CCP4 path, the evaulation mode is not available')
            if self.path:
                print('CCP4 path found, evaulation mode is available')
        else:
            print('CCP4 not detected, the evaulation mode is not available')

        # Pymol settings
        cmd.set('mouse_selection_mode', 1) # set selection mode to residue
        # 2=chain, 4=object
        cmd.set('orthoscopic')
        cmd.set('label_size', '20')
        cmd.set('dash_color', 'brightorange')
        try:
            cmd.undo_disable()
        except Exception:
            cmd.set('suspend_undo', 1)

        # selection criteria possible values
        self.sele_prot_types = list(map(str, sorted(list(set(self.pdb_list.data('prot_abr'))))))
        self.sele_protofilament_types = list(map(str, sorted(list(set(self.pdb_list.data('protofilament'))))))
        self.sele_mainclust_types = ['1', '2', '3']
        self.sele_secclust_types = ['1', '2', '3', '4', '5', '6', '7', '8', '9']
        # selection criteria selection containers
        self.sele_prot = []
        self.sele_protofilament = []
        self.sele_mainclust = []
        self.sele_secclust = []

        # inital wizard parameters
        self.picking_mode_types = ['reset_chains', 'connected', 'Delete', 'OFF']
        self.picking_mode = 'OFF'
        self.need_prompt = 'ON'
        self.generate_chains = 'OFF'
        self.view_count = 0
        self.seq = ''
        self.seq_len = 0
        self.main_axis = 0

        self.colors = ['red', 'blue', 'orange', 'yellow', 'cyan', 'magenta', 'lightblue', 'brown', 'purple'] #for picking
        self.color_range = 'absolute'
        self.color_inter = 'inter'
        self.color_last_method = 'sheets'

        self.manual_connected_sections = []
        self.axis = 'OFF'
        self.fibril_axis = 'OFF'
        self.error_message = ''
        self.busy = 0
        self.pick_count = 0 # number of sheet picked
        self.pick_count_cs = 0 # number of sheet picked in connected section
        self.sheet_count = -1 # number of chain in sheet
        self.pdb_count = 0
        self.cluster_count = 0 #where to start for clusters

        self.window_current = 0

        self.labels = 'OFF'
        self.interface_count = 0 # which interface to display
        self.ext_count = 0 # counter for extended sheets advanced
        self.n_Msheet = 0 # sheet counters for automatic
        self.n_Bsheet = 0
        self.favor_param = None
        self.overwrite_state = 0

        # set up database state
        self.mode = 'database'
        self.calc_state = [0]
        self.main_object = None
        self.reset_counters()
        try:
            self.next_peptide(x=0)
        except IndexError:
            self.mode_change("evaluate")
        self.menu_update()
        cmd.refresh_wizard()


    def mode_change(self, new_mode):
        """Changes between the two main modes."""

        if new_mode == 'database':
            if self.busy:
                print('Wait for completion')
                self.error_message = 'Wait for completion'
                return
            if len(self.pdb_list.lista) == 1:
                print('Database is missing, database mode is not available')
                self.error_message = 'Database is missing, database mode is not available'
                return
            self.picking_mode = 'OFF'
            cmd.delete('all')
            self.delete_all(['everything'])
            self.reset_counters()
            self.mode = new_mode
            try:
                cmd.cd(self.main_path)
                cmd.cd('%s'%self.pdb_list[0].code)
            except Exception:
                print('Database is missing, database mode is not available')
                self.error_message = 'Database is missing, database mode is not available'
                self.menu_update()
                self.mode = 'evaluate'
                return
            self.pdb_count = 0
            self.pdb_list.lista.pop(-1)
            self.next_peptide(x=0) #50
            self.menu_update()

        elif new_mode == 'evaluate':
            if not self.path:
                print('Without CCP4, evaluation mode is not available')
                self.error_message = 'Without CCP4, evaluation mode is not available'
                return

            self.picking_mode = 'OFF'
            cmd.delete('all')
            self.delete_all(['everything'])
            cmd.cd(self.main_path)
            self.reset_counters()
            self.color_last_method = 'sheets'
            self.pdb_list.add_blank()
            self.pdb_count = -1
            self.calc_state = [0]
            self.mode = new_mode
            self.menu_update()


    def next_peptide_cluster(self, n=1):
        """Loads in next fibril in the order of clustering"""

        if not self.sele_prot or len(self.sele_prot) > 1:
            print('Select exactly one protein type.')
            self.error_message = 'Select exactly one protein type.'
            return

        if self.sele_prot[0] in ('misc.', 'SA1', 'OTHER'):
            print('This type has not enough elements to clusterize.')
            self.error_message = 'This type has not enough elements to clusterize.'
            return


        self.sele_protofilament = []
        with open('../%s_nofit_clusters_order_overlap.json'%self.sele_prot[0], "r") as f:
            sort_index = json.load(f)

        if n == -1:
            self.cluster_count -= 2

        if self.cluster_count < 0:
            self.cluster_count = len(sort_index) -1

        if self.cluster_count >= len(sort_index):
            self.cluster_count = 0

        print(self.cluster_count+1, sort_index[self.cluster_count])
        self.next_peptide(goto=sort_index[self.cluster_count][0:4])

        self.cluster_count += 1


    @busy_decorator
    @suspend_decorator
    def next_peptide(self, n=1, x=None, goto=None):
        """Cleans up and loads in next fibril from the selection"""

        if self.pdb_count >= len(self.pdb_list.lista):
            self.pdb_count = -1

        if self.mode == 'evaluate':
            pass
            #cmd.get_wizard().mode_change("database")

        cmd.delete('all')
        self.delete_all(['everything'])
        self.reset_counters()
        self.picking_mode = 'OFF'
        cmd.set('mouse_selection_mode', 1) # set selection mode to residue


        cmd.cd(self.main_path)
        #write_save('acw_personal_database.json', self.pdb_list)
        if self.advanced:
            #write_save('amy_struct_saveout.json', self.pdb_list)
            #shutil.copy(self.main_path+'/acw_personal_database.json', self.main_path+'/%s'%self.cur_pdb().code)
            pass
        cmd.cd('%s'%self.cur_pdb().code)

        while True:
            if x is not None:
                self.pdb_count = x
                break
            if goto is not None:
                if isinstance(goto, int):
                    print('To go to a given structure write its pdb code or its number in the database')
                    self.error_message = 'To go to a given structure write its pdb code or its number in the database'
                    self.pdb_count = 0
                else:
                    for pos, amy in enumerate(self.pdb_list.lista):
                        if amy.code == goto:
                            self.pdb_count = pos
                            break
                    else:
                        try:
                            goto = int(goto)-1
                            if 0 <= goto < len(self.pdb_list.lista):
                                self.pdb_count = goto
                            else:
                                print('Given number or PDB-code not in range, going to first structure')
                                self.error_message = 'Given number or PDB-code not in range, going to first structure'
                                self.pdb_count = 0
                        except ValueError:
                            print('File name not found, going to first')
                            self.error_message = 'File name not found, going to first'
                            self.pdb_count = 0
                    break

            # check next structure
            self.pdb_count += n

            if self.pdb_count == -2:
                self.pdb_count = len(self.pdb_list.lista)-1
            #print(self.pdb_count)
            if self.pdb_count >= len(self.pdb_list.lista) or self.pdb_count == -1:
                self.pick_count = -2
                self.pdb_count = -1
                cmd.cd('../%s'%self.pdb_list[0].code)
                print('No more structures.\nPress next or previous peptide to continue')
                cmd.refresh_wizard()
                return

            if self.advanced:
                print(self.cur_pdb().code)
                #find errors
                """
                if self.cur_pdb().rmsd > 0.5:
                    break
                else:
                    continue
                """

            # criteria to open:
            selection = True
            if self.sele_prot:
                if self.cur_pdb().prot_abr not in self.sele_prot:
                    selection = False
            if self.sele_protofilament:
                if self.cur_pdb().protofilament not in self.sele_protofilament:
                    selection = False
            if self.sele_mainclust and self.sele_secclust:

                for cluster in self.cur_pdb().clusters:
                    if cluster in [[int(x), int(y)] for x in self.sele_mainclust for y in self.sele_secclust]:
                        break
                else:
                    selection = False

            if selection:
                break


        # change directory
        cmd.cd('..')
        if not os.path.isdir(self.cur_pdb().code):
            os.makedirs(self.cur_pdb().code)
        cmd.cd(self.cur_pdb().code)

        try:
            #for analyzed pdb
            self.setup()

            # new files
            #self.read_in() # first time

            #self.test2()
            # self.picking_mode = 'connected'
            #self.picking_mode = 'reset_chains'

            self.menu_update()

            # if len(self.cur_pdb().interface) > 0:
            #     self.next_interface()

            # self.toggle_labels()

        except Exception:
            #traceback.print_exc()
            print('Entry is missing files or corrupted')
            self.error_message = 'Entry is missing files or corrupted'
            return

        finally:
            self.reset()


    def next_interface(self):
        """Unused"""

        self.reset()


    def read_in(self):
        """Open fibril first time."""

        self.main_object = self.cur_pdb().code

        # load
        cmd.load('%s.pdb'%self.cur_pdb().code)

        cmd.select('sele', 'none')
        self.view_count = 2
        self.view()

        #read in protein name
        with open('%s.pdb'%self.cur_pdb().code, 'r+') as f:
            for line in f.read().splitlines():
                if line[0:21] == 'COMPND   2 MOLECULE: ':
                    self.cur_pdb().prot = line.strip()[21:-1]
                    print(self.cur_pdb().prot)
                    break

        # find all sheets in assymetric unit
        cmd.select('sele_main', '%s and polymer.protein'%self.main_object, 0, 1, 0, 1)

        x = 0
        while len(get_data('sele_main', 'chain')) > 0:
            x += 1

            try:
                cmd.select('sele', 'bychain first sele_main', 0, 1, 0, 1)
                atom_name = self.main_object + str(self.pick_count)
                self.pick_sheet(atom_name)
            except IndexError:
                cmd.select('sele_main', 'sele_main and not sele_%s'%(self.pick_count-1), 0, 1, 0, 1)
                cmd.delete('sele_%s'%self.pick_count)
                continue

            cmd.select('sele_main', 'sele_main and not sele_%s'%(self.pick_count-1), 0, 1, 0, 1)

            if x == 500:
                print('Error in locating sheets')
                cmd.refresh_wizard()
                return

        cmd.delete('sele_main')
        self.cur_pdb().protofilament = str(self.pick_count)
        print(self.cur_pdb().protofilament)

        cmd.refresh_wizard()


    def setup(self):
        """Setups and finds the sheets from a open structure."""

        if self.mode == 'evaluate' and self.calc_state[0] != 1:
            print('Please use the methods in order.')
            self.error_message = 'Please use the methods in order.'
            return

        if self.mode == 'evaluate':
            self.calc_state = ['1a']

        # doesnt create new folder!
        #cmd.load('%s.pdb'%self.cur_pdb().code)
        
        if self.mode == 'evaluate':
            self.picking_mode = 'OFF'
            cmd.delete('all')
            self.delete_all(['everything'])
            
        for n in range(int(self.cur_pdb().protofilament)): #avoids states
            cmd.load('monomer_%s.pdb'%n)
            cmd.copy_to(self.cur_pdb().code, 'monomer_%s'%n)
            cmd.delete('monomer_%s'%n)

        self.main_object = cmd.get_names()[0]

        # setup the sheets
        # remove stuff that Sc doesnt like
        # remove water, salt, etc.
        # cmd.remove('not polymer.protein')
        # remove hydrogens
        cmd.remove('elem H')
        # select & remove all non A altlocs
        cmd.remove(' not (alt ''+A)')
        # remove the alt identifier
        cmd.alter('all', "alt=''")

        """
        # renumber chains if needed
        all_chains = list(dict.fromkeys(get_data('polymer.protein', 'chain')))
        self.cur_pdb().first_res = str(get_data('first (chain %s and polymer.protein)'%all_chains[0], 'resi')[0])
        for chain in all_chains:
            offset = int(get_data('first (chain %s and polymer.protein)'%chain, 'resi')[0]) -1

            if offset==0:
                continue
            else:
                cmd.alter('chain %s and polymer.protein'%chain, "resi=str(int(resi)-%s)" % str(offset))
        """

        #check for uneven sidechain
        
        first_chains = pdb_reader('monomer_0.pdb').splice()
        for n in range(int(self.cur_pdb().protofilament)):
            chains = pdb_reader('monomer_%s.pdb'%n).splice()
            if len(chains) <3:
                print('Not enough chains in the fibril.')
                self.error_message = 'Not enough chains in the fibril.'
                if self.mode == 'evaluate':
                    self.calc_state = [5]
                    return
            if len(chains) != len(first_chains):
                print('The number of chains is different in the protofilaments.')
                self.error_message = 'The number of chains is different in the protofilaments.'
                if self.mode == 'evaluate':
                    self.calc_state = [5]
                    return
            select_model_chain([chains[0]], 'sele')
            chain_len = len(get_data('(bychain sele) and name CA ', 'resv'))
            for chain in chains:
                select_model_chain([chain], 'sele')
                #print(len(get_data('(bychain sele) and name CA ', 'resv')))
                if len(get_data('(bychain sele) and name CA ', 'resv')) != chain_len:
                    print('The chains contain a different number of residues.')
                    self.error_message = 'The chains contain a different number of residues.'
                    self.cur_pdb().uneven = 'yes'
                    if self.mode == 'evaluate':
                        self.calc_state = [5]
                        return

        # setup sheets
        cmd.color('green', 'elem C')
        cmd.orient('all')
        cmd.show('wire')

        # color and sheets
        color_methods(self, 'sheets') # creates sele_X!
        cmd.select('sele', 'none')
        #color_methods(self, 'recolor-sc-6')

        #count resvrange
        chains = pdb_reader("monomer_0.pdb").splice()
        self.sheet_count = len(chains)
        x, y, z = chains[0].lista[0].xyz()
        cmd.select('sele', 'bychain sele or ( X>%s and X<%s and Y>%s and Y<%s and Z>%s and Z<%s)'%(x-0.005, x+0.005, y-0.005, y+0.005, z-0.005, z+0.005), 0, 1, 0, 1)
        self.resvrange = get_data('(bychain sele) and name O ', 'resv')
        if self.mode == 'evaluate':
            new_entry_analyser(self)

        """
        #concentric cyclinderlinder_calc
        p, cost, status = lstsq_cylinder_calc_axis(self)
        p_dict= {}
        p_dict['cost'] = cost
        p_dict['status'] = status
        p_dict['rmax'] = max(p[6:])
        p_dict['start'] = p[0:3]
        p_dict['end'] = p[3:6]
        p_dict['radiuses'] = p[6:]
        self.cur_pdb().p_cylinder = {'original_axis':p_dict}          #this resets the save!!!
        """

        #self.test()
        #self.test0()

        #generate more fibrils

        if self.generate_chains == 'ON' or self.mode == 'evaluate':
            rmsd, to_generate = self.extend_filaments()
            self.cur_pdb().rmsd = round(float(rmsd), 3)
        
        if self.mode == 'evaluate':
            self.cur_pdb().to_generate = to_generate
        #self.test()

        """
        p, cost, status = lstsq_cylinder_calc_axis(self)
        p_dict= {}
        p_dict['cost'] = cost
        p_dict['status'] = status
        p_dict['rmax'] = max(p[6:])
        p_dict['start'] = p[0:3]
        p_dict['end'] = p[3:6]
        p_dict['radiuses'] = p[6:]
        self.cur_pdb().p_cylinder['generated_axis'] = p_dict
        """

        if self.mode == 'evaluate':
            self.calc_state = [2]

        #recolor
        if self.color_last_method != 'sheets':
            self.color(method=self.color_last_method)

        self.pick_count = 0
        cmd.select('sele', 'none')
        cmd.unpick()
        cmd.deselect()
        cmd.refresh_wizard()


    def reset(self):
        """Reset everything picking related."""
        cmd.delete('sele')
        cmd.select('none')
        cmd.unpick()
        cmd.refresh_wizard()


    def suspend_refresh(self, wait=0.05):
        """Wait to allow pymol to draw a frame."""
        cmd.set('suspend_updates', 0, quiet=1)
        time.sleep(wait)
        cmd.set('suspend_updates', 1, quiet=1)


    def reset_counters(self):
        """Resets all state and data variables"""

        self.sheet_count = -1 # number of chain in sheet
        self.pick_count = 0 # number of sheet picked
        self.pick_count_cs = 0 # number of sheet picked in connected section
        self.interface_count = 0 # which interface to display
        self.window_current = 0

        self.ext_count = 0 # counter for extended sheets advanced
        self.n_Msheet = 0 # sheet counters for automatic
        self.n_Bsheet = 0
        self.favor_param = None
        self.manual_connected_sections = []

        self.view_count = 0
        self.labels = 'OFF'
        self.axis = 'OFF'
        self.fibril_axis = 'OFF'
        self.seq = ''
        self.main_axis = 0
        self.overwrite_state = 0
        cmd.refresh_wizard()


    @busy_decorator
    def reset_monomers(self):
        
        if self.mode == 'evaluate':
            try:
                self.main_object = cmd.get_names()[0]
            except IndexError:
                print('No open structure found')
                self.error_message = 'No open structure found'
                return
            self.cur_pdb().code = self.main_object

        cmd.deselect()
        cmd.delete('sele*')
        cmd.delete('dist*')
        cmd.delete('obj*')
        cmd.color('green', 'elem C')
        self.picking_mode = 'reset_chains'
        self.pick_count = 0


    def reset_selection(self):
        """Resets selection, also goes to first structure"""
        self.sele_prot = []
        self.sele_protofilament = []
        self.sele_mainclust = []
        self.sele_secclust = []
        self.cluster_count = 0
        self.next_peptide(x=0)
        self.menu_update()
        cmd.refresh_wizard()


    def delete_all(self, to_rm=[]):
        """Delete sheets, files and interfaces"""

        cmd.deselect()
        cmd.delete('sele*')
        cmd.delete('dist*')
        self.pick_count = 0

        if 'everything' in to_rm:
            cmd.delete('all')
            self.reset_counters()

        if 'everything' in to_rm or 'files' in to_rm:
            for name in ('distance*', 'sc*', 'areaimol*', 'Sheet_*', 'sele*'):
                for file in glob.glob(name):
                    os.remove(file)
            if self.main_object:
                for file in glob.glob(self.main_object + '_*.pdb'):
                    os.remove(file)

        if 'interface' in to_rm:
            cmd.util.cbag('all')
            cmd.delete('lengths')
            self.cur_pdb().interface = []
            self.interface_count = 0

        if 'last_amy' in to_rm:
            self.pdb_list.lista.pop(-1)
            self.pdb_list.add_blank()
            self.color_last_method = 'sheets'

        if 'sheets' in to_rm: # unused, deletes sheets and folder
            for file in glob.glob('Bsheet_*'):
                os.remove(file)

        if 'everything' in to_rm:
            self.calc_state = [0]
            self.main_object = None

        cmd.refresh_wizard()


    def cleanup(self):
        """Clean up, executed on quitting"""
        self.delete_all(['everything'])
        self.reset()
        cmd.delete('sele*')

        if self.mode == 'database':
            cmd.cd(self.main_path)
        elif self.mode == 'evaluate':
            self.pdb_list.lista.pop(-1)

        if self.advanced:
            write_save('amy_struct_saveout.json', self.pdb_list)
        if len(self.pdb_list.lista) > 0:
            #write_save('acw_personal_database.json', self.pdb_list)
            pass

        print('Cleanup finished')


    def get_prompt(self):
        """Text to be displayed at top left of display"""
        self.prompt = []

        if self.need_prompt == 'OFF':
            self.prompt = None
        
        else:
            if self.mode == 'database':
                if self.pick_count == -2:
                    self.prompt = ['No more structures.', 'Press next or previous structure to continue']
                else:
                    self.prompt = ['Database: %s %s/%s'%(self.cur_pdb().code, str(self.pdb_count+1), str(len(self.pdb_list.lista)))]
                    self.prompt.append(self.cur_pdb().prot)
                    try:
                        self.prompt.append('%s protofilament(s), %s rmsd'%(self.cur_pdb().protofilament, round(float(self.cur_pdb().rmsd), 3)))
                        if self.advanced:
                            self.prompt.append(''.join([' cluster: %s/%s'%(x[0], x[1]) for x in self.cur_pdb().clusters]))
                        else:
                            self.prompt.append('cluster(s): '+''.join(['%s/%s '%(n+1, x[0]) for n, x in enumerate(self.cur_pdb().clusters)]))
                    except AttributeError:
                        print('missing values')

            elif self.mode == 'evaluate':
                self.prompt.append('Evaluation')
                state = {
                    0:'Open structure and assign protofilaments',
                    1:'If all protofilaments are found, click on setup fibrils',
                    2:'Fibril has been setup, click on assign local interfaces',
                    3:'Local interfaces found, click on sliding window, to analyze them',
                    4:'Finished',
                    5:'Press reset to continue',
                    '1a':''}
                self.prompt.append(state[self.calc_state[0]])

            if self.error_message:
                self.prompt.append(self.error_message)

        return self.prompt


    def toggle_labels_residue(self):
        """Shows and hides the sheet labels"""

        if self.labels == 'residue':
            cmd.delete('Labels')
            self.labels = 'OFF'
        else:
            cmd.delete('Labels')
            try:
                for n in range(int(self.cur_pdb().protofilament)):
                    cmd.save('Sheet_%s.pdb'%n, 'sele_%s'%n)
                    protofilament = pdb_reader('Sheet_%s.pdb'%n).splice()

                    for x, name in enumerate(self.resvrange):
                        # xyz = protofilament[-1].residues()[x].center()
                        xyz = protofilament[-1].select('resi', [str(name)]).center()
                        name = str(name)
                        cmd.pseudoatom('Labels', pos=xyz, name=name, chain=ch_n(n))
                    cmd.label('Labels', 'name')
                    cmd.hide('wire', 'Labels')

            except Exception:
                #traceback.print_exc()
                print('No sheets were detected')
                self.error_message = 'No sheets were detected'

            self.labels = 'residue'


    def toggle_labels_chain(self):
        """Shows and hides the sheet labels"""

        if self.labels == 'chain':
            cmd.delete('Labels')
            self.labels = 'OFF'
        else:
            cmd.delete('Labels')
            try:
                for n in range(int(self.cur_pdb().protofilament)):
                    cmd.save('Sheet_%s.pdb'%n, 'sele_%s'%n)
                    protofilament = pdb_reader('Sheet_%s.pdb'%n).splice()

                    for x, chain in enumerate(protofilament):
                        # xyz = protofilament[-1].residues()[x].center()
                        xyz = chain.lista[0].xyz()
                        chain = chain.lista[0].chain
                        cmd.pseudoatom('Labels', pos=xyz, name=ch_n(n)+str(x), chain=chain)
                    cmd.label('Labels', 'chain')
                    cmd.hide('wire', 'Labels')

            except Exception:
                #traceback.print_exc()
                print('No sheets were detected')
                self.error_message = 'No sheets were detected'

            self.labels = 'chain'


    def toggle_labels_descriptor(self):
        """Shows and hides the sheet labels"""

        if self.labels == 'descriptor':
            cmd.delete('Labels')
            self.labels = 'OFF'
        else:
            cmd.delete('Labels')
            try:
                for n in range(int(self.cur_pdb().protofilament)):
                    cmd.save('Sheet_%s.pdb'%n, 'sele_%s'%n)
                    protofilament = pdb_reader('Sheet_%s.pdb'%n).splice()
                    self.resvrange = get_data('(bychain sele_%s) and name O ' % n, 'resv')
                    for x, name in enumerate(self.resvrange):
                        # xyz = protofilament[-1].residues()[x].center()
                        xyz = protofilament[-1].select('resi', [str(name)]).center()
                        b = protofilament[-1].select('resi', [str(name)]).lista[0].b
                        if b == -1:
                            continue
                        name = str(name)
                        cmd.pseudoatom('Labels', pos=xyz, name=name, chain=ch_n(n), b=b)
                    cmd.label('Labels', 'b')
                    cmd.hide('wire', 'Labels')

            except Exception:
                #traceback.print_exc()
                print('No sheets were detected')
                self.error_message = 'No sheets were detected'

            self.labels = 'descriptor'


    def toggle_fibril_axis(self):
        """Shows and hides the fibril axis"""

        if self.fibril_axis == 'ON':
            cmd.delete('Fibril_axis')
            self.fibril_axis = 'OFF'
        else:
            cmd.delete('Fibril_axis')
            try:

                cmd.save('Sheet_0.pdb', 'sele_0')
                p_filament = pdb_reader("Sheet_0.pdb").splice()
                diff = []
                for axis in ['x', 'y', 'z']:
                    coords = []
                    for chain in p_filament:
                        coords.append(chain.median(axis))

                    diff.append(max(coords)-min(coords))
                main_axis=['x', 'y', 'z'][diff.index(max(diff))]

                start = self.cur_pdb().p_cylinder['generated_axis']['start']
                end = self.cur_pdb().p_cylinder['generated_axis']['end']
                start[{'x':0, 'y':1, 'z':2}[axis]] = min(p_filament[0].data(axis))
                end[{'x':0, 'y':1, 'z':2}[axis]] = max(p_filament[-1].data(axis))

                cmd.pseudoatom('axis', pos=start, name='Start')
                cmd.pseudoatom('axis', pos=end, name='End')
                cmd.distance('Fibril_axis', 'axis and name Start', 'axis and name End')
                cmd.hide('labels', 'Fibril_axis')
                cmd.delete('axis')

            except Exception:
                #traceback.print_exc()
                print('No fibril axis were detected')
                self.error_message = 'No fibril axis were detected'

            self.fibril_axis = 'ON'


    @busy_decorator
    @suspend_decorator
    def do_select(self, name):
        """This activates first when clicking"""
        cmd.get_wizard().launcher("multi", "do_pick", 0)


    def do_pick(self, bondFlag):
        """Commands after clicking"""
        if bondFlag:
            return


        try:
            self.error_message = ''
            if self.picking_mode == 'reset_chains':
                atom_name = self.main_object + str(self.pick_count)
                self.pick_sheet(atom_name)
                self.cur_pdb().protofilament = str(self.pick_count)
                #print(self.pick_count)
                """
                if self.pick_count < 1:
                    self.sheet_count = len(pdb_reader('monomer_0.pdb').splice())
                    print('chains: '+str(self.sheet_count))
                """
                if self.pick_count and self.mode == 'evaluate':

                    self.calc_state = [1]

            elif self.picking_mode == 'connected':
                # get residue number
                try:
                    resn = get_data('sele and sele_0', 'resv')[0]

                    if resn in self.manual_connected_sections:
                        raise IndexError
                    if len(self.manual_connected_sections) >= 2:
                        a, b = self.manual_connected_sections[0:2]
                        if min(a, b) <= resn and resn <= max(a, b):
                            raise IndexError

                except IndexError:
                    cmd.select('sele', 'none')
                    print('Wrong residue picked.')
                    return

                cmd.select('sele', 'none')
                cmd.color('red', 'sele_0 and resi %s'%resn)
                self.manual_connected_sections.append(resn)
                print('Selected residue: %s (%s/4)'%(resn, len(self.manual_connected_sections)))
                print(len(self.manual_connected_sections))
                if len(self.manual_connected_sections) == 4:

                    cmd.save('Sheet_midobj.pdb', 'obj%s'%(self.cur_pdb().to_generate//2))
                    which_monomer = {}
                    for n in range(int(self.cur_pdb().protofilament)):
                        chains = pdb_reader("monomer_%s.pdb"%n).splice()
                        for chain in chains:
                            which_monomer[chain.lista[0].chain] = str(n)

                    # selection
                    cmd.select('sele', 'none')
                    cmd.select('calc_0', 'none')
                    cmd.select('calc_1', 'none')

                    print(self.manual_connected_sections)
                    a, b, c, d = self.manual_connected_sections
                    self.manual_connected_sections = []
                    for aa in range(min(a, b), max(a, b)+1):
                        cmd.select('calc_0', '(calc_0 or sele_0 and resi %s)'%aa, 0, 1, 0, 1)
                    for aa in range(min(c, d), max(c, d)+1):
                        cmd.select('calc_1', '(calc_1 or sele_0 and resi %s)'%aa, 0, 1, 0, 1)
                    cmd.select('calc_0-1', 'calc_0 or calc_1', 0, 1, 0, 1)

                    # color and save
                    cmd.color('grey', 'all')
                    cmd.color('blue', 'calc_0')
                    cmd.color('yellow', 'calc_1')
                    cmd.save('Sheet_0.pdb', 'calc_0')
                    cmd.save('Sheet_1.pdb', 'calc_1')
                    cmd.save('Sheet_0-1.pdb', 'calc_0-1')

                    # calc
                    sc_input(0, self.main_object)
                    sc_input(1, self.main_object)
                    sc, area1, area2 = sc_calc_out(1, self.main_object, self.path)
                    sc = str(sc)
                    scarea = str(round((float(area1)+float(area2))/(2*self.sheet_count), 3))
                    print('Sc = %s'%sc)

                    areaimol_calc_2(0, self.main_object, self.path, True)
                    areaimol_calc_2(1, self.main_object, self.path, True)
                    area, fingerprint = areaimol_out_4(1, self.main_object, which_monomer)
                    area = area / self.sheet_count * self.cur_pdb().to_generate
                    print('Buried area = %s'%(area))

                    #(length, heigth, length_coord) = perf_calc_fast(1, self.sheet_count, 'calc_0', 'calc_1')
                    perf_output = perf_calc_chains(1, self.sheet_count, self.cur_pdb().to_generate)
                    perfi_chains = []
                    for length, heigth, _ in perf_output:
                        if heigth == 0 or length == 0:
                            perfi_chains.append(0)
                        else:
                            perfi_chains.append(round(float(area)/heigth/length, 3))
                            #print(round(float(area)/heigth/length, 3))

                    perfi = round(sum(perfi_chains) / len(perfi_chains), 2)
                    print('SDi = %s'%perfi)

                    if self.pick_count == 0:
                        with open('acwi-out8.txt', 'w+') as f:
                            f.write(str(self.pick_count+1) +' - '+sc +' '+ str(area) +' '+ str(perfi)+' (%s-%s)+(%s-%s)\n'%(min(a,b), max(a,b), min(c,d), max(c,d)))
                    else:
                        with open('acwi-out8.txt', 'a+') as f:
                            f.write(str(self.pick_count+1) +' - '+sc +' '+ str(area) +' '+ str(perfi)+' (%s-%s)+(%s-%s)\n'%(min(a,b), max(a,b), min(c,d), max(c,d)))
                    self.pick_count += 1

            elif self.picking_mode == 'Delete':
                try:
                    self.pick_sheet('to_rm')
                    cmd.select('to_rm', 'bychain to_rm', 0, 1, 0, 1)
                    cmd.remove('to_rm')
                    # cmd.hide('everything', 'sele')

                except IndexError:
                    print('No sheet found')

            else:
                pass
        finally:
            self.busy = 0


    def pick_sheet(self, name, color=False):
        """Extend the pick selection to the whole sheet"""

        resn = get_data('sele', 'resv')[0]

        print('picked:', resn)
        #print(get_data('(bychain sele) and name O ', 'resv'))

        self.resvrange = get_data('(bychain sele) and name O ', 'resv')

        # sheet pick
        bond_dist = 3.4
        res_range = 'resi %s-%s'%(resn-3, resn+3)
        selector_text = 'bychain %s and ((name O within %s of (sele and %s and name N)) or (name N within %s of (sele and %s and name O)))' %(res_range, str(bond_dist), res_range, str(bond_dist), res_range)
        atoms = 0
        for n in range(50):
            cmd.select('sele', selector_text, 0, 1, 0, 1)
            atoms2 = cmd.count_atoms('sele')
            if atoms == atoms2:
                break
            atoms = atoms2

        # coloring
        try:
            cmd.color(self.colors[self.pick_count], 'sele and elem C')
        except IndexError:
            cmd.color('blue', 'sele') # if there are no more defined colors

        # saving
        cmd.select('sele_%s'%self.pick_count, 'sele')
        cmd.remove('sele_%s and not polymer.protein'%self.pick_count)
        cmd.save('monomer_%s.pdb'%self.pick_count, 'sele_%s'%self.pick_count)

        self.pick_count += 1

        cmd.delete('sele')
        cmd.unpick()
        cmd.select('none')
        cmd.refresh_wizard()


    def moving_window(self, atom_name):
        """Full hexa moving window method (found in development)"""

        #moving_window_method(self, atom_name)
        cmd.refresh_wizard()


    @busy_decorator
    @suspend_decorator
    def find_connections(self):
        """For each resiude it checks which other residues are in contact with."""
        for contact_file in glob.glob('contacts*.txt'):
            os.remove(contact_file)

        if self.mode == 'evaluate' and self.calc_state[0] != 2:
            print('Please use the methods in order.')
            self.error_message = 'Please use the methods in order.'
            return

        which_monomer = {}
        for n in range(int(self.cur_pdb().protofilament)):
            chains = pdb_reader("monomer_%s.pdb"%n).splice()
            for chain in chains:
                which_monomer[chain.lista[0].chain] = str(n)
        print('obj%s'%(self.cur_pdb().to_generate//2))
        cmd.save('Sheet_midobj.pdb', 'obj%s'%(self.cur_pdb().to_generate//2))

        cmd.remove('all and not polymer.protein')
        cmd.select('sele', 'none')
        for monomer in range(int(self.cur_pdb().protofilament)):
            cmd.select('sele', 'sele or sele_%s'%monomer)
        cmd.remove('not sele')
        cmd.select('sele', 'none')

        for monomer in range(int(self.cur_pdb().protofilament)):

            # get resv range

            cmd.select('sele', 'none')
            chains = pdb_reader('monomer_%s.pdb'%monomer).splice()
            x, y, z = chains[0].lista[0].xyz()
            cmd.select('sele', 'bychain sele or ( X>%s and X<%s and Y>%s and Y<%s and Z>%s and Z<%s)'%(x-0.005,x+0.005,y-0.005,y+0.005,z-0.005,z+0.005), 0, 1, 0, 1)

            self.sheet_count = len(pdb_reader('monomer_%s.pdb'%monomer).splice())

            # print(get_data('(bychain sele) and name CA ', 'resv'))
            self.resvrange = get_data('(bychain sele) and name CA ', 'resv')

            output = []

            cmd.alter('sele', 'b=0')

            for x, resi in enumerate(self.resvrange):
                print(resi)
                resi = str(resi)

                # check for no side chain = G
                if lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]]=='G':
                    print('G: no side chain ')

                    with open('contacts_%s.txt'%monomer, 'a') as f:
                        line=resi+' '+resi+'_%s\n'%monomer
                        f.write(line)

                    with open('contacts_info_%s.txt'%monomer, 'a') as f:
                        line=resi+' '+resi+'_%s_G_0_0\n'%monomer
                        f.write(line)

                    continue

                cmd.color('gray', 'all')
                # select sheet, color it and save it
                cmd.select('sele', 'sele_%s and resi %s and not bb.'%(monomer, resi), 0, 1, 0, 1)
                cmd.select('calc_1', 'sele', 0, 1, 0, 1)
                # safety check for empty selection
                if len(get_data('calc_1 ', 'resv')) == 0:
                    print('empty selection')
                    with open('contacts_%s.txt'%monomer, 'a') as f:
                        line=resi+' '+resi+'_%s\n'%monomer
                        f.write(line)
                    with open('contacts_info_%s.txt'%monomer, 'a') as f:
                        line=resi+' '+resi+'_%s_X_0_0\n'%monomer
                        f.write(line)
                    continue
                cmd.save('Sheet_1.pdb', 'calc_1')
                cmd.color('yellow', 'calc_1')

                # cmd.select('calc_0', 'sele_0 and not sele', 0, 1, 0, 1)
                # cmd.select('calc_0', '(sele_0 and not resi %s) or (all and not sele_0)'%resi, 0, 1, 0, 1)
                cmd.select('calc_0', 'byres all and ((sele_%s and resi %s) around 12)'%(monomer, resi), 0, 1, 0, 1)

                cmd.save('Sheet_0.pdb', 'calc_0')
                cmd.color('blue', 'calc_0')
                self.suspend_refresh()
                cmd.select('calc_0-1', 'calc_0 or calc_1', 0, 1, 0, 1)
                cmd.save('Sheet_0-1.pdb', 'calc_0-1')
                # gif
                # cmd.save('%s_iw_%s.png'%(self.main_object, x+1))

                # areaimol_input(1, self.main_object)
                areaimol_calc_2(0, self.main_object, self.path, True)
                areaimol_calc_2(1, self.main_object, self.path, True)

                #area1, area1X, areaX, carea1, carea1X, careaX = areaimol_out_3(1, self.main_object, self.sheet_count)
                # print('Buried area = %s'%(round((area1+areaX-area1X)/(2*self.sheet_count), 2)))

                # find and list contacts: residue number, monomer number, name, area aa, area bb
                #c, contacts = atomic_area_prot(1, self.main_object, which_monomer, self.sheet_count)

                print()
                #currently in failed state
                c, contacts = areaimol_out_4_atomic(1, self.main_object, which_monomer, self.sheet_count)

                # print(contacts)

                # split into to sheet
                main_contacts = []
                for contact in contacts:
                    main_contacts.append(contact[0]+' '+contact[1])
                main_contacts = list(dict.fromkeys(main_contacts))
                # print(main_contacts)

                # writing contacts was it like this? also change back G
                with open('contacts_%s.txt'%monomer, 'a') as f:
                    line = resi
                    for contact in main_contacts:
                        line = line + ' ' + contact.replace(' ', '_')
                    line = line + "\n"
                    f.write(line)

                # writing contacts_info
                with open('contacts_info_%s.txt'%monomer, 'a') as f:
                    line = resi+' '
                    for contact in contacts:
                        for n in contact:
                            line = line + str(n) + '_'
                        line = line.rstrip('_') + ' '
                    line = line + "\n"
                    f.write(line)

        if self.mode == 'evaluate':
            self.calc_state = [3]
        cmd.refresh_wizard()
        print('done')


    @busy_decorator
    @suspend_decorator
    def inteligent_window(self):
        """calculate based on connected sidechains"""

        if self.mode == 'evaluate' and self.calc_state[0] != 3:
            print('Please use the methods in order.')
            self.error_message = 'Please use the methods in order.'
            return

        cmd.remove('all and not polymer.protein')
        #self.sheet_count = len(pdb_reader('monomer_0.pdb').splice())

        which_monomer = {}
        for n in range(int(self.cur_pdb().protofilament)):
            chains = pdb_reader("monomer_%s.pdb"%n).splice()
            for chain in chains:
                which_monomer[chain.lista[0].chain] = str(n)

        self.sheet_count = len(pdb_reader('monomer_0.pdb').splice())*self.cur_pdb().to_generate
        cmd.save('Sheet_midobj.pdb', 'obj%s'%(self.cur_pdb().to_generate//2))

        cmd.remove('all and not polymer.protein')
        cmd.select('sele', 'none')
        for monomer in range(int(self.cur_pdb().protofilament)):
            cmd.select('sele', 'sele or sele_%s'%monomer)
        cmd.remove('not sele')
        cmd.select('sele', 'none')

        for monomer in range(int(self.cur_pdb().protofilament)):

            # get resv range
            self.resvrange = get_data('(bychain first sele_%s) and name CA '%monomer, 'resv')
            print(self.resvrange)

            self.window_current = 0

            output = []
            output_lengths = []

            all_contacts = []
            with open('contacts_%s.txt'%monomer, 'r') as f:
                data = f.read().splitlines()
                all_contacts = []
                for d in data:
                    main_contacts = []
                    for n in d.split()[1:]:
                        main_contacts.append(n.replace('_', ' '))
                    all_contacts.append(main_contacts)

            #intra contacts
            while self.window_current < len(self.resvrange):
                try:
                    resi = str(self.resvrange[self.window_current])
                    print(resi, 'intra')

                    """
                    #contact area
                    with open('contacts_info_%s.txt'%monomer, 'r') as f:
                        data = f.read().splitlines()
                        contacts = []
                        for n in data[self.window_current].split()[1:]:
                            contacts.append(n.replace("_", " ").split())
                    contact area debug
                    contacts_area_d = {}
                    for contact in contacts:
                        try:
                            contacts_area_d[contact[0]+' '+contact[1]] += float(contact[3]) + float(contact[4])
                        except:
                            contacts_area_d[contact[0]+' '+contact[1]] = float(contact[3]) + float(contact[4])
                    contacts_area = []
                    for contact in contacts_area_d:
                        contacts_area.append([contact, contacts_area_d[contact]/self.sheet_count])
                    contacts_area = sorted(contacts_area, key=lambda contact : contact[1])
                    print(contacts_area)
                    """

                    # selecting contacts
                    main_contacts = []
                    for contact in all_contacts[self.window_current]:
                        if contact.split()[1] == str(monomer):
                            main_contacts.append(contact)
                    cmd.select('sele', 'none')
                    cmd.select('calc_0', 'none')
                    cmd.select('calc_1', 'none')
                    sheet_0 = []
                    sheet_1 = []
                    for contact in main_contacts:
                        if contact.split()[1] == str(monomer) and abs(int(contact.split()[0])-int(resi)) < 4: #contact near window
                            cmd.select('calc_0', '(calc_0 or sele_%s and resi %s)'%(contact.split()[1], contact.split()[0]), 0, 1, 0, 1)
                            sheet_0.append(contact.replace(' ', '/'))
                        else: #contact elsewhere
                            cmd.select('calc_1', '(calc_1 or sele_%s and resi %s)'%(contact.split()[1], contact.split()[0]), 0, 1, 0, 1)
                            sheet_1.append(contact.replace(' ', '/'))

                    # color and save
                    cmd.color('grey', 'all')
                    cmd.color('blue', 'calc_0')
                    cmd.color('green', 'sele_%s and resi %s'%(monomer, resi))
                    cmd.color('yellow', 'calc_1')
                    cmd.select('calc_0-1', 'calc_0 or calc_1', 0, 1, 0, 1)
                    self.suspend_refresh()

                    print(sheet_0)
                    print(sheet_1)
                    time.sleep(0.1)
                    # escape if no interaction found
                    if len(sheet_0) + len(sheet_1) < 3 or len(sheet_1) < 1:
                        print('no contact')
                        output.append([resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], '-', '-', '-', '-', '-', '-'])
                        output_lengths.append('-')
                        continue

                    # save
                    cmd.save('Sheet_0.pdb', 'calc_0')
                    cmd.save('Sheet_1.pdb', 'calc_1')
                    cmd.save('Sheet_0-1.pdb', 'calc_0-1')

                    # calc
                    sc_input(0, self.main_object)
                    sc_input(1, self.main_object)
                    areaimol_calc_2(0, self.main_object, self.path, True)
                    areaimol_calc_2(1, self.main_object, self.path, True)

                    area, fingerprint = areaimol_out_4(1, self.main_object, which_monomer)
                    area = area / self.sheet_count * self.cur_pdb().to_generate
                    print('Buried area = %s'%(round(area,2)))
                    if len(fingerprint[0]) + len(fingerprint[1]) < 3:
                        print('not enough contact')
                        output.append([resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], '-', '-', '-', '-', '-', '-'])
                        output_lengths.append('-')
                        continue

                    sc, area1, area2 = sc_calc_out(1, self.main_object, self.path)
                    for sc_file in glob.glob('/tmp/raxis/sc*'):
                        os.remove(sc_file)
                    sc = str(sc)
                    scarea = str(round((float(area1)+float(area2))/(2*self.sheet_count), 3))
                    print('Sc = %s'%sc)

                    data_coord = save_coords('all')
                    #(length, heigth, length_coord) = perf_calc_fast(1, self.sheet_count, 'calc_0', 'calc_1')
                    perf_output = perf_calc_chains(1, self.sheet_count, self.cur_pdb().to_generate) #[(0, 0, 0), (0, 0, 0)]
                    reset_coords('all', data_coord)

                    perfi_chains = []
                    for length, heigth, _ in perf_output:
                        if heigth == 0 or length == 0:
                            perfi_chains.append(0)
                        else:
                            perfi_chains.append(round(float(area)/heigth/length, 3))
                            #print(round(float(area)/heigth/length, 3))

                    perfi = round(sum(perfi_chains) / len(perfi_chains), 2)
                    print('SDi = %s'%perfi)

                    output.append([resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], sc, scarea, str(round(area, 2)), str(perfi), '_'.join(sheet_0), '_'.join(sheet_1)])
                    output_lengths.append([resi, perf_output[0][2]]) #first lengthcoord
                    cmd.unpick()
                    cmd.deselect()
                    cmd.select('sele', 'none')

                except Exception:
                    #traceback.print_exc()
                    output.append([resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], '-', '-', '-', '-', '-', 'error'])
                    output_lengths.append('error')
                finally:
                    self.window_current += 1

            # save file
            # resi, name, sc, scarea, ba, sdi, resi_0, resi_1
            with open('acwi-out6i_%s.txt'%monomer, 'w+') as f:
                for out in output:
                    f.write(' '.join(out) + '\n')
            with open('acwi-out6i_%s-lengths.txt'%monomer, 'w+') as f:
                for out in output_lengths:
                    f.write(' '.join(out) + '\n')

            #inter contacts
            self.window_current = 0
            while self.window_current < len(self.resvrange):
                try:
                    resi = str(self.resvrange[self.window_current])
                    print(resi, 'inter')

                    # selecting contacts
                    main_contacts = all_contacts[self.window_current]
                    if all(contact.split()[1] == str(monomer) for contact in main_contacts):
                        continue

                    cmd.select('sele', 'none')
                    cmd.select('calc_0', 'none')
                    cmd.select('calc_1', 'none')
                    sheet_0 = []
                    sheet_1 = []
                    for contact in main_contacts:
                        if contact.split()[1] == str(monomer) and abs(int(contact.split()[0])-int(resi)) < 4: #contact near window
                            cmd.select('calc_0', '(calc_0 or sele_%s and resi %s)'%(contact.split()[1], contact.split()[0]), 0, 1, 0, 1)
                            sheet_0.append(contact.replace(' ', '/'))
                        else: #contact elsewhere
                            cmd.select('calc_1', '(calc_1 or sele_%s and resi %s)'%(contact.split()[1], contact.split()[0]), 0, 1, 0, 1)
                            sheet_1.append(contact.replace(' ', '/'))

                    # color and save
                    cmd.color('grey', 'all')
                    cmd.color('blue', 'calc_0')
                    cmd.color('green', 'sele_%s and resi %s'%(monomer, resi))
                    cmd.color('yellow', 'calc_1')
                    cmd.select('calc_0-1', 'calc_0 or calc_1', 0, 1, 0, 1)
                    self.suspend_refresh()

                    time.sleep(0.1)
                    # escape if no interaction found
                    print(sheet_0)
                    print(sheet_1)

                    if len(sheet_0) + len(sheet_1) < 3 or len(sheet_1) < 1:
                        print('no contact')
                        output[self.window_current] = [resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], '-', '-', '-', '-', '-', '-']
                        output_lengths[self.window_current] = '-'
                        continue

                    # save
                    cmd.save('Sheet_0.pdb', 'calc_0')
                    cmd.save('Sheet_1.pdb', 'calc_1')
                    cmd.save('Sheet_0-1.pdb', 'calc_0-1')

                    # calc
                    sc_input(0, self.main_object)
                    sc_input(1, self.main_object)
                    areaimol_calc_2(0, self.main_object, self.path, True)
                    areaimol_calc_2(1, self.main_object, self.path, True)

                    area, fingerprint = areaimol_out_4(1, self.main_object, which_monomer)
                    area = area / self.sheet_count * self.cur_pdb().to_generate
                    print('Buried area = %s'%(round(area,2)))
                    if len(fingerprint[0]) + len(fingerprint[1]) < 3:
                        print('not enough contact')
                        output[self.window_current] = [resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], '-', '-', '-', '-', '-', '-']
                        output_lengths[self.window_current] = '-'
                        continue

                    sc, area1, area2 = sc_calc_out(1, self.main_object, self.path)
                    for sc_file in glob.glob('/tmp/raxis/sc*'):
                        os.remove(sc_file)
                    sc = str(sc)
                    scarea = str(round((float(area1)+float(area2))/(2*self.sheet_count),3))
                    print('Sc = %s'%sc)

                    data_coord = save_coords('all')
                    #(length, heigth, length_coord) = perf_calc_fast(1, self.sheet_count, 'calc_0', 'calc_1')
                    perf_output = perf_calc_chains(1, self.sheet_count, self.cur_pdb().to_generate) #[(0, 0, 0),(0, 0, 0)]
                    reset_coords('all', data_coord)

                    perfi_chains = []
                    for length, heigth, _ in perf_output:
                        if heigth == 0 or length == 0:
                            perfi_chains.append(0)
                        else:
                            perfi_chains.append(round(float(area)/heigth/length, 3))
                            #print(round(float(area)/heigth/length, 3))

                    perfi = round(sum(perfi_chains) / len(perfi_chains), 2)
                    print('SDi = %s'%perfi)

                    output[self.window_current] = [resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], sc, scarea, str(round(area,2)), str(perfi), '_'.join(sheet_0), '_'.join(sheet_1)]
                    output_lengths[self.window_current] = [resi, perf_output[0][2]]
                    cmd.unpick()
                    cmd.deselect()
                    cmd.select('sele', 'none')

                except Exception:
                    #traceback.print_exc()
                    output[self.window_current] = [resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], '-', '-', '-', '-', '-', 'error']
                    output_lengths[self.window_current] = 'error'
                finally:
                    self.window_current += 1

            # save file
            # resi, name, sc, scarea, ba, sdi, resi_0, resi_1
            with open('acwi-out6_%s.txt'%monomer, 'w+') as f:
                for out in output:
                    f.write(' '.join(out) + '\n')
            with open('acwi-out6_%s-lengths.txt'%monomer, 'w+') as f:
                for out in output_lengths:
                    f.write(' '.join(out) + '\n')

        if self.mode == 'evaluate':
            self.calc_state = [4]
            cmd.delete('lengths')
            cmd.delete('calc_*')
            
        cmd.refresh_wizard()
        print('done')


    @busy_decorator
    @suspend_decorator
    def inteligent_window_recalc(self):
        """recalculate inteligent window, if there was a change"""

        cmd.remove('all and not polymer.protein')
        #self.sheet_count = len(pdb_reader('monomer_0.pdb').splice())

        which_monomer = {}
        for n in range(int(self.cur_pdb().protofilament)):
            chains = pdb_reader("monomer_%s.pdb"%n).splice()
            for chain in chains:
                which_monomer[chain.lista[0].chain] = str(n)

        self.sheet_count = len(pdb_reader('monomer_0.pdb').splice())*self.cur_pdb().to_generate
        cmd.save('Sheet_midobj.pdb', 'obj%s'%(self.cur_pdb().to_generate//2))

        for monomer in range(int(self.cur_pdb().protofilament)):

            # get resv range
            self.resvrange = get_data('(bychain first sele_%s) and name CA '%monomer, 'resv')
            print(self.resvrange)

            self.window_current = 0

            output = []
            output_lengths = []

            all_contacts = []
            with open('contacts_%s.txt'%monomer, 'r') as f:
                data = f.read().splitlines()
                all_contacts = []
                for d in data:
                    main_contacts = []
                    for n in d.split()[1:]:
                        main_contacts.append(n.replace('_', ' '))
                    all_contacts.append(main_contacts)

            #read back save file
            with open('acwi-out6i_%s.txt'%monomer, 'r') as f:
                output = list(map(str.split, f.readlines()))
            with open('acwi-out6i_%s-lengths.txt'%monomer, 'r') as f:
                output_lengths = list(map(str.split, f.readlines()))
            with open('/home/raxis/pythoncodes/ACWF-public/%s/acwi-out6i_%s.txt'%(self.cur_pdb().code, monomer), 'r') as f:
                output_old = list(map(str.split, f.readlines()))

            #intra contacts
            while self.window_current < len(self.resvrange):
                try:
                    #self.window_current += 2
                    resi = str(self.resvrange[self.window_current])
                    print(resi, 'intra')

                    # 0-resi, 1-name, 2-sc, 3-scarea, 4-ba, 5-sdi, 6-resi_0, 7-resi_1
                    if output[self.window_current][6] == '-' and output_old[self.window_current][6] == '-':
                        continue

                    if set(output[self.window_current][6].split('_')) == set(output_old[self.window_current][6].split('_')) and \
                       set(output[self.window_current][7].split('_')) == set(output_old[self.window_current][7].split('_')):

                        print('same')
                        output[self.window_current][2] = output_old[self.window_current][2] #sc
                        output[self.window_current][3] = output_old[self.window_current][3] #sc area
                        #output[self.window_current][4] = output_old[self.window_current][4] #ba
                        #output[self.window_current][5] = output_old[self.window_current][5] #sdi
                        continue

                    print('different')

                    # selecting contacts
                    main_contacts = []
                    for contact in all_contacts[self.window_current]:
                        if contact.split()[1] == str(monomer):
                            main_contacts.append(contact)
                    cmd.select('sele', 'none')
                    cmd.select('calc_0', 'none')
                    cmd.select('calc_1', 'none')
                    sheet_0 = []
                    sheet_1 = []
                    for contact in main_contacts:
                        if contact.split()[1] == str(monomer) and abs(int(contact.split()[0])-int(resi)) < 4: #contact near window
                            cmd.select('calc_0', '(calc_0 or sele_%s and resi %s)'%(contact.split()[1],contact.split()[0]), 0, 1, 0, 1)
                            sheet_0.append(contact.replace(' ', '/'))
                        else: #contact elsewhere
                            cmd.select('calc_1', '(calc_1 or sele_%s and resi %s)'%(contact.split()[1],contact.split()[0]), 0, 1, 0, 1)
                            sheet_1.append(contact.replace(' ', '/'))

                    # color and save
                    cmd.color('grey', 'all')
                    cmd.color('blue', 'calc_0')
                    cmd.color('green', 'sele_%s and resi %s'%(monomer, resi))
                    cmd.color('yellow', 'calc_1')
                    cmd.select('calc_0-1', 'calc_0 or calc_1', 0, 1, 0, 1)
                    self.suspend_refresh()

                    print(sheet_0)
                    print(sheet_1)
                    time.sleep(0.1)
                    # escape if no interaction found
                    if len(sheet_0) + len(sheet_1) < 3 or len(sheet_1) < 1:
                        print('no contact')
                        output[self.window_current] = [resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], '-', '-', '-', '-', '-', '-']
                        output_lengths[self.window_current] = '-'
                        continue

                    # save
                    cmd.save('Sheet_0.pdb', 'calc_0')
                    cmd.save('Sheet_1.pdb', 'calc_1')
                    cmd.save('Sheet_0-1.pdb', 'calc_0-1')

                    # calc
                    sc_input(0, self.main_object)
                    sc_input(1, self.main_object)
                    areaimol_calc_2(0, self.main_object, self.path,True)
                    areaimol_calc_2(1, self.main_object, self.path,True)

                    area, fingerprint = areaimol_out_4(1, self.main_object, which_monomer)
                    area = area / self.sheet_count * self.cur_pdb().to_generate
                    print('Buried area = %s'%(area))
                    if len(fingerprint[0]) + len(fingerprint[1]) < 3:
                        print('not enough contact')
                        output[self.window_current] = [resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], '-', '-', '-', '-', '-', '-']
                        output_lengths[self.window_current] = '-'
                        continue

                    sc, area1, area2 = sc_calc_out(1, self.main_object, self.path)
                    for sc_file in glob.glob('/tmp/raxis/sc*'):
                        os.remove(sc_file)
                    sc = str(sc)
                    scarea = str(round((float(area1)+float(area2))/(2*self.sheet_count),3))
                    print('Sc = %s'%sc)

                    #(length, heigth, length_coord) = perf_calc_fast(1, self.sheet_count, 'calc_0', 'calc_1')
                    perf_output = perf_calc_chains(1, self.sheet_count, self.cur_pdb().to_generate)
                    perfi_chains = []
                    for length, heigth, _ in perf_output:
                        if heigth == 0 or length == 0:
                            perfi_chains.append(0)
                        else:
                            perfi_chains.append(round(float(area)/heigth/length, 3))
                            #print(round(float(area)/heigth/length, 3))

                    perfi = round(sum(perfi_chains) / len(perfi_chains), 2)
                    print('SDi = %s'%perfi)

                    output[self.window_current] = [resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], sc, scarea, str(round(area,2)), str(perfi), '_'.join(sheet_0), '_'.join(sheet_1)]
                    output_lengths[self.window_current] = [resi, perf_output[0][2]] #first lengthcoord
                    cmd.unpick()
                    cmd.deselect()
                    cmd.select('sele', 'none')

                except Exception:
                    #traceback.print_exc()
                    output[self.window_current] = [resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], '-', '-', '-', '-', '-', 'error']
                    output_lengths[self.window_current] = 'error'
                    return
                finally:
                    self.window_current += 1

            # save file
            # resi, name, sc, scarea, ba, sdi, resi_0, resi_1
            with open('acwi-out6i_%s.txt'%monomer, 'w+') as f:
                for out in output:
                    f.write(' '.join(out) + '\n')
            with open('acwi-out6i_%s-lengths.txt'%monomer, 'w+') as f:
                for out in output_lengths:
                    f.write(' '.join(out) + '\n')

            with open('acwi-out6_%s.txt'%monomer, 'r') as f:
                output = list(map(str.split, f.readlines()))
            with open('acwi-out6_%s-lengths.txt'%monomer, 'r') as f:
                output_lengths = list(map(str.split, f.readlines()))
            with open('/home/raxis/pythoncodes/ACWF-public/%s/acwi-out6_%s.txt'%(self.cur_pdb().code, monomer), 'r') as f:
                output_old = list(map(str.split, f.readlines()))

            #inter contacts
            self.window_current = 0
            while self.window_current < len(self.resvrange):
                try:
                    resi = str(self.resvrange[self.window_current])
                    print(resi, 'inter')

                    # 0-resi, 1-name, 2-sc, 3-scarea, 4-ba, 5-sdi, 6-resi_0, 7-resi_1
                    if output[self.window_current][6] == '-' and output_old[self.window_current][6] == '-':
                        continue
                    if set(output[self.window_current][6].split('_')) == set(output_old[self.window_current][6].split('_')) and \
                       set(output[self.window_current][7].split('_')) == set(output_old[self.window_current][7].split('_')):

                        print('same')
                        output[self.window_current][2] = output_old[self.window_current][2] #sc
                        output[self.window_current][3] = output_old[self.window_current][3] #sc area
                        #output[self.window_current][4] = output_old[self.window_current][4] #ba
                        #output[self.window_current][5] = output_old[self.window_current][5] #sdi
                        continue

                    print('different')
                    # selecting contacts
                    main_contacts = all_contacts[self.window_current]
                    if all(contact.split()[1] == str(monomer) for contact in main_contacts):
                        continue

                    cmd.select('sele', 'none')
                    cmd.select('calc_0', 'none')
                    cmd.select('calc_1', 'none')
                    sheet_0 = []
                    sheet_1 = []
                    for contact in main_contacts:
                        if contact.split()[1] == str(monomer) and abs(int(contact.split()[0])-int(resi)) < 4: #contact near window
                            cmd.select('calc_0', '(calc_0 or sele_%s and resi %s)'%(contact.split()[1],contact.split()[0]), 0, 1, 0, 1)
                            sheet_0.append(contact.replace(' ', '/'))
                        else: #contact elsewhere
                            cmd.select('calc_1', '(calc_1 or sele_%s and resi %s)'%(contact.split()[1],contact.split()[0]), 0, 1, 0, 1)
                            sheet_1.append(contact.replace(' ', '/'))

                    # color and save
                    cmd.color('grey', 'all')
                    cmd.color('blue', 'calc_0')
                    cmd.color('green', 'sele_%s and resi %s'%(monomer, resi))
                    cmd.color('yellow', 'calc_1')
                    cmd.select('calc_0-1', 'calc_0 or calc_1', 0, 1, 0, 1)
                    self.suspend_refresh()

                    time.sleep(0.1)
                    # escape if no interaction found
                    print(sheet_0)
                    print(sheet_1)

                    if len(sheet_0) + len(sheet_1) < 3 or len(sheet_1) < 1:
                        print('no contact')
                        output[self.window_current] = [resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], '-', '-', '-', '-', '-', '-']
                        output_lengths[self.window_current] = '-'
                        continue

                    # save
                    cmd.save('Sheet_0.pdb', 'calc_0')
                    cmd.save('Sheet_1.pdb', 'calc_1')
                    cmd.save('Sheet_0-1.pdb', 'calc_0-1')

                    # calc
                    sc_input(0, self.main_object)
                    sc_input(1, self.main_object)
                    areaimol_calc_2(0, self.main_object, self.path,True)
                    areaimol_calc_2(1, self.main_object, self.path,True)

                    area, fingerprint = areaimol_out_4(1, self.main_object, which_monomer)
                    area = area / self.sheet_count * self.cur_pdb().to_generate
                    print('Buried area = %s'%(area))
                    if len(fingerprint[0]) + len(fingerprint[1]) < 3:
                        print('not enough contact')
                        output[self.window_current] = [resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], '-', '-', '-', '-', '-', '-']
                        output_lengths[self.window_current] = '-'
                        continue

                    sc, area1, area2 = sc_calc_out(1, self.main_object, self.path)
                    for sc_file in glob.glob('/tmp/raxis/sc*'):
                        os.remove(sc_file)
                    sc = str(sc)
                    scarea = str(round((float(area1)+float(area2))/(2*self.sheet_count),3))
                    print('Sc = %s'%sc)

                    #(length, heigth, length_coord) = perf_calc_fast(1, self.sheet_count, 'calc_0', 'calc_1')
                    perf_output = perf_calc_chains(1, self.sheet_count, self.cur_pdb().to_generate)
                    perfi_chains = []
                    for length, heigth, _ in perf_output:
                        if heigth == 0 or length == 0:
                            perfi_chains.append(0)
                        else:
                            perfi_chains.append(round(float(area)/heigth/length, 3))
                            #print(round(float(area)/heigth/length, 3))

                    perfi = round(sum(perfi_chains) / len(perfi_chains), 2)
                    print('SDi = %s'%perfi)

                    output[self.window_current] = [resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], sc, scarea, str(round(area,2)), str(perfi), '_'.join(sheet_0), '_'.join(sheet_1)]
                    output_lengths[self.window_current] = [resi, perf_output[0][2]]
                    cmd.unpick()
                    cmd.deselect()
                    cmd.select('sele', 'none')

                except Exception:
                    #traceback.print_exc()
                    output[self.window_current] = [resi, lettercodes[get_data('sele_%s and resi %s'%(monomer, resi), 'resn')[0]], '-', '-', '-', '-', '-', 'error']
                    output_lengths[self.window_current] = 'error'
                finally:
                    self.window_current += 1

            # save file
            # resi, name, sc, scarea, ba, sdi, resi_0, resi_1
            with open('acwi-out6_%s.txt'%monomer, 'w+') as f:
                for out in output:
                    f.write(' '.join(out) + '\n')
            with open('acwi-out6_%s-lengths.txt'%monomer, 'w+') as f:
                for out in output_lengths:
                    f.write(' '.join(out) + '\n')

        print('done')


    @busy_decorator
    @suspend_decorator
    def inteligent_window_gif(self):

        # get resv range
        cmd.remove('all and not polymer.protein')
        cmd.select('sele', 'none')
        chains = pdb_reader("monomer_0.pdb").splice()
        x, y, z = chains[0].lista[0].xyz()
        cmd.select('sele', 'bychain sele or ( X>%s and X<%s and Y>%s and Y<%s and Z>%s and Z<%s)'%(x-0.005,x+0.005,y-0.005,y+0.005,z-0.005,z+0.005), 0, 1, 0, 1)

        # print(get_data('(bychain sele) and name O ', 'resv'))
        self.resvrange = get_data('(bychain sele) and name O ', 'resv')
        self.sheet_count = len(pdb_reader("monomer_0.pdb").splice())

        for file in glob.glob('%s_iw_*.png'%self.main_object):
            os.remove(file)

        for x, resi in enumerate(self.resvrange):

            print(resi)
            resi = str(resi)

            # reading contacts
            with open('contacts.txt', 'r') as f:
                data = f.read().splitlines()
                main_contacts = []
                for n in data[x].split()[1:]:
                    main_contacts.append(n.replace("_", " "))

            with open('contacts_info.txt', 'r') as f:
                data = f.read().splitlines()
                contacts = []
                for n in data[x].split()[1:]:
                    contacts.append(n.replace("_", " ").split())

            contacts_area_d = {}
            for contact in contacts:
                try:
                    contacts_area_d[contact[0]+' '+contact[1]] += float(contact[3]) + float(contact[4])
                except KeyError:
                    contacts_area_d[contact[0]+' '+contact[1]] = float(contact[3]) + float(contact[4])
            contacts_area = []
            for contact in contacts_area_d:
                contacts_area.append([contact, contacts_area_d[contact]/self.sheet_count])
            contacts_area = sorted(contacts_area, key=lambda contact : contact[1])
            # print(contacts_area)

            # selecting contacts
            cmd.select('sele', 'none')
            cmd.select('calc_0', 'none')
            cmd.select('calc_1', 'none')
            sheet_0 = []
            sheet_1 = []
            for contact in main_contacts:
                if contact.split()[1] == '0' and abs(int(contact.split()[0])-int(resi)) <4:
                    cmd.select('calc_0', '(calc_0 or sele_%s and resi %s)'%(contact.split()[1], contact.split()[0]), 0, 1, 0, 1)
                    sheet_0.append(contact.split()[0])
                else:
                    cmd.select('calc_1', '(calc_1 or sele_%s and resi %s)'%(contact.split()[1], contact.split()[0]), 0, 1, 0, 1)
                    sheet_1.append(contact.split()[0])

            # color and save
            cmd.color('grey', 'all')
            cmd.color('blue', 'calc_0')
            cmd.color('green', 'sele_0 and resi %s'%resi)
            cmd.color('yellow', 'calc_1')
            cmd.select('calc_0-1', 'calc_0 or calc_1', 0, 1, 0, 1)
            self.suspend_refresh()

            # screenshot
            cmd.save('%s_iw_%s.png'%(self.main_object, x+1))
            # time.sleep(0.5)


    @busy_decorator
    @suspend_decorator
    def connected_sections(self):
        """Edit, display and calculate the connected sections."""

        # get resv range
        cmd.remove('all and not polymer.protein')
        cmd.select('sele', 'none')
        chains = pdb_reader("monomer_0.pdb").splice()
        chain_count = len(chains)
        x, y, z = chains[0].lista[0].xyz()
        cmd.select('sele', 'bychain sele or ( X>%s and X<%s and Y>%s and Y<%s and Z>%s and Z<%s)'%(x-0.005,x+0.005,y-0.005,y+0.005,z-0.005,z+0.005), 0, 1, 0, 1)

        # print(get_data('(bychain sele) and name O ', 'resv'))
        self.resvrange = get_data('(bychain sele) and name O ', 'resv')
        self.sheet_count = len(pdb_reader("monomer_0.pdb").splice())

        find_connected_sections()

        # reading
        with open('connected_section.txt', 'r') as f:
            data = f.read().splitlines()
            connected_section = []
            for n in data:
                connected_section.append(n.split("_"))

        cmd.select('sele', 'none')
        cmd.select('calc_0', 'none')
        cmd.select('calc_1', 'none')


        if self.pick_count_cs > len(connected_section)-1:
            print('no more connected sections, resetting to first')
            self.pick_count_cs = 0
            self.menu_update()
            self.reset()
            return

        # selection
        next_part = False
        section=[[], []]
        for aa in connected_section[self.pick_count_cs]:
            if aa == '-':
                next_part = True
                continue
            if not next_part:
                section[0].append(aa)
                cmd.select('calc_0', '(calc_0 or sele_0 and resi %s)'%aa, 0, 1, 0, 1)
            else:
                section[1].append(aa)
                cmd.select('calc_1', '(calc_1 or sele_0 and resi %s)'%aa, 0, 1, 0, 1)
        cmd.select('calc_0-1', 'calc_0 or calc_1', 0, 1, 0, 1)
        print(connected_section[self.pick_count_cs])

        # tips for cutting it half manually
        for x in range(2):
            length = len(section[x])
            if length < 4:
                continue
            if length % 2 == 0:
                section[x].pop(int(length/2-1))
                section[x].pop(int(length/2-1))
                section[x].insert(int(length/2-1), '-')
            else:
                section[x].pop(int(length/2-1))
                section[x].pop(int(length/2-1))
                section[x].pop(int(length/2-1))
                section[x].insert(int(length/2-1), '-')
        # print(section[0])
        # print(section[1])

        # color and save
        cmd.color('grey', 'all')
        cmd.color('blue', 'calc_0')
        cmd.color('yellow', 'calc_1')
        cmd.save('Sheet_0.pdb', 'calc_0')
        cmd.save('Sheet_1.pdb', 'calc_1')
        cmd.save('Sheet_0-1.pdb', 'calc_0-1')


        # calc
        sc_input(0, self.main_object)
        sc_input(1, self.main_object)
        sc, area1, area2 = sc_calc_out(1, self.main_object, self.path)
        sc = str(sc)
        scarea = str(round((float(area1)+float(area2))/(2*self.sheet_count), 3))
        print('Sc = %s'%sc)

        areaimol_calc_2(0, self.main_object, self.path, True)
        areaimol_calc_2(1, self.main_object, self.path, True)
        area1, area1X, areaX, carea1, carea1X, careaX = areaimol_out_3(1, self.main_object, self.sheet_count)
        area = str(round((area1+areaX-area1X)/(2*self.sheet_count), 2))
        print('Buried area = %s'%(area))

        if self.pick_count_cs == 0:
            with open('acwi-out7.txt', 'w+') as f:
                f.write(str(self.pick_count_cs) +' - '+sc +' '+ scarea +' '+ area+'\n')
        else:
            with open('acwi-out7.txt', 'a+') as f:
                f.write(str(self.pick_count_cs) +' - '+sc +' '+ scarea +' '+ area+'\n')

        self.pick_count_cs += 1

        print('done')
        cmd.unpick()
        cmd.deselect()
        cmd.select('sele', 'none')
        self.menu_update()
        self.reset()
        print('done')


    def exp_image(self):
        """Exports a png picture."""

        cmd.ray('3000', '3000')
        n_output = 1
        while os.path.isfile('%s_%s.png'%(self.main_object, n_output)):
            n_output += 1
        cmd.save('%s_%s.png'%(self.main_object, n_output))
        print('%s_%s.png saved'%(self.main_object, n_output))

        self.reset()


    def set_value(self, value_type, value, object_type):
        """Change a saved value to an other."""

        if object_type == 'wizard':
            setattr(self, value_type, value)
        elif object_type == 'amyloid':
            setattr(self.pdb_list[self.pdb_count], value_type, value)
        elif object_type == 'interface':
            setattr(self.pdb_list[self.pdb_count].interface[self.interface_count-1], value_type, value)
        if object_type == 'coloring':
            setattr(self, value_type, value)
            self.color(method=self.color_last_method)

        self.menu_update()
        self.reset()


    def toggle_value(self, value_type, values):
        """Toggle a saved value in lists."""
        for value in values:
            if value in getattr(self, value_type + '_types'):
                if value in getattr(self, value_type):
                    getattr(self, value_type).remove(value)
                else:
                    getattr(self, value_type).append(value)

        self.menu_update()
        self.reset()


    def def_changer(self, setting):
        """Toggles default settings."""

        setattr(self, setting, {'ON':'OFF', 'OFF':'ON'}[getattr(self, setting)])

        self.menu_update()
        self.reset()


    def cur_pdb(self):
        """Current pdb structure."""
        return self.pdb_list[self.pdb_count]


    def cur_iface(self, i=0):
        """Current interface of pdb structure."""
        return self.pdb_list[self.pdb_count].interface[self.interface_count+i]


    def automatic_calc(self):
        """Calculates all possible interfaces of the assymetric unit."""

        self.menu_update()


    def save_to_database(self):
        """Save finished evaluation to the personal database."""

        try:
            self.main_object = cmd.get_names()[0]
        except IndexError:
            print('No structure detected')
            self.error_message = 'No structure detected'
            cmd.refresh_wizard()
            return

        if os.path.isdir(self.main_object):
            shutil.rmtree(self.main_object)

        try:
            os.makedirs(self.main_object)
        except OSError:
            print('Cannot create a folder')
            self.error_message = 'Cannot create a folder'
            cmd.refresh_wizard()
            return

        shutil.copyfile('%s.pdb'%self.main_object, '%s/%s.pdb'%(self.main_object, self.main_object))

        print('\nSaving to database')

        self.cur_pdb().code = self.main_object
        self.cur_pdb().prot = '?'
        self.cur_pdb().method = 'ELECTRON_MICROSCOPY'
        self.cur_pdb().seq = '-'
        self.cur_pdb().clas = '-'
        self.cur_pdb().protofilament = '-'
        self.cur_pdb().seq_len = '-'
        self.cur_pdb().sec_class = '-'

        # changing to database
        self.pdb_count = len(self.pdb_list.lista)-1
        cmd.cd(self.main_path)
        write_save('acw_personal_database.json', self.pdb_list)

        self.pdb_list.add_blank()
        self.next_peptide(x=self.pdb_count)
        self.manual_mode()
        self.menu_update()


    def delete_from_database(self, answer=0):
        """Delete a structure from the personal database."""

        if self.overwrite_state == 0:
            self.overwrite_state = 1
            self.menu_update()

        else:

            self.overwrite_state = 0
            self.menu_update()

            if answer:
                print('%s was deleted'%self.cur_pdb().code)
                self.error_message = '%s was deleted'%self.cur_pdb().code
                x = self.pdb_count

                cmd.cd(self.main_path)
                amy_lista = read_save("amy_flength_saveout_10_bad.json")
                amy_lista.add(self.cur_pdb())
                write_save("amy_flength_saveout_10_bad.json", amy_lista)
                cmd.cd(self.cur_pdb().code)

                self.next_peptide(-1)
                self.pdb_list.lista.pop(x)
                self.next_peptide(0)


    def overwrite(self, answer=0):
        """Prompt the user if they want to overwrite a previous output folder."""

        self.overwrite_state = 0
        self.menu_update()
        if answer:
            for n, amy in enumerate(self.pdb_list.lista):
                if amy.code == self.main_object:
                    if amy.source != 'PDB':
                        self.pdb_list.lista.pop(n)
                        shutil.rmtree(self.main_object)
                        break
                    else:
                        print('This name is already used in the PDB database, cannot be overwritten')
                        self.error_message = 'This name is already used in the PDB database, cannot be overwritten'
                        return

            else:
                shutil.rmtree(self.main_object)
            print('Files were deleted')
            self.error_message = 'Files were deleted'
            cmd.get_wizard().launcher("multi", "setup")


    def exporter(self):
        """Export current selection to json and txt."""

        return


    def extend_filaments(self):
        """Extend filaments if needed"""

        cmd.save('Sheet_0.pdb', 'sele_0')
        p_filament = pdb_reader("Sheet_0.pdb").splice()
        chain_count = len(p_filament)
        diff = []
        for axis in ['x', 'y', 'z']:
            coords = []
            for chain in p_filament:
                coords.append(chain.median(axis))

            diff.append(max(coords)-min(coords))
        main_axis=['x', 'y', 'z'][diff.index(max(diff))]
        #print(main_axis)

        chains_mobile = []
        chains_target = []
        h_chain_min = []
        h_chain_max = []
        h_repeat = []

        for n in range(int(self.cur_pdb().protofilament)):
            cmd.save('Sheet_0.pdb', 'sele_%s'%n)
            p_filament = pdb_reader("Sheet_0.pdb").splice()
            chains = [chain.lista[0].chain for chain in p_filament]
            chains_mobile.append(chains[:-1])
            chains_target.append(chains[1:])
            h_chain_min.append(min(p_filament[0].data(main_axis)))
            h_chain_max.append(max(p_filament[-1].data(main_axis)))
            h_repeat.append(abs(p_filament[0].median(main_axis) - p_filament[1].median(main_axis))) #4.8

        #print(chains_mobile)
        #print(chains_target)
        if not all(map(lambda x: len(x) == len(chains_mobile[0]), chains_mobile)):
            print('different number of chains in protofilaments')
            raise ValueError
        #print(h_repeat)
        h_chain = abs(max(h_chain_max)-min(h_chain_min))
        h_repeat = sum(h_repeat)/len(h_repeat)
        #print(h_repeat)
        #print(h_chain)

        cmd.save('Sheet_0.pdb', 'sele_0')
        p_filament = pdb_reader("Sheet_0.pdb").splice()
        chain_count = len(p_filament)
        #print(chain_count)

        to_generate = int(h_chain*2 + h_repeat*2 + 3.8) / h_repeat / chain_count #upper limit
        #print(to_generate)
        to_generate = (int(to_generate)+1)//2*2+1 #rounded up to nearest odd number
        #print(to_generate)
        if to_generate < 3:
            to_generate = 3

        #to_generate = int(90/chain_count)
        rmsd = 0
        for c in range(1, to_generate):
            if c == 1:
                cmd.copy('obj1', self.main_object)
            else:
                cmd.copy('obj%s'%c, 'obj%s'%(c-1))

            for n in range(chain_count):
                to_fit = []
                for x in range(int(self.cur_pdb().protofilament)):
                    for chain_n in range(len(chains_mobile[0])):
                        mobile = 'obj%s and bb. and chain %s and name CA'%(c, chains_mobile[x][chain_n])
                        target = 'obj%s and bb. and chain %s and name CA'%(c, chains_target[x][chain_n])

                        """
                        #test selection
                        xyz1 = zip(get_data(mobile, 'x'), get_data(mobile, 'y'), get_data(mobile, 'z'))
                        xyz2 = zip(get_data(target, 'x'), get_data(target, 'y'), get_data(target, 'z'))
                        for coord1, coord2 in zip(xyz1, xyz2):
                            cmd.distance('d', 'x=%s and y=%s and z=%s'%coord1, 'x=%s and y=%s and z=%s'%coord2)
                        cmd.hide('labels')
                        """
                        to_fit.append(mobile)
                        to_fit.append(target)

                rmsd = cmd.pair_fit(*to_fit)

        for n in range(int(self.cur_pdb().protofilament)):
            cmd.save('Sheet_0.pdb', 'sele_%s'%n)
            chains = [chain.lista[0].chain for chain in pdb_reader("Sheet_0.pdb").splice()]
            cmd.select('sele_%s'%n, 'none')
            cmd.select('sele_%s'%n, 'chain '+'+'.join(chains), 0, 1, 0, 1)

        color_methods(self, 'protofilaments__')
        return str(rmsd), to_generate


    def test(self):

        #extend with rise and turn

        """
        p0 = np.array(self.cur_pdb().p_cylinder['original_axis']['start'])
        q0 = np.array(self.cur_pdb().p_cylinder['original_axis']['end'])

        cmd.save('Sheet_0.pdb', self.main_object + ' and sele_0') #'sele_0'
        p_filament = pdb_reader("Sheet_0.pdb").splice()

        p = get_closest_point_3d(p0, q0, np.array(p_filament[0].center()).reshape(1, 3))[0]
        q = get_closest_point_3d(p0, q0, np.array(p_filament[-1].center()).reshape(1, 3))[0]

        angle = self.cur_pdb().ang_rot_calc
        height = self.cur_pdb().axis_rise_calc
        for c in range(1, self.cur_pdb().to_generate):
            if c == 1:
                cmd.copy('obj2_1', self.main_object)
            else:
                cmd.copy('obj2_%s'%c, 'obj2_%s'%(c-1))

            vector = [a-b for a, b in zip(p, vector_segment(p, q, height*len(p_filament)))]
            cmd.translate(vector, selection='obj2_%s'%c, camera=0, object=None)
            cmd.rotate(list(p-q), angle=angle*len(p_filament), selection='obj2_%s'%c, camera=0, object=None, origin=list(p))
        """

        p0 = np.array([float(self.cur_pdb().map_x)/2, float(self.cur_pdb().map_y)/2, 0])
        q0 = np.array([float(self.cur_pdb().map_x)/2, float(self.cur_pdb().map_y)/2, 100])

        cmd.save('Sheet_0.pdb', self.main_object + ' and sele_0') #'sele_0'
        p_filament = pdb_reader("Sheet_0.pdb").splice()

        p = get_closest_point_3d(p0, q0, np.array(p_filament[0].center()).reshape(1, 3))[0]
        q = get_closest_point_3d(p0, q0, np.array(p_filament[-1].center()).reshape(1, 3))[0]

        angle = float(self.cur_pdb().ang_rot)
        height = float(self.cur_pdb().axis_rise)

        sym = round(4.8 / float(self.cur_pdb().axis_rise))

        for c in range(1, self.cur_pdb().to_generate):
            if c == 1:
                cmd.copy('obj2_1', self.main_object)
            else:
                cmd.copy('obj2_%s'%c, 'obj2_%s'%(c-1))

            vector = [a-b for a, b in zip(p, vector_segment(p, q, height*len(p_filament)*sym))]
            cmd.translate(vector, selection='obj2_%s'%c, camera=0, object=None)
            cmd.rotate(list(p-q), angle=angle*len(p_filament)*sym, selection='obj2_%s'%c, camera=0, object=None, origin=list(p))

        return


    def extend_filaments_manual(self, to_generate):

        chains_mobile = []
        chains_target = []

        for n in range(int(self.cur_pdb().protofilament)):
            cmd.save('Sheet_0.pdb', 'sele_%s'%n)
            p_filament = pdb_reader("Sheet_0.pdb").splice()
            chains = [chain.lista[0].chain for chain in p_filament]
            print(chains)
            chains_mobile.append(chains[:-1])
            chains_target.append(chains[1:])
            residues = [res.lista[0].resi for res in p_filament[0].residues()]

        cmd.save('Sheet_0.pdb', 'sele_0')
        chain_count = len(p_filament)

        if to_generate > 1:
            cmd.copy('obj1', self.main_object)
            for n in range(chain_count):
                to_fit = []
                for x in range(int(self.cur_pdb().protofilament)):
                    for chain_n in range(len(chains_mobile[0])):
                        mobile = 'obj1 and bb. and chain %s'%chains_mobile[x][chain_n]
                        target = 'obj1 and bb. and chain %s'%chains_target[x][chain_n]
                        # for res_n in residues:
                        #     mobile = 'obj1 and name CA and chain %s and resi %s'%(chains_mobile[x][chain_n], res_n)
                        #     target = 'obj1 and name CA and chain %s and resi %s'%(chains_target[x][chain_n], res_n)

                        to_fit.append(mobile)
                        to_fit.append(target)
                cmd.pair_fit(*to_fit)


            if to_generate > 2:
                for c in range(2, to_generate):
                    cmd.copy('obj%s'%c, 'obj%s'%(c-1))
                    for n in range(chain_count):

                        to_fit = []
                        for x in range(int(self.cur_pdb().protofilament)):
                            for chain_n in range(len(chains_mobile[0])):
                                mobile = 'obj%s and bb. and chain %s'%(c, chains_mobile[x][chain_n])
                                target = 'obj%s and bb. and chain %s'%(c, chains_target[x][chain_n])
                                to_fit.append(mobile)
                                to_fit.append(target)
                        cmd.pair_fit(*to_fit)

            for n in range(int(self.cur_pdb().protofilament)):
                cmd.save('Sheet_0.pdb', 'sele_%s'%n)
                chains = [chain.lista[0].chain for chain in pdb_reader("Sheet_0.pdb").splice()]
                cmd.select('sele_%s'%n, 'none')
                cmd.select('sele_%s'%n, 'chain '+'+'.join(chains), 0, 1, 0, 1)

        color_methods(self, 'protofilaments__')


    def create_cylinder_triangles(self, start, end, radius, color=(1, 0, 0), alpha=1.0, segments=72):
        # Convert to numpy arrays
        start = np.array(start)
        end = np.array(end)
        axis = end - start
        axis_len = np.linalg.norm(axis)
        axis_dir = axis / axis_len

        # Create orthogonal vectors u and v (forming the circle plane)
        if np.allclose(axis_dir, [0, 0, 1]):
            ortho = np.array([1, 0, 0])
        else:
            ortho = np.array([0, 0, 1])
        u = np.cross(axis_dir, ortho)
        u /= np.linalg.norm(u)
        v = np.cross(axis_dir, u)

        cgo = [BEGIN, TRIANGLES, COLOR] + list(color) + [ALPHA, alpha]

        angle_step = 2 * math.pi / segments

        for i in range(segments):
            theta1 = i * angle_step
            theta2 = (i + 1) * angle_step

            # Points on start circle
            p1 = start + radius * (math.cos(theta1) * u + math.sin(theta1) * v)
            p2 = start + radius * (math.cos(theta2) * u + math.sin(theta2) * v)

            # Points on end circle
            q1 = end + radius * (math.cos(theta1) * u + math.sin(theta1) * v)
            q2 = end + radius * (math.cos(theta2) * u + math.sin(theta2) * v)

            # Triangle 1
            cgo.extend([VERTEX] + list(p1))
            cgo.extend([VERTEX] + list(p2))
            cgo.extend([VERTEX] + list(q2))

            # Triangle 2
            cgo.extend([VERTEX] + list(p1))
            cgo.extend([VERTEX] + list(q2))
            cgo.extend([VERTEX] + list(q1))

        cgo.append(END)
        return cgo


    def test0(self):
        """draw specific cyclinder"""

        """
        start = self.cur_pdb().p_cylinder['generated_axis']['start']
        end = self.cur_pdb().p_cylinder['generated_axis']['end']
        start[2] = 40
        end[2] = 60
        r = self.cur_pdb().p_cylinder['generated_axis']['radiuses'][36]
        cmd.load_cgo( [9.0] + start + end + [r, 0.6,0.6,0.6, 0.6,0.6,0.6], "cylinder37")
        end[2] = 70
        r = self.cur_pdb().p_cylinder['generated_axis']['radiuses'][16]
        cmd.load_cgo( [9.0] + start + end + [r, 0.3,0.3,0.3, 0.3,0.3,0.3], "cylinder17")
        end[2] = 90
        r = 0.5
        cmd.load_cgo( [9.0] + start + end + [r, 0,0,0, 0,0,0], "cylinder0" )
        """
        """
        start = self.cur_pdb().p_cylinder['generated_axis']['start']
        end = self.cur_pdb().p_cylinder['generated_axis']['end']
        start[2] = 50
        end[2] = 90
        r = self.cur_pdb().p_cylinder['generated_axis']['radiuses'][36]
        cmd.load_cgo( [27.0] + start + end + [r, r, 1,1,0, 1,1,0, 0,0], "cylinder37")
        end[2] = 90
        r = self.cur_pdb().p_cylinder['generated_axis']['radiuses'][26]
        cmd.load_cgo( [ 27.0] + start + end + [r, r, 1,1,0, 1,1,0, 0,0], "cylinder27")
        end[2] = 90
        r = self.cur_pdb().p_cylinder['generated_axis']['radiuses'][16]
        cmd.load_cgo( [ 27.0] + start + end + [r, r, 1,1,0, 1,1,0, 0,0], "cylinder17")
        end[2] = 90
        r = 0.5
        cmd.load_cgo( [ 27.0] + start + end + [r, r, 0,0,0, 0,0,0, 0,0], "cylinder0" )
        """

        """
        start = self.cur_pdb().p_cylinder['generated_axis']['start']
        end = self.cur_pdb().p_cylinder['generated_axis']['end']
        start[2] = 50
        end[2] = 90

        radius = self.cur_pdb().p_cylinder['generated_axis']['radiuses'][16]
        cgo = self.create_cylinder_triangles(start, end, radius, color=(1.0, 1.0, 0.0), alpha=0.7, segments=72)
        cmd.load_cgo(cgo, 'custom_cylinder16')

        radius = self.cur_pdb().p_cylinder['generated_axis']['radiuses'][36]
        cgo = self.create_cylinder_triangles(start, end, radius, color=(1.0, 1.0, 0.0), alpha=0.7, segments=72)
        cmd.load_cgo(cgo, 'custom_cylinder36')

        radius = 0.5
        cgo = self.create_cylinder_triangles(start, end, radius, color=(0.0, 0.0, 0.0), alpha=1, segments=72)
        cmd.load_cgo(cgo, 'custom_cylinder0')
        """


        start = self.cur_pdb().p_cylinder['original_axis']['start']
        end = self.cur_pdb().p_cylinder['original_axis']['end']
        start[2] = 50
        end[2] = 90

        radius = max(self.cur_pdb().p_cylinder['original_axis']['radiuses'])
        cgo = self.create_cylinder_triangles(start, end, radius, color=(0.0, 1.0, 0.0), alpha=0.7, segments=72)
        cmd.load_cgo(cgo, 'cylinder_a_original')

        cmd.pseudoatom('axis', pos=start, name='Start')
        cmd.pseudoatom('axis', pos=end, name='End')
        cmd.distance('Fa_original', 'axis and name Start', 'axis and name End')
        cmd.hide('labels', 'Fa_original')
        cmd.delete('axis')

        start = self.cur_pdb().p_cylinder['generated_axis']['start']
        end = self.cur_pdb().p_cylinder['generated_axis']['end']
        start[2] = 50
        end[2] = 90

        radius = max(self.cur_pdb().p_cylinder['generated_axis']['radiuses'])
        cgo = self.create_cylinder_triangles(start, end, radius, color=(1.0, 1.0, 0.0), alpha=0.7, segments=72)
        cmd.load_cgo(cgo, 'cylinder_a_generated')

        cmd.pseudoatom('axis', pos=start, name='Start')
        cmd.pseudoatom('axis', pos=end, name='End')
        cmd.distance('Fa_generated', 'axis and name Start', 'axis and name End')
        cmd.hide('labels', 'Fa_generated')
        cmd.delete('axis')

        start = self.pdb_list_old[self.pdb_count].p_cylinder['generated']['start']
        end = self.pdb_list_old[self.pdb_count].p_cylinder['generated']['end']


        radius = max(self.pdb_list_old[self.pdb_count].p_cylinder['generated']['radiuses'])
        cgo = self.create_cylinder_triangles(start, end, radius, color=(0.0, 0.0, 1.0), alpha=0.7, segments=72)
        cmd.load_cgo(cgo, 'cylinder_free')

        cmd.pseudoatom('axis', pos=start, name='Start')
        cmd.pseudoatom('axis', pos=end, name='End')
        cmd.distance('Fa_free', 'axis and name Start', 'axis and name End')
        cmd.hide('labels', 'Fa_free')
        cmd.delete('axis')


        start = [float(self.cur_pdb().map_x)/2, float(self.cur_pdb().map_y)/2, 50]
        end = [float(self.cur_pdb().map_x)/2, float(self.cur_pdb().map_y)/2, 90]

        radius = max(self.cur_pdb().p_cylinder['generated_axis']['radiuses'])
        cgo = self.create_cylinder_triangles(start, end, radius, color=(1.0, 0.0, 0.0), alpha=0.7, segments=72)
        cmd.load_cgo(cgo, 'cylinder_new')

        cmd.pseudoatom('axis', pos=start, name='Start')
        cmd.pseudoatom('axis', pos=end, name='End')
        cmd.distance('Fa_new', 'axis and name Start', 'axis and name End')
        cmd.hide('labels', 'Fa_new')
        cmd.delete('axis')

        """
        cmd.delete('all')

        for amy in self.pdb_list_old.lista:
            start = amy.p_cylinder['generated']['start']
            end = amy.p_cylinder['generated']['end']
            start = list(np.array(start)-np.array(end))
            print(start)
            cmd.pseudoatom('axis', pos=start, name='Start')
            cmd.pseudoatom('axis', pos=(0, 0, 0), name='End')
            cmd.distance('Fibril_axis', 'axis and name Start', 'axis and name End')
            cmd.hide('labels', 'Fibril_axis')
            cmd.delete('axis')
        """


    def test1(self):
        """Advanced test method 1."""

        #angle and rise, warp calculator

        p0 = np.array(self.cur_pdb().p_cylinder['original_axis']['start'])
        q0 = np.array(self.cur_pdb().p_cylinder['original_axis']['end'])

        cmd.save('Sheet_0.pdb', 'sele_0 and polymer')
        p_filament = pdb_reader("Sheet_0.pdb").splice()

        p = get_closest_point_3d(p0, q0, np.array(p_filament[0].center()).reshape(1, 3))[0]
        q = get_closest_point_3d(p0, q0, np.array(p_filament[-1].center()).reshape(1, 3))[0]


        #p, q = np.minimum(p0, q0), np.maximum(p0, q0)

        angles = []

        for monomer in range(int(self.cur_pdb().protofilament)):

            cmd.save('Sheet_0.pdb', 'sele_%s and polymer'%monomer)
            p_filament = pdb_reader("Sheet_0.pdb").splice()

            for resi in self.cur_pdb().get_chain('pf', monomer).resi:
                for n in range(len(p_filament)-1):
                    a0 = p_filament[n].select('resi', [str(resi)]).select('name', ['CA']).lista[0].xyz()
                    a3 = p_filament[n+1].select('resi', [str(resi)]).select('name', ['CA']).lista[0].xyz()

                    angles.append(get_dihedral(a0, p, q, a3))
        #print(angles)
        angle = np.degrees(sum(angles) /len(angles))
        self.cur_pdb().ang_rot_calc = round(angle, 4)
        print('original: ', self.cur_pdb().ang_rot)
        print('original: ', self.cur_pdb().axis_sym)
        print('measured: ', angle)

        #rise
        heights = []
        for monomer in range(int(self.cur_pdb().protofilament)):

            cmd.save('Sheet_0.pdb', 'sele_%s and polymer'%monomer)
            p_filament = pdb_reader("Sheet_0.pdb").splice()

            for n in range(len(p_filament)-1):
                a1 = np.array(p_filament[n].select('name', ['CA']).center()).reshape(1, 3)
                a2 = np.array(p_filament[n+1].select('name', ['CA']).center()).reshape(1, 3)

                h = get_closest_point_3d(p, q, a1)-get_closest_point_3d(p, q, a2)
                heights.append(np.linalg.norm(h))

        height = sum(heights) /len(heights)
        print('original: ', self.cur_pdb().axis_rise)
        print('measured: ', height)

        self.cur_pdb().axis_rise_calc = round(height, 3)

        #warping
        h_chain = []
        for monomer in range(int(self.cur_pdb().protofilament)):

            cmd.save('Sheet_0.pdb', 'sele_%s and polymer'%monomer)
            p_filament = pdb_reader("Sheet_0.pdb").splice()

            for chain in range(len(p_filament)):

                #only for z!
                main_axis = 'z'
                h_chain_min = min(p_filament[0].select('name', ['CA']).data(main_axis))
                h_chain_max = max(p_filament[0].select('name', ['CA']).data(main_axis))
                h_chain.append(abs(h_chain_max-h_chain_min))
        print(h_chain)

        self.cur_pdb().warping = round(sum(h_chain)/len(h_chain), 3)
        print(self.cur_pdb().warping)
        return


    @suspend_decorator
    def test2(self):
        #test2_cluster_calc
        #cluster_calc

        cmd.delete('all')
        cmd.cd('..')

        #prot = 'TAU'
        #prots = ['B2M', 'ISLET', 'LLC', 'PRION', 'SA2', 'TDP43', 'TMP106', 'TTR']
        #prots = ['TTR']
        #prots = [ 'ABETA', 'ASYN', 'B2M', 'ISLET', 'LLC', 'PRION', 'SA2', 'TAU', 'TDP43', 'TMP106', 'TTR']
        prots = ['TDP43']

        for prot in prots:
            try:
                data = [] #(x_data, y_data, rmsd, code_pf-code_pf)
                data_names = []
                x_data = 0
                y_data = 0

                for n, amy1 in enumerate(self.pdb_list.lista):
                    if amy1.prot_abr != prot:
                        continue

                    for pf1 in range(int(amy1.protofilament)):
                        cmd.load('%s/monomer_%s.pdb'%(amy1.code, pf1))
                        name1 = '%s_%s'%(amy1.code, pf1)
                        data_names.append(name1)
                        print(name1)
                        cmd.set_name('monomer_%s'%pf1, name1)
                        model1 = pdb_reader('%s/monomer_%s.pdb'%(amy1.code, pf1)).splice()
                        model1.pop(int(len(model1)/2))
                        select_model_chain(model1, 'to_delete')
                        cmd.remove('to_delete')
                        self.suspend_refresh(wait=0.01)
                        for amy2 in self.pdb_list.lista:
                            if amy2.prot_abr != prot:
                                continue

                            for pf2 in range(int(amy2.protofilament)):
                                name2 = '%s_%s'%(amy2.code, pf2)
                                #print(name2)
                                if name1 == name2:
                                    y_data += 1
                                    continue
                                cmd.load('%s/monomer_%s.pdb'%(amy2.code, pf2))

                                cmd.set_name('monomer_%s'%pf2, name2)
                                model2 = pdb_reader('%s/monomer_%s.pdb'%(amy2.code, pf2)).splice()
                                model2.pop(int(len(model2)/2))
                                select_model_chain(model2, 'to_delete')
                                cmd.remove('to_delete')

                                rmsd = cmd.align(name1, name2, cycles=0)
                                #rmsd = cmd.super(name1, name2, cycles=0)
                                #rmsd = cmd.super(name1, name2)

                                self.suspend_refresh(wait=0.001)
                                print(rmsd)
                                data.append((x_data, y_data, rmsd[0], '%s-%s'%(name1, name2), rmsd[2]))
                                cmd.delete(name2)
                                y_data += 1

                        cmd.delete(name1)
                        x_data += 1
                        y_data = 0

                #print(x_data)
                #print(data)
                for line in data:
                    #print(line)
                    pass

                with open('%s_data_nofit.json'%prot, 'w+') as f:
                    json_txt = json.dumps(data, default=lambda obj: obj.__dict__, indent=4)
                    f.write(json_txt)
                with open('%s_name_nofit.json'%prot, 'w+') as f:
                    json_txt = json.dumps(data_names, default=lambda obj: obj.__dict__, indent=4)
                    f.write(json_txt)
            except Exception:
                #traceback.print_exc()
                with open('%s_data_nofit.json'%prot, 'w+') as f:
                    f.write('error')
                continue

    @suspend_decorator
    def test2_cluster_calc_all(self):
        #test2_cluster_calc_all
        #cluster_calc

        cmd.delete('all')
        cmd.cd('..')

        try:
            data = [] #(x_data, y_data, rmsd, code_pf-code_pf)
            data_names = []
            x_data = 0
            y_data = 0

            for n, amy1 in enumerate(self.pdb_list.lista):

                for pf1 in range(int(amy1.protofilament)):
                    cmd.load('%s/monomer_%s.pdb'%(amy1.code, pf1))
                    name1 = '%s_%s'%(amy1.code, pf1)
                    data_names.append(name1)
                    print(name1)
                    cmd.set_name('monomer_%s'%pf1, name1)
                    model1 = pdb_reader('%s/monomer_%s.pdb'%(amy1.code, pf1)).splice()
                    model1.pop(int(len(model1)/2))
                    select_model_chain(model1, 'to_delete')
                    cmd.remove('to_delete')
                    self.suspend_refresh(wait=0.01)
                    for amy2 in self.pdb_list.lista:

                        for pf2 in range(int(amy2.protofilament)):
                            name2 = '%s_%s'%(amy2.code, pf2)
                            #print(name2)
                            if name1 == name2:
                                y_data += 1
                                continue
                            cmd.load('%s/monomer_%s.pdb'%(amy2.code, pf2))

                            cmd.set_name('monomer_%s'%pf2, name2)
                            model2 = pdb_reader('%s/monomer_%s.pdb'%(amy2.code, pf2)).splice()
                            model2.pop(int(len(model2)/2))
                            select_model_chain(model2, 'to_delete')
                            cmd.remove('to_delete')

                            rmsd = cmd.align(name1, name2, cycles=0)
                            #rmsd = cmd.super(name1, name2, cycles=0)
                            #rmsd = cmd.super(name1, name2)

                            self.suspend_refresh(wait=0.001)
                            print(rmsd)
                            data.append((x_data, y_data, rmsd[0], '%s-%s'%(name1, name2), rmsd[2]))
                            cmd.delete(name2)
                            y_data += 1

                    cmd.delete(name1)
                    x_data += 1
                    y_data = 0

            #print(x_data)
            #print(data)

            with open('all_data_nofit.json', 'w+') as f:
                json_txt = json.dumps(data, default=lambda obj: obj.__dict__, indent=4)
                f.write(json_txt)
            with open('all_name_nofit.json', 'w+') as f:
                json_txt = json.dumps(data_names, default=lambda obj: obj.__dict__, indent=4)
                f.write(json_txt)
        except Exception:
            #traceback.print_exc()
            with open('all_data_nofit.json', 'w+') as f:
                f.write('error')


    @suspend_decorator
    def test2_cluster_calc_cryo(self):
        #test2_cluster_calc_cryo
        #cluster_calc
        # timeresoved only

        cmd.delete('all')
        cmd.cd('..')

        #prot = 'TAU'
        prots = ['timeresolved']
        codes = ['8Q9M', '8QJJ', '8Q9K', '8Q9L', '8Q9O', '8Q88', '8Q9F', '8Q8E', '8Q9H', '8Q8F', '8Q9G', '8Q8D', '8Q9I', '8Q8L', '8Q9J', '8QCP', '8Q27', '8Q8U', '8Q8C', '8Q99', '8Q7T', '8Q9A', '8QCR', '8Q9B', '8Q9C', '8Q9D', '8Q7P', '8Q9E', '8Q2J', '8Q2K', '8Q8V', '8Q2L', '8Q8W', '8Q7F', '8Q8X', '8Q7L', '8Q8Y', '8Q7M', '8Q8Z', '8Q98', '8Q97', '8Q8R', '8Q8M', '8Q8S']
        method = '_nofit'

        for prot in prots:
            try:
                data = [] #(x_data, y_data, rmsd, code_pf-code_pf)
                data_names = []
                x_data = 0
                y_data = 0

                for n, amy1 in enumerate(self.pdb_list.lista):
                    if amy1.code not in codes:
                        continue

                    for pf1 in range(int(amy1.protofilament)):
                        cmd.load('%s/monomer_%s.pdb'%(amy1.code, pf1))
                        name1 = '%s_%s'%( amy1.code, pf1)
                        data_names.append(amy1.note+'_'+name1)
                        print(1, name1)
                        cmd.set_name('monomer_%s'%pf1, name1)
                        model1 = pdb_reader('%s/monomer_%s.pdb'%(amy1.code, pf1)).splice()
                        model1.pop(int(len(model1)/2))
                        select_model_chain(model1, 'to_delete')
                        cmd.remove('to_delete')
                        self.suspend_refresh(wait=0.01)
                        for amy2 in self.pdb_list.lista:
                            if amy2.code not in codes:
                                continue

                            for pf2 in range(int(amy2.protofilament)):
                                name2 = '%s_%s'%(amy2.code, pf2)
                                #print(2, name2)
                                if name1 == name2:
                                    y_data += 1
                                    continue
                                cmd.load('%s/monomer_%s.pdb'%(amy2.code, pf2))

                                cmd.set_name('monomer_%s'%pf2, name2)
                                model2 = pdb_reader('%s/monomer_%s.pdb'%(amy2.code, pf2)).splice()
                                model2.pop(int(len(model2)/2))
                                select_model_chain(model2, 'to_delete')
                                cmd.remove('to_delete')

                                rmsd = cmd.align(name1, name2, cycles=0)

                                self.suspend_refresh(wait=0.01)
                                #print(rmsd)
                                data.append((x_data, y_data, rmsd[0], '%s-%s'%(name1, name2)))
                                cmd.delete(name2)
                                y_data += 1

                        cmd.delete(name1)
                        x_data += 1
                        y_data = 0

                with open('%s_data%s_time.json'%(prot, method), 'w+') as f:
                    json_txt = json.dumps(data, default=lambda obj: obj.__dict__, indent=4)
                    f.write(json_txt)
                with open('%s_name%s_time.json'%(prot, method), 'w+') as f:
                    json_txt = json.dumps(data_names, default=lambda obj: obj.__dict__, indent=4)
                    f.write(json_txt)
            except Exception:
                #traceback.print_exc()
                with open('%s_data_nofit_time.json'%prot, 'w+') as f:
                    f.write('error')
                continue


    @suspend_decorator
    def test2_cluster_watch(self):
        #test2_cluster_watch
        #watch clustering result

        with open('TAU_clusters.json', "r") as f:
            clusters = json.load(f)

        for cluster in clusters:

            for amy1 in clusters[cluster]:
                cmd.delete('all')
                code, pf1 = amy1.split('_')
                cmd.load('%s/monomer_%s.pdb'%(code, pf1))
                print(amy1)
                cmd.set_name('monomer_%s'%pf1, amy1)
                select_model_chain(pdb_reader('%s/monomer_%s.pdb'%(code, pf1)).splice()[3:], 'to_delete')
                cmd.remove('to_delete')
                self.suspend_refresh(wait=0.3)

            self.suspend_refresh(wait=2)


    def testx(self):
        """Advanced test method 2."""


        """
        # find axis of rotation with cicrle fitting
        for monomer in range(int(self.cur_pdb().protofilament)):

            cmd.save('Sheet_0.pdb', 'sele_%s'%monomer)
            p_filament = pdb_reader("Sheet_0.pdb").splice()

            diff = []
            for axis in ['x', 'y', 'z']:
                coords = []
                for chain in p_filament:
                    coords.append(chain.median(axis))

                diff.append(max(coords)-min(coords))

            main_axis=['x', 'y', 'z'][diff.index(max(diff))]

            centers = []

            for resi in self.cur_pdb().get_chain('pf', monomer).resi:
                points = []

                #show atoms
                test = mAtom()
                for n in range(len(p_filament)):
                    test.add(p_filament[n].select('resi', [str(resi)]).select('name', ['CA']).lista[0])
                select_model(test, 'test')

                #circle
                for n in range(len(p_filament)):
                    xyz = p_filament[n].select('resi', [str(resi)]).select('name', ['CA']).lista[0].xyz()
                    xyz.pop({'x':0, 'y':1, 'z':2}[main_axis])
                    points.append(xyz)
                #print(points)
                print(get_lc(points[0], points[-1]))
                if get_lc(points[0], points[-1]) < 3:
                    print('too close')
                    continue

                x, y, r, e = cf.least_squares_circle(points)
                if r > 200:
                    continue
                #self.cgoCircle(x, y, p_filament[1].median(main_axis), r)
                cmd.pseudoatom('tmp1', pos=(x, y, p_filament[1].median(main_axis)))
                centers.append((x, y))

        print('done')

        centrum = [mean([c[0] for c in centers]), mean([c[1] for c in centers]), p_filament[1].median(main_axis)]
        cmd.pseudoatom('tmp2', pos=centrum)
        cmd.show('sphere', 'tmp2')
        cmd.color('green', 'tmp2')
        centrum = [median([c[0] for c in centers]), median([c[1] for c in centers]), p_filament[1].median(main_axis)]
        cmd.pseudoatom('tmp3', pos=centrum)
        cmd.show('sphere', 'tmp3')
        cmd.color('yellow', 'tmp3')

        """


        """
        #? some kind of manual extend?
        chains_mobile = []
        chains_target = []

        for n in range(int(self.cur_pdb().protofilament)):
            cmd.save('Sheet_0.pdb', 'sele_%s'%n)
            p_filament = pdb_reader("Sheet_0.pdb").splice()
            chains = [chain.lista[0].chain for chain in p_filament]
            print(chains)
            chains_mobile.append(chains[:-1])
            chains_target.append(chains[1:])

        print(chains_mobile)
        print(chains_target)
        cmd.save('Sheet_0.pdb', 'sele_0')
        chain_count = len(p_filament)

        try:
            cmd.select('new', 'objx')
        except:
            cmd.copy('objx', self.main_object)

        for n in range(chain_count):
            to_fit = []
            for x in range(int(self.cur_pdb().protofilament)):
                for chain_n in range(len(chains_mobile[0])):
                    mobile = 'objx and bb. and chain %s'%chains_mobile[x][chain_n]
                    target = 'objx and bb. and chain %s'%chains_target[x][chain_n]
                    to_fit.append(mobile)
                    to_fit.append(target)
            cmd.pair_fit(*to_fit)
            break
        """

        """
        cmd.set('line_width', '2.5')
        cmd.set('nb_spheres_size', '0.350')
        cmd.set('stick_radius', '0.350')
        """
        pass


    def test3(self, *args):
        """Advanced test method 3."""

        #chain data creator new_entry_analyser
        cmd.remove('not polymer.protein')
        for n in range(int(self.cur_pdb().protofilament)):
            chains = pdb_reader('monomer_%s.pdb'%n).splice()
            for chain in chains:

                chain_data = amy_chain()
                chain_data.pf = n
                chain_data.letter = chain.lista[0].chain
                x, y, z = chains[0].lista[0].xyz()
                cmd.select('sele', 'none')
                cmd.select('sele', 'bychain sele or ( X>%s and X<%s and Y>%s and Y<%s and Z>%s and Z<%s)'%(x-0.005,x+0.005,y-0.005,y+0.005,z-0.005,z+0.005), 0, 1, 0, 1)
                chain_data.resi = get_data('(bychain sele) and name O ', 'resv')
                chain_data.seq = get_data('(bychain sele) and name O ', 'oneletter')

                self.cur_pdb().add_chain(chain_data)

        print('done')


    def test4(self, *args):
        """Advanced test method 4."""
        #warping old!


        v = np.array(self.cur_pdb().p_cylinder['generated']['start']) - np.array(self.cur_pdb().p_cylinder['generated']['end'])
        v = v / np.sqrt(np.sum(v ** 2)) #norm
        angle = - math.degrees(math.acos(np.dot((0, 1, 0), v))) #angle from dot
        v_cross = list(np.cross((0, 1, 0), v)) #axis from cross

        #limit to skip
        if angle > -0.5 or angle < -179.5:
            pass
        else:
            cmd.rotate(v_cross, angle=angle, selection='all', camera=0, object=None, origin=[0, 0, 0])

        cmd.save('Sheet_0.pdb', 'sele_0')
        p_filament = pdb_reader("Sheet_0.pdb").splice()
        chain_count = len(p_filament)
        diff = []
        for axis in ['x', 'y', 'z']:
            coords = []
            for chain in p_filament:
                coords.append(chain.median(axis))

            diff.append(max(coords)-min(coords))
        main_axis=['x', 'y', 'z'][diff.index(max(diff))]
        print(main_axis)

        h_chain = []

        for n in range(int(self.cur_pdb().protofilament)):
            cmd.save('Sheet_0.pdb', 'sele_%s'%n)
            p_filament = pdb_reader("Sheet_0.pdb").splice()
            h_chain_min = min(p_filament[0].select('name', ['CA']).data(main_axis))
            h_chain_max = max(p_filament[0].select('name', ['CA']).data(main_axis))
            print(p_filament[0].select('name', ['CA']).data(main_axis))
            print('min', p_filament[0].select('name', ['CA']).data(main_axis).index(h_chain_min))
            print('max', p_filament[0].select('name', ['CA']).data(main_axis).index(h_chain_max))
            h_chain.append(abs(h_chain_max-h_chain_min))
        print(h_chain)
        print(sum(h_chain)/len(h_chain))
        return sum(h_chain)/len(h_chain)


    def test5(self, *args):
        """Advanced test method 5."""
        #max radius

        start = np.array(self.cur_pdb().p_cylinder['original_axis']['start'])
        end = np.array(self.cur_pdb().p_cylinder['original_axis']['end'])

        maxr = []
        for n in range(int(self.cur_pdb().protofilament)):
            cmd.save('Sheet_0.pdb', 'sele_%s'%n)
            p_filament = pdb_reader("Sheet_0.pdb").splice()
            rs = p_filament[0].coord_list()
            rs = np.reshape(rs, (-1, 3))
            maxr.append(np.max(get_dist_3d_line(start, end, rs)))

        return max(maxr)


    def reloader(self):
        """Advanced test method 3."""
        cmd.set_wizard()

        reload(sys.modules['acwf'])
        print('reloaded')


    def re_calc_method(self, only_this=False):
        """Advanced method for recalculating the database."""

        try:
            while self.pdb_count < len(self.pdb_list.lista):


                self.busy = 0
                try:
                    self.next_peptide(x=self.pdb_count)
                    print(self.cur_pdb().code)

                    #self.find_connections()
                    #self.inteligent_window()
                    """
                    try:
                        self.cur_pdb().error_during_auto
                    except:
                        continue
                    """

                    #self.setup() #its inside next!
                    #new_entry_analyser(self)

                    """
                    if self.cur_pdb().rmsd != self.pdb_list_old[self.pdb_count].rmsd:
                        print('old: ', self.pdb_list_old[self.pdb_count].rmsd)
                        print('new: ', self.cur_pdb().rmsd)

                        self.find_connections()
                        self.inteligent_window()
                    else:
                        print('no change')
                    """

                    #self.inteligent_window_recalc()
                    #self.test3()
                    #self.test()

                    
                    #rmax, warp, angle, rise
                    rmax = self.test5()
                    self.cur_pdb().p_cylinder['original_axis']['rmax'] = rmax
                    print(rmax)

                    self.test1()
                    print(self.cur_pdb().code)
                    print(self.pdb_count)
                    
                    
                    #main axis (z)
                    cmd.save('Sheet_0.pdb', 'sele_0 and polymer')
                    p_filament = pdb_reader("Sheet_0.pdb").splice()
                    diff = []
                    for axis in ['x', 'y', 'z']:
                        coords = []
                        for chain in p_filament:
                            coords.append(chain.median(axis))

                        diff.append(max(coords)-min(coords))
                    main_axis=['x', 'y', 'z'][diff.index(max(diff))]
                    self.cur_pdb().main_axis = main_axis
                    
                    pass
                except Exception:
                    #traceback.print_exc()
                    self.cur_pdb().error_during_auto = 'yes'
                finally:
                    self.suspend_refresh()
                    self.reset()
                    self.pdb_count += 1
                if only_this:
                    break

        finally:
            print("\a")
            cmd.set('suspend_updates', 0, quiet=1)


    def view(self):
        """Reorients the pymol viewer."""
        cmd.orient('all')

        my_view = list(cmd.get_view())
        views = {'x':[0, 0, 1, 1, 0, 0, 0, 1, 0],
                 'y':[0, 1, 0, 0, 0, 1, 1, 0, 0],
                 'z':[1, 0, 0, 0, 1, 0, 0, 0, 1],
                 0:[0, 0, 1, 1, 0, 0, 0, 1, 0]}

        for n in range(9):
            my_view[n] = views[self.main_axis][n]
        cmd.set_view(tuple(my_view))

        if self.view_count == 1:
            cmd.turn('z', 90)
        if self.view_count == 2:
            cmd.turn('x', 90)
        if self.view_count == 3:
            cmd.turn('y', 90)
            self.view_count = 0
            return

        self.view_count += 1


    def color(self, method='sheet'):
        # recolor

        color_methods(self, method)
        self.color_last_method = method
        if self.labels == 'descriptor':
            self.toggle_labels_descriptor()
            self.toggle_labels_descriptor()


    def draw_axis_method(self):
        """Draws the coordinate axis."""
        if self.axis == 'ON':
            for n in ['X_axis', 'Y_axis', 'Z_axis', 'axis']:
                cmd.delete(n)
            self.axis = 'OFF'
        else:
            draw_axis()
            self.axis = 'ON'


    def manual_mode(self):
        """Turns on or of the manual mode."""
        if self.picking_mode == 'OFF':

            if self.path:

                cmd.set('mouse_selection_mode', 1)
                cmd.util.cbaw('all')
                cmd.util.cbag('?sele_0')
                cmd.delete('?lengths')
                self.picking_mode = 'connected'

            else:
                print('CCP4 is not available, calculations cannot be performed')
                self.error_message = 'CCP4 is not available, calculations cannot be performed'
        else:
            self.pick_count = 0
            color_methods(self, 'sheets')
            cmd.delete('?lengths')
            self.picking_mode = 'OFF'

        cmd.refresh_wizard()


    def launcher(self, thread, name, *args, **kwargs):
        """Launch the other motheds with single or multithreads."""

        try:
            self.error_message = ''
            self_fn = getattr(self, name)

            if thread == 'multi':
                t = threading.Thread(target=self_fn, args=args, kwargs=kwargs)
                t.setDaemon(1)
                t.start()

            elif thread == 'single':
                self_fn(*args, **kwargs)

            cmd.refresh_wizard()
        
        except Exception:
            cmd.set('suspend_updates', 0, quiet=1)
            print('Unknown error occured')
            self.error_message = 'Unknown error occured'


    def menu_update(self):
        """Updates the menu buttons."""

        for s in [self.sele_prot, self.sele_protofilament]:
            s.sort()

        selection = [[2, 'List selection:', ''], [1, 'Protein: ' + ','.join(self.sele_prot), []]]
        for t in self.sele_prot_types:
            selection[-1][-1].append([1, t, 'cmd.get_wizard().launcher("single","toggle_value","sele_prot",["%s"])'%t])

        selection.append([1, 'Protofilament: ' + ', '.join(self.sele_protofilament), []])
        for t in self.sele_protofilament_types:
            selection[-1][-1].append([1, t, 'cmd.get_wizard().launcher("single","toggle_value","sele_protofilament",["%s"])'%t])
        """
        selection.append([1, 'Main cluster: ' + ', '.join(self.sele_mainclust), []])
        for t in self.sele_mainclust_types:
            selection[-1][-1].append([1, str(t), 'cmd.get_wizard().launcher("single","toggle_value","sele_mainclust",["%s"])'%t])

        selection.append([1, 'Secondary cluster: ' + ', '.join(self.sele_secclust), []])
        for t in self.sele_secclust_types:
            selection[-1][-1].append([1, str(t), 'cmd.get_wizard().launcher("single","toggle_value","sele_secclust",["%s"])'%t])
        """
        selection.append([1, 'Reset criteria ', 'cmd.get_wizard().launcher("single","reset_selection")'])
        self.menu['selection'] = selection

        # coloring methods
        self.f_coloring_types = ['even-odd__', 'recolor-sc-6', 'recolor-area-6', 'recolor-sdi-6', 'protofilaments__']
        colors = [[ 2, 'Coloring', '' ]]
        for color in self.f_coloring_types:
            colors.append([1, color[:-2], 'cmd.get_wizard().color("%s")'%color])
        color_range_setting = {'absolute':'relative', 'relative':'absolute'}[self.color_range]
        colors.append([1, 'Color range: %s'%self.color_range, 'cmd.get_wizard().set_value("color_range","%s","coloring")'%color_range_setting])
        color_inter_setting = {'inter':'intra', 'intra':'inter'}[self.color_inter]
        colors.append([1, 'Color inter/intra: %s'%self.color_inter, 'cmd.get_wizard().set_value("color_inter","%s","coloring")'%color_inter_setting])
        self.menu['Colors'] = colors

        self.menu['Advanced'] = [[2, 'Advanced', ''],
                                 [1, 'test0', 'cmd.get_wizard().launcher("multi","test0")'],
                                 [1, 'test1', 'cmd.get_wizard().launcher("multi","test")'],
                                 [1, 'test2', 'cmd.get_wizard().launcher("multi","test2")'],
                                 [1, 'test3', 'cmd.get_wizard().launcher("multi","test3")'],
                                 [1, 'test4', 'cmd.get_wizard().launcher("multi","test4")'],
                                 [1, 'test5', 'cmd.get_wizard().launcher("multi","test5")'],
                                 [1, 'Picking: Reset chains', 'cmd.get_wizard().set_value("picking_mode","reset_chains","wizard")'],
                                 [1, 'Picking: Delete', 'cmd.get_wizard().set_value("picking_mode","Delete","wizard")'],
                                 [1, 'Picking: Manual connected', 'cmd.get_wizard().set_value("picking_mode","connected","wizard")'],
                                 [1, 'Moving window', 'cmd.get_wizard().launcher("multi","moving_window","%s")'%self.main_object],
                                 [1, 'Recalculation of database', 'cmd.get_wizard().launcher("multi","re_calc_method")'],
                                 [1, 'Recalculation of this', 'cmd.get_wizard().launcher("multi","re_calc_method",True)']
                                 ]

        self.menu['Toggle'] = [[2, 'Toggle labels and axes:', ''],
                                 [1, 'Residue descriptors', 'cmd.get_wizard().launcher("single","toggle_labels_descriptor")'],
                                 [1, 'Residue ID', 'cmd.get_wizard().launcher("single","toggle_labels_residue")'],
                                 [1, 'Chain ID', 'cmd.get_wizard().launcher("single","toggle_labels_chain")'],
                                 [1, 'Fibril axis', 'cmd.get_wizard().launcher("single","toggle_fibril_axis")'],
                                 [1, 'Cartesian axes', 'cmd.get_wizard().draw_axis_method()'],
                                 ]
        if self.mode == 'evaluate':
            self.menu['Toggle'].pop(4)
        self.menu['Settings'] = [[2, 'Settings:', ''],
                                 [1, 'Show top panel: %s'%self.need_prompt, 'cmd.get_wizard().def_changer("need_prompt")'],
                                 [1, 'Generate chains: %s'%self.generate_chains, 'cmd.get_wizard().def_changer("generate_chains")'],
                                 ]

        cmd.refresh_wizard()


    def get_panel(self):
        """Configures the buttons on rigth panel."""


        if self.pdb_count == -1 and self.mode == 'database':
            view_interface = 'View interfaces (0/0)'
        elif len(self.cur_pdb().interface) == 0:
            view_interface = 'View interfaces (0/0)'
        else:
            view_interface = 'View interfaces (%s/%s)'%(self.interface_count, len(self.cur_pdb().interface))
        if self.picking_mode == 'OFF':
            pick_mode = 'ON'
        else:
            pick_mode = 'OFF'

        if self.mode == 'database':

            try:
                connected_section_list = []
                with open('connected_section.txt', 'r') as f:
                    data = f.read().splitlines()

                    for n in data:
                        connected_section_list.append(n.split("_"))
                con_sec = '(%s/%s)'%(self.pick_count_cs, len(connected_section_list))
            except:
                con_sec = ''


            panel = [
                [1, 'ACWF Wizard: Database', ''],
                [3, 'Selection criteria', 'selection'],
                [2, 'Next structure', 'cmd.get_wizard().launcher("single","next_peptide")'],
                [2, 'Previous structure', 'cmd.get_wizard().launcher("single","next_peptide",-1)'],

                [1, 'Clusters', ''],
                [2, 'Next in cluster', 'cmd.get_wizard().launcher("single","next_peptide_cluster",1)'],
                [2, 'Prev in cluster', 'cmd.get_wizard().launcher("single","next_peptide_cluster",-1)'],

                [1, 'Viewing:', ''],
                [3, 'Coloring', 'Colors'],
                [2, 'Export image', 'cmd.get_wizard().launcher("single","exp_image")'],
                [3, 'Toggle labels and axes', 'Toggle'],
                [2, 'Reset/rotate view', 'cmd.get_wizard().launcher("single","view")'],

                [1, 'Menu:', ''],
                [2, 'Change to ACWF Evaluation', 'cmd.get_wizard().launcher("single","mode_change","evaluate")'],
                [3, 'Settings', 'Settings'],
                [2, 'Quit', 'cmd.set_wizard()']
                ]

            try:
                if self.cur_pdb().source == 'personal':
                    pass
                    # panel.insert(11, [2, 'Delete from database', 'cmd.get_wizard().launcher("single","delete_from_database")'])
            except:
                pass

            if self.overwrite_state:
                panel.insert(0, [1, 'This deletes its folder,', 1])
                panel.insert(1, [1, 'are you sure?', 1])
                panel.insert(2, [2, 'Yes', 'cmd.get_wizard().launcher("single","delete_from_database",1)'])
                panel.insert(3, [2, 'No', 'cmd.get_wizard().launcher("single","delete_from_database",0)'])

        elif self.mode == 'evaluate':

            panel = [
                [1, 'ACWF Wizard: Evaluation', 1],
                [2, '1: Assign protofilaments', 'cmd.get_wizard().launcher("single","reset_monomers")'],
                [2, '2: Setup fibrils', 'cmd.get_wizard().launcher("single","setup")'],
                [2, '3: Assign local interfaces', 'cmd.get_wizard().launcher("multi","find_connections")'],
                [2, '4: Sliding window', 'cmd.get_wizard().launcher("multi","inteligent_window")'],
                [2, '5: Reset', 'cmd.get_wizard().launcher("single","delete_all",["everything","interface","last_amy"])'],

                [1, 'Viewing:', ''],
                #[3, 'Coloring', 'Colors'],
                [2, 'Export image', 'cmd.get_wizard().launcher("single","exp_image")'],
                [3, 'Toggle labels and axes', 'Toggle'],
                [2, 'Reset/rotate view', 'cmd.get_wizard().launcher("single","view")'],
                
                [1, 'Menu:', ''],
                [2, 'Change to ACWF Database', 'cmd.get_wizard().launcher("single","mode_change","database")'],
                [3, 'Settings', 'Settings'],
                [2, 'Quit', 'cmd.set_wizard()']
                ]

            if self.calc_state[0] == 4:
                panel.insert(7, [3, 'Coloring', 'Colors'],)

            if self.overwrite_state:
                panel.insert(0, [1, 'Overwrite existing data?', 1])
                panel.insert(1, [2, 'Yes', 'cmd.get_wizard().launcher("single","overwrite",1)'])
                panel.insert(2, [2, 'No', 'cmd.get_wizard().launcher("single","overwrite",0)'])

        if self.advanced:
            panel.append([3, 'Advanced', 'Advanced'])
            panel.append([2, 'Reloader', 'cmd.get_wizard().launcher("multi","reloader")'])

        return panel


def acw_starter():
    try:
        cmd.set_wizard(AcWizardFibrilMultiSheet(acw_database_path))
    except Exception:
        cmd.set_wizard()

cmd.extend('acwf', acw_starter)


def acwf_goto(goto=0):
    try:
        cmd.get_wizard().cmd.get_wizard().launcher('single', 'next_peptide', goto=goto)
    except AttributeError:
        cmd.set_wizard(AcWizardFibrilMultiSheet(acw_database_path))
        cmd.get_wizard().cmd.get_wizard().launcher('single', 'next_peptide', goto=goto)

cmd.extend('acwfgoto', acwf_goto)

print("Type 'acwf' to start the Amyloid Coordinate Wizard")
cmd.set_wizard(AcWizardFibrilMultiSheet(acw_database_path))
