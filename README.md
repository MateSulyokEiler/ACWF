# ACWF – Amyloid coordinate wizard fibril

ACWF is a Pymol plugin for finding local interfaces in amyloid fibrils, and calculating their descriptors: shape complementarity, buried surface area and surface detail index. It can be used on fibril structurres measured by cryo-EM. The plugin also serves as a database and a viewer for the amyloid fibril structures downloaded from the [PDB](https://www.rcsb.org/). 

## Requirements:
The plugin works with recent Windows and Linux [Pymol](https://www.pymol.org/) editions with python version 2.7 and above. The calculation methods require installation of the [CCP4](https://www.ccp4.ac.uk/) program package (including programs AREAIMOL and SC). 

## Installation:
1. Download or clone this repository. 
The folder contains the plugin file (acw.py and its module folder); a configuration file (config.txt); and the database of available PDB structures of amyloid fibril structures (set of folders and acw_database.json file). The database should be in the same folder to use the database mode.
2. (Optional) Edit the path of CCP4 installation folder in the configuration file (config.txt, inserting the path after the equals sign). More information about system specific syntax is in the config file. Note, for using only the “database mode” (if the plugin will be used only as a database and viewer), no need for pre-installation of CCP4 nor editing the config file.
3. The plugin can be started in PyMol by loading the acwf.py file by selecting ‘file/run script…’ or ‘file/open…’ then typing 'acwf' in the console.

## User guide:
Information about using the wizard can be found in the UserGuide.pdf.

## Reference:
Comming soon.
Please cite our previous work: [Unravelling the Complexity of Amyloid Peptide Core Interfaces](https://pubs.acs.org/doi/10.1021/acs.jcim.4c01479).
