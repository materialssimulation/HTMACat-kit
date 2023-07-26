import requests
import os
import shutil
from os.path import join

def runCRN_net(configfname='ReactInfo'):
    with open(configfname) as f:
        config = f.read()
    post_dict = {'submit':'api', 'ReactInfo':config}
    url = 'http://www.catalysthub.net/CRN-api.php'
    r = requests.post(url, data=post_dict, verify=False)
    if r.status_code != 200:
        print('ERROR: status not 200 but', r.status_code)
        return
    else:
        print('Success...')
        return r.text

"""
Created on July 26 20:37 2023.

@author: RongxinChen
"""

def run_crnconfiggen(logname="CRNGenerator_log.txt"):
    """
        Generate the adsorption configurations for all intermediate species and 
        polar subgraph elementary steps based on CRN log file. 
    """
    if os.path.exists(join(os.path.abspath('.'), 'stems')) or os.path.exists(join(os.path.abspath('.'), 'all_species')):
        raise FileExistsError("The directories 'stems' and 'all_species' already exist in the current directory!\
                              Please delete these directories and try again!")

    smiles_list = crnlog_readspecies(logname)
    create_input(smiles_list)
    
    stem_list = crnlog_getstems(logname)
    stem_value = crnlog_readstems(stem_list)
    create_stem_dir(stem_value)

    if os.path.exists(join(join(join(os.path.abspath('.'), 'stems'), '0_stem'), '0_step')):
        print('Success...')
    else:
        print('Generation failed, please recheck input information!')

def crnlog_readspecies(inputfile):
    '''
    Return a str list for storing SMILES expressions for all reacting species.
    '''
    smiles_list = []
    index_line = 0
    index_smiles = 0
    
    with open(inputfile, 'r') as f:
        f_lines = f.readlines()
        start_line = 999
        for line in f_lines:
            if line[:26] == 'All possible species found':
                # Cerify start_line
                start_line = index_line+2
            if index_line >= start_line:
                if not line[0].isdigit():
                    break
                else:
                    smiles_name = get_smilesname(line)
                    smiles_list.append(smiles_name)
            index_line += 1            
    
    return smiles_list
    
def get_smilesname(line):
    '''
    Return the SMILES expression from one line in log file
    '''
    judge_smlies = 0
    init_index = 0
    final_index = 0
    strjudge_list = [' ', '(', ')']
    
    for i in range(len(line)):
        if line[i] == ' ' and  (not (line[i+1] in strjudge_list)):
            judge_smlies += 1
            
        if judge_smlies == 2 and init_index == 0:
            init_index = i

        final_index = i
    return line[init_index+2:final_index]

def create_input(smiles_list):
    '''
        Create a directory with the same name for each element in smiles_list and generate input files based on the template
    '''
    work_dir = os.path.abspath('.')
    if not os.path.exists(join(work_dir,"config.yaml")):
        raise ValueError("There is no template .yaml file!")
        
    illegal_list = ["\\", "/", "|"]
    for i in range(len(smiles_list)):
        new_dir = join(join(work_dir, 'all_species'),str(i)+'_'+smiles_list[i])        
        for j in illegal_list:
            if j in smiles_list[i] :
                smiles_list[i].replace(j, '?')
        os.makedirs(new_dir)
        shutil.copy(join(work_dir,"config.yaml"), join(new_dir,"config.yaml"))    #Copy between dirs
        change_config(join(new_dir,"config.yaml"), smiles_list[i])

def change_config(filename, smiles_name):
    '''
        Generate config files based on species information
    '''
    with open(filename, 'r+') as f:
        lines = f.readlines()
        ads_line = 0
        for i in range(len(lines)):
            if 'ads' in lines[i][:10]:
                ads_line = i
        wds = lines[ads_line+1].split()
        end_wds = wds[-1]
        for wd in wds:
            if wd[0] == '"' and wd[-2] == '"':
                template_smiles = wd[1:-2]
                
    f1 = open(filename,"r")
    content = f1.read()
    f1.close()

    t1 = content.replace(template_smiles, smiles_name)
    f2 =  open(filename,"w") 
    f2.write(t1)
    f2.close

class Stem:
    '''
    Attributes:
        species: all species for subgraphs
        steps: elementary steps for subgraphs
    '''
    def __init__(self, species, steps):
        self.species = species
        self.steps = steps
    
def crnlog_readstems(stem_list):
    '''
        Return: list
            A list whose elements are Stem objects
    '''
    stem_value = []
    for i in range(len(stem_list)):
        _species, _steps = crnlog_analstem(stem_list[i])
        stem_value.append(Stem(_species, _steps))
    
    return stem_value
        
                
def crnlog_getstems(inputfile):
    '''
        Return: list
            A list consisting of the row of each subgraph in CRNGenerator_log.txt
    '''
    with open(inputfile, 'r') as f:
        f_lines = f.readlines()
        stem_list = []
        for i in range(len(f_lines)):
            if f_lines[i][:3] == '***':
                _startline = i
                while (not not f_lines[i].strip()):
                    i += 1
                _finalline = i
                stem_list.append(f_lines[_startline: _finalline])
        return stem_list
    
def crnlog_analstem(stem):
    '''
        Return species and elementary steps information of a subgraph
    '''
    stem_species = []
    for i in range(len(stem)):
        if stem[i][:2]== '--' and i != 1:
            _finalspeices = i
            
    for i in range(2, _finalspeices):
        stem_species.append(stem[i].split()[-1])
        
    stem_step = []
    for i in range(_finalspeices+1, len(stem)):
        _step = stem[i].split()[-5:]
        stem_step.append(_step)
    return stem_species, stem_step

def create_stem_dir(stem_value):
    '''
     Create directories for all subgraphs
    '''
    work_dir = os.path.abspath('.')
    if not os.path.exists(join(work_dir,"config.yaml")):
        raise ValueError("There is no template .yaml file!")
        
    for i in range(len(stem_value)):
        for j in range(len(stem_value[i].steps)):
            new_dir = join(join(join(work_dir, 'stems'), str(i)+'_stem'), str(j)+'_step')
            os.makedirs(new_dir)
            shutil.copy(join(work_dir, 'config.yaml'), join(new_dir, 'config_old.yaml'))    #用于目录之间的拷贝操作
            change_stemconfig(join(new_dir,'config_old.yaml'), join(new_dir,'config.yaml'), stem_value[i].steps[j])
            os.remove(join(new_dir, 'config_old.yaml'))
    

        
def change_stemconfig(filename, newfilename, stem_step):
    '''
        modify config file based on the information of one of the elementay steps of a Stem.
    '''
    _idxarrow = -1
    if '<-->' in stem_step:
        _idxarrow = stem_step.index('<-->')

    
    coads_list = []
    if _idxarrow == 1:
        ads_name = stem_step[0]
        coads_list.append(stem_step[2])
        coads_list.append(stem_step[4])
    elif _idxarrow == 3:
        ads_name = stem_step[4]
        coads_list.append(stem_step[0])
        coads_list.append(stem_step[2])
    else:
        raise ValueError('Wrong stem format!')
    
    with open(filename, 'r') as f_old:
        lines = f_old.readlines()
        ads_line = 0
        ads_wdgap = 0
        for i in range(len(lines)):
            if 'ads' in lines[i][:10]:
                ads_line = i
                
        for i in range(len(lines[ads_line+1])):
            if lines[ads_line+1][i] != ' ':
                ads_wdgap = i
                break
        wds = lines[ads_line+1].split()
        end_wds = wds[-1]
        for wd in wds:
            if wd[0] == '"' and wd[-2] == '"':
                template_smiles = wd[1:-2]
    
    with open(newfilename, 'w') as f_new:
        idx_line = -1
        for line in lines:
            idx_line += 1
            if idx_line == ads_line+1:
                modified_line = lines[ads_line+1].replace(template_smiles, ads_name)
                f_new.write(modified_line)
                f_new.write(lines[ads_line].replace('ads', 'coads'))
                f_new.write(' '*ads_wdgap+'- [s: "'+coads_list[0]+'", s: "'+coads_list[1]+'", 1, 1]')
                break
            else:
                f_new.write(line)