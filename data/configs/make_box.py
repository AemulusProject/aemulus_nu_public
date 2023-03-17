#!/usr/bin/env python
import os
import sys
import yaml
import numpy as np
from yaml import Loader
from run_zwindstroom import get_growth_expansion_info

def read_yaml(fname):
    with open(fname,'r') as fp:
        config = yaml.load(fp, Loader=Loader)
    return config

def write_template(infile,outfile,tvars):
    # read
    with open(os.path.expandvars(infile),'r') as fp:
        temp = fp.readlines()
    
    # fill template
    temp = "".join(temp)    
    temp = temp.format(**tvars)    
    
    # write to file
    with open(os.path.expandvars(outfile),'w') as fp:
        fp.write(temp)

def read_sigma8(fname):
    for line in open(os.path.expandvars(fname),'r'):
        pass
    items = line.strip().split()
    return float(items[-1])

def read_hzi(fname):
    hz = np.genfromtxt(fname)
    return hz[0,-1]

# code to make a box

sysname = sys.argv[1]
bnum = int(sys.argv[2])
npart = sys.argv[3]

if npart == '1400':
    CosmoFile = './LH_np7_n100_mnu05_randomized.txt'
elif npart == '2800':
    CosmoFile = './LH_np7_n100_mnu05_randomized.txt'

cosmos = np.genfromtxt(CosmoFile,names=True)
SeedFile = './seeds.dat'
seeds = np.genfromtxt(SeedFile)

SystemFile = os.path.join('systems',sysname,'%s.yaml' % sysname)
sysconfig = read_yaml(SystemFile)
sysconfig['Cores'] = sysconfig['Nodes_{}'.format(npart)]*sysconfig['TasksPerNode_{}'.format(npart)]
sysconfig['CoresPerTask'] = int(sysconfig['CoresPerNodeTotal']) // int(sysconfig['TasksPerNode_{}'.format(npart)])

BoxJobPath = os.path.expandvars(os.path.join(sysconfig['JobBase'], 'Box{}_{}'.format(bnum, npart)))
BoxOutputPath = os.path.expandvars(os.path.join(sysconfig['OutputBase'], 'Box{}_{}'.format(bnum, npart)))
RStarOutputPath = os.path.expandvars(os.path.join(sysconfig['OutputBase'], 'Box{}_{}'.format(bnum, npart), 'output', 'rockstar'))

# fill template vars
print(bnum)
print(cosmos['H0'][bnum])

# cosmology
h = cosmos['H0'][bnum]/100
tvars = {}
tvars['nu_mass_ev'] = cosmos['nu_mass_ev'][bnum]
tvars['nu_mass_ev_ov_3'] = tvars['nu_mass_ev'] / 3
tvars['Omeganuh2'] = cosmos['nu_mass_ev'][bnum] / (93.14)
tvars['OmegaNu'] = tvars['Omeganuh2'] / h**2
tvars['OmegaM'] = (cosmos['omch2'][bnum] + cosmos['ombh2'][bnum])/ h**2 + tvars['OmegaNu']
tvars['OmegaB'] = cosmos['ombh2'][bnum]/ h**2
tvars['OmegaCDM'] = cosmos['omch2'][bnum] / h**2
#print(tvars['OmegaM'], tvars['OmegaCDM'], tvars['OmegaB'], tvars['OmegaCDM'] + tvars['OmegaB'])
tvars['OmegaMh2'] = tvars['OmegaM']*h**2
tvars['OmegaBh2'] = tvars['OmegaB']*h**2
tvars['OmegaCDMh2'] = tvars['OmegaCDM']*h**2
tvars['OmegaL'] = 1.0 - tvars['OmegaM']
tvars['As'] = cosmos['As'][bnum]/1e9
print(tvars['As'])
tvars['SpectralIndex'] = cosmos['ns'][bnum]
tvars['w0'] = cosmos['w0'][bnum]
tvars['wa'] = 0.0
tvars['h'] = h
tvars['H0'] = tvars['h']*100.0
tvars['hierarchy'] = 'normal'
tvars['Neff'] = 3.046

# box
#tvars['zini'] = 99.0
tvars['zini'] = 12.0
tvars['aini'] = 1.0/(1.0 + tvars['zini'])
if npart=='1400':
    tvars['BoxL'] = 1050.0
    tvars['NpOneThird'] = 1400
    tvars['NmeshRSD'] = 700
    tvars['NpOneThird_nu'] = 1400
    tvars['NmOneThird'] = 2100
    tvars['NmOneThird_nu'] = 384
    tvars['Soft'] = 0.020
    tvars['MaxMem'] = 6000

elif npart=='2800':
    tvars['BoxL'] = 800.0
    tvars['NpOneThird'] = 2800
    tvars['NpOneThird_nu'] = 1400
    tvars['NmOneThird'] = 4200
    tvars['NmeshRSD'] = 700    
    tvars['NmOneThurd_nu'] = 512
    tvars['Soft'] = 0.0075
    tvars['MaxMem'] = 2850

# snapshots
tvars['NumSnaps'] = 30
tvars['SnapStartRedshift'] = 3.0
tvars['TimeOfFirstSnapshot'] = 1.0/(1.0 + tvars['SnapStartRedshift'])
# exp(timediff in log(a) / (NumSnaps - 1))
# use NumSnaps -1 since code always writes last snapshot
tvars['TimeBetSnapshot'] = np.exp((np.log(1.0)-np.log(tvars['TimeOfFirstSnapshot']))/(tvars['NumSnaps']-1))

# ics
tvars['seed'] = int(seeds[bnum])
tvars['GlassFile'] = 'grid_file_NU_1_CDM_1.dat'
tvars['ICsBaseName'] = 'ics'
tvars['ZAOutputDir'] = os.path.join(BoxOutputPath,'monofonic_za/')
tvars['3LPTOutputDir'] = os.path.join(BoxOutputPath,'monofonic/')
tvars['ICsNio'] = sysconfig['Nio_{}'.format(npart)]
tvars['CAMBHz'] = 'zwind_hubble.txt'

# lgadget
if npart == '1400':
    tvars['ZAFile'] = os.path.join(tvars['ZAOutputDir'],tvars['ICsBaseName'])
    tvars['3LPTFile'] = os.path.join(tvars['3LPTOutputDir'],tvars['ICsBaseName'])
    tvars['ICsFile'] = os.path.join(tvars['3LPTOutputDir'],tvars['ICsBaseName'])    
    
elif npart == '2800':
    tvars['ZAFile'] = os.path.join(tvars['ZAOutputDir'],tvars['ICsBaseName'])
    tvars['3LPTFile'] = os.path.join(tvars['3LPTOutputDir'],tvars['ICsBaseName'])
    tvars['ICsFile'] = os.path.join(tvars['3LPTOutputDir'],tvars['IcsBaseName'])
    
tvars['OutputDir'] = os.path.join(BoxOutputPath,'output')
tvars['SnapBaseName'] = 'snapshot'
tvars['Nio'] = sysconfig['Nio_{}'.format(npart)]
tvars['TimeLimit'] = sysconfig['TimeLimitHours']*60.0*60.0
tvars['Cores'] = sysconfig['Cores']
tvars['AnzuOutputDir'] = os.path.join(BoxOutputPath, 'anzu_fields_2xsmooth/') 
tvars['ValidationOutputDir'] = os.path.join(BoxOutputPath, 'validation', 'pk')

# make paths
for pth in [BoxJobPath,BoxOutputPath,RStarOutputPath,tvars['OutputDir'],
            tvars['ZAOutputDir'],tvars['3LPTOutputDir'],tvars['AnzuOutputDir'],
            tvars['ValidationOutputDir']]:
    os.system('mkdir -p %s' % pth)    

for pth in ['run', 'exec', 'scripts']:
    os.system('mkdir -p %s' % os.path.join(BoxJobPath,pth))

# copy exec and scripts over
os.system('cp -r %s/* %s/exec' % (sysconfig['ExecDir'],os.path.join(BoxJobPath)))
os.system('cp -r * %s' % (os.path.join(BoxJobPath,'scripts/.')))

# do zwind
basepath = '{}/run/'.format(BoxJobPath)
get_growth_expansion_info(tvars, basepath, sysname)
os.system('cp grid_file_NU_1_CDM_1.dat %s' % os.path.join(BoxJobPath,'run','grid_file_NU_1_CDM_1.dat'))

tvars['Hzi'] = read_hzi(os.path.join(BoxJobPath,'run','zwind_hubble.txt'))
print('Hzi: {}'.format(tvars['Hzi']))

write_template(os.path.join('./config_templates','monofonic.param'),
               os.path.join(BoxJobPath,'run','monofonic.param'),tvars)

write_template(os.path.join('./config_templates','monofonic_za.param'),
               os.path.join(BoxJobPath,'run','monofonic_za.param'),tvars)

write_template(os.path.join('./config_templates','gadget3.param'),
               os.path.join(BoxJobPath,'run','gadget3.param'),tvars)
os.system('cp %s %s' % (os.path.join('./config_templates','savelist.txt'), os.path.join(BoxJobPath,'run')))

write_template(os.path.join('./config_templates','rstar.param'),
               os.path.join(BoxJobPath,'run','aemulus_nu_halos.cfg'),tvars)
os.system('cp %s %s' % (os.path.join('./config_templates','bgc2.sh'), os.path.join(BoxJobPath,'run')))

write_template(os.path.join('./config_templates','anzu_fields.param'),
               os.path.join(BoxJobPath,'run','anzu_fields.param'),tvars)

# write job scripts
jobconfig = {}
jobconfig['JobName'] = 'Box{}_{}-i'.format(bnum, npart)
jobconfig['Nodes'] = sysconfig['Nodes_{}'.format(npart)]
jobconfig['TasksPerNode'] = sysconfig['TasksPerNode_{}'.format(npart)]
jobconfig['Cores'] = sysconfig['Cores']
jobconfig['CoresPerTask'] = sysconfig['CoresPerTask']
jobconfig['TimeLimitHours'] = sysconfig['TimeLimitHours']
jobconfig['Cores_anzu'] = int(npart)
jobconfig['Nodes_anzu'] = int(sysconfig['Nodes_{}'.format(npart)]) // 2
jobconfig['CoresPerTask_anzu'] = sysconfig['CoresPerNodeTotal'] * jobconfig['Nodes_anzu'] // jobconfig['Cores_anzu']
jobconfig['TimeLimitHours_anzu'] = '4'
jobconfig['ZACommand'] = '../exec/monofonic/build/monofonIC monofonic_za.param'
jobconfig['3LPTCommand'] = '../exec/monofonic/build/monofonIC monofonic.param'
jobconfig['MergeCommand'] = '../exec/scripts/merge_ics.py ngenic.param monofonic.param'
jobconfig['ICsFile'] = tvars['ICsFile']
jobconfig['BoxJobPath'] = BoxJobPath

write_template(os.path.join('./config_templates','bgc2.sh'),
               os.path.join(BoxJobPath,'run','bgc2.sh'),jobconfig)
os.system('chmod a+x %s' % (os.path.join(BoxJobPath,'run/bgc2.sh')))

jobconfig['LGadgetCommand'] = '../exec/g3_yb_hz_{}/P-Gadget3 gadget3.param'.format(npart)
jobconfig['Email'] = sysconfig['EMAIL']

if sysname == 'stampede2':
    pass
else:
    jobconfig['QOS'] = sysconfig['QOS']
    jobconfig['Repo'] = sysconfig['Repo']

jobconfig['Npart'] = npart

write_template(os.path.join('systems',sysname,'job.first.sh'),
               os.path.join(BoxJobPath,'run','job.first.sh'),jobconfig)

jobconfig['JobName'] = 'Box{}_{}-r'.format(bnum, npart)
write_template(os.path.join('systems',sysname,'job.restart.sh'),
               os.path.join(BoxJobPath,'run','job.restart.sh'),jobconfig)

jobconfig['JobName'] = 'Box{}_{}-rstar'.format(bnum, npart)
jobconfig['OutputDir'] = BoxOutputPath
jobconfig['RstarOutDir'] = RStarOutputPath
jobconfig['BoxNum']    = '%03d' % bnum
jobconfig['snapnum']   = '{snapnum}'
jobconfig['done']      = '{done}'

write_template(os.path.join('systems',sysname,'job.halos.sh'),
               os.path.join(BoxJobPath,'run','job.halos.sh'),jobconfig)

jobconfig['JobName'] = 'Box{}_{}-a'.format(bnum, npart)
jobconfig['Prefix'] = ''
if 'sherlock' not in sysname:
    write_template(os.path.join('systems',sysname,'job.archive.sh'),
                   os.path.join(BoxJobPath,'run','job.archive.sh'),jobconfig)

if 'perlmutter' in sysname:
    write_template(os.path.join('systems',sysname,'job.unarchive.sh'),
                   os.path.join(BoxJobPath,'run','job.unarchive.sh'),jobconfig)
    

jobconfig['JobName'] = 'Box{}_{}-v'.format(bnum, npart)
jobconfig['OutputDir'] = BoxOutputPath
jobconfig['ICsFile']    = tvars['ICsFile']
jobconfig['BoxJobPath'] = BoxJobPath
jobconfig['SnapString'] = '42 70 99'
jobconfig['RedshiftString'] = '{} 2.977 1.0186 0.0'.format(tvars['zini'])
jobconfig['BoxL'] = tvars['BoxL']

write_template(os.path.join('systems',sysname,'job.validate.sh'),
               os.path.join(BoxJobPath,'run','job.validate.sh'),jobconfig)

jobconfig['JobName'] = 'Box{}_{}-anzu'.format(bnum, npart)
write_template(os.path.join('systems',sysname,'job.anzu.sh'),
               os.path.join(BoxJobPath,'run','job.anzu.sh'),jobconfig)
