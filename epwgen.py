#______________________RUN_ENVIRONMENT_______________________#
#prefix for any data files
pf = 'KTO'  
#name of the parent directory created by this script
base_dir = pf 
#name of job on cluster. Be mindful of regex metacharacters
jobname = pf
#location of the pseudopotentials where the calculations are run. Use absolute paths to directories througout.
pps_dir = 'ppsdir' 
#number of processors to be used for scf, bands and nscf calculations
num_of_cpu_scf = 16
#number of processors to be used for the phonon calculations (i.e. if you split the calculation into calculations of
#irreducible modes, each of these calculations will be run with the specified amount of CPUs)
num_of_cpu_ph = 8
#number of processors to be used for any epw calculation
num_of_cpu_epw = 48
#memory per cpu
ram = 1024
#time limit in hours for scf, nscf and bands calculations
scf_t = 4
#time limit in hours for phonon calculations, i.e. for the whole calculation or per irreducbile q-point or irreducible representation
#depending on whether the phonon calculations are split or not
q_t = 4
#time limit in hours for epw calculations 
epw_t = 4
#specify if phonon calculations should be split into calculations of irreducible q-points
split_q = True
#specify if phonon calculations should additionally be split into calculations of irreducible represenations.
#split_q is automatically enabled if split_irr is enabled and split_irr is disabled if split_q is disabled.
split_irr = True
#if the right electron and phonon bands have already been calculated in another parent directory you can set
#ref_bands to True. epw.sh will then run with the corresponding input of the parent directory s.t. the bands
#don't have to be recalculated (start from 2nd submission script).
#If you set ref_bands to False you need to calculate the bands first (start form 1st submission script).
ref_bands = False
#Parent directory containing the already calculated bands (in the respective ELB and PHB directories).
#You only need to specify this if ref_bands = True.
#If ref_bands = True the prefix has to be the same as the one for the phonon calculations
bands_dir = ''
#list of names of necessary modules. Should include: 
#    -Quantum Espresso including EPW package
#    -python
#    -gnuplot (optional)
modules = ['quantum_espresso/6.2.1', 'python/3.6.0', 'gnuplot']

#note: time limits for calculations which should not take long have been set to 4h (shortest cluster limit)
#these include wannier, bands_ip, q2r, matdyn


#__________________________STRUCTURE__________________________#
#lattice vectors in Angstrom
#Syntax:
#'''
#CELL_PARAMETERS angstrom
#$v1(1) $v1(2) $v1(3)
#$v2(1) $v2(2) $v2(3)  
#$v3(1) $v3(2) $v3(3)
#'''
lattice = '''
CELL_PARAMETERS angstrom
   4.006337510   0.000000000   0.000000000
   0.000000000   4.006337510   0.000000000
   0.000000000   0.000000000   4.006337510
'''
#number of atoms in the unit cell
num_of_atoms = 5

#number of distinct atom types in unit cell
num_of_atom_types = 3                 

#atom type specification as required in pw.x input file
#Syntax:
#'''
#ATOMIC_SPECIES
#$atom_label_1 $atom_weight_1 $pps_filename_1
#$atom_label_2 $atom_weight_2 $pps_filename_2
#etc.
#'''
atoms = '''      
ATOMIC_SPECIES
K 39.0983 K_ONCV_PBE_FR_4.0.upf
Ta 180.94788 Ta_ONCV_PBE_FR_4.0.upf
O 15.9994 O_ONCV_PBE_FR_4.0.upf
'''

#Atomic positions in relative positions to the lattice vectors
#Syntax: 
#'''
#ATOMIC_POSITIONS crystal
#   $atom_label_1 $x1 $y1 $z1
#   $atom_label_2 $x2 $y2 $z2
#   $atom_label_3 $x3 $y3 $z3
#   etc.
#'''

atom_positions = '''
ATOMIC_POSITIONS crystal
K  0.000000 0.000000 0.000000 
Ta 0.500000 0.500000 0.500000
O  0.000000 0.500000 0.500000 
O  0.500000 0.000000 0.500000 
O  0.500000 0.500000 0.000000
'''

#____________________CALCULATION_PARAMETERS____________________#
#plane wave cut-off energy in Ry
e_cut = 80  
#coarse k-grid size
k_x = 8
k_y = 8
k_z = 8
#coarse q-grid size. The q and k-grid sizes need to be whole multiples of each other
q_x = 2
q_y = 2
q_z = 2
#fine grid sizes
kf_x = 24
kf_y = 24
kf_z = 24
qf_x = 24
qf_y = 24
qf_z = 24
#for random fine grids, set random_sampling to True and specify the number of respective random points. Note that for
#the Eliashberg functions no random fine grids can be used which is why you still need to specify fine grids above in that case
random_sampling = True
random_nkf = 100000
random_nqf = 50000

#number of bands calculated. Set to 0 if this should be determined automatically.
num_of_bands = 0    
#amount of additional charge per unit cell (electrons negative, holes positive).
nelect = -0.01
#accoustic sum rule type ('simple', 'crystal', 'one-dim', 'zero-dim', 'no' to disable)
asr_type = 'simple'

#scf convergence threshold in Ry
delta_scf = '1.0d-8'
#phonon convergence threshold in Ry^2
delta_ph = '1.0d-14'

#scf diagonalization algorithm ('cg' or 'david'). 
diag_algo = 'cg'

#spin-orbit coupling
soc = False
                    

#____________________HIGH_SYMMETRY_LINES_______________________#
#define an array of high symmetry points.
#Syntax: ['$point_label_1 $k_x1 $k_y1 $k_z1', '$point_label_2 $k_x2 $k_y2 $k_z2', etc]
hisym_points = ['G 0 0 0', 
                'X 0.5 0 0',
                'M 0.5 0.5 0',
                'R 0.5 0.5 0.5'
               ]
#define the high symmetry paths using the point labels defined in hisym_points
#Syntax: ['$point_label_1', '$point_label_2', '$point_label_1', etc. ]
hisym_path = ['G', 'M', 'X', 'G', 'R', 'M', 'G']

#number of interpolation points per high symmetry line
path_prec = 100


#_________________________WANNIERIZATION_______________________#
#The parameters for the wannierization can also be changed in the epw.sh script once you have obtained 
#the electron bands and you want to adjust them. 

#number of bands at and above the Fermi energy that get wannierized. 
#This needs to correspond to the right number implied by your specified projections.
num_of_wan = 3
#array of initial projections for Wannier functions. Check wannier90 documentation for syntax.
#Set to 'random' if you have no initial guess
#Syntax: ['$proj1', '$proj2', etc.] (example: ['random', 'W:l=2,mr=2,3,5'])
wannier_init = ['Ta:l=2,mr=2,3,5']

#automatic wannierization window determination
auto_window = True

#manual wannierization window specification. You only need to specify these parameters
#if auto_window = False or if the automatic determination doesn't work.
#highest band that does not get wannierized,i.e. the $num_of_wan bands after band $wan_min get wannierized
wan_min = 0
#inner wannierization window (in eV)
inner_bottom = 0.0
inner_top = 0.0
#outer wannierization window (n eV)
outer_bottom = 0.0
outer_top = 0.0

#___________________CRITICAL TEMPERATURE_______________________#
#minimum temperature in K for which superconducting gap is calculated
T_min = 0.1
#maximum temperature in K
T_max = 5
#temperature points between T_min and T_max for which superconducting gap is calculated
num_of_T = 10
#___________________________END________________________________#
#____________________NOW_RUN_THE_SCRIPT________________________#

import os
import string
import math

if not os.path.exists(base_dir):
    os.makedirs(base_dir) 
    
if not os.path.exists(os.path.join(base_dir, 'ELB')):
    os.makedirs(os.path.join(base_dir, 'ELB')) 
    
if not os.path.exists(os.path.join(base_dir, 'PHB')):
    os.makedirs(os.path.join(base_dir, 'PHB')) 
    
if not os.path.exists(os.path.join(base_dir, 'EPM')):
    os.makedirs(os.path.join(base_dir, 'EPM')) 

if not os.path.exists(os.path.join(base_dir, 'ISO')):
    os.makedirs(os.path.join(base_dir, 'ISO')) 
    
path_prec = str(path_prec)
num_of_hsp = len(hisym_path)  

#get the k-points for bands.in and matdyn.in
kpoints = ''''''
for hisym_point in hisym_path:
    for i in range(0, len(hisym_points)):
        hisym_point_comp = hisym_points[i].split(' ')[0]
        if hisym_point == hisym_point_comp:
            kpoints += hisym_points[i].split(' ')[1] + '\t' + hisym_points[i].split(' ')[2] + '\t' + \
                                    hisym_points[i].split(' ')[3] + '\t' + path_prec
            kpoints += '\n'            

#get the k-points for wannier90
index = 0
kpoints_wannier = ''''''
for hisym_point in hisym_path:
    for i in range(0, len(hisym_points)):
        hisym_point_comp = hisym_points[i].split(' ')[0]
        if hisym_point == hisym_point_comp:
            #first point in path
            if index == 0:
                kpoints_wannier += 'wdata(' + str(index+3) + ') = ' + '\'' + hisym_points[i].split(' ')[0] + '\t' + \
                                   hisym_points[i].split(' ')[1] + '\t' + \
                                   hisym_points[i].split(' ')[2] + '\t' + hisym_points[i].split(' ')[3] + '\t'
                index += 1
            #last point in path
            elif index == (len(hisym_path) - 1):
                kpoints_wannier += hisym_points[i].split(' ')[0] + '\t' + hisym_points[i].split(' ')[1] + '\t' + \
                                   hisym_points[i].split(' ')[2] + '\t' + hisym_points[i].split(' ')[3] + '\'' 
           #intermediate point
            else:
                kpoints_wannier += hisym_points[i].split(' ')[0] + '\t' + hisym_points[i].split(' ')[1] + '\t' + \
                                   hisym_points[i].split(' ')[2] + '\t' + hisym_points[i].split(' ')[3] + '\''
                kpoints_wannier += '\n'
                kpoints_wannier += '    ' + 'wdata(' + str(index+3) + ') = ' + '\'' + hisym_points[i].split(' ')[0] + '\t' + \
                                   hisym_points[i].split(' ')[1] + '\t' + \
                                   hisym_points[i].split(' ')[2] + '\t' + hisym_points[i].split(' ')[3] + '\t'
                index += 1
                
kpoints_wannier += '\n'
kpoints_wannier += '    ' + 'wdata(' + str(len(hisym_path) + 2) + ') = \'end kpoint_path\''
kpoints_wannier += '\n'
kpoints_wannier += '    ' + 'wdata(' + str(len(hisym_path) + 3) + ') = \'bands_num_points =' +  path_prec + '\''

#get the projections for wannier90
index = 0
projections = ''''''
for p in wannier_init:
    projections += 'proj(' + str(index+1) + ') = \'' + str(p) + '\''+ '\n' + '\t' 
    index += 1
    
#get the projections for epw.sh
projections_epw = '''('''
for p in wannier_init:
    projections_epw +=  '\"' + p + '\" '
projections_epw = projections_epw[:-1]
projections_epw += ')'
    
#function to create job submission files
#LSF
def make_job_sub(_jobname,_num_of_cpu, _ram, _time, _infile, _outfile, _program, _condition, _escape=False):
    job_sub = ''''''
    job_sub += '#!/bin/bash' + '\n'
    job_sub += '#BSUB -n ' + str(_num_of_cpu)  + '\n'
    job_sub += '#BSUB -R \"' + "rusage[mem=" + str(_ram)  + ']\"\n'
    job_sub += '#BSUB -W ' + str(_time)  + ':00' + '\n'
    job_sub += '#BSUB -o ' + _jobname  + '.o' + '\n'
    job_sub += '#BSUB -e ' + _jobname  + '.e' + '\n'
    job_sub += '#BSUB -J ' + _jobname + '\n'
    if not _condition == '':
        job_sub += '#BSUB -w ' + _condition  + '\n'
    if not _program == '':
        job_sub += 'mpirun ' + _program + ' -npool ' + str(_num_of_cpu) + ' -in ' + _infile + ' > ' + _outfile
    if _escape == True:
        job_sub = job_sub.replace("#", "\#")
    return str(job_sub)

#function to check if conditions need to be removed and if submission scripts need to be submitted
#LSF
def check_cond_sub(_index, _scf = False):
    cond_sub = ''''''
    if not _scf:
        cond_sub += 'if (($calc_start == ' + str(_index) + '))' + '\n' + 'then' + '\n'
        cond_sub += 'line=$(grep -n \"#BSUB -w\" job.sh' + ' | cut -d : -f 1)' + '\n'
        cond_sub += 'sed -i ${line}d job.sh' + '\n' + 'fi' + '\n'
    
    cond_sub += 'if ! $no_sub' + '\n' + 'then' + '\n'
    cond_sub += 'bsub<job.sh' + '\n' + 'fi' + '\n'
    return str(cond_sub)
    
    
#function to generate coarse and fine grid specification for epw input
def generate_fine_grids(_random_sampling, _random_nkf, _random_nqf,
                   _kf_x, _kf_y, _kf_z, _qf_x, _qf_y, _qf_z, _file_q = False, _prefix = ""):
    grids = ''''''
    if _random_sampling:
        grids += 'rand_k = .true.' + '\n' + '    '       
        grids += 'rand_nk = ' + str(_random_nkf) + '\n' + '    '
        if not _file_q:
            grids += 'rand_q = .true.' + '\n' + '    '
            grids += 'rand_nq = ' + str(_random_nqf) + '\n' + '    '
        else:
            grids +=  'filqf = \'' + str(_prefix) + '_band.kpt\'' + '\n' + '    '
    else:
        grids += 'nkf1 = ' + str(_kf_x) + '\n' + '    '
        grids += 'nkf2 = ' + str(_kf_y) + '\n' + '    '
        grids += 'nkf3 = ' + str(_kf_z) + '\n' + '    '
        if not _file_q:
            grids += 'nqf1 = ' + str(_qf_x) + '\n' + '    '
            grids += 'nqf2 = ' + str(_qf_y) + '\n' + '    '
            grids += 'nqf3 = ' + str(_qf_z) + '\n' + '    '
        else:
            grids +=  'filqf = \'' + str(_prefix) + '_band.kpt\'' + '\n' + '    '
    return str(grids)

#function to generate module loading commands 
def generate_modules(_modules):
    module_commands = '''module purge \n'''
    for module in _modules:
        module_commands += 'module load ' + module + ' &>/dev/null \n'
    return module_commands
    
asr_enable = '.false.'
if asr_type == 'simple' or asr_type == 'crystal' or asr_type == 'one-dim' or asr_type == 'zero-dim':
    asr_enable = '.true.'

if auto_window == True:
    auto_window = 'true'
else:
    auto_window = 'false'
    
if split_q == True:
    split_q = 'true'
else:
    split_q = 'false'
    
if split_irr == True:
    split_irr = 'true'
else:
    split_irr = 'false'

dvscf_dir = '../PHB/save'
if ref_bands == True:
    ref_bands = 'true'
    dvscf_dir = bands_dir +  '/PHB/save'
else:
    ref_bands = 'false'

noncolin_enable = '.false.'
lspinorb_enable = '.false.'
if soc == True:
    noncolin_enable = '.true.'
    lspinorb_enable = '.true.'
   
        
#_______________________________________INPUT_LISTS_______________________________________#
#_________________________________________________________________________________________#
scf_in = ['''
&CONTROL
    calculation   = 'scf'
    restart_mode  = 'from_scratch'      
    prefix        = '{pf}'              
    pseudo_dir    = '{pps_dir}'         
    outdir        = './'
/
&SYSTEM   
    ibrav       = 0      
    nat         = {num_of_atoms}       
    ntyp        = {num_of_atom_types}   
    ecutwfc     = {e_cut}               
    occupations = 'smearing'           
    degauss     = 7.35d-4               
    tot_charge  = {nelect}   
    noncolin    = {noncolin_enable}
    lspinorb    = {lspinorb_enable}
/
&ELECTRONS
    conv_thr    = {delta_scf}              
    diagonalization = '{diag_algo}'              
/
{lattice}
{atoms}
{atom_positions}
K_POINTS automatic
{k_x} {k_y} {k_z} 0 0 0
'''.format(pf = pf, pps_dir = pps_dir, lattice = lattice,
           num_of_atoms = num_of_atoms, num_of_atom_types = num_of_atom_types, e_cut = e_cut, nelect = nelect,
           noncolin_enable = noncolin_enable, lspinorb_enable = lspinorb_enable, delta_scf = delta_scf,
           diag_algo = diag_algo, atoms = atoms, atom_positions = atom_positions, k_x = k_x, k_y = k_y, k_z = k_z)]

bands_in = ['''
&CONTROL
    calculation   = 'bands'
    restart_mode  = 'from_scratch'      
    prefix        = '{pf}'              
    pseudo_dir    = '{pps_dir}'         
    outdir        = './'
/
&SYSTEM    
    ibrav       = 0            
    nat         = {num_of_atoms}        
    ntyp        = {num_of_atom_types}   
    ecutwfc     = {e_cut}               
    occupations = 'smearing'            
    degauss     = 7.35d-4               
    tot_charge  = {nelect}              
    nbnd        = {num_of_bands}  
    noncolin    = {noncolin_enable}
    lspinorb    = {lspinorb_enable}
/
&ELECTRONS
    conv_thr    = {delta_scf}               
    diagonalization = '{diag_algo}'
/
{lattice}
{atoms}
{atom_positions}
K_POINTS crystal_b
{num_of_hsp}
{kpoints}
'''.format(pf = pf, pps_dir = pps_dir, lattice = lattice,
           num_of_atoms = num_of_atoms, num_of_atom_types = num_of_atom_types, e_cut = e_cut, nelect = nelect,
           noncolin_enable = noncolin_enable, lspinorb_enable = lspinorb_enable, num_of_bands = num_of_bands, 
           delta_scf = delta_scf, diag_algo = diag_algo, atoms = atoms, atom_positions = atom_positions, 
           num_of_hsp = num_of_hsp, kpoints = kpoints)]

bands_ip_in = ['''
&BANDS
 prefix='{pf}'
 outdir='.'
 filband = '{pf}bands.dat'
 lsym=.true.
/
'''.format(pf = pf)]

#the missing "/" at the end in ph_in is intended. Don't put one there.
ph_in = ['''
&INPUTPH
    prefix   = '{pf}'       
    outdir   = './'          
    fildyn   = '{pf}.dyn'    
    fildvscf = 'dvscf'       
    ldisp    = .true.        
    nq1      = {q_x}    
    nq2      = {q_y} 
    nq3      = {q_z} 
    asr      = {asr_enable}        
    tr2_ph   = {delta_ph}
    recover = .true.
'''.format(pf = pf, q_x = q_x, q_y = q_y, q_z = q_z, asr_enable = asr_enable, delta_ph = delta_ph)]

q2r_in = ['''
&INPUT
   fildyn='{pf}.dyn'
   zasr='{asr_type}'
   flfrc='{pf}.fc'           
/
'''.format(pf = pf, asr_type = asr_type)]

matdyn_in = ['''
&INPUT
    asr='{asr_type}'
    flfrc='{pf}.fc'
    flfrq='{pf}.freq'       
    q_in_band_form=.true.
/
{num_of_hsp}
{kpoints}
'''.format(pf = pf, num_of_hsp = num_of_hsp, kpoints = kpoints, asr_type = asr_type)]

#kpoints in nscf_in get appended when epw.sh is executed
nscf_in = ['''
&CONTROL
    calculation   = 'nscf',
    restart_mode  = 'from_scratch'       
    prefix        = '{pf}'              
    pseudo_dir    = '{pps_dir}'         
    outdir        = './'
/
&SYSTEM    
    ibrav       = 0          
    nat         = {num_of_atoms}        
    ntyp        = {num_of_atom_types}  
    ecutwfc     = {e_cut}               
    occupations = 'smearing'            
    degauss     = 7.35d-4               
    tot_charge  = {nelect}              
    nbnd        = {num_of_bands}
    nosym       = .true.     
    noncolin    = {noncolin_enable}
    lspinorb    = {lspinorb_enable}
/
&ELECTRONS
    conv_thr    = {delta_scf}               
    diagonalization = '{diag_algo}'
/
{lattice}
{atoms}
{atom_positions}
'''.format(pf = pf, pps_dir = pps_dir, lattice = lattice,
           num_of_atoms = num_of_atoms, num_of_atom_types = num_of_atom_types, e_cut = e_cut, nelect = nelect,
           noncolin_enable = noncolin_enable, lspinorb_enable = lspinorb_enable, 
           num_of_bands = num_of_bands, delta_scf = delta_scf, diag_algo = diag_algo, atoms = atoms, atom_positions = atom_positions)]

wannier_in = ['''
&inputepw
    prefix      = '{pf}'
    outdir      = './'
    dvscf_dir   = '{dvscf_dir}'
    
    elph        = .false.
    epbwrite    = .true.         
    epwwrite    = .true.         
    
    wannierize  = .true.         
    nbndsub     =  {num_of_wan}             
    nbndskip    =  {wan_min}            
    num_iter    = 10000           
    dis_win_min = 0.0            
    dis_win_max = 0.0
    dis_froz_min = 0.0            
    dis_froz_max = 0.0
    
    {projections}      
    wdata(1) = 'bands_plot = .true.'
    wdata(2) = 'begin kpoint_path'
    {kpoints_wannier}
        
    nk1 = {k_x}
    nk2 = {k_y}
    nk3 = {k_z}
    
    nq1 = {q_x}
    nq2 = {q_y}
    nq3 = {q_z}
    
    {fine_grids}
/
'''.format(pf = pf, dvscf_dir = dvscf_dir, num_of_wan = num_of_wan, wan_min = wan_min, projections = projections,
           kpoints_wannier = kpoints_wannier, asr_type = asr_type, k_x = k_x, k_y = k_y, k_z = k_z,
           q_x = q_x, q_y = q_y, q_z = q_z,
           fine_grids = generate_fine_grids(random_sampling,random_nkf,random_nqf,kf_x,kf_y,kf_z,qf_x,qf_y,qf_z))]

wannier_epm_in = ['''
&inputepw
    prefix      = '{pf}'
    outdir      = './'
    dvscf_dir   = '{dvscf_dir}'
    
    elph        = .true.   
    epbwrite    = .true.
    epwwrite    = .true.
    
    wannierize  = .true.         
    nbndsub     =  {num_of_wan}             
    nbndskip    =  {wan_min}            
    num_iter    = 10000           
    dis_win_min = 0.0            
    dis_win_max = 0.0
    dis_froz_min = 0.0            
    dis_froz_max = 0.0
    
    {projections}      
    wdata(1) = 'bands_plot = .true.'
    wdata(2) = 'begin kpoint_path'
    {kpoints_wannier}
   
    !fsthick     = 10                         
    !degaussw    = 0.1   
    
    asr_typ     = '{asr_type}'
    
    nk1 = {k_x}
    nk2 = {k_y}
    nk3 = {k_z}
    
    nq1 = {q_x}
    nq2 = {q_y}
    nq3 = {q_z}
    
    {fine_grids}
/
'''.format(pf = pf, dvscf_dir = dvscf_dir, num_of_wan = num_of_wan, wan_min = wan_min, projections = projections,
           kpoints_wannier = kpoints_wannier, asr_type = asr_type, k_x = k_x, k_y = k_y, k_z = k_z,
           q_x = q_x, q_y = q_y, q_z = q_z,
           fine_grids = generate_fine_grids(random_sampling,random_nkf,random_nqf,kf_x,kf_y,kf_z,qf_x,qf_y,qf_z))]

ph_lw_in = ['''
&inputepw
    prefix      = '{pf}'
    outdir      = './'
    dvscf_dir   = '{dvscf_dir}'

    epwread     = .true.           
    epwwrite    = .false.           
    kmaps       = .true.

    elph        = .true.          
    
    phonselfen = .true.          
    delta_approx = .true.          
   
    !fsthick     = 10                      
    !degaussw    = 0.1          
    
    efermi_read = .true.
    fermi_energy = 0.0
    
    nk1 = {k_x}
    nk2 = {k_y}
    nk3 = {k_z}
    
    nq1 = {q_x}
    nq2 = {q_y}
    nq3 = {q_z}
    
    {fine_grids}
/
'''.format(pf = pf, dvscf_dir = dvscf_dir,
           kf_x = kf_x, kf_y = kf_y, kf_z = kf_z, k_x = k_x, k_y = k_y, k_z = k_z, q_x = q_x, q_y = q_y, q_z = q_z,
           fine_grids = generate_fine_grids(random_sampling,random_nkf,random_nqf,kf_x,kf_y,kf_z,qf_x,qf_y,qf_z,True,pf))]


a2F_in = ['''
&inputepw
    prefix      = '{pf}'
    outdir      = './'
    dvscf_dir   = '{dvscf_dir}'

    epwread     = .true.           
    epwwrite    = .false.           
    kmaps       = .true.

    elph        = .true.           
    
    !fsthick     = 10                     
    !degaussw    = 0.1              
    
    phonselfen  = .true.
    delta_approx = .true.
    a2f         = .true.           
    
    efermi_read = .true.
    fermi_energy = 0.0

    nk1 = {k_x}
    nk2 = {k_y}
    nk3 = {k_z}
    
    nq1 = {q_x}
    nq2 = {q_y}
    nq3 = {q_z}
    
    {fine_grids}
/
'''.format(pf = pf, dvscf_dir = dvscf_dir, kf_x = kf_x, kf_y = kf_y, kf_z = kf_z, qf_x = qf_x, qf_y = qf_y, qf_z = qf_z,
           k_x = k_x, k_y = k_y, k_z = k_z, q_x = q_x, q_y = q_y, q_z = q_z,
           fine_grids = generate_fine_grids(random_sampling,random_nkf,random_nqf,kf_x,kf_y,kf_z,qf_x,qf_y,qf_z))]

eliashberg_iso = ['''
&inputepw
    prefix      = '{pf}'
    outdir      = './'
    dvscf_dir   = '{dvscf_dir}'

    elph        = .true.          
    eliashberg  = .true.          

    epbwrite    = .true.        
    epwwrite    = .true.         
    ephwrite    = .true.          
    epbread     = .false.
    epwread     = .false.
    kmaps       = .false.
    
    wannierize  = .true.         
    nbndsub     =  {num_of_wan}             
    nbndskip    =  {wan_min}            
    num_iter    = 10000           
    dis_win_min = 0.0            
    dis_win_max = 0.0
    dis_froz_min = 0.0            
    dis_froz_max = 0.0
    
    {projections}      
    wdata(1) = 'bands_plot = .true.'
    wdata(2) = 'begin kpoint_path'
    {kpoints_wannier}

    laniso      = .false.         
    liso        = .true.          

    limag       = .true.
    lpade       = .true.
    lacon       = .true.
    nsiter      = 500             
    conv_thr_iaxis = 1.0d-3       
    conv_thr_racon = 1.0d-3
            
    tempsmin = {T_min}            
    tempsmax = {T_max}
    nstemp = {num_of_T}

    muc     = 0.09                 

    mp_mesh_k = .true.                
    
    efermi_read = .true.
    fermi_energy = 0.0
    
    nk1 = {k_x}
    nk2 = {k_y}
    nk3 = {k_z}
    
    nq1 = {q_x}
    nq2 = {q_y}
    nq3 = {q_z}

    {fine_grids}
/
'''.format(pf = pf, dvscf_dir = dvscf_dir, num_of_wan = num_of_wan, wan_min = wan_min, projections = projections, 
           kpoints_wannier = kpoints_wannier, T_min = T_min, T_max = T_max, num_of_T = num_of_T, kf_x = kf_x, kf_y = kf_y, kf_z = kf_z,
           qf_x = qf_x, qf_y = qf_y, qf_z = qf_z, k_x = k_x, k_y = k_y, k_z = k_z, q_x = q_x, q_y = q_y, q_z = q_z,
           fine_grids = generate_fine_grids(False,random_nkf,random_nqf,kf_x,kf_y,kf_z,qf_x,qf_y,qf_z))]

#_________________________________________________________________________________________#
#_________________________________________________________________________________________#

#__________________________________SUBMISSION_SCRIPTS_____________________________________#
#_________________________________________________________________________________________#

#submission script for all calculations regarding electron and phonon bands
ep_bands_sh = ['''
#!/bin/bash
#Execute from the parent directory (containing ELB, PHB, EPM, ISO).

#Specify the starting and ending point of calculation. E.g. if the phonon calculations stopped because the time limit...
#...was too short, set calc_start=6 and run the submisison script again. The calculations will continue from where they stopped.
#1: electron bands scf calculation
#2: electron bands calculation
#3: electron bands interpolation
#4: phonon bands scf calculation
#5: phonon calculation initialization
#6: phonon calculations for every irreducible q-point/representation
#(If the phonon calculations are split, #6 will set up a manager that automatically splits the phonon calculation...
#... You can kill the manager while it's making the submissions (e.g. if you notice that it generates too many jobs) and...
#...resubmit starting from #6 where it will continue from where it stopped. If you restart, restart from #6 not from #5.)
#7: phonon calculation finalization
#8: phonon q2r run
#9: phonon matdyn run
#(#10: tidy up - deletes everything but the input/output files, the bands (i.e. electron and phonon bands) and the data...
#      ...needed for EPM calculations. This is not necessary but saves space. Usable only as a ...
#      ...seperate (calc_start = calc_end) execution. Only use once you are sure that everything worked.)
calc_start=1
calc_end=9

#If you need a specific submission script (e.g. for the q2r run) set calc_start = calc_end (e.g. = 8) and set no_sub
#to true. You'll find the job script as job.sh in the corresponding directory.
no_sub=false

#load modules 
{module_commands}


#________________________________EXECUTE_THE_SCRIPT________________________________#
#__________________________________________________________________________________#
split_q={split_q}
split_irr={split_irr}

#__________PREPARE_INPUT___________#
base_dir=$(pwd)
cp $base_dir/ELB/scf.in $base_dir/PHB

#__________JOB_SUBMISSION___________#
#electron bands calculation
cd $base_dir/ELB

if (($calc_start < 10)) && (($calc_end == 10))
then
echo "Error: please only tidy up once you are sure that all calculations have finished correctly"
exit
fi

if (($calc_start > $calc_end))
then
echo "Error: calculation order makes no sense"
exit
fi

#scf
if (($calc_start <= 1)) && (($calc_end >= 1))
then

cat > job.sh << EOF
{elb_scf_sub}
#get the number of bands from the scf.out file and put them in the bands.in file
num_of_bands=\$(grep nbnd bands.in | awk '{{print \$3}}')
if [ \$num_of_bands -eq 0 ]
then
    num_of_bands=\$(grep -m1 "number of Kohn-Sham states" scf.out | awk '{{print \$5}}')
    line=\$(grep -n nbnd bands.in | cut -d : -f 1)
    sed -i "\${{line}}s/[^ ]*[^ ]/\${{num_of_bands}}/3" bands.in
fi
EOF

{elb_scf_cond_sub}
fi

#bands
if (($calc_start <= 2)) && (($calc_end >= 2))
then

cat > job.sh << EOF
{elb_sub}
EOF

{elb_cond_sub}
fi

#bands_ip
if (($calc_start <= 3)) && (($calc_end >= 3))
then

cat > job.sh << EOF
{elb_ip_sub}
EOF

{elb_ip_cond_sub}
fi

#phonon calculation
cd $base_dir/PHB

#scf
if (($calc_start <= 4)) && (($calc_end >= 4))
then

cat > job.sh << EOF
{ph_scf_sub}
EOF

{ph_scf_cond_sub}
fi

#phonon calculations
#initialize the ph calculations
if (($calc_start <= 5)) && (($calc_end >= 5))
then

cat > job.sh << EOF
{ph_init_sub}
cp ph.in ph_start.in
echo "    start_irr = 0" >> ph_start.in
echo "    last_irr  = 0" >> ph_start.in
echo "    /" >> ph_start.in
mpirun ph.x -npool 1 -in ph_start.in > ph_start.out
EOF

{ph_init_cond_sub}
fi

#set up a manager that splits the phonon calculation
#IF CALC_START
if (($calc_start <= 6)) && (($calc_end >= 6))
then

#if no parallelization is wanted start the phonon calculation normally
#IF NOSPLIT
if ! $split_q
then
cp ph.in ph_all.in
echo "/" >> ph_all.in

cat > job.sh << EOF
{ph_all_sub}
EOF
#ENDIF NOSPLIT

#IF SPLIT_Q
#if irreducible q-point parallelization is enabled
elif $split_q && ! $split_irr
then

cat > job.sh << EOF
{ph_manager_sub}
#get the number of irreducible q-points
irr_qs=\$(sed "2q;d" {pf}.dyn0 | awk '{{print $1}}')

#run a phonon calculation for every irreducible q-point
for ((q=1; q <= \$irr_qs; q++))
do

#make directories for the irreducible q-point and copy necessary stuff there
if [ ! -d  q\${{q}}/_ph0/{pf}.phsave ]
then
mkdir -p q\${{q}}/_ph0/{pf}.phsave
cp -r {pf}.* q\${{q}}
cp -r _ph0/{pf}.phsave/* q\${{q}}/_ph0/{pf}.phsave
fi

cd q\${{q}}
#prepare the input file
cp ../ph.in ph_q\${{q}}.in
echo "    start_q = \$q" >> ph_q\${{q}}.in
echo "    last_q = \$q" >> ph_q\${{q}}.in
echo "/" >>  ph_q\${{q}}.in

#make the job file
cat > job_temp.sh << EOF1
{ph_q_sub}
\#delete the large wave functions in the q-directory immediately after a calculation is done
if [ -f ph_q\${{q}}.out ] && [ \${{q}} -ne 1 ]
then
status=\\\\\$(tail -n2 ph_q\${{q}}.out | head -n1 | awk '{{print \\\\\$2}}')
if [ "\\\\\$status" == "DONE." ]
then 
rm _ph0/{pf}.q_\${{q}}/*.wfc*
rm {pf}.save/wfc**.dat
rm *.wfc*
fi
fi
EOF1

#remove escapes
sed -i 's/\\\\\#/#/g' job_temp.sh

#in the case of a restart only submit the job if it's not finished yet
if [ -f ph_q\${{q}}.out ]
then
status=\$(tail -n2 ph_q\${{q}}.out | head -n1 | awk '{{print \$2}}')
if [ ! "\$status" == "DONE." ]
then
#LSF
bsub<job_temp.sh
fi
else
#LSF
bsub<job_temp.sh
fi
cd ..
done

#update the conditions once all jobs are submitted
collect_id=\$(bjobs -J {jobname}_ph_collect | tail -n1 | awk '{{print \$1}}')
#LSF
bmod -w "{jobname}_ph_q*" \$collect_id

EOF
#ENDIF SPLIT_Q

#IF SPLIT_IRR
#if irreducible representation parallelization is enabled
elif $split_irr
then
cat > job.sh << EOF
{ph_manager_sub}
#get the number of irreducible q-points
irr_qs=\$(sed "2q;d" {pf}.dyn0 | awk '{{print $1}}')

#get the number of irreducible representations of each q-point and save them in an array
declare -a irreps
for ((i=0; i < irr_qs; i++))
do
    q=\$((i+1))
    irreps_el=\$(grep -A1 "<NUMBER_IRR_REP" _ph0/{pf}.phsave/patterns.\${{q}}.xml | tail -n1)
    irreps[\$i]=\$irreps_el
done
    
#run a phonon calculation for every irreducible representation of every irreducible q-point 
for ((q=1; q <= irr_qs; q++))
do
i=\$((q-1))
for ((r=1; r <= irreps[i]; r++))
do

#make directories for the irreducible representation and copy necessary stuff there
if [ ! -d  q\${{q}}_r\${{r}}/_ph0/{pf}.phsave ]
then
mkdir -p q\${{q}}_r\${{r}}/_ph0/{pf}.phsave
cp -r {pf}.* q\${{q}}_r\${{r}}
cp -r _ph0/{pf}.phsave/* q\${{q}}_r\${{r}}/_ph0/{pf}.phsave
fi

cd q\${{q}}_r\${{r}}
#prepare the input file
cp ../ph.in ph_q\${{q}}_r\${{r}}.in

echo "    start_q = \$q" >> ph_q\${{q}}_r\${{r}}.in
echo "    last_q = \$q" >> ph_q\${{q}}_r\${{r}}.in
echo "    start_irr = \$r" >> ph_q\${{q}}_r\${{r}}.in
echo "    last_irr = \$r" >> ph_q\${{q}}_r\${{r}}.in
echo "    /" >>  ph_q\${{q}}_r\${{r}}.in

#make the job file
cat > job_temp.sh << EOF1
{ph_q_r_sub}
\#delete the wave functions immediately after the calculation is done
if [ -f ph_q\${{q}}_r\${{r}}.out ] && [ \${{q}} -ne 1 ]
then
status=\\\\\$(tail -n2 ph_q\${{q}}_r\${{r}}.out | head -n1 | awk '{{print \\\\\$2}}')
if [ "\\\\\$status" == "DONE." ]
then 
rm _ph0/{pf}.q_\${{q}}/*.wfc*
rm {pf}.save/wfc**.dat
rm *.wfc*
fi
fi
EOF1

#remove escapes
sed -i 's/\\\\\#/#/g' job_temp.sh

#in the case of a restart only submit the job if it's not finished yet
if [ -f ph_q\${{q}}_r\${{r}}.out ]
then
status=\$(tail -n2 ph_q\${{q}}_r\${{r}}.out | head -n1 | awk '{{print \$2}}')
if [ ! "\$status" == "DONE." ]
then
#LSF
bsub<job_temp.sh
fi
else
#LSF
bsub<job_temp.sh
fi
cd ..
done
done

#update the conditions once all jobs are submitted
collect_id=\$(bjobs -J {jobname}_ph_collect | tail -n1 | awk '{{print \$1}}')
#LSF
bmod -w "{jobname}_ph_q*" \$collect_id

EOF
fi
#ENDIF SPLIT_IRR


{ph_manager_cond_sub}
fi
#ENDIF CALC_START


#collect all results
#IF CALC_START
if (($calc_start <= 7)) && (($calc_end >= 7))
then

#do nothing if the calculations are not split 
#IF NOSPLIT
if ! $split_q
then
cat > job.sh << EOF
{ph_collect_sub}
echo "nothing to collect"
EOF
#LSF
line=$(grep -n "#BSUB -w" job.sh | cut -d : -f 1)
sed -i "${{line}}s/manager/all/1" job.sh
#ENDIF NOSPLIT

#IF SPLIT_Q
elif $split_q && ! $split_irr
then
cat > job.sh << EOF
{ph_collect_sub}

#make sure that all phonon calculations have finished without problems
for ph in \$(ls ph_q*/ph_q*.out)
do
status=\$(tail -n2 \$ph | head -n1 | awk '{{print \$2}}')
if [ ! "\$status" == "DONE." ]
then
echo "Error: not all phonon calculations have finished or they ran into problems. Resubmit the phonon calculation (without initialization) before collecting."
exit
fi
done

#collect files
irr_qs=\$(sed "2q;d" {pf}.dyn0 | awk '{{print $1}}')
for ((q=1; q<=irr_qs;q++))
do

#copy the dynmats
cp q\${{q}}/_ph0/{pf}.phsave/dynmat.* _ph0/{pf}.phsave

#copy directory containing dvscf
if ((q==1))
then
cp q1/_ph0/{pf}.dvscf1 _ph0
else
cp -r q\${{q}}/_ph0/{pf}.q_\${{q}} _ph0
fi
done

#prepare input file and execute
cp ph.in ph_end.in
echo "    /" >> ph_end.in
mpirun ph.x -npool 1 -in ph_end.in > ph_end.out
EOF
#ENDIF SPLIT_Q

#IF SPLIT_IRR
elif $split_irr
then
cat > job.sh << EOF
{ph_collect_sub}

#make sure that all phonon calculations have finished without problems
for ph in \$(ls q*_r*/ph_q*_r*.out)
do
status=\$(tail -n2 \$ph | head -n1 | awk '{{print \$2}}')
if [ ! "\$status" == "DONE." ]
then
echo "Error: not all phonon calculations have finished or they ran into problems. Resubmit the phonon calculation (without initialization) before collecting."
exit
fi
done

#collect files
declare -a irreps
irr_qs=\$(sed "2q;d" {pf}.dyn0 | awk '{{print $1}}')
for ((i=0; i < irr_qs; i++))
do
    q=\$((i+1))
    irreps_el=\$(grep -A1 "<NUMBER_IRR_REP" _ph0/{pf}.phsave/patterns.\${{q}}.xml | tail -n1)
    irreps[\$i]=\$irreps_el
done

for ((q=1; q <= irr_qs; q++))
do
i=\$((q-1))
size_new=0
size_old=0
touch dvscf1_old
for ((r=1; r <= irreps[i]; r++))
do

#copy the dynmats
cp q\${{q}}_r\${{r}}/_ph0/{pf}.phsave/dynmat.* _ph0/{pf}.phsave

#combine the dvscf files bytewise
if ((q==1))
then
size_new=\$(ls -l q1_r\${{r}}/_ph0/{pf}.dvscf1 | awk '{{print \$5}}')
_count=\$((size_new - size_old))
_skip=\$((size_new - _count))
size_old=\$size_new
dd if=q1_r\${{r}}/_ph0/{pf}.dvscf1 of=dvscf1_temp skip=\$_skip count=\$_count iflag=skip_bytes,count_bytes
cat dvscf1_old dvscf1_temp > dvscf1_new
mv dvscf1_new dvscf1_old

else
size_new=\$(ls -l q\${{q}}_r\${{r}}/_ph0/{pf}.q_\${{q}}/{pf}.dvscf1 | awk '{{print \$5}}')
_count=\$((size_new - size_old))
_skip=\$((size_new - _count))
size_old=\$size_new
dd if=q\${{q}}_r\${{r}}/_ph0/{pf}.q_\${{q}}/{pf}.dvscf1 of=dvscf1_temp skip=\$_skip count=\$_count iflag=skip_bytes,count_bytes
cat dvscf1_old dvscf1_temp > dvscf1_new
mv dvscf1_new dvscf1_old
fi
done

#move the combined dvscf files to the right location
if ((q==1))
then
mv dvscf1_old _ph0/{pf}.dvscf1
else
if [ ! -d _ph0/{pf}.q_\${{q}} ]
then
mkdir _ph0/{pf}.q_\${{q}}
fi
mv dvscf1_old _ph0/{pf}.q_\${{q}}/{pf}.dvscf1
fi
done

#prepare input file and execute
cp ph.in ph_end.in
echo "    /" >> ph_end.in
mpirun ph.x -npool 1 -in ph_end.in > ph_end.out

EOF
fi
#ENDIF SPLIT_IRR

{ph_collect_cond_sub}
fi
#ENDIF CALC_START

#q2r
if (($calc_start <= 8)) && (($calc_end >= 8))
then

cat > job.sh << EOF
{q2r_sub}
if [ -f {pf}.dyn1.xml ]
then
    cp {pf}.dyn0 {pf}.dyn0.xml
    sed -i "s/{pf}\.dyn/{pf}\.dyn\.xml/g" q2r.in
fi
mpirun q2r.x -npool 1 -in q2r.in > q2r.out
EOF

{q2r_cond_sub}
fi

#matdyn
if (($calc_start <= 9)) && (($calc_end >= 9))
then

cat > job.sh << EOF
{matdyn_sub}
if [ -f {pf}.dyn1.xml ]
then
    sed -i "s/{pf}\.fc/{pf}\.fc\.xml/g" matdyn.in
fi
mpirun matdyn.x -npool 1 -in matdyn.in > matdyn.out
EOF

{matdyn_cond_sub}
fi

#tidying up
if (($calc_start == 10)) && (($calc_end == 10))
then
cat > job.sh << EOF
{tidy_sub}
#phonon directory
cd $base_dir/PHB
rm -r q**/
rm *.xml
rm *.e
rm *.o
rm job*
rm *dvscf*
rm *wfc*
rm matdyn.modes
#electron directory
cd $base_dir/ELB
rm -r *.save/
rm *.xml
rm *.e
rm *.o
rm job*
rm *wfc*
EOF

if ! $no_sub
then
bsub<job.sh
fi
fi
'''.format(elb_scf_sub = make_job_sub(jobname + '_elb_scf',num_of_cpu_scf,ram,scf_t,'scf.in','scf.out','pw.x',''),
           elb_scf_cond_sub = check_cond_sub(1, True),
           elb_sub = make_job_sub(jobname + '_elb',num_of_cpu_scf,ram,scf_t,'bands.in','bands.out','pw.x',jobname + '_elb_scf'),
           elb_cond_sub = check_cond_sub(2),
           elb_ip_sub = make_job_sub(jobname + '_elb_ip',num_of_cpu_scf,ram,4,'bands_ip.in','bands_ip.out','bands.x',jobname + '_elb'),
           elb_ip_cond_sub = check_cond_sub(3),
           ph_scf_sub = make_job_sub(jobname + '_ph_scf',num_of_cpu_scf,ram,scf_t,'scf.in','scf.out','pw.x',''),
           ph_scf_cond_sub = check_cond_sub(4, True),
           ph_init_sub = make_job_sub(jobname + '_ph_init',1,ram,4,'','','',jobname + '_ph_scf'),
           ph_init_cond_sub = check_cond_sub(5),
           ph_all_sub = make_job_sub(jobname + '_ph_all',num_of_cpu_ph,ram,q_t,'ph_all.in','ph_all.out','ph.x',jobname + '_ph_init'),
           ph_manager_sub = make_job_sub(jobname + '_ph_manager',1,ram,4,'','','',jobname + '_ph_init'),
           ph_manager_cond_sub = check_cond_sub(6),
           ph_q_sub = make_job_sub(jobname + '_ph_q\${q}',num_of_cpu_ph,ram,q_t,'ph_q\${q}.in','ph_q\${q}.out','ph.x','',True),
           ph_q_r_sub = make_job_sub(jobname + '_ph_q\${q}_r\${r}',num_of_cpu_ph,ram,q_t,'ph_q\${q}_r\${r}.in','ph_q\${q}_r\${r}.out','ph.x','',True),
           ph_collect_sub = make_job_sub(jobname + '_ph_collect',1,ram,4,'','','',jobname + '_ph_manager'),
           ph_collect_cond_sub = check_cond_sub(7),
           q2r_sub = make_job_sub(jobname + '_q2r',1,ram,4,'','','',jobname + '_ph_collect'),
           q2r_cond_sub = check_cond_sub(8),
           matdyn_sub = make_job_sub(jobname + '_matdyn',1,ram,4,'','','',jobname + '_q2r'),
           matdyn_cond_sub = check_cond_sub(9),
           tidy_sub = make_job_sub(jobname + '_tidy',1,ram,4,'','','',''),
           module_commands = generate_modules(modules),
           split_q = split_q, split_irr = split_irr, pf = pf, scf_t = scf_t, q_t = q_t, jobname = jobname, num_of_cpu_ph = num_of_cpu_ph)]

#submission script for all calculations regarding electron-phonon coupling
epw_sh = ['''
#!/bin/bash
#Execute from the parent directory (containing ELB, PHB, EPM, ISO directories).

#Specify the number of bands to be wannierized, the index of the highest non-wannierized band and the initial guess...
#...for the wannier functions. If you change this to get the windows in #3 right, you have to run #1 (with no_sub = true) again ...
#...in order to reset the corresponding wannierization windows again.
num_of_wan={num_of_wan}
wannier_init={projections_epw}

#Specify if the wannierization energy windows should be automatically determined. If you set it to false, specify the...
#...windows in the next block. The inner window will be the largest possible window that contains only the specified bands...
#...but no parts of other bands. The outer window will be the smallest possible window that contains the whole specified...
#bands including other bands should they fall in this window. You can narrow the inner and outer window with dE which specifies...
#the value in eV that gets added/subtracted to the bottom/top of the inner window (i.e. negative values will broaden it).
auto_window={auto_window}
dE_inner=0.0
dE_outer=0.0


#Set the inner and outer wannierization energy windows manually (in eV). You only need to specify the values of this block...
#... if you disabled auto_window or if the automatic determination didn't work.
wan_min={wan_min}
inner_bottom={inner_bottom}
inner_top={inner_top}
outer_bottom={outer_bottom}
outer_top={outer_top}

#Specify the starting and ending point of calculation. 
#1: input preparation and scf calculation (set no_sub to true for just the input preparation)
#2: nscf calculation for phonon linewidths
#3: wannierization of the electron bands
#Note it is a good idea to stop the calculation here and check the wannierized bands (i.e. compare them to the...
#...Bloch bands). You can use #11 to do that.
#4: calculation of electron-phonon coupling coefficients for phonon linewidths
#5: phonon linewidth calculation along high symmetry lines i.e. bands 
#6: phonon linewidth and a2F calculation in fine grid 
#7: scf calculation for solving the isotropic Eliashberg functions (Tc)
#8: nscf caclulation for Tc
#9: solving the isotropic Eliashberg functions
#(#11: Plot the Wannier and Bloch bands using gnuplot on the cluster (make sure x11 forwarding is enabled)...
#     ...This won't submit a job but simply plot the Wannier bands obtained after #3 and compare them to the...
#     ...Bloch bands obtained through pw.x (located in ELB). Furthermore the energy windows will also be plotted...
#     ...Usable only as a separate (calc_start = calc_end) execution.
#(#12: tidy up - deletes everything but the obtained end results. This is not necessary but saves space...
#     ...Usable only as a seperate (calc_start = calc_end) execution. Only use once you are sure that everything worked.)

calc_start=1
calc_end=10

#If you need a specific submission script (e.g. for wannierization) set calc_start = calc_end (e.g. = 3) and set no_sub
#to true. You'll find the job script as job.sh in the corresponding directory.
no_sub=false

#load modules 
{module_commands}

#________________________________EXECUTE_THE_SCRIPT________________________________#
#__________________________________________________________________________________#

ref_bands={ref_bands}
base_dir=$(pwd)
ref_dir="{bands_dir}"
if ! $ref_bands
then
ref_dir=$base_dir
fi

#__________POSTPROCESSING__________#

if [ ! -d $ref_dir/PHB/save ]
then
cd $ref_dir/PHB
python pp.py
cd $base_dir
fi

#__________PREPARE_INPUT___________#
#IF CALC_START
if (($calc_start == 1))
then
cp $base_dir/ELB/scf.in $base_dir/EPM

cp $base_dir/EPM/nscf.in_orig $base_dir/EPM/nscf.in
cp $base_dir/EPM/wannier.in_orig $base_dir/EPM/wannier.in
cp $base_dir/EPM/wannier_epm.in_orig $base_dir/EPM/wannier_epm.in
cp $base_dir/EPM/ph_lw.in_orig $base_dir/EPM/ph_lw.in
cp $base_dir/EPM/a2F.in_orig $base_dir/EPM/a2F.in
cp $base_dir/ISO/eliashberg_iso.in_orig $base_dir/ISO/eliashberg_iso.in

#generate and append uniform k-grid for nscf.in
python $base_dir/EPM/kmesh.py {k_x} {k_y} {k_z} 1 >> $base_dir/EPM/nscf.in

#append irreducible q-points to epw input files where necessary
num_of_irr_k=$(sed "2q;d" $ref_dir/PHB/*.dyn0 | awk '{{print $1}}')
echo "$num_of_irr_k cartesian" >> $base_dir/EPM/wannier_epm.in
tail -n +3 $ref_dir/PHB/*.dyn0 >> $base_dir/EPM/wannier_epm.in
echo "$num_of_irr_k cartesian" >> $base_dir/EPM/ph_lw.in
tail -n +3 $ref_dir/PHB/*.dyn0 >> $base_dir/EPM/ph_lw.in
echo "$num_of_irr_k cartesian" >> $base_dir/EPM/a2F.in
tail -n +3 $ref_dir/PHB/*.dyn0 >> $base_dir/EPM/a2F.in
echo "$num_of_irr_k cartesian" >> $base_dir/ISO/eliashberg_iso.in
tail -n +3 $ref_dir/PHB/*.dyn0 >> $base_dir/ISO/eliashberg_iso.in


#get and set the Fermi energy where necessary
Ef=$(grep Fermi $ref_dir/ELB/scf.out | awk '{{print $5}}')
line=$(grep -n fermi_energy $base_dir/EPM/ph_lw.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$Ef/3" $base_dir/EPM/ph_lw.in
line=$(grep -n fermi_energy $base_dir/EPM/a2F.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$Ef/3" $base_dir/EPM/a2F.in
line=$(grep -n fermi_energy $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$Ef/3" $base_dir/ISO/eliashberg_iso.in


#copy the updated and necessary input files to the ISO directory
cp $base_dir/EPM/scf.in $base_dir/ISO
cp $base_dir/EPM/nscf.in $base_dir/ISO

#if enabled, get the wannierization energy windows automatically from the electron band calculation
#and set the values in the input files where necessary

#IF AUTO_WINDOW
if $auto_window
then

#get the number of bands
bands=$(head -n1 $ref_dir/ELB/*bands.dat | awk '{{print $3}}' | sed "s/,//")
#get the number of k-points in the total band path
points=$(head -n1 $ref_dir/ELB/*bands.dat | awk '{{print $5}}')

#make a temporary file that stores the band energies sequentially per band
(cat $ref_dir/ELB/*bands.dat.gnu | sed 's/[^ ]*[^ ]//1' | sed 's/\ //g' | sed '/^$/d') > $ref_dir/ELB/bands_temp

#run the py script to find the bands
windows=$(python $ref_dir/ELB/wannier_windows.py $num_of_wan $bands $points $Ef $ref_dir/ELB/bands_temp)
smallest_inner=$(echo $windows | awk '{{print $1}}')
largest_inner=$(echo $windows | awk '{{print $2}}')
smallest_outer=$(echo $windows | awk '{{print $3}}')
largest_outer=$(echo $windows | awk '{{print $4}}')
wan_min=$(echo $windows | awk '{{print $5}}')

#narrow the bands
smallest_inner=$(echo "$smallest_inner + $dE_inner" | bc -l)
largest_inner=$(echo "$largest_inner - $dE_inner" | bc -l)
smallest_outer=$(echo "$smallest_outer + $dE_outer" | bc -l)
largest_outer=$(echo "$largest_outer - $dE_outer" | bc -l)

inner_bottom=$smallest_inner
inner_top=$largest_inner
outer_bottom=$smallest_outer
outer_top=$largest_outer

fi
#ENDIF AUTO_WINDOW

#set the wannierization parameters
line=$(grep -n nbndsub $base_dir/EPM/wannier.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$num_of_wan/3" $base_dir/EPM/wannier.in
line=$(grep -n nbndskip $base_dir/EPM/wannier.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$wan_min/3" $base_dir/EPM/wannier.in

index=0
for p in "${{wannier_init[@]}}"
do
    ((index++))
    line=$(grep -n "proj(${{index}})" $base_dir/EPM/wannier.in | cut -d : -f 1)
    sed -i -e "${{line}}s/[^ ]*[^ ]/'$p'/3" $base_dir/EPM/wannier.in
done

line=$(grep -n nbndsub $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$num_of_wan/3" $base_dir/EPM/wannier_epm.in
line=$(grep -n nbndskip $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$wan_min/3" $base_dir/EPM/wannier_epm.in

index=0
for p in "${{wannier_init[@]}}"
do
    ((index++))
    line=$(grep -n "proj(${{index}})" $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
    sed -i -e "${{line}}s/[^ ]*[^ ]/'$p'/3" $base_dir/EPM/wannier_epm.in
done

line=$(grep -n nbndsub $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$num_of_wan/3" $base_dir/ISO/eliashberg_iso.in
line=$(grep -n nbndskip $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$wan_min/3" $base_dir/ISO/eliashberg_iso.in

index=0
for p in "${{wannier_init[@]}}"
do
    ((index++))
    line=$(grep -n "proj(${{index}})" $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
    sed -i -e "${{line}}s/[^ ]*[^ ]/'$p'/3" $base_dir/ISO/eliashberg_iso.in
done


line=$(grep -n dis_win_min $base_dir/EPM/wannier.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$outer_bottom/3" $base_dir/EPM/wannier.in
line=$(grep -n dis_win_max $base_dir/EPM/wannier.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$outer_top/3" $base_dir/EPM/wannier.in
line=$(grep -n dis_froz_min $base_dir/EPM/wannier.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$inner_bottom/3" $base_dir/EPM/wannier.in
line=$(grep -n dis_froz_max $base_dir/EPM/wannier.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$inner_top/3" $base_dir/EPM/wannier.in

line=$(grep -n dis_win_min $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$outer_bottom/3" $base_dir/EPM/wannier_epm.in
line=$(grep -n dis_win_max $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$outer_top/3" $base_dir/EPM/wannier_epm.in
line=$(grep -n dis_froz_min $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$inner_bottom/3" $base_dir/EPM/wannier_epm.in
line=$(grep -n dis_froz_max $base_dir/EPM/wannier_epm.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$inner_top/3" $base_dir/EPM/wannier_epm.in

line=$(grep -n dis_win_min $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$outer_bottom/3" $base_dir/ISO/eliashberg_iso.in
line=$(grep -n dis_win_max $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$outer_top/3" $base_dir/ISO/eliashberg_iso.in
line=$(grep -n dis_froz_min $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$inner_bottom/3" $base_dir/ISO/eliashberg_iso.in
line=$(grep -n dis_froz_max $base_dir/ISO/eliashberg_iso.in | cut -d : -f 1)
sed -i -e "${{line}}s/[^ ]*[^ ]/$inner_top/3" $base_dir/ISO/eliashberg_iso.in

fi
#ENDIF CALC_START


#__________JOB_SUBMISSION___________#
if ( (($calc_start == 11)) && ((! $calc_end == 11)) ) || ( ((! $calc_start == 11)) && (($calc_end == 11)) )
then
echo "Error: plotting the Bloch and Wannier bands is only possible as a separate execution"
exit
elif (($calc_start == 11)) && (($calc_start == $calc_end))
then

#get the positions of the last k-points of the bands
last_QE=$(tail -n2 $ref_dir/ELB/*bands.dat.gnu | head -n1 | awk '{{print $1}}')
last_EPM=$(tail -n2 $base_dir/EPM/{pf}_band.dat | head -n1 | awk '{{print $1}}')

#get the windows
smallest_inner=$(grep dis_froz_min $base_dir/EPM/wannier_epm.in | awk '{{print $3}}')
largest_inner=$(grep dis_froz_max $base_dir/EPM/wannier_epm.in | awk '{{print $3}}')
smallest_outer=$(grep dis_win_min $base_dir/EPM/wannier_epm.in | awk '{{print $3}}')
largest_outer=$(grep dis_win_max $base_dir/EPM/wannier_epm.in | awk '{{print $3}}')

#make the gnuplot script
cat > $base_dir/ELB/check_bands.gnu << EOF
reset
set terminal x11 dashed title "{pf}"
set grid
set linestyle 1 lc rgb "black" lw 2
set linestyle 2 lc rgb "red" lw 1
set linestyle 3 lc rgb "blue" lw 2
set linestyle 4 lc rgb "green" lw 2

y_low = $smallest_outer - 1
y_high = $largest_outer + 1
set yrange [y_low:y_high]

set ylabel "[eV]"

set arrow from graph 0, first $smallest_outer to graph 1, first $smallest_outer ls 3 nohead
set arrow from graph 0, first $largest_outer to graph 1, first $largest_outer ls 3 nohead
set arrow from graph 0, first $smallest_inner to graph 1, first $smallest_inner ls 4 nohead
set arrow from graph 0, first $largest_inner to graph 1, first $largest_inner ls 4 nohead

set arrow from graph 0.95, first $smallest_outer to graph 0.95, first $largest_outer ls 3 heads
set arrow from graph 0.9, first $smallest_inner to graph 0.9, first $largest_inner ls 4 heads

set label "outer window" at graph 0.96, graph 0.5 rotate by 90 tc lt 3  font "Helvetica, 12"
set label "inner window" at graph 0.91, graph 0.5 rotate by 90 tc lt 2 font "Helvetica, 12"
blochf=system("ls $ref_dir/ELB/*bands.dat.gnu")
plot blochf u (\$1/$last_QE):2 w l ls 1 t "Bloch", \\
     "$base_dir/EPM/{pf}_band.dat" u (\$1/$last_EPM):2 w l ls 2 t "Wannier"
EOF
gnuplot -p "$base_dir/ELB/check_bands.gnu" 
fi

if ( (($calc_start == 12)) && ((! $calc_end == 12)) ) || ( ((! $calc_start == 12)) && (($calc_end == 12)) )
then
echo "Error: please only tidy up once you are sure that all calculations have finished correctly"
exit
fi

if (($calc_start > $calc_end))
then
echo "Error: calculation order makes no sense"
exit
fi
#EPM calculations
cd $base_dir/EPM

#scf
if (($calc_start <= 1)) && (($calc_end >= 1))
then
cat > job.sh << EOF
{epm_scf_sub}
#get the number of bands from the scf.out file and put them in the nscf.in file
num_of_bands=\$(grep nbnd nscf.in | awk '{{print \$3}}')
if [ \$num_of_bands -eq 0 ]
then
    num_of_bands=\$(grep "number of Kohn-Sham states" scf.out | awk '{{print \$5}}')
    line=\$(grep -n nbnd nscf.in | cut -d : -f 1)
    sed -i "\${{line}}s/[^ ]*[^ ]/\${{num_of_bands}}/3" nscf.in
fi
EOF

{epm_scf_cond_sub}
fi

#nscf
if (($calc_start <= 2)) && (($calc_end >= 2))
then

cat > job.sh << EOF
{epm_nscf_sub}
EOF

{epm_nscf_cond_sub}
fi

#wannierization
if (($calc_start <= 3)) && (($calc_end >= 3))
then

cat > job.sh << EOF
{wannier_sub}
EOF

{wannier_cond_sub}
fi

#electron phonon matrix
if (($calc_start <= 4)) && (($calc_end >= 4))
then

cat > job.sh << EOF
{epm_sub}
EOF

{epm_cond_sub}
fi

#phonon linewidth
if (($calc_start <= 5)) && (($calc_end >= 5))
then

cat > job.sh << EOF
{ph_lw_sub}
mv linewidth.phself linewidth.phself_bands
EOF

{ph_lw_cond_sub}
fi

#a2F
if (($calc_start <= 6)) && (($calc_end >= 6))
then

cat > job.sh << EOF
{a2F_sub}
mv linewidth.phself linewidth.phself_grid
EOF

{a2F_cond_sub}
fi

#ISO calculations
cd $base_dir/ISO
#scf
if (($calc_start <= 7)) && (($calc_end >= 7))
then

cat > job.sh << EOF
{iso_scf_sub}
#get the number of bands from the scf.out file and put them in the nscf.in file
num_of_bands=\$(grep nbnd nscf.in | awk '{{print \$3}}')
if [ \$num_of_bands -eq 0 ]
then
    num_of_bands=\$(grep "number of Kohn-Sham states" scf.out | awk '{{print \$5}}')
    line=\$(grep -n nbnd nscf.in | cut -d : -f 1)
    sed -i "\${{line}}s/[^ ]*[^ ]/\${{num_of_bands}}/3" nscf.in
fi
EOF

{iso_scf_cond_sub}
fi

#nscf
if (($calc_start <= 8)) && (($calc_end >= 8))
then

cat > job.sh << EOF
{iso_nscf_sub}
EOF

{iso_nscf_cond_sub}
fi

#isotropic eliashberg
if (($calc_start <= 9)) && (($calc_end >= 9))
then

cat > job.sh << EOF
{iso_sub}
EOF

{iso_cond_sub}
fi

#tidying up
if (($calc_start == 12)) && (($calc_end == 12))
then
cat > job.sh << EOF
{tidy_sub}
#phonon directory
cd $base_dir/PHB
rm -r _ph0 *.save
rm *.dyn* *.e *.o *.fc
#EPM directory
cd $base_dir/EPM
rm -r *save
rm *.e *.o
rm *.epb* *wfc* *.xml *.wout *.win *.ukk *.nnkp *.kmap *.epmatwe1 *.chk decay* *.fmt

EOF
if ! $no_sub
then
bsub<job.sh
fi

fi
'''.format(epm_scf_sub = make_job_sub(jobname + '_epm_scf',num_of_cpu_scf,ram,scf_t,'scf.in','scf.out','pw.x',''),
           epm_scf_cond_sub = check_cond_sub(1, True),
           epm_nscf_sub = make_job_sub(jobname + '_epm_nscf',num_of_cpu_scf,ram,scf_t,'nscf.in','nscf.out','pw.x',jobname + '_epm_scf'),
           epm_nscf_cond_sub = check_cond_sub(2),
           wannier_sub = make_job_sub(jobname + '_wannier',num_of_cpu_scf,ram,4,'wannier.in','wannier.out','epw.x',jobname + '_epm_nscf'),
           wannier_cond_sub = check_cond_sub(3),
           epm_sub = make_job_sub(jobname + '_epm',num_of_cpu_epw,ram,epw_t,'wannier_epm.in','wannier_epm.out','epw.x',jobname + '_wannier'),
           epm_cond_sub = check_cond_sub(4),
           ph_lw_sub = make_job_sub(jobname + '_ph_lw',num_of_cpu_epw,ram,epw_t,'ph_lw.in','ph_lw.out','epw.x',jobname + '_epm'),
           ph_lw_cond_sub = check_cond_sub(5),
           a2F_sub = make_job_sub(jobname + '_a2F',num_of_cpu_epw,ram,epw_t,'a2F.in','a2F.out','epw.x',jobname + '_ph_lw'),
           a2F_cond_sub = check_cond_sub(6),
           iso_scf_sub = make_job_sub(jobname + '_iso_scf',num_of_cpu_scf,ram,scf_t,'scf.in','scf.out','pw.x',''),
           iso_scf_cond_sub = check_cond_sub(7, True),
           iso_nscf_sub = make_job_sub(jobname + '_iso_nscf',num_of_cpu_scf,ram,scf_t,'nscf.in','nscf.out','pw.x',jobname + '_iso_scf'),
           iso_nscf_cond_sub = check_cond_sub(8),
           iso_sub = make_job_sub(jobname + '_iso',num_of_cpu_epw,ram,epw_t,'eliashberg_iso.in','eliashberg_iso.out','epw.x',jobname + '_iso_nscf'),
           iso_cond_sub = check_cond_sub(9),
           tidy_sub = make_job_sub(jobname + '_tidy',1,ram,4,'','','',''),
           module_commands = generate_modules(modules),
           num_of_wan = num_of_wan, wan_min = wan_min, projections_epw = projections_epw, auto_window = auto_window,
           inner_bottom = inner_bottom, inner_top = inner_top, outer_bottom = outer_bottom, outer_top = outer_top,
           ref_bands = ref_bands, bands_dir = bands_dir, k_x = k_x, k_y = k_y, k_z = k_z, pf = pf)]

#_________________________________________________________________________________________#
#_________________________________________________________________________________________#

#___________________________________AUXILIARY_SCRIPTS_____________________________________#
#_________________________________________________________________________________________#

#generating uniform k-mesh. Taken from the wannier90 repository on github (https://github.com/wannier-developers/wannier90/tree/develop/utility)
#and rewritten in python
kmesh_py = ['''
from argparse import ArgumentParser
import sys
parser = ArgumentParser()
parser.add_argument('grid_size', type=int, nargs=3)
parser.add_argument('weighted', type=int, default=0, nargs='?')
args = parser.parse_args()

n1 = args.grid_size[0]
n2 = args.grid_size[1]
n3 = args.grid_size[2]
weighted = args.weighted

if n1 <= 0:
    sys.exit("n1 needs to be > 0")
if n2 <= 0:
    sys.exit("n2 needs to be > 0")
if n3 <= 0:
    sys.exit("n3 needs to be > 0")
  
totpts = n1 * n2 * n3

if not weighted:
    print("K_POINTS crystal")
    print("%d" % (totpts))
    for i in range(0, n1):
        for j in range(0, n2):
            for k in range(0, n3):
                print("%12.8f%12.8f%12.8f" % (i/n1,j/n2,k/n3))

else:
    print("K_POINTS crystal")
    print("%d" % (totpts))
    for i in range(0, n1):
        for j in range(0, n2):
            for k in range(0, n3):
                print("%12.8f%12.8f%12.8f%14.6e" % (i/n1,j/n2,k/n3,1/totpts))
''']

#postprocessing of the phonon calculations for EPW calculations. Taken from the EPW repository on github
#(https://github.com/QEF/q-e/tree/master/EPW/bin) and slightly changed
pp_py = ['''    
#!/usr/bin/python
#
# Post-processing script from of PH data in format used by EPW
# 14/07/2015 - Creation of the script - Samuel Ponce
# 14/03/2018 - Automatically reads the number of q-points - Michael Waters
# 14/03/2018 - Detect if SOC is included in the calculation - Samuel Ponce 
# 13/11/2018 - Write dyn files in xml format for SOC case - Shunhong Zhang (USTC)
# 
import numpy as np
import os
from xml.dom import minidom

# Convert the dyn files to the xml form, for SOC case - Shunhong Zhang (USTC)
def dyn2xml(prefix):
    ndyn=int(os.popen('head -2 {{0}}.dyn0|tail -1'.format(prefix)).read())
    for idyn in range(1,ndyn+1):
        print('{{0}}.dyn{{1}} to {{0}}.dyn_q{{1}}.xml'.format(prefix,idyn))
        dynmat=dyn(prefix,idyn)
        dynmat._write_xml()
def get_geom_info():
    phfile = 'ph_end.out'
    if os.path.isfile(phfile)==False:
       print('cannot extract geometry info from ' + phfile)
       return 1
    else:
       volm=float(os.popen('grep -a volume ' + phfile + ' 2>/dev/null|tail -1').readline().split()[-2])
       get_at=os.popen('grep -a -A 3 "crystal axes" ' + phfile + ' 2>/dev/null|tail -3').readlines()
       at=np.array([[float(item) for item in line.split()[3:6]] for line in get_at])
       get_bg=os.popen('grep -a -A 3 "reciprocal axes" ' + phfile + ' 2>/dev/null|tail -3').readlines()
       bg=np.array([[float(item) for item in line.split()[3:6]] for line in get_bg])
       return volm,at,bg

class dyn(object):
    def __init__(self,prefix,idyn):
        self._prefix=prefix
        self._idyn=idyn
        fil='{{0}}.dyn{{1}}'.format(prefix,idyn)
        f=open(fil)
        self._comment=f.readline()
        f.readline()
        line=f.readline().split()
        self._ntype=int(line[0])
        self._natom=int(line[1])
        self._ibrav=int(line[2])
        self._nspin=1
        self._cell_dim=np.array([float(ii) for ii in line[3:]])
        self._volm=0
        self._at=np.zeros((3,3),float)
        self._bg=np.zeros((3,3),float)
        try: self._volm,self._at,self._bg = get_geom_info()
        except: print('warning: lattice info not found')
        for i in range(0,4):
            f.readline()
        self._species=[];
        self._mass=[]
        for i in range(self._ntype):
            line=f.readline().split()
            self._species.append(line[1].strip("'"))
            self._mass.append(float(line[-1])/911.4442)  # normalize to atomic mass
        self._atom_type=np.zeros(self._natom,int)
        self._pos=np.zeros((self._natom,3),float)
        for i in range(self._natom):
            line=f.readline().split()
            self._atom_type[i]=int(line[1])
            for j in range(3): self._pos[i,j]=float(line[j+2])
        self._nqpt=int(os.popen('grep -c "Dynamical  Matrix" {{0}}'.format(fil)).read().split()[0])
        self._qpt=[]
        self._dynmat=np.zeros((self._nqpt,self._natom,self._natom,3,3,2),float)
        f.readline()
        for iqpt in range(self._nqpt):
            f.readline();
            f.readline()
            line=f.readline().split()
            self._qpt.append(np.array([float(item) for item in line[3:6]]))
            f.readline()
            for i in range(self._natom):
                for j in range(self._natom):
                    f.readline()
                    data=np.fromfile(f,sep=' ',count=18,dtype=float).reshape(3,3,2)
                    self._dynmat[iqpt,i,j]=data
        self._qpt=np.array(self._qpt)
        for i in range(5): f.readline()
        self._freq=np.zeros((self._natom*3,2),float)
        self._disp=np.zeros((self._natom*3,self._natom,3,2),float)
        for i in range(self._natom*3):
            line=f.readline().split()
            self._freq[i,0]=float(line[4])
            self._freq[i,1]=float(line[7])
            for j in range(self._natom):
                line=f.readline().split()[1:-1]
                data=np.array([float(item) for item in line]).reshape(3,2)
                self._disp[i,j]=data

    def _write_xml(self):
        doc=minidom.Document()
        root = doc.createElement('Root')
        doc.appendChild(root)
        geom_info=doc.createElement('GEOMETRY_INFO')
        tags=('NUMBER_OF_TYPES','NUMBER_OF_ATOMS','BRAVAIS_LATTICE_INDEX','SPIN_COMPONENTS')
        numbers=(self._ntype,self._natom,self._ibrav,self._nspin)
        for i,(tag,num) in enumerate(zip(tags,numbers)):
            inode=doc.createElement(tag)
            inode.setAttribute('type','integer')
            inode.setAttribute('size','1')
            inode.text=num
            inode.appendChild(doc.createTextNode(str(num)))
            geom_info.appendChild(inode)
        cell_dim=doc.createElement('CELL_DIMENSIONS')
        cell_dim.setAttribute('type','real')
        cell_dim.setAttribute('size','6')
        for i in range(6):
            cell_dim.appendChild(doc.createTextNode('{{0:16.10f}}'.format(self._cell_dim[i])))
        geom_info.appendChild(cell_dim)
        tags=['AT','BG']
        for tag,lat in zip(tags,(self._at,self._bg)):
            inode=doc.createElement(tag)
            inode.setAttribute('type','real')
            inode.setAttribute('size','9')
            inode.setAttribute('columns','3')
            for i in range(3):
                text=' '.join(['{{0:16.10f}}'.format(item) for item in lat[i]])
                inode.appendChild(doc.createTextNode(text))
            geom_info.appendChild(inode)
        volm=doc.createElement('UNIT_CELL_VOLUME_AU')
        volm.setAttribute('type','real')
        volm.setAttribute('size','1')
        volm.appendChild(doc.createTextNode('{{0:16.10f}}'.format(self._volm)))
        geom_info.appendChild(volm)
        for itype in range(self._ntype):
            nt=doc.createElement('TYPE_NAME.{{0}}'.format(itype+1))
            nt.setAttribute('type','character')
            nt.setAttribute('size','1')
            nt.setAttribute('len','3')
            nt.appendChild(doc.createTextNode('{{0}}'.format(self._species[itype])))
            na=doc.createElement('MASS.{{0}}'.format(itype+1))
            na.setAttribute('type','real')
            na.setAttribute('size','1')
            na.appendChild(doc.createTextNode('{{0:16.10f}}'.format(self._mass[itype])))
            geom_info.appendChild(nt)
            geom_info.appendChild(na)
        for iat in range(self._natom):
            at=doc.createElement('ATOM.{{0}}'.format(iat+1))
            at.setAttribute('SPECIES','{{0}}'.format(self._species[self._atom_type[iat]-1]))
            at.setAttribute('INDEX',str(iat+1))
            pos=' '.join(['{{0:16.10f}}'.format(item) for item in self._pos[iat]])
            at.setAttribute('TAU',pos)
            geom_info.appendChild(at)
        nqpt=doc.createElement('NUMBER_OF_Q')
        nqpt.setAttribute('type','integer')
        nqpt.setAttribute('size','1')
        nqpt.appendChild(doc.createTextNode(str(self._nqpt)))
        geom_info.appendChild(nqpt)
        root.appendChild(geom_info)
        for iqpt in range(self._nqpt):
            dynmat=doc.createElement('DYNAMICAL_MAT_.{{0}}'.format(iqpt+1))
            qpt=doc.createElement('Q_POINT')
            qpt.setAttribute('type','real')
            qpt.setAttribute('size','3')
            qpt.setAttribute('columns','3')
            tnode=doc.createTextNode(' '.join(['{{0:16.10f}}'.format(item) for item in self._qpt[iqpt]]))
            qpt.appendChild(tnode)
            dynmat.appendChild(qpt)
            for iat in range(self._natom):
                for jat in range(self._natom):
                    ph=doc.createElement('PHI.{{0}}.{{1}}'.format(iat+1,jat+1))
                    ph.setAttribute('type','complex')
                    ph.setAttribute('size','9')
                    ph.setAttribute('columns','3')
                    for i in range(3):
                        for j in range(3):
                            text='{{0:16.10f}} {{1:16.10f}}'.format(self._dynmat[iqpt,iat,jat,i,j,0],self._dynmat[iqpt,iat,jat,i,j,1])
                            ph.appendChild(doc.createTextNode(text))
                    dynmat.appendChild(ph)
            root.appendChild(dynmat)
        mode=doc.createElement('FREQUENCIES_THZ_CMM1')
        for iomega in range(self._natom*3):
            inode=doc.createElement('OMEGA.{{0}}'.format(iomega+1))
            inode.setAttribute('type','real')
            inode.setAttribute('size','2')
            inode.setAttribute('columns','2')
            inode.appendChild(doc.createTextNode('{{0:16.10f}} {{1:16.10f}}'.format(self._freq[iomega,0],self._freq[iomega,1])))
            idisp=doc.createElement('DISPLACEMENT.{{0}}'.format(iomega+1))
            idisp.setAttribute('tpye','complex')
            idisp.setAttribute('size','3')
            for iat in range(self._natom):
                for j in range(3):
                    tnode=doc.createTextNode('{{0:16.10f}} {{1:16.10f}}'.format(self._disp[iomega,iat,j,0],self._disp[iomega,iat,j,1]))
                    idisp.appendChild(tnode)
            mode.appendChild(inode)
            mode.appendChild(idisp)
        root.appendChild(mode)
        fp = open('{{0}}.dyn{{1}}.xml'.format(self._prefix,self._idyn), 'w')
        doc.writexml(fp, addindent='  ', newl='\\n')

# Return the number of q-points in the IBZ
def get_nqpt(prefix):
  fname = '_ph0/' +prefix+'.phsave/control_ph.xml'

  fid = open(fname,'r')
  lines = fid.readlines() # these files are relatively small so reading the whole thing shouldn't be an issue
  fid.close()

  line_number_of_nqpt = 0
  while 'NUMBER_OF_Q_POINTS' not in lines[line_number_of_nqpt]: # increment to line of interest
    line_number_of_nqpt +=1
  line_number_of_nqpt +=1 # its on the next line after that text

  nqpt = int(lines[line_number_of_nqpt])

  return nqpt

# Check if the calculation include SOC
def hasSOC(prefix):
  fname = prefix+'.save/data-file-schema.xml'

  xmldoc = minidom.parse(fname)
  item = xmldoc.getElementsByTagName('spinorbit')[0]
  lSOC = item.childNodes[0].data
  
  return lSOC

# Check if the calculation was done in sequential
def isSEQ(prefix):
  fname = '_ph0/'+str(prefix)+'.dvscf'
  if (os.path.isfile(fname)):
    lseq = True
  else:
    lseq = False
 
  return lseq
    
# Enter the number of irr. q-points
prefix = '{pf}' #get the prefix directly from this generation script

# Test if SOC
SOC = hasSOC(prefix)

# If SOC detected, but dyn is not in XML always convert it 
if SOC=='true' and not os.path.isfile(prefix + '.dyn1.xml'):
    dyn2xml(prefix)

# If no SOC detected, do nothing 

# Test if seq. or parallel run
SEQ = isSEQ(prefix)

nqpt =  get_nqpt(prefix)

if not os.path.isdir('save'): 
    os.system('mkdir save')

for iqpt in np.arange(1,nqpt+1):
  label = str(iqpt)

  # Case calculation in seq.
  if SEQ:
    # Case with SOC
    if SOC == 'true':
      os.system('cp '+prefix+'.dyn0 '+prefix+'.dyn0.xml')
      os.system('cp '+prefix+'.dyn'+str(iqpt)+'.xml save/'+prefix+'.dyn_q'+label+'.xml')
      if (iqpt == 1):
        os.system('cp _ph0/'+prefix+'.dvscf* save/'+prefix+'.dvscf_q'+label)
        os.system('cp -r _ph0/'+prefix+'.phsave save/')
        #if the q2r run was made with non xml files copy the non xml ifc file
        if not os.path.isfile(prefix + '.fc.xml'):
            os.system('cp '+prefix+'.fc save/ifc.q2r')
        else:
            os.system('cp '+prefix+'.fc.xml save/ifc.q2r.xml')
      else:
        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf* save/'+prefix+'.dvscf_q'+label)
        os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )
    # Case without SOC
    if SOC == 'false':
      os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q'+label)
      if (iqpt == 1):
        os.system('cp _ph0/'+prefix+'.dvscf save/'+prefix+'.dvscf_q'+label)
        os.system('cp -r _ph0/'+prefix+'.phsave save/')
        os.system('cp '+prefix+'.fc save/ifc.q2r')
      else:
        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf save/'+prefix+'.dvscf_q'+label)
        os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )
 
  else:
    # Case with SOC
    if SOC == 'true':
      os.system('cp '+prefix+'.dyn0 '+prefix+'.dyn0.xml')
      os.system('cp '+prefix+'.dyn'+str(iqpt)+'.xml save/'+prefix+'.dyn_q'+label+'.xml')
      if (iqpt == 1):
        os.system('cp _ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
        os.system('cp -r _ph0/'+prefix+'.phsave save/')
        #if the q2r run was made with non xml files copy the non xml ifc file 
        if not os.path.isfile(prefix + '.fc.xml'):
            os.system('cp '+prefix+'.fc save/ifc.q2r')
        else:
            os.system('cp '+prefix+'.fc.xml save/ifc.q2r.xml')
      else:
        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
        os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )
    # Case without SOC
    if SOC == 'false':
      os.system('cp '+prefix+'.dyn'+str(iqpt)+' save/'+prefix+'.dyn_q'+label)
      if (iqpt == 1):
        os.system('cp _ph0/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
        os.system('cp -r _ph0/'+prefix+'.phsave save/')
        os.system('cp '+prefix+'.fc save/ifc.q2r')
      else:
        os.system('cp _ph0/'+prefix+'.q_'+str(iqpt)+'/'+prefix+'.dvscf1 save/'+prefix+'.dvscf_q'+label)
        os.system('rm _ph0/'+prefix+'.q_'+str(iqpt)+'/*wfc*' )
'''.format(pf = pf)]


#script to automatically determine the disentanglement windows.
wannier_windows_py = [''' 
from argparse import ArgumentParser
import sys

#read command line
parser = ArgumentParser()
parser.add_argument('num_of_wan', type=int)
parser.add_argument('bands', type=int)
parser.add_argument('points', type=int)
parser.add_argument('Ef', type=float)
parser.add_argument('bands_file', type=str)
args = parser.parse_args()

num_of_wan = args.num_of_wan
bands = args.bands
points = args.points
Ef = args.Ef
bands_file = args.bands_file

#read the band energies
e_vec = []

ifs = open(bands_file, "r")
for e in ifs:
    e_vec.append(float(e))

#get the index (wan_min) of the first band that gets wannierized i.e. crosses the fermi energy (indices starting at 0)
wan_min = 0
for i in range(0, bands):
    for j in range(0, points):
        index = i*points + j
        if e_vec[index] >= Ef:
            wan_min = i
            break
    if not wan_min == 0:
        break

if wan_min + num_of_wan > bands :
    sys.exit("You chose too many bands to wannierize")
    
#get the outer window by finding the smallest and largest values of the specified bands
smallest_outer = e_vec[wan_min*points]
largest_outer = smallest_outer
for i in range(wan_min*points, (wan_min+num_of_wan)*points):
    if e_vec[i] < smallest_outer:
        smallest_outer = e_vec[i]
    elif e_vec[i] > largest_outer:
        largest_outer = e_vec[i]

#find the "disturbing" bands that fall into the outer window and do not belong to the specified bands
dist_bands = []
for i in range(0, bands):
    if i > wan_min and i <= wan_min + num_of_wan:
        continue
    else:
        for j in range(0, points):
            index = i*points + j
            if e_vec[index] > smallest_outer and e_vec[index] < largest_outer:
                dist_bands.append(i)
                break

#find the inner window using the disturbing bands
smallest_inner = smallest_outer;
largest_inner = largest_outer;
for i in range(0, len(dist_bands)):
    smallest_local = e_vec[dist_bands[i]*points]
    largest_local = smallest_local
    for j in range(0, points):
        index = dist_bands[i]*points + j
        #get the largest and smallest values of the current band
        if e_vec[index] < smallest_local:
            smallest_local = e_vec[index]
        elif e_vec[index] > largest_local:
            largest_local = e_vec[index]

    #adjust the inner band using these maximum points of the disturbing band
    if largest_local > smallest_inner and largest_local < largest_inner:
        smallest_inner = largest_local
    if smallest_local < largest_inner and  smallest_local > smallest_inner:
        largest_inner = smallest_local

#increment band index by 1 as they start at 1 in qe
wan_min += 1
print("%12.8f    %12.8f          %12.8f          %12.8f          %d" %(smallest_inner, largest_inner, smallest_outer, largest_outer, wan_min))

''']

#_________________________________________________________________________________________#
#_________________________________________________________________________________________#

with open(os.path.join(base_dir, 'ELB/scf.in'), 'w') as txt:
    txt.writelines(scf_in)
with open(os.path.join(base_dir, 'ELB/bands.in'), 'w') as txt:
    txt.writelines(bands_in)
with open(os.path.join(base_dir, 'ELB/bands_ip.in'), 'w') as txt:
    txt.writelines(bands_ip_in)
with open(os.path.join(base_dir, 'ELB/wannier_windows.py'), 'w') as txt:
    txt.writelines(wannier_windows_py)
with open(os.path.join(base_dir, 'PHB/ph.in'), 'w') as txt:
    txt.writelines(ph_in)
with open(os.path.join(base_dir, 'PHB/q2r.in'), 'w') as txt:
    txt.writelines(q2r_in)
with open(os.path.join(base_dir, 'PHB/matdyn.in'), 'w') as txt:
    txt.writelines(matdyn_in)
with open(os.path.join(base_dir, 'PHB/pp.py'), 'w') as txt:
    txt.writelines(pp_py)
with open(os.path.join(base_dir, 'EPM/nscf.in_orig'), 'w') as txt:
    txt.writelines(nscf_in)
with open(os.path.join(base_dir, 'EPM/wannier.in_orig'), 'w') as txt:
    txt.writelines(wannier_in)
with open(os.path.join(base_dir, 'EPM/wannier_epm.in_orig'), 'w') as txt:
    txt.writelines(wannier_epm_in)
with open(os.path.join(base_dir, 'EPM/ph_lw.in_orig'), 'w') as txt:
    txt.writelines(ph_lw_in)
with open(os.path.join(base_dir, 'EPM/a2F.in_orig'), 'w') as txt:
    txt.writelines(a2F_in)
with open(os.path.join(base_dir, 'EPM/kmesh.py'), 'w') as txt:
    txt.writelines(kmesh_py)
with open(os.path.join(base_dir, 'ISO/eliashberg_iso.in_orig'), 'w') as txt:
    txt.writelines(eliashberg_iso)
with open(os.path.join(base_dir, 'ep_bands.sh'), 'w') as txt:
    txt.writelines(ep_bands_sh)
with open(os.path.join(base_dir, 'epw.sh'), 'w') as txt:
    txt.writelines(epw_sh)

print("done")
