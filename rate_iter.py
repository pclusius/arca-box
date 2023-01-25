import numpy as np
import sys,os,time
from subprocess import Popen, PIPE, STDOUT
import shutil
from pdb import set_trace as bp

nr_of_iterations = 500

arcaroot_dir = './'
tempfile = 'test'

chem = os.path.join(arcaroot_dir,'src/chemistry/','CH4')
outpath = os.path.join(arcaroot_dir,'INOUT/MCMC_0000/BASENOX/')
# if intermediate output is saved
save_all = True

# should not change but why not
exe_name = './arcabox.exe'
# to measure time
start_time = time.time()
# change to where the program is called
os.chdir(arcaroot_dir)

# read the reactions that will be changed
order, indices = np.genfromtxt(os.path.join(chem,'reactions.txt'), usecols=(0,1), skip_header=1, unpack=True, dtype=int)
indices = indices[order>0]
order   = order[order>0]

# Backup original RATES.dat if not already existing
if not os.path.exists(os.path.join(chem,'ORIG_RATES.dat')) and os.path.exists(os.path.join(chem,'RATES.dat')):
    shutil.copyfile(os.path.join(chem,'RATES.dat'), os.path.join(chem,'ORIG_RATES.dat'))

# Start the iteration loop
iteration = 0

def run_shit(ones=False):
    string_to_print = '&NML_RATES\n'
    for i,ii in enumerate(indices):
        if ones:
            string_to_print += '  RCONST(%d) = %f \n' %(ii, 1.0)
        else:
            # random reaction rate constant multipliers - here should be the MCMC guesses or something
            string_to_print += '  RCONST(%d) = %f \n' %(ii, 10**(np.random.random()*8 -4.0))

    # Footer for the F namelist
    string_to_print += '/\n'
    # write F namelist
    with open(os.path.join(chem,'RATES.dat'), 'w+') as file:
        file.write(string_to_print)

    # Call ARCA box
    box = Popen([exe_name, "%s"%tempfile], stdout=PIPE,stderr=STDOUT,stdin=None)

    # poll the status
    status = box.poll()

    while status == None:
        # poll the status as long as it runs
        status = box.poll()

    # put a bullet in the neck to be sure
    box.kill()
    # check no pulse
    status = box.poll()
    if ones:
        shutil.move(os.path.join(outpath,'CHEM_FINAL.txt'),
            os.path.join(outpath,'CHEM_FINAL_BASE.txt'))

    return

# create base case, the TRUTH
run_shit(True)


while iteration < nr_of_iterations:
    # header for the produced namelist

    print('iteration', iteration)
    run_shit()
    # If each iteration final concetration and reaction rate coefficients are saved
    if save_all:
        shutil.move(os.path.join(outpath,'CHEM_FINAL.txt'),
            os.path.join(outpath,'CHEM_FINAL_%03d.txt'%iteration))

        shutil.copyfile(os.path.join(chem,'RATES.dat'),os.path.join(chem,'RATES%03d.dat'%iteration))

    iteration += 1

print("--- runtime %s seconds / iteration---" % ((time.time() - start_time)/iteration))
