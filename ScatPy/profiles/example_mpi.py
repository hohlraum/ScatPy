# ScatPy Configuration file
# For running DDSCAT on landau parallel

# The profile name
name = 'landau_mpi'

# The operating system
os = 'unix'

# The absolute path to the folder containing materials properties
mat_library = '/home/guests/mark/mat_prop'

def write_script(job):
    """A script to generate a SGE submission script."""
    import posixpath
    import os.path

    # The number of slots to use 
    num_slots = (16,32)
    
    with open(os.path.join(job.folder, 'submit.sge'), 'wb') as f:
 
        f.write('#!/bin/csh\n' )
        f.write('#\n#\n#\n')
        f.write('# ---------------------------\n')
        f.write('# our name \n')
 
        f.write('#$ -N ddscat_mpi_PRE_\n#\n')
        f.write('# pe request\n')                    
        f.write('#$ -pe openmpi %d-%d\n' % tuple(num_slots))
 
        f.write('#\n')
 
        f.write('# stderr >& stdout\n')
        f.write('#$ -j y\n')
        f.write('#\n')
        f.write('# ---------------------------\n')
                  
        f.write('set hostname=`/bin/hostname`\n')
 
        f.write('echo beginning `pwd`\n')
        f.write('date\n')
        mpi = '/usr/mpi/gcc/openmpi-1.4.2/bin/mpirun')
        f.write('time %s -np $NSLOTS -machinefile $TMPDIR/machines /cluster/bin/ddscat_openmpi\n' % (mpi))
        f.write('echo completed `pwd`\n')
        f.write('echo \'------------------------------------------\'\n')
        f.write('date\n')
         
        f.write('foreach old (ddscat_mpi_PRE_.*)\nmv $old "$old:gas/_PRE_//"".txt"\nend\n')
