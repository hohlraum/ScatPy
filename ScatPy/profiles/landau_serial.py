# ScatPy Configuration file
# For running DDSCAT serial on landau

# The profile name
name = 'landau'

# The operating system
os = 'unix'

# The absolute path to the folder containing materials properties
mat_library = '/home/guests/mark/mat_prop'

# The absolute path to the folder containing the DDSCAT executable
ddscat_path = '/opt/local/lib/ddscat/src'

# The absolute path to the folder containing the mpi executable
mpi_path = '/usr/mpi/gcc/openmpi-1.4.2/bin/'


def write_script(job):
    """
    Write a script file for submitting the job to SGE via qsub
    """
    import os.path
    with open(os.path.join(job.folder, 'submit.sge'), 'wb') as f:

        f.write('#!/bin/csh\n' )
        f.write('#\n#\n#\n')
        f.write('# ---------------------------\n')
        f.write('# our name \n')

        f.write('#$ -N ddscat_ser_PRE_\n')
        f.write('#\n')
            
        f.write('# stderr >& stdout\n')
        f.write('#$ -j y\n')
        f.write('#\n')
        f.write('# ---------------------------\n')
        
        
        f.write('set hostname=`/bin/hostname`\n')

        f.write('echo beginning `pwd`\n')
        f.write('date\n')
        f.write('time /cluster/bin/ddscat\n')
        f.write('echo completed `pwd`\n')
        f.write('echo \'------------------------------------------\'\n')
        f.write('date\n')
        
        f.write('foreach old (ddscat_ser_PRE_.*)\nmv $old "$old:gas/_PRE_//"".txt"\nend\n')
