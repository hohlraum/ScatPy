# ScatPy Configuration file
# For running DDSCAT serial on landau

# The profile name
name = 'landau'

# The absolute path to the folder containing materials properties
mat_library = '/home/guests/mark/mat_prop'

# The absolute path to the folder containing the DDSCAT executable
ddscat_path = '/opt/local/lib/ddscat/src'

# The absolute path to the folder containing the mpi executable
mpi_path = '/usr/mpi/gcc/openmpi-1.4.2/bin/'

# The number of slots to use 
num_slots = (16,32)

# The operating system
os = 'unix'


def write_script(self, write_sge=False):
    import posixpath
    
    with open(os.path.join(self.folder, 'submit.sge'), 'wb') as f:
 
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
        mpi=posixpath.join(path, 'mpirun')
        f.write('time %s -np $NSLOTS -machinefile $TMPDIR/machines /cluster/bin/ddscat_openmpi\n' % (mpi))
        f.write('echo completed `pwd`\n')
        f.write('echo \'------------------------------------------\'\n')
        f.write('date\n')
         
        f.write('foreach old (ddscat_mpi_PRE_.*)\nmv $old "$old:gas/_PRE_//"".txt"\nend\n')
