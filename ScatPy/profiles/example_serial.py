# ScatPy Configuration file
# For running luna in serial on SGE via qsub

# The profile name
name = 'luna_serial'

# The operating system
os = 'unix'

# The absolute path to the folder containing materials properties
mat_library = '/home/guests/mark/mat_prop'

def write_script(job):
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
        f.write('time /cluster/bin/ddscat\n' % ddscat_path)
        f.write('echo completed `pwd`\n')
        f.write('echo \'------------------------------------------\'\n')
        f.write('date\n')
         
        f.write('foreach old (ddscat_ser_PRE_.*)\nmv $old "$old:gas/_PRE_//"".txt"\nend\n')
