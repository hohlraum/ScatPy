# ScatPy Configuration file

# Paths should be specififed using unix forward slashes / even on 
# Windows systems

# The profile name
name = 'system_default'

# The absolute path to the folder containing materials properties
mat_library = '~/Documents/Analysis/ddscat/mp/'

# The absolute path to the folder containing the DDSCAT executable
ddscat_path = '/usr/local/bin/'

# The absolute path to the folder containing the mpi executable
#mpi_path = '/usr/mpi/gcc/openmpi-1.4.3/bin/'

# The operating system ('unix', 'mac', 'windows')
os = 'unix'

# An optional script that is run when the job is written

#def write_script(job):
#   """Useful postflight processing"""
#    pass