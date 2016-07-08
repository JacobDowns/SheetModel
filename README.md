This repository is an implementation of a model of subglacial drainage through linked cavities. To run the code, you need FEniCS 1.6, which can be installed on Ubuntu by running the following commands.

```bash
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install fenics
sudo apt-get dist-upgrade
```

You will also need to install the Python colored library.

```bash
sudo pip install colored
```

To install, add a line to your .bashrc file:
```bash
export PYTHONPATH="<insert SheetModel directory here>:$PYTHONPATH"
```
You can test the installation by opening an ``ipython`` terminal and typing:

```python
from sheet_model import *
```
The current repository contains only the model source code. A plethora of model experiments can be found here: https://github.com/JacobDowns/SheetExperiments

Once you have the experiments, the easist way to get started is to run the two commands

```bash
python gen_inputs_high.py
python gen_inputs_low.py
```
in the SheetExperiments/gen\_inputs directory. These scripts will create some default hdf5 inputs files in the inputs\_sheet/inputs directory called inputs\_high.hdf5 and inputs\_low.hdf5. These input files contain fields including bed elevation, ice thickness, initial sliding speed, initial melt rate, and boundary conditions for use in synthetic simulations.  

Frome here you will be able run the sim_ref_steady.py script in SheetExperiments/ref\_steady. 

```python
"""
Reference simulation steady state.
"""

from dolfin import *
from dolfin import MPI, mpi_comm_world
from sheet_model import *
from constants import *


### Model inputs

# Process number
MPI_rank = MPI.rank(mpi_comm_world())

model_inputs = {}
pcs['k'] = 5e-3
pcs['l_r'] = 0.1
model_inputs['input_file'] = '../inputs_sheet/inputs/inputs_high.hdf5'
model_inputs['out_dir'] = 'out_ref_steady/'
model_inputs['constants'] = pcs
model_inputs['newton_params'] = prm

# Create the sheet model
model = SheetModel(model_inputs)


### Run the simulation

# Seconds per day
spd = pcs['spd']
# End time
T = 90.0 * spd
# Time step
dt = spd / 2.0

while model.t < T:  
  if MPI_rank == 0: 
    current_time = model.t / spd
    print ('%sCurrent time: %s %s' % (fg(1), current_time, attr(0)))
  
  model.step(dt)
  model.write_pvds(['pfo', 'h'])
  
model.write_steady_file('../inputs_sheet/steady/ref_steady')
```


