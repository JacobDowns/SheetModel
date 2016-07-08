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

To install, add a line to you .bashrc file:
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
in the gen\_inputs directory. These will create some default hdf5 inputs files in the inputs_sheet/inputs directoryr called


