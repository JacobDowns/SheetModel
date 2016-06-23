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
You can test everything is working by opening an ``ipython`` terminal and typing:

```python
from sheet_model import *
```


