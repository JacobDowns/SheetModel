import abc
from dolfin import *
from constants import *
from colored import fg, attr

""" Abstract model class. """

class Model(object):
  __metaclass__ = abc.ABCMeta

  def __init__(self, model_inputs):

    # Load a dictionary of model inputs
    self.model_inputs = model_inputs

    # Check if an input file is specified 
    if 'input_file' in self.model_inputs:
      self.input_file = HDF5File(mpi_comm_world(), self.model_inputs['input_file'], 'r') 
    else :
      # Otherwise yell at the user (that would be myself)
      raise Exception("Please specify a valid input file.")
      sys.exit(1)
      
    # Output directory for writing model output
    if 'out_dir' in self.model_inputs:
      self.out_dir = model_inputs['out_dir'] + '/'
    else :
      # Otherwise yell at the user (that would be myself)
      raise Exception("Please specify an output directory.")
      sys.exit(1)
    
    # Load the mesh        
    self.mesh = Mesh()
    self.input_file.read(self.mesh, "mesh", False)  
      
    # Model time
    self.t = 0.0
      

  # Loads model inputs and initializes input and output files
  def init_model(self):
    # If there is a dictionary of physical constants specified in model_inputs, 
    # use it. Otherwise use the defaults. 
    if 'constants' in self.model_inputs :
      self.pcs = self.model_inputs['constants']
    else :
      self.pcs = pcs
    
    if self.to_continue():
      # Continue an existing simulation
    
      # Close the input file and re-open in append mode
      self.input_file.close()
      self.input_file = HDF5File(mpi_comm_world(), self.model_inputs['input_file'], 'a') 
      
      # Checkpoint files will be appended to the existing output file rather than 
      # making a new one 
      self.output_file = self.input_file
      
      self.continue_simulation()
    else :
      # Start a new simulation 
      
      # Check if there's a given name for the checkpoint file
      checkpoint_file = 'out'
      if 'checkpoint_file' in self.model_inputs:
        checkpoint_file = self.model_inputs['checkpoint_file']
      
      # Create a new ouput file for checkpointing.
      self.output_file = HDF5File(mpi_comm_world(), self.out_dir + checkpoint_file + ".hdf5", 'w')
      
      self.start_simulation()
      

  @abc.abstractmethod
  def to_continue(self):
      """Check the input file to check if we're starting or continuing a simulation."""
      return


  @abc.abstractmethod
  def load_inputs(self):
      """Load model inputs from the input file."""
      return
      
      
  @abc.abstractmethod
  def continue_simulation(self):
      """Continue an existing simulation."""
      return
      
      
  @abc.abstractmethod
  def start_simulation(self):
      """Start a new simulation."""
      return
      
  
  @abc.abstractmethod
  def step(self, dt):
      """Step the model forward by dt."""
      return
        
  # Assigns a value to a model variable. This function first looks in the 
  # model_inputs dictionary then the input_file. 
  def assign_var(self, name, var, var_func):
    # Look for the variable in the model_inputs dictionary. It could be either
    # a function or expression.
    if 'name' in self.model_inputs:
      input_var = self.model_inputs[name]
      
      if isinstance(input_var, dolfin.Expression):
        # If the input variable is an expression, we need to project it to 
        # get a function
        var = input_var
        var_func.assign(project(input_var, self.V_cg))
      else :
        # If the input variable is a function, then we'll copy it rather than
        # use the original function to prevent unwanted metaphysical linkage
        var.assign(input_var)
        var_func.assign(input_var)
          
    else :
      # If we didn't find anything in the model_inputs dictionary, check the input_file.
      if self.read_var(name, var) :
        var_func.assign(var)
      else :
        # If the input isn't there either raise an exception
        raise Exception("Could not load model input: " + str(name))


  # Assigns a value to a function. Looks in the model inputs dictionary then
  # the inputs file        
  def assign_func(self, f, name):
    if name in self.model_inputs:
      f.assign(self.model_inputs[name])
    else :
      self.input_file.read(f, name)
    
  
  # Reads in a variable from the input file
  def read_var(self, name, var):
    try :
      self.input_file.read(var, name)
      return True
    except :
      return False