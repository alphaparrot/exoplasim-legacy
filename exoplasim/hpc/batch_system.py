#This is just an importation wrapper--various modules import this, and this imports 
#various predefined quantities from different possible submission systems.


#Select a batch submission system (ONLY UNCOMMENT ONE)
from torque import *
#from slurm import *