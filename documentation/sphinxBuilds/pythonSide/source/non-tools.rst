Simulation Starter Functions
****************************

SMODT simulations are initialized with a Python 'starter' function which 
checks the settings in the settings file are valid, there is data in the data
files, creates an output directory and puts a copy of the SMODT code into it.
The MCMC Process Manager is then called to start the individual Simulated 
Annealing and MCMC chains.  As this stage of the simulation is very 
computational heavy, it is taken care of with C++ code to take advantage of 
its speed.  Once the chains are complete, other C++ functions will conduct 
the heavy stages of the post-completion statistical calculations.  Following 
those, the MCMC Process Manager will then take care of a few other statistical
calculations, produce plots summarizing the results and write the final output
data to disk.  The user can then simply review the results and plots in the 
output directory, and even write scripts of their own to make further plots or 
statistical calculations on the output data.  For such tasks, please review the 
tools already available to make writting these personalized scripts 
faster/easier.

SMODT Starter
=============================

.. automodule:: SMODTstarter
   :members: main

Simulation Manager
====================================

.. automodule:: SimulationManager
   :members: 

   

   
