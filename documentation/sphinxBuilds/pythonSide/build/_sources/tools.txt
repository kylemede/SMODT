Tools
*****

"Tools" are the SMODT nomenclature for utilities or library functions and they 
are contained in "toolboxes".  These were created as the simulator was built 
to not only provide code clarity through abstraction, but also to minimize 
code redundancy.  These toolboxes have grown to contain many useful functions 
for not just performing the standard simulation runs and post-completion plotting, 
but others that can be used directly by the user's own scripts to create more 
specific plots or perform other statistical calculations on the results.

The C++ side of the simulator performs the computational heavy tasks
involved in making a MCMC or Simulated Annealing chains and some of the 
post-completion statistical calculations.  As such, many of the tools on the 
Python side have been translated into twins on the C++ side.  Testing was done
to ensure the both versions of those tools provide the exact same results to 
within the required precision.
   
General Tools
=============================

.. automodule:: toolboxes.generalToolbox
   :members:

DI Tools
=============================

.. automodule:: toolboxes.DItoolbox
   :members:
   
RV Tools
=============================

.. automodule:: toolboxes.RVtoolbox
   :members:
   
Plotting Tools
===============================

.. automodule:: toolboxes.plotToolbox
   :members:

RAM tracking Tool
===============================

.. automodule:: toolboxes.ramTracker
   :members:

Artificial Data Making Tool
===============================

.. automodule:: toolboxes.artificialDataMaker
   :members:

