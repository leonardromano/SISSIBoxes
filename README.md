# SISSIBoxes

This script is designed to make the user familiar with how to obtain uniform fluid boxes from the SISSI simulation snapshot data on OP3PGN.   
You can use the functions introduced in this script to build a more automated analysis (i.e. loops over snapshots, bubbles, etc.).   
Paths to data files are hard-coded into the source code. This has the advantage that it allows the user to run these scripts with minimal knowledge about the data management,   
but it should be kept in mind, if one tries to run it on a different system. That simply won't work unless the code is modifed significantly.   
Currently only Leonard Romano, Andreas Burkert, Manuel Behrendt and Mingyu Hu have access to these files on OP3PGN (and system admins).  

**Parallelization Note**: The core functions allow for threaded use.   
Execution can be heavily accelerated by using multiple threads, i.e. by calling the script with multiple threads      
```julia -t nthreads get_SISSI_boxes.jl```, where _nthreads_ is the number of parallel threads. 
