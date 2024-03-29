# Explanation

The ***path_enumeration.py*** file contains all functions needed for full path enumeration and enumeration with multi-bin constraint.\
The ***pairedbin_enumeration.py*** file contains all functions needed for converting paired-bins into multi-bins and for post-filtering of transcripts.\
The ***main.py*** file is used to read the graph data, create the corresponding graphs, extract the paths and solve the optimization problem. \
The ***optimize.py*** file is used to solve the optimization problem with gurobi. \
The ***flowProblem.py*** file is used to solve the problem using a flow network.


# Usage

***$ python main.py [input.graph] [options]***  \
***python main.py -help*** lists all available options

# Example
Full path enumeration and optimization with gurobi using L1 norm \
***$ python main.py human_geuvadis_simulated_5sets.graph -full -opt -norm1***
