# topology_reconstruction
Python codes for the topology reconstruction after polymer condensation. The codes process files output from LAMMPS only. Namely, a bond list bond.* is read as input. Example input LAMMPS script (run.lmp) is also provided in this repository. 

==================== input and output files =======================================

input files: bond.* lammps bond lists.

output files: hist_linear_*.txt, hist_ring_*.txt: histograms of the polymer lengths in linear and ring state vs time. 

average_length.txt: the average polymer length in time (average over linear and rings)

Stat_polys.txt: first column = time, second column = the number of polymers in linear state, third column = the number of polymers in ring state

===================== technical details ============================================

To run: with python 3 interpreter run the main.py file. The input parameters are Nmols = initial number of molecules in the system before the condensation process begins, Nmon = initial number of beads per polymer, t_0 = the first time frame to read (i.e., bonds.0 if t_0 = 0), t_max = the last time frame to read (i.e., bonds.1000 if t_max = 1000), dt = the interval between two time frames (how often the bond lists are saved from LAMMPS), sim_dt = the timestep used during the simulation. 

==================== information on the codes ======================================

main.py: The main script that reads the input parameters and calls the readoutput_lmp.topo_reconstruct.

readoutput_lmp.py: Topology reconstruction code. The bond list are read here and the lookup subrutine is called to find the connectivity matrix between all the connected beads. 

properties.py: some properties, such as the polydispersity index, are calculated here. More properties can be added in this script later on by the user. 
