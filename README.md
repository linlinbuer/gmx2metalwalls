# gmx2metalwalls
Generate metalwalls input files (data.inpt and runtime.inpt) based on .gro file (without electrodes) and .itp file of solute.

# Usage
- Required files: .gro (system without electrodes), .itp (topology file of molecular solute)
- copy input.dat from database folder to your working folder
- modify input.dat based on your system
- python3 gmx2metalwalls -c *.gro  -t *.itp -f input.dat
