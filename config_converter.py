import MDAnalysis as mda
import numpy as np
from commands import parse_args
from input_reader import read_input_file

A_TO_BOHR = 1.88972613  # Conversion factor from angstrom to Bohr
args = parse_args()  # Get command-line arguments
# read input parameters

if args.input:
    params = read_input_file(args.input)
    solute_name = params["solute_type"]
    ion_names = params["ion_types"]


def load_system(gro_file, xyz_file):
    """Load GROMACS (.gro) and gold slab (.xyz) files into MDAnalysis Universes."""
    u = mda.Universe(gro_file)
    u_gold = mda.Universe(xyz_file)
    return u, u_gold

def convert_to_bohr(atoms):
    """Convert atom coordinates to Bohr."""
    atoms.positions *= A_TO_BOHR

def translate_system(water, solute, ions, au1, au2):
    """Translate the system to align with the gold slab."""
    z_sys_high = max(np.max(water.positions[:,2]), np.max(ions.positions[:,2]))
    z_sys_low = min(np.min(water.positions[:,2]), np.min(ions.positions[:,2]))
    z_au1, z_au2 = np.max(au1.positions[:,2]), np.min(au2.positions[:,2])
    dz_gold = z_au2 - z_au1
    dz_system = z_sys_high - z_sys_low
    
    if dz_system > dz_gold:
        print(dz_system, dz_gold)
        raise ValueError("Water and solute do not fit within the gold slab.")
    
    z_trans = z_sys_high - z_au2 + 5
    print(f"Translating system by {z_trans} Bohr")
    water.translate([0, 0, -z_trans])
    solute.translate([0, 0, -z_trans])
    ions.translate([0, 0, -z_trans])

def sort_atoms_by_name(atoms, reference_resid):
    """Sort atoms in the order of the topology file."""
    order = {name: index for index, name in enumerate(reference_resid.names)}
    sorted_indices = np.argsort([order[name] for name in atoms.names], kind='stable')
    return atoms[sorted_indices]


def process_conversion(gro_file, xyz_file, output_file):
    """Main function to process conversion from GROMACS to MetalWalls format."""
    u, u_gold = load_system(gro_file, xyz_file)
    water = u.select_atoms('resname SOL')
    if solute_name:
        solute = u.select_atoms(f'resname {solute_name}')
    else:
        solute = u.select_atoms("resid 999999")
    
    ion_atomgroups = {}  # Dict to hold each ion's AtomGroup
    ions_combined = u.select_atoms('resname DUMMY_NOTHING')  # start with empty AtomGroup

    for name in ion_names:
        atomgroup = u.select_atoms(f'resname {name}')
        # print(f'done with {name}, count: {len(atomgroup)}')
        ion_atomgroups[name] = atomgroup
        ions_combined += atomgroup

    au1 = u_gold.select_atoms("name Au1")
    au2 = u_gold.select_atoms("name Au2")

    
    
    
    convert_to_bohr(water)

    convert_to_bohr(solute)
    convert_to_bohr(ions_combined)

    translate_system(water, solute, ions_combined, au1, au2)
    sorted_water = sort_atoms_by_name(water, u.select_atoms(f'resid {min(water.resids)}'))
    if solute_name:
        sorted_solute = sort_atoms_by_name(solute, u.select_atoms(f'resid {min(solute.resids)}'))
    
    if solute_name:
        combined_atoms = sorted_water + sorted_solute + ions_combined
    else:
        combined_atoms = sorted_water + ions_combined
        # Rename specific atom names
    for atom in combined_atoms:
        if atom.name == 'OW':
            atom.name = 'O'
        elif atom.name == 'HW1':
            atom.name = 'H1'
        elif atom.name == 'HW2':
            atom.name = 'H2'
    combined_gold_slab = au1 + au2
    
    # Write the MetalWALL input file
    output_file = 'data.inpt'  # Output file name
    with open(output_file, 'w') as f:
        # # Write the number of atoms
        # f.write(f"{len(combined_atoms) + len(combined_gold_slab)}\n")
        # f.write(f"daje\n")
        # Write the header
        Natoms = len(combined_atoms) + len(combined_gold_slab)
        f.write("# header\n")
        f.write("step                           1\n")
        f.write(f"num_atoms                      {Natoms}\n")
        f.write(f"num_electrode_atoms             {len(combined_gold_slab)}\n")
        f.write("# box\n")
        f.write(" 6.922066797550001E+01  6.922066797550001E+01  1.915000000000000E+02\n")
        f.write(f"# coordinates :    {Natoms} atoms - step  2.790584176380000E+08\n")
        
        
        # Write the atom names and coordinates
        for atom in combined_atoms:
            name = atom.name
            x, y, z = atom.position
            f.write(f"{name} {x:.6f} {y:.6f} {z:.6f}\n")
            
        # Write the atom names and coordinates
        for atom in combined_gold_slab:
            name = atom.name
            x, y, z = atom.position
            f.write(f"{name} {x:.6f} {y:.6f} {z:.6f}\n")

        print(f"{output_file} has been created.")



def count(gro_file, xyz_file):
    """Main function to process conversion from GROMACS to MetalWalls format."""
    u, u_gold = load_system(gro_file, xyz_file)
    water = u.select_atoms('resname SOL')
    solute = u.select_atoms(f'resname {solute_name}')
    gold = u_gold.select_atoms('name Au1 Au2')

    ion_atomgroups = {}  # Dict to hold each ion's AtomGroup
    ions_combined = u.select_atoms('resname DUMMY_NOTHING')  # start with empty AtomGroup

    for name in ion_names:
        atomgroup = u.select_atoms(f'resname {name}')
        ion_atomgroups[name] = atomgroup
        ions_combined += atomgroup

    
    n_water = water.residues.n_residues
    n_solute = solute.residues.n_residues
    n_ions = {}
    n_gold = gold.n_atoms

    for name, ag in ion_atomgroups.items():
        n_ions[name] = ag.residues.n_residues
        print(f"{name}: {n_ions[name]} ions")
    print(f"number of electrodes: {n_gold}")
    
    return n_water, n_solute, n_ions, n_gold