from pathlib import Path
import numpy as np
from commands import parse_args 
from input_reader import read_input_file
from config_converter import count

# Constants for LJ pairs
lj_CL = "lj_pair CL CL  0.05349 4.830"
lj_NA = "lj_pair NA NA 1.475 2.160"
lj_CS = "lj_pair CS CS 0.376 3.601"


# Path to the script's directory
script_dir = Path(__file__).resolve().parent
data_dir = script_dir / 'database'


# Helper functions
def format(r):
    return "%.7f" % r

def ff(r):
    return float(format(r))

def str2Value(line):
    items = line.strip('\t').split()
    converted_items = []
    for item in items:
        try:
            converted_item = int(item)
        except ValueError:
            try:
                converted_item = float(item)
            except ValueError:
                converted_item = str(item)
        converted_items.append(converted_item)
    return converted_items

def extract_section(section_name):
    """Extracts all lines under a given section until the next section starts."""
    extracted_data = []
    inside_section = False

    with open(args.topology, "r") as file:
        for line in file:
            line = line.strip()

            # Detect a new section (stops if a new [] block appears)
            if line.startswith("[") and line.endswith("]"):
                if inside_section:
                    break  # Stop reading when a new section starts
                inside_section = line.strip("[]").strip() == section_name
                continue  # Skip section header itself
            
            # Collect data if inside the target section
            if inside_section and line and not line.startswith(";"):  # Skip empty lines and comments
                extracted_data.append(str2Value(line))

    return extracted_data

def check_if_impropers(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    return any("impropers" in line.lower() for line in lines)

def convert_bonds(bonds, atomType):
    N = len(bonds)
    bonds_converted = []
    constraints = []
    A1 = []
    A2 = []

    # Extract indices and parameters
    index1 = [bonds[i][0:2] for i in range(N)]
    r0_gro = np.array([bonds[i][3] for i in range(N)])
    k0_gro = np.array([bonds[i][4] for i in range(N)])

    # Map atom indices to atom types
    for i in range(len(index1)):
        A1.append(atomType[index1[i][0] -1])# -1 because python starts from 0, but the index of atom starts from 1 in gromacs
        A2.append(atomType[index1[i][1] -1])

    # Convert units
    r0_MW = r0_gro * 18.8973  # nm to bohr
    k0_MW = (k0_gro * 0.00038 / (18.8973 * 18.8973)) / 2  # gro to MW units

    # Separate bonds and constraints
    for i in range(N):
        if A1[i].startswith("H") or A2[i].startswith("H"):
            # Add to constraints
            constraints.append(f"constraint {A1[i]} {A2[i]} {ff(r0_MW[i])}")
        else:
            # Add to bonds
            bonds_converted.append(f"harmonic_bond {A1[i]} {A2[i]} {ff(k0_MW[i])} {ff(r0_MW[i])}")


    return bonds_converted, constraints

def convert_angles(angles, atomType):
    N = len(angles)
    angles_converted = []
    A1 = [atomType[angles[i][0] - 1] for i in range(N)]
    A2 = [atomType[angles[i][1] - 1] for i in range(N)]
    A3 = [atomType[angles[i][2] - 1] for i in range(N)]
    theta0_MW = [ff(np.deg2rad(angles[i][4])) for i in range(N)]  # deg to rad
    k0_MW = [ff(angles[i][5] * 0.00038 / 2) for i in range(N)]  # gro to MW units

    for i in range(N):
        angles_converted.append(f"harmonic_angle {A2[i]} {A1[i]} {A3[i]} {k0_MW[i]} {theta0_MW[i]}")

    return angles_converted

def convert_dihedrals(dihedrals, atomType):
    N = len(dihedrals)
    di_converted = []
    A1 = [atomType[dihedrals[i][0] - 1] for i in range(N)]
    A2 = [atomType[dihedrals[i][1] - 1] for i in range(N)]
    A3 = [atomType[dihedrals[i][2] - 1] for i in range(N)]
    A4 = [atomType[dihedrals[i][3] - 1] for i in range(N)]
    v1_MW = [ff(dihedrals[i][5] * 0.00038) for i in range(N)]
    v2_MW = [ff(dihedrals[i][6] * 0.00038) for i in range(N)]
    v3_MW = [ff(dihedrals[i][7] * 0.00038) for i in range(N)]
    v4_MW = [ff(dihedrals[i][8] * 0.00038) for i in range(N)]

    for i in range(N):
        di_converted.append(f"dihedral {A1[i]} {A2[i]} {A3[i]} {A4[i]} {v1_MW[i]} {v2_MW[i]} {v3_MW[i]} {v4_MW[i]}")

    return di_converted

def convert_improper(impropers, atomType):
    N = len(impropers) - 1
    di_converted = []
    A1 = [atomType[impropers[i][0] - 1] for i in range(N)]
    A2 = [atomType[impropers[i][1] - 1] for i in range(N)]
    A3 = [atomType[impropers[i][2] - 1] for i in range(N)]
    A4 = [atomType[impropers[i][3] - 1] for i in range(N)]
    v1_MW = [ff(impropers[i][6] * 0.00038 * 2) for i in range(N)]
    v2_MW = [ff(0) for i in range(N)]
    v3_MW = [ff(0) for i in range(N)]
    v4_MW = [ff(0) for i in range(N)]

    for i in range(N):
        di_converted.append(f"improper {A1[i]} {A2[i]} {A3[i]} {A4[i]} {v1_MW[i]} {v2_MW[i]} {v3_MW[i]} {v4_MW[i]}")

    return di_converted

def convert_LJ(LJ, atom_dict):
    N = len(LJ)
    LJ_converted =[]
    inx = []
    epsilon = []
    sigma = []
    for i in range(len(LJ)):
        inx.append(LJ[i][0])
        sigma.append(LJ[i][5]*10)
        epsilon.append(LJ[i][6])
    atom_names = [atom_dict[key] for key in inx]
    
    for i in range(N):
        LJ_converted.append(f"lj_pair {atom_names[i]} {atom_names[i]} {epsilon[i]} {sigma[i]}")

    return LJ_converted
    
    
# Parse command-line arguments
args = parse_args()
if args.input:
    params = read_input_file(args.input)
    charge_rescale = params["charge_rescale"]
    solute_name = params["solute_type"]
    ion_names = params["ion_types"]
    
if args.config:
    n_water, n_solute, n_ions, n_gold = count(args.config, data_dir/'gold_au.xyz')


# Generate species section
def generate_species_section(atoms, charge_rescale, n_solute):
    N = len(atoms)
    atomType = [atoms[i][4] for i in range(N)]
    charge = [ff(atoms[i][6])*float(charge_rescale) for i in range(N)]
    mass = [atoms[i][7] for i in range(N)]

    species_section = []
    for i in range(N):
        species_section.append(f"""
    species_type
      name        {atomType[i]}
      count      {n_solute}
      charge  point {charge[i]}
      mobile True
      mass   {mass[i]}
""")
    return "\n".join(species_section)


# Main function
def main():
    print('****Grabbing info to write into runtime.inpt****')
    if solute_name:
        print(f'***Getting solute info for {solute_name}****')
        # Extract data from topology file
        atoms = extract_section('atoms')
        bonds = extract_section('bonds')
        angles = extract_section('angles')
        # angles = extractIntra('angles', 'dihedrals')
        impropers = extract_section('impropers')
        dihedrals = extract_section('dihedrals')
        LJ = extract_section('atomtypes')


        atomType = [atoms[i][4] for i in range(len(atoms))]
        atomInd = [atoms[i][1] for i in range(len(atoms))]
        atom_dict = {k: v for k, v in zip(atomInd, atomType)}

        # Convert sections
        bonds_converted, constraints = convert_bonds(bonds, atomType)
        angles_converted = convert_angles(angles, atomType)
        dihedrals_converted = convert_dihedrals(dihedrals, atomType)
        impropers_converted = convert_improper(impropers, atomType) if impropers else []
        LJ_converted = convert_LJ(LJ, atom_dict)

        # Generate species section
        species_solute_section = generate_species_section(atoms, charge_rescale, n_solute)
    
    ion_types = list(n_ions.keys())
    print(f'***Getting ions info for {ion_types}****')
    
    lj_ion_section = []
    species_ion_data = []
    
    if 'NA' in ion_types: 
        with open(data_dir/"species_Na_template.dat", "r") as f:
            species_ion = f.read()
        lj_ion_section.append(lj_NA)
        species_ion_data.append(species_ion.format(n_ion=n_ions['NA']))
        print(f"**writing NA in species section**")
    if 'CL' in ion_types:
        with open(data_dir/"species_Cl_template.dat", "r") as f:
            species_ion = f.read()
        lj_ion_section.append(lj_CL)
        species_ion_data.append(species_ion.format(n_ion=n_ions['CL']))
        print(f"**writing CL in species section**")

    species_ion_data_str = "\n\n".join(species_ion_data)
    lj_ion_section_str = "\n".join(lj_ion_section)

    with open(data_dir/"runtime_template.dat", "r") as f:
        RUNTIME_TEMPLATE = f.read()

    # Generate runtime content
    
    if solute_name:
        print('***write solutes and ions***')
        cons_para = f"\t constraints_algorithm rattle 1.0e-7 100"
        # Combine intra section
        intra_section = bonds_converted + constraints + [cons_para] + angles_converted + dihedrals_converted
    
            
        runtime_content = RUNTIME_TEMPLATE.format(
            n_water=n_water,  # Placeholder, update as needed
            species_solute_section=species_solute_section,
            species_ion_section=species_ion_data_str,  # Placeholder, update as needed
            n_Au1 = n_gold//2,
            n_Au2 = n_gold//2,
            molecule_solute_section=f"molecule_type\n  \t name {solute_name}\n  \t count {n_solute}\n \t  sites {' '.join(atomType)}",
            intra_section= "\n".join("\t" + line for line in intra_section),
            lj_solute_section="\n".join("\t" + line for line in LJ_converted),
            lj_ion_section= lj_ion_section_str
        )
    else:
        print('***No solute, write ions***')
        runtime_content = RUNTIME_TEMPLATE.format(
            n_water=n_water,  # Placeholder, update as needed
            species_solute_section=" " ,
            species_ion_section=species_ion_data_str,  # Placeholder, update as needed
            n_Au1 = n_gold//2,
            n_Au2 = n_gold//2,
            molecule_solute_section=" ",
            intra_section= " ",
            lj_solute_section=" ",
            lj_ion_section= lj_ion_section_str
        )
        

    # Write to output file
    with open('runtime.inpt', 'w') as f:
        f.write(runtime_content)

    print(f"runtime.inpt has been created.")

if __name__ == "__main__":
    main()



