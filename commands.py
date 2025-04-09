import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Convert GROMACS input files to MetalWalls format.")

    parser.add_argument("-c", "--config", help="Convert a GROMACS .gro configuration file to MetalWalls format.")
    parser.add_argument("-t", "--topology", help="Convert a GROMACS .itp topology file to MetalWalls format.")
    parser.add_argument("-f", "--input", help="Input parameters to the program to specify the system") 

    return parser.parse_args()




# # Example of how to use this function
# input_file = 'input_parameters.txt'
# params = read_input_file(input_file)
# print(params)