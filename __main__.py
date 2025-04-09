#!/usr/bin/env python
# coding: utf-8
# Created by Wanlin Chen, 08.04.2025, Bochum

from pathlib import Path
import shutil

from config_converter import process_conversion
import topology_converter
from commands import parse_args 
from input_reader import read_input_file




def main():
    
    # Path to the script's directory
    script_dir = Path(__file__).resolve().parent
    data_dir = script_dir / 'database'

    # src_input = script_dir / 'database' / 'input.dat'
    # dst_input = Path.cwd() / 'input.dat'

    # # Copy the file
    # shutil.copy2(src_input, dst_input)
    
    args = parse_args()  # Get command-line arguments
    
    if args.config:
        process_conversion(args.config, data_dir/'gold_au.xyz', 'data.inpt')

    topology_converter.main()


if __name__ == "__main__":
    main()
