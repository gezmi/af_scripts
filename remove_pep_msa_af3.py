#!/usr/bin/env python3
import argparse
import json
import os
from pathlib import Path
from typing import List, Union

def modify_json_file(input_path: Path, output_path: Path, chains=2, remove_templates=True) -> None:
    """
    Modify a JSON file by setting the 'unpairedMsa' field to an empty string
    in the specified nested structure.
    
    Args:
        input_path: Path to the input JSON file
        output_path: Path where the modified JSON will be saved
    """
    try:
        with open(input_path, 'r') as f:
            data = json.load(f)
        
        # Modify the specific nested field
        try:
            for ch in chains:
                data['name'] = f"{data['name']}_nomsa"
                data['sequences'][ch-1]['protein']['unpairedMsa'] = ""
                data['sequences'][ch-1]['protein']['pairedMsa'] = ""
                if remove_templates:
                    data['sequences'][ch-1]['protein']['templates'] = ""
                    data['name'] = f"{data['name']}_nomsa_notemp"
        except (KeyError, IndexError) as e:
            print(f"Warning: Could not find expected structure in {input_path}")
            print(f"Error details: {str(e)}")
            print('Typically, there was a problem with chain numbering')
            return
        
        # Create output directory if it doesn't exist
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Write modified data to output file
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
            
        print(f"Successfully processed: {input_path} -> {output_path}")
            
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in {input_path}")
        print(f"Error details: {str(e)}")
    except Exception as e:
        print(f"Error processing {input_path}")
        print(f"Error details: {str(e)}")

def process_files(input_files: List[Path], output_dir: Path, chains=[2], remove_templates=True) -> None:
    """
    Process multiple JSON files and save modified versions to the output directory.
    
    Args:
        input_files: List of paths to input JSON files
        output_dir: Directory where modified files will be saved
    """
    for input_path in input_files:
        output_name = input_path.stem.replace('_data', '') + "_no_pep_msa.json"
        output_path = output_dir / output_name
        modify_json_file(input_path, output_path, chains, remove_templates)

def main():
    parser = argparse.ArgumentParser(
        description="Modify JSON files by setting 'sequences[1].protein.unpairedMsa' to empty string"
    )
    parser.add_argument(
        'input_files',
        nargs='+',
        type=Path,
        help='Input JSON file(s) to process'
    )
    parser.add_argument(
        '-c', '--chains',
        type=str,
        default=2,
        help='Chain number(s) to remove MSA from, starting from 1. Concatenate with ,'
    )
    parser.add_argument(
        '-r', '--remove_templates',
        action='store_true',
        help='Also remove templates from selected chains, not just MSA.'
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        default=Path.cwd(),
        help='Output directory for modified files (default: current directory)'
    )
    
    args = parser.parse_args()
    
    # Validate input files
    valid_files = []
    for file_path in args.input_files:
        if not file_path.exists():
            print(f"Warning: Input file does not exist: {file_path}")
            continue
        if not file_path.is_file():
            print(f"Warning: Not a file: {file_path}")
            continue
        if file_path.suffix.lower() != '.json':
            print(f"Warning: Not a JSON file: {file_path}")
            continue
        valid_files.append(file_path)
    
    if not valid_files:
        print("Error: No valid input files to process")
        return
    
    chains = [int(x) for x in args.chains.split(',')]
    process_files(valid_files, args.output_dir, chains, args.remove_templates)

if __name__ == "__main__":
    main()
