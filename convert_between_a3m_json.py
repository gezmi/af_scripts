#!/usr/bin/env python
# coding: utf-8

"""
Convert A3M files to JSON format and vice versa.

This script provides functionality to convert A3M files into a structured JSON format for input at as AF3,
as well as to convert JSON files back into A3M format, for ColabFold input.

Author: Julia K. Varga
Email: julia.varga@mail.huji.ac.il

Hebrew University of Jerusalem
2025
"""

import os
import re
import json
import argparse

from tqdm.utils import disp_len


###### A3M TO JSON ######

def split_sequence_by_chains(sequence, lengths):
	"""
	Split a sequence into chains based on specified lengths.

	Args:
		sequence (str): The full protein sequence
		lengths (list): List of integer lengths for each chain

	Returns:
		list: List of split chain sequences
	"""
	result = []
	current_pos = 0
	
	for target_length in lengths:
		chain_seq = ""
		actual_length = 0
		
		# Keep adding characters until we reach the target length for this chain
		while actual_length < target_length and current_pos < len(sequence):
			char = sequence[current_pos]
			chain_seq += char
			
			# Only count uppercase letters toward the actual length
			if not char.islower():
				actual_length += 1
			
			current_pos += 1
		
		result.append(chain_seq)
	
	return result


def get_data_from_a3m(a3m_lines):
	"""
	Extract main metadata and sequences from A3M file lines.

	Args:
		a3m_lines (list): Lines from the A3M file

	Returns:
		tuple: (chains, dummy_chain_names, sequences)
	"""
	header = a3m_lines[0].strip().strip('#').split('\t')
	chains = list(map(int, header[0].split(',')))
	stoch = header[1].split(',')  # Stochasticity parameters
	
	dummy_chain_names = a3m_lines[1].strip().strip('>').split('\t')
	
	# Calculate total length as sum of all chain lengths
	total_length = sum(chains)
	
	full_seq = a3m_lines[2].strip()
	
	if len(full_seq) != total_length:
		raise ValueError(f"Length of full sequence {len(full_seq)} does not match total length {total_length}")
	
	sequences = split_sequence_by_chains(full_seq, chains)
	
	return chains, dummy_chain_names, sequences, header[0], stoch


def split_msa_vertically(chains, dummy_chain_names, a3m_lines):
	"""
	Split MSA entries vertically based on chain definitions.

	Args:
		chains (list): List of chain lengths
		dummy_chain_names (list): List of chain identifiers
		a3m_lines (list): Lines from the A3M file (excluding header)

	Returns:
		dict: Dictionary mapping chain names to their MSA entries
	"""
	# Create dictionary for each chain
	dict_a3m_per_chain = {chain: [] for chain in dummy_chain_names}
	
	for line in a3m_lines:
		if line.startswith('>'):
			# Add the header line to all chains
			for chain in dummy_chain_names:
				dict_a3m_per_chain[chain].append(line)
		else:
			# Split the line into corresponding chains
			split_sequences = split_sequence_by_chains(line.strip(), chains)
			
			for i, chain in enumerate(dummy_chain_names):
				# Add the sequence to the corresponding chain if it's not only gaps
				dict_a3m_per_chain[chain].append(split_sequences[i])
	
	return dict_a3m_per_chain


def create_one_job_server(name, sequences, stoch, seeds=[1], use_templates=True, max_temp_date='3000-01-01'):
	return  {
		"modelSeeds": seeds,
		"name": name,
		"dialect": 'alphafoldserver',
		"version": 1,
		"sequences": [
			{
				'proteinChain': {
					"sequence": seq,
					"useStructureTemplate": use_templates,
					"maxTemplateDate": max_temp_date,
					"count": int(stoch[i]) if stoch and i < len(stoch) else 1
				}
			}
			for i, seq in enumerate(sequences)
		]
	}


def create_local_job(name, sequences, msas, stoch, seeds=[1]):
	
	return  {
		"modelSeeds": seeds,
		"name": name,
		"dialect": 'alphafold3',
		"version": "2",
		"sequences": [
			{
				'protein': {
					"id": [chr(65 + i)],  # Single letter for each chain (A, B, C...)
					"sequence": seq,
					"unpairedMsa": ''.join(msas[i]) if i < len(msas) else [],
					"pairedMsa": '',
					"copies": int(stoch[i]) if stoch and i < len(stoch) else 1
				}
			}
			for i, seq in enumerate(sequences)
		]
	}


def create_json_structure(name, sequences, stoch=None, msas=None, add_path=False, server=False,
                          seeds=[1], use_templates=True, max_temp_date='3000-01-01'):
	"""
	Create the JSON structure with separate chains.

	Args:
		name (str): Name for the header
		sequences (list): List of chain sequences
		msas (list, optional): List of MSAs for each chain

	Returns:
		dict: JSON-compatible dictionary with chain information
	"""
	msas = [] if msas is None else msas
	
	if server:
		json_data = [create_one_job_server(name, sequences, stoch, seeds, use_templates, max_temp_date)]
	else:
		json_data = create_local_job(name, sequences, msas, stoch, seeds)

	return json_data


def write_chainwise_a3m(dict_a3m_per_chain, output_dir='.'):
	files = {}
	for chain, msa in dict_a3m_per_chain.items():
		if not os.path.exists(output_dir):
			os.makedirs(output_dir)
		
		# Write each chain's MSA to a separate file
		output_file = os.path.join(output_dir, f"{chain}.a3m")
		with open(output_file, 'w') as f:
			f.write('\n'.join(msa))
		
		# create full path
		output_file = os.path.abspath(output_file)
		files[chain] = output_file
	
	return files


def process_a3m_file(input_file, output_dir='.', add_path=False, suffix='', server=False,
                     seeds=[1], use_templates=True, max_temp_date='3000-01-01'):
	"""
	Process an A3M file and convert it to a structured JSON format.

	Args:
		input_file (str): Path to the input A3M file
		output_file (str, optional): Path to the output JSON file
									 (defaults to input filename with .json extension)

	Returns:
		dict: The generated JSON data structure
	"""
	
	basename = os.path.basename(os.path.splitext(input_file)[0])
	
	# Read the A3M file
	with open(input_file, 'r') as f:
		lines = f.readlines()
	
	# Strip lines and remove empty lines
	lines = [line.strip() for line in lines if line.strip()]
	
	# Extract data from the A3M file
	chains, dummy_chain_names, sequences, header_name, stoch = get_data_from_a3m(lines)
	
	# Split MSA entries vertically by chain
	dict_a3m_per_chain = split_msa_vertically(chains, dummy_chain_names, lines[1:])
	
	if add_path:
		# if we want smaller JSON files, we just write out the MSA-s, and only provide their path
		dict_a3m_per_chain = write_chainwise_a3m(dict_a3m_per_chain, output_dir)
	
	# Create output filename
	if suffix:
		basename = f'{basename}_{suffix}'
	
	# Create the JSON structure
	json_data = create_json_structure(basename, sequences, stoch, list(dict_a3m_per_chain.values()),
	                                  add_path, server, seeds, use_templates, max_temp_date)

	output_file = f"{basename}.json"
	output_file = os.path.join(output_dir, os.path.basename(output_file))
	
	# Write the JSON file
	with open(output_file, 'w') as f:
		json.dump(json_data, f, indent=4)
	
	print(f"Processed {input_file} -> {output_file}")


###### JSON TO A3M ######

def build_header(seqlens, nums_copies):
	header = f"#{','.join(map(str, seqlens))}\t{','.join(map(str, nums_copies))}"

	return header


def get_paired_headers(seqlens, dict_keys):
	# Get paired sequences: in their headers (MSA dict keys), there are tabs, 1 less than the len of seqlens
	paired_headers = [key for key in dict_keys if key.count('\t') == len(seqlens) - 1]

	return paired_headers


def process_chain(json_chain):
	num_copies = json_chain['protein']['copies'] if 'copies' in json_chain['protein'].keys() else len(json_chain['protein']['id'])
	sequence = json_chain['protein']['sequence']

	seqlen = len(sequence)

	# get the MSA
	unpaired_msa = json_chain['protein']['unpairedMsa'].split('\n')
	if len(unpaired_msa) == 1:
		if os.path.isfile(json_chain['protein']['unpairedMsa']):
			# read without new line
			with open(json_chain['protein']['unpairedMsa'], 'r') as f:
				unpaired_msa = f.readlines()
				unpaired_msa = [line.strip() for line in unpaired_msa[1:] if line.strip()]

	# split into paired and unpaired dictionaries based on 'E+' or 'E-' in header
	unpaired_dict = {}
	paired_dict = {}

	for i in range(0, len(unpaired_msa), 2):
		if i+1 < len(unpaired_msa):
			header = unpaired_msa[i]
			seq = unpaired_msa[i+1]
			if ('E+' in header or 'E-' in header) or '\t' not in header:
				unpaired_dict[header] = seq
			else:
				paired_dict[header] = seq

	return num_copies, sequence, seqlen, unpaired_dict, paired_dict


def build_paired_msa(paired_list, seqlens):
	all_keys = sorted(list(set(key for d in paired_list for key in d.keys())))

	# loop through keys, and search them in each of the dict of paired_list. if you can find it, extract the sequence.
	# if you cannot find it, pad it with the corresponding seqlens

	paired_msa = []
	for i, key in enumerate(all_keys):
		# get the sequence from each of the paired_dicts
		seqs = []
		for j, d in enumerate(paired_list):
			if key in d.keys():
				seqs.append(d[key])
			else:
				# pad with gaps
				seqs.append('-' * seqlens[j])

		# add the sequence to the unpaired_list
		paired_msa.append(key)
		paired_msa.append(''.join(seqs))

	return paired_msa


def build_unpaired_msa(seqlens, unpaired_list):
	unpaired_padded_msa = []

	for i, chain_msa in enumerate(unpaired_list):
		msa = chain_msa

		# sum chains before i, calculate padding with this length
		chain_len_before = sum(seqlens[:i])
		chain_len_after = sum(seqlens[i+1:])

		for header in msa.keys():
			seq = msa[header]

			# pad with gaps before and after
			padded_seq = '-' * chain_len_before + seq + '-' * chain_len_after

			unpaired_padded_msa.append(header)
			unpaired_padded_msa.append(padded_seq)

	return unpaired_padded_msa


def process_json_file(input_file, output_dir='.', suffix=''):
	with open(input_file, 'r') as f:
		json_data = json.load(f)
	
	basename = json_data['name']
	
	# extract necessary data from json file
	nums_copies, sequences, seqlens, unpaired_list, paired_list = zip(
		*[process_chain(chain) for chain in json_data['sequences']])
	
	# first line
	header = build_header(seqlens, nums_copies)
	
	# process paired and unpaired MSA separately
	paired_msa = build_paired_msa(paired_list, seqlens)
	unpaired_msa = build_unpaired_msa(seqlens, unpaired_list)
	
	# merge all the info
	full_msa = [header] + paired_msa + unpaired_msa
	
	# write to file, join with new lines
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	
	if suffix:
		basename = f'{basename}_{suffix}'
	
	output_file = f"{basename}.a3m"
	output_file = os.path.join(output_dir, os.path.basename(output_file))
	with open(output_file, 'w') as f:
		f.write('\n'.join(full_msa))
	
	print(f"Processed {input_file} -> {output_file}")


#####################

def main():
	"""
	Main function to parse command-line arguments and process the input files based on their extension
	"""
	parser = argparse.ArgumentParser(description='Convert A3M files to a JSON structure with separated chains.')
	parser.add_argument('input_file', type=str, help='Path to the input file(s)')
	parser.add_argument('--output_dir', '-o', help='Path to the output JSON file (default: input_name.json)')
	parser.add_argument('--add_path', '-p', action='store_true', help='Do not add the whole MSA, just the path, when the output is JSON file')
	parser.add_argument('--suffix', '-s', default='', help='append string to the output file name, before the extension')
	parser.add_argument('--server', '-r', action='store_true', help='aWhen converting to json, use the AF3 server format-')
	parser.add_argument('--seeds', '-e', default='1', type=str, help='Number of seeds for the AF3 server')
	parser.add_argument('--no_templates', '-n', action='store_true', help='Do not use templates on the AF3 server (default: use)')
	parser.add_argument('--max_temp_date', '-m', default='3000-01-01', help='Max template date when using server')
	args = parser.parse_args()
	
	# write to file, join with new lines
	if not os.path.exists(args.output_dir):
		os.makedirs(args.output_dir)
	
	if args.input_file.endswith('.a3m'):
		seeds = [int(seed) for seed in args.seeds.split(',')]
		
		if args.server and len(seeds) > 1:
			print('The server only supports one seeds. Aborting.')
			sys.exit(1)
		
		use_templates = not args.no_templates
		process_a3m_file(args.input_file, args.output_dir, args.add_path, args.suffix, args.server,
		                 seeds, use_templates, args.max_temp_date)
	elif args.input_file.endswith('.json'):
		process_json_file(args.input_file, args.output_dir, args.suffix)
	else:
		print("Error: Input file must be either an A3M or JSON file.")
		return


if __name__ == "__main__":
	main()
