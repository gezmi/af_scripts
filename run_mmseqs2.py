import logging
logger = logging.getLogger(__name__)
import os
from pathlib import Path
import argparse
import glob
import json
import types, sys

# Otherwise on some CPU-s, the code crashes with not being able to import jax
class MockModule(types.ModuleType):
    def __getattr__(self, name):
        return None

sys.modules['colabfold.alphafold.extra_ptm'] = MockModule('colabfold.alphafold.extra_ptm')
sys.modules['jax'] = MockModule('jax')
sys.modules['jax.numpy'] = MockModule('jax.numpy')

# This was the problem before the patch above
from colabfold.batch import *


def split_input_by_chains(sequence, lengths):
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
	header = a3m_lines[0].strip().strip('#').split('\t')
	chains = list(map(int, header[0].split(',')))
	
	dummy_chain_names = a3m_lines[1].strip().strip('>').split('\t')
	
	# sum chains for total length, with map
	total_length = sum(map(int, chains))
	
	full_seq = a3m_lines[2].strip()
	
	if len(full_seq) != total_length:
		raise ValueError(f"Length of full sequence {len(full_seq)} does not match total length {total_length}")
	
	sequences = split_input_by_chains(full_seq, chains)
	
	return chains, dummy_chain_names, sequences


def split_msa_vertically(chains, dummy_chain_names, a3m_lines):
	# create dictionary for each chain
	dict_a3m_per_chain = {chain: [] for i, chain in enumerate(dummy_chain_names)}
	
	for line in a3m_lines:
		if line.startswith('>'):
			[dict_a3m_per_chain[chain].append(line) for chain in dummy_chain_names]
		else:
			# split the line into the corresponding chains
			split_sequences = split_input_by_chains(line.strip(), chains)
			
			for i, chain in enumerate(dummy_chain_names):
				dict_a3m_per_chain[chain].append(split_sequences[i])
	
	return dict_a3m_per_chain

			
def create_a3m_mmseqs(sequences, jobname, use_templates):
	version = importlib_metadata.version("colabfold")
	commit = get_commit()
	if commit:
		version += f" ({commit})"
	user_agent = f"colabfold/{version}"
	
	(unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_features) = \
		get_msa_and_templates(query_sequences=sequences, jobname=jobname, pair_mode="unpaired_paired",
		                      msa_mode='mmseqs2_uniref_env',
		                      a3m_lines=None, result_dir=result_dir, use_templates=use_templates, custom_template_path=None,
		                      user_agent=user_agent)
	
	msa = msa_to_str(unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality)
	result_dir.joinpath(f"{jobname}.a3m").write_text(msa)  # write out the msa in the common format
	
	return msa


def add_to_json(json_data, dict_a3m_per_chain):
	template_paths = glob.glob(f"*/templates_*/")
	
	i = 101
	for l, chain in enumerate(json_data['sequences']):
		matching_paths = [path for path in template_paths if f'templates_{i}' in path]
		if matching_paths:  # Check if any paths matched
			json_data['sequences'][l]['protein']['templates'] = f"{os.getcwd()}/{matching_paths[0]}"
		else:
			json_data['sequences'][l]['protein']['templates'] = ''
		
		if 'protein' in chain.keys():
			json_data['sequences'][l]['protein']['unpairedMsa'] = '\n'.join(dict_a3m_per_chain[str(i)])
			i += 1
		
		# we need to null this, otherwise it attempts pairing and concatenates weirdly
		json_data['sequences'][l]['protein']['pairedMsa'] = ''
		
		
	return json_data

if __name__ == "__main__":
	
	# add arguments
	parser = argparse.ArgumentParser(description="Run MMseqs2 for an input JSON file of AF3")
	parser.add_argument("-i", "--input", type=str, required=True, help="Input JSON file")
	parser.add_argument("--use_templates", action="store_true", help="Use tempaltes")
	parser.add_argument("--pair_mode", type=str, help="Pairing mode: paired, unpaired, unpaired_paired (default)")
	
	args = parser.parse_args()
	
	result_dir=Path('.')
	
	# check if the input file exists
	if not os.path.exists(args.input):
		raise FileNotFoundError(f"Input file {args.input} does not exist.")
	
	# move json file to _org.json instead of .json
	json_file = args.input.replace('.json', '_org.json')
	if os.path.exists(json_file):
		os.remove(json_file)
	os.rename(args.input, json_file)
	
	# load the JSON file
	with open(json_file, 'r') as f:
		data = json.load(f)
	
	# get info from json
	jobname = data['name']
	sequences = [x['protein']['sequence'] for x in data['sequences'] if 'protein' in x.keys()]
	
	a3m_lines = create_a3m_mmseqs(sequences, jobname, args.use_templates)
	
	# get chains, dummy chain names and sequences
	a3m_lines = a3m_lines.split('\n')
	chains, dummy_chain_names, sequences = get_data_from_a3m(a3m_lines)
	
	# also put into unpairedMSA - this will not be paired in the pipeline, keep it as is
	dict_a3m_per_chain = split_msa_vertically(chains, dummy_chain_names, a3m_lines[1:])
	
	# add to json
	json_data = add_to_json(data, dict_a3m_per_chain)
	
	# write out the json, in place of the original
	json_file = args.input.replace('.json', '_paired.json')
	with open(args.input, 'w') as f:
		json.dump(json_data, f, indent=4)
