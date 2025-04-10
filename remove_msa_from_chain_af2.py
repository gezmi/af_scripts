#!/usr/bin/env python3
import argparse
import os
import sys
from pathlib import Path


def main():
	parser = argparse.ArgumentParser(
		description="Remove a chain's MSA from an a3m file"
	)
	parser.add_argument(
		'input_msa',
		type=Path,
		help='Input a3m file to process'
	)
	parser.add_argument(
		'--chain', '-c',
		type=int,
		help='Chain number to remove (starting from 1)'
	)
	parser.add_argument(
		'-o', '--output-dir',
		type=Path,
		help='Output directory for modified files (default: same as input file)'
	)
	
	args = parser.parse_args()
	
	input_msa = args.input_msa
	peptide = args.chain
	outdir = args.output_dir
	
	# open a file and iterate over its lines
	new_lines = ''
	lengths = [1]
	with open(input_msa) as f:
		lines_iter = iter(f)
		skip = False
		edit_msa = False
		for line in lines_iter:
			# print(skip, edit_msa, '\t' not in line, line)
			if line.startswith('#'):
				# print('header')
				# remove #, split on whitespace
				lengths, _ = line[1:].split()
				lengths = lengths.split(',')
				
				if args.chain > len(lengths) or args.chain == 0:
					print(f'Chain number {args.chain} is larger than the number of chains in the file or 0')
					sys.exit(1)
					
				new_lines += line
			elif '>101\t102' in line:
				# print('first paired')
				if len(lengths) > 2:  # we need to keep the rest, just pad the peptide
					edit_msa = True
				else:  # we can omit the paired msa
					skip = True
				new_lines += line
				
				# add also next line, safely using the iterator
				try:
					next_line = next(lines_iter)
					new_lines += next_line
				except StopIteration:
					# Handle case where there's no next line
					pass
			elif line.startswith(f'>10{peptide}'):  # remove peptide entry
				# print('peptide')
				skip = True
				new_lines += line
				# add also next line, safely using the iterator
				try:
					next_line = next(lines_iter)
					new_lines += next_line
				except StopIteration:
					# Handle case where there's no next line
					pass
			elif line.startswith('>') and '\t' not in line:  # keep MSA for the rest
				# print('other chains')
				skip = False
				edit_msa = False
				new_lines += line
			elif edit_msa and not line.startswith('>'):
				# print('edit msa')
				line_to_add = ''
				for i, l in enumerate(lengths):
					if i != peptide - 1:
						# loop through the line and add everything, until you didnt count l capital and - characters
						# then remove the same amount of characters from the line
						capital_dash = 0
						for c in line:
							line_to_add += c
							if c.isupper() or c == '-':
								capital_dash += 1
							if capital_dash == int(l):
								break
						line = line[capital_dash:]
					else:
						line_to_add += '-' * int(l)
				new_lines += f'{line_to_add}\n'
			elif not skip:
				# print('keep')
				new_lines += line
			else:
				# print('skip')
				pass
	
	# Generate output path
	if outdir:
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		
		basename = os.path.basename(str(input_msa))
		outpath = os.path.join(outdir, f'{basename.replace(".a3m", "_nomsa.a3m")}')
	else:
		outpath = str(input_msa).replace('.a3m', '_nomsa.a3m')
	
	with open(outpath, 'w') as f:
		f.write(new_lines)
	
	print(f'New a3m file written to {outpath}')


if __name__ == "__main__":
	main()