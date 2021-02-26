#/bin/python3

import xml.etree.ElementTree as ET
from sys import argv
from os import path
import re
import argparse


parser = argparse.ArgumentParser(description='Parse one pepXML file(s) and convet to tsv, optionally applying an FDR filter.')
parser.add_argument('-f', '--file', type=str, required=True, help='Path to one pepXML files.')
parser.add_argument('-q', '--qvalue', action='store_true', help='Caculcate q-values.')
parser.add_argument('-c', '--fdr_cutoff', type=float, help='Apply a specified FDR cutoff.')
parser.add_argument('-d', '--decoy_prefix', type=str, default='rev_', help='The decoy prefix, used for calculating q-values.')
parser.add_argument('-o', '--output_file', type=str, required=True, help='The output filename.')

args = parser.parse_args()

def get_namespace(element):
	m = re.match(r'\{.*\}', element.tag)
	return m.group(0) if m else ''

pepxml_file = args.file

print('Loading pepXML file')
ET.register_namespace('', 'http://regis-web.systemsbiology.net/pepXML')
tree = ET.parse(pepxml_file)
root = tree.getroot()
ns = get_namespace(root)

f_out = args.output_file

print('Parsing pepXML file')

with open(f_out, 'w') as f:

	f.write('Spectrum\tScan\tRT\tPeptide\tSpectraST_Peptide\tProtein\tLabel\tiProphet_prob\n')

	for spectrum_query in root.iter('{}spectrum_query'.format(ns)):
		spectrum = spectrum_query.attrib['spectrum']
		scan = spectrum_query.attrib['start_scan']
		rt = spectrum_query.attrib['retention_time_sec']
		
		hit = spectrum_query.find('{}search_result'.format(ns)).find('{}search_hit'.format(ns))
		mod_info = hit.find('{}modification_info'.format(ns))
		pep = hit.attrib['peptide']
		protein = hit.attrib['protein']

		if protein.startswith(args.decoy_prefix):
			label = 'decoy'
		else:
			label = 'target'

		if mod_info is not None:
			spectrast_pep = mod_info.attrib['modified_peptide'] + '/' + spectrum_query.attrib['assumed_charge']
		else:
			spectrast_pep = hit.attrib['peptide'] + '/' + spectrum_query.attrib['assumed_charge']

		for analysis in hit.findall('{}analysis_result'.format(ns)):
			if analysis.attrib['analysis'] == 'peptideprophet':
				peptideprophet_prob = analysis.find('{}peptideprophet_result'.format(ns)).attrib['probability']
			else:
				iprophet_prob = analysis.find('{}interprophet_result'.format(ns)).attrib['probability']

		f.write('\t'.join([spectrum, scan, rt, pep, spectrast_pep, protein, label, iprophet_prob]) + '\n')

if args.qvalue:
	print('Adding q-values')

	with open(f_out, 'r') as f:
		header = f.readline().strip().split()
		contents = [x.strip().split() for x in f.readlines()]

	prob_index = header.index('iProphet_prob')
	label_index = header.index('Label')

	target_probs = [float(x[prob_index]) for x in contents if x[label_index] == 'target']
	decoy_probs = [float(x[prob_index]) for x in contents if x[label_index] == 'decoy']

	for spec in contents:
		prob = float(spec[prob_index])
		n_target = len([1 for x in target_probs if x >= prob])
		n_decoy = len([1 for x in decoy_probs if x >= prob])
		q = float(n_decoy) / (n_decoy + n_target)
		spec.append(str(q))
	with open(f_out, 'w') as f:
		header.append('q_value')
		f.write('\t'.join(header) + '\n')
		if args.fdr_cutoff:
			for line in [x for x in contents if float(x[-1]) <= args.fdr_cutoff]:
				f.write('\t'.join(line) + '\n')
		else:
			for line in contents:
				f.write('\t'.join(line) + '\n')
