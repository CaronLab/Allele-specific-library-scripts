#/usr/bin/python3

import subprocess
import os
import re
from pathlib import Path
from typing import List, Union
from multiprocessing import Pool
from datetime import datetime
from configparser import ConfigParser
from uuid import uuid4
import random
from itertools import islice


ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
common_aa = "ARNDCQEGHILKMFPSTWYV"

def chunk_list(it, size):
    it = iter(it)
    return iter(lambda: tuple(islice(it, size)), ())

class Job:
    def __init__(self,
                 command: Union[str, List[str]],
                 working_directory: Union[str, Path, None],
                 sample=None):

        self.command = command
        self.working_directory = working_directory
        self.returncode = None
        self.time_start = str(datetime.now()).replace(' ', '')
        self.time_end = ''
        self.stdout = ''
        self.stderr = ''
        self.sample = sample

    def run(self):
        if self.working_directory is not None:
            os.chdir(self.working_directory)

        command = self.command.split(' ') if isinstance(self.command, str) else self.command
        p = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        self.stdout, self.stderr = p.communicate()
        self.time_end = str(datetime.now()).replace(' ', '')
        self.returncode = p.returncode


def run(job: Job):
    job.run()
    return job


def _run_multiple_processes(jobs: List[Job], n_processes: int):
    pool = Pool(n_processes)
    returns = pool.map(run, jobs)
    pool.close()
    return returns


def remove_modifications(peptides: Union[List[str], str]):
    if isinstance(peptides, str):
        return ''.join(re.findall('[a-zA-Z]+', peptides))
    unmodified_peps = []
    for pep in peptides:
        pep = ''.join(re.findall('[a-zA-Z]+', pep))
        unmodified_peps.append(pep)
    return unmodified_peps


def remove_previous_and_next_aa(peptides: Union[List[str], str]):
    return_one = False
    if isinstance(peptides, str):
        peptides = [peptides]
        return_one = True
    for i in range(len(peptides)):
        if peptides[i][1] == '.':
            peptides[i] = peptides[i][2:]
        if peptides[i][-2] == '.':
            peptides[i] = peptides[i][:-2]
    if return_one:
        return peptides[0]
    return peptides


def replace_uncommon_aas(peptide):
    pep = peptide
    for aa in peptide:
        if aa not in common_aa:
            pep = pep.replace(aa, 'X')
    return pep


def create_netmhcpan_peptide_index(peptide_list):
    netmhcpan_peps = {}
    for i in range(len(peptide_list)):
        if len(peptide_list[i]) < 1:
            continue
        netmhc_pep = replace_uncommon_aas(peptide_list[i])
        netmhcpan_peps[peptide_list[i]] = netmhc_pep
    return netmhcpan_peps


class Helper:
    """
    example usage:
    cl_tools.make_binding_prediction_jobs()
    cl_tools.run_jubs()
    cl_tools.aggregate_netmhcpan_results()
    cl_tools.clear_jobs()
    """
    def __init__(self,
                 peptides: List[str] = None,
                 alleles: List[str] = ('HLA-A03:02', 'HLA-A02:02'),
                 n_threads: int = 0,
                 tmp_dir: str = '/tmp',
                 output_dir: str = None):
        """
        Helper class to run NetMHCpan on multiple CPUs from Python. Can annotated a file with peptides in it.
        """

        self.NETMHCPAN = 'netMHCpan'

        if isinstance(alleles, str):
            if ',' in alleles:
                alleles = alleles.split(',')
            elif ' ' in alleles:
                alleles = alleles.split(' ')
            else:
                alleles = [alleles]
        self.alleles = alleles
        self.min_length = 8
        self.peptides = peptides
        self.netmhcpan_peptides = dict()
        self.predictions = dict()
        self.wd = Path(output_dir) if output_dir else Path(os.getcwd())
        self.temp_dir = Path(tmp_dir) / 'PyNetMHCpan'
        if not self.wd.exists():
            self.wd.mkdir(parents=True)
        if not self.temp_dir.exists():
            self.temp_dir.mkdir(parents=True)
        self.predictions_made = False
        self.not_enough_peptides = []
        if n_threads < 1 or n_threads > os.cpu_count():
            self.n_threads = os.cpu_count()
        else:
            self.n_threads = n_threads
        self.jobs = []

    def add_peptides(self, peptides: List[str]):
        if not self.peptides:
            self.peptides = []
        peptides = remove_previous_and_next_aa(peptides)
        peptides = remove_modifications(peptides)
        use_peptides = [x for x in peptides if len(x) >= self.min_length]
        self.peptides += use_peptides

        if len(use_peptides) < len(peptides):
            print(f'{len(peptides)-len(use_peptides)} peptides were outside the length restrictions and '
                  f'were removed from the peptide list.')

        self.netmhcpan_peptides = create_netmhcpan_peptide_index(self.peptides)

        self.predictions = {pep: {} for pep in self.peptides}

    def _make_binding_prediction_jobs(self):
        if not self.peptides:
            print("ERROR: You need to add some peptides first!")
            return
        self.jobs = []

        # split peptide list into chunks
        peptides = list(self.netmhcpan_peptides.values())
        random.shuffle(peptides)  # we need to shuffle them so we don't end up with files filled with peptide lengths that take a LONG time to compute (this actually is a very significant speed up)

        if len(peptides) > 100:
            chunks = chunk_list(peptides, int(len(peptides)/self.n_threads))
        else:
            chunks = [peptides]
        job_number = 1

        for chunk in chunks:
            if len(chunk) < 1:
                continue
            fname = Path(self.temp_dir, f'peplist_{job_number}.csv')
            # save the new peptide list, this will be given to netMHCpan
            with open(str(fname), 'w') as f:
                f.write('\n'.join(chunk))
            # run netMHCpan
            command = f'{self.NETMHCPAN} -p -f {fname} -a {",".join(self.alleles)}'.split(' ')
            job = Job(command=command,
                      working_directory=self.temp_dir)
            self.jobs.append(job)
            job_number += 1

    def _run_jobs(self):
        self.jobs = _run_multiple_processes(self.jobs, n_processes=self.n_threads)

    def _clear_jobs(self):
        self.jobs = []

    def _aggregate_netmhcpan_results(self):
        for job in self.jobs:
            if job.returncode != 0:
                print(job.stdout)
                print(job.stderr)
                print('ERROR: There was a problem in NetMHCpan. See the above about for possible information.')
                exit(1)
            self._parse_netmhc_output(job.stdout.decode())

        #self.predictions.to_csv(str(Path(self.temp_dir) / f'netMHCpan_predictions.csv'))

    def _parse_netmhc_output(self, stdout: str):
        lines = stdout.split('\n')
        # works for NetMHCpan 4.0 and 4.1, will need to keep an eye on future releases
        allele_idx = 1
        peptide_idx = 2
        rank_idx = 12
        for line in lines:
            line = line.strip()
            line = line.split()
            if not line or line[0] == '#' or not line[0].isnumeric():
                continue
            allele = line[allele_idx].replace('*', '')
            peptide = line[peptide_idx]
            rank = line[rank_idx]

            if float(rank) <= 0.5:
                binder = 'Strong'
            elif float(rank) <= 2.0:
                binder = 'Weak'
            else:
                binder = 'Non-binder'

            self.predictions[peptide][allele] = {'rank': rank, 'binder': binder }

    def make_predictions(self):
        self.temp_dir = self.temp_dir / str(uuid4())
        self.temp_dir.mkdir(parents=True)
        self._make_binding_prediction_jobs()
        self._run_jobs()
        self._aggregate_netmhcpan_results()
        self._clear_jobs()

    def annotate_file(self, filename: str,
                      peptide_column: str = 'Peptide',
                      delimiter: str = '\t'):

        # clear the peptide list
        self.peptides = []

        with open(filename, 'r') as f:
            header = f.readline().strip().split(delimiter)
            content = [x.strip().split(delimiter) for x in f.readlines()]
        pep_index = header.index(peptide_column)
        peptides = [x[pep_index] for x in content]
        self.add_peptides(peptides)
        self.make_predictions()

        new_content = []

        for allele in self.alleles:
            header.append(f'{allele}_rank')
        for line in content:
            pep = replace_uncommon_aas(remove_modifications(remove_previous_and_next_aa(line[pep_index])))
            # pep = self.netmhcpan_peptides[pep]
            if len(pep) < self.min_length:
                continue
            for allele in self.alleles:
                rank = float(self.predictions[pep][allele]['rank'])
                #line.insert(pep_index-1, str(np.log(rank)))
                line.append(str(rank))
            new_content.append(line)

        f_out = self.wd / (Path(filename).stem + '_annotated.tsv')
        with open(f_out, 'w') as f:
            f.write('\t'.join(header) + '\n')
            for line in new_content:
                f.write('\t'.join(line) + '\n')

if __name__ == '__main__':
    from sys import argv
    from os import path
    pep_file = argv[1]
    alleles = argv[2:]
    output_dir = path.split(pep_file)[0]

    predictor = Helper(alleles=alleles, output_dir=output_dir)
    predictor.annotate_file(pep_file)
