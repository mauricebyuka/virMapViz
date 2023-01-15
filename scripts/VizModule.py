#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 21:28:06 2022

@author: Maurice
"""
import os
import sys
import csv
import subprocess
from Bio import SeqIO


def check_dependancy(programs=[], rScript='', logfile=None):
    log = open(logfile, 'a')
    # the argument programs is a list of dependancies
    check = 'Checking dependancies:\n'
    log.write(check +'\n')
    print(check)
    missingPrograms = []

    # Check Biopython
    if 'Biopython' in programs:
        try:
            from Bio import SeqIO
            check = f' Biopython: OK ---'
            log.write(check +'\n')
            print(check)
        except ModuleNotFoundError:
            missingPrograms.append('Biopython')
            check =' Biopython: Module not found'
            log.write(check +'\n')
            print(check)

    # Check for system programs
    for program in programs:
        if not program == 'Biopython':
            track = subprocess.run(['which', program], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            #path = track.stdout.strip().decode('utf-8')
            if track.returncode != 0:
                check = f' {program}: Executable not found'
                log.write(check +'\n')
                print(check)
                missingPrograms.append(program)
            else:
                check = f' {program}: OK ---'
                log.write(check +'\n')
                print(check)

    # check for R libraries
    if not rScript == '':
        results = subprocess.run(['Rscript', rScript], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in results.stdout.decode('utf-8').strip().split('\n'):
            log.write(line + '\n')
            print(line)
            if 'Library not found' in line:
                lib = line.split(':')[0]
                missingPrograms.append(lib)
    log.close()
    return missingPrograms


class genViz:
    def __init__(self, outDir='None', outFileBase='None', fastq1='None',
                 fastq2='None', unpairedFastq='None', ref='None', gff='None', cpus=2, zs='None', ze='None', logfile='None'):
        self.outDir = outDir
        self.outFileBase = outFileBase
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.unpairedFastq = unpairedFastq
        self.ref = ref
        self.gff = gff
        self.cpus = cpus
        self.zs = zs
        self.ze = ze
        seqRecord = SeqIO.read(self.ref, 'fasta')
        self.refLen = len(seqRecord.seq)
        self.log = logfile

    def __str__(self):
        
        # params = [
        # f'\nParameters to be used:\n',
        # f' Basename: {self.outFileBase}',
        # f' OutDir: {os.path.abspath(self.outDir)}',
        # f' First read file: {os.path.abspath(self.fastq1)}',
        # f' Second read file: {os.path.abspath(self.fastq2)}',
        # f' Unpared reads file: {os.path.abspath(self.unpairedFastq)}',
        # f' Reference genome: {os.path.abspath(self.ref)}',
        # f' Interval file: {os.path.abspath(self.gff)}',
        # f' Zoom will start at: {self.zs}',
        # f' Zoom will end at: {self.ze}',
        # f' Threads to be used: {self.cpus}'
        # ]

        params = [
        f'\nParameters to be used:\n',
        f' Basename: {self.outFileBase}',
        f' OutDir: {os.path.abspath(self.outDir)}',
        f' First read file: {self.fastq1}',
        f' Second read file: {self.fastq2}',
        f' Unpared reads file: {self.unpairedFastq}',
        f' Reference genome: {os.path.abspath(self.ref)}',
        f' Interval file: {os.path.abspath(self.gff)}',
        f' Zoom will start at: {self.zs}',
        f' Zoom will end at: {self.ze}',
        f' Threads to be used: {self.cpus}'
        ]

        with open (self.log, 'a') as log:
            for param in params: 
                log.write(param + '\n')
                print(param)

        return ''

    def trimmer(self):
        with open (self.log, 'a') as log:
            message = '\n====== Trimming raw reads for quality control ======\n'
            log.write(message +'\n')
            print(message)
            newReadsFolder = os.path.join(self.outDir, 'trimmed_reads')
            if not os.path.exists(newReadsFolder):
                os.makedirs(newReadsFolder)
            
            if self.unpairedFastq != 'None':
                base = os.path.basename(self.unpairedFastq)
                newRead = os.path.join(newReadsFolder, f'trimmed_{base}')
                cmd = f'fastp -i {self.unpairedFastq} -o {newRead} -w {self.cpus} -j /dev/null -h /dev/null -q 20 -u 40 -l 15 -z 4'
                newRead1 = newRead
                newRead2 = 'None'
            else:
                base1 = os.path.basename(self.fastq1)
                base2 = os.path.basename(self.fastq2)
                newRead1 = os.path.join(newReadsFolder, f'trimmed_{base1}')
                newRead2 = os.path.join(newReadsFolder, f'trimmed_{base2}')
                cmd = f'fastp -i {self.fastq1} -I {self.fastq2} -o {newRead1} -O {newRead2} -w {self.cpus} -j /dev/null -h /dev/null -q 20 -u 40 -l 15 -z 4'
            if not os.path.exists(newRead1) and not os.path.exists(newRead2):
                execute = subprocess.run(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                for line in execute.stdout.decode('utf-8').strip().split('\n'):
                    log.write(line +'\n')

                if execute.returncode != 0:
                    print("Failed to trimming reads for quality. Check the log file for more details.")
                    sys.exit()

            else:
                message = '\n=============== TRIMMED READS EXIST ALREADY ========== SKIPPING!!!'
                log.write(message +'\n')
                print(message)
            outmessage = [
            '\nDone Trimming:\n',
            f'\t---> Trimmed read1: {newRead1}',
            f'\t---> Trimmed read2: {newRead2}\n',
            ]
            for message in outmessage:
                log.write(message + '\n')
                print(message)
            return [newRead1, newRead2]

    def mapper(self, read1=None, read2=None, mapper='bwa', rmdup=True):
        with open (self.log, 'a') as log:
            message1 = '\n====== Aligning reads to the reference genome ======\n'
            message2 = f'\t---> Reference length: {self.refLen:,}'
            log.write(message1 + '\n')
            print(message1)
            log.write(message2 + '\n')
            print(message2)
            rawBam = os.path.join(self.outDir, self.outFileBase + '_raw.bam')
            outBam = os.path.join(self.outDir, self.outFileBase + '_dedupped.bam')  # PCR duplicates removed
            
            if mapper == 'bwa':  
                cmd1 = f'bwa index {self.ref}'
                if self.unpairedFastq != 'None':
                    cmd2 = f'bwa mem -t {self.cpus} {self.ref} {read1} | samtools view -b | samtools sort -@{self.cpus} > {rawBam}'
                else:
                    cmd2 = f'bwa mem -t {self.cpus} {self.ref} {read1} {read2} | samtools view -b | samtools sort -@{self.cpus} > {rawBam}'
                
                

                if not os.path.exists(self.ref + '.pac'):
                    message = '\t---> Indexing the reference genome'
                    log.write(message + '\n')
                    print(message)

                    execute = subprocess.run(cmd1.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                    for line in execute.stdout.decode('utf-8').strip().split('\n'):
                        log.write(line + '\n')

                    if execute.returncode != 0:
                        print("Failed to index the reference. Check the log file for more details.")
                        sys.exit()

            elif mapper == 'bowtie2':
                cmd1 = f'bowtie2-build {self.ref} {self.ref}'
                if self.unpairedFastq != 'None':
                    cmd2 = f'bowtie2 -p {self.cpus} -x {self.ref} -U {read1} | samtools view -b | samtools sort -@{self.cpus} > {rawBam}'
                else:
                    cmd2 = f'bowtie2 -p {self.cpus} -x {self.ref} -1 {read1} -2 {read2} | samtools view -b | samtools sort -@{self.cpus} > {rawBam}'

                if not os.path.exists(self.ref + '.1.bt2'):
                    message = '\t---> Indexing the reference genome'
                    log.write(message + '\n')
                    print(message)

                    execute = subprocess.run(cmd1.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                    for line in execute.stdout.decode('utf-8').strip().split('\n'):
                        log.write(line + '\n')
                        
                    if execute.returncode != 0:
                        print("Failed to index the reference. Check the log file for more details.")
                        sys.exit()

            elif mapper == 'bbmap.sh':
                rawSam = os.path.join(self.outDir, self.outFileBase + '_raw.sam')
                if self.unpairedFastq != None:
                    cmd2a = f'bbmap.sh threads={self.cpus} ref={self.ref} in={read1} out={rawSam} nodisk'
                else:
                    cmd2a = f'bbmap.sh threads={self.cpus} ref={self.ref} in={read1} in2={read2} out={rawSam} nodisk'
                cmd2b = f'samtools view -b {rawSam} | samtools sort -@{self.cpus} > {rawBam}'

            cmd3 = f'samtools rmdup {rawBam} {outBam}'


            if not os.path.exists(outBam):
                message1 = f'\t---> Aligning reads with {mapper}'
                log.write(message1 + '\n')
                print(message1)

                if mapper == 'bbmap.sh':

                    execute2a = subprocess.run(cmd2a, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                    for line in execute2a.stdout.decode('utf-8').strip().split('\n'):
                        log.write(line + '\n')

                    if execute2a.returncode != 0:
                        print("bbmap failed. Check the log file for more details.")
                        sys.exit()

                    execute2b = subprocess.run(cmd2b, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                    for line in execute2b.stdout.decode('utf-8').strip().split('\n'):
                        log.write(line + '\n')

                    if execute2b.returncode != 0:
                        print("Samtools failed. Check the log file for more details.")
                        sys.exit()
                    os.remove(rawSam)
                else:

                    execute2 = subprocess.run(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                    for line in execute2.stdout.decode('utf-8').strip().split('\n'):
                        log.write(line + '\n')

                    if execute2.returncode != 0:
                        print(f"{mapper} failed. Check the log file for more details.")
                        sys.exit()

                message2 = '\t---> Removing PCR and optical duplicates'
                log.write(message2 + '\n')
                print(message2)               

                execute3 = subprocess.run(cmd3.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                for line in execute3.stdout.decode('utf-8').strip().split('\n'):
                    log.write(line + '\n')

                if execute3.returncode != 0:
                    print("Samtools failed. Check the log file for more details.")
                    sys.exit()

            else:
                message1 = '\n====== ALIGNMENT FILE EXISTS ALREADY ======== SKIPPING!!!'
                log.write(message1 + '\n')
                print(message1)
 
            if rmdup == False:
                dupInfo = '\nThe option to remove duplicates was turned off. The raw alignment file will be used. \n'
                log.write(dupInfo + '\n')
                print(dupInfo)
                
                outBam = rawBam

            outinfo = ['\nDone with mapping:\n', f'\t---> Alignment file: {outBam}']

            for info in outinfo:
                log.write(info + '\n')
                print(info)

        return outBam


    def makeBed(self, inBam):
        with open (self.log, 'a') as log:
            inmessage = '\n====== Making a coverage bedGraph file from the alignment file ======\n'
            log.write(inmessage + '\n')
            print(inmessage)

            bedGraph = os.path.join(self.outDir, f'{self.outFileBase}.bedGraph')
            if not os.path.exists(bedGraph):
                cmd = f'bedtools genomecov -ibam {inBam} -bga'
                with open(bedGraph, 'w') as bed:
                    execute = subprocess.run(cmd.split(), stdout=bed, stderr=subprocess.PIPE)
                    for line in execute.stderr.decode('utf-8').strip().split('\n'):
                        log.write(line + '\n')

                    if execute.returncode != 0:
                        print("Failed to make the coverage file from thre bam file. Check the log file for more details.")
                        sys.exit()
            else:
                message = '\n====== BEDGRAPH FILE EXISTS ALREADY ======== SKIPPING!!!\n'
                log.write(message + '\n')
                print(message)

            outmessage = f'\t---> bedGraph file: {bedGraph}'
            log.write(outmessage + '\n')
            print(outmessage)

            with open(bedGraph, 'r') as bed:
                chrom = bed.readline().split('\t')[0]
            return [chrom, os.path.abspath(bedGraph)]


    def collectVarNames(self, values):  # values is a list of values
        outFile = os.path.join(self.outDir, self.outFileBase + '_values.csv')
        with open(outFile, 'w') as vals:
            writer = csv.writer(vals, dialect=csv.excel_tab)
            writer.writerow(['sample', 'out_dir', 'virus', 're_fasta', 'gff', 'zoomStart', 'zoomEnd', 'seqEnd', 'bam', 'bedfile', 'chromosome'])
            writer.writerow(values)
        return(outFile)


    def generateGraphs(self, rScript, valueFile):
        with open (self.log, 'a') as log:
            inmessage = '\n====== Generating coverage graphs ======\n'
            log.write(inmessage + '\n')
            print(inmessage)

            execute = subprocess.run(['Rscript', rScript, valueFile], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            for line in execute.stdout.decode('utf-8').strip().split('\n'):
                log.write(line + '\n')

            if execute.returncode != 0:
                print("Failed to produce the coverage images. Check the log file for more details.")
                sys.exit()
            else:
                message = 'The following alignment graphs were generated:\n'
                log.write(message + '\n')
                print(message)

                file1 = os.path.join(self.outDir, f'graphs/{self.outFileBase}_full_coverage.png')
                # file2 = os.path.join(self.outDir, f'graphs/{self.outFileBase}_zoom_{self.zs}-{self.ze}_coverage.png')
                # file3 = os.path.join(self.outDir, f'graphs/{self.outFileBase}_full_alignment.png')
                # file4 = os.path.join(self.outDir, f'graphs/{self.outFileBase}_zoom_{self.zs}-{self.ze}_alignment.png')
                for file in [file1]:#, file2]:#, file4]:
                    if os.path.exists(file):
                        outmessage = f'\t---> {file}'
                        log.write(outmessage + '\n')
                        print(outmessage)
                    else:
                        log.write(f'Graph {os.path.basename(file)} was not generated! Something is wrong with input files\n')
                        print(f'Graph {os.path.basename(file)} was not generated! Check log file for more details')
                print('\n')
