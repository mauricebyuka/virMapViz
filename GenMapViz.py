#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 13:24:25 2022

@author: Maurice
"""

import os
import argparse
from Bio import SeqIO
from scripts import VizModule

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sample",
                    dest="sample",
                    help="Sample Name. Will be the base name for result files.",
                    required=True)
parser.add_argument("-v", "--virus",
                    help="Name of the virus",
                    default='Some Virus')
parser.add_argument("-1", "--fastq1",
                    help="First pair of raw reads.",
                    default=None)
parser.add_argument("-2", "--fastq2",
                    help="Second pair of raw reads.",
                    default=None)
parser.add_argument("-u", "--unpaired",
                    help="Unpaired reads")
parser.add_argument("-r", "--refFile",
                    help="Path to the reference genome file (in fasta format).",
                    default=None)
parser.add_argument("-g", "--gff",
                    help="Path to the gene file (in GFF3 format). Can also be a tab delimited file with 3 columns: start, end, symbol.",
                    default=None)
parser.add_argument("-o", "--outdir",
                    help="Path to output directory.",
                    required=True)

parser.add_argument("-a", "--aligner",
                    help="Short read aligner to be used. options = bwa, bowtie2, bbmap. Default: bwa (bwa mem will be used).",
                    default = "bwa")

parser.add_argument("-d", "--dupl",
                    help="Indicate whether optical duplicated should be remove or kept. By default, duplicates are removed.",
                    choices=['remove', 'keep'],
                    default= 'remove')

parser.add_argument("-t", "--threads",
                    help="Number of threads. Default: 1 or half of available cpus.",
                    type=int,
                    default=None)

parser.add_argument("-qc", "--qualCntr",
                    help="Perfom quality control prior to alignment (fastp will be used).",
                    action='store_true')

parser.add_argument("-zs", "--zoom_start",
                    help="The start of the area to be zoomed. Default: 1/3 of the reference length.",
                    type=int,
                    default=None)

parser.add_argument("-ze", "--zoom_end",
                    help="The end of the area to be zoomed. Default: 2/3 of the reference length.",
                    type=int,
                    default=None)



args = parser.parse_args()


def main():
    sampleName = args.sample
    resultsFolder = args.outdir
    if not os.path.exists(resultsFolder):
        os.makedirs(resultsFolder)

    logFile = os.path.join(resultsFolder, sampleName + '_log.txt')
    if os.path.exists(logFile):
        os.remove(logFile)

    script_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'scripts')

    # Check for dependancy executables
    executables = ['Biopython', 'samtools', 'bcftools', 'bedtools']
    qual = args.qualCntr
    if qual:
        executables.append('fastp')

    mapTool = args.aligner
    if mapTool == 'bbmap':
        mapTool = 'bbmap.sh'
    executables.append(mapTool)

    missingDepend = VizModule.check_dependancy(executables, os.path.join(script_dir, 'checkLib.R'), logFile)

    if len(missingDepend) > 0:
        print('\nThe following programs were not found. You may want to install them or check that your environment is correctly set.\n')
        for program in missingDepend:
            print(f'\t{program}')
        print('\n')
    else:  # Run the program

        outValues = [sampleName, resultsFolder]  # Values to be used for plotting.

        # Set input files

        virusName = args.virus
        outValues.append(virusName)

        refFasta = args.refFile
        outValues.append(os.path.abspath(refFasta))

        gffFile = args.gff
        outValues.append(os.path.abspath(gffFile))

        unpairedReads = args.unpaired
        if unpairedReads:
            rawFastq1 = 'None'
            rawFastq2 = 'None'
        else:
            rawFastq1 = args.fastq1
            rawFastq2 = args.fastq2
            unpairedReads = 'None'
        

        # Set cpus to be used
        cpus = args.threads
        cpusAvail = int(os.cpu_count())
        if cpus != None:
            cpusUsed = cpus
        elif cpusAvail == 1:
            cpusUsed = cpusAvail
        else:
            cpusUsed = int(cpusAvail / 2)

        # Set bounds of the area to be zoomed out
        seqRecord = SeqIO.read(refFasta, 'fasta')
        refLen = len(seqRecord.seq)

        if args.zoom_start == None:
            zoomStart = int(refLen / 3)
        else:
            zoomStart = args.zoom_start

        if args.zoom_end == None:
            zoomEnd = int(refLen * 2 / 3)
        else:
            zoomEnd = int(args.zoom_end)

        # Show parameters to be used
        vz = VizModule.genViz(outDir=resultsFolder, outFileBase=sampleName,
                              fastq1=rawFastq1, fastq2=rawFastq2, unpairedFastq=unpairedReads, ref=refFasta,
                              gff=gffFile, cpus=cpusUsed, zs=zoomStart, ze=zoomEnd, logfile=logFile)
        print(vz)

        seqEnd = vz.refLen
        outValues.extend([zoomStart, zoomEnd, seqEnd])

        # Quality control & Alignement to the reference

        duplicateStatus = args.dupl
        if duplicateStatus == 'remove':
            if qual:
                trimmed = vz.trimmer()
                bamFile = vz.mapper(read1=trimmed[0], read2=trimmed[1], mapper=mapTool, rmdup=True)
            elif unpairedReads:
                bamFile = vz.mapper(read1=unpairedReads, read2='None', mapper=mapTool, rmdup=True)
            else:
                bamFile = vz.mapper(read1=rawFastq1, read2=rawFastq2, mapper=mapTool, rmdup=True)

        elif  duplicateStatus == 'keep':
            if qual:
                trimmed = vz.trimmer()
                bamFile = vz.mapper(read1=trimmed[0], read2=trimmed[1], mapper=mapTool, rmdup=False)
            elif unpairedReads:
                bamFile = vz.mapper(read1=unpairedReads, read2='None', mapper=mapTool, rmdup=False)
            else:
                bamFile = vz.mapper(read1=rawFastq1, read2=rawFastq2, mapper=mapTool, rmdup=False)

        outValues.append(os.path.abspath(bamFile))

        # Make the coverage bedGraph
        bedResults = vz.makeBed(bamFile)

        chromosome = bedResults[0]
        bedGraph = bedResults[1]
        outValues.append(os.path.abspath(bedGraph))
        outValues.append(chromosome)

        # Input to R:
        ValuesHandler = vz.collectVarNames(outValues)

        # Generating coverage files
        vz.generateGraphs(os.path.join(script_dir, 'GenViz.R'), ValuesHandler)


if __name__ == '__main__':
    main()
