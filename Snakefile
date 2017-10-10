#!/usr/bin/env python3


import os
import re


#########
#GLOBALS#
#########

read_dir = 'data/reads'


###########
#Functions#
###########


#########
#Setup###
#########


#get a list of fastq files
read_dir_files = list((dirpath, filenames)
    for (dirpath, dirnames, filenames) 
    in os.walk(read_dir))

all_fastq_files = []

for dirpath, filenames in read_dir_files:
    for filename in filenames:
        if 'fastq.gz' in filename:
            all_fastq_files.append(os.path.join (dirpath, filename))

#list of unique samples
fastq_basenames = [os.path.basenames(x).split('.')[0] for x in all_fastq_files]

all_samples = list(set(re.sub('^MA-(?P<id>\w+)_\d+$', '\g<id>', x)
    for x in fastq_basenames)) 

#########
#Rules###
#########

#target
rule all:
    input:
        expand("output/trim_decon/Ma-{strain}.fastq.gz',
            strain=all_samples)

#trim & decontaminate read files

rule trim_decon:
    input:
        r1 = 'data/reads/Ma-{strain}_1.fastq.gz',
        r2 = 'data/reads/Ma-{strain}_2.fastq.gz'
    output:
        fq = 'output/trim_decon/Ma-{strain}.fastq.gz',
        f_stats = 'output/trim_decon/Ma-{strain}_filter-stats.txt',
        t_stats = 'output/trim_decon/Ma-{strain}_trim-stats.txt'
    log:
        filter = 'output/trim_decon/Ma-{strain}_filter.log',  
        trim = 'output/trim_decon/Ma-{strain}_trim.log'
    threads:
        10
    shell:
        'bin/bbmap/bbduk.sh '
        'threads={threads}'
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        'ref=bin/bbmap/resources/phix174_ill.ref.fa.gz '
        'hdist=1'
        'stats={output.f_stats}        
        '2> {log.filter}'
        '|'
        'bin/bbmap/bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'out={output.fq} '
        'ref=bin/bbmap/resources/adapters.fa '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '


