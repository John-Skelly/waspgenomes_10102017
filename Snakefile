#!/usr/bin/env python3


import os
import re
import pathlib

#########
#GLOBALS#
#########

read_dir = 'data/reads'
meraculous_config_file = 'src/meraculous_config.txt'
read_set = ['norm', 'trim_decon']
k = ['31', '71', '127']
diploid_mode = ['0', '1']

###########
#Functions#
###########

def resolve_path(x):
    mypath=pathlib.Path(x).resolve()
    return str(mypath)

def find_completed_assemblies():
    my_files = list((dirpath, filenames)
                    for (dirpath, dirnames, filenames) 
                    in os.walk('output/meraculous'))
    my_fasta_files = []
    for dirpath, filenames in my_files:
        for filename in filenames:
            if ('final.scaffolds.fa' in filename 
                    and 'meraculous_final_results' in dirpath):
                my_path=os.path.join(dirpath, filename)
                my_fasta_files.append(resolve_path(my_path))
    return(my_fasta_files)

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
fastq_basenames = [os.path.basename(x).split('.')[0] for x in all_fastq_files]

all_samples = list(set(re.sub('^Ma-(?P<id>\w+)_\d+$', '\g<id>', x)
    for x in fastq_basenames)) 

# read the meraculous config
with open(meraculous_config_file, 'rt') as f:
    meraculous_config_string = ''.join(f.readlines())

#########
#Rules###
#########

#target
rule all:
    input:
        expand(('output/meraculous/{strain}/{read_set}/k_{k}/diplo_{diploid_mode}/'
                'meraculous_final_results/final.scaffolds.fa'),
               strain=all_samples, read_set=read_set, k=k, diploid_mode=diploid_mode)

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
        trim = 'output/trim_decon/Ma-{strain}_trim.log',
        repair = 'output/trim_decon/Ma-{strain}_repair.log'
    threads:
        10
    shell:
        'bin/bbmap/bbduk.sh '
        'threads={threads} '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        'ref=bin/bbmap/resources/phix174_ill.ref.fa.gz '
        'hdist=1 '
        'stats={output.f_stats} '       
        '2> {log.filter}'
        ' | '
        'bin/bbmap/bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'out=stdout.fastq '
        'ref=bin/bbmap/resources/adapters.fa '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '
        '| '
        'bin/bbmap/repair.sh '
        'in=stdin.fastq ' 
        'out={output.fq} '
        '2> {log.repair} '
        

# merge overlaping reads

rule merge:
    input:
        r1 = 'output/trim_decon/Ma-{strain}.fastq.gz'
    output:
        fq_merged = 'output/merge/Ma-{strain}_merged.fastq.gz',
        fq_unmerged = 'output/merge/Ma-{strain}_unmerged.fastq.gz',
        ihist = 'output/merge/Ma-{strain}_hist.txt'
    log:
        merge = 'output/merge/Ma-{strain}_merged.log'
    threads:
        10
    shell:
        'bin/bbmap/bbmerge.sh '
        'threads={threads} '
        'in={input.r1} '
        'verystrict=t '
        'out={output.fq_merged} '
        'outu={output.fq_unmerged} '
        'ihist={output.ihist} '
        '2> {log.merge} '

# normalise

rule norm:
    input:
        r1 = 'output/trim_decon/Ma-{strain}.fastq.gz'
    output:
        fq_norm = 'output/norm/Ma-{strain}.fastq.gz',
        fq_toss = 'output/norm/Ma-{strain}_toss.fastq.gz',
        hist = 'output/norm/Ma-{strain}_hist.txt',
        hist_out = 'output/norm/Ma-{strain}_hist_out.txt',
        peaks = 'output/norm/Ma-{strain}_peaks.txt'
    log:
        norm = 'output/norm/Ma-{strain}_norm.log'
    threads:
        50
    shell:
        'bin/bbmap/bbnorm.sh '
        'in={input.r1} '
        'threads={threads} '
        'out={output.fq_norm} '
        'outt={output.fq_toss} '
        'hist={output.hist} '
        'histout={output.hist_out} '
        'target=50 '
        'min=5 '
        'peaks={output.peaks} '
        '2> {log.norm} '
        
# run meraculous
rule meraculous:
    input:
        fastq = 'output/{read_set}/Ma-{strain}.fastq.gz'
    threads:
        50
    params:
        outdir = 'output/meraculous/{strain}/{read_set}/k_{k}/diplo_{diploid_mode}/'
    output:
        config = ('output/meraculous/{strain}/{read_set}/k_{k}/diplo_{diploid_mode}/'
                'config.txt'),
        contigs = ('output/meraculous/{strain}/{read_set}/k_{k}/diplo_{diploid_mode}/'
                'meraculous_final_results/final.scaffolds.fa')
    log:
        'output/meraculous/{strain}/{read_set}/meraculous.log'
    run:
        my_fastq = resolve_path(input.fastq)
        if wildcards.strain == 'MA3':
            if wildcards.k == '31':
                my_dmin = '46'
            elif wildcards.k == '71':
                my_dmin = '38'
            elif wildcards.k == '127':
                my_dmin = '26'
            else:
                raise ValueError('kmer not matched or whatever it does not really matter')
            if wildcards.read_set == 'norm':
                my_dmin = str(int(my_dmin)//2)
        else:
            my_dmin = '0'
        my_conf = meraculous_config_string.format(
            my_fastq, wildcards.k, wildcards.diploid_mode, my_dmin, threads)
        with open(output.config, 'wt') as f:
            f.write(my_conf)
        shell(
            'bin/meraculous/run_meraculous.sh '
            '-dir {params.outdir} '
            '-config {output.config} '
            '&> {log}')

#assembly stats
rule assembly_stats:
    input:
        fa = find_completed_assemblies
    output:
        'output/assembly_stats/stats.txt'
    run:
        my_inputfiles = ','.join(input.fa)
        shell('bin/bbmap/statswrapper.sh '
              'in={my_inputfiles} '
              'minscaf=1000 '
              'format=2 '
              '>{output} '
        
