#!/usr/bin/env python3


import os
import re
import pathlib

###########
#Functions#
###########

def resolve_path(x):
    mypath=pathlib.Path(x).resolve()
    return str(mypath)

def find_completed_assemblies(wildcards):
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

def parse_fasta_path(fasta_path):
    stripped_path = re.sub('^.+/output/meraculous/(?P<id>.+)/meraculous_.+$', '\g<id>', fasta_path)
    path_elements = stripped_path.split('/')
    return {'fasta': fasta_path, 
            'strain': path_elements[0], 
            'read_set': path_elements[1],
            'k': path_elements[2],
            'diploid_mode': path_elements[3]}

def dmin_writer(directory):
    stripped_path = re.sub('^.*output/meraculous/(?P<id>.*)/.*$', '\g<id>', directory)
    path_elements = stripped_path.split('/')
    if path_elements[0] == 'MA3':
        dmin_filepath = 'output/meraculous/' + stripped_path + '/meraculous_mercount/dmin.txt'
        if os.path.exists(dmin_filepath):
            with open(dmin_filepath) as f:
                return f.read()
        else:
            return '0'                
    else:
        return '0'

#########
#GLOBALS#
#########

read_dir = 'data/reads'
meraculous_config_file = 'src/meraculous_config.txt'
read_set = ['norm', 'trim_decon']
k = [ '63', '67', '71', '75', '79']
# k = ['31', '63', '67', '71', '75', '79', '127']
diploid_mode = ['0', '1']
augustus_config_dir = resolve_path('bin/augustus/config')
hymenoptera_odb = resolve_path('data/hymenoptera_data')
meraculous_threads = 50
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'
busco_container = 'shub://TomHarrop/singularity-containers:busco_3.0.2'
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

# generate a list of Busco targets
completed_assemblies = [parse_fasta_path(x) for x in find_completed_assemblies('')]
busco_targets = [('output/busco/'
                  '{strain}/{read_set}/{k}/{diploid_mode}/'
                  'run_busco/full_table_busco.tsv').format(**x) 
                 for x in completed_assemblies]
#########
#Rules###
#########

#target
rule all:
    input:
        expand(('output/meraculous/{strain}/{read_set}/k_{k}/diplo_{diploid_mode}/'
                'meraculous_final_results/final.scaffolds.fa'),
               strain=all_samples, read_set=read_set, k=k, diploid_mode=diploid_mode)

rule dmin_finder_targets:
    input:
        expand(('output/meraculous/{strain}/{read_set}/k_{k}/diplo_{diploid_mode}/'
                'meraculous_mercount/dmin.txt'),
               strain=all_samples, read_set=read_set, k=k, diploid_mode=diploid_mode)

rule kmer_coverage_targets:
    input:
        expand(('output/kmer_plot/Ma-{strain}_plot.pdf'),
               strain=all_samples)
                
rule meraculous_config_targets:
    input:
        expand(('output/meraculous/{strain}/{read_set}/k_{k}/diplo_{diploid_mode}/'
                'config.txt'),
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

#meraculous config
rule meraculous_config:
    input:
        fastq = 'output/{read_set}/Ma-{strain}.fastq.gz'
    threads:
        1
    params:
        outdir = 'output/meraculous/{strain}/{read_set}/k_{k}/diplo_{diploid_mode}/'
    output:
        config = ('output/meraculous/{strain}/{read_set}/k_{k}/diplo_{diploid_mode}/'
                'config.txt'),
    run:
        my_fastq = resolve_path(input.fastq)
        my_dmin = '0'
        my_conf = meraculous_config_string.format(
            my_fastq, wildcards.k, wildcards.diploid_mode, my_dmin, meraculous_threads)
        with open(output.config, 'wt') as f:
            f.write(my_conf)


# run meraculous
rule meraculous:
    input:
        fastq = 'output/{read_set}/Ma-{strain}.fastq.gz',
        config = ('output/meraculous/{strain}/{read_set}/k_{k}/diplo_{diploid_mode}/'
                'config.txt')
    threads:
        meraculous_threads
    params:
        outdir = 'output/meraculous/{strain}/{read_set}/k_{k}/diplo_{diploid_mode}'
    output:
        contigs = ('output/meraculous/{strain}/{read_set}/k_{k}/diplo_{diploid_mode}/'
                'meraculous_final_results/final.scaffolds.fa')
    log:
        'output/meraculous/{strain}/{read_set}/k_{k}/diplo_{diploid_mode}/meraculous.log'
    run:
        shell(
            'bin/meraculous/run_meraculous.sh '
            '-dir {params.outdir} '
            '-config {input.config} '
            '&> {log}')

#find dmin
rule dmin_finder:
    input:
        mercount_file = ('output/meraculous/{strain}/{read_set}/k_{k}'
                       '/diplo_{diploid_mode}/meraculous_mercount/mercount.hist')
    output:
        dmin_out = ('output/meraculous/{strain}/{read_set}/k_{k}/'
                  'diplo_{diploid_mode}/meraculous_mercount/dmin.txt'),
        dmin_plot = ('output/meraculous/{strain}/{read_set}/k_{k}/'
                   'diplo_{diploid_mode}/meraculous_mercount/dmin_plot.pdf')
    log:
        log = ('output/meraculous/{strain}/{read_set}/k_{k}/'
             'diplo_{diploid_mode}/meraculous_mercount/dmin_finder.log')
    singularity:
        r_container
    script:
        'src/dmin_finder.R'

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
              'format=3 '
              '>{output} ')

#kmer coverage analysis
rule plot_kmer_coverage:
    input:
        hist_before = 'output/norm/Ma-{strain}_hist.txt',
        hist_after = 'output/norm/Ma-{strain}_hist_out.txt',
        peaks = 'output/norm/Ma-{strain}_peaks.txt'
    output:
        plot = 'output/kmer_plot/Ma-{strain}_plot.pdf'
    threads:
        1
    log:
        log = 'output/kmer_plot/Ma-{strain}.log'
    script:
        'src/plot_kmer_coverage.R'

#busco analysis
rule busco_targets:
    input: busco_targets

rule busco:
    input:
        fasta = ('output/meraculous/{strain}/{read_set}/{k}/{diploid_mode}/'
                 'meraculous_final_results/final.scaffolds.fa')
    output:
        tsv = ('output/busco/'
               '{strain}/{read_set}/{k}/{diploid_mode}/'
               'run_busco/full_table_busco.tsv')
    params:
        wd = 'output/busco/{strain}/{read_set}/{k}/{diploid_mode}/'
    threads: 10
    run:
        my_fasta = resolve_path(input.fasta)
        shell(
              #'cd {params.wd} || exit 1 ; '
              './TomHarrop-singularity-containers-master-busco_3.0.2.simg '
              'run_BUSCO.py '
              '-i {my_fasta} '
              '-c {threads} '
              '-o busco '
              '-m geno '
              '-l {hymenoptera_odb} '
              '-s nasonia '
              '-f '                                                             #force
              '&> busco.log')

#Combine busco results
rule combine_busco_results:
    input: busco_targets = busco_targets
    output: 
        rds = 'output/busco/full_table_combine.Rds',
        plot = 'output/busco/stats_plot.pdf'
    script:
        'src/combine_busco_results.R'
