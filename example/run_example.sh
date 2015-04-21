#!/bin/sh

export basedir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
export bin_dir=$(dirname ${basedir})/bin

${bin_dir}/ISmapper.pl --genome UH02_contigs_annotation.fa --is ISAba1.fasta --reference ACICU.fasta --out output 
