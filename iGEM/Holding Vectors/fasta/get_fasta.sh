#!/bin/bash

roots=(3 5)
suffices=(A B C D)

for root in "${roots[@]}"; do
    for suffix in "${suffices[@]}"; do
        plasmid_name="pSB1C${root}${suffix}"

        wget -O "${plasmid_name}.fasta" https://parts.igem.org/cgi/partsdb/composite_edit/putseq.cgi\?part\=${plasmid_name}

        plasmid_name="pSB1C${root}S${suffix}"
        wget -O "${plasmid_name}.fasta" https://parts.igem.org/cgi/partsdb/composite_edit/putseq.cgi\?part\=${plasmid_name}
    done
done
