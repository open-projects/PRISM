#    'frameshift': 'CDS,UTR5,OffFrame,UTR3,ncRNA,Frameshift,Intronic,Intergenic',
#    'prio1': 'Extra,CDS,UTR5,OffFrame,UTR3,ncRNA,Intronic,Intergenic',
#    'prio2': 'CDS,Extra,UTR5,OffFrame,UTR3,ncRNA,Intronic,Intergenic',
#    'prio3': 'CDS,UTR5,OffFrame,UTR3,ncRNA,Extra,Intronic,Intergenic',

# if HLA files are placed in the sample directories:
#./prism.batch.py -g h.ens90 -extra input/Mel15OP1_mut_only_new_header.fasta -i input -o output -cat frameshift prio2 prio3

./prism.batch.py -g h.ens90 -hla 17T_input/17T/17T.hla -extra 17T_input/17T/17T.avi.final.fasta -i 17T_input -o 17T_output -cat frameshift prio2 prio3 -r 17T.run_prism.sh
