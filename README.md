# chimera_assembler.pl

Usage:
    chimera_assembler.pl -master 16S_20120506.fas -master COI_20120506.fas
    -optional Opsine1_20120506.fas -unit station

Options:
    -master <gene fasta file>, sequences that must be present

    -optional
            <gene fasta file>, optional sequences

    -min    minimum number of sequences from the optional genes [0]

    -unit   <indiv|station|species|motu>, defines the unit level at which
            the assembly [indiv]

    -criteria
            <bestindiv|bestchimera>, for any assembly unit higher than the
            individual, either choose the the best individual or the best
            chimera (ie. best set of sequences irrespective of the
            individual they belong to). The decision is made based on ACGT
            counts [bestindiv]

    -maxambiguity
            maximum percentage of ambiguity code a seq can have to be
            selected [5]

    -useW   use the sequences that have the working flag (W_). By default
            theses sequences are not used

    -outtab
             name of the output assembly [assembly_table.tsv]

    -selectedunits
            name of a file listing the units to keep, the others will
            discarded. The format should be the same as for the sequences
            (Species, Species|Station, ...)

    -exportseq
            export fasta files with the sequences selected in the assembly
            table

    -prefix prefix given to the output sequence files in addition to the
            original name [chimera_assembler_]

    -motu   <filename>, directs to a file that has at least two columns: 1)
            sequence name in E3S format (core or full), 2) MOTU name. If a
            single MOTU is found in a sp|location, the MOTU propagated to
            the whole sp|location.

    -noNCBI Do not use NCBI sequences

    -noTranscript
            Do not use RNA-seq contigs


