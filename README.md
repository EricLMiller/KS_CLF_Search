# KS_CLF_Search — getORFs.py

This script finds all possible ORFs, defined as a start codon followed by a stop codon, for a region within a DNA sequence and returns the single best ORF that matches the query. As an example, BLASTp results can be fed into this program and the best full-length ORF will be returned. This then allows full genes to be pulled out of BLASTp results. 


What is an ORF, in this context? Any set of amino acids that:
- is bounded by an initial in-frame stop codon
- followed by any potential start codon
- terminating in a stop codon or the end of a sequence (in the case of draft genomes or small DNA fragments)
- minimum length of ORF (default = 300AAs)
If no inital in-frame stop codon is possible, the first possible start codon for any given frame is used instead. 


There are important caveats in this script. It specifically conciders alternative start codons (specifically, XTG, ATT, and ATC), which many, many times leads to a ORF longer than the actual gene on the 5' end. Additionally, only the local region surrounding the given start and end positions are scanned. By default, start and stop codons within 500 nucleotides of these given positions are used, but this default distance can be changed. 


This script also relies on BioPython for a small section of the code related to translating a sequence; as this is a base function of BioPython, the specific BioPython version number should be unimportant.


## Input files
The BLAST input file has no header line. It should consist of four tab-separated columns:
- Identifier
- Blast_start_nucleotide_position
- Blast_end_nucleotide_position
- Blast_protein_sequence

Example:
A43592	1117	1317	IVITGMGLVSVFGNDIDTFYNKLLEGESGISIIDRFDASSFSVRFGGQIRDFSSKGYIDGKNDRRLD

If the end_nucleotide_position is less than the start_nucleotide_position, then the gene is assumed to be running right-to-left and the reverse complement will be used. 

The DNA sequence input file is a FASTA, with the identifiers matching perfectly between the BLAST input file and the DNA input file. The DNA sequences are only read into memory if they are listed in the BLAST input file, a lower amount of RAM usage and any size of DNA sequence file. The script assumes that there is only one header line per sequence and that all headers are unique.




## Command line arguments
+ '-f', '--fasta', FASTA file with entire DNA sequence; potentially used with tblastn to obtain BLAST file
+ '-b', '--blast', protein BLAST file, with no inital header row
+ '-o', '--output', output file for this program
+ '-s', '--minSize', minimum size of the BLAST hit, in nucleotides. default of 300 nucleotides.
+ '-d', '--distance', distance before and after BLAST hit that should be scanned for full ORF. default of 500 nucleotides.
+ '-c', '--coverage', proportion of an ORF that must be covered by the BLAST hit. default of 0.5.



## Bio.Seq complaints
Note: Bio.Seq will complain, but not quit, if given a DNA sequence to translate that is of length not divisible by three.

This error message is: /usr/local/lib/python3.6/dist-packages/Bio/Seq.py:2338: BiopythonWarning: Partial codon, len(sequence) not a multiple of three. Explicitly trim the sequence or add trailing N before translation. This may become an error in future.

It's not great, as this is inserted at the end of a line, right after one of our status messages below.


## Possible error messages
(none of which quit the program automatically)
+ Line 59:	Error with line:  + line (Issue with the BLAST file at this line, but continue with the next line in the file)
+	Line 67: Rejected ORF due to very short original BLAST hit
+	Line 81: Rejected ORF due to stop codon within original BLAST hit
+ Line 247: Blast hit not fully covered (Print a warning that the best ORF does not cover the entire BLAST hit.)
+	Line 293: Cannot find sequence... Refering to a BLAST hit not found in the FASTA file.
+	Line 313: No viable ORF found for ...
+	Line 319: Rejected ORF with xxx Blast coverage for ...
+ The standard output will also state: XX sequences found; XX sequences not found. (This refers to which BLAST sequences were not found in the provided FASTA file.)


## Output file
Output file will have these tab-delimited columns:
- Identifier
- Blast_start_nucleotide_position
- Blast_end_nucleotide_position
- Proportion of BLAST sequence covered by this ORF
- Calculated ORF start position
- Calculated ORF end position
- ORF protein length
- ORF protein sequence





# KS_CLF_Search — MatchKS_CLS.py

Given a file with the format below (called "KS_CLF_Cutoff.txt" - see line 33 to change this):

\[first line as a header; all columns below separated by one tab]
- KS_Number
- KS_Subject
- Seq_ID
- Subject_Start
- Subject_End
- \[empty column]
- CLF_Number
- CLF_Subject
- Seq_ID
- Subject_Start
- Subject_End

The remainder of this file after this first line are the KS's and CLF's. A KS and a CLF on one line does not imply that these 2 genes have any relationship.



This program finds KS's and CLF's that share an identical 'Seq_IDs'. These KS's and CLF's must be a maximum of distance_Cutoff from each other to be concidered part of the same cluster (see line 34). For matching pairs of KS's and CLF's, the directionality of each gene in relation to each other is determined as either:

- KS first - with both genes in same direction
- CLF first - with both genes in same direction
- Diverge - genes are in different directions, pointing away from each other
- Converge - genes are in different diretions, pointing towards each other

This program only considers pairs of KS's and CLF's. Given KS_1...KS_2...CLF_1...CLF_2 all with the same Seq_ID, the following pairs will be reported, assuming all genes are within distance_Cutoff from each other:
- KS_1	CLF_1
- KS_2	CLF_1
- KS_1	CLF_2
- KS_2	CLF_2


Output is printed directly to standard output


# Example usage:

python matchKS_CLS.py > KS_CLF_Cutoff_Output.txt 

This re-directs the output to the output file named above.

Output is: KS_ID (from column 1 of input); CLF_ID (column 5); minimum distance between genes; directionality; start of KS gene; end of KS gene; start of CLF gene; end of CLF gene
