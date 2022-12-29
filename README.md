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
