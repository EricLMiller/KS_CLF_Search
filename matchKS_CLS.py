#!/usr/bin/python

# EL Miller - 29 June 2021
# Given a file with this format (called "KS_CLF_Cutoff.txt" - see line 33):
# [[first line as a header]] KS Number	KS Subject Seq ID	Subject Start	Subject End	[empty column]	CLF Number	CLF Subject Seq ID	Subject Start	Subject End
# ... [[rest of lines are these KS's and CLF's]]
# [[a KS and a CLF on one line does not imply that these 2 genes have any relationship]]
#
# This program finds KS's and CLF's that share an identical 'Seq IDs'.
# These KS's and CLF's must be a maximum of distance_Cutoff from each other to be concidered part of the same cluster (see line 34)
# For matching pairs of KS's and CLF's, the directionality of each gene in relation to each other is determined as either:
#	KS first - with both genes in same direction
#	CLF first - with both genes in same direction
#	Diverge - genes are in different directions, pointing away from each other
#	Converge - genes are in different diretions, pointing towards each other
# This program only considers pairs. Given KS_1...KS_2...CLF_1...CLF_2 all with the same Seq ID, the following pairs will be reported,
# assuming all genes are within distance_Cutoff from each other:
#	KS_1	CLF_1
#	KS_2	CLF_1
#	KS_1	CLF_2
#	KS_2	CLF_2
#
# Output is printed directly to standard output
# Example usage:
#		python matchKS_CLS.py > KS_CLF_Cutoff_Output.txt 
# This re-directs the output to the output file named above.
#
# Output is: KS_ID (from column 1 of input); CLF_ID (column 5); minimum distance between genes; directionality; 
# (cont.) start of KS gene; end of KS gene; start of CLF gene; end of CLF gene
#


fileName = "KS_CLF_Cutoff.txt"		## input file, with fileName hardcoded into this program
distance_Cutoff = 2000				## maximum distance between any part (start or ends) of genes to still be considered a pair

inputFile = open(fileName, 'r')			## read in data a list of lines
lines = inputFile.readlines()

KS_to_Position = dict()				## set up a series of dictionaries to be modified. KS_id_to_Subject: key is column 1 in the input
CLF_to_Position = dict()			##		file (assumped to be unique) with an added 'KS_' in front of it. Value (called Subject) is then the 'Seq id'.	
KS_id_to_Subject = dict()			##		KS_to_Postion: key is this 'Seq id'. Value is a list of tuples: (KS_id, gene start, gene end)
CLF_id_to_Subject = dict()

for i in range(1, len(lines)):		## for every line, skipping the header line
	tmpArray = lines[i].split("\t")	## assume columns are separated by tabs
	if tmpArray[0]:					## examine KS Number
		KS_id = "KS_" + tmpArray[0]	## mark it as a KS_ number
		KS_Subject = tmpArray[1]
		KS_id_to_Subject[KS_id] = KS_Subject ## set this dict
		if KS_Subject not in KS_to_Position:	## if this Seq id is not in the dict, then add it, with this value being an empty list
			KS_to_Position[KS_Subject] = []
		KS_to_Position[KS_Subject].append((KS_id, int(tmpArray[2]), int(tmpArray[3])))	## add the tuple of (KS_id, gene start, gene end) to this list.

	if tmpArray[5]:					## do the same, if there is a CLF also on this line.
		CLF_id = "CLF_" + tmpArray[5]
		CLF_Subject = tmpArray[6]
		CLF_id_to_Subject[CLF_id] = CLF_Subject 
		if CLF_Subject not in CLF_to_Position:
			CLF_to_Position[CLF_Subject] = []
		CLF_to_Position[CLF_Subject].append((CLF_id ,int(tmpArray[7]), int(tmpArray[8])))

inputFile.close()					## close up the input file	


## now that the relevant data is stored in these 4 dictionaries, start analysing the data.

for key in sorted(KS_to_Position):	## for every KS_id. There's no analog to the CLF_id, as we are looking for matches between the two gene types.
	if key in CLF_to_Position:		## if the Seq id does NOT have a CLF, continue to the next KS. There is no match. Otherwise, continue in this loop.
		for i in KS_to_Position[key]:	## for all KS's within this Seq id...
			KS_id = i[0]
			KS_Start = i[1]
			KS_End = i[2]

			for j in CLF_to_Position[key]:	## for all CLF's within this same Seq id...
				CLF_id = j[0]
				CLF_Start = j[1]
				CLF_End = j[2]

				if min(KS_Start, KS_End) < min(CLF_Start, CLF_End):			## This checks that the KS is first
					distance = min(CLF_Start, CLF_End) - max(KS_Start, KS_End)	## Minimum distance between genes, again knowing that KS is first.
					if distance <= distance_Cutoff:							## Are we within the distance cutoff? If not, this is not a match
						direction = ""										## We now know that we have a match. What direction is it?
						if KS_Start < KS_End:								## Set of if-elses to check directionality
							if CLF_Start < CLF_End:
								direction = "KS_First"
							else:
								direction = "Converge"
						else:
							if CLF_Start < CLF_End:
								direction = "Diverge"
							else:
								direction = "CLF_First"
						## Print to standard output, separated into columns by tabs
						print("\t".join([KS_id, CLF_id, key, str(distance), direction, str(KS_Start), str(KS_End), str(CLF_Start), str(CLF_End)]))


				else:															## If the KS is not first, then the CLF must be first.
					distance = min(KS_Start, KS_End) - max(CLF_Start, CLF_End)	## Minimum distance between genes, again knowing that CLF is first.
					if distance <= distance_Cutoff:								## Are we within the distance cutoff? If not, this is not a match
						direction = ""											## We now know that we have a match. What direction is it?
						if KS_Start < KS_End:									## Set of if-elses to check directionality
							if CLF_Start < CLF_End:
								direction = "CLF_First"
							else:
								direction = "Diverge"
						else:
							if CLF_Start < CLF_End:
								direction = "Converge"
							else:
								direction = "KS_First" 
						## Print to standard output, separated into columns by tabs
						print("\t".join([KS_id, CLF_id, key, str(distance), direction, str(KS_Start), str(KS_End), str(CLF_Start), str(CLF_End)]))

