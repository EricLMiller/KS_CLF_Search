#!/usr/bin/python3

import sys, os, re, string, argparse
from Bio.Seq import Seq
from difflib import SequenceMatcher


## Returns the reverse compliment of a sequence
## Will not translate lowercase, or many non-ATCG letters.
def reverseCompliment(myString):
	complement = {
		'A': 'T',
		'C': 'G',
		'G': 'C',
		'T': 'A',
		'R': 'Y',
		'Y': 'R',
		'S': 'S',
		'W': 'W',
		'K': 'M',
		'M': 'K',
		'B': 'V',
		'D': 'H',
		'H': 'D',
		'V': 'B',
		'N': 'N',
		'X': 'X',
		'-': '-',
		'.': '-'}
	return ''.join([complement[base] for base in myString[::-1]])


## Translates DNA to protein, using Bio.seq funtion
## However, many bacteria have alternate start codons.
## Therefore, if the first codon ends in 'TG', or is ATT, or is ATC, then treat this as an 'M'.  
def translate(sequence):
	proteinSequence = Seq(sequence).translate()
	if sequence[1:3] == 'TG' or sequence[0:3] == 'ATT' or sequence[0:3] == 'ATC':   # specific cases for alternate codons only
		proteinSeqArray = list(proteinSequence)	# Break protein into list
		proteinSeqArray[0] = 'M'				# Replace the first residue
		proteinSequence = "".join(proteinSeqArray)	# and re-join together the protein
	return(proteinSequence)




## Read the input file, which has the following structure:
## Identifier	Blast_start_nucleotide_position	Blast_end_nucleotide_position	Blast_protein_sequence
## Example"
## A43592	1117	1317	IVITGMGLVSVFGNDIDTFYNKLLEGESGISIIDRFDASSFSVRFGGQIRDFSSKGYIDGKNDRRLD
## This file has no header line
def readInput(blastFileName, minSizeBlastHit, permittedBlastCoverage):
	idToInterval = {}		## Key is ID; value is list of tuples, each as (blastStart, blastEnd)

	input = open(blastFileName, 'r')
	for line in input:
		tmpList = line.strip().split("\t")		## Strip off whitespace, and split based on tabs
		if len(tmpList) != 4 or not tmpList[1].isdigit() or not tmpList[2].isdigit():	## lines must be exactly 4 entries long, and have two numbers in the propper positions
			print("Error with line: " + line)	## Print error, but continue with the next line in the file
			continue

		[seqID, blastStart, blastEnd, AAseq] = [tmpList[0], int(tmpList[1]), int(tmpList[2]), tmpList[3][1:-1].replace('-', '')]		## temporarily store these variables.
		# We are getting rid of any gaps in the AA sequence returned by blast, as well as getting rid of any start codons and stop codons at the edges of the hit.
		# Start codons are removed, as translation can give a non-methionine residue if it is not an ATG start codon. 

		if (abs(blastStart - blastEnd) + 1) < minSizeBlastHit:
			print("Rejected ORF due to very short original BLAST hit: {0}\t{1}\t{2}".format(seqID, blastStart, blastEnd))
			continue 		# Blast hits need to have a minimum length.

		if '*' in AAseq:	# If there is a stop codon in the middle (as last residue was chopped off), take the first ORF that is longer than the permitted blast coverage proportion.
							# As this proportion is expected to be greater than 0.5, there will be only one such ORF.
							# This is how to handle BLAST hits such as ...TGGWED*SCD, where a stop codon is present within the BLAST hit.
			AAseqSegments = AAseq.split('*')
			longEnoughSegment = ""
			for segment in AAseqSegments:
				if len(segment) / float(len(AAseq)) >= permittedBlastCoverage:
					longEnoughSegment = segment
					break

			if not longEnoughSegment: # What if no segment is long enough, though?
				print("Rejected ORF due to stop codon within original BLAST hit: {0}\t{1}\t{2}".format(seqID, blastStart, blastEnd))
				continue

			sys.stdout.write("We trimmed entry due to stop codon within original BLAST hit. Old: {0}\t{1}\t{2}\t\tNew:".format(seqID, blastStart, blastEnd)) ## Report that we trimmed a BLAST hit
			segIndex = AAseq.find(longEnoughSegment) # Where is this shorten sequence?
			AAseq = longEnoughSegment
			if blastStart < blastEnd:	# Adjust the BLAST hit start and stop to relect this shortened segment
				blastStart += (3 * segIndex)
				blastEnd = blastStart + (len(AAseq) * 3)
			else:
				blastStart -= (3* segIndex)
				blastEnd = blastStart - (len(AAseq) * 3)
			sys.stdout.write("{0}\t{1}\t{2}\n".format(seqID, blastStart, blastEnd)) # Last bit of reporting that we trimmed a BLAST hit.

		if seqID not in idToInterval:	# sequence ID's can have more than one blast hit, so create a list to store the blast tuples
			idToInterval[seqID] = []

		idToInterval[seqID].append((blastStart, blastEnd, AAseq))	# Actually append this hit that we have been checking, as this is the variable to return

	input.close()
	return(idToInterval)



## Read in DNA sequences based on command line input file, but
## only if they are listed the blast input list. 
## This cuts down on RAM space needed (as the DNA sequence files can then be of any size), but does take more time to scan through the files
## Assumes fasta file organization with one header line per sequence
## idToInterval: passed in variable, keys are which sequence headers we want to retain.
def readSeqs(fastaFileName, idToInterval):
	idToSequence = {}			# Main variable that stores the DNA sequences
	idsOnly = set(idToInterval)	# Create new 'set' variable with all ids.

	fastaInput = open(fastaFileName, 'r')	# FASTA file with the DNA sequences	
	seqID = ""				# temp variable needed to remember id names
	seq = ""				# temp variable that stores sequences, given that sequences are more than 1 line long

	for line in fastaInput:	# Work through the file, line by line (whole file is not stored in memory)
		if line[0] == '>':	# We have hit a header line
			if seqID:		# Is this the end of a sequence that we are interested in?
				idToSequence[seqID] = seq	# If so, save the sequence,
				[seqID, seq] = ["", ""]		# reset these temporary variables

			possibleID = re.match(r'>([a-zA-Z0-9_]+)', line)	# Note that IDs from DNA must only have these characters; dots are not included. This prevents issues with mismatches between the BLAST hit results and the downloaded sequences, but this may cause problems in the future. 
			if possibleID and possibleID.group(1) in idsOnly:	# Are we looking for this header ID?
				seqID = possibleID.group(1)
				idsOnly.remove(seqID)		# Remove the ID from the list that we are looking for. This should result in a slow increase in speed if there are multiple instances of the header ID in the fasta file

		elif seqID:							# Not a header line, and we are looking for the DNA sequence
			seq += line.strip()				# Add this line of nucleotides to the sequence variable

	if seqID:								# Notice that this is a repeat of lines above. This is needed, as it is out of the main loop; this catches literally the last sequence in the input fasta file, if needed.  
		idToSequence[seqID] = seq 		

	fastaInput.close()						# Close file
	print("{0} sequences found; {1} sequences not found.".format(len(idToSequence), len(idsOnly)))		# Can be used for debugging
	return(idToSequence)					# Return the main variable. 


## This function returns all of the possible ORFs in a given sequence
## given an inital start and stop location to examine.
## The reverse compliment is used, if specified (as a boolean)
## What is an ORF, in this context? Any set of amino acids that:
## 		- is bounded by an initial in-frame stop codon (!!!)
##		- followed by any potential start codon
##		- terminating in a stop codon or the end of a sequence
## if no inital in-fram stop codon is possible, the first possible start codon
## for any given frame is used instead. 
def getORFs(sequence, initialOffset, endingOffset, reverse):
	startCodons = set(['ATG', 'GTG', 'TTG', 'ATT', 'CTG'])	# All of the potential start/stop codons
	stopCodons = set(['TAG', 'TAA', 'TGA'])
	ORFs = []		## found ORFs; list of tuples (absolute start position, absolute stop position, protein sequence)

	sequence = sequence.upper()		# Convert to uppercase, to standardize sequences
	if reverse:						# With the reverse boolean, flip and compliment the DNA sequence
		sequence = reverseCompliment(sequence)

	for frame in [0, 1, 2]:			# For each frame, one by one:
		startPositions = []			# List all start and stop codons, as in which nucleotide of the DNA sequence they are
		stopPositions = []
		ORFbyCodon = []				# temporary storage of ORFs, using the relative position of codons in the given sequence
									# list of tuples (start position, stop position, protein sequence)
		codons = [sequence[i:i+3] for i in range(frame, len(sequence), 3)]		# Split DNA sequence into codons
		for i in range(0, len(codons)):		# Search all codons for start and stop codons; record location
			if codons[i] in startCodons:	
				startPositions.append(i)
			elif codons[i] in stopCodons:
				stopPositions.append(i)

		if len(startPositions) == 0 or startPositions[0] != 0 :	# Start with an inital start codon
			startPositions.insert(0, 0)							# This then catches any ORF fragments at the start of the DNA sequence
		if len(stopPositions) == 0 or stopPositions[-1] != len(codons) - 1:		# End with a terminal stop codon for the same reasons
			stopPositions.append(len(codons) - 1)

		startPositions.reverse() 	# For the pop() algorithm below to work correctly, we need these arrays to be high numbers to low numbers.
		stopPositions.reverse()

		while len(startPositions) > 0 and len(stopPositions) > 0:	# Get all start/stop position pairs that are possible in this frame
			proteinSeq = translate(sequence[startPositions[-1] * 3 + frame:(stopPositions[-1] + 1) * 3 + frame])
			ORFbyCodon.append((startPositions.pop(), stopPositions.pop(), proteinSeq))

			while len(startPositions) > 0 and startPositions[-1] < ORFbyCodon[-1][1]:
				startPositions.pop()
			while len(stopPositions) > 0 and len(startPositions) > 0 and stopPositions[-1] < startPositions[-1]:
				stopPositions.pop()

		if reverse:		# Transfer ORFbyCodon into ORFs; 
						# Moves from relative position in a given DNA sequence to its absolute position
			for [codonStart, codonEnd, proteinSeq] in ORFbyCodon:
				ORFs.append((endingOffset - ((codonStart)* 3 + frame) , endingOffset - ((codonEnd + 1) * 3 + frame) + 1, proteinSeq))

		else:	# Transfer ORFbyCodon into ORFs; 
				# Moves from relative position in a given DNA sequence to its absolute position
				# +1 is needed below because offsets start with 0, and not 1. 
			for [codonStart, codonEnd, proteinSeq] in ORFbyCodon:
				ORFs.append((codonStart * 3 + initialOffset + frame + 1, (codonEnd + 1) * 3 + initialOffset + frame, proteinSeq))		

	return(ORFs)


## Out of a list of ORFs, which one is the 'best':
## Which one covers the returned blast sequence, in both position 
## and amino acid contents?
def findBestORF(ORFpositions, seqID, blastStart, blastEnd, AAseq):
	bestORF = ""		# Storage for the 'best' ORF

	# First, return an ORF that covers the BLAST hit entirely, 
	# and the entire AA sequence of the BLAST hit
	for (geneStart, geneEnd, proteinSeq) in sorted(ORFpositions):		
		if geneStart < geneEnd:
			if geneStart <= blastStart and geneEnd >= blastEnd and AAseq in proteinSeq:		# BLAST position is covered entirely; our ORF might be longer than the BLAST returned amino acid sequence
				bestORF = (geneStart, geneEnd, proteinSeq, 1)		# 1 indicating complete blast coverage
		else:	## Dependant on which direction the gene is
			if geneStart >= blastStart and geneEnd <= blastEnd and AAseq in proteinSeq:
				bestORF = (geneStart, geneEnd, proteinSeq, 1)		# 1 indicating complete blast coverage
	if bestORF:		## If there is an ORF that entirely covers the BLAST hit and the returned amino acid sequence, then we've found the best ORF.
		return(bestORF)


	# What if the BLAST hit is not contained within a viable ORF?
	# This can be the case if the BLAST hit goes beyond the start codon, or beyond the stop codon.
	# Both do occur, especailly while using a high e-value cutoff.
	# In this case, we want the ORF that covers the BLAST hit the most, both via position and in amino acid content.
	coverage = 0 		# Keep track of the coverage in the best BLAST hit was covered
						# This is a minimum value for amino acid sequence and BLAST position.
	for (geneStart, geneEnd, proteinSeq) in sorted(ORFpositions):	# For each potential ORF:
		[geneCoverage, seqCoverage] = [0, 0]	# These are two different ways of measuring coverage

		# This section focuses on BLAST position coverage
		if geneStart < geneEnd:
			geneCoverage = min(blastEnd, geneEnd) - max(blastStart, geneStart) + 1		# Number of bases of the BLAST positions that are covered
		else:
			geneCoverage = min(blastStart, geneStart) - max(blastEnd, geneEnd) + 1		# Same as above, but if the gene is in the reverse direction.
		geneCoverage /= float(abs(blastStart - blastEnd) + 1)	# Convert to a proportion of the BLAST hit covered

		# This section focuses on BLAST amino acid coverage
		seqMatch = SequenceMatcher(None, AAseq, proteinSeq ,autojunk=False)		# Find longest shared substring between BLAST return result and the ORF
		longestMatch = seqMatch.find_longest_match(0, len(AAseq), 0, len(proteinSeq))
		seqCoverage = longestMatch.size / float(len(AAseq))		# Convert to a proportion

		minCoverage = min(geneCoverage, seqCoverage)			# Take the minimum between the positional coverage and the amino acid coverage

		if minCoverage > coverage:		# If the coverage is greater than the previous 'best' ORF, replace it.
			bestORF = (geneStart, geneEnd, proteinSeq, minCoverage)
			coverage = minCoverage

	print("Blast hit not fully covered:\t{0}\t{1}\t{2}\t{3}".format(coverage, seqID, blastStart, blastEnd))	# Print a warning that the best ORF does not cover the entire BLAST hit.
	return(bestORF)		## This is then the best 


## Parse all command line arguments from user, or use default values
def getArguments():
	myArgParser = argparse.ArgumentParser(description='getORFs: get longest possible ORF that covers a protein BLAST hit', formatter_class=argparse.RawTextHelpFormatter)
	myArgParser.add_argument('-f', '--fasta', required=True, help='FASTA file with entire DNA sequence; potentially used with tblastn to obtain BLAST file')
	myArgParser.add_argument('-b', '--blast', required=True, help='protein BLAST file, with no inital header row')
	myArgParser.add_argument('-o', '--output', required=True, help='output file for this program')
	myArgParser.add_argument('-s', '--minSize', type=int, default=300, help='minimum size of the BLAST hit, in nucleotides. default of 300 nucleotides.')
	myArgParser.add_argument('-d', '--distance', type=int, default=500, help='distance before and after BLAST hit that should be scanned for full ORF. default of 500 nucleotides.')
	myArgParser.add_argument('-c', '--coverage', type=float, default=0.5, help='proportion of an ORF that must be covered by the BLAST hit. default of 0.5.')

	return(myArgParser.parse_args())


##################
## Main program start here
##################
if __name__=='__main__':
	myArgs = getArguments()
	fastaFileName = myArgs.fasta
	blastFileName = myArgs.blast
	outputFileName = myArgs.output
	minSizeBlastHit = myArgs.minSize
	distanceBeforeAfterBlast = myArgs.distance
	permittedBlastCoverage = myArgs.coverage

	# Read in BLAST hits; based on this, get the sequences required from the FASTA file specified.
	idToInterval = readInput(blastFileName, minSizeBlastHit, permittedBlastCoverage)
	idToSequence = readSeqs(fastaFileName, idToInterval)

	# If the output file already exists, exit out as to not modify it.
	if os.path.isfile("outputFileName"):
		sys.exit("Output file {0} already exists. Exiting.".format(outputFileName))

	# Create header row for the results.
	output = open("tester_ORF.txt", "w")
	output.write("\t".join(["seqID", "blastStart", "blastEnd", "blastCoverage", "geneStart", "geneEnd", "proteinLength", "proteinSeq"]))
	output.write("\n")

	# For each sequence in the BLAST file...
	for seqID in idToInterval:
		for [blastStart, blastEnd, AAseq] in idToInterval[seqID]:
			if seqID not in idToSequence:
				print("Cannot find sequence {0} in {1}".format(seqID, fastaFileName))
				break

			# Calculate the nucleotide interval for the sequence, from the blast results.
			# Notice that this is extended in both directions by the variable 'distanceBeforeAfterBlast'
			initalOffset = max(0, min(blastStart, blastEnd) - distanceBeforeAfterBlast - 1)
			endingOffset = min(len(idToSequence[seqID]), max(blastStart, blastEnd) + distanceBeforeAfterBlast)

			# And get the DNA sequence
			sequence = idToSequence[seqID][initalOffset:endingOffset]

			# Find the possible ORFs, while concidering which direction the BLAST hit is in.
			if blastStart < blastEnd: # one direction for the BLAST hit
				ORFpositions = getORFs(sequence, initalOffset, endingOffset, False)
			else: # and the other direction for the BLAST hit
				ORFpositions = getORFs(sequence, initalOffset, endingOffset, True)

			bestORF = findBestORF(ORFpositions, seqID, blastStart, blastEnd, AAseq)
			# If there is no ORF that is valid, report out.
			if not bestORF:
				print("No viable ORF found for {0}\t{1}\t{2}".format(seqID, blastStart, blastEnd))
				continue

			[geneStart, geneEnd, proteinSeq, blastCoverage] = bestORF[0:4]
			# If there is no ORF that is valid (that meats the minimum coverage, for example), report out.
			if blastCoverage < permittedBlastCoverage:
				print("Rejected ORF with {0:f} Blast coverage for {1}\t{2}\t{3}".format(blastCoverage, seqID, blastStart, blastEnd))
				continue

			# Write to output
			output.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(seqID, blastStart, blastEnd, blastCoverage, geneStart, geneEnd, len(proteinSeq), proteinSeq))

	output.close()


# Note: Bio.Seq will complain, but not quit, if given a DNA sequence to translate that is of length not divisible by three.
# 		This error message is: /usr/local/lib/python3.6/dist-packages/Bio/Seq.py:2338: BiopythonWarning: Partial codon, len(sequence) not a multiple of three. Explicitly trim the sequence or add trailing N before translation. This may become an error in future.
#  BiopythonWarning,
#		It's not great, as this is inserted at the end of a line, right after one of our status messages.

# Possible error messages (none of which quit the program automatically)
# 	Line 59:	Error with line:  + line
#					Issue with the BLAST file at this line, but continue with the next line in the file
#	Line 67: Rejected ORF due to very short original BLAST hit
#	Line 81: Rejected ORF due to stop codon within original BLAST hit
# 	Line 247: Blast hit not fully covered
# 					Print a warning that the best ORF does not cover the entire BLAST hit.
#	Line 293: Cannot find sequence... Refering to a BLAST hit not found in the FASTA file.
#	Line 313: No viable ORF found for ...
#	Line 319: Rejected ORF with xxx Blast coverage for ...

# Output file will state this:
# 	XX sequences found; XX sequences not found.
# This refers to which BLAST sequences were not found in the provided FASTA file.
