import Bio as bp # import the biopython module
from Bio import SeqIO, SeqUtils

#read the fasta sequence from file
handle =  open("sequence.fasta")
genome = SeqIO.read(handle,"fasta")
print("The analysis genome ID is:",genome.name,"\nThe description is:",genome.description)
genome_text = str(genome.seq) # get the sequence as a string from the genome SeqRecord object

handle.close()


# AnalyzeBacterialGenomeis our highest level function for analyzing an input bacterial genome.
# input: a string (sense) representing a DNA sequence from one strand of the genome, 5' to 3'
# output: a list (senseORFs) containing all the ORFs found within the genome, a list (antisenseORFs), a dictionary (proteinsSenseDict) whose keys are integers from 1 to the number of proteins found, and values are a 2 element list containing the protein aa sequence and the molecular mass of that sequence, and a dictionary (proteinsAntisenseDict)
def AnalyzeBacterialGenome(sense = ""):
    antisense5to3 = GetReverse(GetComplement(sense))

    senseORFs = GetORFs(sense) # user defined function that scans over an input DNA sequence and returns all the ORFs it contains
    antisenseORFs = GetORFs(antisense5to3)

    over50_senseORFs = FilterByLength(senseORFs,50) # user defined function that returns ORFs with length over n (the input integer)
    over50_antisenseORFs = FilterByLength(antisenseORFs,50)

    # get the proteins
    proteins50aa_sense = PerformTranslations(over50_senseORFs)
    proteins50aa_antisense = PerformTranslations(over50_antisenseORFs)

    # get the dictionarys to return
    proteinsAndMW_sense = GetProteinDictionary(proteins50aa_sense)
    proteinsAndMW_antisense = GetProteinDictionary(proteins50aa_antisense)

    return proteinsAndMW_sense, proteinsAndMW_antisense
    



# translate each element of a list where each element is a string representing a DNA sequence
def PerformTranslations(DNAseqList):
    translatedList = []
    for i in range(0,len(DNAseqList)):
        translatedList.append(translate(DNAseqList[i]))
    return translatedList

# generate the MW dictionary for a list of proteins
# use Biopython utility to get molecular weights of all the proteins
# https://biopython.org/docs/1.75/api/Bio.SeqUtils.html
# check some calculated MW yourself: https://www.bioinformatics.org/sms/prot_mw.html
# Biopython module for molecular weights is in Da, so divide by 1000 to get kDa
def GetProteinDictionary(proteinList):
    # initialize dictionary with Keys "Protein x" pointing to a list whose 0th element is the final protein sequence, and 1th element is the MW of the protein in kDa
    key = "protein "
    proteinDictionary = dict() # initialize the dictionary to return

    # range over all the proteins make a key for each of them called "Protein i" and assign value as the sequence for the protein and the MW of the final protein in kDa
    for i in range(0,len(proteinList)): # for all the proteins in the list
        i_num = str(i + 1)
        keyI = key + i_num # the ith key for the protein
        proteinDictionary[keyI] = [proteinList[i], bp.SeqUtils.molecular_weight(proteinList[i], "protein")/1000]

    return proteinDictionary


def GetComplement(someDNA): 
    compDNA_Arr = [] # initialize a list aka an array
    for i in range(0,len(someDNA)): # loop over all the indices of someString, starting at the last index and decrementing down
        cur_nt = someDNA[i]
        if cur_nt == "A":
            compDNA_Arr.append("T")
        elif cur_nt == "T":
            compDNA_Arr.append("A")
        elif cur_nt == "C":
            compDNA_Arr.append("G")
        elif cur_nt == "G":
            compDNA_Arr.append("C")
        else:
            compDNA_Arr.append(cur_nt) # this means some character besides A, T, C, G given - return itself.
    
    compDNA = "".join(compDNA_Arr)
    return compDNA

def GetReverse(someString): # you can do this in one line with someString[::-1]
    reverseStringArr = [] # initialize a list aka an array
    for i in range(len(someString)-1,-1,-1): # loop over all the indices of someString, starting at the last index and decrementing down
        reverseStringArr.append(someString[i])
    
    reverseString = "".join(reverseStringArr)
    return reverseString


def GetORFs(sequence5to3):
    ORFs = [] # initialize list of ORFs
    ranOnce = False
    
    # range over the whole sequence looking at each possible ORF. scan the sequence starting at some position i, then check if the the 3 nucleotides from sequence5to3[i:i+3] define a start codon. if so, continue, otherwise don't. if so, range over the rest of the sequence until a stop codon is found. note: this transcript is what will be read during translation.
    # how do we know when to stop the for loop? assuming the sequence is 3 or longer, we know that for any length of sequence that the last two nucleotides wont have a full frame. e.g. ATCG will have ATC, then TCG, but then nothing else for C or G. therefore, range over the sequence up to length of sequence minus 2
    seqLength = len(sequence5to3)
    for i in range(seqLength - 2):
        firstCodon = sequence5to3[i:i+3] # get the first potential codon at this position

        if firstCodon == "ATG": # if this a start codon in the coding strand, then range over the rest of the sequence until a stop is reached. append this.
            j = i+3 # initialize an inner variable for ranging over the rest of the sequence starting at position i
            nextCodon = sequence5to3[i+3:i+6] # initalize nextCodon as second codon in reading frame

            # range over the rest of the sequence until a stop codon is found OR the end of the sequence has been reached
            while CheckIfStop(nextCodon) and j < seqLength:
                
                j += 3 # get next codon index frame
                nextCodon = sequence5to3[j:j+3] # get the next codon
                
            
            if j >=seqLength and ranOnce == False:
                ranOnce = True # no need to run this again
                print("\nsome ORFs found where start codon is present but no stop codon is present, so the end of the genome was reached! These were excluded from the potential ORF collection\n")
            
            if j <seqLength:
                ORFs.append(sequence5to3[i:j+3]) # add the new ORF to list ORFs. note, the while loop sbotS once a stop codon has reached so we have to add the stop codon to the end

    return ORFs


# CheckIfStop checks if some 3 length sttring representing a DNA codon (someCodon) equals one of the stop codons (in terms of DNA, so NO uracil!)
def CheckIfStop(someCodon):
    check = True
    if someCodon == "TAA" or someCodon == "TAG" or someCodon == "TGA":
        check = False
    return check

# iterates over some list of iterables (someList, a list which contains lists and/or strings) and removes elements whose length are less than input length (exclusionFactor) and returns a new list (filteredList) with all the elements whose lengths are greater than or equal to exclusionFactor.
def FilterByLength(someList, exclusionFactor):
    filteredList = []
    for element in someList:
        if len(element) >= exclusionFactor:
            filteredList.append(element)

    return filteredList

# this is taken from https://www.geeksforgeeks.org/dna-protein-python-3/
# this will translate  DNA sequences into their corresponding amino acid based on the codon
# Note: the transcription step is ignored, since the sequence is already 5' to 3' and mRNA would just have the same 5' to 3' sequence as the coding strand (when reading the complement from the template) with T's subbed with U's
# Note 2: updated so that stop codons are mapped to empty string so nothing is appended when the stop codon is reached
def translate(seq):
       
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
        'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein

# Input: a dictionary someDict with values of integer (keys can be anything), an integer n, optionally listVal
# default n = 1
# default listVal = 0
# listVal is optional, and tells the function what to do when the value for the key is a list
# Output: a dictionary outputDict with n keys/values whose keys and values are are those keys/value pairs in someDict with the greatest values
# if value are equal, then the first key/value in the dict is the first one that is returned. Increase n to capture more values which may duplicate.

def GetTopORBotNDict(someDict,n=1,listVal=0,bot=False):
    copySomeDict = someDict.copy()
    outputDict = {} # see output of func defn
    
    # keep running this loop while the values of outputDict are being populated (see Output in func defn)
    while len(outputDict) < n:
        topKey = ""
        topVal = 0
        
        # range over all items in the copy dictionary
        for key in copySomeDict:
            # check if current item has a value greater than someKey, if listVal is non-zero then this means the function call intended for the dict value to be a list
            curVal = copySomeDict[key]
            if listVal != 0:
                curVal = curVal[listVal] # get the specified list element from the list defined by the value associated to the key
            
            
            if curVal > topVal:
                topKey = key
                topVal = curVal

            if bot == True:
                if curVal < topVal:
                    topKey = key
                    topVal = curVal
        
        # remove the most recently found topKey/topValue pair from the input dictionary        
        copySomeDict.pop(topKey) 
        # collect topKey/topValue pair into output dictionary
        outputDict[topKey] = outputDict.get(topKey,topVal) 
        
    return outputDict


  
senseDict, antisenseDict = AnalyzeBacterialGenome(genome_text)
print("\nGenome analysis completed! Getting the smallest proteins now.")


from Bio import Entrez # connects us to https://www.ncbi.nlm.nih.gov/
from Bio.Blast import NCBIWWW # allows us to run blast
Entrez.email = "egaskin@andrew.cmu.edu"

topNum = 100
botSense = GetTopORBotNDict(senseDict,topNum,1,True)
botSense = list(botSense.keys())
botAntisense = GetTopORBotNDict(antisenseDict,topNum,1,True)
botAntisense = list(botAntisense.keys())

print("\n{0} smallest proteins found! Running BLASTp now.".format(topNum))

# print the things that will be blasted 
""" 
print("Number of potential proteins:",len(senseDict))
    for i in range(-5,0):
        desiredKey_sense = botSense[i] # get the key from the min list
        sequenceToBlast_sense = senseDict[desiredKey_sense][0] # use the key to get the sequence from the dict from analysis output
        print("\nSense protein",desiredKey_sense,"has length",len(sequenceToBlast_sense),"amino acids")
        print("sequence",sequenceToBlast_sense)

        desiredKey_antisense = botAntisense[i]
        sequenceToBlast_antisense = antisenseDict[desiredKey_antisense][0]
        print("\nAntisense protein",desiredKey_antisense,"has length",len(sequenceToBlast_antisense),"amino acids")
        print("sequence",sequenceToBlast_antisense)
 """

 # open the save file in append mode. so we can  compile all the results
save_file_sense = open("my_SenseBlast.xml","a")
save_file_antisense = open("my_AntisenseBlast.xml","a")

# run BLASTp on the last 5 of the final lists from the sense and antisense dictonaries
# https://biopython.org/docs/1.75/api/Bio.Blast.NCBIWWW.html
par1 = "blastp" # parameter 1: program
par2 = "np"
par3 = 10 # parameter 3: alignments
par4 = 10 # parameter 4: hitlist_size
par5 = True
par6 = "XML"
# ranging from -5 to 0 gives -5, -4... -1. we can use this to reference the last elements of the min list, which will be larger than ones toward the beginning
for i in range(-5,0):
    desiredKey_sense = botSense[i] # get the key from the min list
    sequenceToBlast_sense = senseDict[desiredKey_sense][0] # use the key to get the sequence from the dict from analysis output
    print("\nrunning BLASTp for sense",desiredKey_sense,"which has length",len(sequenceToBlast_sense),"amino acids")
    print("sequence",sequenceToBlast_sense)
    result_handle_sense = NCBIWWW.qblast(program = par1, sequence = sequenceToBlast_sense, database=par2,alignments=par3,hitlist_size=par4,short_query=par5,format_type=par6)
    additionalInfoSense = "\nSense Protein:" + desiredKey_sense + "\n"
    save_file_sense.write(additionalInfoSense)
    save_file_sense.write(result_handle_sense.read())

    desiredKey_antisense = botAntisense[i]
    sequenceToBlast_antisense = antisenseDict[desiredKey_antisense][0]
    print("\nrunning BLASTp for antisense",desiredKey_antisense,"which has length",len(sequenceToBlast_antisense),"amino acids")
    print("sequence",sequenceToBlast_antisense)
    # result_handle_antisense = NCBIWWW.qblast("blastp","nr", sequenceToBlast_antisense, hitlist_size = 10, short_query = True)
    result_handle_antisense = NCBIWWW.qblast(program = par1, sequence = sequenceToBlast_antisense, database=par2,alignments=par3,hitlist_size=par4,short_query=par5,format_type=par6)
    additionalInfoAnti = "\nAntisense Protein:" + desiredKey_antisense + "\n"
    save_file_antisense.write(additionalInfoAnti)
    save_file_antisense.write(result_handle_antisense.read())

save_file_sense.close()
save_file_antisense.close()
result_handle_sense.close()
result_handle_antisense.close()

