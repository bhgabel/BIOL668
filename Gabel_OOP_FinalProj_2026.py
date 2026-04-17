## RENAME this file YourLastName_OOP_FinalProject_2026.py

##Assignment: Add to the constructor and methods of a parent class and child classes
##            which inherit the base class properties. NOTE: You are not allowed
##            to import any specialized libraries for this project (e.g., no Biopython)
##            The idea is for you to write these methods from scratch.

## Begin with the parent Seq class and the child DNA class we created in lecture below.
## 

### Seq Class
#
#  Constructor:
#  (1) Use the string functions upper and strip to clean up self.sequence.
#  (2) Add a variable self.kmers to the constructor and make it equal to an empty list.

#  Methods:
#  (1) Add a method called make_kmers that makes overlapping kmers of a given length from self.sequence
#      appends these to self.kmers. Default kmer parameter=3.
#  (2) Add a method called fasta that returns a fasta formatted string like this:
#      >species gene
#      AGATTGATAGATAGATAT


### DNA Class: INHERITS Seq class
#   
#  Constructor:
#  Use re.sub to change any non nucleotide characters in self.sequence into an 'N'.
#      re.sub('[^ATGCU]','N',sequence) will change any character that is not a
#      capital A, T, G, C or U into an N. (Seq already uppercases and strips.)

#  Methods:
#  (1) Add a method called print_info that is like print_record, but adds geneid and an
#      empty space to the beginning of the string.
#  (2) Add a method called reverse_complement that returns the reverse complement of
#      self.sequence
#  (3) Add a method called six_frames that returns all 6 frames of self.sequence
#      This include the 3 forward frames, and the 3 reverse complement frames

### RNA Class:  INHERITS DNA class
#  
#  Construtor:
#  Use the super() function (see DNA Class example).
#  (1) Automatically change all Ts to Us in self.sequence. 
#  (2) Add self.codons equals to an empty list

#  Methods:
#  (1) Add make_codons which breaks the self.sequence into 3 letter codons
#      and appends these codons to self.codons unless they are less than 3 letters long.
#  (2) Add translate which uses the Global Variable standard_code below to
#      translate the codons in self.codons and returns a protein sequence.

### Protein Class: INHERITS Seq class
#
#  Construtor:
#  Use the super() function (see DNA Class example).
#  Use re.sub to change any non LETTER characters in self.sequence into an 'X'.

#  Methods:
#  The next 2 methods use a kyte_doolittle and the aa_mol_weights dictionaries.
#  (2) Add total_hydro, which return the sum of the total hydrophobicity of a self.sequence
#  (3) Add mol_weight, which returns the total molecular weight of the protein
#      sequence assigned to the protein object. 



# Methods added by Francis:
# 1. Override on len() to return length of the sequence attribute of Seq class
# 2. A charge() method to find the total charge of a protein based on AA charges at physiological pH


import re

standard_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

kyte_doolittle={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,
                'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'X':0,'Y':-1.3}

aa_mol_weights={'A':89.09,'C':121.15,'D':133.1,'E':147.13,'F':165.19,
                'G':75.07,'H':155.16,'I':131.17,'K':146.19,'L':131.17,
                'M':149.21,'N':132.12,'P':115.13,'Q':146.15,'R':174.2,
                'S':105.09,'T':119.12,'V':117.15,'W':204.23,'X':0,'Y':181.19}


class Seq:

    def __init__(self, sequence, gene, species):
        self.sequence = sequence.strip().upper()
        self.gene = gene
        self.species = species
        self.kmers = []

    def __str__(self):
        return self.sequence

    def print_record(self):
        print(self.species + " " + self.gene + ": " + self.sequence)

    def make_kmers(self, k=3):
        self.kmers=[]
        for i in range(0, len(self.sequence)-k+1):
            kmer = self.sequence[i:i+k]
            if len(kmer) != k: pass

            if kmer not in self.kmers:
                self.kmers.append(kmer)

    def fasta(self):
        return ">" + self.species + " " + self.gene + "\n" + self.sequence
    
    #Count occurances of provided motif
    #@param motif [optional]: pattern to find in sequence
    #   if not provided, will use self.kmers list
    def kmerCount(self, motif=""):
        """Count occruances of provided motif or of object's kmers

        #create Seq object
        >>> seq = Seq("ACGTACGT", "gene", "species")

        #check with empty self.kmers list
        >>> seq.kmerCount()
        0

        >>> seq.make_kmers(k=3)
        >>> seq.kmerCount("ACGT")
        2
        >>> seq.kmerCount()
        {'ACG': 2, 'CGT': 2, 'GTA': 1, 'TAC': 1}

        #check with motif not in sequence
        >>> seq.kmerCount("XYZ")
        0
        """
        if motif == "":
            #check that self.kemers is non-empty
            if len(self.kmers) == 0:
                return 0
            
            counts = {kmer: 0 for kmer in self.kmers}
            for kmer in counts.keys():
                for i in range(0, len(self.sequence)):
                    if kmer == self.sequence[i:i+len(kmer)]:
                        counts[kmer] += 1
            return counts

        else:
            motif = motif.strip().upper()
            count = 0
            for i in range(0, len(self.sequence)):
                if motif == self.sequence[i:i+len(motif)]:
                    count += 1
            return count

    #Override equal and not equal operators
    def __eq__(self, other):
        """
        >>> seq1 = Seq("ACGTACGT", "gene", "species")
        >>> seq2 = Seq("ACGTACGT", "x", "y")
        >>> seq3 = Seq("ACGT", "gene", "species")
        >>> seq1 == seq2
        True
        >>> seq1 == seq3
        False

        #compare to non- Seq object -> false
        >>> seq1 == "ACGTACGT"
        False
        """

        if isinstance(other, Seq):
            return self.sequence == other.sequence
        else:
            return False

    def __ne__(self, other):
        return not self == other
    
    # FL: Override len to return length of the sequence not the class
    def __len__(self):
        return len(self.sequence)
    
class DNA(Seq):

    def __init__(self, sequence, gene, species, geneid, **kwargs):
        super().__init__(sequence, gene, species)
        self.sequence = re.sub('[^ATGCU]','N', self.sequence)
        self.geneid = geneid
 
    def analysis(self):
        """ Returns the GC count of the DNA sequence

        >>> seq = DNA("ACGTACGT", "gene", "species", "geneid")
        >>> seq.analysis()
        4
        >>> DNA("GGGGCC", "gene", "species", "geneid").analysis()
        6
        >>> DNA("XYZ", "gene", "species", "geneid").analysis()
        0
        """
        gc = len(re.findall('G',self.sequence) + re.findall('C',self.sequence))
        return gc

    def print_info(self):
        print(self.species + " " + self.gene + " " + self.geneid + ": " + self.sequence)

    def reverse_complement(self):
        """Returns the reverse complement of the DNA sequence

        >>> DNA("ACGTACGT", "gene", "species", "geneid").reverse_complement()
        'ACGTACGT'

        >>> DNA("ACGTAAA", "gene", "species", "geneid").reverse_complement()
        'TTTACGT'

        >>> DNA("ACGTXYZ", "gene", "species", "geneid").reverse_complement()
        'NNNACGT'
        """
        reverse = self.sequence[::-1]
        complement = ""
        for base in reverse:
            if base == 'A':
                complement += 'T'
            elif base == 'C':
                complement += 'G'
            elif base == 'G':
                complement += 'C'
            elif base == 'T':
                complement += 'A'
            else:
                complement += 'N'
        return complement

    #returns 6 lists, indexes 0,1,2 are forward; 3,4,5 are reverse
    def six_frames(self):
        frames = []

        #forward
        frames.append(self.sequence[0:])
        frames.append(self.sequence[1:])
        frames.append(self.sequence[2:])

        reverse = self.reverse_complement()
        frames.append(reverse[0:])
        frames.append(reverse[1:])
        frames.append(reverse[2:])

        return frames



class RNA(DNA):

    def __init__(self, sequence, gene, species, geneid, **kwargs):
        super().__init__(sequence, gene, species, geneid)
        self.sequence = re.sub('T', 'U', self.sequence)
        self.codons = []
        
    def make_codons(self):
        self.codons = []
        for i in range(0, len(self.sequence), 3):
            codon = self.sequence[i:i+3]
            if len(codon) < 3: pass
            else:
                self.codons.append(codon)
 
    def translate(self):
        protein = ""
        for codon in self.codons:
            protein += standard_code.get(codon, 'X')
        return protein

class Protein(Seq):

    def __init__(self, sequence, gene, species, *args, **kwargs,):
        super().__init__(sequence, gene, species)
        self.sequence = re.sub('[^A-Z]', 'X', self.sequence)

    def total_hydro(self):
        hydro = 0
        for aa in self.sequence:
            hydro += kyte_doolittle[aa]
        return hydro

    def mol_weight(self):
        weight = 0
        for aa in self.sequence:
            weight += aa_mol_weights[aa]
        return weight
    
    # FL: estimate total charge of the protein by subtracting all - charges (acid) from all + charges (base)
    # based on amino acid charge at physiological pH
    def charge(self):
        base = len(re.findall('[RHL]', self.sequence))
        acid = len(re.findall('[DE]', self.sequence))
        total_charge = base-acid
        return total_charge




#Testing section
if __name__ == "__main__":

    import doctest
    doctest.testmod()

