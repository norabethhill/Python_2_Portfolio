# Python_2_Portfolio
This is the portfolio for the Python code I learned this quarter! 

## Sequence Objects - Parts 1, 2, 3, & 4

In this analysis, we created biologically relevant strings, utilized manipulations for seq objects, and added sequences together. We transcribed and translated proteins and mRNA.



```python
from Bio.Seq import Seq
```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index,letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# We can also print the length of each sequence
print(len(my_seq))
```

    5



```python
# We can access different positions with seq objects
print(my_seq[0])
```

    G



```python
# We are looking at the position of the last object
print(my_seq[4])
```

    G



```python
# We get an error if we try to access position 5 because our sequence only goes up to 4
print(my_seq[2])
```

    T



```python
Seq("AAAA").count("AA")
```




    2




```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGGAAAAATCGAGT")
```


```python
# This tells us how long my sequence is
len(my_seq)
```




    36




```python
# We can count the specific number of letters
my_seq.count("G")
```




    11




```python
# This can be important for designing primers
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    44.44444444444444




```python
from Bio.SeqUtils import gc_fraction 
```


```python
my_seq = Seq("GATCGATGGGCCTATATAAGGAAAGAAATGCAT")
```


```python
gc_fraction(my_seq)
```




    0.3939393939393939




```python
# We can slice sequences into multiple parts
my_seq[4:12]
```




    Seq('GATGGGCC')




```python
my_seq[0::3]
```




    Seq('GCTGTAAAGAC')




```python
my_seq[1::3]
```




    Seq('AGGCATGAATA')




```python
# This is important for selecting our start position within the codons
my_seq[2:3]
```




    Seq('T')




```python
# We can use negatives to start from the opposite end
my_seq[::-1]
```




    Seq('TACGTAAAGAAAGGAATATATCCGGGTAGCTAG')




```python
str(my_seq)
```




    'GATCGATGGGCCTATATAAGGAAAGAAATGCAT'




```python
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    >Name
    GATCGATGGGCCTATATAAGGAAAGAAATGCAT
    



```python
# The name tells us what our fasta string is, like a label
```


```python
# We can add 2 sequences together 
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
```


```python
spacer = Seq("N" *10) 
```


```python
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
# Sometimes we have case sensitive sequences 
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
dna_seq.lower()
```




    Seq('acgtacgt')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
# We can search for certain codons 
"gtac" in dna_seq
```




    False




```python
"GTAC" in dna_seq
```




    False




```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GATCGATGGGGCCTATATATAGGATCGAAATCG")
```


```python
# We can base pair the sequence 
my_seq.complement()
```




    Seq('CTAGCTACCCCGGATATATATCCTAGCTTTAGC')




```python
my_seq.reverse_complement()
```




    Seq('CGATTTCGATCCTATATATAGGCCCCATCGATC')




```python
protein_seq = Seq("EVRNAK")
protein_seq.complement() 
```




    Seq('EBYNTM')




```python
# These are ambiguity codes
```


```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')




```python
# Now we are translating it into a protein
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
messenger_rna.translate()
```




    Seq('MAIVMGR*KGAR*')




```python
# That is our translated protein
# The asterisk represent stop codons 
```


```python
# We can specify codons
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('MAIVMGRWKGAR*')




```python
coding_dna.translate(table = 2)
```




    Seq('MAIVMGRWKGAR*')




```python
coding_dna.translate(table =2, to_stop=True)
```




    Seq('MAIVMGRWKGAR')




```python
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('MAIVMGRWKGAR!')




```python
# This is for irregular stop codons, like ones found in bacteria
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HR*')




```python
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
gene.translate(table = "Bacterial", cds = True)
```




    Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
# Vertebrate Mitochondrial is table 2
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
# These are the stop codons for the mito table
mito_table.stop_codons 
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
# These are the start codons for the mito table
mito_table.start_codons 
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
# We used the mito table because this is mitochondrial DNA and these are our start codons
```


```python
# Now lets talk about sequence comparison 
seq = Seq("ACGT")
```


```python
"ACGT" == seq1
```




    True




```python
seq1 == "ACGT"
```




    True




```python
# This can be used if know the length of a sequence but not the contents
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
# This is useful for comparing gene sizes 
```


```python
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
seq[1000:1020]
```




    Seq(None, length=20)




```python
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length =10)
```


```python
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


```python
# We can create a mutable sequence 
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# We are replacing the letter in the 5th position with a C
mutable_seq[5] = "C"
```


```python
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# That just removes the first T
```


```python
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
# This gives us the reverse of our sequence 
```


```python
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GCTTGTTAGGGTCGTTGGAAGGGTGGAATTTCGTCGTCGTGGTTA"
```


```python
reverse_complement(my_string)
```




    'TAACCACGACGACGAAATTCCACCCTTCCAACGACCCTAACAAGC'




```python
back_transcribe(my_string)
```




    'GCTTGTTAGGGTCGTTGGAAGGGTGGAATTTCGTCGTCGTGGTTA'




```python
# We can translate the string to a protein
translate(my_string)
```




    'AC*GRWKGGISSSWL'


    
## Sequence Annotations - Parts 1, 2, 3, & 4 

In this analysis, we worked with Sequence Annotations over four parts. 


```python
from Bio.SeqRecord import SeqRecord
```


```python
from Bio.Seq import Seq
```


```python
simple_seq = Seq("GATC")
```


```python
simple_seq_r = SeqRecord(simple_seq)
```


```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
# We have created a record but have not provided any information yet
```


```python
# We can also pass the ID to it
simple_seq_r.id = "AC12345"
```


```python
simple_seq_r.description = "Made up sequence for the VDB Computational Biology Class"
```


```python
print(simple_seq_r.description)
```

    Made up sequence for the VDB Computational Biology Class



```python
simple_seq_r.seq
```




    Seq('GATC')




```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
simple_seq_r.annotations["evidence"] = "None. This is just an example"
```


```python
print(simple_seq_r.annotations["evidence"])
```

    None. This is just an example



```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
# We can also do per letter annotations 
simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
```


```python
print(simple_seq_r.letter_annotations)
```

    {'phred_quality': [40, 40, 38, 30]}



```python
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.fna.txt", "fasta")
```


```python
# This gives us the information for the fasta file
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
# We can look at individual parts 
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.dbxrefs
```




    []




```python
record.annotations
```




    {}




```python
record.features
```




    []




```python
# We are now working with a larger file
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'NC_005816.1'




```python
record.name
```




    'NC_005816'




```python
record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.letter_annotations 
```




    {}




```python
len(record.annotations)
```




    13




```python
record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
# This source is where the DNA is found
```


```python
# This will give us the project number
record.dbxrefs
```




    ['Project:58037']




```python
len(record.features)
```




    41




```python
from Bio import SeqFeature
```


```python
start_pos = SeqFeature.AfterPosition(5)
# The feature starts after position 5 in the sequence
```


```python
end_pos = SeqFeature.BetweenPosition(9, left=8, right=9)
```


```python
# Now we can create a location
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
print(my_location)
```

    [>5:(8^9)]



```python
my_location.start
```




    AfterPosition(5)




```python
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
int(my_location.end)
# This will give you just the number of the location
```




    9




```python
int(my_location.start)
```




    5




```python
exact_location = SeqFeature.SimpleLocation(5,9)
```


```python
print(exact_location)
```

    [5:9]



```python
exact_location.start
```




    ExactPosition(5)




```python
# Now we are going into pt 3
```


```python
from Bio.SeqRecord import SeqRecord 
```


```python
# Now we are going to create a seq record 
```


```python
record = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
...         "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
...         "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
...         "SSAC"),
                   id="gi|14150838|gb|AAK54648.1|AF376133",
...       description="chalcone synthase [Cucumis sativus]",
... )
```


```python
# Now we are going to print our record as a fasta file
```


```python
print(record.format("fasta"))
```

    >gi|14150838|gb|AAK54648.1|AF376133 chalcone synthase [Cucumis sativus]
    MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD
    GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK
    NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM
    SSAC
    



```python
# This gives us a different format (not fasta)
```


```python
print(record)
```

    ID: gi|14150838|gb|AAK54648.1|AF376133
    Name: <unknown name>
    Description: chalcone synthase [Cucumis sativus]
    Number of features: 0
    Seq('MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVG...SAC')



```python
from Bio import SeqIO
```


```python
# now we are going to read in a record 
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# This gives the length of the record 
len(record)
```




    9609




```python
# We can look at the length of the features
len(record.features)
```




    41




```python
# There are 41 features within this record 
```


```python
# This prints the 20th feature
print(record.features[20])
```

    type: gene
    location: [4342:4780](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(record.features[21])
```

    type: CDS
    location: [4342:4780](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# this feature is a CDS
```


```python
# We can subdivide the record 
sub_record = record[4300:4800]
```


```python
len(sub_record)
```




    500




```python
# This tells us how many features we have within the sub record
len(sub_record.features)
```




    2




```python
# there are 2 features within the subrecord 
```


```python
sub_record.features[0]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='gene', qualifiers=...)




```python
sub_record.features[1]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='CDS', qualifiers=...)




```python
print(sub_record.features[0])
```

    type: gene
    location: [42:480](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
# the location changes bc we are looking at the subrecord 
```


```python
print(sub_record.features[1])
```

    type: CDS
    location: [42:480](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# we can also look at the annotations 
sub_record.annotations 
```




    {'molecule_type': 'DNA'}




```python
# We are looking at a DNA molecule 
```


```python
sub_record.dbxrefs
```




    []




```python
sub_record.annotations["topology"] = "linear"
```


```python
# We can add to our annotations for the record 
```


```python
sub_record.annotations 
```




    {'molecule_type': 'DNA', 'topology': 'linear'}




```python
sub_record.id 
```




    'NC_005816.1'




```python
# this gives the record id
```


```python
sub_record.name 
```




    'NC_005816'




```python
# the name and ID are going to be similar in this case 
```


```python
sub_record.description 
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
sub_record.description = 'Yersinia pestis Microtus str. 91001 plasmid pPCP1, partial sequence'
```


```python
sub_record.description 
```




    'Yersinia pestis Microtus str. 91001 plasmid pPCP1, partial sequence'




```python
# Now it says partial sequence instead of complete sequence 
```


```python
print(sub_record.format("genbank")[:200] + "...")
```

    LOCUS       NC_005816                500 bp    DNA     linear   UNK 01-JAN-1980
    DEFINITION  Yersinia pestis Microtus str. 91001 plasmid pPCP1, partial sequence.
    ACCESSION   NC_005816
    VERSION     NC_00...



```python
# Now onto sequence annotation pt 4 
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
len(record)
```




    9609




```python
len(record.features)
```




    41




```python
record.dbxrefs
```




    ['Project:58037']




```python
# We can also look up the record annotations and keys 
record.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
# this lets us see all the different things we can oull up
```


```python
shifted = record[2000:] + record[:2000]
```


```python
shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
len(shifted)
```




    9609




```python
# the length is the same bc we didnt change anything about the information 
```


```python
len(shifted.features)
```




    40




```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
shifted.dbxrefs
```




    []




```python
shifted.dbxrefs = record.dbxrefs[:]
```


```python
shifted.dbxrefs
```




    ['Project:58037']




```python
shifted.annotations = record.annotations.copy()
```


```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
# now we can look at reverse complement 
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13



```python
# We printed the first value as a string
# The string is our record id
```


```python
# the length of the record, features, dbxrefs, and annotations are all printed as integers 
```


```python
rc = record.reverse_complement(id ="Testing")
```


```python
rc
```




    SeqRecord(seq=Seq('CAGGGGTCGGGGTACGCATTCCCTCATGCGTCAATATTATCTGGCATTGCGATG...ACA'), id='Testing', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
# now we have the same thing just backwards 
```


```python
print("%s %i %i %i %i" % (rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
```

    Testing 9609 41 0 0


## Sequence I/O - Parts 1, 2, & 3

In this analysis, we worked on input/output with Seq files throughout three parts.


```python
from Bio import SeqIO
```


```python
# We can manipulate these long databases to analyze certain parts or entries
```


```python
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    gi|2765658|emb|Z78533.1|CIZ78533
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    gi|2765657|emb|Z78532.1|CCZ78532
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    gi|2765656|emb|Z78531.1|CFZ78531
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    gi|2765655|emb|Z78530.1|CMZ78530
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    gi|2765654|emb|Z78529.1|CLZ78529
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    gi|2765652|emb|Z78527.1|CYZ78527
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    gi|2765651|emb|Z78526.1|CGZ78526
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    gi|2765650|emb|Z78525.1|CAZ78525
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    gi|2765649|emb|Z78524.1|CFZ78524
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    gi|2765648|emb|Z78523.1|CHZ78523
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    gi|2765647|emb|Z78522.1|CMZ78522
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    gi|2765646|emb|Z78521.1|CCZ78521
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    gi|2765645|emb|Z78520.1|CSZ78520
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    gi|2765644|emb|Z78519.1|CPZ78519
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    gi|2765643|emb|Z78518.1|CRZ78518
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    gi|2765642|emb|Z78517.1|CFZ78517
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    gi|2765641|emb|Z78516.1|CPZ78516
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    gi|2765640|emb|Z78515.1|MXZ78515
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    gi|2765639|emb|Z78514.1|PSZ78514
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    gi|2765638|emb|Z78513.1|PBZ78513
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    gi|2765637|emb|Z78512.1|PWZ78512
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    gi|2765636|emb|Z78511.1|PEZ78511
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    gi|2765635|emb|Z78510.1|PCZ78510
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    gi|2765634|emb|Z78509.1|PPZ78509
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    gi|2765633|emb|Z78508.1|PLZ78508
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    gi|2765632|emb|Z78507.1|PLZ78507
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    gi|2765631|emb|Z78506.1|PLZ78506
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    gi|2765630|emb|Z78505.1|PSZ78505
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    gi|2765629|emb|Z78504.1|PKZ78504
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    gi|2765628|emb|Z78503.1|PCZ78503
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    gi|2765627|emb|Z78502.1|PBZ78502
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    gi|2765626|emb|Z78501.1|PCZ78501
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    gi|2765625|emb|Z78500.1|PWZ78500
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    gi|2765624|emb|Z78499.1|PMZ78499
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    gi|2765623|emb|Z78498.1|PMZ78498
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    gi|2765622|emb|Z78497.1|PDZ78497
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    gi|2765621|emb|Z78496.1|PAZ78496
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    gi|2765620|emb|Z78495.1|PEZ78495
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    gi|2765619|emb|Z78494.1|PNZ78494
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    gi|2765618|emb|Z78493.1|PGZ78493
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    gi|2765617|emb|Z78492.1|PBZ78492
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    gi|2765616|emb|Z78491.1|PCZ78491
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    gi|2765615|emb|Z78490.1|PFZ78490
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    gi|2765614|emb|Z78489.1|PDZ78489
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    gi|2765613|emb|Z78488.1|PTZ78488
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    gi|2765612|emb|Z78487.1|PHZ78487
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765611|emb|Z78486.1|PBZ78486
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    gi|2765610|emb|Z78485.1|PHZ78485
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    gi|2765609|emb|Z78484.1|PCZ78484
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    gi|2765608|emb|Z78483.1|PVZ78483
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    gi|2765607|emb|Z78482.1|PEZ78482
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    gi|2765606|emb|Z78481.1|PIZ78481
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    gi|2765605|emb|Z78480.1|PGZ78480
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    gi|2765604|emb|Z78479.1|PPZ78479
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    gi|2765603|emb|Z78478.1|PVZ78478
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    gi|2765602|emb|Z78477.1|PVZ78477
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    gi|2765601|emb|Z78476.1|PGZ78476
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    gi|2765600|emb|Z78475.1|PSZ78475
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    gi|2765599|emb|Z78474.1|PKZ78474
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    gi|2765598|emb|Z78473.1|PSZ78473
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    gi|2765597|emb|Z78472.1|PLZ78472
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    gi|2765596|emb|Z78471.1|PDZ78471
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765595|emb|Z78470.1|PPZ78470
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    gi|2765594|emb|Z78469.1|PHZ78469
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    gi|2765593|emb|Z78468.1|PAZ78468
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    gi|2765592|emb|Z78467.1|PSZ78467
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    gi|2765591|emb|Z78466.1|PPZ78466
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    gi|2765590|emb|Z78465.1|PRZ78465
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    gi|2765589|emb|Z78464.1|PGZ78464
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    gi|2765588|emb|Z78463.1|PGZ78463
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    gi|2765587|emb|Z78462.1|PSZ78462
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    gi|2765586|emb|Z78461.1|PWZ78461
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765585|emb|Z78460.1|PCZ78460
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    gi|2765584|emb|Z78459.1|PDZ78459
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    gi|2765583|emb|Z78458.1|PHZ78458
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    gi|2765582|emb|Z78457.1|PCZ78457
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    gi|2765581|emb|Z78456.1|PTZ78456
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765580|emb|Z78455.1|PJZ78455
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    gi|2765579|emb|Z78454.1|PFZ78454
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    gi|2765578|emb|Z78453.1|PSZ78453
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    gi|2765577|emb|Z78452.1|PBZ78452
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    gi|2765576|emb|Z78451.1|PHZ78451
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    gi|2765575|emb|Z78450.1|PPZ78450
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    gi|2765574|emb|Z78449.1|PMZ78449
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    gi|2765573|emb|Z78448.1|PAZ78448
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    gi|2765572|emb|Z78447.1|PVZ78447
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    gi|2765571|emb|Z78446.1|PAZ78446
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    gi|2765570|emb|Z78445.1|PUZ78445
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    gi|2765569|emb|Z78444.1|PAZ78444
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    gi|2765568|emb|Z78443.1|PLZ78443
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    gi|2765567|emb|Z78442.1|PBZ78442
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    gi|2765566|emb|Z78441.1|PSZ78441
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    gi|2765565|emb|Z78440.1|PPZ78440
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    gi|2765564|emb|Z78439.1|PBZ78439
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
# Now doing the same thing but for the other file
```


```python
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    Z78532.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    Z78531.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    Z78530.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    Z78529.1
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    Z78527.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    Z78526.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    Z78525.1
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    Z78524.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    Z78523.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    Z78522.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    Z78521.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    Z78520.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    Z78519.1
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    Z78518.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    Z78517.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    Z78516.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    Z78515.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    Z78514.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    Z78513.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    Z78512.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    Z78511.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    Z78510.1
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    Z78509.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    Z78508.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    Z78507.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    Z78506.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    Z78505.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    Z78504.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    Z78503.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    Z78502.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    Z78501.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    Z78500.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    Z78499.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    Z78498.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    Z78497.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    Z78496.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    Z78495.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    Z78494.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    Z78493.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    Z78492.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    Z78491.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    Z78490.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    Z78489.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    Z78488.1
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    Z78487.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78486.1
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    Z78485.1
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    Z78484.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    Z78483.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    Z78482.1
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    Z78481.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    Z78480.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    Z78479.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    Z78478.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    Z78477.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    Z78476.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    Z78475.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    Z78474.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    Z78473.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    Z78472.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    Z78471.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78470.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    Z78469.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    Z78468.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    Z78467.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    Z78466.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    Z78465.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    Z78464.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    Z78463.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    Z78462.1
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    Z78461.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78460.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    Z78459.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    Z78458.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    Z78457.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    Z78456.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78455.1
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    Z78454.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    Z78453.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    Z78452.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    Z78451.1
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    Z78450.1
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    Z78449.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    Z78448.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    Z78447.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    Z78446.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    Z78445.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    Z78444.1
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    Z78443.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    Z78442.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    Z78441.1
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    Z78440.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
# Now we are extracting things as a list
identifiers = [seq_record.id for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")]
```


```python
identifiers 

```




    ['Z78533.1',
     'Z78532.1',
     'Z78531.1',
     'Z78530.1',
     'Z78529.1',
     'Z78527.1',
     'Z78526.1',
     'Z78525.1',
     'Z78524.1',
     'Z78523.1',
     'Z78522.1',
     'Z78521.1',
     'Z78520.1',
     'Z78519.1',
     'Z78518.1',
     'Z78517.1',
     'Z78516.1',
     'Z78515.1',
     'Z78514.1',
     'Z78513.1',
     'Z78512.1',
     'Z78511.1',
     'Z78510.1',
     'Z78509.1',
     'Z78508.1',
     'Z78507.1',
     'Z78506.1',
     'Z78505.1',
     'Z78504.1',
     'Z78503.1',
     'Z78502.1',
     'Z78501.1',
     'Z78500.1',
     'Z78499.1',
     'Z78498.1',
     'Z78497.1',
     'Z78496.1',
     'Z78495.1',
     'Z78494.1',
     'Z78493.1',
     'Z78492.1',
     'Z78491.1',
     'Z78490.1',
     'Z78489.1',
     'Z78488.1',
     'Z78487.1',
     'Z78486.1',
     'Z78485.1',
     'Z78484.1',
     'Z78483.1',
     'Z78482.1',
     'Z78481.1',
     'Z78480.1',
     'Z78479.1',
     'Z78478.1',
     'Z78477.1',
     'Z78476.1',
     'Z78475.1',
     'Z78474.1',
     'Z78473.1',
     'Z78472.1',
     'Z78471.1',
     'Z78470.1',
     'Z78469.1',
     'Z78468.1',
     'Z78467.1',
     'Z78466.1',
     'Z78465.1',
     'Z78464.1',
     'Z78463.1',
     'Z78462.1',
     'Z78461.1',
     'Z78460.1',
     'Z78459.1',
     'Z78458.1',
     'Z78457.1',
     'Z78456.1',
     'Z78455.1',
     'Z78454.1',
     'Z78453.1',
     'Z78452.1',
     'Z78451.1',
     'Z78450.1',
     'Z78449.1',
     'Z78448.1',
     'Z78447.1',
     'Z78446.1',
     'Z78445.1',
     'Z78444.1',
     'Z78443.1',
     'Z78442.1',
     'Z78441.1',
     'Z78440.1',
     'Z78439.1']




```python
# Now all the identifiers are pulled out 
```


```python
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
print(first_record.id)
```

    gi|2765658|emb|Z78533.1|CIZ78533



```python
print(first_record.description)
```

    gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
second_record = next(record_iterator)
```


```python
print(second_record.id)
```

    gi|2765657|emb|Z78532.1|CCZ78532



```python
print(second_record.description)
```

    gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
# This is how to access the records in any order 
records = list(SeqIO.parse("ls_orchid.gbk.txt", "genbank"))
```


```python
print("Found %i records" % len(records))
```

    Found 94 records



```python
print("The last record")
last_record = records [-1]
print(last_record.id)
print(repr(last_record.seq))
print(len(last_record))
```

    The last record
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
print("The first record")
first_record = records[-3]
print(first_record.id)
print(repr(first_record.seq))
print(len(first_record))
```

    The first record
    Z78441.1
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703



```python
from Bio import SeqIO
```


```python
record_iterator = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
# Remember this record has all our organisms
```


```python
first_record = next(record_iterator)
```


```python
print(first_record)
```

    ID: Z78533.1
    Name: Z78533
    Description: C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Number of features: 5
    /molecule_type=DNA
    /topology=linear
    /data_file_division=PLN
    /date=30-NOV-2006
    /accessions=['Z78533']
    /sequence_version=1
    /gi=2765658
    /keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2']
    /source=Cypripedium irapeanum
    /organism=Cypripedium irapeanum
    /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium']
    /references=[Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')



```python
print(first_record.annotations)
```

    {'molecule_type': 'DNA', 'topology': 'linear', 'data_file_division': 'PLN', 'date': '30-NOV-2006', 'accessions': ['Z78533'], 'sequence_version': 1, 'gi': '2765658', 'keywords': ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'source': 'Cypripedium irapeanum', 'organism': 'Cypripedium irapeanum', 'taxonomy': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], 'references': [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]}



```python
print(first_record.annotations.keys())
```

    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references'])



```python
print(first_record.annotations.values())
```

    dict_values(['DNA', 'linear', 'PLN', '30-NOV-2006', ['Z78533'], 1, '2765658', ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'Cypripedium irapeanum', 'Cypripedium irapeanum', ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]])



```python
print(first_record.annotations["source"])
```

    Cypripedium irapeanum



```python
print(first_record.annotations["organism"])
```

    Cypripedium irapeanum



```python
# The organism is an orchid
# Usually the source will give you the common name
```


```python
all_species = []
```


```python
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    all_species.append(seq_record.annotations["organism"])
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
# This is a list of all 94 species in the document we are analyzing
```


```python
all_species = [
    seq_record.annotations["organism"]
    for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")
]
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
# This gives us a list annotation; same list as the other one
```


```python
all_species = []
```


```python
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    all_species.append(seq_record.description.split()[1])
# This is the other file now
```


```python
print(all_species)
```

    ['C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', 'C.lichiangense', 'C.yatabeanum', 'C.guttatum', 'C.acaule', 'C.formosanum', 'C.himalaicum', 'C.macranthum', 'C.calceolus', 'C.segawai', 'C.pubescens', 'C.reginae', 'C.flavum', 'C.passerinum', 'M.xerophyticum', 'P.schlimii', 'P.besseae', 'P.wallisii', 'P.exstaminodium', 'P.caricinum', 'P.pearcei', 'P.longifolium', 'P.lindenii', 'P.lindleyanum', 'P.sargentianum', 'P.kaiteurum', 'P.czerwiakowianum', 'P.boissierianum', 'P.caudatum', 'P.warszewiczianum', 'P.micranthum', 'P.malipoense', 'P.delenatii', 'P.armeniacum', 'P.emersonii', 'P.niveum', 'P.godefroyae', 'P.bellatulum', 'P.concolor', 'P.fairrieanum', 'P.druryi', 'P.tigrinum', 'P.hirsutissimum', 'P.barbigerum', 'P.henryanum', 'P.charlesworthii', 'P.villosum', 'P.exul', 'P.insigne', 'P.gratrixianum', 'P.primulinum', 'P.victoria', 'P.victoria', 'P.glaucophyllum', 'P.supardii', 'P.kolopakingii', 'P.sanderianum', 'P.lowii', 'P.dianthum', 'P.parishii', 'P.haynaldianum', 'P.adductum', 'P.stonei', 'P.philippinense', 'P.rothschildianum', 'P.glanduliferum', 'P.glanduliferum', 'P.sukhakulii', 'P.wardii', 'P.ciliolare', 'P.dayanum', 'P.hennisianum', 'P.callosum', 'P.tonsum', 'P.javanicum', 'P.fowliei', 'P.schoseri', 'P.bougainvilleanum', 'P.hookerae', 'P.papuanum', 'P.mastersianum', 'P.argus', 'P.venustum', 'P.acmodontum', 'P.urbanianum', 'P.appletonianum', 'P.lawrenceanum', 'P.bullenianum', 'P.superbiens', 'P.purpuratum', 'P.barbatum', 'C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', 'C.lichiangense', 'C.yatabeanum', 'C.guttatum', 'C.acaule', 'C.formosanum', 'C.himalaicum', 'C.macranthum', 'C.calceolus', 'C.segawai', 'C.pubescens', 'C.reginae', 'C.flavum', 'C.passerinum', 'M.xerophyticum', 'P.schlimii', 'P.besseae', 'P.wallisii', 'P.exstaminodium', 'P.caricinum', 'P.pearcei', 'P.longifolium', 'P.lindenii', 'P.lindleyanum', 'P.sargentianum', 'P.kaiteurum', 'P.czerwiakowianum', 'P.boissierianum', 'P.caudatum', 'P.warszewiczianum', 'P.micranthum', 'P.malipoense', 'P.delenatii', 'P.armeniacum', 'P.emersonii', 'P.niveum', 'P.godefroyae', 'P.bellatulum', 'P.concolor', 'P.fairrieanum', 'P.druryi', 'P.tigrinum', 'P.hirsutissimum', 'P.barbigerum', 'P.henryanum', 'P.charlesworthii', 'P.villosum', 'P.exul', 'P.insigne', 'P.gratrixianum', 'P.primulinum', 'P.victoria', 'P.victoria', 'P.glaucophyllum', 'P.supardii', 'P.kolopakingii', 'P.sanderianum', 'P.lowii', 'P.dianthum', 'P.parishii', 'P.haynaldianum', 'P.adductum', 'P.stonei', 'P.philippinense', 'P.rothschildianum', 'P.glanduliferum', 'P.glanduliferum', 'P.sukhakulii', 'P.wardii', 'P.ciliolare', 'P.dayanum', 'P.hennisianum', 'P.callosum', 'P.tonsum', 'P.javanicum', 'P.fowliei', 'P.schoseri', 'P.bougainvilleanum', 'P.hookerae', 'P.papuanum', 'P.mastersianum', 'P.argus', 'P.venustum', 'P.acmodontum', 'P.urbanianum', 'P.appletonianum', 'P.lawrenceanum', 'P.bullenianum', 'P.superbiens', 'P.purpuratum', 'P.barbatum']



```python
# This is abbreveated differently because it is a fasta file and not a genbank file 
```


```python
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
first_record.id
```




    'gi|2765658|emb|Z78533.1|CIZ78533'




```python
# This gives us the first line/ID of the fasta file
```


```python
# We can rename it
first_record.id = "new id"
```


```python
first_record.id
```




    'new id'




```python
first_record
```




    SeqRecord(seq=Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC'), id='new id', name='gi|2765658|emb|Z78533.1|CIZ78533', description='gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', dbxrefs=[])




```python
first_record.description = first_record.id + " " + "Mutations induced randomly"
```


```python
print(first_record.format("fasta"[:200]))
```

    >new id new id Mutations induced randomly
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAA
    CGATCGAGTGAATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGT
    GACCCTGATTTGTTGTTGGGCCGCCTCGGGAGCGTCCATGGCGGGTTTGAACCTCTAGCC
    CGGCGCAGTTTGGGCGCCAAGCCATATGAAAGCATCACCGGCGAATGGCATTGTCTTCCC
    CAAAACCCGGAGCGGCGGCGTGCTGTCGCGTGCCCAATGAATTTTGATGACTCTCGCAAA
    CGGGAATCTTGGCTCTTTGCATCGGATGGAAGGACGCAGCGAAATGCGATAAGTGGTGTG
    AATTGCAAGATCCCGTGAACCATCGAGTCTTTTGAACGCAAGTTGCGCCCGAGGCCATCA
    GGCTAAGGGCACGCCTGCTTGGGCGTCGCGCTTCGTCTCTCTCCTGCCAATGCTTGCCCG
    GCATACAGCCAGGCCGGCGTGGTGCGGATGTGAAAGATTGGCCCCTTGTGCCTAGGTGCG
    GCGGGTCCAAGAGCTGGTGTTTTGATGGCCCGGAACCCGGCAAGAGGTGGACGGATGCTG
    GCAGCAGCTGCCGTGCGAATCCCCCATGTTGTCGTGCTTGTCGGACAGGCAGGAGAACCC
    TTCCGAACCCCAATGGAGGGCGGTTGACCGCCATTCGGATGTGACCCCAGGTCAGGCGGG
    GGCACCCGCTGAGTTTACGC
    



```python
from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
rec1 = SeqRecord(
    Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
        "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
        "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
        "SSAC",
),
    id = "gi|14150838|gb|AAK54658.1|AF376133_1",
    description = "chalcone synthase [Cucumis sativus]",
)
```


```python
rec2 = SeqRecord(
Seq("YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEYMAPSLDARQ"
        "DMVVEIPKLGKEAAVKAIKEWGQ",
   ),
    id = "gi | 13919613 |gb| AAK33142.1",
    description = "chalcone synthatse[Fragaria vesca subsp. bracteata]",
) 
```


```python
rec3 = SeqRecord(
Seq("MVTVEEFRRAQCAEGPATVMAIGTATPSNCVDQSTYPDYYFRITNSEHKVELKEKFKRMC"
        "EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARQDIVVVEVPKLGKEAAQKAIKEWGQP"
        "KSKITHLVFCTTSGVDMPGCDYQLTKLLGLRPSVKRFMMYQQGCFAGGTVLRMAKDLAEN"
        "NKGARVLVVCSEITAVTFRGPNDTHLDSLVGQALFGDAAAVIIGSDPIPEVERPLFELV"
        "SAAQTLLPDSEGAIDGHLREVGLTFHLLKDVPGLIISKNIEKSLVEAFQPLGISDWNSLFW"
        "IAHPGGPAILDQVELKLGLKQEKLKATRKVLSNYGNMSSACVLFILDEMRKAASAKEGLGT"
        "TGEGLEWGVLFGFGPGLTVETVVLHSVAT",
    ),
    id="gi|13925890|gb|AAK49457.1",
    description="chalcone synthase [Nicotiana tabacum]",
)
```


```python
my_records = [rec1, rec2, rec3]
```


```python
from Bio import SeqIO
SeqIO.write(my_records, "my_example.faa", "fasta")
```




    3




```python
records = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
count= SeqIO.write(records, "my_example.fasta", "fasta")
print("Converted %i records" % count)
```

    Converted 94 records



```python
# Now my example is written as a fasta file; we converted them from peptides to a fasta file
```


```python
for record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(record.id)
    print(record.seq.reverse_complement())
```

    Z78533.1
    GCGTAAACTCAGCGGGTGCCCCCGCCTGACCTGGGGTCACATCCGAATGGCGGTCAACCGCCCTCCATTGGGGTTCGGAAGGGTTCTCCTGCCTGTCCGACAAGCACGACAACATGGGGGATTCGCACGGCAGCTGCTGCCAGCATCCGTCCACCTCTTGCCGGGTTCCGGGCCATCAAAACACCAGCTCTTGGACCCGCCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACCACGCCGGCCTGGCTGTATGCCGGGCAAGCATTGGCAGGAGAGAGACGAAGCGCGACGCCCAAGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCAAGATTCCCGTTTGCGAGAGTCATCAAAATTCATTGGGCACGCGACAGCACGCCGCCGCTCCGGGTTTTGGGGAAGACAATGCCATTCGCCGGTGATGCTTTCATATGGCTTGGCGCCCAAACTGCGCCGGGCTAGAGGTTCAAACCCGCCATGGACGCTCCCGAGGCGGCCCAACAACAAATCAGGGTCACCACGGGAGCAATGCCCCCGGTGAGCTGAGTACACCGGTCCTCCGGATTCACTCGATCGTTTATTCCACGGTCTCATCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78532.1
    GCCTCAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTGAATGGAAATCAACTGCCCAATGGTTATTTTAGCTCCATTGGGGTTCAATTAGGTTCTTGTGTAGGTTCGAAAAAATACAACAACATGGGGGATTCAAATAGCAGCCTTATGACTGTTAGCATTCTCCACCTCGTGCCACATTCCTACCCATCAAAGCAACAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCATTCACATCCGTATAATGCCAGCTTAGCGATATGCCAAGCAAGCATTGGTAGGAGAGACGCAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTAGCTGCGTTCTTCATCGATGTTAGAGCCAAGATATCCGTTGCGAGAGTCATAAAAATTCACTGGCATGCTCAGTAGCATACTGCCCCTCTGATTTTTTCTGACAATAATGTCATTCATCAGTGATGCTTTATATGACTTGGCGCAAACTGCACCGTACTAAAGTTCAAACCTGCCATGAAAGCTCTTGAGGAGGCCCAACAACAAAGCAGGGTCACGACAAAAGCAGTGCCACGACGAGCTGAGTTACCACAGGTCCTCCAGATTCACTCGATCATATATTCTGTTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78531.1
    TTACTCACGGGTGGCCCGCCTGACTGGGGTCGCATCTGAATGGAAATCACCGCCCGAGGGCGGTTTTGCGTCCACTGGGGGTTCAAACAGGTTCTTCTGTAGGCCTGACGAGTACGACAACACGGGGGATTCGAACGGCAGCCTTGCGGCTGCCAGCATCCTCCACCTACTGCCAAGTTCCGGGCCATCAGAACACCGATGCTTAGACCCGCCGCACCTAGGCGCAAGGGGCCAATCTTTCACGTCCTCACAATGACAGACTAGCCGCATGCCAATCAAGCATTATCAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAGATATCCCGTTGCCGAGAGTCGTCATAAATTCATTGGCACGCGCAACAGACGCCGCCCCTCCGAGTTTTTCTTGACAAAAATGCCATCCATCGGTGACGCTCCATATGACTTGGCGCAAACTGCGCCGCGCTAGACGTTCAAACCCGCCATGAAAGTTCCCGAGGCGGCCCGGTCGCAAACCGGGTTCACCACGAGAGCAAAGCCACGGTGAGCCGTGTAACCACGGGTCCTCCGGATTCACTCGATCGTATGTTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78530.1
    ATGGGGTGGCCCCCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAACATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGCGTTCTTCATCGATGCAAGAGTCAAGATATCCGTCTGCCGAGAGTCATCAAAATTCATTGGAATGAACATTAACACACCACACCTCCGACGTTGTGTGGGGCAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGGACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGATATCCCTAGCGAGCCAAATTACCACAAGTCCTCCAGATTCACTCAATCGTTTATTATGTTGTTTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78529.1
    TTTTAGCTCAGCTGGTGGCCCGCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAATATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGGGATTCTTCATCGATGCAAGAGCCAAGATATCCCTCGGCGAGAGTCATCAACAATTCATTGGCATGACCAATAACACACCACACCTCCGACTTTTTTGACAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGCACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGAAATCCCTAGCGAGCCAAATAACCACAAGTCCTCCAGATTCACTCAATCGTATATTCTGCTGTCTCAACAATGTCCTTCGGCAGCTCGCCGT
    Z78527.1
    GGGTGGCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACAATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACCCAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAAGAGCCAGATATCCGTTGGCCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAACACACTGCCACTCCAACTTTTGGTGACAGTAGTGCCATTTGCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGCCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGCGAGCCGAGTAACCACAGGTCATCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78526.1
    ACATAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGGCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAGCACACTGCCACTCCAACTTTTGGTGACAATAATGCCATTGTTCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGGCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGGGAGCCGAGTAACCACAGGTCCTCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78525.1
    TGCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAGTCAACCATCCAAGGGTCATTTTGCCTCCACTGGGGTTCAAACAAGTTATTATATAGGTCCGATAAATACGATAACATGGGGGATTTAAATGACAACCTTGTGGCAGCCAATGTCCTCCACCTCCTGCCAAGTTTCGGGCCATCAAAACATTCCTTAGACCCAACGCGCCAAGGCACAAGGGGCCAATCTTTCGCATCCGCACAATGCAAGCCTAGATATATGGACAGGCATTGGCAAGAGAGACACAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCCGATGCAAGAGCCAAGATATCCGCTTGTCGAGAGTCATCAAAATTATATTCGGCATGCACAGGAGGACACCGCCCCTCCAACTTTTTTGACAATAATGTCATTTGTCGGTGATGCTTCATATGACTTGGCACAAACTGTGCCGTGCTAGAGTTTGAAACCTGCCATGAAAGCTCCCGAGGAGGCCCAACGACAAATTTGGGTCACCACAAAAGCAAAGCCCTCGGCAAGCCGAATAACCACAGGTCCTCCGGATTCACTCGATGTATATTCTGCTATCTCAACA
    Z78524.1
    GCTTAATCTCAGGGTGGCCACCTGACCTGGGGTTGCATCTGAATGAAAATCAACCGCCCAAGGGTCATTTTTCCTTCATTGGGGTTCAAACAAGTTCTTTTATAGGTTCAACAAATACGACAACATGGGGAATTCAAATGGCAGCCTTGCAGCTACCAACATTCTCCACCTCCTGCCGAGTTCCGAGCTATCAAAATACTAATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATATTTCATCCGCACAATGCTAACCTAGCTATATGCTGGGCAAGCATTGGNAGGAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGGATGTAAAGAGCCAAGTTCCGTTGCGAGAGTCATCAAAAAATTCATTAGCATGCACAACAGCATGCCGCCCCTCCGACTTTTTTAACAACAATGCCATTCATCAGGGATGCTTATATGACTTGGCACAAACTGCGCCGTGCTAAAGGTTTGAACCCGCCATGAAAGCTCCTAAGGAGGCCCAATTAAGGTCACCACAAAAGCAAAGCACTGACGAGCCAAGTAACCACATGTCCTCCATATTCACTCAATCATATATTCTACTATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78523.1
    CTTAAACTCAGCGGTGGCCCCGCCTGACCTGGGTCGCAATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGGTTCAAACGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCGCATCCACGCAATGCAAACCAATCTATATGATAGGAAAGCATTGATAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACATATGTCGAGAGTCAATAAAAAATTCATTATGCATGCATGCAAGAGCACGCTTCCACTCTAGCTTTTGGGACAATAATGTCATTCCTCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCACACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAGTTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTGCCGCAGGTTCACCTACGGAAACCTGGTTACG
    Z78522.1
    CTCAGCGGGTGGCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCCCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACGCAATGCAAACCTGTCTATATGACGTATAACCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTTCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCGATCACACCACGATATGTCGAGAGTCATAAAAAATTCATTTGCATGTATGCAATAGCACGCTTCCACTCCAACTTTTTGGACAATAATGTCATTCATCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGGGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78521.1
    GGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGAGTTGTTTGCCTCCATTGGGATTCAAACATGTTTTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGTACCTAGGCACAAGGGGCCATTCTTTCACATCCACGCAATGCAAACCTATCTATATGACATGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTTGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAATTCATTTGCATGCATGCAATAGCACGCTTCCACTCCAACTTTTTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGGCAAATTAGGGTCACCGCCAAAGCCAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78520.1
    AAACTCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCCAGGTCCGACAAACACGACAATATGGGGGATTCACATGACAACCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCGATATGACAGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCACATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTTCGCTGCGTTCTTCATCCGATGCAAGAGCCAAGATATCCCTTGTCGAGAGTCATAAAATATTCTATTGCATGCATGCAATAGCACGCTTCTACTCCAACTTTTTGGACAATAATGTCATTCATCGGTGATGATTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78519.1
    TAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCAACAAACGCGACAATATGGGGGATTCAAATGATAGCCTTATATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACATTAAGTTCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCTATATGATGGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAAATTCATTTGCATGCATGCAGTAGCACGCTTCCACTCCAACTTTTGGACAATAATGTCATTCATCGGCAATGCTTCATATGACTTGGTGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCAGATTCACTCGATCATAT
    Z78518.1
    GGAAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCAAGGGTTGTTTTGCCTCCATAGGGGTTCAAATAGGTTCTTCTGTGGGTCCGGCAAATACGACAACATGTGGGATTTGAACGACAGCCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTTGGGCCATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGAGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGTATGCGTAACACACACCGTCCATCCGAGTTTGGGGGACAATAATGTAATTCGTCGGCGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGCTCAAACCCGCCATAAAAGCTCTCGAGGATGCCCAACTACAAATCAGGGTCACCACAAAAGTTAAGCCTTCGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCGAATATTCTACTATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78517.1
    GCTTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAATGACAGTCTTGCGACTATCAACATGATCCACCTCCTGCCGGGGTTCGGGCCACCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGTCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGGAAGAGCAAGATATCCGTTGGCCGAGAGTCATCAAAAATTTATTGGCATGCACAACAGCACACCGCCCATCCGACTTTTGTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCAAACTGCACCGTGATAGAGGTTCAAACCCGCCATAAAATCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCTGGATTCACTCGATCGTATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78516.1
    TTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCATATCTGAATGGCAATCAATCACTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAACGATAGTCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTCGGGACATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTAAGATATCCGTTGACGAGAGTCATCAAAAAATTATTGGTATGCACAATAGCACACCGCCCATCCGACTTTGGTGACAATAATGTCATTCGTCGGTGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGTTCAAACCCGCCATAAAAGCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCATATATTATACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78515.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCTTCCAAATGGCAATCAACCGCCCATGGGTTGGGTTGATGTGGCACCAATGGGGTTCCGAAGGGTCCTACGAATTATGATGGCCCATCAAATCGGACAACATACGGCATTCGCACAACAGCTTTGGCTCAGCCATTGCCCATCCACCTCTTGTCGAATTTTGGACCATCAAAAGCCCGAGGTCTTAGACCCATTGAACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGACCGAACGGTATTCTTGGACAACCATTGGTGGGAGAGACACAGTACACAACGCCCAGGCAGGCGTGCCCTTAGCCTGACGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTAAAAATGTAACTTGGGGGGCACGAATCAAGTGCCACCCCCACCGGCTACAAGGGAACACCCCTATGCCGATTCTTGTCTCAAAAGACTTGGCGCGAAGCTGCGCCGGGCTAGAAGTTTGATCGCCTTCAAGAACACTCCCAAGGAGGCCCGATGGCAGATCAAAGGCCACCACAGTAGCGAAAGCCCCTGGGAGATCGAAGTACACCGGTCCTCCAGATTAACTCGATCGTTTATTCTACGGTCTCAGCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78514.1
    TAGCGGGTAGTCTCACCTGACCTGGGGTTGCATCCTAATGGCCGTCAACCGCCCATGGGTTGGGTCGACGTGCCTCCGACGGGGTTCGAAAGGGATCATTCCGGCCCATCAAATACGACAACATAGGGCATTCGCACAACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGACGAGAGTTGTCAGACGAATCTTAAGGGGGCACAGCAGCACGCCGGCCCTCCCCGGTCCATCCATGGAGGACGAAGACCCTTGTCGGTTCTATCTCATTTGACTTGGAGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATAGAGAACTCTCCCAAGGAGGTTCGATGGGAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTATCTCGATCGTATATTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78513.1
    CTCAGCGGTAGCTCACTGACTGGGTTGCATCCTAATGGCGTCACCGCCCATGGGTTGGGTCGACGTGCCTCCGAAGGGGTTCGAAAGGGATCGTTCTGGCACATCAAATACGACAACATAGGGCATTCGCACAAAAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACAGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCGAAAATTCTTTCGGGCACAGCAGCATGCCGGCCTCCCCGGCTCCATGGAGGACGATGCCCTTGCCGGTTCTATCTCATTTGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATGGAGAAATCTCCCAATGAGGCTCGATGGCAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTAACTCGATCGTGTATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78512.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGGATCCAAATGACCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAGCTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCGAGAGTTGTCAAAAAATAATTTCATTGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGGGGACGATGCCCTCCGCCGGTTCCATCTCATTTGATTTGGGGCGAAACTGCGCCGGGCCAAGGGTTTCATTTGCCATCAAGAATTCCTCCCAGAGGGTTCGACGGAAATTCACGGTCACCACAGTAGCCAAGGCCCCTGGGAGACCAAACTACACCGGCCCTCCGGATTAATTCGATCGTTTTTTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78511.1
    TCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGCTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTTCATCGGATGCAAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAAAAATTCATGGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAATTCTCCCAAGGAGGTTCGACGGAAATTCACGGCCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78510.1
    TAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAAAACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAATCTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCCATCCGATGCAAAGAGGCGAGATATCCGTTGCGAGAGTTGTCAAAAAATCCGAGGGGGGAACGGAAGAATGCCGGCCCTCCCCGGTTCCATGGGGGGGGACGCCCCCTGCCGGTCCTATCTCATTTGATTGGGGGGGAAATTGCGCCGGGCCAAGGGTCCATTGGCCACCAAGAATTCTCCCAGGGGGGTTCGATGGAAATTCACGGCCACCAAGGGGGGAAAGCCCCTGGGGAGACCAAATTAAACCGGTCCTCCGGTTTAATTCGATCGTTTTTTCGGGGGCTTAAAAAGGAATCCTCCCGAAGGTCACCTCGGAACCCTGGTTAG
    Z78509.1
    TCCTAAAATTAACGGGTGCCTTAACTTGCCGGGGTTAACCAAATGCCCGTCACCCCCCATTGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAAAACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAACCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCTGCTTATCACATTTAAACCTTGTTACGTCCGTTGCCGAGAGTTGTCAGAAAAATTCATGGGGGGGCACGGAGCATGCCGGCCCTCCCCGCTCCATGGAGGACGACGCCCTCCGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78508.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAGTGGCCGTCACCGCCCATGGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAATTCATTGGGGGCACGGCAGCATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78507.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAATGGCCGTCGACCGCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGTGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAAATTCATTGGGGGACGGAAGCACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGAAATTCACGGTCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78506.1
    TCAGCGGGTGGCCTCACCTGACTGGGTTGGCATCCAATGGCCGTCAAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTATGAGCTGAGGGCTCCTTGAACGAGAGTTGTCAGCCCGATTCAAAGGGGGTACGGCGGCATGCCGGTCCTCCCCGGTTCCATGGAGGACGAAGCCCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78505.1
    AAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGACCGAACAGTATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTGACGAGAGTTGTCACAAAGATTCATTGGGGGTATGGTGGTATGCCGGGCCTCCCCGGTTCCATGGAGGACGAAGACCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCACGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAATATCATGGTCACCGCAGTAGCCACCGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78504.1
    TTACTCAGCGGTGGCTCACTGACTGGGGTTGCATCCAATGGCCGTCACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGTTTTGTCGCAGCCACCGTCCGTCCACCTCTTGTCGGATTTGGTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGTCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAGATTCATTGGGGGCACGGCGGGCATGCCGGCCCTCCCCGGCTCCATGGAGGGACGATGCCCTCTGACGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78503.1
    TTATCTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAAAGAGCTGAGACTCCGTTGCCGAGAGTTGTCAAAATAATTCATTGGGGGTACGGTAGCATGCCGGCCCTCCCCGGTTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGTAAATCACGGACACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAAGGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78502.1
    GCGTAAACTCACGGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTACGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTTCGCTGCGTCCTTCCATCGATGCAAGAGCCGAGATTCCGGCCGAGAGTTGTCAATATAAAATTCATTGGGGGGACGGTAGTATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGACCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGTCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGGGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78501.1
    TCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCTCGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACAACACTTATCACAATTTCGCTGCGTCCTTCCATCCGATGCAAGGAGCCGAGATTTCCCGTTTGGCCGAGAGTTGCCCAAAAAAAAATTCATTGGGGGCACGGCAACACGCCGGCCCTCCCCGGCTCCATGGAGGACGATTCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTTCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78500.1
    CTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCGAGATTCCCGTTGGCCGAGAAGTTTTCCCAAAGAAAAATTCATTGGGGGCACGGCAGGACGCCGGCCCTCCCCGGCTCCATGGAGGGCGATGCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCAATTATTTTGCGGTCTCAACAATGAGCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78499.1
    GGTTGCTCACCTGACCTGGGGTCGCATCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACACTTATCGCATTTCGCTGCGTCCTCCATCCGATGCAAAGAGCCGAGATATCCGTTGGCTGAGAGTTGTCAAAAAAATAATTGGGAGAGGTCAAGGCCACATGCCGTCCCCTCCTCGAATTCAACGAAGGGGGTATGTCCTTCACCATTTATGTGTCATATGACTTGGTGCAAAACTGCGACGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGAAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGCGATCTCAACAATGATCCCTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78498.1
    GCTTAAACTCAGCGGGTTGCTCACTGACTGGGGTCGCAACAAATGGCCATCAACTGATCACGGGTTGATGTGCCTCCAGTGGGGTTCACAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAGTTTTGTTGTAATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGTACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCCAAAAATAATTTGGAGAGGGGGAGTGATTGTAGGGTCAAGGGCAAAATGCCGCCCCCTCCTTCGCCATTATTGTGTCATATGACTTGGGGCAAAACTGCCCCGGGCAAGGGTAAGATCGCCGGCAAGAAAACTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCACCCATGGGTGACCGGAGTAAACTGATCCTCTGGAGTAACTCGATCAATTATTATGTGATCTCAACAATGACCTTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78497.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTGTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCGCGTTGCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACACCCGCACCGGCTGAACAGTATGCCTGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAGAGCTGAGATCTCCGTTGCTGAGAGTTGTTTCACAAAAAATAATTTGGAGAGGGGGAGGGAGTGTAGGCTCAAGGCCACATGCCCTTCACCTCCTTGAATTCACAAAGGGCCGTGCCCTTCACCATTTATGTGTCATGTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCGCCCATGGGCGACCAAAGTAAACTGATCCTCTGGATTCACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78496.1
    GCTTAAACTCAGCTGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGGTGCGTCCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTCAAAAAATTAATTGGAGAGGTCAAGGCCACATGCCGGCCCCTCCTCGAATTCACGAAGGGATATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGCAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78495.1
    CACTCAGCGGGTTGCTCACTGACCTGGGGTCGCAACCGAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTGGATTATGCCGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAACAGTGTTGTTGTAGCCTGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCTGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTGAGGGAGAGAGATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTCATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGGATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTCAGAAAAATACTTTGGAGAGGGGGAGAGCGAGTGCAGTGTCAAGGCCACATGCCGCCCCCCCTCCTTGAATACACGAAGGGCTATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACCGCGCCGGACAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCTCCATGGCAACTCTAGGTCACCGCAATAGCAAGCGCCCATGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78494.1
    CTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTACGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGACAGGCGTTCCCTTGGCCTGTTGGCCTCGGGCCCAACTTGCGTTCAAAGACTCGATGTTAACGGGATTCTGCAATTCACACCTGAGATATCCGAAGTTGAGAGTGTTACCTCAATAATGGGGAGTTGGTCAAGGCAGTATGCCGTCCCCCTTCATTAATTATGGGTCATTTGATTTGGCGCATTACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGTTCCCAAGGAGGCCCGATGGTAAATCTAGGTCACTGCAACAGTAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGACCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78493.1
    GGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCACACCGGATGACCAGTACGTTTGGACAGCCTCATTAAGGGAGAGATATGAGTCGCAATGTCCAGGCAGGCGTCCCCTTGGGCTGATGGCCTCGCGCGCAACTTGCGTTCAAAGATTCGATGGTTCACGGGATTCTGCAAGGCACACCATTTATCGAATCTCGCTGCGTTCCTTCATCGATGCATGAGCCGAGATATCCGTTGCTGAGAGTTGTTACCTAAATACTTGGGGGCTGGTCAAGGAAGAATGCCGCCCCCCTTCACCACTTATGTAAAATTTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGCTCCCAGGGAGGCCCGATGGAAATTCTAGGTCACTGCAACAGAAAATGCCCATGGGTGACCAAAGTAAACCGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78492.1
    TATTAAACTCAGCGGTGTGCCTCACCTGACTGGGATCGCAACCAAATGGCCAATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCTGATTAAGCTGGCCCATGTAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAATTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCGTGCAGCTCTTAGACCCCCCGTCCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTCCACGGGATTCTGCAAATCACACCAGTATCGCATTTCGCTGGGTTCTACATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGAGTGGGTCAAGGCAGTACGCTCCCCCCCTTCACCAATTATGTGACATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78491.1
    GCTTATCTCAGCGGGTTGCTCACCTGACTGGGATCGCAACAAATGGCCATTCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTCCAAAGGGTCTTCTGATTATGCTGGTCCATGTAACACAACAACCCGGGGCATTCGCACAACATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAATATAATTTGGAGCGGGTCAAGGCAGCACGCCTCCCCCCAACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCTAAGGAGGCCCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78490.1
    TCAGCTGGTTGCTCACCTGACTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGACGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAAATTATGCTGGCCCATCTAATACGACAACGCGGGGCATTCGCACAACACTTGTTGCAGAACAGTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCGTTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCGAGGCAGCATGCCGCCCCCTTCACCAATTATGTGCCATATGACTTGGCGCAAAACTGCGCCGGACTAGGGTTCAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78489.1
    GCCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCATGGGGTTCAAAAGGGTCTTCAGATATGCTGGGCCCATCTAATACGACAACTCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAATCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACCCACATCCCCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATCTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTCGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCTCGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78488.1
    AGCCGGTTGCTCACCTGACTGGGGTCGCAACAAATGTCCATCAACTGATCATGGGTTGATGGGCCTCCGATGGGGTTCAAAAGGGTCTTCAGATTATGATGGCCCATCTAATACGACAACCCGGGGGATCCGCACAACAGTTTGGTTGTGGCGTTGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGGCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTCCATCCCGTTGATGAGAGTTGTTAGAAAATAATTAGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCACCAATTATGTGTCGTATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTAAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAACAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTGCGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACAG
    Z78487.1
    TTACTCAGCGGTTGCTCACCTGACTGGGGTCGCAACAAGTGGCCCCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTACGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGCATTTAGGACCATCGATAGCCCATGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCGAACTCACATCCGCACCGGCTGAACAATATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTAAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGAGTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78486.1
    TCAGCGGGTTGCCTCACCTGACCTGGGGTCGCTACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGGCACTGCAATAGCAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATAATCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGATTCACCTACGGAAACCTCGTGACG
    Z78485.1
    TACTAAATCAGGGGTTGCTCAGCTGACTGGGGTCGCAACACATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAGGGGTCTTTAGATTGTGCTGGCCCGTCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGGCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGCCCAAACTCACATCCGCACCGGCTGCACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCCCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGATATCTCAACAATGATCCATCCGCAGATTCACCTTCGGACACCAGGTTCAG
    Z78484.1
    AAAACTCAGGGAGTTGCTTCACCTGACTTGGGGTCGCAACCCAATGGACCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTGTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGGGAGGGTAGGAAACATGCCGCCCCTTCACAATTATGTGTCATATGACTTGGGCCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGGAGGCTCGATGGTAAACTCTCGGTCACTGCAATAGCCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCCCAGGTTCACCTACGGAAACCTTGTTACG
    Z78483.1
    TGCTAAACTCAGGGGTTGCTTCACTTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGATATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATTGTTCACGGGATTCTGCAATTCACACCATTTATGATTTCGGTTTGGTCTTCATCCATGCAAGAGGCCAGATTTCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78482.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACGCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGNTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGAAGAGGCGAGATATCCGTTGCNGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTACGGTTTAGATCGCCAGCAAGAAAGCTCCCAGGAGGCTCGATGGCAAATCTCGGTCACTGCAGTAGA
    Z78481.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTCCGTTGCTGAGAGTTGTTAAAAAAATGATTTGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78480.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCGCACAACATTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78479.1
    ACTCAGAGGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCGTCAACTGATCTTGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAGACTCGATGGTTCACGGGATCTCTGTATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGGCAAGGCAAGGCAGCATCCCTCCCCTTCCAGTAATGTCTCATATAAATTGGCGCAAATCTGCGCCGGGCAAGGGTTTAGTTCGGCTCGATGGCAAATCTAGGTCACTGAAACAGCAAATGCCCATGGGTGACCAAAGTAACCTGATCCTCCTGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78478.1
    GCCTAAACTCAGCGGGTGCCTCATCTGACCTGGGTCGCAACCAAATGTCCTCAACTGATAATGGGTTGATGGGCCTCCAATGGGGGTCAAAAGGGTCTTTAGATTATGCTGACCCATCTAAACGACAACCGGGGCATTCGCACAACAGTTCTGTTGTCACATAGCTCGTCCACCTCTTGCCGTATTTAGGACCATCCAACCCATACAGCTCTTAGACCCCCGTACCAAAGGACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGCAGCCTCATTAAGGGAGAGATATGTCTCACAATTTCCAGGCAGGCGTCCCCTTGACCTGATGGCCTCGGACGCAACATGCGTTCAAAACTCGATGTTCAGGGATCTGCTATCACGCCATTATCAATATTTGGGGGAGCCATGGCAGCATGCCGCCCCTTCCAGTTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAACTCTAGGCCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCCCCAGATTAACTCGATCAATTATTATGTGATCTCAACACTGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78477.1
    GCATAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCCGCAACTTGCGTTCCAAAGACTCGATGGTTCAGGGATTCTTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAGAGCCGAGATACTCGTTGCTGAGAGTTGGTTAACAAATATTTGGGGAGGGCAAGGCAGCATGCCGCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGCAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78476.1
    GGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGATTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCGCAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTGGATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCGGATGGCCTAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGCCAAGGCAGCATGCCGCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78475.1
    ACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCATTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGACGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78474.1
    AAGCTCAACGGGTTGCTTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAGAAGGGTCTGTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGTGTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTCAAAAATATAATTTGGGGAGGGTCAAGGCAGGATGCCGCCCTTCCAATTATGTGTCATATGATTGGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGAAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTACGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78473.1
    CCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGGCAAGAACAAGGGGCCAAACTCACATCCGCACCGACTGAACAGTATGTATGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAACAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCACCAAGAAAGCTCCCATGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78472.1
    GCTTAAACTCAGCGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATTGGCCTCTAAAGGGGTTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGAACTACAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAAGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATAGACAGCCTCATTAAAGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATTCCGTTGCTGAGAGTTATCAAAAAATAATTTGGGGAGGGTCTTAGACTGCATGCCGCCCCCTTCCAATTATGTGTGAAATGACTTGGCGCAAGACTGCGCCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78471.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTACTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGATGTTGAGAGTTGTCAAAAAATTAATTTGGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCATTATGTGTCGGGGTGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGACAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACCGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78470.1
    AACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTAGTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTAAAAAAATAATTTGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78469.1
    AACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGCACTATCATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCTAGACCGCATGCCGCCCCCTTCCAATTATGTGTGATATGACTTGGCGCAAGACTGCGCCGGGCAACGAATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78468.1
    AACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCACACAATAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTCCATCAAAAGCCCAGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGAGTTCCGGAGCTGAGAGTCGTTAAAAAAATGATATGGGGAGGGTCAAGGCAGCATGCTGCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTCAGATCGCCATCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATTCCCATGAGTGACCAAAGTAAACTGATCCTCCAGATTAACTCAATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78467.1
    TCAGCGGGTTGCCTCACCTGACCTGGGTCGCAAACATATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAGAGGGTCTTTAGATTATGCTGGCCCAGCTAATACGACAACCCGGGTCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCACCAAAGGCACATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTACCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78466.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACGCAAGAACAAGGGGCCAAACTCACATCCGCACAGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCAAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTACAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGAAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78465.1
    GCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACTCGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTTCGTTCAAAGACTCGATGGTTTACGGGATTCTTCAATTCACACCAATTATCGCATTTCGCTGCGTTCCTCATCGATTCAAGAACCGAGAAATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78464.1
    GCTTAAACTCAGCGGGTTGCTCACCTGACCTGGGGTCGCACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATTACGAAAGCCCGGGGAATCCGAACAGAAATTTGGTTGAAGAATAGTCCGCCCACCTCTTGCCGTTTTTAGGACCACCAAAAGCCCATGAAGTTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGAACCGGTTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATTTGATTCGAAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGAACTTGGCGTTCAAAGACTCGATGGTTCACGGGATTCTGAAATTCACACCATTTATCGGATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCTCCGGGCAAGGGTTTAGATCTCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGATCAATTATTATGTGATCTCAACAATGACCCTTCCGCTCACCTACGGAAACCTTGTTACG
    Z78463.1
    GCTTAAACTCAGCGGGTTGCTCACCTGATCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGTATTCGCACAGCAATTTTGTTGCAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAACATGCCGCCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTCCCCCGGGCGAGGGTTTAGATCCCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGAACAATGATTATGTGATCTCAACAATGAACCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78462.1
    ATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTGCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGGTCACATCCGGAGACCTCGTGACG
    Z78461.1
    TTAACTCAGCGGCTTGCTCACTGACTGGGTCGCAACAATGGCATCAACTGATCATGGGTTGATGGGCTCCAATGGGGTTCAAAGGGTCTTCAGATTACGGTGGCCCATCTTATACGACAACCCGGGGGCTTTGGACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGGATGTATGGACAGGCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGGAATTCACACCATTTATCGGATTTCGGTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTGAAAAAATTATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGTCATTTGACTGGGGGAAAAATTGCGCCGGGTAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTCCAGCAGCAACTGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78460.1
    TAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTGGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78459.1
    AAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAATTTATCGCAATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78458.1
    CAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCGTACGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78457.1
    CTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACGATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGGTCCTTCATCCGATGAAAGAGCCGAGATATCCGTTGTTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGCCATTTGACTGGGGGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTACAACAGAAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78456.1
    GCTAAACTCAGGGTTGCCTCACCTGACCTGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCCAGACCGGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTCCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78455.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGAAGCTCTTAGACCCCCCGTGCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCAACCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGGATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGAAAGGGTTTAGTTCTCGCCAAAAAGAAAGTTCCCAAGGGGGTTCGATGGAAATTCTAGGTCACTCCAAAAGAAATTGTTCATGGGTGACCAAAGTAAAAAGTTCCTCCAGATTAACTCGATCAATTTTTATGTGATCTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTGGTTACG
    Z78454.1
    GTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGNGGTTGATGGGCCTCCAATGGGTTCAAAAGGGGCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAGAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGTCTGATGGCCTCGGCCGCAACTTGCGTTCACAGACTCGATGGTTCACGGGATTCTGCAAATCACACCATTTATCGCATGAGTTCCGTTGATGAGAGTTGTTAACAATAGTGGGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCGATTATGTGTCAAATGACTTGGAGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGTAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78453.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAGTTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78452.1
    TGCTAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACAGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCATCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATTCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAAAAAGAAAGCTCCCAAGGAGGCTTGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAATAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78451.1
    GCTCAGCTGGTTGCTCACTGACTGGTGTCGCATCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCTATAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGGCTTCGCACAACAATTTTATTGTAGGGATAGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAATCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGGTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATTTGACTCGCGATGCCCAGGGAGGGGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGATTTCGCTGGTTCTTCATCGATGAAAGAGCGAGATTTCCGTTGTTGAGAGTTGCAAAAAATACTTTGGGGAGGGTCAGGGCAGGATGCCCGCCCCCTACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCACGGGTTTAGATCTCGCCAGCAAGAAAGCTCCCAGGGGGGCTCCATGGAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTACACCTACGGAAACCTTGTTACG
    Z78450.1
    CTCAGGGGTTGCTCAGCTGATCTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGTCCTCCAATGGGGTTCAACAGGGTCTACAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGCACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCAGTTTTGAGGACCATCATATGCCCATGCAGTTCTTAGACCCCCCGGACCAAAGAACAAGGGGCCACACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAGGGGAGAGATATGACTCGGAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGGCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGTAAATGCTCATGGATGACCAAAGTAAACAGATCCTCCAGCTTAACTCGATCAATTATTATGTGATATCAGCAATGATCCTTCC
    Z78449.1
    GCATAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTAGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGGAGGGTCAGGGCAGGATGCCACCCTTCACCAATTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78448.1
    CCTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCANCCAAATGGCCATCAACTGATCATGGGCTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTNTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGCGTTCTTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAAGGGAAGGATGCCACCCCCCTTCACCAATTATGTGTCATACGACTTGGCGCAAAACTGCGCCGGGAAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78447.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATAGGATGCCACCCCCTTCACCACTTATGTGTCATTTGACTGGGGGCAAAATTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAAAGCTCCCAGGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGGAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78446.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGAGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGGCCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCTGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTATCGCAATTTCGCTGGTCCTCCATCGATGAAGAGCCGAGTATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGGAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78445.1
    ACAGCTGTTGCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTGCACGGGATTCTGCAATTCACACCATTTATCGCATTTGCCGAAGACTGAGTCTCCGTTGCTGAGAGTTGTTAAAAAAATAGTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCAAACGACTTGGAGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78444.1
    AATGGCCATCAACTGACATGGGTTGATGGGCTCCAATGGGGTCCAAAAGGGTCTTCAGATATCGGTGGCCCATCTTATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCATCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTAAGGGAAGAATGCCACCCCCTCCACCAATTTTGTGTAATTTGATTGGGGGAAAAATTGAGGGTTTAGATATCGCCAACAAGGATGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGTTCACCCTACGGAAACCTTGTTACG
    Z78443.1
    CCTCACCTGACTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTATGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGGTGTGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGTGGCCCAAAACTCTGCCGGGCAAGGGTTTAGATATCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78442.1
    ACTCAGCGGGTTGCTCAGCTGACTGGGGTCGCAACAAATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTTGTTGTGGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCGGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTTTGGATAGCCTCGTTAGGGGAGAGATATGATTCCCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCAAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78441.1
    CTCAGCGCGTTGCTCAGCTGACTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCACAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAACTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATATCGGCAATGACCTTCC
    Z78440.1
    TGCTAAACTCAGGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGCTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCAACGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGATATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAATGCTCCCAAGGAGGGTCGATGGAAATTCTAGGTCACTACAACAGAAATTGCCCATGGGTGACCAAAATATGCAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGGAGGTCCACCTACGGAAACCTTGTTACG
    Z78439.1
    GGCCCAACTAAACGGCCAACCGGGGCATTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTTCCGTATTTAGGACCATCCAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATTACCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCCAACTCCGCCGGGCAGGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATG



```python
# This prints 94 reverse complements 
```


```python
records = [
    rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
                          for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
]
```


```python
# This adds rc before every ID so that we know it is a reverse complement
```


```python
len(records)
```




    94




```python
records = [
    rec.reverse_complement(id = "rc" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
]
```


```python
len(records)
```




    18




```python
records = (
rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
)

SeqIO.write(records, "rev_comp.fasta", "fasta")
```




    18


## Multiple Sequence Alignment - Parts 1, 2, 3, 4, & 5

In this analysis, we worked with taking multiple versions of something and aligning the sequences. We were able to look at how closely related DNA strands were through pairwise alignment. 

```python
#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/PF05371_seed.sth
# This is the file we are using
```


```python
from Bio import AlignIO
```


```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
# This is our alignment 
print(alignment)
```

    Alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
print("Alignment length %i" % alignment.get_alignment_length())
```

    Alignment length 52



```python
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
for record in alignment:
    if record.dbxrefs:
        print("%s %s" % (record.id, record.dbxrefs))
```

    COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
    COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
    Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
    COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']



```python
# We can see all the annotations for all the alignments 
for record in alignment:
    print(record)
```

    ID: COATB_BPIKE/30-81
    Name: COATB_BPIKE
    Description: COATB_BPIKE/30-81
    Database cross-references: PDB; 1ifl ; 1-52;
    Number of features: 0
    /accession=P03620.1
    /start=30
    /end=81
    Per letter annotation for: secondary_structure
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA')
    ID: Q9T0Q8_BPIKE/1-52
    Name: Q9T0Q8_BPIKE
    Description: Q9T0Q8_BPIKE/1-52
    Number of features: 0
    /accession=Q9T0Q8.1
    /start=1
    /end=52
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA')
    ID: COATB_BPI22/32-83
    Name: COATB_BPI22
    Description: COATB_BPI22/32-83
    Number of features: 0
    /accession=P15416.1
    /start=32
    /end=83
    Seq('DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA')
    ID: COATB_BPM13/24-72
    Name: COATB_BPM13
    Description: COATB_BPM13/24-72
    Database cross-references: PDB; 2cpb ; 1-49;, PDB; 2cps ; 1-49;
    Number of features: 0
    /accession=P69541.1
    /start=24
    /end=72
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPZJ2/1-49
    Name: COATB_BPZJ2
    Description: COATB_BPZJ2/1-49
    Number of features: 0
    /accession=P03618.1
    /start=1
    /end=49
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA')
    ID: Q9T0Q9_BPFD/1-49
    Name: Q9T0Q9_BPFD
    Description: Q9T0Q9_BPFD/1-49
    Database cross-references: PDB; 1nh4 A; 1-49;
    Number of features: 0
    /accession=Q9T0Q9.1
    /start=1
    /end=49
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPIF1/22-73
    Name: COATB_BPIF1
    Description: COATB_BPIF1/22-73
    Database cross-references: PDB; 1ifk ; 1-50;
    Number of features: 0
    /accession=P03619.2
    /start=22
    /end=73
    Per letter annotation for: secondary_structure
    Seq('FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA')



```python
# This is how to write an alignment file
```


```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment 
```


```python
align1 = MultipleSeqAlignment(
...    [
...        SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha"),
...        SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta"),
...        SeqRecord(Seq("ACTGCTAGDTAG"), id="Gamma"),
...    ]
... )
>>> align2 = MultipleSeqAlignment(
...    [
...        SeqRecord(Seq("GTCAGC-AG"), id="Delta"),
...        SeqRecord(Seq("GACAGCTAG"), id="Epsilon"),
...        SeqRecord(Seq("GTCAGCTAG"), id="Zeta"),
...    ]
... )
>>> align3 = MultipleSeqAlignment(
...    [
...        SeqRecord(Seq("ACTAGTACAGCTG"), id="Eta"),
...        SeqRecord(Seq("ACTAGTACAGCT-"), id="Theta"),
...        SeqRecord(Seq("-CTACTACAGGTG"), id="Iota"),
...    ]
... )
```


```python
my_alignments = [align1, align2, align3]
```


```python
my_alignments
```




    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7f6900e9ad50>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7f6900e9ae50>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7f6900e9af50>]




```python
print(my_alignments)
```

    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7f6900e9ad50>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7f6900e9ae50>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7f6900e9af50>]



```python
from Bio import AlignIO
AlignIO.write(my_alignments, "my_example.phy", "phylip")
```




    3




```python
alignments = AlignIO.parse("my_example.phy", "phylip")
```


```python
for alignment in alignments:
    print(alignment)
    print()
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma
    
    Alignment with 3 rows and 9 columns
    GTCAGC-AG Delta
    GACAGCTAG Epsilon
    GTCAGCTAG Zeta
    
    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota
    



```python
alignments = list(AlignIO.parse("my_example.phy", "phylip"))
```


```python
last_align = alignments[-1]
```


```python
print(last_align)
```

    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota



```python
# That is all align 3
```


```python
first_align = alignments[0]
```


```python
print(first_align)
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma



```python
# That is all from align 1
```


```python
from Bio import AlignIO
```


```python
count = AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.aln", "clustal")
```


```python
# We are converting our file to a clustal file
```


```python
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")
```


```python
count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal")
```


```python
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
# We can also make a phylip format 
AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip")
```




    1




```python
AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip-relaxed")
```




    1




```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
name_mapping = {}
for i, record in enumerate(alignment):
    name_mapping[i] = record.id
    record.id = "seq%i" % i 
```


```python
print(name_mapping)
```

    {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', 3: 'COATB_BPM13/24-72', 4: 'COATB_BPZJ2/1-49', 5: 'Q9T0Q9_BPFD/1-49', 6: 'COATB_BPIF1/22-73'}



```python
AlignIO.write([alignment], "PF05371_seed.phy", "phylip")
```




    1




```python
# We can also do slicing alignments 
```


```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
print("Number of rows: %i" % len(alignment))
```

    Number of rows: 7



```python
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
print(alignment[3:7])
```

    Alignment with 4 rows and 52 columns
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
print(alignment[2, 6])
```

    T



```python
# T is the sixth character in the second row
# You start counting the rows from 0 (not 1) 
```


```python
# This is another way to do it
print(alignment[2].seq[6])
```

    T



```python
# This is how to pull out a single column 
print(alignment[:, 6])
```

    TTT---T



```python
# This is how to select a range 
print(alignment[3:6, :6])
```

    Alignment with 3 rows and 6 columns
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49



```python
# This is all the alignments through 6
print(alignment[:, :6])
```

    Alignment with 7 rows and 6 columns
    AEPNAA COATB_BPIKE/30-81
    AEPNAA Q9T0Q8_BPIKE/1-52
    DGTSTA COATB_BPI22/32-83
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49
    FAADDA COATB_BPIF1/22-73



```python
print(alignment[:, 6:9])
```

    Alignment with 7 rows and 3 columns
    TNY COATB_BPIKE/30-81
    TNY Q9T0Q8_BPIKE/1-52
    TSY COATB_BPI22/32-83
    --- COATB_BPM13/24-72
    --- COATB_BPZJ2/1-49
    --- Q9T0Q9_BPFD/1-49
    TSQ COATB_BPIF1/22-73



```python
# We can pull out certain parts of sequences 
print(alignment[:, 9:])
```

    Alignment with 7 rows and 43 columns
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
edited = alignment[:, :6] + alignment[:, 9:]
```


```python
print(edited)
```

    Alignment with 7 rows and 49 columns
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
# We got rid of NYA 
```


```python
edited.sort()
```


```python
print(edited)
```

    Alignment with 7 rows and 49 columns
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49



```python
# This sorts based on the name of the ID descending alphabetically
```


```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment 
# We have to load everything back up 
```


```python
alignment = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTCCTA"), id="seq1"),
        SeqRecord(Seq("AAT-CTA"), id="seq2"),
        SeqRecord(Seq("CCTACT-"), id="seq3"),
        SeqRecord(Seq("TCTCCTC"), id="seq4"),
    ]
)
```


```python
# This gives us our alignment 
print(alignment)
```

    Alignment with 4 rows and 7 columns
    ACTCCTA seq1
    AAT-CTA seq2
    CCTACT- seq3
    TCTCCTC seq4



```python
# Now we can are working with substitutions 
substitutions = alignment.substitutions 
```


```python
# This shows our substitutions 
print(substitutions)
```

        A    C    T
    A 2.0  4.5  1.0
    C 4.5 10.0  0.5
    T 1.0  0.5 12.0
    



```python
m = substitutions.select("ATCG")
```


```python
print(m)
```

        A    T    C   G
    A 2.0  1.0  4.5 0.0
    T 1.0 12.0  0.5 0.0
    C 4.5  0.5 10.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
# Now we have added G; it was initially left off bc there were no Gs in our sequence
```


```python
m = substitutions.select("ACTG")
#This is changing the order
```


```python
print(m)
```

        A    C    T   G
    A 2.0  4.5  1.0 0.0
    C 4.5 10.0  0.5 0.0
    T 1.0  0.5 12.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
import Bio.Align.Applications 
```


```python
dir(Bio.Align.Applications)
```




    ['ClustalOmegaCommandline',
     'ClustalwCommandline',
     'DialignCommandline',
     'MSAProbsCommandline',
     'MafftCommandline',
     'MuscleCommandline',
     'PrankCommandline',
     'ProbconsCommandline',
     'TCoffeeCommandline',
     '_ClustalOmega',
     '_Clustalw',
     '_Dialign',
     '_MSAProbs',
     '_Mafft',
     '_Muscle',
     '_Prank',
     '_Probcons',
     '_TCoffee',
     '__all__',
     '__builtins__',
     '__cached__',
     '__doc__',
     '__file__',
     '__loader__',
     '__name__',
     '__package__',
     '__path__',
     '__spec__']




```python
from Bio.Align.Applications import ClustalwCommandline
```


```python
#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.aln
#This is where our file is coming from
```


```python
from Bio import AlignIO
```


```python
align = AlignIO.read("opuntia.aln.txt", "clustal")
```


```python
# This gives us our alignment for the opuntia file
print(align)
```

    Alignment with 7 rows and 906 columns
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191



```python
from Bio import Phylo
```


```python
#https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.dnd
#We will use this file now
```


```python
tree = Phylo.read("opuntia.dnd.txt", "newick")
```


```python
# Here is our tree
Phylo.draw_ascii(tree)
```

                                 _______________ gi|6273291|gb|AF191665.1|AF191665
      __________________________|
     |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
     |                          |__|
     |                             |_____ gi|6273289|gb|AF191663.1|AF191663
     |
    _|_________________ gi|6273287|gb|AF191661.1|AF191661
     |
     |__________ gi|6273286|gb|AF191660.1|AF191660
     |
     |    __ gi|6273285|gb|AF191659.1|AF191659
     |___|
         | gi|6273284|gb|AF191658.1|AF191658


```python
# This is sequence alignment pt 4 video
```


```python
from Bio import Align
```


```python
aligner = Align.PairwiseAligner()
```


```python
aligner = Align.PairwiseAligner(match_score = 1.0)
#We are setting our scores to one point each
```


```python
target = "GAACT"
```


```python
query = "GAT"
```


```python
score = aligner.score(target , query)
```


```python
score
```




    3.0




```python
# This is because we have 3 overlaps
```


```python
alignments = aligner.align(target, query)
```


```python
for alignment in alignments:
    print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    
    target            0 GAACT 5
                      0 |-|-| 5
    query             0 G-A-T 3
    



```python
aligner.mode = "local"
```


```python
target = "AGAACTC"
```


```python
query = "GAACT"
```


```python
score = aligner.score(target, query)
```


```python
score
```




    5.0




```python
alignments = aligner.align(target, query)
```


```python
for alignment in alignments:
    print(alignment)
```

    target            1 GAACT 6
                      0 ||||| 5
    query             0 GAACT 5
    



```python
# it is ignoring the first and last nucleotide and counting GAACT bc those are the same; they overlap
```


```python
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: 0.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: local
    



```python
aligner.algorithm
```




    'Smith-Waterman'




```python
# We are using Smith-Waterman
```


```python
# We can also set a significance 
aligner.epsilon
```




    1e-06




```python
from Bio import Align
```


```python
aligner = Align.PairwiseAligner()
```


```python
target = "GAACT"
```


```python
query = "GAT"
```


```python
alignments = aligner.align(target, query)
```


```python
alignment = alignments[0]
```


```python
alignment
```




    <Alignment object (2 rows x 5 columns) at 0x7f22e4527310>




```python
alignment.score
```




    3.0




```python
# We can pull out the target sequence
alignment.target 
```




    'GAACT'




```python
# We can pull out the query sequence
alignment.query
```




    'GAT'




```python
print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    



```python
# We can also get coordinates 
alignment.coordinates 
```




    array([[0, 2, 4, 5],
           [0, 2, 2, 3]])




```python
# We can also get the length of our alignment
len(alignment)
```




    2




```python
alignment.shape
```




    (2, 5)




```python
# The alignment is 2 rows and 5 columns 
```


```python
aligner.mode = "local"
```


```python
local_alignments = aligner.align("TGAACT", "GAC")
```


```python
local_alignment = local_alignments[0]
```


```python
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
# This drops of the two end letters since its a local alignment 
```


```python
local_alignment.shape
```




    (2, 4)




```python
# The alignment is 2 rows and 4 columns 
```


```python
# Now lets look at global
```


```python
aligner.mode = "global"
```


```python
aligner = Align.PairwiseAligner(match = 1.0, mismatch_score = -10)
```


```python
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: -10.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: global
    



```python
alignments = aligner.align("AAACAAA", "AAAGAAA")
```


```python
# We can look at the length
len(alignments)
```




    2




```python
print(alignments[0])
```

    target            0 AAAC-AAA 7
                      0 |||--||| 8
    query             0 AAA-GAAA 7
    



```python
# It would rather put gaps in than have a mismatch; there are 6 positives and 2 zero points 
```


```python
print(alignments[1])
```

    target            0 AAA-CAAA 7
                      0 |||--||| 8
    query             0 AAAG-AAA 7
    



```python
# This just reverses the G and the C
```


```python
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
local_alignment.sort()
```


```python
print(local_alignment)
```

    target            0 GA-C 3
                      0 ||-| 4
    query             1 GAAC 5
    



```python
# This is sequence alignment video Pt 5
```


```python
from Bio import Align
```


```python
from Bio.Seq import reverse_complement 
```


```python
target = "AAACCC"
```


```python
query = "AACC"
```


```python
aligner = Align.PairwiseAligner(mismatch_score= -1, internal_gap_score = -1)
```


```python
aligner.score(target, query)
```




    4.0




```python
aligner.score(target, reverse_complement(query))
```




    0.0




```python
aligner.score(target, reverse_complement(query), strand = "-")
```




    4.0




```python
aligner.score(target, query, strand = "-")
```




    0.0




```python
alignments = aligner.align(target, query)
```


```python
len(alignments)
```




    1




```python
# This means that there is 2 alignments
```


```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
# This is the first alignment 
```


```python
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	+	1	5	0	1	4,	0,
    



```python
# A bed file is an alignment file 
```


```python
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	-	1	5	0	1	4,	0,
    



```python
# This gives us the same thing as our normal alignment 
```


```python
alignments = aligner.align(target, query, strand = "-")
```


```python
len(alignments)
```




    2




```python
print(alignments[0])
```

    target            0 AAACCC----  6
                      0 ---------- 10
    query             4 ------GGTT  0
    



```python
print(alignments[1])
# This will be the reverse of alignment 0
```

    target            0 ----AAACCC  6
                      0 ---------- 10
    query             4 GGTT------  0
    



```python
aligner.left_gap_score = -0.5
```


```python
aligner.right_gap_score = -0.2
```


```python
aligner.score(target, query)
```




    3.3




```python
alignments = aligner.align(target, query)
```


```python
len(alignments)
```




    1




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
# This is our alignment score compared to our reverse complement 
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
print(alignments)
```

    <Bio.Align.PairwiseAlignments object at 0x7f22e4534f10>



```python
aligner.score(target, reverse_complement(query), strand = "-")
```




    3.3




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             4 -AACC- 0
    



```python
# This is for the positive strand 
aligner.score(target, query, strand = "+")
```




    3.3




```python
# mapping our DNA to the reverse complement gives us a worse alignment than mapping to our standard string
# If we map it to the reverse complement we get an alignment score of -2 
```


## Blast 

In this section, we worked with Blast. 

```python
from Bio.Blast import NCBIWWW
```


```python
NCBIWWW.email = "neh1023@gmail.com"
```


```python
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
```


```python
# This is the file we are using
#https://github.com/biopython/biopython/blob/master/Doc/examples/m_cold.fasta
# This is a fasta file for cold aclimation protein 
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("m_cold.fasta", format = "fasta")
```


```python
print(record)
```

    ID: gi|8332116|gb|BE037100.1|BE037100
    Name: gi|8332116|gb|BE037100.1|BE037100
    Description: gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence
    Number of features: 0
    Seq('CACTAGTACTCGAGCGTNCTGCACCAATTCGGCACGAGCAAGTGACTACGTTNT...TTC')



```python
# This is our fasta record
```


```python
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-8-eebaa3d54861> in <module>
    ----> 1 result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
    

    NameError: name 'NCBIWWW' is not defined



```python
with open("m_cold_2.fasta", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
# We are going to read in our result handle results and write our result handle
# Then we are going to close the connection 
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-9-a5ee02ca1283> in <module>
          1 with open("m_cold_2.fasta", "w") as out_handle:
    ----> 2     out_handle.write(result_handle.read())
          3 result_handle.close()
          4 # We are going to read in our result handle results and write our result handle
          5 # Then we are going to close the connection


    NameError: name 'result_handle' is not defined



```python
from Bio.Blast import NCBIXML
# Now we are going to change it to a XML file
```


```python
result_handle = open("my_blast.xml")
```


    ---------------------------------------------------------------------------

    FileNotFoundError                         Traceback (most recent call last)

    <ipython-input-11-5d5d42105f61> in <module>
    ----> 1 result_handle = open("my_blast.xml")
    

    FileNotFoundError: [Errno 2] No such file or directory: 'my_blast.xml'



```python
blast_record = NCBIXML.read(result_handle)
```


```python
# Now we can parse through it
```


```python
E_VALUE_THRESH = 0.04
```


```python
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGHTMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    <ipython-input-1-b6a7171bc552> in <module>
    ----> 1 for alignment in blast_record.alignments:
          2     for hsp in alignment.hsps:
          3         if hsp.expect < E_VALUE_THRESH:
          4             print("****ALIGHTMENT****")
          5             print("sequence:", alignment.title)


    NameError: name 'blast_record' is not defined



```python
# This is the issue that I came to your office about and you told me to leave this notebook as it is bc the code cannot be run twice
```


## Challenge 1 



```python
from Bio.Blast import NCBIWWW
```


```python
NCBIWWW.email = "neh1023@gmail.com"
```


```python
result_handle = NCBIWWW.qblast("blastn", "nt", "4436")
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("MSH2.fasta", format = "fasta")
```


```python
print(record)
```

    ID: NC_000002.12:47403067-47709830
    Name: NC_000002.12:47403067-47709830
    Description: NC_000002.12:47403067-47709830 Homo sapiens chromosome 2, GRCh38.p14 Primary Assembly
    Number of features: 0
    Seq('AAGCTGATTGGGTGTGGTCGCCGTGGCCGGACGCCGCTCGGGGGACGTGGGAGG...TAA')



```python
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
```


```python
# we are finding the record sequence 
```


```python
with open("MSH2.fasta", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()

# now we are reading in the results 
```


```python
from Bio.Blast import NCBIXML
```


```python
result_handle = open("my_blast.xml")
```


```python
# now we are reading in the handle
blast_record = NCBIXML.read(result_handle)
```


```python
# now we are looking at E values of 0.04 and up
E_VALUE_THRESH = 0.04
```


```python
# This gives us the alignment 
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGNMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0.75] + "...")
```


```python
# i was working with the MSH2 gene
```











## Open CV - Parts 1, 2, & 3 

In this section, we looked at computer vision for phenotyping. Computer Vision helps us process images with Python to write a script for measurements. 


```python
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
import cv2
```


```python
# Open CV is another package we will use for phenomics 
# Computer vision helps us process imqges with python to write a script for measurements 
```


```python
img = cv2.imread("Puppy.jpg")
```


```python
type(img)
```




    numpy.ndarray




```python
# This is a path that does not exist 
img_wrong = cv2.imread('wrong/path/doesnot/abcdegh.jpg')
```


```python
type(img_wrong)
```




    NoneType




```python
# NoneType means we have loaded an incorrect image; this is important to check for
```


```python
# This will give us our image
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fac54c09cd0>




![output_8_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/a249d795-ba2f-4f5c-91fa-b276b32cae17)



```python
# The picture turn blue bc the package we are using expects a different order of colors, but we can change this
fix_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(fix_img)
```




    <matplotlib.image.AxesImage at 0x7fac51b3d250>




![output_10_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/a3a8e522-c7a5-4d2a-bc7f-582808cfb5d3)




```python
# Now the image is fixed 
```


```python
# We can also make gray images
img_gray = cv2.imread("Puppy.jpg", cv2.IMREAD_GRAYSCALE)
img_gray.shape
```




    (3888, 5184)




```python
plt.imshow(img_gray)
```




    <matplotlib.image.AxesImage at 0x7fac51ab2590>




![output_13_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/4729a5d2-42aa-4e85-8172-2691b73b4baa)




```python
plt.imshow(img_gray, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fac51a22550>




![output_14_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/847845f7-c0ab-4d3b-ad26-6ea196c3c0f9)




```python
# Gray scale images have less information in them 
```


```python
# We can also resize our images 
fix_img.shape 
```




    (3888, 5184, 3)




```python
# Those are the number of pixels 
```


```python
new_img = cv2.resize(fix_img, (1000,400))
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7fac51a14290>




![output_18_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/9c233204-3007-419b-9450-95303c2fc7e9)




```python
# Now we have stretched our image out 
new_img.shape
```




    (400, 1000, 3)




```python
w_ratio = 0.5
h_ratio = 0.5

new_img = cv2.resize(fix_img, (0,0), fix_img, w_ratio, h_ratio)
```


```python
# Now the scale will be different
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7fac51973210>




![output_21_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/a981a910-719c-44d5-8b63-9f5dbff336fd)




```python
new_img.shape
```




    (1944, 2592, 3)




```python
# The 0.5 ratio cuts our image in half
```


```python
# We can also flip images 
flip_img = cv2.flip(fix_img, 0)
plt.imshow(flip_img)
```




    <matplotlib.image.AxesImage at 0x7fac518e2290>




![output_24_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/1d3f1148-49db-4b94-a8cc-aa584023ed61)




```python
# Now the image has been flipped vertically
```


```python
flip_img2 = cv2.flip(fix_img, -1)
plt.imshow(flip_img)
```




    <matplotlib.image.AxesImage at 0x7fac518d2090>




![output_26_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/9c79d0b2-815b-45b3-81f3-c4df4877b552)




```python
type(fix_img)
```




    numpy.ndarray




```python
cv2.imwrite('Puppy_fixed_image.jpg', flip_img)
```




    True




```python
# True means that it is right 
```


```python
# It does not save the color convert, only the image manipulation 
```


```python
import cv2
import matplotlib.pyplot as plt
%matplotlib inline 
```


```python
img = cv2.imread("Puppy.jpg")
```


```python
plt.imshow(img)
# The image will be blue again
```




    <matplotlib.image.AxesImage at 0x7fac51843090>




![output_33_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/9fd087a0-dbca-4a84-9ee4-092886e5b6f3)




```python
# Now we have to convert the colors 
img1 = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7fac4ddfae50>




![output_35_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/6cab6542-694d-47d1-8604-dd5490dc19ba)




```python
# Now we have our colors fixed 
img2 = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
# HSV formats are a bit older 
```


```python
plt.imshow(img2)
# HSV stands for Hue, Saturation, and Value
```




    <matplotlib.image.AxesImage at 0x7fac4dde7b10>




![output_37_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/79480b2b-7cf6-4541-99a2-ae060f0fd3e1)




```python
img3 = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)
```


```python
plt.imshow(img3)
```




    <matplotlib.image.AxesImage at 0x7fac4dd27f10>




![output_39_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/fea78973-9c84-489e-933a-83c1780c7839)




```python
# This is the HLS image 
```


```python
# Now we are using another image 
img1 = cv2.imread('DNC.jpg')
img2 = cv2.imread("Puppy.jpg")
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7fac4dc34d90>




![output_42_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/95047b70-ca10-4c5d-b5ed-57008206aae8)




```python
img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7fac4dc16e90>




![output_44_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/29db5bcd-b746-48f3-9224-93472f45c1a5)




```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7fac4db7af90>




![output_45_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/fbae1e5e-832a-4050-bc58-0e36a83ab01c)




```python
# Now we are resizing the images
img1 = cv2.resize(img1,(1200,1200))
img2 = cv2.resize(img2,(1200, 1200))
```


```python
alpha = 0.5
beta = 0.5
```


```python
blended = cv2.addWeighted(img1, alpha, img2, beta, gamma=0)
```


```python
plt.imshow(blended)
```




    <matplotlib.image.AxesImage at 0x7fac4a1a8e50>




![output_49_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/272d5f62-349e-4390-b047-ff5c23a443ec)




```python
# We have over ran the 2 images 
# This is very important for microscopy and fluoresence image overlay 
```


```python
alpha = 0.8
beta = 0.2

blended1 = cv2.addWeighted(img1, alpha, img2, beta, 0)
plt.imshow(blended1)
```




    <matplotlib.image.AxesImage at 0x7fac4a1174d0>




![output_51_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/51c76143-a8e2-4cb1-8d0a-5e585846880e)




```python
# Now the Do Not Copy is much darker and the puppies are much lighter 
```


```python
# We can also flip the alpha and beta values to reverse which is darker 
```


```python
# We can try blending with different image sizes 
img1 = cv2.imread('DNC.jpg')
img2 = cv2.imread('Puppy.jpg')

img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)

img1 = cv2.resize(img1, (200,200))
```


```python
large_img = img2
small_img = img1

x_offset = 0
y_offset = 0

x_end = x_offset + small_img.shape[1]
y_end = y_offset + small_img.shape[0]

large_img[y_offset:y_end, x_offset:x_end] = small_img

plt.imshow(large_img)
```




    <matplotlib.image.AxesImage at 0x7fac4a08b090>




![output_55_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/80becdcc-8767-4292-8313-ad2cd601e3d6)




```python
# https://github.com/worklifesg/Python-for-Computer-Vision-with-OpenCV-and-Deep-Learning
```


```python
# We will be working with new images
```


```python
import cv2
import matplotlib.pyplot as plt
%matplotlib inline 
```


```python
img = cv2.imread('rainbow.jpg')
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fac49ff3550>




![output_60_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/e22d94f1-e3d3-4f4b-9efe-d707619d1bb9)




```python
# Now we are going to do some thresholding 
img = cv2.imread('rainbow.jpg', 0)
```


```python
plt.imshow(img, cmap = 'gray')
# This will make the rainbow gray 
```




    <matplotlib.image.AxesImage at 0x7fac49fdb190>




![output_62_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/48b010a5-da43-462d-87c3-a85d1c6113df)




```python
ret1, thresh1 = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY)
```


```python
ret1
```




    127.0




```python
plt.imshow(thresh1, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fac49f35e10>




![output_65_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/a0fa9683-3fda-4466-b8b3-48737379939a)




```python
# Everything lighter turns white and anything darker turns black bc of our threshold
```


```python
# We can also inverse the limits 
img2 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img2, 127, 255, cv2.THRESH_TRUNC)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fac49f19d90>




![output_67_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/b6d5bc60-4391-47a8-9d2e-0f7f008a40c8)




```python
img3 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img3, 127, 255, cv2.THRESH_TOZERO)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fac49e7c8d0>




![output_68_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/81de9b78-6ec5-4290-942e-97c7bea8a476)




```python
img_r = cv2.imread('crossword.jpg', 0)
plt.imshow(img_r, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fac49e60150>




![output_69_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/3d6cb9f9-b0fc-4b0c-a366-7c8208f0265c)




```python
def show_pic(img):
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
show_pic(img_r)
```


![output_71_0](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/7b9e755d-471a-44e4-af1e-a8223d960e24)




```python
ret, th1 = cv2.threshold(img_r, 127, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![output_72_0](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/31f7d9e5-269b-4f45-ae96-59fa454ed07e)




```python
# We are turning the picture to just black on white
```


```python
ret, th1 = cv2.threshold(img_r, 200, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![output_74_0](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/85b0aceb-ae4b-4265-813c-ba17947838d4)




```python
th2 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)
```


```python
show_pic(th2)
```


![output_76_0](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/d696bd95-d68c-45b1-8991-a1044d7fbcb9)




```python
# We can blend this to do different thresholds 
blended = cv2.addWeighted(src1 = th1, alpha = 0.6,
                          src2 = th2, beta = 0.4, gamma = 0)

show_pic(blended)
```


![output_77_0](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/16a1af93-4461-4920-b157-962d62d85e03)




```python
# We layered th2 under th1 
```


```python
th3 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)

blended = cv2.addWeighted(src1 = th1, alpha = 0.6,
                          src2 = th2, beta = 0.4, gamma=0)

show_pic(blended)
```


![output_79_0](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/a8e9b5b6-0356-4c54-b033-4fc2cde93146)




```python
# Essentially we have removed background
```


## Corner Detection 

In this section, we analyzed methods for determining corners in images. 




```python
# We are using computer vision to detect features 
```


```python
import cv2
import numpy as np
import matplotlib.pyplot as plt 
%matplotlib inline 
```


```python
# This is for the green chessboard
flat_chess = cv2.imread('Chessboard_green.png')
flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2RGB)
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fed243884d0>




![output_2_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/f5afad89-8ab1-4b9e-a364-4043d95f0df3)




```python
# Now we will turn the green chessboard to gray 
gray_flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_flat_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fed242d0350>




![output_3_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/8fb45e31-0828-4374-b68b-6aeeb2631e01)




```python
# This gray scale chessboard will be our template 
```


```python
# Now we will plot our real chessboard
```


```python
real_chess = cv2.imread("Chessboard.jpg")
real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fed1fb559d0>




![output_7_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/cdb72d58-f9b7-42d7-a89e-263cb790b645)




```python
# Now we will convert the real chessbpard to grayscale
```


```python
gray_real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_real_chess, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fed1fb43ed0>




![output_9_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/1ba7fcb7-17cc-4e3b-85a3-6d3894c0742e)




```python
gray = np.float32(gray_flat_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)

dst = cv2.dilate(dst, None)
```


```python
# Above is the corner Harris method
```


```python
flat_chess[dst>0.01*dst.max()] = [255,0,0]

plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fed1fab8410>




![output_12_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/7b59d665-5af9-4dfe-85a4-349a18490f28)




```python
# Now anytime a corner is detected, it will become red
```


```python
# Now we will use the corner detection on a real chessboard
```


```python
gray = np.float32(gray_real_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k=0.04)
dst = cv2.dilate(dst, None)

real_chess[dst>0.01*dst.max()] = [255, 0, 0]

plt.imshow(real_chess)

```




    <matplotlib.image.AxesImage at 0x7fed1fa19c90>




![output_15_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/9c025dd3-789b-4123-96e4-1e0b7f6a2117)




```python
# Now we will use the Shi-Tomasi Corner detection method 
```


```python
corners = cv2.goodFeaturesToTrack(gray_flat_chess, 64, 0.01, 10)

```


```python
corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(flat_chess, (x,y),3,(255,0,0), -1)
    
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fed1f985e10>




![output_18_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/ef3aded4-4d45-4fd9-a6d9-fcc2f6c56196)




```python
# 49 corners have been detected 
```


```python
# Now we will do the Shi-Tomasi method on the real chessboard

corners = cv2.goodFeaturesToTrack(gray_real_chess, 100, 0.01, 10)

corners = np.int0(corners)

for i in corners:
    x, y = i.ravel()
    cv2.circle(real_chess, (x,y), 3, (0,255,0), -1)
    
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fed1f979210>




![output_20_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/54f2e1f6-7425-4ff5-b4e0-54475bb0f6d0)




```python
# This is an overlap of the two chessboards 
```

## Edge Detection 

In this analysis, we looked at ways to define and detect edges in biological data. 


```python
import cv2
```


```python
import numpy as np
```


```python
import matplotlib.pyplot as plt
%matplotlib inline 
```


```python
# I will be using the puppy image
img = cv2.imread("Puppy.jpg")
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fc7ccf70b10>




![output_3_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/5424811f-8da2-4728-b552-af33984fc786)




```python
# For edge detection, we do not need to convert the colors
```


```python
edges = cv2.Canny(image =img, threshold1 = 127, threshold2 = 127)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fc7cc63cf10>




![output_5_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/203054af-c5ae-4b5f-b27d-ad840548a3a2)




```python
med_value = np.median(img)
med_value
```




    107.0




```python
# 107 is our median value 
```


```python
# Now we are changing the threshold limits 

lower = int(max(0, 0.7*med_value))
upper = int(min(255,1.3*med_value))

edges = cv2.Canny(img, threshold1 = lower, threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fc7cc5bc7d0>




![output_8_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/f86dd182-d81e-49eb-ac2f-1c8de8be08b5)




```python
edges = cv2.Canny(image = img, threshold1 = lower, threshold2 = upper +100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fc7cc5310d0>




![output_9_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/f747d652-d94d-4952-b24f-2d2a03a6210f)




```python
# Now we will work with blurring 

blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fc7cc512fd0>




![output_10_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/f4a33786-090d-4a43-a9a3-f88628770af7)




```python
# We can also try increasing the kernel size; we are changing it from 5 to 7

blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fc7cc486d10>




![output_11_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/31b81b0f-2e47-46f2-bd85-48449d605401)




```python
# Now we are increasing the upper threshold limit by 50

blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 50)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fc7cc3f4ad0>




![output_12_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/45b12504-c00b-483a-bc58-47c26c44697a)




```python
# We are now increasing the upper limit by 100

blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fc7cc3e3810>




![output_13_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/8d882f1a-8d57-4bdc-b569-ccae8e43d1f4)




```python
# Now we are changing kernel size and the upper limit 

blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fc7cc357550>




![output_14_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/d1a0dad8-038d-434b-9bd7-3e7516a0d79b)




```python
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 75)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fc7cc2c62d0>




![output_15_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/7e2937b3-9419-4a73-be4b-aa4e4191b9c3)




```python
blurred_img = cv2.blur(img, ksize = (3,3))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fc7cc2aaf90>




![output_16_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/a725e108-3f25-421b-a3fa-015105eef118)




```python
blurred_img = cv2.blur(img, ksize = (2,2))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper - 100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fc7cc219d50>




![output_17_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/ba1ca37a-b009-4413-81de-af698222acdf)




```python
# You have to play around with kernel size and threshold values to find what works best for your image 
```


```python
# On my image smaller kernel sizes and subtracting from the upper limit gave me better results 
```


## Object Detection 

In this section, we used object detection via Open CV for phenotyping purposes. 



```python
import cv2
```


```python
import numpy as np
```


```python
import matplotlib.pyplot as plt 
```


```python
%matplotlib inline 
```


```python
# The training sunflower is one singular flower and the testing sunflower is a field of sunflowers
```


```python
# We are loading our training flower
full = cv2.imread('Training_sunflower.jpeg')
```


```python
# We are converting the color for training 

full = cv2.cvtColor(full, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(full)
# This will show the photo in the right colors
```




    <matplotlib.image.AxesImage at 0x7f13e1e18690>




![output_7_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/dffc1272-9476-448f-8afd-9e1e0db00a8c)




```python
# Now we will upload the testing sunflower image 

test = cv2.imread('Searching_sunflower.jpeg')
```


```python
test = cv2.cvtColor(test, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(test)
```




    <matplotlib.image.AxesImage at 0x7f13e15af410>




![output_10_1](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/fdc8f728-eecd-4bac-9637-2e9db46865b6)




```python
# This is our searching/testing image with multiple sunflowers 
```


```python
print("Test image shape:", full.shape)
print("Training image shape", test.shape)
```

    Test image shape: (1329, 1200, 3)
    Training image shape (413, 612, 3)



```python
methods = ['cv2.TM_CCOEFF', 'cv2.TM_CCOEFF_NORMED', 'cv2.TM_CCORR', 'cv2.TM_CCORR_NORMED', 'cv2.TM_SQDIFF', 'cv2.TM_SQDIFF_NORMED']
# These are all of the different methods we will try 
```


```python
for m in methods:
    
    test_copy = test.copy()
    method = eval(m)
    
    res = cv2.matchTemplate(test_copy, full, method)
    
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
    
    if method in [cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED]:
        top_left = min_loc
    else:
        top_left = max_loc 
        
    height, width, channels = full.shape
    bottom_right = (top_left[0] + width, top_left[1] + height)
    
    cv2.rectangle(test_copy, top_left, bottom_right, (255,0,0),10)
    
    plt.subplot(121)
    plt.imshow(res)
    plt.title("Heatmap of template matching")
    plt.subplot(122)
    plt.imshow(test_copy)
    plt.title('Detection of template')
    
    plt.suptitle(m)
    
    plt.show()
    print('\n')
    print('\n')
```


![output_14_0](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/b57b0850-6bb2-4d94-b299-42b5aa26fae1)



    
    
    
    



![output_14_2](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/dac5a763-9837-43cb-a4b4-a67334efef26)



    
    
    
    



![output_14_4](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/c700dcdf-4285-4f6f-a2dc-65654594be42)



    
    
    
    



![output_14_6](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/b7c13978-4551-4b03-935c-19a0ab8b77eb)



    
    
    
    



![output_14_8](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/f744f78b-2510-496b-b0a2-08008cb778d0)



    
    
    
    



![output_14_10](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/7a149c1f-aa2b-40e6-af40-e51f39d21776)



    
    
    
    



```python
# It is making a heatmap of the input image to detect the sunflower heads 
```


## Feature Matching 


Feature method is useful for sorting through data to find specific objects. 


```python
import cv2
import numpy as np 
import matplotlib.pyplot as plt
%matplotlib inline 
```


```python
def display(img, cmap = 'gray'):
    fig = plt.figure(figsize = (12,10))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
# We are going to be identifying cereals from a collage 
```


```python
apple_jacks = cv2.imread("Apple_jacks.jpeg", 0)
display(apple_jacks)

# This is for just the Apple Jacks box
```


![output_3_0](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/310db89f-f3f8-441d-9aac-e53d6dc1633e)




```python
cereals = cv2.imread('All_cereals.jpeg', 0)
display(cereals)

# This is for the cereal collage 
```


![output_4_0](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/06d7e3d3-e680-4625-acb8-fde97c3b0ff5)




```python
# The first method we will do is brute force matching 

orb = cv2.ORB_create()

kp1,des1 = orb.detectAndCompute(apple_jacks, mask=None)
kp2,des2 = orb.detectAndCompute(cereals, mask=None)
```


```python
# Now we will run the bruce force matching 

bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck = True)
matches = bf.match(des1, des2)
```


```python
matches = sorted(matches, key = lambda x:x.distance)
```


```python
apple_jacks_matches = cv2.drawMatches(apple_jacks, kp1, cereals, kp2, matches[:25], None, flags = 2)
```


```python
display(apple_jacks_matches)
```


![output_9_0](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/95f2fdb9-80bc-4b60-afe6-ada07b367a91)




```python
# This method is not the best for matching just apple jacks 
```


```python
sift = cv2.SIFT_create()
```


```python
kp1, des1 = sift.detectAndCompute(apple_jacks, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
# Now we will use brute force with the sift method 

bf = cv2.BFMatcher()
matches = bf.knnMatch(des1, des2, k=2)
```


```python
# Now we will apply a ratio test 

good = []

for match1, match2 in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
# The good matches abide by the 75% rule 
```


```python
print('Length of total matches:', len(matches))
print('Length of good matches:', len(good))
```

    Length of total matches: 2136
    Length of good matches: 39



```python
sift_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, good, None, flags =2)
display(sift_matches)
```


![output_17_0](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/3df90164-ad67-41f2-9b7d-5f6f104545da)




```python
# Now we will try flan base matching 
```


```python
sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(apple_jacks, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
flann_index_KDtree = 0
index_params = dict(algorithm=flann_index_KDtree, trees = 5)
search_params = dict(checks=50)
```


```python
# Now we can set parameters 

flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k=2)

good = []

for match1, match2, in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
flann_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, good, None, flags = 0)
display(flann_matches)
```


![output_22_0](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/0adf191e-f92f-4db4-b1e3-e6fb4288973e)




```python
# Now we are adding a mask to it 

sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(apple_jacks, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
flann_index_KDtree = 0
index_params = dict(algorithm= flann_index_KDtree, trees = 5)
search_param = dict(checks =50)
```


```python
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k = 2)
```


```python
matchesMask = [[0,0] for i in range(len(matches))]
```


```python
# 0,0 means pure black will be added 
```


```python
# Now we can reduce the distance between features to get better matches 

for i, (match1, match2) in enumerate(matches):
    if match1.distance <0.75*match2.distance:
        matchesMask[i] = [1,0]
        
draw_params = dict(matchColor = (0, 255, 0),
                   singlePointColor = (255,0,0),
                   matchesMask = matchesMask,
                   flags = 0) 
```


```python
flann_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, matches, None, **draw_params)

display(flann_matches)
```


![output_29_0](https://github.com/norabethhill/Python_2_Portfolio/assets/132499791/e1d0adfc-bd8b-448b-b17a-a0bba1463138)




```python
# now we are showing all unique features we identified on both the apple jacks box and on all the cereal boxes 
# we changed the line color so we can easily see the best box to match 
```










