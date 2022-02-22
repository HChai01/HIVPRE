# Feature encoding
This is the Java implementation for encoding some genome-, proteome- and gene ontology-based features that require data processing.

Here is the list and description of the compiled Java programs:

For genome-based features:
- Nucleotide_composition.java -- calculate 23 nucleotide composition in the coding sequence.
- Codon_usage_CDS.java -- calculate 64 codon usage in the coding sequence.
- Codon_usage_mRNA.java -- calculate 64 codon usage in the mRNA.

For proteome-based features:
- Amino_acid_composition.java -- calculate 37 amino acid composition in the coding sequence.

For gene ontology-based features:
- GO_root_biological_process.java -- map the Gene Ontology term to 29 child term of [biological process](https://www.ebi.ac.uk/QuickGO/term/GO:0008150) through the derivation tree.
- GO_root_molecular_function.java -- map the Gene Ontology term to 16 child term of [molecular function](https://www.ebi.ac.uk/QuickGO/term/GO:0003674) through the derivation tree.
- GO_root_cellular_component.java -- map the Gene Ontology term to 21 child term of [cellular component](https://www.ebi.ac.uk/QuickGO/term/GO:0005575) through the derivation tree.

# Instruction
Our codes require support from [Java Platform and Standard Edition Development Kit (JDK)](https://www.oracle.com/java/technologies/downloads/#jdk17-linux). 

They are well supported in Linux or MAX OS.

There are only five simple steps to use our codes:
- Step 1: install a JDK environment on your computer;
- Step 2: download 'Compiled_code.zip' and uncompress it;
- Step 3: add your nucleotide, protein sequences or gene ontology profiles in the example files (e.g., /compiled_code/data/Protein_example.txt);
- Step 4: ues 'Terminal' to the directory of 'Compiled_code';
- Step 5: use 'java [program_name]' to encode features for your sequences or gene ontology profiles (e.g., `java Nucleotide_composition`).

# Input and output examples
The example input files are placed at /HIVPRE/data/ while outputs for them can be found at /HIVPRE/Example_output/.

Here is the list and description for the input examples:
- DNA_CDs_example.txt -- 3 coding sequences in [FASTA format](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp).
- mRNA_example.txt -- 3 mRNA sequences in FASTA format.
- Protein_example.txt -- 3 protein sequences in FASTA format.
- Gene_Ontology_example.txt -- 6 genes ontology annotaion data. 

# Web server
Our web server for the HIVPRE project is accessible at: http://hivpre.cvr.gla.ac.uk/.

# Citation
Chai H, Gu Q, Hughes J, Robertson DL (2022) In silico prediction of HIV-1-host molecular interactions and their directionality. PLoS Comput Biol 18(2): e1009720. https://doi.org/10.1371/journal.pcbi.1009720


