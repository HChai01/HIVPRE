# Feature encoding
This is the Java implementation for encoding some genome-, proteome- and gene ontology-based features that require data processing.

Here is the list and description:
For genome-based features:
- Nucleotide_composition.class -- calculate 23 nucleotide composition in the coding sequence.
- Codon_usage_CDS.class -- calculate 64 codon usage in the coding sequence.
- Codon_usage_mRNA.class -- calculate 64 codon usage in the mRNA.

For proteome-based features:
- Amino_acid_composition.class -- calculate 37 amino acid composition in the coding sequence.

For gene ontology-based features:
- GO_root_biological_process.class -- map the Gene Ontology term to 29 child term of [biological process](https://www.ebi.ac.uk/QuickGO/term/GO:0008150) through the derivation tree.
- GO_root_molecular_function.class -- map the Gene Ontology term to 16 child term of [molecular function](https://www.ebi.ac.uk/QuickGO/term/GO:0003674) through the derivation tree.
- GO_root_cellular_component.class -- map the Gene Ontology term to 21 child term of [cellular component](https://www.ebi.ac.uk/QuickGO/term/GO:0005575) through the derivation tree.

# Environment and settings
Our codes require support from [Java Platform and Standard Edition Development Kit (JDK)](https://www.oracle.com/java/technologies/downloads/#jdk17-linux). 

We recommend to install [Eclipse IDE](https://www.eclipse.org/downloads/packages/release/kepler/sr2/eclipse-ide-java-ee-developers) to run our codes.


# Input and output examples
The example input files are placed at /HIVPRE/data/ while outputs for them can be found at /HIVPRE/data/Example_output/.

Here is the list and description for the input examples:
- DNA_CDs_example.txt -- 3 coding sequences in [FASTA format](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp).
- mRNA_example.txt -- 3 mRNA sequences in FASTA format.
- Protein_example.txt -- 3 protein sequences in FASTA format.
- Gene_Ontology_example.txt -- 6 genes ontology annotaion data. 

# Web server
Our web server for the HIVPRE project is accessible at: http://hivpre.cvr.gla.ac.uk/.

# Citation
Haiting Chai, Quan Gu, Joseph Hughes and David L. Robertson (202X), In silico prediction of HIV-1-host molecular interactions and their directionality, PLOS Computational Biology, XX(X): XXX-XXX, PMID: XXXXXXXX.


