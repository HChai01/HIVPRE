# Feature encoding
This is the Java implementation for encoding some genome-, proteome- and gene ontology-based features that require data processing.

- Here is the list and description:
  For genome-based features:
    Nucleotide_composition.java -- calculate 23 nucleotide composition in the coding sequence.
    Codon_usage_CDS.java -- calculate 64 codon usage in the coding sequence.
    Codon_usage_mRNA.java -- calculate 64 codon usage in the mRNA.
  For proteome-based features:
  
  For gene ontology-based features:
    GO_root_biological_process.java -- map the Gene Ontology term to 29 child term of biological process (GO:0008150) through the derivation tree.
    GO_root_molecular_function.java -- map the Gene Ontology term to 16 child term of molecular function (GO:0003674) through the derivation tree.
    GO_root_cellular_component.java -- map the Gene Ontology term to 21 child term of cellular component (GO:0005575) through the derivation tree.

