[![Build Status](https://travis-ci.org/clinical-genomics-uppsala/fumic.svg?branch=master)](https://travis-ci.org/clinical-genomics-uppsala/fumic)

# fumic

FUMIC is a python-based program for the identification and classification of FFPE-artefacts in UMI-tagged sequence data. Using a VCF and BAM-file as input, FUMIC is able to successfully identify group and collapse all reads aligning to a position called by the VCF, generating consensus sequences for each UMI as well as identifying their string of origin before amplification. These consensus sequences are then compared with their mate, and thus FUMIC is able to not only identify C:G$>$T:A artefacts left by hydrolytic deamination, but also identify true mutations, deletions, unknowns and any other type of mismatch. 

FUMIC requires the user to have basic understanding of the data they wish to study, namely the location of the UMI-tag and how it is structurally stored within the  BAM-file. More specifically it requires the user to define both the location of the tag (query name or the RX-field) as well as the character the UMI is separated by if such a character exists. 

From this input, FUMIC generates a modified VCF-file as output. The output VCF is a copy  of the input VCF but has a modified "FILTER" field where any classified FFPE-artefact will display  "FFPE". Furthermore, the output VCF with also have a modified "FORMAT" field where the molecular support for the variant position having no mutation a true mutation, an FFPE-artefact, an unknown, or a deletion will be displayed. This field also contains the molecular support for the reference genome nucleotide as well as the called variant nucleotide for paired reads on str1, str2, as well as the support on single reads belonging to string 1 and string 2.
