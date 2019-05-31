[![Build Status](https://travis-ci.org/clinical-genomics-uppsala/fumic.svg?branch=master)](https://travis-ci.org/clinical-genomics-uppsala/fumic)

# FUMIC 
FFPE-artefact UMI-based Mapper for Imputation in Cancer-sample tissue data

## Getting Started
FUMIC is a python-based program for the identification and classification of FFPE-artefacts in UMI-tagged sequence data. Using a VCF and BAM-file as input, FUMIC is able to successfully identify group and collapse all reads aligning to a position called by the VCF, generating consensus sequences for each UMI as well as identifying their string of origin before amplification. These consensus sequences are then compared with their mate, and thus FUMIC is able to not only identify C:G$>$T:A artefacts left by hydrolytic deamination, but also identify true mutations, deletions, unknowns and any other type of mismatch. 

FUMIC requires the user to have basic understanding of the data they wish to study, namely the location of the UMI-tag and how it is structurally stored within the  BAM-file. More specifically it requires the user to define both the location of the tag (query name or the RX-field) as well as the character the UMI is separated by if such a character exists. 

From this input, FUMIC generates a modified VCF-file as output. The output VCF is a copy  of the input VCF but has a modified "FILTER" field where any classified FFPE-artefact will display  "FFPE". Furthermore, the output VCF with also have a modified "FORMAT" field where the molecular support for the variant position having no mutation a true mutation, an FFPE-artefact, an unknown, or a deletion will be displayed. This field also contains the molecular support for the reference genome nucleotide as well as the called variant nucleotide for paired reads on str1, str2, as well as the support on single reads belonging to string 1 and string 2.

### Prerequisites
FUMIC is based on the python module Pysam, and thus requires this to be installed. Furthermore, FUMIC requires the module pandas to be run. Both modules can be obtained for free through their respective githubs, or easily installed through pip.

```
sudo pip install pandas
sudo pip install pysam
```

There are four fundamental assumptions made by FUMIC to yield the desired output 
(1) The input sequence data is UMI-tagged 
(2) The sequence data alignment includes gapped bases 
(3) The UMI-tag is located either in the query-name or RX-tag fields in the BAM (translated to SAM) file. 
(4) If the UMI-tag is located in the query-name, it is located as the last entry and separated from the rest of the query name through  

The first assumption is necessary due to FUMIC's algorithm working in a classifying manner. To identify all reads stemming from a source molecule, a common identifier in the form of the UMI is vital for collapsing reads into a consensus sequence. The second assumption is necessary to properly locate the called variant within each subsequent read belonging to a UMI of interest. FUMIC uses the reference genome to identify the correct position for each subsequent read. And thus, If gapped bases are not included, this position will be incorrect, thus yielding an incorrect comparison. The third and fourth assumption are both necessary to ensure that the UMI-tagged data can be properly extracted from each read. 

### Running the tests
FUMIC comes supplemented with a test_fumic.py function which contains unit-tests for FUMIC. These can be easily run through Travis, or manually through the terminal.

```
python test_fumic.py
```
