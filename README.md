[![Build Status](https://travis-ci.org/clinical-genomics-uppsala/fumic.svg?branch=master)](https://travis-ci.org/clinical-genomics-uppsala/fumic)

# FUMIC 
FFPE-artefact UMI-based Mapper for Imputation in Cancer-sample tissue data

## Getting Started
FUMIC is a python-based program for the identification and classification of FFPE-artefacts in UMI-tagged sequence data. Using a VCF and BAM-file + BAI-file as input, FUMIC is able to successfully identify group and collapse all reads aligning to a position called by the VCF. From this FUMIC generates consensus sequences for each UMI as well as identifying their string of origin before amplification. These consensus sequences are then compared with their mate, and thus FUMIC is able to not only identify C:G$>$T:A artefacts left by hydrolytic deamination, but also identify true mutations, deletions, unknowns and any other type of mismatch. 

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

### Quickstart
Required input arguments for running FUMIC are -b and -v,  which are the respective paths to the .bam and .vcf file. Furthermore, an indexed BAM (.bai) file is required for extracting desired segments of the BAM-file. The other input flags are not required, but should be changed if the default value is not representative of the desired output. To minimize run-time and CPU-load FUMIC can run on multiple threads. Unfortunately, as pickling cannot deal with open filehandles, multiprocessing is not a viable option as this would require the file to be opened for every read aligning to the variant-position. Instead, FUMIC uses the python "threading" module with a producer-consumer approach, where the producer generates and populates a queue, and the consumer thread extracts the inhabitants of this queue for analysis. To control this threading process, the arguments threads (-t) and queueSize (-qs) determine the number of threads to be run and the size of the threading queue respectively.
The default values for threads and queueSize respectively are one active thread and an infinite queue, but can be set to any integer value desired. 

The default FFPE-classification mode focuses solely on C:G$>$T:A artefacts, however if desired the program can also identify any mismatching consensus nucleotides using the input flag ffpeBases (-fb) with the option "all". Lastly, FUMIC is entirely dependent on the UMI-tag being properly extracted to ensure that reads are assigned to String 1 or String 2 as origin. Therefore, the user can specify through the umiPosition (-up) tag if the UMI-tag is located in the query-name ("qrn") or the RX-tag respectively ("rx"). Furthermore, the UMI-tag needs to be split in half to be rearranged correctly, which can be done using the input splitCharacter (-sc) which represents the character on which to split the tag. For reads where the UMI-tag is not separated by a tag, the input "" should be used to split the tag in half. 

The final input to consider is csvFile (-cf) which controls whether or not FUMIC generates an output CSV file based on the fumic-output. This CSV generates a separate row for each variant-record with columns for the molecular support for the reference genome nucletoide, the variant-call nucleotide, the number of FFPE-calls, the overall frequency of FFPE-artefacts for each variant-record, and the type of mismatch for the variant-record. The default setting is to generate the CSV, but if this is not required the function can be turned off using the input  "no".

```
| Flag  | Name | Function | Required |  Default | Alternative|
| --- | --- | --- | --- | --- | --- |
| -b | inputBAM | Input BAM file path | Yes | N/A | Any |
| -v | inputVCF | Input VCF file path | Yes | N/A | Any |
| --- | --- | --- | --- | --- | --- |
| --- | --- | --- | --- | --- | --- |
| --- | --- | --- | --- | --- | --- |
```

## Functions

### 

### Tests
FUMIC comes supplemented with a test_fumic.py function which contains unit-tests for FUMIC. These can be easily run through Travis, or manually through the terminal.

```
python test_fumic.py
```

More specifically, the python test_fumic.py function contain unit-tests for most functions within FUMIC. These are as follows:

```
test_qrn_ext
tests the qrn_ext function for three variations of a UMI-tag stored in the query name
```

```
test_cha_splt
tests the cha_splt function for a case with a + character, a _ character, a - character and a X character
```

```
test_umi_maker
tests the umi_maker function for a case with a forward read1 read, a forward read 2 read, a reverse read 1 read and a reverse read 2 read. 
```

```
test_pos_hits
tests the pos_hits function for a case with no mutation, a case with a mutation, a case with a ffpe-artefact, a case with an unknown, a a case with a deletion, and finally for a case with a mismatched pair of consensus nucleotides
```

```
test_ffpe_finder
tests the ffpe_finder function for a case with no mutation, a case with a mutation, a case with a ffpe-artefact, a case with an unknown, and a case with a deletion
```

```
test_ffpe_finder
tests the var_extract function for a case with no mutation, a case with a mutation, a case with a ffpe-artefact, a case with an unknown, and a case with a deletion
```


