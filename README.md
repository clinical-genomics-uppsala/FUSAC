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

| Flag | Name | Function | Required | Default | Alternative |
| --- | --- | --- | --- | --- | --- |
| -b | inputBAM | Input BAM file path | Yes | N/A | Any |
| -v | inputVCF | Input VCF file path | Yes | N/A | Any |
| -t | threads | No. threads to run the program | No | 1 | Any integer |
| -qs | queueSize | Threading deque size | No | Infinite | Any |
| -fb | ffpeBases | Bases used for FFPE classification | No | C:T>G:A | all |
| -up | umiPosition | Location of the UMI-tag in a read | No | Query-name (qrn) | Rx-tag (rx) |
| -sc | splitCharacter | Split character for the UMI-tag | No | + | Any |
| -cf | csvFile | Generate an output CSV file | No | yes | no |

#### Example 1
We wish classify all mismatches belonging to the file example_bam using the example_vcf file. The Reads in the example\_bam file have their UMI-tag stored in the query-name, which is separated by the character "_". The program is being run on a laptop with 4 cores, and we wish to limit the queue to 9 variant-records. 

```
python fumic.py -b example_bam.bam -v example_vcf.vcf -t 3 -qs 9 -b all -sc _
```

#### Example 2
We wish classify only C:T>G:A artefacts belonging to the file example\_bam using the example_vcf file. The Reads in the example\_bam file have their UMI-tag stored in the RX-tag, and are not separated by any character. The program is being run on a cluster with 16 cores. We do not wish to limit the queue, but rather have it infinite. Furthermore, we do not wish to generate a CSV file.

```
python fumic.py -b example_bam.bam -v example_vcf.vcf -t 15 -up rx -sc "" -cf no
```

### Interpreting FUMIC's output
The output from FUMIC can be easily extracted for analysis, or read directly to form a quick opinion. In the following example one modified variant-record is presented. The example given below is a finctional FFPE-classed variant record generated by FUMIC. The field INFO have been cut out for sake of clarity and are marked as "...".

| Field | Value |
| --- | --- |
| CHROM | chr1 |
| POS | 4367323 |
| ID | rs1490413 |
| REF | G |
| ALT | A |
| QUAL | . |
| FILTER | FFPE |
| INFO | ... |
| FORMAT | GT:AD:AF:DP:UMI:SUMI |
| GT:AD:AF:DP: | 0/1:571,632:0.527:1203: |
| UMI | 235;313;15;0;0;313;328;250;235;272;322;306;256: |
| SUMI | 0;563;0;0;0;563;563;0;0;570;571;4;7 |

In the example variant-record above, we see that the variant-caller initially identified the record as a mutation from G to A based on the "REF" and "ALT" fields. We then see in the "Filter" field that FUMIC have deemed this to be a record with molecular support for a FFPE artefact. 

Further investigation in the "FORMAT" field, more specifically the two fields UMI and SUMI added by FUMIC shows us that the molecular support for no mutation is 235, the support for a true mutation is 313 and the support for a FFPE-artefacts is 15. The molecular support for a deletion or an unknown is deemed to be 0. 

We can then further identify the molecular support for the reference nucleotide G on string 1 for paired reads as 313, and on string 2 328. The support for the variant nucleotide A on string 1 for paired reads is deemed to be 250, and on string 2 235. We can quickly see that this score is consistent with the 15 FFPE-artefact classifications that were found by FUMIC. For paired reads without their mate the molecular support for the reference nucleotide is  272 on string 1 and 322 on string 2. And the support for the variant nucleotide is 306 on string 1, and 256 on string 2.

To summarize, we can see that for position 4367323 15 molecules support an FFPE artefact. However, as 313 molecules support a true mutation, this specific position is unlikely to be a true FFPE-artefact, and thus these results are more likely to be caused by some other factor. 

## Reference manual
The reference manual covers all functions belonging to FUMIC, describing their purpose, mehodology and input/output. Use to get a better understanding of FUMIC or if you have any questions. 

### Main
After the user have supplemented their desired input arguments using the flags covered in the section "Quickstart", the main function will use the vcf-file as an input argument to populate a newly generated double-ended-queue (deque) of size -qs, with the variant-calls found within the VCF-file. Based on the flag "threads", the function will then create -t consumer threads, their results appended to a result-queue. Once the VCF has been iterated through entirely, the output will be written to an output vcf-file.

### QueueThread
The QueueThread function is a producer which takes the variant-calls found within the VCF-file and populates a deque (while not full) to be used by the consumer function.

| Input | Function |
| --- | --- |
| vcf_file | VCF filehandle |
| thr_que | Deque to be populated |
| ID | rs1490413 |

### ResultThread
While the deque is not empty, the ResultThread function calls upon the vcf_extract function using the variant-record extracted as input. The subsequent results are then stored in a separate thread while not None.

| Input | Function |
| --- | --- |
| bam_path | Path to BAM-file |
| thr_que | Populated deque to be used as input for vcf_extract |
| res_que | Queue to be populated with the results from vcf_extract |
| ffpe_b | Optional input argument controlling which mismatches to consider  for FFPE-classification |
| ext_fun | Function for extracting the UMI-tag from a read |
| spl_fun | Function used for splitting the UMI-tag in a read |
| spl_cha | Character used for splitting the UMI-tag |

### vcf_extract
The vcf_extract function uses the supplemented variant-record to extract all reads in the BAM-file overlapping with its position. This newly generated list is used for the var_extract function to return molecular data. The output from var_extract is then subsequently used in the inf_builder function. Finally, the output from inf_builder is added to the copied input record and returned.

| Input | Function |
| --- | --- |
| record | Variant-record of interest
| bam_file | BAM-file filehandle |
| ffpe_b | Optional input argument controlling which mismatches to consider  for FFPE-classification |
| ext_fun | Function for extracting the UMI-tag from a read |
| spl_fun | Function used for splitting the UMI-tag in a read |
| spl_cha | Character used for splitting the UMI-tag |

### var_extract
he var\_extract function first calls the umi_maker function which creates a dict based on the directionality and umi-tags of the supplemented reads in the bam_lst. This dict is then used to call the pos_hits function which will return a dict of consensus nucleotides for each UMI. For paired reads (ie: exists on both string 1 and string 2 for a UMI) the output from the pos_hits function is then used to call the ffpe_finder function which returns a dict containing the variant type for eachUMI along with the molecular support for each variant type. The pos_checker function returns a dict with data regarding the support for each variant type on the variant-record position, as well as the nucleotides present for each variant with respect to their UMI.

| Input | Function |
| --- | --- |
| record | Variant-record of interest
| bam_lst |Input list of BAM-reads aligning to the variant call |
| ffpe_b | Optional input argument controlling which mismatches to consider  for FFPE-classification |
| rec_pos | The position of the variant in the reference genome |
| var_nuc | The nucleotide called in the variant-record |
| ref_nuc | The nucleotide found in the reference genome at the variant-call position |
| ffpe_b | Optional input argument controlling which mismatches to consider  for FFPE-classification |
| ext_fun | Function for extracting the UMI-tag from a read |
| spl_fun | Function used for splitting the UMI-tag in a read |
| spl_cha | Character used for splitting the UMI-tag |

| Dict | Structure |
| --- | --- |
| Input | Example list: bam_lst = [read_1, read_2 ... ] |
| Output | Example dict: mpd_res[umi_key] = {"Single_Hits": Str1_Hits: {}, Str2_Hits:{C,T}, "Mate_Hits": Mutation_Hits": {}, "FFPE_Hits": {"String_1": C, "String_2": T}, "N_Hits": {}, "Del_Hits": {}, "Reference_Support": 0, "Mutation_Support": 0, \\ "FFPE_Support": 1, "N_Support": 0, "Del_Support": 0 |

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


