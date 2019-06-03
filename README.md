[![Build Status](https://travis-ci.org/clinical-genomics-uppsala/fumic.svg?branch=master)](https://travis-ci.org/clinical-genomics-uppsala/fumic)

# FUMIC 
FFPE-artefact UMI-based Mapper for Imputation in Cancer-sample tissue data

## Getting Started
FUMIC is a python-based program for the identification and classification of FFPE-artefacts in UMI-tagged sequence data. Using a VCF and BAM-file + BAI-file as input, FUMIC is able to successfully identify group and collapse all reads aligning to a position called by the VCF. From this FUMIC generates consensus sequences for each UMI as well as identifying their string of origin before amplification. These consensus sequences are then compared with their mate, and thus FUMIC is able to not only identify C:G>T:A artefacts left by hydrolytic deamination, but also identify true mutations, deletions, unknowns and any other type of mismatch. 

FUMIC requires the user to have basic understanding of the data they wish to study, namely the location of the UMI-tag and how it is structurally stored within the  BAM-file. More specifically it requires the user to define both the location of the tag (query name or the RX-field) as well as the character the UMI is separated by if such a character exists. 

From this input, FUMIC generates a modified VCF-file as output. The output VCF is a copy  of the input VCF but has a modified "FILTER" field where any classified FFPE-artefact will display  "FFPE". Furthermore, the output VCF with also have a modified "FORMAT" field where the molecular support for the variant position having no mutation a true mutation, an FFPE-artefact, an unknown, or a deletion will be displayed. This field also contains the molecular support for the reference genome nucleotide as well as the called variant nucleotide for paired reads on str1, str2, as well as the support on single reads belonging to string 1 and string 2.

### Prerequisites
FUMIC is based on the python module Pysam, and thus requires this to be installed. Furthermore, FUMIC requires the module pandas to be run. Both modules can be obtained for free through their respective github-pages, or easily installed through pip.

```
sudo pip install pandas
sudo pip install pysam
```

There are four fundamental assumptions made by FUMIC to yield the desired output 
  1. The input sequence data is UMI-tagged 
  2. The sequence data alignment includes gapped bases 
  3. The UMI-tag is located either in the query-name or RX-tag fields in the BAM (translated to SAM) file. 
  4. If the UMI-tag is located in the query-name, it is located as the last entry and separated from the rest of the query name through  

The first assumption is necessary due to FUMIC's algorithm working in a classifying manner. To identify all reads stemming from a source molecule, a common identifier in the form of the UMI is vital for collapsing reads into a consensus sequence. The second assumption is necessary to properly locate the called variant within each subsequent read belonging to a UMI of interest. FUMIC uses the reference genome to identify the correct position for each subsequent read. And thus, If gapped bases are not included, this position will be incorrect, thus yielding an incorrect comparison. The third and fourth assumption are both necessary to ensure that the UMI-tagged data can be properly extracted from each read. 

### Quickstart
Required input arguments for running FUMIC are -b and -v,  which are the respective paths to the .bam and .vcf file. Furthermore, an indexed BAM (.bai) file is required for extracting desired segments of the BAM-file. The other input flags are not required, but should be changed if the default value is not representative of the desired output. To minimize run-time and CPU-load FUMIC can run on multiple threads. Unfortunately, as pickling cannot deal with open filehandles, multiprocessing is not a viable option as this would require the file to be opened for every read aligning to the variant-position. Instead, FUMIC uses the python "threading" module with a producer-consumer approach, where the producer generates and populates a queue, and the consumer thread extracts the inhabitants of this queue for analysis. To control this threading process, the arguments threads (-t) and queueSize (-qs) determine the number of threads to be run and the size of the threading queue respectively.
The default values for threads and queueSize respectively are one active thread and an infinite queue, but can be set to any integer value desired. 

The default FFPE-classification mode focuses solely on C:G>T:A artefacts, however if desired the program can also identify any mismatching consensus nucleotides using the input flag ffpeBases (-fb) with the option "all". Lastly, FUMIC is entirely dependent on the UMI-tag being properly extracted to ensure that reads are assigned to String 1 or String 2 as origin. Therefore, the user can specify through the umiPosition (-up) tag if the UMI-tag is located in the query-name ("qrn") or the RX-tag respectively ("rx"). Furthermore, the UMI-tag needs to be split in half to be rearranged correctly, which can be done using the input splitCharacter (-sc) which represents the character on which to split the tag. For reads where the UMI-tag is not separated by a tag, the input "" should be used to split the tag in half. 

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
The var_extract function first calls the umi_maker function which creates a dict based on the directionality and umi-tags of the supplemented reads in the bam_lst. This dict is then used to call the pos_hits function which will return a dict of consensus nucleotides for each UMI. For paired reads (ie: exists on both string 1 and string 2 for a UMI) the output from the pos_hits function is then used to call the ffpe_finder function which returns a dict containing the variant type for eachUMI along with the molecular support for each variant type. The pos_checker function returns a dict with data regarding the support for each variant type on the variant-record position, as well as the nucleotides present for each variant with respect to their UMI.

| Input | Function |
| --- | --- |
| bam_lst |Input list of BAM-reads aligning to the variant call |
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
| Output | Example dict: mpd_res[umi_key] = {"Single_Hits": Str1_Hits: {}, Str2_Hits:{C,T}, "Mate_Hits": Mutation_Hits": {}, "FFPE_Hits": {"String_1": C, "String_2": T}, "N_Hits": {}, "Del_Hits": {}, "Reference_Support": 0, "Mutation_Support": 0, "FFPE_Support": 1, "N_Support": 0, "Del_Support": 0 |

### inf_builder
The inf_builder function uses the output from var_extract to generate a list containing strings representing the data found for each record, more specifically support for each variant-type, as well as the support for the reference and variant call for string 1 and string 2.

The inp_dict is meant to be the output from the var_extract function,and is divided into two dicts named Single Hits and Mate Hits. The Single-dict contains the molecular support for the reference genome nucleotide and the variant nucleotide based on all reads without a mate. The Mate-Hits dicts instead contains data regarding the variant-classification, the support for each variant type,and the molecular support for the reference genome nucleotide and the variant nucleotide based on reads with a mate.

Returns a list containing the support for each variant-type, as well as the support for the reference and variant call for string 1 and string 2

| Input | Function |
| --- | --- |
| read | Read of interest
| inp_dict |Input dict dict for mapped and unmapped reads. Each of these dicts containing a single-hits and a mate-hits dict. The mate-hits dict in turn contains data regarding if the variant is a mutation,  no mutation, FFPE-artfefact deletion or N-call. Whereas the single-hits dict contains positional data for reads with no mate |
| ref_nuc | Nucleotide in reference genome for the variant-record variant position |
| var_nuc | Variant nucleotide for the variant-record variant-call |

| Dict | Structure |
| --- | --- |
| Input | Example dict for a FFPE artefact: inp_dict = umi_key: {"Single_Hits": Str1_Hits: {}, Str2_Hits:{C,T}, "Mate_Hits": Mutation_Hits": {},"FFPE_Hits": {"String_1": C, "String_2": T}, "N_Hits": {}, "Del_Hits": {}, "Reference_Support": 0, "Mutation_Support": 0, "FFPE_Support": 1, "N_Support": 0, "Del_Support": 0}} |
| Output | Example output list for one FFPE artefact: [0;0;1;0;0;1;0;1;0;0;0;0;0] |

### csv_maker
The csv_maker function generates an output CSV-file based on the fumic output containing data for each
variant-record. More specifically regarding the molecular support for the reference genome nucleotide, the variant-call nucleotide,the number of FFPE-calls, the overall frequency of FFPE-artefacts for each variant-record, and the type of mismatch for the variant-record. Through extracting this info and calling the csv_record_maker function, it populates a series of list that are then written to the new CSV.

| Input | Function |
| --- | --- |
| read | Read of interest |
| vcf_file | The output VCF file generated by fumic |
| ref_nuc | Nucleotide in reference genome for the variant-record variant position |
| ffpe_b | Optional input argument controlling which mismatches to consider for FFPE-classification |

### csv_record_maker
The csv_record_maker extracts information from a FFPE-artefact tagged record extracted the output from fumic.py, more specifically from the "samples field. Generates data regarding the molecular support for the nucleotide in the reference genome for the variant-call position, the variant nucleotide in the variant-record, the number of FFPE-artefacts found, the number of unknowns found, the number of deletions found, as well as the percentile ratio of FFPE artefacts in the variant-record. The generated data is then used to populate the lists used as input.

| Input | Function |
| --- | --- |
| read | Read of interest
| pos_lst| List containing the position for records |
| change_lst | List containing the mismatch for studied variant-records |
| ref_lst | List containing the molecular support for the reference genome nucleotide \\ for each studied variant call position respectively |
| var_lst | List containing the molecular support for the variant-record \\ variant-call position for each studied variant call position respectively |
| ffpe_lst | List containing the molcular support for FFPE-artefacts on the variant-record variant-call position for each studied variant call position respectively |
| perc_lst | List containing the percentual ratio of molecules support FFPE-artefacts to the total number of molecules studied for each studied variant call position respectively |
| record | The variant-record of interest |

### umi_maker
The umi_maker function rearranges the UMI-tag belonging to a read, based on if the read is read 1 or read 2 in combination with its directionality. To extract the UMI from the read the ext_fun function is used call either qrn_ext or rx_ext based on user input. The UMI is then transformed into a string and used as input for the spl_fun function. Returns the query-name of the read, the strand it belongs to, and the adjusted UMI-sequence.

| Input | Function |
| --- | --- |
| read | Read of interest
| splt_umi| UMI-tag for the read split into two strings |

### qrn_ext
The qrn_ext function extracts a umi-tag from a read, based on the key being present as the last item in the query-name. Returns the umi-tag.

| Input | Function |
| --- | --- |
| read | Read from which the umi-tag is to be extracted from |

### rx_ext
The rx_ext function extracts a umi-tag from a read, based on the key being present in the RX-tag. Returns the umi-tag.

| Input | Function |
| --- | --- |
| read | Read from which the umi-tag is to be extracted from |

### cha_splt
The cha_splt function splits the umi_string based on the split-character argument. Returns a list containing the umi-tag split into two components.

| Input | Function |
| --- | --- |
| umi_str | A string representing the umi-tag to be split |
| char | Character to split the umi-string by |

### hlf_splt
The hlf_splt function splits the umi_string in half based on its length. Returns a list containing the umi-tag split into two components. 

| Input | Function |
| --- | --- |
| umi_str | A string representing the umi-tag to be split |
| char | Character to split the umi-string by (not used but required by the function call) |

### pos_hits
The pos_hits function selects the most prominent base for a UMI of interest. The function works through iterating through all query-names in the input list and determines if the query-name has a mate or not. The function then calls the base_check function to retrieve the base matching the variant-record position for each read. In the next step the retrieved base is matched against a dict, and depending on the outcome adds to a counter representative of the base. This process is repeated for each query name and its subsequent read, and the resulting dict is then used to determine the most prominent nucleotide for the UMI, effectively collapsing all reads belonging to a UMI. Returns the consensus nucleotide

| Input | Function |
| --- | --- |
| inp_dict | Input list of reads categorized by their query-name |
| rec_pos | The position of the called variant in the reference genome |

| Dict | Structure |
| --- | --- |
| Input | Example dict: input_dict = {example_name_UMI_ACTGCA+ACTGCA: {read1, read2}, example_2_name_UMI_TGACGT+TGACGT: {read2}} |
| Output | Example list: cons_lst = [UMI_tag_1: C, UMI_tag_2: T] |

### base_check
The base_check function checks the variant-record position against the supplemented read, and then extracts the nucleotide belonging to this position in the read. Returns the nucleotide in the read mapping against the variant-record position.

| Input | Function |
| --- | --- |
| read | Read of interest |
| rec_pos | The position of the called variant in the reference genome |

### ffpe_finder
The ffpe_finder function is made to classify the variant type for paired UMI-reads. All-together the UMI and its variant-record position can be classified as: No mutation, Mutation, FFPE-artefact, Unknown (N) or Deletion (-).

The function uses a dict of paired reads containing their consensus nucleotides categorized through their UMI-tag, which is then iterated through for every UMI. It then uses the consensus nucleotide originating from string 1 and string 2 to classify the UMI through comparing these to one another. If the two consensus nucleotides are equal to one another and furthermore equal to the base in the reference genome, the UMI is determined to be "No mutation".If the two consensus nucleotides are equal to one another and furthermore equal to the variant in the variant-record, they are instead deemed to be a "True mutation". In default mode, a FFPE classification only occurs if there is a mismatch between the two consensus nucleotides, if one of the consensus nucleotides is equal to the variant in the variant-record, and finally if the mismatch is of a C:T, T:C or G:A, A:G type. Alternatively, if the flag "ffpe_b" has been called with the input "all", the function instead classifies any mismatch between the  consensus nucleotides for string 1 and string 2 as a FFPE-artefact. If any of the consensus nucleotides are equal to N or -. the UMI is instead deemed to be "Unknown" or "Deletion" respectively.

After a UMI is classified, a counter is added to, and the UMI is stored within a dict named after the variant type. Once the algorithm has iterated through every UMI within the cons\_dict, it creates a new dict containing all variant-type dicts as well as their molecular support, which is then returned

| Input | Function |
| --- | --- |
| cons_dict | Dict containing the consensus nucleotides for String 1 and String 2 classified through their UMI-tag |
| var_nuc | The nucleotide called in the variant-record |
| ref_nuc | The nucleotide found in the reference genome at the variant-call position |
| ffpe_b | Optional input argument controlling which mismatches to consider  for FFPE-classification |

| Dict | Structure |
| --- | --- |
| Input | Example dict for a FFPE artefact: cons_dict = {String_1_Hits: C, String_2_Hits: T} |
| Output | Example dict for a FFPE-artefact: var_dict = {"Single_Hits": Str1_Hits: {}, Str2_Hits:{C,T}, "Mate_Hits": Mutation_Hits": {},"FFPE_Hits": {"String_1": C, "String_2": T}, "N_Hits": {}, "Del_Hits": {}, "Reference_Support": 0, "Mutation_Support": 0, "FFPE_Support": 1, "N_Support": 0, "Del_Support": 0}} |

### mol_count
The mol_count function uses the output generated by the var\_extract function, more specifically support for each variant-type as well as the support for the reference and variant call for str1 and str2. Returns a list consisting of the support for each variant type

| Input | Function |
| --- | --- |
| inp_dict | The output dict from the var_extract function |

| Dict | Structure |
| --- | --- |
| Input | Example dict for a FFPE-artefact: inp_dict = umi_key : {"Single_Hits": Str1_Hits: {}, Str2_Hits:{C,T}, "Mate_Hits": Mutation_Hits": {},"FFPE_Hits": {"String_1": C, "String_2": T}, "N_Hits": {}, "Del_Hits": {}, "Reference_Support": 0, "Mutation_Support": 0, "FFPE_Support": 1, "N_Support": 0, "Del_Support": 0}} |
| Output | Example dict for a FFPE-artefact: var_dict = {"Single_Hits": Str1_Hits: {}, Str2_Hits:{C,T}, "Mate_Hits": Mutation_Hits": {},"FFPE_Hits": {"String_1": C, "String_2": T}, "N_Hits": {}, "Del_Hits": {}, "Reference_Support": 0, "Mutation_Support": 0, "FFPE_Support": 1, "N_Support": 0, "Del_Support": 0}} |
| Output | Example list for a FFPE-artefact: [0,0,1,0,0] |

### nuc_count
The nuc_count function uses the output generated by the var_extract function, more specifically the support for a given nucleotide of interest. Returns a dict containing subsequent dicts with the support for the nucleotide for paired reads on str1, str2, as well as the support for the nucleotide on single reads belonging to string 1 and string 2.

| Input | Function |
| --- | --- |
| inp_dict | The output dict from the var_extract function |
| nuc | The nucleotide of interest |

| Dict | Structure |
| --- | --- |
| Input | Example dict for a FFPE-artefact: inp_dict = umi_key : {"Single_Hits": Str1_Hits: {}, Str2_Hits:{C,T}, "Mate_Hits": Mutation_Hits": {},"FFPE_Hits": {"String_1": C, "String_2": T}, "N_Hits": {}, "Del_Hits": {}, "Reference_Support": 0, "Mutation_Support": 0, "FFPE_Support": 1, "N_Support": 0, "Del_Support": 0}} |
| Output | Example dict for a FFPE-artefact: var_dict = {"Single_Hits": Str1_Hits: {}, Str2_Hits:{C,T}, "Mate_Hits": Mutation_Hits": {},"FFPE_Hits": {"String_1": C, "String_2": T}, "N_Hits": {}, "Del_Hits": {}, "Reference_Support": 0, "Mutation_Support": 0, "FFPE_Support": 1, "N_Support": 0, "Del_Support": 0}} |
| Output | Example dict for the nucleotide C for a imaginary input_dict: n_sup = {"Paired": {"String_1": {C: {12}}, "String_2": C: {11},  "String_1_Single": C: {5}, "String_2_Single": C: {2}} |

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

## FAQ
The FAQ aims to answer questions the reader may have regarding FUMIC and its use.

### Why do i need a VCF-file to use FUMIC? 
The VCF is required to identify and target known variant calls. As a FFPE-artefact will be initially identified as a variant in the single-reads, they are then used to evaluate if they are in fact FFPE-artefacts or not. We as developers felt no need to construct a custom variant-caller due to the abundance of excellent and sophisticated variant-calling software already existing. 

### Why is there a need for a BAI file, is the BAM not enough?
A BAM file is stored in binary format, and thus has no internal structure to use for purposes of fetching. The BAI file provides an indexed form of the BAM, allowing us to fetch specific reads, which is a function FUMIC requires. 

### The molecular support for the variant-call variants in my output seems comparatively low to the read-depth, why is this? 
FUMIC can only classify variant types if the UMI it is studying has both a string 1 and string 2 consensus sequence. Thus, if the studied data has many singletons or if only one of the strings align to the variant-call position, there will be a limited amount of classifications that can be made. The single paired reads and the singletons instead provide molecular support for the variant call nucleotide as well as the reference genome nucleotide for the variant-call position. 

### What are the recommended number of threads to run FUMIC on?:
The optimal number of threads depends entirely on the system FUMIC is to be run on. Ideally, one should use n-1 threads, where n is the total number of available processors on the system. 

### How do i know where the UMI-tag is stored, and how do i know what character is separating it (if there is any)? 
Due to the high variability of sequence data, this is not an automated process by FUMIC, but instead requires the user to manually inspect one of their reads. We recommend the user to utilize the free software SAMtools \cite{SAMtools} for this purpose.

### Where can i dowload FUMIC, and is FUMIC free?
FUMIC is a publically available software, and can be acquired free of charge through its github page: https://github.com/clinical-genomics-uppsala/fumic

### I have identified a bug or have suggestions for improving the program, where can i contact you regarding this?
All inquires regarding updating the software should be done through FUMIC's github: https://github.com/clinical-genomics-uppsala/fumic

### I have a question not covered by the FAQ, where can i contact you regarding this?:
All inquires regarding questions about the program should be done through FUMIC's github: https://github.com/clinical-genomics-uppsala/fumic

## Special thanks and credits
FUMIC were developed by Hugo Swenson and Clinical Genomics, Uppsala University in 2019. Special thanks and credit goes to Claes Ladenvall & Patrik Smeds at Clinical Genomics Uppsala for their support and expertise in regards to developing this program. Furthermore i would like Adam Ameur at Department of Immunology, Genetics and Pathology and Uppsala Genome Center for creative feedback and suggestions on how to improve the program.

