# PIRCHE Anonymization Client

The purpose of this script is to be able to send anonymized HLA typing data of patient and donor pairs to the PIRCHE web service to calculate the PIRCHE scores. 
This scripts runs locally and only the anonymized data is transferred to the PIRCHE web service.

__Important note: The current version is for testing and aligning on the input and output formats only. Anonymization is supported but not final therefore disabled by default. It can be enabled by the -a parameter. No low resolution data support yet.__ 

## Requirements
1. Working Python 3 environment with the modules needed by the script (see imports)
2. HLA typing data stored in a csv file according to the structure [importTemplate.csv](importTemplate.csv)
3. An account on the PIRCHE web service (https://www.pirche.com) with API access enabled and API access token configured

## Running the script
Just run the script locally and provide all parameters needed. The results will be stored in a file containing PIRCHE I and II scores separated for each loci as well as added up to one PIRCHE I and II score. 

Script Parameters:

| Short | Long              |Required| Description                                                                               |
|:------|:-------------     |:------:|:-----                                                                                     |
| -v    | --verbose         |        |Verbose mode                                                                               |
| -url  | --url             | x      |URL to the PIRCHE web service                                                              |
| -u    | --user            | o      |PIRCHE web service user*                                                                   |
| -p    | --password        | o      |PIRCHE web service user password*                                                          |
| -k    | --apikey          | o      |PIRCHE web service user API Key*                                                           |
| -i    | --input           | x      |HLA typing data input file                                                                 |
| -pp   | --population      |        |Population for HLA typing data provided (needed for low res high res conversion)           |
| -hp   | --haplotypes      | x      |NMDP haplotype table (either 2007 or 2011 or equally formatted)                            |
| -t    | --threshold       | x      |Frequency threshold for haplotypes generation 0.0 to 1.0                                   |
| -ps   | --population_short| x      |Population short code as used in the NMDP haplotype table header row                       |
| -o    | --output          | x      |Output file name                                                                           |
| -a    | --anonymization   |        |Enable anonymization. Default - False - no anonymization. To enable it set parameter to True                             |
| -s    | --salt            | o      |Salt (password) used to anonymize input data. Use identical password when submitting same HLA data set multiple times.** |
| -k    | --kanonymization  |        |Number of smoke hla data sets (genotypes) generated per patient and per donor. Default - 5.                              |

*only either of user & password or api key credentials are needed
**required if anonymization is enabled

Example call (with anonymization enabled):<br>
`-v -url https://www.pirche.com -key apikey -i importTemplate.csv -hp NMDP2011_5_locus.xls -t 0.9 -ps EURCAU -o PIRCHE_results.csv -a True -s xyz123pass`
