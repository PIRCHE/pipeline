# PIRCHE Anonymization Client

The purpose of this script is to be able to send anonymized HLA typing data of patient and donor pairs to the PIRCHE web service to calculate the PIRCHE scores. 
This scripts runs locally and only the anonymized data is transferred to the PIRCHE web service.

__Important note: The current version is feature complete but still in testing. Anonymization is supported and enabled by default.__ 

## Requirements
1. Working Python 3 environment with the modules needed by the script (see imports)
2. HLA typing data stored in a csv file according to the structure [importTemplate.csv](importTemplate.csv). So one line per patient donor pair. All values are separated by a comma and structured according to the header. Currently only the locus A, B, C, DRB1, DQB1 are supported due to limited haplotype data available for generating the "smoke/fake" data needed for anonymization. 
   
   The following HLA resolutions are supported:
   1. molecular high resolution: C*05:01 (more than 4 digits can be provided BUT this is not considered when calculating the PIRCHE score)
   2. molecular low resolution: C*05 (2 digits, with asterisk) 
   3. serological equivalent (low resolution): C5 (so if it is a serological value you should provided it as such when using the anonymization client - without asterisk)

3. An account on the PIRCHE web service (https://www.pirche.com) with API access enabled and API access token configured
4. The following files must be downloaded and present in the script execution folder/location:<br>
4.1 IMGTHLA g-groups reference table (rel_dna_ser.txt) - [Link](https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_g.txt) (Right click on Link and select "Save Link As ..")<br>
4.2 IMGTHLA dna-ser reference table (hla_nom_g.txt) - [Link](https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/rel_dna_ser.txt) (Right click on Link and select "Save Link As ..")<br>
4.3 IMGTHLA hla alleles reference table (Allelelist.txt) - [Link](https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt) (Right click on Link and select "Save Link As ..")<br>
4.4 NMDP haplotype table file - [Link](https://bioinformatics.bethematchclinical.org/WorkArea/DownloadAsset.aspx?id=6383) --> file must be renamed to <code>2007_haplotypes.xls</code> after download<br> 

## Running the script
Just run the script locally and provide all parameters needed. The results will be stored in a file containing PIRCHE I and II scores separated for each loci as well as added up to one PIRCHE I and II score. 

Parameters:

| Short | Long              |Needed<sup>1</sup>| Description                                                                                                                        |
|:------|:-------------     |:------:  |:-----                                                                                                                                      |
| -v    | --verbose         |          |Verbose mode                                                                                                                                |
| -url  | --url             | x        |URL to the PIRCHE web service                                                                                                               |
| -u    | --user            | o        |PIRCHE web service user <sup>2</sup>                                                                                                        |
| -p    | --password        | o        |PIRCHE web service user password <sup>2</sup>                                                                                               |
| -k    | --apikey          | o        |PIRCHE web service user API Key <sup>2</sup>                                                                                                |
| -i    | --input           | x        |HLA typing data input file                                                                                                                  |
| -o    | --output          | x        |Output file name                                                                                                                            |
| -a    | --anonymization   | d        |Enable anonymization. Default - True - Anonymization enabled                                                                                |
| -s    | --salt            | x        |Salt (password) used to anonymize input data. Use identical password when submitting same HLA data set multiple times. <sup>3</sup>         |
| -k    | --kanonymization  | d        |Number of smoke hla data sets (genotypes) generated per patient and per donor. Default - 5.                                                 |
| -pp   | --population      | d        |Population for HLA typing data provided (needed for low res high res conversion). Default - NMDP EUR haplotypes (2007).                     |
| -gg   | --ggroups         | d        |HLA g-groups reference table file name and path. Default - hla_nom_g.txt                                                                    |
| -ds   | --dstable         | d        |HLA dna ser reference table file and path. Default - rel_dna_ser.txt                                                                        |
| -al   | --allelelist      | d        |HLA alleles reference table file. Default - Allelelist.txt                                                                                  |
| -ht   | --haplotypes      | d        |NMDP haplotype table file name and path (either 2007 or 2011 or equally formatted). Default - 2007_haplotypes.xls <sup>4</sup>              |
| -hf   | --haplofileformat | d        |NMDP haplotype table file alleles format (either alleles XXXX (2007) or locus + alleles L*XX:XX (2011)). Default - XXXX <sup>4</sup>        |
| -hp   | --haplofilepop    | d        |NMDP haplotype table file population short code as used in the header row (e.g. EUR (2007) or EURCAU (2011)) . Default - EUR <sup>4</sup>   |
| -hth  | --haplothreshold  | d        |Frequency threshold for haplotypes generation 0.0 to 1.0. Default - 0.8 <sup>4</sup>                                                        |
| -prx  | --proxyhttp		| o        |HTTP proxy - full http proxy url (ip or dns) with protocol and port (http(s)://proxy:port)                                                  |
| -prxs | --proxyhttps		| o        |HTTPS proxy - full https proxy url (ip or dns) with protocol and port (http(s)://proxy:port)                                                |

<sup>1</sup> x - value must be given | o - value is optional | d - value has a default but can be overwritten<br>
<sup>2</sup> only either of user & password or api key credentials are needed<br>
<sup>3</sup> not required if anonymization is disabled<br>
<sup>4</sup> if haplotype table file is changed format and population must be changed accordingly; threshold might be adjusted as well<br>

Example call (with anonymization enabled and using default haplotype table values):<br>
`-v -url https://www.pirche.com -key apikey -i importTemplate.csv -o PIRCHE_results.csv -s xyz123pass`

Example call (with anonymization enabled and using custom haplotype table values):<br>
`-v -url https://www.pirche.com -key apikey -i importTemplate.csv -o PIRCHE_results.csv -s xyz123pass -ht 2011_haplotypes.xls -hf L*XX:XX -hp EURCAU -hth 0.7`
