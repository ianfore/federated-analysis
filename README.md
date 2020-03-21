# federated-analysis

HIPAA requires that medical records be kept private, yet scientists need this data in order to conduct their research.  So instead of sharing data, hospitals and clinics can run software on that data and share the summarized, aggregated, anonymized results with the interested scientists.  In this way, people's privacy is protected, and science can move forward on important research.

This software repository contains scripts that achieve the above goal.  

To run a co-occurrence analysis, perform the following steps:
1. Change directory to the top-level directory of the repository.

cd federated-analysis/

2. Run the runMe.sh script as follows:

./runMe.sh BreastCancer.shuffle.vcf 37 75 [13,17] False

where:
* BreastCancer.shuffle.vcf is the name of the VCF file in the federated/analysis/data directory

* 37 is the version of the human genome for the coordinates in the VCF file

* 75 is the build of ensembl for retrieving gene names for the coordinates

* [13,17] is a list of chromosomes which contain the SNPs of interest

* False is a boolean value for whether the VCF data is phased (True) or not (False)

3. This will generate a report in federated-analysis/data called vusFinalData.json which contains a list of VUS, each in the following format:

"(13, 32914461, 'A', 'C')": [ 	<-- VUS
        0.0005026755310523755,  <-- P(VUS is benign with a pathogenic variant in trans) 
        0.0001,			<-- P(VUS is pathogenic with another pathogenic variant in trans)
        82,			<-- n = number of times this VUS occurred in cohort
        2,			<-- k = number of times this VUS co-occurred with pathogenic variants
        0.04087136171331834,	<-- likelihood this variant is pathogenic = (p2^k)*(1-p2)^(n-k) / (p1^k)*(1-p1)^(n-k)
        [
            [			<-- co-occurring pathogenic variant #1
                13,
                32914936,
                "TATTAA",
                "T"
            ],
            [			<-- co-occurring pathogenic variant #2
                13,
                32920978,
                "C",
                "T"
            ]
        ]
    ],




Additionally, this software allows users to run a Docker container which has the necessary code to perform basic statistical analysis and validity checking.  There's a configuration file that the cooperating owner of the data must fill out in conjunction with the scientist to define the fields of interest in the data set.  
There are three Python modules in this repository: one which performs default analysis (dataAnalyzer.py), one that performs custom analysis (customDataAnalyzer.py), and one that creates a table (supplementaryTable4.py).  

The default data analysis outputs the following information:
1. data file name 
2. data file header present?
3. data file field delimiter
4. total records read from data file
5. for each field of interest
    - name of field
    - type of field (numerical, categorical, or free-form)
    - total counts for each value
    - for numerical data, the min, max, mean, and median values
6. bad values (those that don't conform to the type)
7. missing values (fields that are not populated)

Note that the default data analysis is generic -- it's completely devoid of any application or context.  If the scientist wishes to perform specific analyses on the data, then they must implement the custom data analyzer.  The custom data analyzer is provided an object that encapsulates all the default data analysis.  The custom code can then perform application-specific analyses on the data. 


In order to use this solution, perform the following steps.

1. Change directory to the top-level directory of the repository.

cd federated-analysis/

2. Optionally, edit the config/conf.json file to reflect the metadata regarding the data file (file name, header line, field delimiter) as well as the correct fields of interest.

vi config/conf.json

3. Run the runMe.sh script as follows:

./runMe.sh


To run the unit tests, perform the following steps:

1. change directory to the top-level directory of the repository

cd federated-analysis/

2. run the following command:

python -m unittest tests.test_dataAnalyzer
