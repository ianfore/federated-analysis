# Summary 

HIPAA requires that medical records be kept private, yet scientists need this data in order to conduct their research.  So instead of sharing data, hospitals and clinics can run software on that data and share the summarized, aggregated, anonymized results with the interested scientists.  In this way, people's privacy is protected, and science can move forward on important research.

This software repository contains scripts that achieve the above goal.  

# Cooccurrence analysis
Run the following steps to find variants of uncertain significance (VUS) in a VCF file that co-occur with known pathogenic variants in the BRCA Exchange database.

## Prepare for co-occurrence 
To prepare for a co-occurrence analysis, perform the following steps:

1. Clone this github repository to your local system where the VCF file resides.

```console
$ **git clone https://github.com/BRCAChallenge/federated-analysis **
```

2. Put a copy of the VCF file in the federated-analysis/data directory.

```console
$ cp <your-vcf-file> federated-analysis/data
```

3. Make sure the VCF file has read permissions for the world.

```console
$ chmod a+r federated-analysis/data/*.vcf
```
 

## Run co-occurrence analysis container from command line
To run the co-occurrence analysis, perform the following steps:

1. Change directory to the top-level directory of the repository.

```console
$ cd federated-analysis/
```

2. Run the runMe.sh script as follows:

```console
$ ./runMe_nontopmed.sh BreastCancer.shuffle.vcf 37 75 13 False BRCA2 

$ ./runMe_nontopmed.sh BreastCancer.shuffle.vcf 37 75 17 False BRCA1
```

where:
* BreastCancer.shuffle.vcf is the name of the VCF file in the federated-analysis/data directory

* 38 is the version of the human genome for the coordinates in the VCF file (37 is also supported)

* 99 is the build of ensembl for retrieving gene names for the coordinates (use 75 for hg37)

* 13 (or 17) is the chromosome to filter in the VCF file

* True is a boolean value for whether the VCF data is phased (True) or not (False)

* BRCA2 (or BRCA1) is the name of the gene of interest on the chromosome of interest

3. This will generate a report in federated-analysis/data called output.json which contains a list of VUS, each in the following format:

```json
"(13, 32317399, 'T', 'G')": {
            "likelihood data": {
                "p1": 0.25,
                "p2": 0.001,
                "n": 2,
                "k": 2,
                "likelihood": 1.6e-05
            },
            "allele frequencies": {
                "maxPop": "Allele_frequency_genome_AMR_GnomAD",
                "maxPopFreq": 0.001179,
                "cohortFreq": 1.0
            },
            "pathogenic variants": [
                [
                    13,
                    32338277,
                    "G",
                    "T"
                ]
            ],
            "RARE": true

```


# Pathology statistics

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

```console
cd federated-analysis/
```

2. Edit the config/conf.json file to reflect the metadata regarding the data file (file name, header line, field delimiter) as well as the correct fields of interest.

```console
vi config/conf.json
```

3. Run the runMe.sh script as follows:

```console
./runMe.sh analyze
```

# Software unit testing 

To run the unit tests, perform the following steps:

1. change directory to the top-level directory of the repository

```console
cd federated-analysis/
```

2. run the following command:

```console
python -m unittest tests.test_dataAnalyzer
```

# Run co-occurrence analysis container from WDL
Define the environment variables for the workflow, such as the following example shows:
```
PYTHON_SCRIPT=/home/jcasalet/src/cooccurrenceFinder.py
VCF_FILE=/data/chr13_brca2.vcf 
BRCA_FILE=/data/brca-variants.tsv 
OUTPUT_FILENAME=13-out.json 
HG_VERSION=38 
ENSEMBL_RELEASE=99 
PHASED=True 
SAVE_VARS=True 
CHROM=13 
GENE=BRCA2 
NUM_CORES=16
INDIVIDUALS_FILENAME=13-vpi.json
VARIANTS_FILENAME=13-ipv.json
```

Then run the WDL workflow using the values of the environment variables defined above.

