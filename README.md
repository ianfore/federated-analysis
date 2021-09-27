# Summary 

HIPAA requires that medical records be kept private, yet scientists need this data in order to conduct their research.  So instead of sharing data, hospitals and clinics can run software on that data and share the summarized, aggregated, anonymized results with the interested scientists.  In this way, people's privacy is protected, and science can move forward on important research.

This software repository contains scripts that achieve the above goal.  

# Cooccurrence analysis
Run the following steps to find variants of uncertain significance (VUS) in a VCF file that co-occur with known pathogenic variants in the BRCA Exchange database.

## Prepare for the co-occurrence analysis
To prepare for a co-occurrence analysis, perform the following steps:

1. Clone this github repository to your local system where the VCF file resides.

```console
$ git clone https://github.com/BRCAChallenge/federated-analysis
```

2. Change directory to the federated-analysis directory. 

```console
$ cd federated-analysis 
```

3. Examine the VCF file in the examples/data directory.

```console
$ cat examples/BRCA/data/brca2.vcf
```

4. Examine the phenotype data file in the examples/data directory.

```console
$ cat examples/BRCA/data/brca2-pathology.tsv
```

5. Examine the variant pathogenicity file in the examples/data directory.

```console
$ cat examples/BRCA/data/brca2-pathogenicity.tsv
```

6. Examine the data quality report configuration file in the examples/config directory.

```console
$ cat examples/BRCA/config/brca2-report-config.tsv
```

## Run the report workflow.

1. Run the runMe.sh script.

```console
$ ./runMe.sh -c 13 -p True -g BRCA2 -dd $(pwd)/examples/BRCA/data -cd $(pwd)/examples/config -vf brca2.vcf -vpf brca2-pathogenicity.tsv -rc brca2-report-config.json -spf brca2-pathology.tsv
```

where:
* -c is the chromosome number

* -p is the Boolean phased value (True or False)

* -g is the gene name

* -dd is the absolute path to the data directory

* -cd is the absolute path to the configuration directory

* -vf is the VCF file in the data directory

* -sp is the sample pathology file in the data directory 

* -rc is the data quality report configuration file in the config directory

* -spf is the sample phenotype data file in the data directory 


2. This will generate a report in the federated-analysis/examples/BRCA/data directory called `BRCA2-cooccurrence.json` which contain a list of VUS, each in the following format:

```json
"(13, 32911164, 'T', 'A')": {
            "likelihood data": {
                "p1": 0.0015891032917139615,
                "p2": 0.001,
                "n": 25,
                "k": 1,
                "likelihood": 0.6382577479687377
            },
            "allele frequencies": {
                "maxPop": null,
                "maxPopFreq": 0.0,
                "minPop": null,
                "minPopFreq": 0.0,
                "cohortFreq": 0.0008107669855683476
            },
            "pathogenic variants": [
                [
                    13,
                    32911297,
                    "TAAAC",
                    "T"
                ]
            ]

```

3. This will also create a JSON file called `BRCA2-intersection.json` in the federated-analysis/examples/BRCA/data directory which intersect the phenotype data with the genotype data. 

4. Finally, this will also create a data quality report called `brca2-data-quality-report.txt` in the examples/BRCA/data directory.  Note that the default data analysis is generic -- it's completely devoid of any application or context.  If the scientist wishes to perform specific analyses on the data, then they must implement the custom data analyzer.  The custom data analyzer is provided an object that encapsulates all the default data analysis.  The custom code can then perform application-specific analyses on the data. 


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
