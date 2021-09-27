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
{
    "cooccurring vus": {
        "(13, 32353519, 'A', 'G')": {
            "likelihood data": {
                "p1": 0.375,
                "p2": 0.001,
                "n": 1,
                "k": 1,
                "likelihood": 0.0026666666666666666
            },
            "allele frequencies": {
                "maxPop": null,
                "maxPopFreq": 0.0,
                "minPop": null,
                "minPopFreq": 1.0,
                "cohortFreq": 0.25
            },
            "pathogenic variants": [
                [
                    13,
                    32338749,
                    "AATTAC",
                    "A"
                ]
            ]
        },
        "(13, 32355250, 'T', 'C')": {
            "likelihood data": {
                "p1": 0.375,
                "p2": 0.001,
                "n": 2,
                "k": 1,
                "likelihood": 0.0042624
            },
            "allele frequencies": {
                "maxPop": null,
                "maxPopFreq": 0.0,
                "minPop": null,
                "minPopFreq": 1.0,
                "cohortFreq": 0.5
            },
            "pathogenic variants": [
                [
                    13,
                    32316508,
                    "GAC",
                    "G"
                ]
            ]
        },
        "(13, 32353470, 'A', 'C')": {
            "likelihood data": {
                "p1": 0.375,
                "p2": 0.001,
                "n": 1,
                "k": 1,
                "likelihood": 0.0026666666666666666
            },
            "allele frequencies": {
                "maxPop": null,
                "maxPopFreq": 0.0,
                "minPop": null,
                "minPopFreq": 1.0,
                "cohortFreq": 0.25
            },
            "pathogenic variants": [
                [
                    13,
                    32340836,
                    "GACAA",
                    "G"
                ]
            ]
        }
    },
    "homozygous vus": {
        "(13, 32355250, 'T', 'C')": {
            "count": 1,
            "maxPop": null,
            "maxPopFreq": 0.0,
            "minPop": null,
            "minPopFreq": 1.0,
            "cohortFreq": 0.25
        }
    }
}

```

3. This will also create a JSON file called `BRCA2-intersection.json` in the federated-analysis/examples/BRCA/data directory which intersect the phenotype data with the genotype data. 

````json
{
    "cooccurring": {
        "(13, 32353470, 'A', 'C')": {
            "phenotype": [
                {
                    "Age at onset": [
                        51
                    ],
                    "CarrierGene": [
                        "BRCA2"
                    ],
                    "ER": [
                        "Negative"
                    ],
                    "Family history / breast cancer": [
                        0
                    ],
                    "HER2": [
                        "1+"
                    ],
                    "ID": [
                        2
                    ],
                    "PgR": [
                        "Positive"
                    ]
                }
            ]
        },
        "(13, 32353519, 'A', 'G')": {
            "phenotype": [
                {
                    "Age at onset": [
                        66
                    ],
                    "CarrierGene": [
                        "BRCA2"
                    ],
                    "ER": [
                        "Negative"
                    ],
                    "Family history / breast cancer": [
                        0
                    ],
                    "HER2": [
                        "Negative"
                    ],
                    "ID": [
                        3
                    ],
                    "PgR": [
                        "Negative"
                    ]
                }
            ]
        },
        "(13, 32355250, 'T', 'C')": {
            "phenotype": [
                {
                    "Age at onset": [
                        57
                    ],
                    "CarrierGene": [
                        "BRCA2"
                    ],
                    "ER": [
                        "Positive"
                    ],
                    "Family history / breast cancer": [
                        1
                    ],
                    "HER2": [
                        "3+"
                    ],
                    "ID": [
                        1
                    ],
                    "PgR": [
                        "Positive"
                    ]
                }
            ]
        }
    },
    "homozygous": {
        "(13, 32355250, 'T', 'C')": {
            "phenotype": [
                {
                    "Age at onset": [
                        0
                    ],
                    "CarrierGene": [
                        "NonCarrier"
                    ],
                    "ER": [
                        NaN
                    ],
                    "Family history / breast cancer": [
                        0
                    ],
                    "HER2": [
                        NaN
                    ],
                    "ID": [
                        4
                    ],
                    "PgR": [
                        NaN
                    ]
                }
            ]
        }
    },
    "numMissing": 0
}

````

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
