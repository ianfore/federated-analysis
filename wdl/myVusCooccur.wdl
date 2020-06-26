version 1.0

workflow cooccurrence {
	input {
	File PYTHON_SCRIPT
	File VCF_FILE
	File BRCA_FILE
	String OUTPUT_FILENAME
	String VPI_FILENAME
	String IPV_FILENAME
	String HG_VERSION
	String ENSEMBL_RELEASE
	String PHASED
	String CHROM
	String GENE
	String NUM_CORES
	}

	call run_cooccurrence {
	input: 
		python_script=PYTHON_SCRIPT,
		vcf_file=VCF_FILE,
		brca_file=BRCA_FILE,
		hg_version=HG_VERSION,
		ensembl_release=ENSEMBL_RELEASE,
		phased=PHASED,
		chrom=CHROM,
		gene=GENE,
		num_cores=NUM_CORES,
		output_filename=OUTPUT_FILENAME,
		individuals_filename=VPI_FILENAME,
		variants_filename=IPV_FILENAME,
	}
		
}

task run_cooccurrence {
	input {
		File python_script
		File vcf_file
		File brca_file
		String output_filename
		String individuals_filename
		String variants_filename
		String hg_version
		String ensembl_release
		String phased
		String chrom
		String gene
		String num_cores
	}

	command <<<
		export PYTHONPATH=/ 
		export PYTHONIOENCODING=UTF-8 
		/usr/bin/python3  ~{python_script} --v ~{vcf_file} --o ~{output_filename} --h ~{hg_version} --e ~{ensembl_release} --c ~{chrom} --g ~{gene} --p ~{phased} --n ~{num_cores} --b ~{brca_file} --vpi ~{individuals_filename} --ipv ~{variants_filename}
	>>>

	runtime {
		docker: 'brcachallenge/federated-analysis:cooccurrence'
    		memory: "8192 MB"
    		disk: "local-disk 20 HDD"   ## hardcoded disk size (20) and type (HDD)
	} 
}


task setup_perms {
	input {
		File output_file
	}
	command <<<
		if [ $(uname) == "Darwin" ]
		then
			PREV_PERMS=$(stat -f "%OLp" ${DATA_PATH})
		else
			PREV_PERMS=$(stat -c "%a" ${DATA_PATH})
		fi
		chmod 1777 ${DATA_PATH}
	>>>

}
