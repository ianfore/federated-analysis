version 1.0

workflow cooccurrence {
	input {
	File PYTHON_SCRIPT
	File VCF_FILE
	File VARIANT_PATHOGENICITY_FILE
	File? ANNO_FILE
	File  GNOMAD_FILE
	String OUTPUT_FILENAME
	String VPI_FILENAME
	String IPV_FILENAME
	String ALL_FILENAME
	String TOUT_FILENAME
	String HG_VERSION
	String ENSEMBL_RELEASE
	String PHASED
	String P2
	String CHROM
	String GENE
	String NUM_CORES
	String SAVE_FILES
	}

	call run_cooccurrence {
	input: 
		python_script=PYTHON_SCRIPT,
		vcf_file=VCF_FILE,
		variant_pathogenicity_file=VARIANT_PATHOGENICITY_FILE,
		anno_file=ANNO_FILE,
		gnomad_file=GNOMAD_FILE,
		hg_version=HG_VERSION,
		ensembl_release=ENSEMBL_RELEASE,
		phased=PHASED,
		p2=P2,
		chrom=CHROM,
		gene=GENE,
		num_cores=NUM_CORES,
		output_filename=OUTPUT_FILENAME,
		individuals_filename=VPI_FILENAME,
		all_filename=ALL_FILENAME,
		variants_filename=IPV_FILENAME,
		tout_filename=TOUT_FILENAME,
		save_files=SAVE_FILES
	}
	output { 
		File ipv_file = IPV_FILENAME
                File vpi_file = VPI_FILENAME 
                File out_file = OUTPUT_FILENAME 
		File all_file = ALL_FILENAME
		File tout_file = TOUT_FILENAME
	}

		
}

task run_cooccurrence {
	input {
		File python_script
		File vcf_file
		File variant_pathogenicity_file
		File? anno_file
		File gnomad_file
		String hg_version
		String ensembl_release
		String phased
		String p2
		String chrom
		String gene
		String num_cores
		String output_filename
		String individuals_filename
		String variants_filename
		String tout_filename
		String all_filename
		String save_files
	}

	command <<<
		export PYTHONPATH=/ 
		export PYTHONIOENCODING=UTF-8 
		/usr/bin/python3  ~{python_script} --vcf ~{vcf_file} --h ~{hg_version} --e ~{ensembl_release} --c ~{chrom} --g ~{gene} --p ~{phased} --p2 ~{p2} --n ~{num_cores} --vpf ~{variant_pathogenicity_file} --anno "~{anno_file}" --save ~{save_files} --gf ~{gnomad_file}
	>>>
	
	output {
		File actual_ipv_file = "~{variants_filename}"
		File actual_vpi_file = "~{individuals_filename}"
		File actual_out_file = "~{output_filename}"
		File actual_tout_file = "~{tout_filename}"
		File actual_all_file = "~{all_filename}"
		

	}

	runtime {
		docker: "brcachallenge/federated-analysis:cooccurrence"
		memory: "8192 MB"
		bootDiskSizeGb: 50
    		disk: "local-disk 100 HDD"   ## hardcoded disk size (20) and type (HDD)
    		#disk: "local-disk 200"   ## hardcoded disk size (20) and type (HDD)
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
