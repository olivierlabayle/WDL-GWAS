version 1.0

import "structs.wdl"


task regenie_step_1 {
    input {
        String docker_image
        String group_name
        File bed_file
        File bim_file
        File fam_file
        File sample_list
        File covariates_file
        Array[String] covariates_list
        String cv_folds = "5"
        String bsize = "1000"
        String maf = "0.01"
        String mac = "10"
    }

    # String group_name = sub(basename(sample_list, ".txt"), "gwas.individuals.", "")

    command <<<
        genotypes_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)

        # In principle this is already taken care of in `make_group_bed_qced`, but we ensure here that the input is correct because it is cheap.
        # Only retain bi-allelic (REGENIE can't handle non-biallelic) and frequent variants within the sample list
        plink2 \
            --bfile ${genotypes_prefix} \
            --keep ~{sample_list} \
            --mac ~{mac} \
            --maf ~{maf} \
            --min-alleles 2 \
            --max-alleles 2 \
            --write-snplist \
            --out biallelic_frequent

        # Parse cross-validation option
        cv_option="--cv ~{cv_folds}"
        if [[ ~{cv_folds} == "loocv" ]]; then
            cv_option="--loocv"
        fi

        # phenotype from group_name
        phenotype=$(echo ~{group_name} | cut -d'.' -f2)

        # Find the type of the phenotype (quantitative or binary)
        phenotype_col_idx=$(head -1 ~{covariates_file} | tr '\t' '\n' | grep -n ${phenotype} | cut -d: -f1)
        uniq_vals_count=$(cut -f"${phenotype_col_idx}" ~{covariates_file} | sort -u | grep -v "NA" | wc -l)
        trait_type="--bt"
        if [ "${uniq_vals_count}" -gt 3 ]; then # two binary values + header
            trait_type="--qt"
        fi

        conda run -n regenie_env regenie \
            --step 1 \
            --bed ${genotypes_prefix} \
            --keep ~{sample_list} \
            --extract biallelic_frequent.snplist \
            --phenoFile ~{covariates_file} \
            --phenoColList ${phenotype} \
            --covarFile ~{covariates_file} \
            --covarColList ~{sep="," covariates_list} \
            --minMAC ~{mac} \
            ${cv_option} \
            ${trait_type} \
            --bsize ~{bsize} \
            --lowmem \
            --out ~{group_name}.step1
        awk '{sub(".*/", "", $2); print $1, $2}' ~{group_name}.step1_pred.list > ~{group_name}.step1_pred.listrelative
    >>>

    output {
        Array[File] phenotypes_loco = glob("${group_name}.step1_*.loco")
        File list = "${group_name}.step1_pred.listrelative"
        # RegenieStep1Files step1_files = object {
        #     phenotypes_loco: glob("${group_name}.step1_*.loco"),
        #     list: "${group_name}.step1_pred.listrelative"
        # }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }
}

workflow gwas_step_1 {
    input {
        String docker_image
        String gwas_software
        String group_name
        File bed_file
        File bim_file
        File fam_file
        File sample_list
        File covariates_file
        Array[String] covariates_list
        String cv_folds
        String bsize
        String maf
        String mac
    }

    if (gwas_software == "regenie") {
        call regenie_step_1 {
            input:
                docker_image = docker_image,
                group_name = group_name,
                bed_file = bed_file,
                bim_file = bim_file,
                fam_file = fam_file,
                sample_list = sample_list,
                covariates_file = covariates_file,
                covariates_list = covariates_list,
                cv_folds = cv_folds,
                bsize = bsize,
                maf = maf,
                mac = mac
        }
    }

    output {
        Array[File]? regenie_loco_preds = regenie_step_1.phenotypes_loco
        File? regenie_list = regenie_step_1.list
        # RegenieStep1Files? regenie_step_1_files = regenie_step_1.step1_files
    }
}