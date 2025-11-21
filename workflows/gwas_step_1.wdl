version 1.0

import "structs.wdl"

task make_fake_regenie_step_1 {
    input {
        String group_name
    }

    command <<<
        touch "~{group_name}.step1_1.loco"
        touch "~{group_name}.step1_pred.listrelative"
    >>>

    output {
        Array[File] phenotypes_loco = glob("${group_name}.step1_*.loco")
        File list = "${group_name}.step1_pred.listrelative"
    }
}

task make_fake_saige_step_1 {
    input {
        String group_name
    }

    command <<<
        touch "~{group_name}.rda"
        touch "~{group_name}.varianceRatio.txt"
    >>>

    output {
        File model_file = "${group_name}.rda"
        File variance_ratio_file = "${group_name}.varianceRatio.txt"
    }
}

task saige_step_1 {
    input {
        String docker_image
        String group_name
        File bed_file
        File bim_file
        File fam_file
        File sample_list
        File covariates_file
        Array[String] covariates_list
        String maf = "0.01"
        String mac = "10"
        String npcs = "10"
    }

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
            --make-bed \
            --out biallelic_frequent

        # phenotype from group_name
        phenotype=$(echo ~{group_name} | cut -d'.' -f2)

        # Make covariates list
        pc_list=$(printf "CHR0_OUT_PC%s," {1..~{npcs}} | sed 's/,$//')
        full_covariates_list="~{sep="," covariates_list},${pc_list}"

        # Find the type of the phenotype (quantitative or binary)
        phenotype_col_idx=$(head -1 ~{covariates_file} | tr '\t' '\n' | grep -n ${phenotype} | cut -d: -f1)
        uniq_vals_count=$(cut -f"${phenotype_col_idx}" ~{covariates_file} | sort -u | grep -v "NA" | wc -l)
        trait_type="binary"
        if [ "${uniq_vals_count}" -gt 3 ]; then # two binary values + header
            trait_type="quantitative --invNormalize=TRUE"
        fi

        step1_fitNULLGLMM.R     \
            --plinkFile=biallelic_frequent \
            --phenoFile=~{covariates_file} \
            --phenoCol=${phenotype} \
            --covarColList=${full_covariates_list} \
            --sampleIDColinphenoFile=IID \
            --traitType=${trait_type} \
            --outputPrefix=~{group_name} \
            --nThreads=$(nproc)	\
            --IsOverwriteVarianceRatioFile=TRUE
    >>>

    output {
        File model_file = "${group_name}.rda"
        File variance_ratio_file = "${group_name}.varianceRatio.txt"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }
}

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

        call make_fake_saige_step_1 {
            input:
                group_name = group_name
        }
    }

    if (gwas_software == "saige") {
        call saige_step_1 {
            input:
                docker_image = docker_image,
                group_name = group_name,
                bed_file = bed_file,
                bim_file = bim_file,
                fam_file = fam_file,
                sample_list = sample_list,
                covariates_file = covariates_file,
                covariates_list = covariates_list,
                maf = maf,
                mac = mac
        }

        call make_fake_regenie_step_1 {
            input:
                group_name = group_name
        }
    }

    output {
        Array[File]? regenie_loco_preds = if (gwas_software == "regenie") then regenie_step_1.phenotypes_loco else make_fake_regenie_step_1.phenotypes_loco
        File? regenie_list = if (gwas_software == "regenie") then regenie_step_1.list else make_fake_regenie_step_1.list
        File? saige_model_file = if (gwas_software == "saige") then saige_step_1.model_file else make_fake_saige_step_1.model_file
        File? saige_variance_ratio_file = if (gwas_software == "saige") then saige_step_1.variance_ratio_file else make_fake_saige_step_1.variance_ratio_file
    }
}