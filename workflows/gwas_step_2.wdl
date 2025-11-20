version 1.0

import "structs.wdl"


task regenie_step_2 {
    input {
        String docker_image
        String group_name
        String chr
        File pgen_file
        File pvar_file
        File psam_file
        File sample_list
        File covariates_file
        Array[File]? regenie_loco
        File? regenie_list
        Array[String] covariates_list
        String npcs = "10"
        String bsize = "1000"
        String mac = "10"
        String loco_pca
    }
    
    String chr_out = if (loco_pca == "false") then "0" else "~{chr}"

    command <<<

        for file in  ~{sep=" " regenie_loco}; do
            ln -s "$file" .
        done

        input_prefix=$(dirname "~{pgen_file}")/$(basename "~{pgen_file}" .pgen)

        # Only retain bi-allelic (REGENIE can't handle non-biallelic)
        plink2 \
            --pfile ${input_prefix} \
            --keep ~{sample_list} \
            --max-alleles 2 \
            --min-alleles 2 \
            --rm-dup exclude-all list \
            --make-pgen \
            --out ${input_prefix}.biallelic_frequent
            
        # Make covariates list
        pc_list=$(printf "CHR~{chr_out}_OUT_PC%s," {1..~{npcs}} | sed 's/,$//')
        full_covariates_list="~{sep="," covariates_list},${pc_list}"

        # phenotype from group_name
        phenotype=$(echo ~{group_name} | cut -d'.' -f2)
        echo $phenotype > phenotype.txt

        # Find the type of the phenotype (quantitative or binary)
        phenotype_col_idx=$(head -1 ~{covariates_file} | tr '\t' '\n' | grep -n ${phenotype} | cut -d: -f1)
        uniq_vals_count=$(cut -f"${phenotype_col_idx}" ~{covariates_file} | sort -u | grep -v "NA" | wc -l)
        trait_type="--bt"
        if [ "${uniq_vals_count}" -gt 3 ]; then # two binary values + header
            trait_type="--qt"
        fi

        conda run -n regenie_env regenie \
            --step 2 \
            --pgen ${input_prefix}.biallelic_frequent \
            --keep ~{sample_list} \
            --phenoFile ~{covariates_file} \
            --phenoColList ${phenotype} \
            --write-samples \
            --covarFile ~{covariates_file} \
            --covarColList ${full_covariates_list} \
            ${trait_type} \
            --firth --approx --pThresh 0.01 \
            --minMAC ~{mac} \
            --pred ~{regenie_list} \
            --bsize ~{bsize} \
            --out ~{group_name}.chr~{chr}
    >>>

    output {
        File summary_stats = "${group_name}.chr${chr}_" + read_string("phenotype.txt") + ".regenie"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }
}

workflow gwas_step_2 {
    input {
        String docker_image
        String gwas_software
        String group_name
        PGENFileset imputed_chr_fileset
        File sample_list
        File covariates_file
        Array[String] covariates_list
        Array[File]? regenie_loco_preds
        File? regenie_list
        String regenie_bsize
        String mac
        String npcs
        String loco_pca
    }

    if (gwas_software == "regenie") {
        call regenie_step_2 {
            input:
                docker_image = docker_image,
                group_name = group_name,
                chr = imputed_chr_fileset.chr,
                pgen_file = imputed_chr_fileset.pgen,
                pvar_file = imputed_chr_fileset.pvar,
                psam_file = imputed_chr_fileset.psam,
                sample_list = sample_list,
                covariates_file = covariates_file,
                regenie_loco = regenie_loco_preds,
                regenie_list = regenie_list,
                covariates_list = covariates_list,
                npcs = npcs,
                bsize = regenie_bsize,
                mac = mac,
                loco_pca = loco_pca
        }
    }

    output {
        File? gwas_output = regenie_step_2.summary_stats
    }

}