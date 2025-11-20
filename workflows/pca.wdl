version 1.0

import "structs.wdl"

task pca {
    input {
        String docker_image
        File bed_file
        File bim_file
        File fam_file
        String npcs = 10
        String approx = "true"
        String not_chr
    }

    String output_prefix = "pca." + basename(bed_file, ".ldpruned.bed") + ".chr~{not_chr}_out"

    command <<<
        genotypes_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)

        # Approx PCA option
        approx_option=""
        if [[ "~{approx}" == "true" ]]; then
            approx_option=" approx"
        fi

        # LOCO-PCA option
        loco_option="--not-chr ~{not_chr} "
        if [[ "~{not_chr}" == "0" ]]; then
            loco_option=""
        fi

        plink2 \
            --bfile ${genotypes_prefix} ${loco_option}\
            --pca ~{npcs}${approx_option} \
            --out ~{output_prefix}
    >>>

    output {
        File eigenvec = "${output_prefix}.eigenvec"
        File eigenval = "${output_prefix}.eigenval"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }
}

workflow compute_pcs {
    input {
        String docker_image
        Array[PGENFileset] imputed_genotypes
        File bed_file
        File bim_file
        File fam_file
        String loco_pca
        String npcs
        String approx_pca
    }

    if (loco_pca == "true") {
        scatter (imputed_chr_fileset in imputed_genotypes) {
            call pca as pca_loco {
                input:
                    docker_image = docker_image,
                    not_chr = imputed_chr_fileset.chr,
                    bed_file = bed_file,
                    bim_file = bim_file,
                    fam_file = fam_file,
                    npcs = npcs,
                    approx = approx_pca
            }
        }
    }

    if (loco_pca != "true") {
        # Fake loop to produce an array wven though a single file is created
        scatter (dummy in [1]) {
            call pca as pca_nonloco {
                input:
                    docker_image = docker_image,
                    not_chr = "0",
                    bed_file = bed_file,
                    bim_file = bim_file,
                    fam_file = fam_file,
                    npcs = npcs,
                    approx = approx_pca
            }
        }
    }

    output {
        Array[File]? eigenvec = if (loco_pca == "true") then pca_loco.eigenvec else pca_nonloco.eigenvec
    }
}