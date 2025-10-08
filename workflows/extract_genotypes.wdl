version 1.0

struct PGENFileset {
    String chr
    File pgen
    File psam
    File pvar
}

struct BGENFileset {
    String chr
    File bgen
    File bgi
    File sample
    File vcf_info
    File vcf_info_index
}

struct PLINKFileset {
    String chr
    File bed
    File bim
    File fam
}

workflow extract_ukb_genotypes {
    # Inputs
    input {
        String docker_image = "olivierlabayle/wdl-gwas:main"

        Array[BGENFileset] bgen_filesets

        String qc_genotype_missing_rate = "0.02"
        String qc_individual_missing_rate = "0.02"
        String r2_threshold = "0.9"
    }

    # Filter UKB chromosomes with R2 and critical samples
    scatter (bgen_fileset in bgen_filesets) {
        call make_pgen_and_bed {
            input:
                docker_image = docker_image,
                chr = bgen_fileset.chr,
                bgen_file = bgen_fileset.bgen,
                bgen_bgi_file = bgen_fileset.bgi,
                bgen_sample_file = bgen_fileset.sample,
                vcf_info_file = bgen_fileset.vcf_info,
                vcf_info_file_index = bgen_fileset.vcf_info_index,
                qc_genotype_missing_rate = qc_genotype_missing_rate,
                qc_individual_missing_rate = qc_individual_missing_rate,
                r2_threshold = r2_threshold
        }
    }

    # Merge typed PLINK files across chromosomes
    scatter (set in make_pgen_and_bed.plink_fileset) {
        File bed_files = set.bed
    }

    call merge_ukb_chrs {
        input:
            docker_image = docker_image,
            plink_filesets = make_pgen_and_bed.plink_fileset,
            bed_files = bed_files
    }

    # Outputs

    output {
        Array[PGENFileset] imputed_genotypes = make_pgen_and_bed.pgen_fileset
        PLINKFileset typed_genotypes = merge_ukb_chrs.plink_fileset
    }
}

task merge_ukb_chrs {
    input {
        String docker_image
        Array[PLINKFileset] plink_filesets
        Array[File] bed_files
    }

    command <<<
        # Merge filesets
        for f in ~{sep=" " bed_files}; do
            echo "${f%.bed}"
        done > merge_list.txt

        plink \
            --biallelic-only \
            --merge-list merge_list.txt \
            --make-bed \
            --out typed
    >>>

    output {
        PLINKFileset plink_fileset = object {
            chr: "all",
            bed: "typed.bed",
            bim: "typed.bim",
            fam: "typed.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd2_v2_x4"
    }
}

task make_pgen_and_bed {
    input {
        String docker_image
        String chr
        File bgen_file
        File bgen_bgi_file
        File bgen_sample_file
        File vcf_info_file
        File vcf_info_file_index
        String qc_genotype_missing_rate
        String qc_individual_missing_rate
        String r2_threshold = "0.9"
    }

    String input_prefix = basename(bgen_file, ".bgen")

    command <<<
        # 1. Make PGEN file
        ## Get variants with R2 > r2_threshold
        bcftools view -i 'INFO/R2>~{r2_threshold}' ~{vcf_info_file} | bcftools query -f '%CHROM\t%POS\t%ID\n' > imputed_variants.tsv
        awk 'BEGIN { OFS="\t" } { print $1, $2, $2 }' imputed_variants.tsv > imputed_extract_list.txt
        ## Filter BGEN file and write PGEN with defined filters
        plink2 \
            --bgen ~{bgen_file} ref-first \
            --sample ~{bgen_sample_file} \
            --extract range imputed_extract_list.txt \
            --max-alleles 2 \
            --geno ~{qc_genotype_missing_rate} \
            --mind ~{qc_individual_missing_rate} \
            --make-pgen \
            --out "~{input_prefix}.imputed"
        ## Update variant IDs in pvar
        awk 'BEGIN {OFS="\t"} NR==FNR {id[$2]=$3; next} {if (($2) in id) $3=id[$2]; print}' imputed_variants.tsv "~{input_prefix}.imputed.pvar" > temp_pvar.txt
        mv temp_pvar.txt "~{input_prefix}.imputed.pvar"
        
        # 2. Make a bed file with typed variants
        ## Extract typed variants
        bcftools view -i 'INFO/TYPED=1' ~{vcf_info_file} | bcftools query -f '%CHROM\t%POS\t%ID\n' > typed_variants.tsv
        awk 'BEGIN { OFS="\t" } { print $1, $2, $2 }' typed_variants.tsv > typed_extract_list.txt
        ## Create PLINK bed file
        plink2 \
            --bgen ~{bgen_file} ref-first \
            --sample ~{bgen_sample_file} \
            --extract range typed_extract_list.txt \
            --max-alleles 2 \
            --geno ~{qc_genotype_missing_rate} \
            --mind ~{qc_individual_missing_rate} \
            --make-bed \
            --out "~{input_prefix}.typed"
        ## Update variant IDs in bim
        awk 'BEGIN {OFS="\t"} NR==FNR {id[$2]=$3; next} {if (($4) in id) $2=id[$4]; print}' typed_variants.tsv "~{input_prefix}.typed.bim" > temp_bim.txt
        mv temp_bim.txt "~{input_prefix}.typed.bim"
    >>>

    output {
        PGENFileset pgen_fileset = object {
            chr: chr,
            pgen: "${input_prefix}.imputed.pgen",
            pvar: "${input_prefix}.imputed.pvar",
            psam: "${input_prefix}.imputed.psam"
        }
        PLINKFileset plink_fileset = object {
            chr: chr,
            bed: "${input_prefix}.typed.bed",
            bim: "${input_prefix}.typed.bim",
            fam: "${input_prefix}.typed.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }

}