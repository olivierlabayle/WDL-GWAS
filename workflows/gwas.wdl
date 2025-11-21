version 1.0

import "structs.wdl"
import "pca.wdl"
import "gwas_step_1.wdl"
import "gwas_step_2.wdl"

workflow gwas {
    input {
        String docker_image = "olivierlabayle/wdl-gwas:main"
        File covariates_file
        PLINKFileset genotypes
        Array[PGENFileset]+ imputed_genotypes
        Array[String] groupby = []
        Array[String] filterby = []
        Array[String] covariates = ["AGE", "SEX", "AGE_x_AGE", "AGE_x_SEX"]
        Array[String] phenotypes = ["SEVERE_COVID_19"]
        String julia_use_sysimage = "true"
        String julia_threads = "auto"
        # QC parameters
        String min_cases_controls = "10"
        String npcs = "10"
        String approx_pca = "true"
        String maf = "0.01"
        String mac = "10"
        String ip_values = "1000 50 0.05"
        # PCA parameters
        String loco_pca = "false"
        # GWAS software
        String gwas_software = "regenie"
        # Regenie parameters
        String regenie_cv_folds = "loocv" # or an integer
        String regenie_bsize = "1000"
        # Finemapping parameters
        String min_sig_clump_size = "10"
        String lead_pvalue = "5e-8"
        String p2_pvalue = "5e-5"
        String r2_threshold = "0.1"
        String clump_kb = "250"
        String n_causal = "10"
        String finemap_strategy = "rss"
        String susie_max_iter = "1000"
        # Meta analysis
        Boolean meta_analysis = length(groupby) > 0
        Array[String] meta_exclude = ["ADMIXED"]
        String meta_method = "STDERR"
    }

    # Get generic Julia command
    call get_julia_cmd as get_julia_cmd {
        input:
            use_sysimage = julia_use_sysimage,
            threads = julia_threads
    }

    # Create groups and update covariates
    call make_groups_and_covariates {
        input:
            docker_image=docker_image,
            covariates_file=covariates_file,
            groupby=groupby,
            filters=filterby,
            covariates=covariates,
            phenotypes_list=phenotypes,
            min_cases_controls=min_cases_controls,
            julia_cmd=get_julia_cmd.julia_cmd
    }


    # A GWAS is performed for each group defined by a combination of the user provided `groupby` and a phenotype
    scatter (sample_list in make_groups_and_covariates.groups_individuals) {
        String group_name = sub(basename(sample_list, ".txt"), "gwas.individuals.", "")

        # Make plink bed files for the group 
        call make_group_bed_qced {
            input:
                docker_image = docker_image,
                chr = genotypes.chr,
                genotypes_bed = genotypes.bed,
                genotypes_bim = genotypes.bim,
                genotypes_fam = genotypes.fam,
                sample_list = sample_list,
                maf = maf,
                mac = mac
        }

        # LD prune the bed file for LOCO PCA
        call ld_prune as groups_ld_prune {
            input:
                docker_image = docker_image,
                chr = make_group_bed_qced.plink_fileset.chr,
                bed_file = make_group_bed_qced.plink_fileset.bed,
                bim_file = make_group_bed_qced.plink_fileset.bim,
                fam_file = make_group_bed_qced.plink_fileset.fam,
                output_prefix = basename(make_group_bed_qced.plink_fileset.bed, ".bed") + ".ldpruned",
                ip_values = ip_values,
                maf = maf
        }

        # Perform LOCO PCA
        call pca.compute_pcs as compute_pcs {
            input:
                docker_image = docker_image,
                imputed_genotypes = imputed_genotypes,
                bed_file = groups_ld_prune.ld_pruned_fileset.bed,
                bim_file = groups_ld_prune.ld_pruned_fileset.bim,
                fam_file = groups_ld_prune.ld_pruned_fileset.fam,
                loco_pca = loco_pca,
                npcs = npcs,
                approx_pca = approx_pca
        }

        # Merge covariates and PCs
        call merge_covariates_and_pcs {
            input:
                docker_image = docker_image,
                group_name = group_name,
                covariates_file = make_groups_and_covariates.updated_covariates,
                pcs_files = select_first([compute_pcs.eigenvec]),
                julia_cmd = get_julia_cmd.julia_cmd 
        }

        # First run regenie step 1 using plink genotypes
        call gwas_step_1.gwas_step_1 as run_gwas_step_1 {
            input:
                docker_image = docker_image,
                gwas_software = gwas_software,
                group_name = group_name,
                bed_file = make_group_bed_qced.plink_fileset.bed,
                bim_file = make_group_bed_qced.plink_fileset.bim,
                fam_file = make_group_bed_qced.plink_fileset.fam,
                sample_list = sample_list,
                covariates_file = merge_covariates_and_pcs.covariates_and_pcs,
                covariates_list = make_groups_and_covariates.covariates_list,
                cv_folds = regenie_cv_folds,
                bsize = regenie_bsize,
                maf = maf,
                mac = mac
        }

        # Second run regenie step 2 across imputed chromosomes filesets
        scatter (imputed_chr_fileset in imputed_genotypes) {
            # Perform Regenie Step 2
            call gwas_step_2.gwas_step_2 as run_gwas_step_2 {
                input:
                    docker_image = docker_image,
                    julia_cmd = get_julia_cmd.julia_cmd,
                    gwas_software = gwas_software,
                    group_name = group_name,
                    imputed_chr_fileset = imputed_chr_fileset,
                    sample_list = sample_list,
                    covariates_file = merge_covariates_and_pcs.covariates_and_pcs,
                    covariates_list = make_groups_and_covariates.covariates_list,
                    regenie_loco_preds = select_first([run_gwas_step_1.regenie_loco_preds]),
                    regenie_list = run_gwas_step_1.regenie_list,
                    saige_model_file = run_gwas_step_1.saige_model_file,
                    saige_variance_ratio_file = run_gwas_step_1.saige_variance_ratio_file,
                    regenie_bsize = regenie_bsize,
                    mac = mac,
                    npcs = npcs,
                    loco_pca = loco_pca
            }

            # Finemap results for the chromosome
            call finemapping {
                input:
                    docker_image = docker_image,
                    julia_cmd = get_julia_cmd.julia_cmd,
                    gwas_results = run_gwas_step_2.gwas_output,
                    pgen_file = imputed_chr_fileset.pgen,
                    pvar_file = imputed_chr_fileset.pvar,
                    psam_file = imputed_chr_fileset.psam,
                    chr = imputed_chr_fileset.chr,
                    covariates_file = merge_covariates_and_pcs.covariates_and_pcs,
                    sample_file = sample_list,
                    group_name = group_name,
                    min_sig_clump_size = min_sig_clump_size,
                    lead_pvalue = lead_pvalue,
                    p2_pvalue = p2_pvalue,
                    r2_threshold = r2_threshold,
                    clump_kb = clump_kb,
                    n_causal = n_causal,
                    susie_max_iter = susie_max_iter,
                    finemap_strategy = finemap_strategy
            }
        }

        # Merge GWAS results across chromosomes
        call merge_chr_results as merge_gwas_group_chr_results {
            input:
                docker_image = docker_image,
                output_prefix = group_name + ".gwas",
                julia_cmd = get_julia_cmd.julia_cmd,
                results_files = select_all(run_gwas_step_2.gwas_output)
        }

        # Merge Finemapping results across chromosomes
        call merge_chr_results as merge_fp_group_chr_results {
            input:
                docker_image = docker_image,
                output_prefix = group_name + ".finemapping",
                julia_cmd = get_julia_cmd.julia_cmd,
                results_files = finemapping.finemapping_results
        }

        # Generate GWAS plots
        call make_plots as gwas_group_plots {
            input:
                docker_image = docker_image,
                julia_cmd = get_julia_cmd.julia_cmd,
                gwas_results = merge_gwas_group_chr_results.merged_results,
                finemapping_results = merge_fp_group_chr_results.merged_results,
                maf = maf,
                output_prefix = group_name
        }
    }

    # Meta Analyse
    if (meta_analysis) {
        call meta_analyse {
            input:
                docker_image = docker_image,
                julia_cmd = get_julia_cmd.julia_cmd,
                gwas_results = merge_gwas_group_chr_results.merged_results,
                exclude = meta_exclude,
                method = meta_method
        }

        # Finemap meta-analysed results: this effectively loops through phenotypes
        scatter (meta_gwas_result in meta_analyse.meta_gwas_results) {
            String phenotype = sub(basename(meta_gwas_result, ".gwas.tsv"), "META_ANALYSIS.", "")

            scatter (imputed_chr_fileset in imputed_genotypes) {
                call finemapping_summary_stats {
                    input:
                        docker_image = docker_image,
                        julia_cmd = get_julia_cmd.julia_cmd,
                        gwas_results = meta_gwas_result,
                        covariates_file = make_groups_and_covariates.updated_covariates,
                        sample_files = make_groups_and_covariates.groups_individuals,
                        pgen_file = imputed_chr_fileset.pgen,
                        pvar_file = imputed_chr_fileset.pvar,
                        psam_file = imputed_chr_fileset.psam,
                        chr = imputed_chr_fileset.chr,
                        min_sig_clump_size = min_sig_clump_size,
                        lead_pvalue = lead_pvalue,
                        p2_pvalue = p2_pvalue,
                        r2_threshold = r2_threshold,
                        clump_kb = clump_kb,
                        n_causal = n_causal,
                        exclude = meta_exclude,
                        susie_max_iter = susie_max_iter,
                        phenotype=phenotype
                }
            }

            call merge_chr_results as merge_fp_meta_chr_results {
                input:
                    docker_image = docker_image,
                    output_prefix = "META_ANALYSIS." + phenotype + ".finemapping" ,
                    julia_cmd = get_julia_cmd.julia_cmd,
                    results_files = finemapping_summary_stats.finemapping_results
            }

            call make_plots as gwas_meta_plots {
                input:
                    docker_image = docker_image,
                    julia_cmd = get_julia_cmd.julia_cmd,
                    gwas_results = meta_gwas_result,
                    finemapping_results = merge_fp_meta_chr_results.merged_results,
                    maf = maf,
                    output_prefix = "META_ANALYSIS." + phenotype
            }
        }

        Array[File] phenotypes_meta_plots = flatten(gwas_meta_plots.plots)
    }

    output {
        Array[File] gwas_group_results = merge_gwas_group_chr_results.merged_results
        Array[File] finemapping_group_results = merge_fp_group_chr_results.merged_results
        Array[File] group_plots = flatten(gwas_group_plots.plots)
        Array[File]? meta_gwas_results = meta_analyse.meta_gwas_results
        Array[File]? meta_finemapping_results = merge_fp_meta_chr_results.merged_results
        Array[File]? meta_plots = phenotypes_meta_plots
    }
}


task get_julia_cmd {
    input {
        String use_sysimage = "true"
        String threads = "auto"
    }
    command <<<
        julia_cmd_string="julia --project=/opt/PopGen --startup-file=no"
        if [[ "~{use_sysimage}" == "true" ]]; then
            julia_cmd_string+=" --sysimage=/opt/PopGen/sysimage.so"
        fi
        if [[ "~{threads}" == "auto" ]]; then
            julia_cmd_string+=" --threads=auto"
        fi
        julia_cmd_string+=" /opt/PopGen/bin/wdl-gwas.jl"
        echo "$julia_cmd_string"
    >>>

    output {
        String julia_cmd = read_string(stdout())
    }
}

task ld_prune {
    input {
        String docker_image
        String chr
        File bed_file
        File bim_file
        File fam_file
        String output_prefix = "ld_pruned"
        String ip_values = "1000 50 0.05"
        String maf = "0.01"
    }

    command <<<
        bed_prefix=$(dirname "~{bed_file}")/$(basename "~{bed_file}" .bed)

        plink2 \
            --bfile ${bed_prefix} \
            --indep-pairwise ~{ip_values}
        # Always exclude high LD regions stored in the docker image at /opt/PopGen/assets/exclude_b38.txt
        plink2 \
            --bfile ${bed_prefix} \
            --extract plink2.prune.in \
            --maf ~{maf} \
            --make-bed \
            --exclude range /opt/PopGen/assets/exclude_b38.txt \
            --out ~{output_prefix}
    >>>

    output {
        PLINKFileset ld_pruned_fileset = object {
            chr: "all",
            bed: "${output_prefix}.bed",
            bim: "${output_prefix}.bim",
            fam: "${output_prefix}.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem1_ssd1_v2_x8"
    }
}

task make_plots {
    input {
        String docker_image
        String julia_cmd
        File gwas_results
        File finemapping_results
        String maf = "0.01"
        String output_prefix
    }

    command <<<
        ~{julia_cmd} make-plots \
            ~{gwas_results} \
            ~{finemapping_results} \
            --maf=~{maf} \
            --output-prefix=~{output_prefix}
    >>>

    output {
        Array[File] plots = glob("*.png")
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x8"
    }
}

task meta_analyse {
    input {
        String docker_image
        String julia_cmd
        Array[File] gwas_results
        Array[String] exclude
        String method = "STDERR"
    }

    command <<<
        for f in ~{sep=" " gwas_results}; do
            echo "${f}"
        done > gwas_meta_list.txt

        ~{julia_cmd} meta-analyse \
            gwas_meta_list.txt \
            --exclude=~{sep="," exclude} \
            --method=~{method} \
            --output-prefix=META_ANALYSIS
    >>>

    output {
        Array[File] meta_gwas_results = glob("META_ANALYSIS.*.tsv")
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x8"
    }
}

task merge_chr_results {
    input {
        String docker_image
        String output_prefix = "results.all_chr"
        String julia_cmd
        Array[File] results_files
    }

    command <<<
        for f in ~{sep=" " results_files}; do
            echo "${f}"
        done > merge_list.txt

        ~{julia_cmd} merge-chr-results \
            merge_list.txt \
            --output-prefix=~{output_prefix}
    >>>

    output {
        File merged_results = "${output_prefix}.tsv"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x8"
    }
}

task finemapping_summary_stats {
    input {
        String docker_image
        String julia_cmd
        File gwas_results
        File pgen_file
        File pvar_file
        File psam_file
        String chr
        File covariates_file
        Array[File] sample_files
        String min_sig_clump_size = "3"
        String lead_pvalue = "5e-8"
        String p2_pvalue = "1e-5"
        String r2_threshold = "0.1"
        String clump_kb = "1000"
        String n_causal = "10"
        String susie_max_iter = "1000"
        Array[String] exclude = []
        String phenotype
    }

    command <<<
        # Make samples files list
        for f in ~{sep=" " sample_files}; do
            echo "${f}"
        done > samples_file.txt
        # Get PGEN prefix
        pgen_prefix=$(dirname "~{pgen_file}")/$(basename "~{pgen_file}" .pgen)
        # Run finemapping
        ~{julia_cmd} finemap \
            ~{gwas_results} \
            ${pgen_prefix} \
            ~{covariates_file} \
            samples_file.txt \
            --output-prefix=finemapping.meta_analysis.~{phenotype}.chr~{chr} \
            --min-sig-clump-size=~{min_sig_clump_size} \
            --lead-pvalue=~{lead_pvalue} \
            --p2-pvalue=~{p2_pvalue} \
            --r2-threshold=~{r2_threshold} \
            --clump-kb=~{clump_kb} \
            --n-causal=~{n_causal} \
            --phenotype=~{phenotype} \
            --susie-max-iter=~{susie_max_iter} \
            --rss \
            --exclude=~{sep="," exclude}
    >>>

    output {
        File finemapping_results = "finemapping.meta_analysis.${phenotype}.chr${chr}.tsv"
        File clumping_results = "finemapping.meta_analysis.${phenotype}.chr${chr}.clumps.tsv"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }
}


task finemapping {
    input {
        String docker_image
        String julia_cmd
        File? gwas_results
        File pgen_file
        File pvar_file
        File psam_file
        String chr
        File covariates_file
        File sample_file
        String group_name
        String min_sig_clump_size = "3"
        String lead_pvalue = "5e-8"
        String p2_pvalue = "1e-5"
        String r2_threshold = "0.1"
        String clump_kb = "1000"
        String n_causal = "10"
        String finemap_strategy = "rss" # or "full"
        String susie_max_iter = "1000"
    }

    command <<<
        phenotype=$(echo ~{group_name} | cut -d'.' -f2)
        pgen_prefix=$(dirname "~{pgen_file}")/$(basename "~{pgen_file}" .pgen)

        finemap_opt=""
        if [[ "~{finemap_strategy}" == "rss" ]]; then
            finemap_opt="--rss"
        fi

        ~{julia_cmd} finemap \
            ~{gwas_results} \
            ${pgen_prefix} \
            ~{covariates_file} \
            ~{sample_file} \
            --output-prefix=finemapping.~{group_name}.chr~{chr} \
            --min-sig-clump-size=~{min_sig_clump_size} \
            --lead-pvalue=~{lead_pvalue} \
            --p2-pvalue=~{p2_pvalue} \
            --r2-threshold=~{r2_threshold} \
            --clump-kb=~{clump_kb} \
            --n-causal=~{n_causal} \
            --susie-max-iter=~{susie_max_iter} \
            --phenotype=${phenotype} ${finemap_opt}
    >>>

    output {
        File finemapping_results = "finemapping.${group_name}.chr${chr}.tsv"
        File clumping_results = "finemapping.${group_name}.chr${chr}.clumps.tsv"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }
}

task merge_covariates_and_pcs {
    input {
        String docker_image
        String group_name
        File covariates_file
        Array[File] pcs_files
        String julia_cmd
    }

    String outfile = group_name + ".merged_covariates_and_pcs.tsv"

    command <<<
        for file in  ~{sep=" " pcs_files}; do
            ln -s "$file" .
        done

        ~{julia_cmd} \
            merge-covariates-pcs \
            ~{covariates_file} \
            pca \
            --output=~{outfile}
    >>>

    output {
        File covariates_and_pcs = "${outfile}"
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x16"
    }
}

task make_groups_and_covariates {
    input {
        String docker_image
        File covariates_file
        Array[String] groupby = []
        Array[String] filters = []
        Array[String] covariates = ["SEX", "AGE"]
        Array[String] phenotypes_list = ["SEVERE_COVID_19"]
        String min_cases_controls = "10"
        String julia_cmd
    }

    command <<<
        groupby_string='~{sep="," groupby}'
        groupby_string_opt=""
        if [[ -n "${groupby_string}" ]]; then
            groupby_string_opt="--groupby=${groupby_string}"
        fi

        filters_string='~{sep="," filters}'
        filters_string_opt=""
        if [[ -n "${filters_string}" ]]; then
            filters_string_opt="--filters=${filters_string}"
        fi

        covariates_string='~{sep="," covariates}'

        ~{julia_cmd} \
            make-groups-and-covariates \
            ~{covariates_file} \
            --covariates=${covariates_string} \
            --phenotypes=~{sep="," phenotypes_list} \
            --output-prefix=gwas \
            --min-cases-controls=~{min_cases_controls} ${groupby_string_opt} ${filters_string_opt}
    >>>

    output {
        File updated_covariates = "gwas.covariates.csv"
        Array[String] covariates_list = read_lines("gwas.covariates_list.txt")
        Array[File]+ groups_individuals = glob("gwas.individuals.*")
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x8"
    }
}

task make_group_bed_qced {
    input {
        String docker_image
        String chr
        File genotypes_bed
        File genotypes_bim
        File genotypes_fam
        File sample_list
        String maf = "0.01"
        String mac = "10"
    }

    String group_name = sub(basename(sample_list, ".txt"), "gwas.individuals.", "")

    command <<<
        genotypes_prefix=$(dirname "~{genotypes_bed}")/$(basename "~{genotypes_bed}" .bed)

        plink2 \
            --bfile ${genotypes_prefix} \
            --keep ~{sample_list} \
            --maf ~{maf} \
            --mac ~{mac} \
            --make-bed \
            --out ~{group_name}
    >>>

    output {
        PLINKFileset plink_fileset = object {
            chr: chr,
            bed: "${group_name}.bed",
            bim: "${group_name}.bim",
            fam: "${group_name}.fam"
        }
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x8"
    }
}