version 1.0

task make_gwas_outputs {
    input {
        String docker_image
        String output_prefix = "results.all_chr"
        String julia_cmd
        String maf
        Array[File] results_files
    }

    command <<<
        for f in ~{sep=" " results_files}; do
            echo "${f}"
        done > merge_list.txt

        ~{julia_cmd} make-gwas-outputs \
            merge_list.txt \
            --maf=~{maf} \
            --output-prefix=~{output_prefix}
    >>>

    output {
        File merged_results = "${output_prefix}.tsv"
        Array[File] plots = glob("*.png")
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x8"
        cpu: "8"
        memory: "32G"
        disks: "local-disk 100 SSD"
    }
}

task make_finemapping_outputs {
    input {
        String docker_image
        String output_prefix = "results.all_chr"
        String julia_cmd
        File gwas_results
        Array[File] results_files
    }

    command <<<
        for f in ~{sep=" " results_files}; do
            echo "${f}"
        done > merge_list.txt

        ~{julia_cmd} make-finemapping-outputs \
            merge_list.txt \
            ~{gwas_results} \
            --output-prefix=~{output_prefix}
    >>>

    output {
        File merged_results = "${output_prefix}.tsv"
        Array[File] plots = glob("*.png")
    }

    runtime {
        docker: docker_image
        dx_instance_type: "mem2_ssd1_v2_x8"
        cpu: "8"
        memory: "32G"
        disks: "local-disk 100 SSD"
    }
}