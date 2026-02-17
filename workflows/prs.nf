#! /appl/nextflow-23.04.1.5866/bin/nextflow

// enable DSL2 syntax***
nextflow.enable.dsl = 2

workflow {
    /*
    this section is called first if:
        Apply PGS is the first workflow called
        Apply PGS is called independently of other workflows
        Another PRScsx workflow is not being stitched to Apply PGS workflow
        Channels from another workflow are not inputs to the Apply PGS workflow
    */
    score_weights = PGS_setup()
    indiv_polygenic_scores = PGS(score_weights)
}

workflow PGS_setup {
    main:
        // make summary table file name into a string type
        descriptor_file = new File(params.input_descriptor_table_filename.toString())
        // reads summary table file into groovy making it into a list
        descriptor_lines = descriptor_file.readLines()
        if (descriptor_file.text.tokenize('\n').size() <= 2) {
            descriptor_lines[1..-1].each {
                line ->
                /* groovylint-disable-next-line GStringExpressionWithinString */
                lineG = line.replace('${launchDir}', "${launchDir}")
                line_parts = lineG.toString().trim().replace('[', '').replace(']', '').split(',') as List
                descriptor_tuple = line_parts
            }
        } else {
            descriptor_tuple = []
            descriptor_lines[1..-1].each {
                line ->
                /* groovylint-disable-next-line GStringExpressionWithinString */
                lineG = line.replace('${launchDir}', "${launchDir}")
                line_parts = lineG.toString().trim().replace('[', '').replace(']', '').split(',') as List
                descriptor_tuple.add(line_parts)
            }
        }
        // creates a channel out of the summary table list
        descriptor_channel = Channel.fromList(descriptor_tuple)
        // parse header and retrieve indices of essential columns based on their column names
        descriptor_header = descriptor_file.readLines().get(0).split(',')
        descriptor_header.eachWithIndex { item, index ->
            if (item == params['input_descriptor_table_colnames']['id_colname']) {
                id_index = index
            }
            if (item == params['input_descriptor_table_colnames']['pgs_weights_full_filpath_colname']) {
                score_file_index = index
            }
        }
        // find out indices
        // makes each item (table row) in the list a separate tuple
        if (descriptor_file.text.tokenize('\n').size() <= 2) {
            descriptor_channel_collect = descriptor_channel.collect()
            score_weights = descriptor_channel_collect.map { arr -> new Tuple(
                arr[id_index],
                arr[score_file_index]) }
        } else {
            score_weights = descriptor_channel.map { arr -> new Tuple(
                arr[id_index],
                arr[score_file_index]) }
        }
    emit:
        score_weights
}
workflow PGS {
    take:
        score_weights
    main:
        // create validation cohort list
        validation_population_list = params.validation_populations.keySet().toList()

        // create validation_population, and chromomosome channels to parallelize by
        validation_population = Channel.fromList(validation_population_list)
        chromosome = Channel.fromList(params.chromosome_list)

        // clean vcf process
        // define plink2 object
        my_plink2 = "${params.my_plink2}"
        // combine val pop channel with chromosome channel
        val_pop_vcf = validation_population.combine(chromosome).map { val_pop, chr ->
        new Tuple(
            val_pop,
            chr,
            "${params.validation_populations[val_pop]['vcf_prefix']}${chr}${params.validation_populations[val_pop]['vcf_suffix']}.vcf"
        ) }
        // call process
        val_pop_plink = clean_vcf(val_pop_vcf, my_plink)

        // make column names channels
        score_chr_colname = "${params['score_file_colnames']['chr_colname']}"
        score_pos_colname = "${params['score_file_colnames']['pos_colname']}"
        score_a1_colname = "${params['score_file_colnames']['a1_colname']}"
        score_a2_colname = "${params['score_file_colnames']['a2_colname']}"
        score_pgs_colname = "${params['score_file_colnames']['pgs_colname']}"
        // make script a channel
        plink2_score_input_script = "${launchDir}/scripts/pgs_weight_reformat.py"
        // define python path channel
        my_python = "${params.my_python}"
        // call process
        score_reformat = make_plink2_score_input(score_weights,
                                                plink2_score_input_script,
                                                score_chr_colname,
                                                score_pos_colname,
                                                score_a1_colname,
                                                score_a2_colname,
                                                score_pgs_colname,
                                                my_python)

        // join score files process
        // collect separate input files into one channel
        all_score_inputs = score_reformat.collect()
        // define script channel
        join_score_files_script = "${launchDir}/scripts/join_score_files.py"
        // call process
        combined_score_input = join_score_files(all_score_inputs, my_python, join_score_files_script)

        // intersect variants process
        // define script
        intersect_variants_script = "${launchDir}/scripts/intersect_variants.py"
        // extract pvar file from clean vcf process output
        val_pop_pvar = val_pop_plink.map { val_pop, chr, files ->
        def psam = files.find { it.name.endsWith(".pvar") }
        tuple(val_pop, psam)}.groupTuple(by: [0])
        // define ref panel pvar file
        ref_panel_pvar = "${params.ref_panel_plink_prefix}.pvar"
        // call process
        intersected_variants = intersect_variants(val_pop_pvar, ref_panel_pvar, combined_score_input, my_python, intersect_variants_script)

        // intersect variants target process
        // extract variant list from interesected variants output
        variant_list = intersected_variants.map{ val_pop, variant_list, score_file -> variant_list }
        // call process
        target_intersected = intersect_variants_target(val_pop_plink, variant_list, my_plink2)

        // clean_sample_list process
        // extract psam files from clean vcf process output
        val_pop_psam = val_pop_plink.map { val_pop, chr, files ->
        def psam = files.find { it.name.endsWith(".psam") }
        tuple(val_pop, psam)}.groupTuple(by: [0])
        // create tuple with validation population, sample list, and psam file
        val_pop_sample_list = validation_population.combine(chromosome).map { val_pop, chr ->
        new Tuple(
            val_pop,
            params.validation_populations[val_pop]['population_subset_file'],
            params.validation_populations[val_pop]['population_subset_file_id_col'],
            params.validation_populations[val_pop]['population_subset_file_delim'],
        ) }.groupTuple(by: [0])
        // combine sample list tuple with psam tuple
        clean_sample_list_input = val_pop_psam.join(val_pop_sample_list)
        // define script
        clean_sample_list_script = "${launchDir}/scripts/clean_sample_list.py"
        // call process
        cleaned_sample_list = clean_sample_list(clean_sample_list_input, my_python, clean_sample_list_script)

        // compute_scores process
        // join sample list tuple with plink files tuple
        compute_scores_input = target_intersected.first().join(cleaned_sample_list)
        // create --read-freq plink line
        plink_read_freq_line = params['read_freq'] == 'auto' ? '' : "--read-freq ${params.read_freq}"
        // make channels from other PLINK score parameters
        dosage_transformation = "${params.dosage_transformation}"
        xchr_model = "${params.xchr_model}"
        no_mean_imputation = "${params.no_mean_imputation}"
        independent_se = "${params.independent_se}"
        // call process
        plink_score_output = compute_scores(compute_scores_input,
                                            combine_apply_pgs_input,
                                            plink_read_freq_line,
                                            dosage_transformation,
                                            xchr_model,
                                            no_mean_imputation,
                                            independent_se,
                                            my_plink2)

        // concatenate_plink_score_outputs process
        // group chromosome separated files together
        plink_score_output_grouped = plink_score_output[0].groupTuple(by: 0, size: params.chromosome_list.size())
        // make cat outputs script channel
        cat_outputs_script = "${launchDir}/scripts/cat_plink_score_outputs.py"
        // call process
        cat_scores = concatenate_plink_score_outputs(plink_score_output_grouped,
                                                        cat_outputs_script,
                                                        my_python)

        // pca qc process
        // define ref psam and pgen and ref key
        ref_panel_pgen = "${params.ref_panel_plink_prefix}.pgen"
        ref_panel_psam = "${params.ref_panel_plink_prefix}.psam"
        ref_panel_key = "${params.ref_panel_key}"
        // create tuple
        ref_panel_tuple = tuple(ref_panel_key, ref_panel_pgen, ref_panel_psam, ref_panel_pvar).groupTuple(by : [0])
        // define qc thresholds
        ref_hwe = "${params.ref_hwe_threshold}"
        ref_maf = "${params.ref_maf_threshold}"
        ref_geno = "${params.ref_geno_threshold}"
        ref_mind = "${params.ref_mind_threshold}"
        ref_ld_prune = "${params.ref_ld_prune_params}"
        // call process
        ref_qc_output = pca_qc(ref_panel_tuple,
                                variant_list,
                                my_plink2,
                                ref_hwe,
                                ref_maf,
                                ref_geno,
                                ref_mind,
                                ref_ld_prune)

        // target pca qc process
        // extract plink files from intersect variants process
        target_pca_qc_input = target_intersected.first()
        // extract variant list from ref qc process
        ref_qc_variant_list = ref_qc_output.map { plink1_files, plink2_files, prunein, pruneout, log -> prunein}
        // call process
        target_pca_qc_output = target_pca_qc(target_pca_qc_input, ref_qc_variant_list, my_plink)

        // merge target files for pca process
        // make merge list
        merge_list = target_pca_qc_output.first().groupTuple(by: [0]).map { val_pop, plink_files ->
        def prefixes = records.collect { it[2] }
        tuple(pop, prefixes)}.map { val_pop, prefixes ->
        def list_file = file("${val_pop}.merge_list.txt")
        list_file.text = prefixes.join("\n")
        tuple(pop, list_file)}
        // join merge list and plink files
        target_merge_input = target_pca_qc_output.first().groupTuple(by: [0]).join(merge_list)
        // call process
        target_merge_output = target_merge_pca_qc(target_merge_input, my_plink2)

        // ref pca process
        // extract pruned plink1 files from qc output tuple
        ref_qc_plink1 = ref_qc.first()
        // define fraposa path
        my_fraposa = "${params.my_fraposa}"
        // call process
        ref_pca_output = ref_pca(ref_qc_plink1, my_fraposa)

        // project pca process
        // extract plink files from target_merge_output
        target_merge_output_plink = target_merge_output.first()
        // call process
        pca_projection = project_pca(ref_qc_plink1, target_merge_output_plink, my_fraposa)

        // ancestry adjustment process
        // define relateds
        related_indiv = validation_population.map { val_pop -> 
                        new Tuple(val_pop,
                        params.validation_populations[val_pop]['population_subset_file_delim'])}
        // define pgscatalog path
        my_pgscatalog = "${params.my_pgscatalog}"
        // extract ref plink psam
        ref_pca_clean_psam = ref_qc_plink1.map { ref_key, files ->
        def psam = files.find { it.name.endsWith(".psam") }
        tuple(ref_key, psam)}
        // extract concatenated score file
        cat_scores_filt = cat_scores.map{ val_pop, score, vars -> new Tuple(val_pop, score)}
        // extract pcs
        pcs = pca_projection.map{ val_pop, ref_key, pca.projection.pcs, pca.projection.ref.pcs, dat ->
                                    new Tuple(val_pop, ref_key, pca.projection.pcs, pca.projection.ref.pcs)}
        // join tuples
        ancestry_input = target_merge_output_plink.join(cat_scores_filt).join(related_indiv)
                tuple val(validation_population), path(ref_plink_psam), path(target_plink_files), path(ref_pcs), path(target_pcs)
        // call process
        ancestry_output = ancestry_adjustment(ancestry_input, ref_pca_clean_psam, pcs, my_pgscatalog)

        // pgs percentile process
        // define percentile cutoff from params
        percentile_cutoff = "${params.percentile_cutoff}"
        // extract adjusted score files from ancestry process
        adjusted_score_files = ancestry_output.first()
        // define percentile script
        percentile_script = "${launchDir}/scripts/pgs_percentile.py"
        // call process
        pgs_percentile_output = pgs_percentile(adjusted_score_files,
                                                percentile_script
                                                my_python,
                                                percentile_cutoff)

        // reformat output for emit
        indiv_polygenic_scores = pgs_percentile_output

        // emit params to json
        json_params = dump_params_to_json(params)

    emit:
       indiv_polygenic_scores
}

process clean_vcf {
    publishDir "${launchDir}/plink2/initial_cleaning/${validation_population}/"
    input:
        tuple val(validation_population), val(chromosome), path(vcf_file)
        val(my_plink2)
    output:
        tuple val(validation_population), val(chromosome), path("${validation_population}.cleaned.chr${chromosome}.{pgen, psam, pvar}")
        path(${validation_population}.cleaned.chr${chromosome}.log)
    shell:
        """
        # steps:
        # 1. update variant ID format
        # 2. convert vcf to PLINK2 format

        ${my_plink2} --vcf ${vcf_file} \
        --new-id-max-allele-len 500 \
        --set-all-var-ids chr@:#:\$r:\$a \
        --make-pgen \
        --out ${validation_population}.cleaned.chr${chromosome} > ${validation_population}.cleaned.chr${chromosome}.log
        """
    stub:
        """
        touch ${validation_population}.cleaned.chr${chromosome}.pgen
        touch ${validation_population}.cleaned.chr${chromosome}.psam
        touch ${validation_population}.cleaned.chr${chromosome}.pvar
        path("${validation_population}.cleaned.chr${chromosome}.log")
        """
}

process make_plink2_score_input {
    publishDir "${launchDir}/pgs_weights/reformat/"
    input:
        tuple val(pgs_id), path(pgs_weight_file)
        path(plink2_score_input_script)
        val(score_chr_colname)
        val(score_pos_colname)
        val(score_a1_colname)
        val(score_a2_colname)
        val(score_id_colname)
        val(score_pgs_colname)
        val(my_python)
    output:
        path("${pgs_id}.reformatted.txt")
    shell:
        """
        ${my_python} ${plink2_score_input_script} \
        --pgsId '${pgs_id}' \
        --scoreFile '${pgs_weight_file}' \
        --scoreChrCol '${score_chr_colname}' \
        --scorePosCol '${score_pos_colname}' \
        --scoreA1Col '${score_a1_colname}' \
        --scoreA2Col '${score_a2_colname}' \
        --scorePGSCol '${score_pgs_colname}' \
        """
    stub:
        """
        touch ${pgs_id}.reformatted.txt
        """
}

process join_score_files {
    publishDir "${launchDir}/pgs_weights/join/"
    input:
        path(all_score_files)
        val(my_python)
        path(join_score_files_script)
    output:
        path('combined_score_file.txt')
    script:
        """
        ${my_python} ${join_score_files_script} \
        --score_files ${all_score_files}
        """
    stub:
        """
        touch combined_score_file.txt
        """
}

process intersect_variants {
    publishDir "${launchDir}/intersect_variants/"
    input:
        tuple val(validation_population), path(pvar_file)
        path(ref_panel_file)
        path(combined_score_file)
        val(my_python)
        path(intersect_variants_pgs_script)
    output:
        tuple val(validation_population), path("${validation_population}.variant_list.intersect.txt"), path("${validation_population}.combined_score_file.intersect.txt")
    script:
        """
        ${my_python} ${intersect_variants_pgs_script} \
        --valPop ${validation_population} \
        --targetPvarFileList ${pvar_file} \
        --scoreFile ${combined_score_file} \
        --refPanelPvarFile ${ref_panel_file}
        """
    stub:
        """
        touch ${validation_population}.variant_list.intersect.txt
        touch ${validation_population}.combined_score_file.intersect.txt
        """
}

process intersect_variants_target {
    publishDir "${launchDir}/intersect_variants/target_plink/${validation_population}/"
    input:
        tuple val(validation_population), val(chromosome), path(plink_files)
        path(intersect_variants_list)
        val(my_plink2)
    output:
        tuple val(validation_population), val(chromosome), path("${validation_population}.intersect.chr${chromosome}.{pgen, psam, pvar}")
        path("${validation_population}.intersect.chr${chromosome}.log")
    shell:
        """
        # steps:
        # 1. get plink prefix
        # 2. filter to intersect variants

        # get plink prefix
        echo ${plink_files} | sed 's/[][]//g' | cut -d',' -f1 | sed -e 's/.psam//g' -e 's/.pvar//g' -e 's/.pgen//g' | sed 's/,//g' | sed 's/ //g' > temp
        plink_prefix=\$(cat "temp")

        # make pgen files for plink score
        ${my_plink2} \$plink_prefix \
        --extract ${intersect_variants_list} \
        --make-pgen \
        --out ${validation_population}.intersect.chr${chromosome} > ${validation_population}.intersect.chr${chromosome}.log
        """
    stub:
        """
        touch ${validation_population}.intersect.chr${chromosome}.pgen
        touch ${validation_population}.intersect.chr${chromosome}.psam
        touch ${validation_population}.intersect.chr${chromosome}.pvar
        touch ${validation_population}.intersect.chr${chromosome}.bed
        touch ${validation_population}.intersect.chr${chromosome}.bim
        touch ${validation_population}.intersect.chr${chromosome}.fam
        touch ${validation_population}.intersect.chr${chromosome}.log
        """
}

process clean_sample_list {
    publishDir "${launchDir}/sample_lists/"
    input:
        tuple val(validation_population), path(sample_list), val(id_col), val(sample_delim), path(first_psam_file)
        val(my_python)
        path(clean_sample_list_script)
    output:
        tuple val(validation_population), path("${validation_population}.fam_format_sample_list.txt")
    script:
        """
        ${my_python} ${clean_sample_list_script} \
        --valPop ${validation_population} \
        --sample_list ${sample_list} \
        --sample_delim ${sample_delim} \
        --sample_id_col ${sample_id_col} \
        --psam_file ${first_psam_file}
        """
    stub:
        """
        touch ${validation_population}.fam_format_sample_list.txt
        """
}

process compute_scores {
    publishDir "${launchDir}/plink_score_output/chromosome_separated_outputs/"
    input:
        tuple val(validation_population), val(chromosome), path(plink_files), path(sample_list_file)
        path(combined_score_input)
        val(plink_read_freq_line)
        val(dosage_transformation)
        val(xchr_model)
        val(no_mean_imputation)
        val(independent_se)
        val(my_plink2)
    output:
        tuple val(validation_population), path("${validation_population}.PLINK2_score_output.chr${chromosome}.sscore"), path("${validation_population}.PLINK2_score_output.chr${chromosome}.sscore.vars")
        path("${validation_population}.PLINK2_score_output.chr${chromosome}.log")
    shell:
        """
        # get column numbers of scores in score file
        sed 's/\t/,/g' ${combined_score_input} | head -n1 | sed 's/[^,]//g' | wc -c > temp1
        total_colnums=\$(cat "temp1")
        echo \$total_colnums
        seq 3 \$total_colnums | tr '\n' , | sed 's/,\$//g' > temp2
        plink_score_colnums=\$(cat "temp2")
        echo \$plink_score_colnums
        rm temp1
        rm temp2

        # get plink prefix
        echo ${plink_files} | sed 's/[][]//g' | cut -d',' -f1 | sed -e 's/.psam//g' -e 's/.pvar//g' -e 's/.pgen//g' | sed 's/,//g' | sed 's/ //g' > temp
        plink_prefix=\$(cat "temp")

        # plink score command
            # first line includes score input, and opts to read the score input header
            # first line also includes dosage transformation, no mean imputation, and independent SE parameters
            # first line also specifies the output files to include score sums, not score averages
            # finally, first line directs creation of an output file with the variants plink score was computed with
            # second line specifies the X chromosome model parameter
            # third line specifies the column numbers of the score columns
            # fourth line only keeps individuals specified in the sample list
            # fourth line include --read-freq flag is user opted to use it
            # fifth line specifies the input plink files
            # sixth line removes unequal duplicates- duplicate variant IDs were causing an error
            # seventh line specifies the output file format
        ${my_plink2} --score ${combined_score_input} header-read ${dosage_transformation} ${no_mean_imputation} ${independent_se} cols=+scoresums,-scoreavgs list-variants \
        --xchr-model ${xchr_model} \
        --score-col-nums \$plink_score_colnums \
        --keep ${sample_list_file} ${plink_read_freq_line} \
        ${plink_flag} \$plink_prefix \
        --rm-dup exclude-mismatch \
        --out ${validation_population}.PLINK2_score_output.chr${chromosome} > ${validation_population}.PLINK2_score_output.chr${chromosome}.log
        """
    stub:
        """
        touch ${validation_population}.PLINK2_score_output.chr${chromosome}.sscore
        touch ${validation_population}.PLINK2_score_output.chr${chromosome}.sscore.vars
        touch ${validation_population}.PLINK2_score_output.chr${chromosome}.log
        """
}

process concatenate_plink_score_outputs {
    publishDir "${launchDir}/plink_score_output/concatenated_outputs/", mode:'copy'
    input:
        tuple val(validation_population), path(plink_score_output_files), path(plink_score_var_lists_files)
        path(cat_outputs_script)
        val(my_python)
    output:
        tuple val(validation_population), path("${validation_population}.all_computed_PGS.txt"), path("${validation_population}.all_computed_PGS.variant_list.txt")
    shell:
        """
        ${my_python} ${cat_outputs_script} \
        --valPop '${validation_population}' \
        --score_file_list '${plink_score_output_files}' \
        --variant_list '${plink_score_var_lists_files}'
        """
    stub:
        """
        touch ${validation_population}.all_computed_PGS.txt
        touch ${validation_population}.all_computed_PGS.variant_list.txt
        """
}

process ref_pca_qc {
    publishDir "${launchDir}/pca/ref_qc/"

    input:
        tuple val(ref_key), path(ref_plink_files), path(intersect_variant_list)
        val(my_plink)
        val(hwe_threshold)
        val(maf_threshold)
        val(geno_threshold)
        val(mind_threshold)
        val(ld_prune_params)
    output:
        tuple(validation_population), path("${ref_key}.reference.pruned.{bed, bim, fam}")
        path("${ref_key}.reference.cleaned.{pgen, pvar, psam, log}")
        path(${ref_key}.reference.pruned.prune.in)
        path(${ref_key}.reference.pruned.prune.out)
        path(${ref_key}.reference.pruned.log)
    shell:
        """
        # get plink prefix
        echo ${ref_plink_files} | sed 's/[][]//g' | cut -d',' -f1 | sed -e 's/.psam//g' -e 's/.pvar//g' -e 's/.pgen//g' | sed 's/,//g' | sed 's/ //g' > temp
        plink_prefix=\$(cat "temp")

        # apply initial qc
        ${my_plink} --pfile  \$plink_prefix \
        --hwe ${hwe_threshold} \
        --maf ${maf_threshold} \
        --max-alleles 2 \
        --snps-only just-acgt \
        --rm-dup exclude-all \
        --geno ${geno_threshold} \
        --mind ${mind_threshold}
        --extract ${intersect_variant_list} \
        --allow-extra-chr \
        --autosome \
        --make-pgen \
        --out ${ref_key}.reference.cleaned

        # ld prune
        ${my_plink} --pfile  ${ref_keyn}.reference.cleaned \
        --indep-pairwise ${ld_prune_params}
        --out ${ref_key}.reference.pruned

        # extract ld pruned variants
        ${my_plink} --pfile  ${ref_key}.reference.cleaned \
        --extract ${ref_key}.reference.pruned.prune.in \
        --make-bed \
        --out ${ref_key}.reference.pruned
        """
    stub:
        """
        touch ${ref_key}.reference.pruned.bed
        touch ${ref_key}.reference.pruned.bim
        touch ${ref_key}.reference.pruned.fam
        touch ${ref_key}.reference.pruned.log
        touch ${ref_key}.reference.cleaned.pgen
        touch ${ref_key}.reference.cleaned.pvar
        touch ${ref_key}.reference.cleaned.psam
        touch ${ref_key}.reference.cleaned.log
        touch ${ref_key}.reference.prune.in
        touch ${ref_key}.reference.prune.out
        """
}

process target_pca_qc {
    publishDir "${launchDir}/pca/target_qc/"

    input:
        tuple val(validation_population), val(chromosome), path(plink_files)
        path(ref_qc_variant_list)
        val(my_plink)
    output:
        tuple(validation_population), path("${validation_population}.chr${chromosome}.reference.pruned.{pgen, pvar, psam}")
        path(${validation_population}.chr${chromosome}.reference.pruned.log)
    shell:
        """
        # get plink prefix
        echo ${plink_files} | sed 's/[][]//g' | cut -d',' -f1 | sed -e 's/.psam//g' -e 's/.pvar//g' -e 's/.pgen//g' | sed 's/,//g' | sed 's/ //g' > temp
        plink_prefix=\$(cat "temp")

        # apply initial qc
        ${my_plink} --pfile  \$plink_prefix \
        --extract ${ref_qc_variant_list} \
        --make-pgen \
        --out ${validation_population}.chr${chromosome}.reference.pruned
        """
    stub:
        """
        touch ${validation_population}.chr${chromosome}.reference.pruned.pgen
        touch ${validation_population}.chr${chromosome}.reference.pruned.psam
        touch ${validation_population}.chr${chromosome}.reference.pruned.pvar
        touch ${validation_population}.chr${chromosome}.reference.pruned.log
        """
}

process target_merge_pca_qc {
    publishDir "${launchDir}/pca/target_qc/"

    input:
        tuple val(validation_population), path(plink_files), path(merge_list)
        val(my_plink)
    output:
        tuple(validation_population), path("${validation_population}.all_chr.reference.pruned.{bed, bim, fam}")
        path(${validation_population}.chr${chromosome}.reference.pruned.log)
    shell:
        """
        # merge
        ${my_plink} --pmerge-list  ${merge_list} \
        --make-bed \
        --out ${validation_population}.all_chr.reference.pruned
        """
    stub:
        """
        touch ${validation_population}.all_chr.reference.pruned.bed
        touch ${validation_population}.all_chr.reference.pruned.bim
        touch ${validation_population}.all_chr.reference.pruned.fam
        touch ${validation_population}.all_chr.reference.pruned.log
        """
}

process ref_pca {
    publishDir "${launchDir}/pca/ref/"

    input:
        tuple val(ref_key), path(ref_plink_files)
        val(my_fraposa)

    output:
        tuple(ref_key), path(${ref_key}.reference.pruned.pcs)
    shell:
        """
        # get plink prefix
        echo ${ref_plink_files} | sed 's/[][]//g' | cut -d',' -f1 | sed -e 's/.bed//g' -e 's/.bim//g' -e 's/.fam//g' | sed 's/,//g' | sed 's/ //g' > temp
        plink_prefix=\$(cat "temp")

        # run fraposa pca
        ${my_fraposa} \$plink_prefix \
        --method shrinkage \
        --dim_ref 10
        """
    stub:
        """
        touch ${ref_key}.reference.pruned.pcs
        """
}

process project_pca {
    publishDir "${launchDir}/pca/projection/"

    input:
        tuple val(ref_key), path(ref_plink_files)
        tuple val(validation_population), path(target_plink_files)
        val(my_fraposa)

    output:
        tuple val(validation_population), val(ref_key), path(${validation_population}.${ref_key}.pca.projection.pcs), path(${validation_population}.${ref_key}.pca.projection.ref.pcs), path(${validation_population}.${ref_key}.pca.projection.dat)
    shell:
        """
        # get plink prefixes
        echo ${ref_plink_files} | sed 's/[][]//g' | cut -d',' -f1 | sed -e 's/.bed//g' -e 's/.bim//g' -e 's/.fam//g' | sed 's/,//g' | sed 's/ //g' > temp
        ref_plink_prefix=\$(cat "temp")

        echo ${target_plink_files} | sed 's/[][]//g' | cut -d',' -f1 | sed -e 's/.bed//g' -e 's/.bim//g' -e 's/.fam//g' | sed 's/,//g' | sed 's/ //g' > temp
        target_plink_prefix=\$(cat "temp")

        # run fraposa pca
        ${my_fraposa} \$ref_plink_prefix \
        --method shrinkage \
        --dim_ref 10 \
        --stu_filepref \$target_plink_prefix \
        --out ${validation_population}.${ref_key}.pca.projection
        """
    stub:
        """
        touch ${validation_population}.${ref_key}.pca.projection.pcs
        touch ${validation_population}.${ref_key}.pca.projection.ref.pcs
        touch ${validation_population}.${ref_key}.pca.projection.dat
        """
}

process ancestry_adjustment {
    publishDir "${launchDir}/pca/ancestry"

    input:
        tuple val(validation_population), path(target_plink_files), path(score_files), path(related)
        tuple val(ref_key), path(ref_plink_psam)
        tuple val(validation_population), val(ref_key), path(target_pcs), path(ref_pcs)
        val(my_pgscatalog)

    output:
        tuple(validation_population), val(ref_key), path(${validation_population}.${ref_key}.adjusted_scores.txt)
        path(${validation_population}.${ref_key}.scaled_scores.txt)
        path(${validation_population}.${ref_key}.sample_info.txt)
        path(${validation_population}.${ref_key}.norm_params.json)
        path(${validation_population}.${ref_key}.pop_assignments.txt)
        path(${validation_population}.${ref_key}.logs.txt)
    shell:
        """
        pgscatalog-ancestry-adjust -d ${validation_population} \
        -r reference \
        --psam ${ref_plink_psam} \
        --ref_pcs ${ref_pcs} \
        --target_pcs ${target_pcs} \
        -x ${related} \
        -p ${ref_key} \
        -s ${score_files} \
        -a random_forest \
        --n_popcomp 5 \
        -n zscore \
        --n_normalization 4 \
        --outdir . \
        -v
        """
    stub:
        """
        touch ${validation_population}.${ref_key}.pca.projection.pcs
        touch ${validation_population}.${ref_key}.pca.projection.ref.pcs
        touch ${validation_population}.${ref_key}.pca.projection.dat
        """
}

process pgs_percentile {
    publishDir "${launchDir}/pgs_percentile/"
    input:
        tuple val(validation_population), val(ref_key), path(adjusted_score)
        path(pgs_percentile_script)
        val(my_python)
        val(percentile)
    output:
        path("${validation_population}.${ref_key}.pgs_percentile.txt")
    shell:
        """
        ${my_python} ${pgs_percentile_input_script} \
        --scoreFile '${adjusted_score}' \
        --valPop "${validation_population}" \
        --refKey "${ref_key}" \
        --percentile "${percentile}"
        """
    stub:
        """
        ${validation_population}.${ref_key}.pgs_percentile.txt
        """
}

import groovy.json.JsonBuilder
process dump_params_to_json {
    publishDir "${launchDir}/Summary", mode: 'copy'

    input:
        val params_dict
    output:
        path('plink_score_params.json')
    shell:
        """
        echo '${new JsonBuilder(params_dict).toPrettyString().replace(';', '|')}' > plink_score_params.json
        """
}
