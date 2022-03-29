nextflow.enable.dsl=2

/*
* check params
*/
if (params.help) {
    exit 0, help_message()
}

if (params.pasa_use_mysql) {
    if (params.pasa_mysql_host==''
    ||params.pasa_mysql_dbname==''
    ||params.pasa_mysql_username==''
    ||pasa_mysql_password=='') {
        exit 1,"pasa_mysql_host pasa_mysql_dbname pasa_mysql_username and pasa_mysql_password must be spcified"
    }
}

if (params.rm_species==null){
    params.rm_species= params.species
}


/*
* channels
*/

genome_file = Channel.fromPath( params.genome, checkIfExists: true)
reads_file = Channel.fromPath( params.reads+'*.{fastq,fq}.gz', checkIfExists: true)
proteins_ch = Channel.fromPath(params.proteins+'*.{fasta,fa,faa}',checkIfExists: true)

univec_file = Channel.fromPath( params.univec, checkIfExists: true)
species_ch = Channel.value( params.species )
augustus_config_tarball = Channel.fromPath(params.augustus_config_tarball, checkIfExists: true)
busco_db_ch = Channel.fromPath( params.busco_db, checkIfExists: true)
weights_ch = Channel.fromPath( params.evm_weights, checkIfExists: true)

pasa_config_template_ch = Channel.fromPath( workflow.projectDir + '/assets/pasa/pasa.CONFIG.template', checkIfExists: true)
pasa_align_config_template_ch = Channel.fromPath ( workflow.projectDir + '/assets/pasa/pasa.alignAssembly.Template.txt',checkIfExists: true)
repeatmasker_speices_ch = Channel.value( params.rm_species )
repeatmodeler_dbname_ch = Channel.value( params.rmdb_name )
search_engine = Channel.value( params.rm_search_engine)


/*
scripts
*/
exonerate_format_script = Channel.fromPath( workflow.projectDir + '/scripts/convert_format_exonerate.pl', checkIfExists: true)
splite_script = Channel.fromPath( workflow.projectDir + '/scripts/split_genome.pl', checkIfExists: true)
relocate_script = Channel.fromPath( workflow.projectDir + '/scripts/relocate_augustus.awk', checkIfExists: true)
evm_merge_script = Channel.fromPath( workflow.projectDir + '/scripts/merge_evm_gff.pl', checkIfExists: true)


/*
initialize
*/

//get file name
fileName = file( params.genome ).simpleName
genome_name = Channel.value( fileName )

//check if output_dir already exists
// result_dir = file(params.output,type: 'dir')
// dir_exists = result_dir.exists
// if (dir_exists){
//     exists = dir_exists
//     count = 0
//     while (exists) {
//         count++
//         exists = file("${params.output}_${count}",type: 'dir').exists()
//     }
//     result_dir.renameTo("${params.output}_${count}")
// }


/*
* modules
*/

include {seqkit_length_filter ;seqkit_to_uppercase; seqkit_to_singleline; seqkit_sliding_trimming} from './modules/seqkit.nf'
include {fastp} from './modules/fastp.nf'
include {cdhit} from './modules/cdhit.nf'

workflow preprocess {
    take:
        genome_raw
        reads_raw
        proteins
    main:
        seqkit_to_uppercase(genome_raw)
        seqkit_length_filter(seqkit_to_uppercase.out)
        seqkit_to_singleline(seqkit_length_filter.out)
        genome_preprocessed = seqkit_to_singleline.out

        fastp(genome_name ,reads_raw.toSortedList())
        reads_trimmed = fastp.out.sample_trimmed

        num_of_proteins = 0
        proteins.toList().count().subscribe{
            num_of_proteins = it
        }
        if ( num_of_proteins > 1 && !params.skip_cdhit){
            cdhit(proteins.collect())
            protein_merged = cdhit.out
        } else {
            protein_merged = proteins
        }

    emit:
        genome_preprocessed
        reads_trimmed
        protein_merged
}

include {repeatmasker; build_database; repeatmodeler} from './modules/repeatmasker.nf'
include {bedtools_mask} from './modules/bedtools.nf'

workflow repeat_annotation {
    take:
        genome_file
        species
    
    main:
         
        if (params.rm_lib){
            repeatmodeler_model_ch = Channel.fromPath(params.rm_lib,checkIfExists: true)
            model = repeatmodeler_model_ch
        } 
        else {
            model = build_database(genome_file,repeatmodeler_dbname_ch,search_engine).flatten().take(1)
        }
        repeatmasker(genome_file,species,search_engine)
        repeatmodeler(genome_file,model,search_engine)
        masker = repeatmasker.out.mask_gff
        modeler = repeatmodeler.out.mask_gff
        gffs = Channel.empty().concat(masker,modeler)
        bedtools_mask(genome_file,gffs.collect())
        genome_softmasked = bedtools_mask.out.softmasked_file
        genome_hardmasked = bedtools_mask.out.hardmasked_file
    
    emit:
        genome_softmasked
        genome_hardmasked

}

include { augustus; augustus_partition; create_augustus_config;copy_augustus_model; merge_result_augustus } from './modules/augustus.nf'
include { augustus_to_evm } from './modules/evidencemodeler.nf'
include  {busco; busco_get_model_name } from './modules/busco.nf'
include { train_glimmerhmm; glimmerhmm } from './modules/glimmerhmm.nf'

workflow de_novo {

    take:
        genome_file
        genome_raw
    main:
        seqkit_sliding_trimming(genome_file)
        splited_genomes = augustus_partition(seqkit_sliding_trimming.out,splite_script)

        if (params.augustus_config){
            augustus_config = Channel.fromPath(params.augustus_config,checkIfExists: true, type: 'dir')
            model_species = params.augustus_species
        }
        else {
            augustus_config = create_augustus_config(augustus_config_tarball)
            if (!params.augustus_train_model) {
                busco(genome_raw,busco_db_ch,augustus_config,params.species,params.augustus_species)
                augustus_model = busco.out.busco_model
                model_species = busco_get_model_name(augustus_model)
            } else {
                augustus_model = Channel.fromPath(params.augustus_train_model,checkIfExists: true, type: 'dir')
                model_species = file(params.augustus_train_model).getName()
            }
            copy_augustus_model(augustus_model,augustus_config)
        }

        augustus_out = augustus(splited_genomes.flatten(),augustus_config,model_species)
        merged_augustus_result = merge_result_augustus(augustus_out.collect(),genome_name)
        converted_annotation = augustus_to_evm(merged_augustus_result,relocate_script)

        // if (params.use_glimmerhmm) {
        //     model = train_glimmerhmm(genome_file,exon_file)
        //     glimmerhmm_result = glimmerhmm(genome,model)
        //     result = Channel.empty().concat(converted_annotation,glimmerhmm_result)
        // } else {
        //     result = converted_annotation
        // }
        result = converted_annotation

    emit:
        result

}

include { exonerate; exonerate_partition; merge_result_exonerate; convert_format_exonerate } from './modules/exonerate.nf'
workflow homology_pred {
    take:
        genome_file_ch
        uniprot_ch
    main:
        seqkit_sliding_trimming(genome_file_ch)
        splited_genome = exonerate_partition(seqkit_sliding_trimming.out,splite_script)
        exonerate_out = exonerate(splited_genome.flatten(),uniprot_ch)
        merged_exonerate_out = merge_result_exonerate(exonerate_out.collect(),genome_name)
        converted_annotation = convert_format_exonerate(exonerate_format_script,relocate_script,merged_exonerate_out)
        result = converted_annotation
    emit:
        result

}

include { hisat2index; hisat2 } from './modules/hisat2.nf'
include { trinity_de_novo_assembly; tririty_genome_guided_assembly } from './modules/trinity.nf'
include { pasa ;pasa_concat; pasa_create_tdn; pasa_seq_clean; pasa_mysql_config; pasa_sqlite_config;pasa_assemblies_to_orf } from './modules/pasa.nf'
workflow transcriptome_pred {

    take:
        ref_genome
        masked_genome
        reads_trimmed
    main:
        index = hisat2index(ref_genome)
        hisat2(genome_name,reads_trimmed,ref_genome,index)
        bam = hisat2.out.sample_bam
        assembly_gg=tririty_genome_guided_assembly(bam,600000)
        assembly_denovo=trinity_de_novo_assembly(reads_trimmed)
        assemblies=pasa_concat(assembly_gg,assembly_denovo)
        tdn = pasa_create_tdn(assembly_denovo)
        pasa_seq_clean(assemblies,univec_file)
        assemblies_clean = pasa_seq_clean.out.clean
        assemblies_cln = pasa_seq_clean.out.cln
        if ( params.pasa_config != '' && params.pasa_align_config != '' ){
            pasa_config = Channel.fromPath( params.pasa_config, checkIfExists: true)
            pasa_align_config = Channel.fromPath (params.pasa_align_config,checkIfExists: true)
        }
        else if ( params.pasa_use_mysql ) {
            pasa_mysql_config(pasa_config_template_ch,
            pasa_align_config_template_ch,
            params.pasa_mysql_host,
            params.pasa_mysql_dbname,
            params.pasa_mysql_username,
            params.pasa_mysql_password)

            pasa_config = pasa_mysql_config.out.pasa_conf
            pasa_align_config = pasa_mysql_config.out.pasa_alignassembly_conf
        } else {
            pasa_sqlite_config(pasa_config_template_ch,
            pasa_align_config_template_ch,
            params.pasa_sqlite_path)
            pasa_config = pasa_sqlite_config.out.pasa_conf
            pasa_align_config = pasa_sqlite_config.out.pasa_alignassembly_conf
        }

        pasa(pasa_align_config,pasa_config,masked_genome,assemblies,assemblies_clean,assemblies_cln,tdn)

        assemblies_gff = pasa.out.assemblies_gff 
        assemblies_fasta = pasa.out.assemblies_fasta
        transdecoder_gff = pasa_assemblies_to_orf(assemblies_fasta,assemblies_gff)

    emit:
        assemblies_gff
        transdecoder_gff

}

include { evm_partition; run_evm; evm_convert_and_merge_result; evm_merge_input } from './modules/evidencemodeler.nf'
workflow evidence_modeler_annotation {
    take:
        genome_file_ch
        de_novo_gff_ch
        proteins_gff_ch
        transcript_gff_ch
        transdecoder_gff
        weights_ch
    main:
        gene_pred_gff = evm_merge_input(de_novo_gff_ch,transdecoder_gff)
        evm_partition(genome_file_ch, gene_pred_gff, proteins_gff_ch, 
        transcript_gff_ch, weights_ch, params.evm_segmentsize, params.evm_overlapsize)
        partition = evm_partition.out.partition_list
        commands = evm_partition.out.commands_list.splitText(by:params.evm_batchsize,file:true)
        run_evm(commands)
        evm_convert_and_merge_result(genome_file_ch, partition,run_evm.out.collect(),evm_merge_script)

}

/*
WORKFLOW ENTRY POINT
*/
workflow {
    preprocess(genome_file,reads_file,proteins_ch)
    genome_preprocessed = preprocess.out.genome_preprocessed
    reads_trimmed = preprocess.out.reads_trimmed
    protein_merged = preprocess.out.protein_merged

    repeat_annotation(genome_preprocessed,repeatmasker_speices_ch)
    genome_file_softmasked = repeat_annotation.out.genome_softmasked
    genome_file_hardmasked = repeat_annotation.out.genome_hardmasked

    de_novo(genome_file_softmasked,genome_file)
    de_novo_gff = de_novo.out.result
    
    homology_pred(genome_file_softmasked,protein_merged)
    protein_gff = homology_pred.out.result

    transcriptome_pred(genome_preprocessed,genome_file_hardmasked,reads_trimmed)
    assemblies_gff = transcriptome_pred.out.assemblies_gff
    transdecoder_gff = transcriptome_pred.out.transdecoder_gff

    evidence_modeler_annotation(genome_preprocessed,de_novo_gff,protein_gff,assemblies_gff,transdecoder_gff,weights_ch)
}

def help_message() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________

    ${c_yellow}Usage examples:${c_reset}
    nextflow run dakehero/lgaflow -profile test,local,conda
    nextflow run dakehero/lgaflow --cores 8 --genome c_eleagans.fa --reads input.csv --annotation gtf_virus.csv --autodownload hsa --pathway hsa

    ${c_yellow}Input:${c_reset}
    ${c_green}--genome${c_reset}                 CSV file with genome reference FASTA files (one path in each line)
    ${c_green}--reads${c_reset}                  A CSV file following the pattern: Sample,R,Condition,Source for single-end or Sample,R1,R2,Condition,Source for paired-end
                                        ${c_dim}(check terminal output if correctly assigned)
                                        Per default, all possible comparisons of conditions in one direction are made. Use --deg to change.${c_reset}
    ${c_green}--proteins${c_reset}
    ${c_green}--species${c_reset}       Specifies the species identifier for downstream path analysis. (DEPRECATED)
                                        If `--include_species` is set, reference genome and annotation are added and automatically downloaded. [default: $params.species]
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly | Homo_sapiens.GRCh38.98]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel | Escherichia_coli_k_12.ASM80076v1.45]
                                        - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly | Mus_musculus.GRCm38.99.gtf]
                                        - mau [Ensembl: Mesocricetus_auratus.MesAur1.0.dna.toplevel | Mesocricetus_auratus.MesAur1.0.100]${c_reset}


    ${c_yellow}Preprocessing options:${c_reset}
    --mode                             Either 'single' (single-end) or 'paired' (paired-end) sequencing [default: $params.mode]
    --fastp_additional_params          additional parameters for fastp [default: $params.fastp_additional_params]
    --skip_sortmerna                   Skip rRNA removal via SortMeRNA [default: $params.skip_sortmerna] 
    --hisat2_additional_params        additional parameters for HISAT2 [default: $params.hisat2_additional_params]
    --featurecounts_additional_params  additional parameters for FeatureCounts [default: $params.featurecounts_additional_params]

    ${c_yellow}DEG analysis options:${c_reset}
    --strand                 0 (unstranded), 1 (stranded) and 2 (reversely stranded) [default: $params.strand]
    --tpm                    Threshold for TPM (transcripts per million) filter. A feature is discared, if for all conditions the mean TPM value of all 
                             corresponding samples in this condition is below the threshold. [default: $params.tpm]
    --deg                    A CSV file following the pattern: conditionX,conditionY
                             Each line stands for one differential gene expression comparison.  
                             Must match the 'Condition' labels defined in the CSV file provided via --reads.  
    --pathway                Perform different downstream pathway analysis for the species. [default: $params.pathway]
                             ${c_dim}Currently supported are:
                                 - hsa | Homo sapiens
                                 - mmu | Mus musculus
                                 - mau | Mesocricetus auratus${c_reset}
    --feature_id_type        ID type for downstream analysis [default: $params.feature_id_type]

    ${c_yellow}Transcriptome assembly options:${c_reset}
    --assembly               Perform de novo and reference-based transcriptome assembly instead of DEG analysis [default: $params.assembly]
    --busco_db               The database used with BUSCO [default: $params.busco_db]
                             ${c_dim}Full list of available data sets at https://busco.ezlab.org/v2/frame_wget.html ${c_reset}
    --dammit_uniref90        Add UniRef90 to the dammit databases (time consuming!) [default: $params.dammit_uniref90]

    ${c_yellow}Computing options:${c_reset}
    --cores                  Max cores per process for local use [default: $params.cores]
    --max_cores              Max cores used on the machine for local use [default: $params.max_cores]
    --memory                 Max memory in GB for local use [default: $params.memory]
    --output                 Name of the result folder [default: $params.output]

    ${c_yellow}Caching:${c_reset}
    --permanentCacheDir      Location for auto-download data like databases [default: $params.permanentCacheDir]
    --condaCacheDir          Location for storing the conda environments [default: $params.condaCacheDir]
    --singularityCacheDir    Location for storing the singularity images [default: $params.singularityCacheDir]
    ${c_dim}--workdir                Working directory for all intermediate results [default: $params.workdir] (DEPRECATED: use `-w your/workdir` instead)${c_reset}
    --softlink_results       Softlink result files instead of copying.

    ${c_yellow}Execution/Engine profiles:${c_reset}
    The pipeline supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.: -profile ${c_green}local${c_reset},${c_blue}conda${c_reset}
    
    ${c_green}Executer${c_reset} (choose one):
      local
      slurm
      pbs
    
    ${c_blue}Engines${c_reset} (choose one):
      conda
      docker
      singularity
    
    Per default: -profile local,conda is executed. 

    ${c_dim}For a test run (~ 15 min), add "test" to the profile, e.g. -profile test,local,conda.
    The command will create all conda environments and download and run test data.

    We also provide some pre-configured profiles for certain HPC environments:    
      ara (slurm, conda and parameter customization)
    ${c_reset}
    """.stripIndent()
}
