nextflow.enable.dsl=2

/*
* check params
*/
// if (params.help) {
//     exit 0, help_message()
// }
if (params.pasa_use_mysql){
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
        proteins.count().subscribe{
            num_of_proteins = it
        }
        if ( num_of_proteins > 1){
            cdhit(proteins)
            merged_protein = cdhit.out
        }else {
            merged_protein = proteins
        }

    emit:
        genome_preprocessed
        reads_trimmed
        merged_protein
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

include { augustus; convert_format_augustus; augustus_partition; create_augustus_config;copy_augustus_model; merge_result_augustus } from './modules/augustus.nf'
include { augustus_to_evm } from './modules/evidencemodeler.nf'
include {busco} from './modules/busco.nf'

workflow de_novo {

    take:
        genome_file
        genome_raw
    main:
        seqkit_sliding_trimming(genome_file)
        splited_genomes = augustus_partition(seqkit_sliding_trimming.out,splite_script)
        //get or train busco model
        if (params.augustus_config){
            augustus_config = Channel.fromPath(params.augustus_config,checkIfExists: true, type: 'dir')
        }
        else {
            augustus_config = create_augustus_config(augustus_config_tarball)
            if (!params.augustus_train_model) {
                augustus_model = busco(genome_raw,busco_db_ch,augustus_config,species,params.augustus_species)                
                model_speices = "BUSCO_${params.speices}"
            } else {
                augustus_model = Channel.fromPath(params.augustus_train_model,checkIfExists: true, type: 'dir')
            }
            copy_augustus_model(augustus_model,augustus_config)
            model_speices = params.species
        }
        augustus_out = augustus(splited_genomes.flatten(),augustus_config,model_speices)
        result = merge_result_augustus(augustus_out.collect(),genome_name)
        converted_annotation = augustus_to_evm(result,relocate_script)
    emit:
        result
        converted_annotation

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
        result = merge_result_exonerate(exonerate_out.collect(),genome_name)
        converted_annotation = convert_format_exonerate(exonerate_format_script,relocate_script,result)
    emit:
        result
        converted_annotation

}

include { hisat2index; hisat2 } from './modules/hisat2.nf'
include { trinity_de_novo_assembly; tririty_genome_guided_assembly } from './modules/trinity.nf'
include { pasa ;pasa_concat; pasa_create_tdn; pasa_seq_clean; pasa_mysql_config; pasa_sqlite_config } from './modules/pasa.nf'
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

    emit:
        assemblies_gff

}

include { evm_partition; run_evm; evm_merge_result; evm_convert_to_gff } from './modules/evidencemodeler.nf'
workflow evidence_modeler {
    take:
        genome_file_ch
        de_novo_gff_ch
        proteins_gff_ch
        transcript_gff_ch
        weights_ch
    main:
        evm_partition(genome_file_ch, de_novo_gff_ch, proteins_gff_ch, 
        transcript_gff_ch, weights_ch, params.evm_segmentsize, params.evm_overlapsize)
        partition = evm_partition.out.partition_list
        commands = evm_partition.out.commands_list.splitText(by:params.evm_batchsize,file:true)
        run_evm(commands)
        evm_merge_result(genome_file_ch, partition,run_evm.out.collect())

}
/*
WORKFLOW ENTRY POINT
*/
workflow {
    preprocess(genome_file,reads_file,proteins_ch)

    genome_preprocessed = preprocess.out.genome_preprocessed
    reads_trimmed = preprocess.out.reads_trimmed
    protein_merged = preprocess.out.merged_protein

    repeat_annotation(genome_preprocessed,repeatmasker_speices_ch)

    genome_file_softmasked = repeat_annotation.out.genome_softmasked
    genome_file_hardmasked = repeat_annotation.out.genome_hardmasked

    de_novo(genome_file_softmasked,genome_file)

    homology_pred(genome_file_softmasked,protein_merged)

    de_novo_gff = de_novo.out.converted_annotation
    protein_gff = homology_pred.out.converted_annotation

    assemblies_gff = transcriptome_pred(genome_preprocessed,genome_file_hardmasked,reads_trimmed)

    evidence_modeler(genome_preprocessed,de_novo_gff,protein_gff,assemblies_gff,weights_ch)
}
