nextflow.enable.dsl=2

/*
* channels
*/

genome_file = Channel.fromPath( params.genome, checkIfExists: true)
reads_file = Channel.fromPath( params.reads+'/*.fq.gz', checkIfExists: true)
proteins_ch = Channel.fromPath(params.proteins+'/*.fa',checkIfExists: true)
univec_file = Channel.fromPath( params.univec, checkIfExists: true)
species_ch = Channel.value( params.species )
busco_db_ch = Channel.value( params.busco_db )
weights_ch = Channel.fromPath( params.evm_weights, checkIfExists: true)
pasa_config_ch = Channel.fromPath( params.pasa_config, checkIfExists: true)
repeatmasker_speices_ch = Channel.value( params.rm_species )
repeatmodeler_dbname_ch = Channel.value( params.rmdb_name )
search_engine = Channel.value( params.rm_search_engine)

//get file name
fileName = file( params.genome ).simpleName
genome_name = Channel.value( fileName )


//scripts
//augustus_format_script = Channel.fromPath( workflow.projectDir + '/scripts/convert_format_augustus.pl', checkIfExists: true)
exonerate_format_script = Channel.fromPath( workflow.projectDir + '/scripts/convert_format_exonerate.pl', checkIfExists: true)
splite_script = Channel.fromPath( workflow.projectDir + '/scripts/split_genome.pl', checkIfExists: true)
//merge_evm_gff_script = Channel.fromPath( workflow.projectDir + '/scripts/merge_evm_gff.pl', checkIfExists: true)
relocate_script = Channel.fromPath( workflow.projectDir + '/scripts/relocate_augustus.awk', checkIfExists: true)

/*
* modules
*/

include {seqkit_length_filter ;seqkit_to_uppercase; seqkit_to_singleline; seqkit_sliding_trimming} from './modules/seqkit.nf'
include {busco} from './modules/busco.nf'
include {fastp} from './modules/fastp.nf'

workflow preprocess {
    take:
        genome_raw
        reads_raw
    main:
        seqkit_to_uppercase(genome_raw)
        seqkit_length_filter(seqkit_to_uppercase.out)
        seqkit_to_singleline(seqkit_length_filter.out)
        genome_preprocessed = seqkit_to_singleline.out
        fastp(genome_name ,reads_raw.toSortedList())
        reads_trimmed = fastp.out.sample_trimmed

    emit:
        genome_preprocessed
        reads_trimmed
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
            model = build_database(genome_file,repeatmodeler_dbname_ch,search_engine)
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

include { augustus; convert_format_augustus; split_for_augustus; copy_busco_model; merge_result_augustus } from './modules/augustus.nf'
include { augustus_to_evm } from './modules/evidencemodeler.nf'

workflow de_novo {

    take:
        genome_file
        species
    main:
        seqkit_sliding_trimming(genome_file)
        splited_genomes = split_for_augustus(seqkit_sliding_trimming.out,splite_script)
        //get or train busco model
        if (params.augustus_config){
            augustus_config_ch = Channel.fromPath(params.augustus_config,checkIfExists: true, type: 'dir')
            busco_model = augustus_config_ch
        }
        else {
            busco_db_dir = busco_db_ch
            busco_model = busco(genome_raw, species,busco_db_ch,species)
        }
        copy_busco_model(augustus_config_ch)
        augustus_out = augustus(splited_genomes.flatten(),species)
        result = merge_result_augustus(augustus_out.collect(),genome_name)
        converted_annotation = augustus_to_evm(result,relocate_script)
    emit:
        result
        converted_annotation

}

include { exonerate; split_for_exonerate; merge_result_exonerate; convert_format_exonerate } from './modules/exonerate.nf'
workflow homology_pred {
    take:
        genome_file_ch
        uniprot_ch
    main:
        seqkit_sliding_trimming(genome_file_ch)
        splited_genome = split_for_exonerate(seqkit_sliding_trimming.out,splite_script)
        exonerate_out = exonerate(splited_genome.flatten(),uniprot_ch)
        result = merge_result_exonerate(exonerate_out.collect(),genome_name)
        converted_annotation = convert_format_exonerate(exonerate_format_script,relocate_script,result)
    emit:
        result
        converted_annotation

}

include { hisat2index; hisat2 } from './modules/hisat2.nf'
include { trinity_de_novo_assembly; tririty_genome_guided_assembly } from './modules/trinity.nf'
include { pasa ;pasa_concat; pasa_create_tdn; pasa_seq_clean } from './modules/pasa.nf'
workflow transcriptome_pred {

    take:
        ref_genome
        masked_genome
        reads_trimmed
        pasa_config
    main:
        index = hisat2index(ref_genome)
        hisat2(genome_name,reads_trimmed,ref_genome,index)
        bam = hisat2.out.sample_bam
        assembly_gg=tririty_genome_guided_assembly(bam,600000)
        assembly_denovo=trinity_de_novo_assembly(reads_trimmed)
        assemblies=pasa_concat(assembly_gg,assembly_denovo)
        tdn = pasa_create_tdn(assembly_gg)
        pasa_seq_clean(assemblies,univec_file)
        assemblies_clean = pasa_seq_clean.out.clean
        assemblies_cln = pasa_seq_clean.out.cln
        pasa(pasa_config,masked_genome,assemblies,assemblies_clean,assemblies_cln,tdn)
        trans_gff =pasa.out.trans_gff 

    emit:
        trans_gff

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
        evm_partition(genome_file_ch, de_novo_gff_ch, proteins_gff_ch, transcript_gff_ch,weights_ch,100000,10000)
        partition = evm_partition.out.partition_list
        commands = evm_partition.out.commands_list.splitText().flatten()
        run_evm(commands)
        evm_merge_result(genome_file_ch, partition,run_evm.out.collect())

}
/*
WORKFLOW ENTRY POINT
*/
workflow {
    preprocess(genome_file,reads_file)

    genome_preprocessed = preprocess.out.genome_preprocessed
    reads_trimmed = preprocess.out.reads_trimmed
    repeat_annotation(genome_preprocessed,repeatmasker_speices_ch)

    genome_file_softmasked = repeat_annotation.out.genome_softmasked
    genome_file_hardmasked = repeat_annotation.out.genome_hardmasked

    de_novo(genome_file_softmasked, species_ch)

    homology_pred(genome_file_softmasked,proteins_ch)

    de_novo_gff = de_novo.out.converted_annotation
    protein_gff = homology_pred.out.converted_annotation

    trans_gff = transcriptome_pred(genome_preprocessed,genome_file_hardmasked,reads_trimmed,pasa_config_ch)

    evidence_modeler(genome_preprocessed,de_novo_gff,protein_gff,trans_gff,weights_ch)
}
