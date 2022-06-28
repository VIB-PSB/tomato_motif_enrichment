nextflow.enable.dsl=2

log.info """\
          NEXTFLOW PIPELINE FOR MOTIF ENRICHMENT IN TOMATO
         ==================================================
         Feature file : ${params.Feature_file}
         Directory containing the set files : ${params.Set_directory}
         Set file : ${params.Set_file}
         Output directory : ${params.OutDir}
         Are set files in a directory? : ${params.Sets_in_dir}
         Print help : ${params.help}
         """
         .stripIndent()


def helpMessage() {
    log.info """\
            NEXTFLOW PIPELINE FOR MOTIF ENRICHMENT IN TOMATO
            ==================================================
            Feature file : ${params.Feature_file}
            Directory containing the set files : ${params.Set_directory}
            Set file : ${params.Set_file}
            Output directory : ${params.OutDir}
            are set files in a directory? : ${params.Sets_in_dir}


            Mandatory arguments:
            --Feature_file                 Feature file
            --Set_file                     Set file
            --OutDir                       Output directory

            Optional arguments:
            --Sets_in_dir                  Are set files in a directory?
            --Set_directory                Directory containing the set files
            --help                         Print help

            """
            .stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

process runEnricher{
    publishDir "$params.OutDir", pattern: "*_output_info.txt", mode: 'copy'
    module 'python/x86_64/3.6.5'

    input:
    tuple path(set_file), path(feature_file)
    path enricher_script
    path outfile_script
    path metadata_file

    output:
    path "${set_file.baseName}_enrichmentOutput.txt"
    path "${set_file.baseName}_output_info.txt"

    script:

    """
    ./$enricher_script $feature_file $set_file -f 0.05 -p -o ${set_file.baseName}_enrichmentOutput.txt
    python3 $outfile_script ${set_file.baseName}_enrichmentOutput.txt ${set_file.baseName}_output_info.txt $metadata_file

    """
}

workflow {

    if (params.Sets_in_dir == true) {

        set_files = Channel.fromPath("${params.Set_directory}/*.txt")
                                .combine(Channel.fromPath(params.Feature_file))


    }

    else {

        set_files = Channel.fromPath(params.Set_file)
                                .combine(Channel.fromPath(params.Feature_file))        
    }

   script_enricher = "$baseDir/scripts/enricherv2.4"
   script_add_info = "$baseDir/scripts/addMetadata_mot_enr.py"

   runEnricher(set_files, script_enricher, script_add_info, "$baseDir/data/tomato_motif_info.txt")
}
