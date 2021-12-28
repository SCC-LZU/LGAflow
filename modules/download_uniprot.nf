process download_uniprot{
    label 'basic_tools'
    label 'smalljobs'

    storeDir "${params.permanentCacheDir}/databases/uniprot/${params.busco_db}"
}