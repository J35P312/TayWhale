if(params.help){
    println "TayWhale: a RNA-seq workflow"
    println "Usage: nextflow TayWhale.nf --r1 read1.fq --r2 --read2.fq --sample sampleID--output output_directory --ref STAR_reference_folder -c config"
    println ""
    println "Optional parameters:"
    println ""
    println "Readgroup parameters:"
    println "--rgid         readgroup id"
    println "--rglb         library"
    println "--rgpl         platform"
    println "--rgpu         unit id"

}else{

    if(!file(params.r1)) exit 1, "Error Missing read1 (--r1), type --help for a help message"

    process STAR_Aln{
        input:
           file(params.r1)

        output:


        """


        """
    }


}
