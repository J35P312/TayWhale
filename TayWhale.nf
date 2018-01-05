params.help=false
params.r1="none"
params.r2="none"

if(params.help){
    println "TayWhale: a RNA-seq workflow"
    println "Usage: nextflow TayWhale.nf --r1 read1.fq --r2 --read2.fq --sample sampleID --output output_directory --ref STAR_reference_folder -c config"
    println ""
    println "Optional parameters:"
    println ""
    println "Readgroup parameters:"
    println "--rglb         library"
    println "--rgpl         platform"
    println "--rgpu         unit id"

}else{

    r1=file(params.r1)
    if(!r1.exists()) exit 1, "Error Missing read1 (--r1), type --help for a help message"

    r2=file(params.r2)   
    if(!r2.exists()) exit 1, "Error Missing read2 (--r2), type --help for a help message"



    process STAR_Aln{
        publishDir "${params.working_dir}", mode: 'copy', overwrite: true

        input:

           file r1
           file r2
           
        output:
            file "${params.sample}.Chimeric.out.junction" into junctions
            file "${params.sample}.RG.Aligned.sortedByCoord.out.bam" into bam
            file "${params.sample}.RG.Aligned.sortedByCoord.out.bam" into bai

        """
        
        STAR --genomeDir ${params.ref} --readFilesIn ${r1} ${r2}  --twopassMode Basic --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --chimSegmentReadGapMax parameter 3 --alignSJstitchMismatchNmax 5 -1 5 5 --runThreadN 16 --limitBAMsortRAM 31532137230 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${params.sample}.

        picard AddOrReplaceReadGroups I= ${params.sample}.Aligned.sortedByCoord.out.bam  O= ${params.sample}.RG.Aligned.sortedByCoord.out.bam RGLB=${params.rglb} RGPL=${params.rgpl} RGPU=${params.rgpu} RGSM=${params.sample}
        rm ${params.sample}.Aligned.sortedByCoord.out.bam

        samtools index ${params.sample}.RG.Aligned.sortedByCoord.out.bam
        """
    }

    process STAR_Fusion{
        publishDir "${params.working_dir}", mode: 'copy', overwrite: true
        
        input:
            file junctions

        output:
            file "${params.sample}" into Fusion_dir

        """
            star-Fusion --genome_lib_dir ${params.ref} -J ${junctions} --output_dir ${params.sample}
        """
    }

}
