params.help=false
params.r1="none"
params.r2="none"
params.libtype="A"

if(params.help){
    println "TayWhale: a RNA-seq workflow"
    println "Usage: nextflow TayWhale.nf --r1 read1.fq.gz --r2 --read2.fq.gz --sample sampleID --output output_directory -c TayWhale.conf"
    println ""
    println "Optional parameters:"
    println ""
    println "Readgroup parameters:"
    println "--rglb         library"
    println "--rgpl         platform"
    println "--rgpu         unit id"
    println ""
    println "library-info"
    println "--libtype IU - reads poining Inward Unstranded library"
    println "--libtype ISF - reads pointing Inward, Stranded library, first read on forward strand (default)"
    println ""
    println "I = inward"
    println "O = outward"
    println "M = matching"
    println ""
    println "S = stranded"
    println "U = unstranded"
    println ""
    println "F = read 1 comes from the forward strand"
    println "R = read 1 comes from the reverse strand"
    

}else{

    r1=file(params.r1)
    if(!r1.exists()) exit 1, "Error Missing read1 (--r1), type --help for a help message"

    r2=file(params.r2)   
    if(!r2.exists()) exit 1, "Error Missing read2 (--r2), type --help for a help message"

    r1_salmon=file(params.r1)
    r2_salmon=file(params.r2)

    process salmon{
        publishDir "${params.output}", mode: 'copy', overwrite: true
        cpus 16

        input:

           file r1_salmon
           file r2_salmon
           
        output:
            file "${params.sample}.sf" into salmon_quant

        """
        singularity exec ${params.taywhale} salmon quant --index ${params.salmon_index} --libType ${params.libtype} --output salmon_output -1 <( zcat ${r1_salmon} ) -2 <( zcat ${r2_salmon} )
        mv salmon_output/quant.sf ${params.sample}.sf
        """
    } 

    process STAR_Aln{
        publishDir "${params.output}", mode: 'copy', overwrite: true
        cpus 16

        input:

           file r1
           file r2
           
        output:
            file "${params.sample}.Chimeric.out.junction" into junctions
            file "${params.sample}.RG.Aligned.sortedByCoord.out.bam" into star_bam
            file "${params.sample}.RG.Aligned.sortedByCoord.out.bam.bai" into bai
            file "${params.sample}.ReadsPerGene.out.tab" into geneCounts

        """
        
        singularity exec ${params.taywhale} STAR --genomeDir ${params.STAR_ref_dir} --readFilesIn ${r1} ${r2}  --twopassMode Basic --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --chimSegmentReadGapMax parameter 3 --alignSJstitchMismatchNmax 5 -1 5 5 --runThreadN 16 --limitBAMsortRAM 31532137230 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${params.sample}. --quantMode GeneCounts --outSAMstrandField intronMotif --readFilesCommand gunzip -c
      
        singularity exec ${params.taywhale} picard AddOrReplaceReadGroups I= ${params.sample}.Aligned.sortedByCoord.out.bam  O= ${params.sample}.RG.Aligned.sortedByCoord.out.bam RGLB=${params.rglb} RGPL=${params.rgpl} RGPU=${params.rgpu} RGSM=${params.sample}
        rm ${params.sample}.Aligned.sortedByCoord.out.bam

        singularity exec ${params.taywhale} samtools index ${params.sample}.RG.Aligned.sortedByCoord.out.bam
        """
    }


    
    process STAR_Fusion{
        publishDir "${params.output}", mode: 'copy', overwrite: true
        
        input:
            file junctions

        output:
            file "${params.sample}" into Fusion_dir

        """
            singularity exec ${params.taywhale} STAR-Fusion --genome_lib_dir ${params.ctat_folder} -J ${junctions} --output_dir ${params.sample}
        """
    }

    bam=Channel.create()
    stringtie_bam=Channel.create()

    Channel
            .from star_bam
            .separate( bam, stringtie_bam){ a -> [a,a,a] }

    process GATK_Split{
        publishDir "${params.output}", mode: 'copy', overwrite: true
        
        input:
            file bam
            file bai

        output:
            file "${params.sample}.RG.split.Aligned.sortedByCoord.out.bam" into GATK_bam
            file "${params.sample}.RG.split.Aligned.sortedByCoord.out.bam.bai" into GATK_bai

        """
        java -jar ${params.GATK} -T SplitNCigarReads -R ${params.ref} -I ${bam}  -o ${params.sample}.RG.split.Aligned.sortedByCoord.out.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
        singularity exec ${params.taywhale} samtools index ${params.sample}.RG.split.Aligned.sortedByCoord.out.bam
        """

    }


    process StringTie{
         publishDir "${params.output}", mode: 'copy', overwrite: true

         cpus 16

         input:
            file stringtie_bam

        output:

           file "${stringtie_bam}.stringtie.gff" into stringtie_gff
           file "${params.sample}.stats" into gffstats
           file "${params.sample}.annotated.gtf" into annotated_gff

        """
        singularity exec ${params.taywhale} stringtie ${stringtie_bam} -G ${params.gff} > ${stringtie_bam}.stringtie.gff
        ${params.gffcompare} -r ${params.gff} -o ${params.sample} ${stringtie_bam}.stringtie.gff
        """

    }

    process GATK_ASE{
        publishDir "${params.output}", mode: 'copy', overwrite: true

	cpus 16

        input:
            file GATK_bam
            file GATK_bai

        output:
            file "${params.sample}.vcf" into GATK_haplotype_vcf
            file "${params.sample}.annotated.GATKASE.vcf" into GATK_haplotype_ASE_vcf
            file "${params.sample}.GATKASE.csv" into GATK_ASE_CSV

        """
        java -jar ${params.GATK} -R ${params.ref} -T HaplotypeCaller -I ${GATK_bam} -stand_call_conf 10 -o ${params.sample}.vcf -dontUseSoftClippedBases --min_mapping_quality_score 10
        java -jar ${params.GATK} -R ${params.ref} -T ASEReadCounter -o ${params.sample}.GATKASE.csv -I ${GATK_bam} -sites ${params.sample}.vcf
        python ${params.BootstrapAnn} --vcf ${params.sample}.vcf --ase ${params.sample}.GATKASE.csv >${params.sample}.GATKASE.vcf
        ${params.VEP_exec_file} -i ${params.sample}.GATKASE.vcf  -o ${params.sample}.annotated.GATKASE.vcf --fork 6 ${params.vep_args}
        """
    }


}

