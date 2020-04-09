/*
 * Align reads against virus references
 */

params.fastaRef = ""
params.fastqPath = ""
params.cores = 1


// Select containers
container_hisat = 'docker://makaho/hisat2-zstd'
container_samtools = 'docker://hpobiolab/samtools'

fastaRef = file(params.fastaRef)
genesRef = file(params.genesRef)
virus = fastaRef.simpleName
fastqs = Channel.fromPath("${params.fastqPath}/*.fastq.gz")
              .map { file -> tuple(file.simpleName, file) }
cores = params.cores

process createIndex {

  container container_hisat
  input:
  file(fasta) from fastaRef
  val(virus) from virus

  output:
  file("${virus}*") into genomeIdx


  //Add path to binary as it is not in path
  """
  hisat2-build $fasta $virus
  """

}

process mapReads {

  container container_hisat

  publishDir "results/alignments", pattern: '*log',  mode: 'copy'

  input:
  set sampName, file(fastq) from fastqs
  file(genome) from genomeIdx
  val(cores) from cores
  val(virus) from virus

  output:
  set sampName, virus, file("temp.sam") into alignment
  file("${sampName}.log") into log

  """
  hisat2 --dta -p $cores -x $virus -U $fastq -S temp.sam 2> ${sampName}.log
  """
}

// Sort bam
process sortBam{
  container container_samtools

  input:
  set sampName, virus, file(sam) from alignment

  output:
  file("${sampName}_${virus}.bam") into bamsort

  """
  samtools view -bS $sam > temp.bam
  samtools sort -o "${sampName}_${virus}.bam" temp.bam
  """
}

// Index bam
process indexBams {
  container container_samtools

  publishDir "results/alignments", mode: 'copy'

  input:
  file(bam) from bamsort

  output:
  file("${bam}.bai") into bamsidx
  file("${bam}") into bamsout

  """
  samtools index -b $bam
  """
}
