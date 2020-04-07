/*
 * Align reads against virus references
 */

params.fastaRef = ""
params.fastqPath = ""
params.cores = 1

fastaRef = file(params.fastaRef)
virus = fastaRef.simpleName
fastqs = Channel.fromPath("${params.fastqPath}/*.fastq.gz")
              .map { file -> tuple(file.simpleName, file) }
cores = params.cores

process createIndex {

  input:
  file(fasta) from fastaRef
  val(virus) from virus

  output:
  file("${virus}*") into genomeIdx


  //Add path to binary as it is not in path
  """
  /app/bowtie2-2.4.1-linux-x86_64/bowtie2-build $fasta $virus
  """

}

process mapReads {

  input:
  set sampName, file(fastq) from fastqs
  file(genome) from genomeIdx
  val(cores) from cores
  val(virus) from virus

  output:
  set sampName, virus, file("temp.bam") into alignment

  """
  bowtie2 -p $cores -x $virus -U $fastq -S temp.sam
  samtools view -bS temp.sam > temp.bam
  """
}

// Sort bam
process sortBam{

  input:
  set sampName, virus, file(tmp) from alignment

  output:
  file("${sampName}_${virus}.bam") into bamsort

  """
  samtools sort -o "${sampName}_${virus}.bam" $tmp
  """
}

// Index bam
process indexBams {

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
