/*
 * Remove reads mapping to virus genome
 */

params.alignmentPath = ""
params.virus = "COV_SARS2"

bams = Channel.fromFilePairs("${params.alignmentPath}/*{hg38,${params.virus}}.bam", flat: true)

process makeSharedList {

  input:
  set sampID, file(human), file(virus) from bams

  output:
  file("shared.list") into sharedReads
  set sampID, file(human) into humanBams
  set sampID, file(virus) into virusBams

  """
  samtools view -F4 $human | awk '{print $1}' | sort | uniq > human.list
  samtools view -F4 $virus | awk '{print $1}' | sort | uniq > virus.list
  cat human.list virus.list | sort | uniq -c | sort -nr | awk '{if($1==2) {print $2}}' > shared.list
  ​"""
}

process filterHuman {

  input:
  file(sharedReads) from sharedReads
  set sampID, file(human) from humanBams

  output:
  file("${sampID}_human.uniq.bam") into humanFinal

  """
  picard FilterSamReads I=$human O="${sampID}_human.uniq.bam" READ_LIST_FILE=$sharedReads FILTER=excludeReadList SORT_ORDER=coordinate
  """
}

process filterVirus {

  input:
  file(sharedReads) from sharedReads
  set sampID, file(virus) from virusBams

  output:
  file("${sampID}_virus.uniq.bam") into virusFinal

  """
  picard FilterSamReads I=$virus O="${sampID}_virus.uniq.bam" READ_LIST_FILE=$sharedReads FILTER=excludeReadList SORT_ORDER=coordinate
  """
}
