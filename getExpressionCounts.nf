/*
 * Get viral transcripts and expression counts
 */
params.genesRef = ""
params.inputBams = ""
params.outFolder = ""

container_stringtie = 'docker://bmennis/stringtie'
genesRef = file(params.genesRef)

bamsout = Channel.fromPath("${params.inputBams}")
bamsidx = Channel.fromPath("${params.inputBams}.bai")

bamsidx.into{ bamsidxST1; bamsidxST2}
bamsout.into{ bamsST1; bamsST2}

 // Derive transcripts from stringTie
 process deriveTranscripts {

   container container_stringtie

   input:
   file(bam) from bamsST1
   file(bamidx) from bamsidxST1
   file(generef) from genesRef

   output:
   file("transcripts.gff") into transcripts

   """
   stringtie $bam -o transcripts.gff -G ${generef}
   """
 }

 // Merge transcripts from different stringTie runs
 process mergeTranscripts {

   container container_stringtie

   input:
   file("transcript*.gff") from transcripts.toList()
   file(generef) from genesRef

   output:
   file("newRef.gff") into newref
   file("goodFiles")

   """
   wc -l transcript*.gff | awk -v OFS="\\t" '\$1=\$1' | grep -v '2\\s' | head -n -1 | cut -f2 > goodFiles
   stringtie --merge goodFiles -o newRef.gff -G ${generef}
   """

 }


 // Quantify transcripts from stringTie reference
 process quantifyExpression {

   publishDir "${params.outFolder}", mode: 'copy'

   container container_stringtie

   input:
   file(bam) from bamsST2
   file(bamidx) from bamsidxST2
   file(generef) from newref

   output:
   file("${bam.simpleName}_geneAbundance.tab") into counts
   file(generef)

   """
   stringtie $bam -o out.gtf -eB -G $generef -A ${bam.simpleName}_geneAbundance.tab
   """
 }
