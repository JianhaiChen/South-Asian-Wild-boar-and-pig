cat *1.clean.fq.gz  >all.wild1.fq.gz
cat *2.clean.fq.gz >all.wild2.fq.gz


fastp -w 10 -i *w*1.*q.gz -I *w*2.*q.gz -o cl_1.cl.fq.gz -O cl_2.cl.fq.gz

rdir=/nas_data/04.pig/siberia_wildboar/D1;ref=Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa;t=10;export TMPDIR=./;bwa mem -B 1 -t $t -Ma -Y -R "@RG\tID:${PWD##*/}\tPL:ILLUMINA\tPU:illumina\tLB:${PWD##*/}\tSM:${PWD##*/}" $rdir/$ref  cl_1.cl.fq.gz cl_2.cl.fq.gz |samtools view -@ $t -Sb -o bwa.bam; samtools fixmate -m -@ 10 bwa.bam bwa.fm.bam; samtools sort -@ 10 -m 10G -o sort.bam bwa.fm.bam; samtools markdup -s -@ 10 sort.bam mkdup.bam;samtools index mkdup.bam
