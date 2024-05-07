### This pipeline is to get the one to one orthologs between Sscrofa 11.1 and query assembly

lastdb -P8 -uNEAR db ref


last-train -P0 --revsym -E0.05 -C2 db query > train

export TMPDIR=./; parallel-fasta -j8 "lastal -E0.05 -C2 -p train db" < query | last-split -fMAF+ > 1.maf

export TMPDIR=./; maf-swap 1.maf |awk -v c=${PWD##*/} '/^s/ {$2 = (++s % 2 ? c"." : "duroc.") $2} 1' |last-split -m1 |maf-swap > 2.maf

export TMPDIR=./;last-split -r -m1e-5 1.maf |awk -v c="{PWD##*/}" '/^s/ {$2=(++s%2 ? "duroc.":c".")$2} 1' | last-postmask > 2.maf

grep -v '#' *sing.maf|awk '{$NF=NULL;print}'|awk '{print $2,$NF}'|awk '$2>0'|grep -v scor|grep -v dur|sort -k1,1 |uniq|awk '{sum+=$2}END{print sum}' >alllen

 awk '$7>$8 {print $1,$8,$7};$7<$8{print $1,$7,$8}' blasttab |awk '{print $3-$2}'|awk '{sum+=$1}END{print sum}' > one2onelen


