sourcedir=/data/home/shota/SH05/171117_NB501415_0177_AH32JTBGX5/Data/Intensities/BaseCalls/Results/04a.mapping/*.bam
destdir=/data/projects/sh-xxx-swetansu-hota-mm9-atac-jan-18/Macs2PeakCalls/

for f in $sourcedir
do
  fbase=$(basename "$f")
  echo "$f"
  echo "Inside $fbase"
  
  fsortbase=${fbase/.bam/_MACS2}
  
  echo $fsortbase

  macs2 callpeak -t $f -f BAM -n "$fsortbase" -g 1870000000.0 -p 0.01 --nomodel --shift 100 --extsize 200 -B --SPMR --call-summits --outdir $destdir

done
