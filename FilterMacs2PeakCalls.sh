sourcedir=/data/projects/sh-xxx-swetansu-hota-mm9-atac-jan-18/Jan19_Macs2PeakCalls/*_MACS2_peaks.narrowPeak
destdir=/data/projects/sh-xxx-swetansu-hota-mm9-atac-jan-18/FilteredMacs2PeakCalls/

for f in $sourcedir
do
  fbase=$(basename "$f")
  echo "$f"
  echo "Inside $fbase"
  
  fsortbase=${fbase/_MACS2_peaks.narrowPeak/_Filtered_MACS_peaks.bed}
  
  echo $fsortbase

  awk '{ if($9 >= 1.30103) { print $1"\t"$2"\t"$3"\t"$4}}' $f | sortBed -i > $destdir/"$fsortbase"

done
