# Hota et al. 2021
These scripts describe the estimation of differential chromatin accessibility using the Atac-seq data obtained under Wild-type, BRM knockout conditions, at day 4, Cardiac Progenitor and Cardiomyocte stages and under normal and high BMP4 treatment conditions.

1. The bam files (with the aligned reads) produced from the MonkeyPipeline is used to make the MACS2-based peak calls - <CallMacs2Peaks_SH05.sh>. 
2. The called peaks in the .narrowPeak file are filtered using a q-value threshold of 0.05 - <FilterMacs2PeakCalls.sh>. 
3. Consensus peak calls are then made using the bedops program - <BrmSpecificAnalyses_W_D4_GetConsensusPeaks.sh>
4. The counts of reads mapping to each of the consensus peak calls in each of the samples are then made using the featureCounts program - <GenerateCounts_BrmSpecificAnalyses_w_D4.sh>
5. The count matrix is loaded into edgeR for differential chromatin accessibility - <Differential_Atac_Peak_Analyses.R>
