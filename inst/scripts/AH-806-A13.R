# an example from Julia

qsub -I -l walltime=8:00:00 -l select=1:ncpus=8:mem=8GB:model=c3
module use /mod/scgc
module load anaconda R


indir=/mnt/scgc/scgc_raw/results/nextseq/2021_01_22_BrownJ_SCGC_AI-547_Test/AH-806-A13
fasta=${indir}/AH-806-A13_contigs.fasta
outdir=/mnt/scgc/scgc_nfs/lab/julia/test_new_tetpca

# if prototyping
devel_script=/mnt/scgc/scgc_nfs/lab/ben/packages/tetramers/inst/scripts/tetramer_pipeline.Rscript
Rscript --vanilla $devel_script --input $fasta --output $outdir --num_threads 8 --window 1600 --step 200

# if testing installed
installed_script=/mnt/scgc_nfs/opt/common/R/3.2.2/lib64/R/library/tetramers/scripts/tetramer_pipeline.Rscript
Rscript --vanilla $installed_script --input $fasta --output $outdir --num_threads 8 --window 1600 --step 200
