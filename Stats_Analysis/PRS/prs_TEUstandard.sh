#!/bin/bash

# Specify a job name
#$ -N prs_TEUstandard.sh

# --- Parameters for the Queue Master ---
# Project name and target queue
#$ -P hunter.prjc
#$ -q short.qc
#$ -pe shmem 2

# Specify the working directory
#$ -wd /well/hunter/users/ryo297/ukb/prs/projects/htn-evangelou2018

# Log file locations - we want all logs to go to /well/hunter/projects/ukb
#$ -o /well/hunter/projects/ukb/prs/projects/htn-evangelou2018/logs/
#$ -e /well/hunter/projects/ukb/prs/projects/htn-evangelou2018/logs/

# Print some useful data about the job to help with debugging
echo "------------------------------------------------"
echo "Job ID: $JOB_ID"
echo "SGE Job ID: $SGE_JOB_ID"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

# Paths to important files
UKBDATAPATH=/well/ukbb-wtchg/v3/imputation
TEUSAMPLE=/well/hunter/shared/ukb33952/sample/ukb33952_imp_chr19_v3_s487296.sample
BETAS=inputs/Betas.csv

# Paths to software tools
BGENIXPATH=/apps/well/bgenix/1.1.1/bin/bgenix
CATBGENPATH=/apps/well/bgenix/1.1.1/bin/cat-bgen

# Paths to QC files
IMPUTEDB=/well/hunter/shared/ukb_mfi/imputeinfo.db
IMPUTETSV=/well/hunter/shared/ukb_mfi/ukb_mfi_v3.tsv
USEDINPCA=/well/hunter/shared/ukb33952/sample/usedinpca.txt
WHITEBRIT=/well/hunter/shared/ukb33952/sample/whitebrit.txt

# Load software modules
module load PLINK/2.00a2.3_x86_64
module load SQLite/3.29.0-GCCcore-8.3.0

# Name of the resulting PRS
SBP_PRS=sbp_prs_20211111
DBP_PRS=dbp_prs_20211111

###################################

awk -F, '{ if (NR>1) { print $1 }}' ${BETAS} > outputs/rsidlist.txt

awk -F, '{ if (NR>1) { print sprintf("%02d", $2)":"$3"-"$3 }}' ${BETAS} > outputs/chrposlist.txt

for j in {0..5}
do
  START=$((4*j + 1))
  END=$(( 4*(j+1) < 22 ? 4*(j+1) : 22))
  for ((i=START;i<=END;i++))
  do
    ${BGENIXPATH} -g ${UKBDATAPATH}/ukb_imp_chr${i}_v3.bgen \
-incl-rsids outputs/rsidlist.txt \
-incl-range outputs/chrposlist.txt > outputs/chr_${i}.bgen &
  done
  wait
done

# Combine the .bgen files for each chromosome into one
cmd=""
for i in {1..22}
do
  cmd=$cmd"outputs/chr_${i}.bgen "
done
${CATBGENPATH} -g  $cmd -og outputs/initial_chr.bgen -clobber
${BGENIXPATH} -g outputs/initial_chr.bgen -index -clobber

# Remove the individual chromosome files
for i in {1..22}
do
  rm outputs/chr_${i}.bgen
done


########################################################

# Import the betas into the sqlite database as a table called Betas
sqlite3 outputs/initial_chr.bgen.bgi "DROP TABLE IF EXISTS Betas;"
sqlite3 -separator "," outputs/initial_chr.bgen.bgi ".import ${BETAS} Betas"

sqlite3 outputs/initial_chr.bgen.bgi "DROP TABLE IF EXISTS Joined;"
# And inner join it to the index table (Variants), making a new table (Joined)
# By joining on alleles as well as chromosome and position 
# we can ensure only the relevant alleles from any multi-allelic SNPs are retained
sqlite3 -header -csv outputs/initial_chr.bgen.bgi "CREATE TABLE Joined AS SELECT Variant.* FROM Variant INNER JOIN Betas ON Variant.chromosome = printf('%02d', Betas.chr_name) AND Variant.position = Betas.chr_position AND Variant.allele1 = Betas.effect_allele AND Variant.allele2 = Betas.reference_allele UNION SELECT Variant.* FROM Variant INNER JOIN Betas ON Variant.chromosome = printf('%02d', Betas.chr_name) AND Variant.position = Betas.chr_position AND Variant.allele1 = Betas.reference_allele AND Variant.allele2 = Betas.effect_allele;"

# Filter the .bgen file to include only the alleles specified in the Betas for each SNP 
${BGENIXPATH} -g outputs/initial_chr.bgen -table Joined  > outputs/single_allelic.bgen
${BGENIXPATH} -g outputs/single_allelic.bgen -index -clobber
${BGENIXPATH} -g outputs/single_allelic.bgen -list > outputs/single_allelic.txt


########################################################

sqlite3 "" <<EndOfSqlite3Commands
ATTACH 'outputs/single_allelic.bgen.bgi' AS db;
ATTACH '${IMPUTEDB}' AS Impute;
DROP TABLE IF EXISTS db.Info;
CREATE TABLE db.Info AS 
SELECT Variant.* FROM db.Variant INNER JOIN Impute.Info ON (
Variant.chromosome = Impute.Info.CHROMOSOME AND 
Variant.position = Impute.Info.POSITION AND 
Variant.allele1 = Impute.Info.ALLELE1 AND 
Variant.allele2 = Impute.Info.allele2
) WHERE Impute.Info.INFO >= 0.4;
EndOfSqlite3Commands

${BGENIXPATH} -g outputs/single_allelic.bgen -table Info  > outputs/impute_info.bgen
${BGENIXPATH} -g outputs/impute_info.bgen -index -clobber
${BGENIXPATH} -g outputs/impute_info.bgen -list > outputs/impute_info.txt

########################################################

plink2 --bgen outputs/impute_info.bgen ref-first \
--hard-call-threshold 0.1 \
--sample ${TEUSAMPLE} \
--memory 15000 \
--set-all-var-ids @:# \
--freq \
--make-pgen \
--out outputs/raw 


########################################################

awk '/^[^#]/ { if( $5>0.49 && $5<0.51 && ( ($3=="A" && $4=="T") || ($4=="T" && $3=="A") || ($3=="C" && $4=="G") || ($4=="G" && $3=="C") ) ) { print $0 }}' outputs/raw.afreq > outputs/exclrsIDs_ambiguous.txt

plink2 --pfile outputs/raw \
--memory 15000 \
--exclude outputs/exclrsIDs_ambiguous.txt \
--maf 0.005 \
--write-snplist \
--make-pgen \
--out outputs/snpQC 


plink2 --pfile outputs/snpQC \
--memory 15000 \
--indep-pairwise 500 'kb' 0.3 \
--out outputs/AfterQC

########################################################

# Extract SBP weights
awk -F, 'NR > 1 { print $2":"$3" "$4" "$6 }' ${BETAS} > inputs/score.txt
# Remove NA row for SBP weights
awk '{ if( $3!="NA") {print }}' inputs/score.txt > inputs/SBP_score.txt

# Extract DBP weights
awk -F, 'NR > 1 { print $2":"$3" "$4" "$9 }' ${BETAS} > inputs/DBP_score.txt

plink2 --pfile outputs/raw \
--memory 15000 \
--extract outputs/snpQC.snplist \
--exclude outputs/AfterQC.prune.out \
--keep-fam ${USEDINPCA} \
--score inputs/SBP_score.txt no-mean-imputation cols=fid,denom,dosagesum,scoreavgs,scoresums \
--out outputs/${SBP_PRS}


plink2 --pfile outputs/raw \
--memory 15000 \
--extract outputs/snpQC.snplist \
--exclude outputs/AfterQC.prune.out \
--keep-fam ${USEDINPCA} \
--score inputs/DBP_score.txt no-mean-imputation cols=fid,denom,dosagesum,scoreavgs,scoresums \
--out outputs/${DBP_PRS}

