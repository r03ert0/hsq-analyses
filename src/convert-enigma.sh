hostname

na=("SHIP" "TREND" "ADNI1" "ADNI2" "IMAGEN1" "IMAGEN2" "IMAGEN3" "OSAKA" "LOTHIAN")

d[0]="2015_09SHIP"
d[1]="2015_09SHIP"
d[2]="2015_10ADNI"
d[3]="2015_10ADNI"
d[4]="2016_01Imagen"
d[5]="2016_01Imagen"
d[6]="2016_01Imagen"
d[7]="2016_01Osaka"
d[8]="2016_02Lothian"

g[0]="_Data_Germany-SHIP/ENIGMA2-GCTA-SHIP-SHIP_TREND/ENIGMA2-GCTA-SHIP"
g[1]="_Data_Germany-SHIP/ENIGMA2-GCTA-SHIP-SHIP_TREND/ENIGMA2-GCTA-SHIP-TREND"
g[2]="ADNI_1_GWAS_Plink/ADNI_cluster_01_forward_757LONI"
g[3]="ADNI_GO_2_OmniExpress/ADNI_GO_2_Forward_Bin"
g[4]="all_subjects_all_snps_wave1/all_subjects_all_snps_wave1"
g[5]="all_subjects_all_snps_wave2/all_subjects_all_snps_wave2"
g[6]="all_subjects_all_snps_wave3/all_subjects_all_snps_wave3"
g[7]="Handai_GCTA/Handai_GCTA"
g[8]="ENIGMA_GCTA_LBC1936/LBC36_clean_231009"

p[0]="_Data_Germany-SHIP/ENIGMA2-GCTA-SHIP-SHIP_TREND/ENIGMA2-GCTA-SHIP.pheno"
p[1]="_Data_Germany-SHIP/ENIGMA2-GCTA-SHIP-SHIP_TREND/ENIGMA2-GCTA-SHIP-TREND.pheno"
p[2]="ADNI_1_GWAS_Plink/ADNI_cluster_01_forward_757LONI.pheno"
p[3]="ADNI_GO_2_OmniExpress/ADNI_GO_2_Forward_Bin.pheno"
p[4]="all_subjects_all_snps_wave1/all_subjects_all_snps_wave1.pheno"
p[5]="all_subjects_all_snps_wave2/all_subjects_all_snps_wave2.pheno"
p[6]="all_subjects_all_snps_wave3/all_subjects_all_snps_wave3.pheno"
p[7]="Handai_GCTA/Handai_GCTA.pheno"
p[8]="ENIGMA_GCTA_LBC1936/LBC36_LandRvolumes.csv"

raw='../data/raw/'
o="../data/derived/"
plink='/local/gensoft2/exe/plink/1.90b2m/scripts/plink'
ldir='../bin/lift'
log="../data/derived/log"

mkdir -p $o

if [ true ]; then
# Copy bed/bim files
    if [ true ]; then
    # 1. Convert ped/fam/map to bed/bim
    src=("2016_01Osaka/Handai_GCTA/Handai_GCTA")
    dst=("OSAKA")
    for ((i=0;i<1;i++)); do
        mkdir -p $o/${dst[$i]}/01.genotype/
        $plink --file $raw/${src[$i]} --make-bed --out $o/${dst[$i]}/01.genotype/all --noweb
    done
    fi
    
    # 2. Copy the remaining files
    if [ true ]; then
    src=("2015_09SHIP/_Data_Germany-SHIP/ENIGMA2-GCTA-SHIP-SHIP_TREND/ENIGMA2-GCTA-SHIP" \
         "2015_09SHIP/_Data_Germany-SHIP/ENIGMA2-GCTA-SHIP-SHIP_TREND/ENIGMA2-GCTA-SHIP-TREND" \
         "2015_10ADNI/ADNI_1_GWAS_Plink/ADNI_cluster_01_forward_757LONI" \
         "2015_10ADNI/ADNI_GO_2_OmniExpress/ADNI_GO_2_Forward_Bin" \
         "2016_01Imagen/all_subjects_all_snps_wave1/all_subjects_all_snps_wave1" \
         "2016_01Imagen/all_subjects_all_snps_wave2/all_subjects_all_snps_wave2" \
         "2016_01Imagen/all_subjects_all_snps_wave3/all_subjects_all_snps_wave3" \
         "2016_02Lothian/ENIGMA_GCTA_LBC1936/LBC36_clean_231009")
    dst=("SHIP" "TREND" "ADNI1" "ADNI2" "IMAGEN1" "IMAGEN2" "IMAGEN3" "LOTHIAN");
    for ((i=0;i<8;i++)); do
        mkdir -p $o/${dst[$i]}/01.genotype
        cp $raw/${src[$i]}.bed $o/${dst[$i]}/01.genotype/all.bed
        cp $raw/${src[$i]}.bim $o/${dst[$i]}/01.genotype/all.bim
        cp $raw/${src[$i]}.fam $o/${dst[$i]}/01.genotype/all.fam
    done
    fi
fi

if [ true ]; then
# Create subjects.txt lists
    for i in ${na[@]}; do
        cat $o/$i/01.genotype/all.fam|cut -d' ' -f 2 > $o/$i/subjects.txt
    done
fi

if [  ]; then
# Display summary information for genotypes
    mkdir -p $log
    for ((i=0;i<9;i++)); do
        echo ">> "${na[$i]}
        $plink \
            --bfile $o/${na[$i]}/01.genotype/all \
            --missing \
            --hardy \
            --freq \
            --check-sex \
            --noweb --allow-no-sex --out $log/plink_${na[$i]}
        echo
    done
fi

if [  ]; then
# Display SNP position in bim, hg18 and hg19
find $o -name *.bim|while read bim; do
    echo $bim;
    cat $bim|awk "/1[ \t]rs/{print $1}"|head -n 10|while read row;do
        rs=$(echo $row|awk '{print $2}');
        a=$(echo $row|awk '{printf "%s\t%s\t%s",$2,$1,$4}');
        b=$(mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -B -N -D hg18 -e 'SELECT name,chrom,chromEnd FROM snp130 WHERE name = "'$rs'"');
        c=$(mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -B -N -D hg19 -e 'SELECT name,chromEnd,alleles,alleleFreqs FROM snp138Common WHERE name = "'$rs'"');
        echo "bim:" $a "hg18:" $b "hg19:" $c;
    done;
    echo;
done
fi

if [ true ]; then
# Move all hg18 data to hg19
    if [ true ]; then
    # Modify ped files to lift their SNPs from hg18 to hg19
        sub_na=("SHIP" "TREND" "ADNI1" "IMAGEN1" "IMAGEN2" "IMAGEN3" "LOTHIAN")
        for da in ${sub_na[@]}; do

            # convert plink bim into plink map
            cut -f1-4 $o/$da/01.genotype/all.bim > $o/$da/01.genotype/all.hg18.map

            # rename original file to all.hg18
            mv $o/$da/01.genotype/all.bim $o/$da/01.genotype/all.hg18.bim
            mv $o/$da/01.genotype/all.bed $o/$da/01.genotype/all.hg18.bed
            mv $o/$da/01.genotype/all.fam $o/$da/01.genotype/all.hg18.fam

            python $ldir/LiftMap.py \
                -m $o/$da/01.genotype/all.hg18.map \
                -o $o/$da/01.genotype/bed-lift

            awk '{print $2,$1+0}' $o/$da/01.genotype/bed-lift.map > $o/$da/01.genotype/bed-lift.chrlist
            awk '{print $2,$4+0}' $o/$da/01.genotype/bed-lift.map > $o/$da/01.genotype/bed-lift.poslist

            $plink \
                --bfile $o/$da/01.genotype/all.hg18 \
                --update-chr $o/$da/01.genotype/bed-lift.chrlist \
                --update-map $o/$da/01.genotype/bed-lift.poslist \
                --extract $o/$da/01.genotype/bed-lift.chrlist \
                --make-bed \
                --out $o/$da/01.genotype/all

            rm $o/$da/01.genotype/bed-lift*
        done
    fi
fi

if [ true ]; then
# Move LOTHIAN from 1234 alleles to ACGT
    sub_na=("LOTHIAN")
    for da in ${sub_na[@]}; do
        # rename original file to all.hg18
        mv $o/$da/01.genotype/all.bim $o/$da/01.genotype/all.1234.bim
        mv $o/$da/01.genotype/all.bed $o/$da/01.genotype/all.1234.bed
        mv $o/$da/01.genotype/all.fam $o/$da/01.genotype/all.1234.fam

        # recode 1234 as ACGT
        $plink --bfile $o/$da/01.genotype/all.1234 --alleleACGT --make-bed --out $o/$da/01.genotype/all
    done
fi

if [  ]; then
# Convert phenotype files
    source extract_ADNI_phenotypes.sh
    source extract_SHIP_TREND_phenotypes.sh
fi

