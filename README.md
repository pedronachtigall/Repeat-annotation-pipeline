# Repeat-annotation-pipeline
Bioinformatics pipeline and tutorial for performing repeat annotation in genome assemblies using RepeatModeler and RepeatMasker.

## Dependencies
 - [Python](https://www.python.org/) and [biopython](https://biopython.org/)
 - [RepeatModeler2](https://github.com/Dfam-consortium/RepeatModeler)
 - [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker)
 - [DeepTE](https://github.com/LiLabAtVT/DeepTE)
 - [TERL](https://github.com/muriloHoracio/TERL)
 - [bedtools](https://github.com/arq5x/bedtools2)
 - [SeqKit](https://github.com/shenwei356/seqkit)

Ensure that all dependencies are working properly.

## Summary
 - [Run RepeatModeler](#run-repeatmodeler)
 - [Classify Unknown](#classify-unknown)
 - [Run RepeatMasker](#run-repeatmasker)
	- [Species-specific library](#species-specific-library)
 	- [Species-specific and curated libraries](#species-specific-and-curated-libraries)

## Model species
We will use the Eastern diamondback rattlesnake (*Crotalus adamanteus*) as a model for this tutorial.

The genome is linked to the manuscript ["A Segregating Structural Variant Defines Novel Venom Phenotypes in the Eastern Diamondback Rattlesnake"]() published in *Molecular Biology and Evolution* and it is available in a [figshare repository](https://figshare.com/projects/Eastern_diamondback_rattlesnake_Crotalus_adamanteus_-_haplotype-resolved_genome_assembly/200614).

Download and uncompress the genome data:
```
wget https://figshare.com/ndownloader/files/45450430 -O Cadamanteus_genome.fasta.gz
gzip -d Cadamanteus_genome.fasta.gz
```

##  Run RepeatModeler

##  Classify unknown

##  Run RepeatMasker

### Species-specific library

### Species-specific and curated libraries

:construction:	**Under construction!** :construction:	

<!--

#Run RepeatModeler

mkdir RepeatModeler && cd RepeatModeler

ln -s ../Cadam_primary_chromosomes.fasta Cadam_chr.fa

BuildDatabase -name Crotalus_Adamanteus -engine ncbi Cadam_chr.fa

RepeatModeler -pa 40 -engine ncbi -database Crotalus_Adamanteus

cat Crotalus_Adamanteus-families.fa | seqkit fx2tab | awk '{ print "Cadamanteus_"$0 }' | seqkit tab2fx >Crotalus_Adamanteus-families.prefix.fa

#separate files into known and Unknown
cat Crotalus_Adamanteus-families.prefix.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > Crotalus_Adamanteus-families.prefix.fa.known
cat Crotalus_Adamanteus-families.prefix.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > Crotalus_Adamanteus-families.prefix.fa.unknown


awk '/^>/{print $1; next}{print}' Crotalus_Adamanteus-families.prefix.fa.known > Crotalus_Adamanteus-families.prefix.fa.known.FINAL


 #conda activate deepTE
python /home/pgn22a/programs/DeepTE/DeepTE.py -d deepTE_out -o deepTE_out -i Crotalus_Adamanteus-families.prefix.fa.unknown -sp M -m_dir /home/pgn22a/programs/DeepTE/models/Metazoans_model/
python ~/programs/DeepTE/CleanDeepTEheader.py deepTE_out/opt_DeepTE.fasta Crotalus_Adamanteus-families.prefix.fa.unknown.DeepTE

python /home/pgn22a/programs/TERL/terl_test.py -m /home/pgn22a/programs/TERL/Models/DS3/ -f Crotalus_Adamanteus-families.prefix.fa.unknown.DeepTE

mv TERL* Crotalus_Adamanteus-families.prefix.fa.unknown.DeepTE.TERL

python /home/pgn22a/programs/TERL/FilterTERL.py Crotalus_Adamanteus-families.prefix.fa.unknown.DeepTE.TERL Crotalus_Adamanteus-families.prefix.fa.unknown.DeepTE Crotalus_Adamanteus-families.prefix.fa.unknown.FINAL

#generate the RM2_CuratedSnakes

cat Crotalus_Adamanteus-families.prefix.fa.known.FINAL Crotalus_Adamanteus-families.prefix.fa.unknown.FINAL ~/programs/TElibs/BovB_CR1_TElib.fasta ~/programs/TElibs/Boa_Snakes_Known_TElib_CLEAN.fasta > Cadam_RM2Curated_TElib.fasta

#run RepeatMasker with the RM2Curated_TElib

RepeatMasker -pa 20 -gff -s -a -inv -no_is -norna -xsmall -nolow -div 40 -lib Cadam_RM2Curated_TElib.fasta -cutoff 225 Cadam_chr.fa
export PERL5LIB=/home/pgn22a/miniconda3/envs/EDTA/share/RepeatMasker/
calcDivergenceFromAlign.pl -s Cadam_chr.fa.align.divsum Cadam_chr.fa.align
tail -n 72 Cadam_chr.fa.align.divsum > Cadam_chr.fa.Kimura.distance

python ~/programs/TElibs/AdjustingGFF_RM.py Cadam_chr.fa.out.gff Cadam_RM2Curated_TElib.fasta Cadam_chr.fa.out.adjusted.gff


#running the seriated pipeline (simple+snakes+known+unknown)

mkdir RL_FINAL && cd RL_FINAL

#ln -s ../Cadam_primary_chromosomes.fasta Cadam_chr.fa

mkdir 01_simple_out 02_snakes_out 03_known_out 04_unknown_out

#round 1: annotate/mask simple repeats
RepeatMasker -pa 20 -a -e ncbi -dir 01_simple_out -noint -xsmall Cadam_chr.fa

 #adjust names of outputs

#round 2: annotate/mask using a curated lib and the output from round 1
RepeatMasker -pa 20 -a -e ncbi -dir 02_snakes_out -nolow -lib ~/programs/TElibs/Snakes_Curated_TElib.fasta 01_simple_out/Cadam_chr.simple_mask.masked.fasta

  #nohup ProcessRepeats -xsmall -nolow -gff -a -species tetrapoda Cadam_chr.snakes_mask.cat.gz &
  #bedtools maskfasta -soft -fi ../01_simple_out/Cadam_chr.simple_mask.masked.fasta -bed Cadam_chr.snakes_mask.out.gff -fo Cadam_chr.snakes_mask.masked.fasta

#round 3: annotate/mask known RM2 elements sourced from species-specific de novo repeat library using output from 2nd round of RepeatMasker
RepeatMasker -pa 20 -a -e ncbi -dir 03_known_out -nolow -lib Crotalus_Adamanteus-families.prefix.fa.known.FINAL 02_snakes_out/Cadam_chr.snakes_mask.masked.fasta

#round 4 - annotate/mask unkoen RM2 (classified using DeepTE and TERL) using the output from round 3
RepeatMasker -pa 20 -a -e ncbi -dir 04_unknown_out -nolow -lib Crotalus_Adamanteus-families.prefix.fa.unknown.FINAL 03_known_out/Cadam_chr.known_mask.masked.fasta

#round 5: merge round 1, 2, 3, and 4 annotations

mkdir -p 05_full_out

cat 01_simple_out/Cadam_chr.simple_mask.cat.gz \
02_snakes_out/Cadam_chr.snakes_mask.cat.gz \
03_known_out/Cadam_chr.known_mask.cat.gz \
04_unknown_out/Cadam_chr.unknown_mask.cat.gz \
> 05_full_out/Cadam_chr.full_mask.cat.gz


# resummarize repeat compositions from combined analysis of all RepeatMasker rounds
cd 05_full_out/
ProcessRepeats -gff -a -species tetrapoda Cadam_chr.full_mask.cat.gz

export PERL5LIB=/home/pgn22a/miniconda3/envs/EDTA/share/RepeatMasker/
calcDivergenceFromAlign.pl -s Cadam_chr.full_mask.align.divsum Cadam_chr.full_mask.align
tail -n 72 Cadam_chr.full_mask.align.divsum > Cadam_chr.full_mask.Kimura.distance

python ~/programs/TElibs/AdjustingGFF_RM.py Cadam_chr.full_mask.out.gff ../../RepeatModeler/Cadam_RM2Curated_TElib.fasta Cadam_chr.full_mask.out.adjusted.gff

#soft mask the genome

bedtools maskfasta -soft -fi ../../Cadam_primary_chromosomes.fasta -bed Cadam_chr.full_mask.out.adjusted.gff -fo Cadam_primary_chromosomes.soft.masked.fasta


-->
