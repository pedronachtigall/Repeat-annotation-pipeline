![repeat_annotation_pipeline](/repeat_annotation_workflow.png)

# Repeat-annotation-pipeline
Bioinformatics pipeline and tutorial for performing repeat annotation in genome assemblies using RepeatModeler and RepeatMasker. It also includes the commands needed to calculate the Repeat landscape (i.e., the CpG-adjusted kimura divergence versus the percent of the genome).

## Dependencies
 - [Python](https://www.python.org/) and [biopython](https://biopython.org/)
 - [RepeatModeler2](https://github.com/Dfam-consortium/RepeatModeler)
 - [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker)
 - [DeepTE](https://github.com/LiLabAtVT/DeepTE)
 - [TERL](https://github.com/muriloHoracio/TERL)
 - [bedtools](https://github.com/arq5x/bedtools2)
 - [SeqKit](https://github.com/shenwei356/seqkit)

Ensure that all dependencies are installed and working properly.

## Summary
 - [Run RepeatModeler](#run-repeatmodeler)
 - [Classify Unknown](#classify-unknown)
 - [Run RepeatMasker](#run-repeatmasker)
	- [Species-specific library](#species-specific-library)
 	- [Species-specific and curated libraries](#species-specific-and-curated-libraries)
 		- [One-round annotation](#one-round-annotation)
 		- [Multi-round annotation](#multi-round-annotation)

## Model species
We will use the Eastern diamondback rattlesnake (*Crotalus adamanteus*) as a model for this tutorial.

The genome is linked to the manuscript ["A Segregating Structural Variant Defines Novel Venom Phenotypes in the Eastern Diamondback Rattlesnake"]() published in *Molecular Biology and Evolution* and it is available in a [figshare repository](https://figshare.com/projects/Eastern_diamondback_rattlesnake_Crotalus_adamanteus_-_haplotype-resolved_genome_assembly/200614).

Download and decompress the genome data:
```
wget https://figshare.com/ndownloader/files/45450430 -O Cadamanteus_genome.fasta.gz
gzip -d Cadamanteus_genome.fasta.gz
```

##  Run RepeatModeler
Let's create a species-specific library using RepeatModeler.
```
mkdir RepeatModeler && cd RepeatModeler

ln -s ../Cadamanteus_genome.fasta Cadam_chr.fa

BuildDatabase -name Crotalus_Adamanteus -engine ncbi Cadam_chr.fa

RepeatModeler -pa 20 -engine ncbi -database Crotalus_Adamanteus
```
 - ```-pa``` is the number of threads. Adjust it properly.

:warning: This may take a while!

Let's add a prefix to sequence names and separate files into known and unknown based on RepeatModeler classification.
```
#add a prefix
cat Crotalus_Adamanteus-families.fa | seqkit fx2tab | awk '{ print "Cadamanteus_"$0 }' | seqkit tab2fx > Crotalus_Adamanteus-families.prefix.fa

#separate files into known and unknown
cat Crotalus_Adamanteus-families.prefix.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > Crotalus_Adamanteus-families.prefix.fa.known
cat Crotalus_Adamanteus-families.prefix.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > Crotalus_Adamanteus-families.prefix.fa.unknown

#known set is ready to go
awk '/^>/{print $1; next}{print}' Crotalus_Adamanteus-families.prefix.fa.known > Crotalus_Adamanteus-families.prefix.fa.known.FINAL
```

##  Classify unknown
Let's classify the unknown sequences output by RepeatModeler using two machine learning approaches: DeepTE and TERL.

We will use DeepTE to classify seqeunces using the metazoan models (this model is suitable for the snake species of this tutorial; modify it properly based on your target species).

The TERL will be used to remove false-positives. We will use the model DS3 bacause it is suitable for the snake species (Modify the model properly based on your target species).

Scripts ```CleanDeepTEheader.py``` and ```FilterTERL.py``` are available in this repository. Clone the repsitory or download each file separately.
```
wget https://raw.githubusercontent.com/pedronachtigall/Repeat-annotation-pipeline/refs/heads/main/scripts/CleanDeepTEheader.py
wget https://raw.githubusercontent.com/pedronachtigall/Repeat-annotation-pipeline/refs/heads/main/scripts/FilterTERL.py

python /path/to/DeepTE/DeepTE.py -d deepTE_out -o deepTE_out -i Crotalus_Adamanteus-families.prefix.fa.unknown -sp M -m_dir /path/to/DeepTE/models/Metazoans_model/
python CleanDeepTEheader.py deepTE_out/opt_DeepTE.fasta Crotalus_Adamanteus-families.prefix.fa.unknown.DeepTE

python /path/to/TERL/terl_test.py -m /path/to/TERL/Models/DS3/ -f Crotalus_Adamanteus-families.prefix.fa.unknown.DeepTE
mv TERL* Crotalus_Adamanteus-families.prefix.fa.unknown.DeepTE.TERL
python FilterTERL.py Crotalus_Adamanteus-families.prefix.fa.unknown.DeepTE.TERL Crotalus_Adamanteus-families.prefix.fa.unknown.DeepTE Crotalus_Adamanteus-families.prefix.fa.unknown.FINAL
```

##  Run RepeatMasker
Now that we have the species-specific library ready, we can perform the repeat annotation using RepeatMasker.

### Species-specific library
You can perform the repeat annotation using the species-specific library only.
```
#concatenate known and unknown
cat Crotalus_Adamanteus-families.prefix.fa.known.FINAL Crotalus_Adamanteus-families.prefix.fa.unknown.FINAL > Crotalus_Adamanteus-families.prefix.fa.known.unknown.FINAL

#run repeatmasker
RepeatMasker -pa 20 -gff -s -a -inv -no_is -norna -xsmall -nolow -div 40 -lib Crotalus_Adamanteus-families.prefix.fa.known.unknown.FINAL -cutoff 225 Cadam_chr.fa

#OPTIONAL - calculate the kimura distance
calcDivergenceFromAlign.pl -s Cadam_chr.fa.align.divsum Cadam_chr.fa.align
tail -n 72 Cadam_chr.fa.align.divsum > Cadam_chr.fa.Kimura.distance

#OPTIONAL - generate a modified GFF file, which I like more than the default generated by RepeatMasker.
wget https://raw.githubusercontent.com/pedronachtigall/Repeat-annotation-pipeline/refs/heads/main/scripts/AdjustingGFF_RM.py
python AdjustingGFF_RM.py Cadam_chr.fa.out.gff Crotalus_Adamanteus-families.prefix.fa.known.unknown.FINAL Cadam_chr.fa.out.adjusted.gff
```
 - ```-pa``` is the number of threads. Adjust it properly.
 - ```Cadam_chr.fa.out.adjusted.gff``` is your final annotation file.
 - ```Cadam_chr.fa.Kimura.distance``` is the CpG-adjusted Kimura distance file.

### Species-specific and curated libraries
If there are a curated set of repeat sequences for your target lineage available, it may be useful to include them in the annotation.

For the model species being used here, there are a curated set available from previous studies<sup>[Castoe et al., 2013](https://doi.org/10.1073/pnas.1314475110)</sup>.

Download and merge them to the species-specific library:
```
wget https://raw.githubusercontent.com/pedronachtigall/Repeat-annotation-pipeline/refs/heads/main/TElibs/BovB_CR1_TElib.fasta
wget https://raw.githubusercontent.com/pedronachtigall/Repeat-annotation-pipeline/refs/heads/main/TElibs/Boa_Snakes_Known_TElib_CLEAN.fasta
cat BovB_CR1_TElib.fasta Boa_Snakes_Known_TElib_CLEAN.fasta > Snakes_Curated_TElib.fasta
cat Crotalus_Adamanteus-families.prefix.fa.known.FINAL Crotalus_Adamanteus-families.prefix.fa.unknown.FINAL BovB_CR1_TElib.fasta Boa_Snakes_Known_TElib_CLEAN.fasta > Cadam_RM2Curated_TElib.fasta
```

#### One-round annotation
You can do the annotation in one round (i.e., using the merged file):
```
#run repeatmasker
RepeatMasker -pa 20 -gff -s -a -inv -no_is -norna -xsmall -nolow -div 40 -lib Cadam_RM2Curated_TElib.fasta -cutoff 225 Cadam_chr.fa

#OPTIONAL - calculate the kimura distance
calcDivergenceFromAlign.pl -s Cadam_chr.fa.align.divsum Cadam_chr.fa.align
tail -n 72 Cadam_chr.fa.align.divsum > Cadam_chr.fa.Kimura.distance

#OPTIONAL - generate a modified GFF file, which I like more than the default generated by RepeatMasker.
wget https://raw.githubusercontent.com/pedronachtigall/Repeat-annotation-pipeline/refs/heads/main/scripts/AdjustingGFF_RM.py
python AdjustingGFF_RM.py Cadam_chr.fa.out.gff Cadam_RM2Curated_TElib.fasta Cadam_chr.fa.out.adjusted.gff```
```
 - ```-pa``` is the number of threads. Adjust it properly.
 - ```Cadam_chr.fa.out.adjusted.gff``` is your final annotation file.
 - ```Cadam_chr.fa.Kimura.distance``` is the CpG-adjusted Kimura distance file.

#### Multi-round annotation
OR, you can do the annotation in multiple rounds, which you use one file at a time considering the masked genome of the previous run:
```
mkdir RL_FINAL && cd RL_FINAL

ln -s ../../Cadam_genome.fasta Cadam_chr.fa

mkdir 01_simple_out 02_snakes_out 03_known_out 04_unknown_out

#round 1: annotate/mask simple repeats
RepeatMasker -pa 20 -a -e ncbi -dir 01_simple_out -noint -xsmall Cadam_chr.fa

#round 2: annotate/mask using a curated lib and the masked genome from round 1
RepeatMasker -pa 20 -a -e ncbi -dir 02_snakes_out -nolow -lib ../Snakes_Curated_TElib.fasta 01_simple_out/Cadam_chr.simple_mask.masked.fasta

#round 3: annotate/mask known RM2 elements sourced from species-specific de novo repeat library using masked genome from 2nd round
RepeatMasker -pa 20 -a -e ncbi -dir 03_known_out -nolow -lib ../Crotalus_Adamanteus-families.prefix.fa.known.FINAL 02_snakes_out/Cadam_chr.snakes_mask.masked.fasta

#round 4 - annotate/mask unkown RM2 (classified using DeepTE and filtered using TERL) using the masked genome from round 3
RepeatMasker -pa 20 -a -e ncbi -dir 04_unknown_out -nolow -lib ../Crotalus_Adamanteus-families.prefix.fa.unknown.FINAL 03_known_out/Cadam_chr.known_mask.masked.fasta

#round 5: merge round 1, 2, 3, and 4 annotations

mkdir -p 05_full_out

cat 01_simple_out/Cadam_chr.simple_mask.cat.gz \
02_snakes_out/Cadam_chr.snakes_mask.cat.gz \
03_known_out/Cadam_chr.known_mask.cat.gz \
04_unknown_out/Cadam_chr.unknown_mask.cat.gz \
> 05_full_out/Cadam_chr.full_mask.cat.gz

#resummarize repeat compositions from combined analysis of all RepeatMasker rounds
cd 05_full_out/
ProcessRepeats -gff -a -species tetrapoda Cadam_chr.full_mask.cat.gz

#OPTIONAL - calculate the kimura distance
calcDivergenceFromAlign.pl -s Cadam_chr.full_mask.align.divsum Cadam_chr.full_mask.align
tail -n 72 Cadam_chr.full_mask.align.divsum > Cadam_chr.full_mask.Kimura.distance

#OPTIONAL - generate a modified GFF file, which I like more than the default generated by RepeatMasker.
wget https://raw.githubusercontent.com/pedronachtigall/Repeat-annotation-pipeline/refs/heads/main/scripts/AdjustingGFF_RM.py
python ~/programs/TElibs/AdjustingGFF_RM.py Cadam_chr.full_mask.out.gff ../../RepeatModeler/Cadam_RM2Curated_TElib.fasta Cadam_chr.full_mask.out.adjusted.gff

#soft mask the genome
bedtools maskfasta -soft -fi ../../Cadam_primary_chromosomes.fasta -bed Cadam_chr.full_mask.out.adjusted.gff -fo Cadam_primary_chromosomes.soft.masked.fasta

#hard mask the genome
bedtools maskfasta -fi ../../Cadam_primary_chromosomes.fasta -bed Cadam_chr.full_mask.out.adjusted.gff -fo Cadam_primary_chromosomes.hard.masked.fasta
```
 - ```-pa``` is the number of threads. Adjust it properly.
 - ```Cadam_chr.fa.out.adjusted.gff``` is your final annotation file.
 - ```Cadam_chr.fa.Kimura.distance``` is the CpG-adjusted Kimura distance file.
 - ```Cadam_primary_chromosomes.soft.masked.fasta``` is the soft-masked genome file.
 - ```Cadam_primary_chromosomes.hard.masked.fasta``` is the hard-masked genome file.

## References
If you use this tutorial or any of the resources/scripts, please consider citing: [Nachtigall et al., in prep]().

If you use the curated set of repeat sequences from snakes, cite the original manuscript: [Castoe et al., 2013](https://doi.org/10.1073/pnas.1314475110).

Please, cite the original manuscript of each tool used in this tutorial.
