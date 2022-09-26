SCRIPTS:Venom variation and evolution in populations of montane
pitvipers (*Viperidae* :*Cerrophidion*).
================
Ramses A. Rosales-Garcia, Rhett M. Rautsaw , Erich P. Hofmann, Christoph
I. Grunwald, Jason M. Jones, Hector Franz-Chavez, Ivan T.
Ahumada-Carrillo, Ricardo Ramirez-Chaparro, Miguel Angel De la
Torre-Loranca, Jason L. Strickland, Andrew J. Mason, Matthew L. Holding,
Miguel Borja, Gamaliel Castaneda-Gaytan, Darin R. Rokyta, Tristan D.
Schramer, N. Jade Mellor, Edward A. Myers, Christopher Parkinson
2022 September 26

-   [Assembly](#assembly)
    -   [Pre Assembly](#pre-assembly)
    -   [Assembly](#assembly-1)
        -   [Extender](#extender)
        -   [Ngen](#ngen)
        -   [Trinity](#trinity)
    -   [Post assembly modifications](#post-assembly-modifications)
        -   [Change names and merge](#change-names-and-merge)
-   [Annotation](#annotation)
    -   [Toxcodan](#toxcodan)
        -   [Toxins](#toxins)
        -   [Nontoxins](#nontoxins)
    -   [Manual annotation](#manual-annotation)
        -   [Blast](#blast)
        -   [Annotation with scripts](#annotation-with-scripts)
        -   [Annotation in Geneious](#annotation-in-geneious)
        -   [move back to palmetto](#move-back-to-palmetto)
-   [Post Annotation Filters](#post-annotation-filters)
    -   [Concatenate and Clean](#concatenate-and-clean)
    -   [Chimera Killer](#chimera-killer)
-   [Consense Transcriptomes](#consense-transcriptomes)
    -   [Cluster contigs](#cluster-contigs)
    -   [Filtering Transcriptomes (CDS\_filter and
        RSEM\_combiner)](#filtering-transcriptomes-cds_filter-and-rsem_combiner)
    -   [Clean consense](#clean-consense)
    -   [RSEM, Plots, and Differential
        expression](#rsem-plots-and-differential-expression)
-   [PLA2s Phylogeny](#pla2s-phylogeny)
    -   [Get PLA2s secuences](#get-pla2s-secuences)
-   [Sequence Diversity 1](#sequence-diversity-1)
    -   [Busco and Phylogenetics](#busco-and-phylogenetics)
        -   [Busco](#busco)
        -   [Prepare sequences](#prepare-sequences)
        -   [All genes](#all-genes)
            -   [Alignment, Cleaning and
                Triminig](#alignment-cleaning-and-triminig)
            -   [Gene trees and Species
                trees](#gene-trees-and-species-trees)
        -   [Core-Genes](#core-genes)
            -   [Subset core-genes](#subset-core-genes)
            -   [Concatenated sequences](#concatenated-sequences)
        -   [Scaling tree with IQtree](#scaling-tree-with-iqtree)
-   [Sequence Diversity 2](#sequence-diversity-2)
    -   [Variant Calling](#variant-calling)
    -   [Selection Analysis](#selection-analysis)
        -   [filter coverage](#filter-coverage)
        -   [Run Analysis(Tajimas’ D, *F*<sub>*S**T*</sub>, *π* and
            snpEff)](#run-analysistajimas-d-f_st-pi-and-snpeff)
            -   [snpEff reference](#snpeff-reference)
            -   [Run Analysis](#run-analysis)
        -   [Run Analysis(HyPhy BUSTED)](#run-analysishyphy-busted)
            -   [Clip tree](#clip-tree)
            -   [Run Analysis](#run-analysis-1)
    -   [Prepare Final Data](#prepare-final-data)
        -   [Correct tables](#correct-tables)
        -   [Merge tables](#merge-tables)
-   [Final](#final)

# Assembly

## Pre Assembly

I started with the raw sequences and a list of the ID of the files named
list1

loop to generate the concatenated files R1 forward read and R2 reverse
read

``` bash
for i in `cat list1`; do mkdir -p $i/00raw;mv $i* $i/00raw/;
mkdir $i/01concat;
cat $i/00raw/*_R1_*.fastq.gz > $i/01concat/${i}_R1.fastq.gz;
cat $i/00raw/*_R2_*.fastq.gz > $i/01concat/${i}_R2.fastq.gz;
done
```

Make a “bio” conda environment for next steps

``` bash
conda create -n bio
conda activate bio
conda install parallel trim-galore fastqc biopython
conda install -y biopython bamtools bedtools blast bowtie2 bwa cd-hit emboss fastqc gatk4 jellyfish parallel pear picard pigz rsem samtools sra-tools trim-galore
```

run fastqc to check for quality of the sequences run trim galore to trim
the sequences with low quality and the idex sequences run fastqc again
to check the trimed sequences

``` bash
parallel -a list1 -j 6 --verbose "cd {}
cd 01concat 
fastqc *.fastq.gz
cd ..
mkdir 02trim
trim_galore --paired --phred33 --length 75 -q 5 --stringency 1 -e 0.1 -o 02trim 01concat/{}_R1.fastq.gz 01concat/{}_R2.fastq.gz &> 02trim/{}_tg.log
cd 02trim 
mv {}_R1_val_1.fq.gz {}_R1_trim.fastq.gz  
mv {}_R2_val_2.fq.gz {}_R2_trim.fastq.gz
fastqc *.fastq.gz"
```

Use pear to merge forward and reverse reads

``` bash
parallel -a list1 -j 3 --verbose "echo {}
cd {}
mkdir 03merged
pear -j 4 -f 02trim/{}_R1_trim.fastq.gz -r 02trim/{}_R2_trim.fastq.gz -o 03merged/{}_pear > 03merged/{}_pear.log
cd 03merged
pigz *.fastq
pigz -d {}_pear.assembled.fastq.gz"
```

## Assembly

We used 3 assembly softwares, the seamples were in palmetto HPC cluster
which uses PVS protocol.

``` bash
qsub <pbs_script>
```

### Extender

``` bash
#PBS -N Extender
#PBS -l select=13:ncpus=16:mem=60gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j ose

source .bash_profile

module load anaconda3/5.1.0-gcc/8.3.1
source activate bio

cd $PBS_O_WORKDIR

parallel -a list1 --sshloginfile $PBS_NODEFILE -j1 "source .bash_profile;
 module load anaconda3/5.1.0-gcc/8.3.1;
 source activate bio;
 cd /zfs/venom/Ramses/Cerrophidion;
 cd {};
    cd 03merged
    pigz -d *fastq.gz
    cd ..
 mkdir 04extender;
 cd 04extender;
    Extender3.py -r ../03merged/{}_pear.assembled.fastq -s 1000 -o 120 -msq 30 -mrq 20 -reps 20 -p 0.20 -e 2 -np 20 > {}_extender.log"
```

### Ngen

``` bash
#!/bin/bash
#
#PBS -N Ngen
#PBS -l select=13:ncpus=16:mem=60gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

source .bash_profile

module load anaconda3/5.1.0-gcc/8.3.1
source activate bio

cd $PBS_O_WORKDIR

parallel -a list1 --sshloginfile $PBS_NODEFILE -j1 "source .bash_profile
    module load anaconda3/5.1.0-gcc/8.3.1
    source activate bio 
    cd /zfs/venom/Ramses/Cerrophidion
    cd {}
    mkdir 05Ngen
    cd 05Ngen
    cp /zfs/venom/Ramses/bin/NGen/05Ngen.script .
    sed 's/NAME/{}/g' 05Ngen.script > 05Ngen.script2
    rm 05Ngen.script
    mv 05Ngen.script2 05Ngen.script
    xng 05Ngen.script &> {}_05Ngen.log"
```

the control script is this, it has to be modified with the right path
&lt;05Ngen.script&gt;

``` bash
#!/zfs/venom/Ramses/bin/NGen/usr/bin/xng
set $project_name: "NAME_NGen14"
project kind: transcriptome
workflow kind: de_novo
diskPath path: { "/zfs/venom/Ramses/Cerrophidion/NAME/05Ngen/" }
set $TempDisk: "/zfs/venom/Ramses/Cerrophidion/NAME/05Ngen/TMP"
set $ResultDisk: "/zfs/venom/Ramses/Cerrophidion/NAME/05Ngen/"
set $mersize: 21
set $matchQuery: 50
set $matchTemplate: 80
set $trim: false
set $linkerFile: ""
set $scan: false
set $contaminateFile: ""
set $deleteIntermediates: true
; set up directory structure
set $project_root:       "${ResultDisk}/${project_name}.Transcriptome"
set $intermediatePath:     "${project_root}/Intermediate_Assembly_Results"
set $linkerFile_default:    ""
set $trim_default:      false
; no contam
set $scan_default:      false
set $clusterPath:      "${intermediatePath}/cluster"
set $intermediateFiles:     "${TempDisk}/intermediateFiles/${project_name}"
set $assemblyPath:      "${project_root}/Assemblies"
set $maxSeqs_default:     -1
set $minClusterSizeToOutput_default:  100
set $maxClusterSizeToOutput_default:  100000
set $forceCluster_default:    true
set $deleteIntermediates_default:  true
set $ignorePolyMers_default:   true
setDefaultDirectory defaultDirectory:   "${project_root}"
; directory structure:
;------------------------------------------
;${project_name}.Transcriptome
; Assemblies
;  ${project_name}_AllUnassembled.fastq
;  sub_0
; Reports
;  ${project_name}_AllTranscripts.searchresults
; Transcripts
;  ${project_name}_NGen14.fasta
; Intermediate Assembly Results
;  intermediateFiles
;  cluster       (cluster all query reads and assemble them)
;   rawStacks
;   scripts
;   seqs
;   unassembled
;   unclustered
;   SNGCluster_Stats.txt
;------------------------------------------
;------------------------------------------
; cluster the seqs and use SNG to assemble them
@stackCluster
rnaAssemble
  query: {
            {
                file : "/zfs/venom/Ramses/Cerrophidion/NAME/03merged/NAME_pear.assembled.fastq"
                SeqTech: IlluminaLongReads
            }
   trim: ${trim}
   sngTrim: {
    ; Vector trimming
    setVectorParam  EndCutOff:  130
        MatchSize:  11
        MinTrimLength: 15

    ;; xng will append the reads and destination params to this script
    ; xng will set $queryTrimOutput
    TrimVector
     LinkerFile: "${linkerFile}"
     ;
     destination: "$${queryTrimOutput}"
     ;reads: { }
   }
   scan: ${scan}
   contaminateScan: {
    assembleTemplate template: {

     }
     directoryQueryMer:  "intermediateFiles"
     directoryTemplateMer:  "${ResultDisk}/templateMers/MerSize_${mersize}"
     hits:     "intermediateFiles/$${project}.hits"
     layout:     "intermediateFiles/$${project}.layout"
     output:     "$${queryTrimOutput}"
     unassembled:   "$${queryTrimOutput}_unassembled.fastq"
     results:    "$${queryTrimOutput}_results.txt"
     format:     none

     deleteIntermediates: true
     ignorePolyMers:   ${ignorePolyMers}
   }
  }
;-----------  ;
  maxSeqs: ${maxSeqs}
  merSize: ${mersize}
  merSkipQuery: 0
  deleteIntermediates: ${deleteIntermediates}
  forceMake: ${forceCluster}
  directoryQueryMer: "${intermediateFiles}"
  output:    "${clusterPath}"
  assemblyFolder:  "${assemblyPath}"

  clusterParam: {

   minNewClusterSize: 5
   minSingleMergeClusterSize: 7
   minMultiMergeClusterSize : 7
   minMultiMergeIgnoreFactor: 3.0
   maxMergeSize: 50
   minClusterSizeToOutput: ${minClusterSizeToOutput}
   maxClusterSizeToOutput: ${maxClusterSizeToOutput}
   ignorePolyMers:   ${ignorePolyMers}
  }
  ; how to run sng for all of the clusters that are built
  ; xng will set $seqfile, $rootName, $clusterID
  sngScript: {
     setDefaultDirectory defaultDirectory: "."
     setParam merLength: ${mersize}
        minmatchPercent: 97
        useRepeatHandling: false
        minContigSeqs: 101
     ;setPairSpecifier pairs:{ { forward:".*_([0-9]+)_f" reverse:".*_([0-9]+)_r" min:500 max:10000 } }
     setAssemblyReport file: "$${resultName}" name:"cl_$${clusterID}"
     loadSeq file: "$${seqfile}"
     assemble

     ; assemble the existing contigs, todo remove unassembled seqs
     setParam minmatchPercent:85
     assemble
     WriteUnassembledSeqs file: "$${unassembledName}-unassembled.fastq"
     nameContigs name: "cl_$${clusterID}" contigs: all

     saveProject  file: "$${rootName}.sqd" openInSeqMan: false
     saveProject  file: "$${rootName}.fas" format: fastaNoQual openInSeqMan: false

     closeProject
    }
  sngAlign: true

;------------------------------------------
; generate the final transcript reports
rnaClusterReports
 output:   "Reports"
 assemblyFolder: "${assemblyPath}"

 trimReport:  "${stackCluster.trimResults}"
 contamReport:   "${stackCluster.contaminateResults}"
 clusterReport: "${clusterPath}/results.txt"
quit
```

### Trinity

``` bash
#!/bin/bash
#
#PBS -N Trinity
#PBS -l select=1:ncpus=20:interconnect=fdr:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

source .bash_profile

module load anaconda3/5.1.0-gcc/8.3.1
source activate trinity_env

cd $PBS_O_WORKDIR

parallel -a list1 --sshloginfile $PBS_NODEFILE -j1 "source .bash_profile
    module load anaconda3/5.1.0-gcc/8.3.1
    source activate trinity_env
    cd /zfs/venom/Ramses/Cerrophidion
    cd {}
    mkdir 06trinity
    cd 06trinity
    Trinity --seqType fq --CPU 18 --min_contig_length 200 --max_memory 95G --full_cleanup --output {}_trinity --single ../03merged/{}_pear.assembled.fastq"
```

## Post assembly modifications

### Change names and merge

change contig names with awk merge all contigs remove redundant contigs
with cd-hit

``` bash
#PBS -N Concatenate
#PBS -l select=1:ncpus=8:mem=8gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe


source .bash_profile

anaconda
bio 

cd $PBS_0_WORKDIR

for i in `cat list1`
do cd $i
mkdir 07assembly
awk '/^>/{print ">extenderContig" ++i; next}{print}' < 04extender/Final_Extender_contigs.fasta >  07assembly/${i}_extender.fasta
awk '/^>/{print ">ngenContig" ++i; next}{print}' < 05Ngen/${i}_NGen14.Transcriptome/Transcripts/${i}_NGen14_NovelTranscripts.fas>  07assembly/${i}_ngen.fasta
awk '/^>/{print ">trinityContig" ++i; next}{print}' < 06trinity/${i}_trinity/*.fasta >  07assembly/${i}_trinity.fasta
cd 07assembly
cat ${i}_extender.fasta ${i}_ngen.fasta ${i}_trinity.fasta > ${i}_assembly.fasta
cd-hit-est -i ${i}_assembly.fasta -o ${i}_assembly_reduced.fasta -d 0 -c 1.0
done
```

we had to increase memory use in cd-hit for some of the files, -M

for one sequence we had to recompile cd-hit and increase the memory
allowed

Warning: Some seqs are too long, please rebuild the program with make
parameter MAX\_SEQ=new-maximum-length (e.g. make MAX\_SEQ=10000000)

to do this, download cd-hit from github and follow the instructions to
compile and change maximun lenght.

# Annotation

## Toxcodan

Toxcodan requires Signalp-4.1 and biopython=1.76

### Toxins

Annotation of Toxins

``` bash
#PBS -N venomancer
#PBS -l select=12:ncpus=20:interconnect=fdr:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.com
#PBS -m abe
#PBS -j oe

source .bash_profile

module load anaconda3/5.1.0-gcc/8.3.1
source activate venomancer_env

cd $PBS_O_WORKDIR

parallel -a list1 --sshloginfile $PBS_NODEFILE -j1 " source .bash_profile
module load anaconda3/5.1.0-gcc/8.3.1
source activate venomancer_env
cd /zfs/venom/Ramses/Cerrophidion
cd {}
mkdir 12venomancer
toxcodan.py -s {} -t 07assembly/{}_assembly_reduced.fasta -o 12venomancer -m /zfs/venom/Ramses/bin/ToxCodAn/models -c 18
cd 12venomancer
cat {}_Toxins_cds_Redundancyfiltered.fasta {}_PutativeToxins_cds_SPfiltered.fasta > {}_Toxins.fasta
perl -p -i -e 's/>/>TOXIN_/g' {}_Toxins.fasta"
```

### Nontoxins

you need to download the busco odb and pfam models before run this
chunk, for more information see toxodan guide in github

run hmmpress in Pfam moddel, this is not mentioned in toxcodan guide

``` bash
cd /zfs/venom/Ramses/db/Pfam-A.hmm
hmmpress Pfam-A.hmm
```

``` bash
#PBS -N venomancerNT
#PBS -l select=13:ncpus=20:interconnect=fdr:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.com
#PBS -m abe
#PBS -j oe

source .bash_profile


module load anaconda3/5.1.0-gcc/8.3.1
source activate venomancer_env

cd $PBS_O_WORKDIR

parallel -a list1 --sshloginfile $PBS_NODEFILE -j1 " source .bash_profile
module load anaconda3/5.1.0-gcc/8.3.1
source activate venomancer_env
cd /zfs/venom/Ramses/Cerrophidion
cd {}
mkdir 13nontoxins
cp 12venomancer/{}_NonToxins_contigs.fasta ./13nontoxins/
cd 13nontoxins
codan.py -t {}_NonToxins_contigs.fasta -m /zfs/venom/Ramses/bin/CodAn/models/VERT_full -o {}_NonToxins_codan -18
mv {}_NonToxins_codan/ORF_sequences.fasta {}_NonToxins_CDS.fasta
NonToxinsAnnotation.py -t {}_NonToxins_CDS.fasta -d /zfs/venom/Ramses/db/swissprot -b /zfs/venom/Ramses/db/tetrapoda_odb10 -p /zfs/venom/Ramses/db/Pfam-A.hmm -c 18
mv Annotation_output/annotated.fa {}_NonToxins_annotated.fasta"
```

## Manual annotation

### Blast

get unit prot data base

``` bash
#download swiss prot data base
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

#use pigz to decompres file
anaconda
bio
pigz -d *.fasta.gz

#then make db, is AA fasta, so is protein
makeblastdb -in uniprot_sprot.fasta -dbtype prot
```

run blast

``` bash
#PBS -N blast
#PBS -l select=1:ncpus=16:mem=60gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

source .bash_profile

module load anaconda3/5.1.0-gcc/8.3.1
source activate bio

cd $PBS_O_WORKDIR

parallel -a list1 --sshloginfile $PBS_NODEFILE -j1 "source .bash_profile
    module load anaconda3/5.1.0-gcc/8.3.1
    source activate bio
    cd /zfs/venom/Ramses/Cerrophidion
    cd {}
    mkdir 08blast
    cd 08blast
    blastx -query ../07assembly/{}_assembly_reduced.fasta -db /zfs/venom/Ramses/db/uniprot_sprot.fasta -outfmt 5 -num_threads 13 -evalue 0.00001 -out {}_blastUPT.xml"
```

### Annotation with scripts

Blast parser

``` bash
#PBS -N blast_parse
#PBS -l select=13:ncpus=16:mem=60gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

source .bash_profile

module load anaconda3/5.1.0-gcc/8.3.1
source activate bio

cd $PBS_O_WORKDIR

parallel -a list1 --sshloginfile $PBS_NODEFILE -j1 "source .bash_profile
    module load anaconda3/5.1.0-gcc/8.3.1
    source activate bio
    cd /zfs/venom/Ramses/Cerrophidion
    cd {}
    cd 08blast
    /zfs/venom/Ramses/bin/RokytaScripts/Blast_parse_v6.py ../07assembly/{}_assembly_reduced.fasta {}_blastUPT.xml 90
    mv keywords.fasta {}_Toxins.fasta; mv nontoxins.fasta {}_Nontoxins.fasta"
```

Autoannotator

``` bash
#PBS -N autoannotate
#PBS -l select=13:ncpus=16:mem=60gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

source .bash_profile

module load anaconda3/5.1.0-gcc/8.3.1
source activate bio

cd $PBS_O_WORKDIR

parallel -a list1 --sshloginfile $PBS_NODEFILE -j1 "source .bash_profile
    module load anaconda3/5.1.0-gcc/8.3.1
    source activate bio
    cd /zfs/venom/Ramses/Cerrophidion
    cd {}
    mkdir 09autoannotate
    cd 09autoannotate
    /zfs/venom/Ramses/bin/RokytaScripts/AutoAnnotator.py -i ../08blast/{}_Toxins.fasta -r zfs/venom/Ramses/bin/RokytaScripts/Annotation_Databases/Snake_toxins_CDS.fasta -c cd-hit.est -ri Toxin -pf {} -p 0.80 -s signalp
    /zfs/venom/Ramses/bin/RokytaScripts/Autoannotator.py -i ../08blast/{}_Nontoxins.fasta -r zfs/venom/Ramses/bin/RokytaScripts/Annotation_Databases/ChorrA_nontoxins_CDS.fasta -c cd-hit.est -ri Nontoxin -pf {} -p 0.80 -s signalp
    python /zfs/venom/Ramses/bin/RokytaScripts/AutoAnnotatorFormat.py -i Annotated_{}_Toxins.fasta -f seqman2
    python /zfs/venom/Ramses/bin/RokytaScripts/AutoAnnotatorFormat.py -i Annotatted_{}_Nontoxins.fasta -f seqman2"
```

``` bash
#problem with the nontoxins with Autoannotator.
#the .pbs job did not work, however it worked with a interactive job and using for loop
interactive
anaconda
bio
Cerrophidion

for i in `cat list1`
do echo $i
cd $i/09autoannotate
python /zfs/venom/Ramses/bin/RokytaScripts/AutoAnnotatorFormat.py -i Annotated_*_Toxins.fasta -f seqman2
python /zfs/venom/Ramses/bin/RokytaScripts/AutoAnnotatorFormat.py -i Annotated_*_Nontoxins.fasta -f seqman2
cd ../..
done 
```

NextAnnotate

``` bash
#PBS -N nexannotate
#PBS -l select=13:ncpus=16:mem=60gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

source .bash_script

module load anaconda3/5.1.0-gcc/8.3.1
source activate bio

cd $PBS_O_WORKDIR

parallel -a list1 --sshloginfile $PBS_NODEFILE -j1 "source .bash_profile
    module load anaconda3/5.1.0-gcc/8.3.1
    source activate bio
    cd /zfs/venom/Ramses/Cerrophidion
    cd {}
    mkdir 10nextannotate
    cd 10nextannotate
    /zfs/venom/Ramses/bin/AndrewScripts/NextAnnotate_v0.3.py -i ../09autoannotate/Unannotated_{}_Toxins.fasta -x ../08blast/{}_blastUPT.xml -pf {}"
```

``` bash
interactive 
anaconda
bio
Cerrophidion

for i in `cat list1`
do cd $i/10nextannotate/NextAnnotator
python /zfs/venom/Ramses/bin/RokytaScripts/AutoAnnotatorFormat.py -i Nextannotated_toxins.fasta -f seqman2
cd ../../..
done
```

### Annotation in Geneious

move the sequences to a Desktop with Geneious

make the directory structure suggested by Rokyta for Geneious and import
the sequences

Geneious/ \| Cgodm-CLPxxxx/ \| - \| Toxins/ \| - \| Nontoxins/ \| - \|
Unannotated/ \| - - \| Toxins/ \| - - \| Nontoxins/ \| - - \| z\_Trash/
\| - \| Combined \| - - \| CDS

for the toxins in unannotated folder

make geneious to show all the CDS, use the xml files from blast to check
if any of the CDS match with a toxin following the length and the
percentace of similarity if there is a match annotate the toxin family
and the name of the contig (TOXIN\_<Toxin Family>\_<contig name>) in the
CDS then move the sequence to the annotated folder if there is not a
match move the sequence to the trash folder

to make this esier trim the xml files with the next code

``` bash
#PBS -N autoannotate
#PBS -l select=13:ncpus=16:mem=60gb:interconnect=fdr,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

source .bash_profile

module load anaconda3/5.1.0-gcc/8.3.1
source activate bio

cd $PBS_O_WORKDIR

parallel -a list1 --sshloginfile $PBS_NODEFILE -j1 "source ~/.bash_profile
module load anaconda3/5.1.0-gcc/8.3.1
source activate bio
cd /zfs/venom/Ramses/Cerrophidion
cd {}/08blast
XML_Barber_v3.py -f ../10nextannotate/NextAnnotator/Still_unannotated_toxins.fasta -x {}_blastUPT.xml -o {}_blastUPT_reduced.xml
grep '>' ../10nextannotate/NextAnnotator/Still_unannotated_toxins.fasta > Still_unannotated_list
perl -p -i -e 's/\>//g' Still_unannotated_list
    for i in `cat Still_unannotated_list`
    do echo $i
    mkdir Bv
    python /zfs/venom/Ramses/bin/RokytaScripts/BV.py ${i}_blastUPT_reduced2.xml ${i} > Bv/${i}_Bv 
    done
cd ../.."
```

### move back to palmetto

export the CDS to a directory named 11finalannotation

``` bash
#continuation afterr finish work in geneious and export the files as name_CDS.toxin in the directory 11finalannotation

#create the folder with this code
for i in `cat list1`
do echo $i
cd $i
mkdir 11finalannotation
cd ..
done

#then I can go back to cerrophidion and  run the next script
conda activate bio
#needs bioconda
#wasn't working, I had to comment the next line
##from sets import Set

parallel -a list1 -j 2 --verbose "cd {}/11finalannotation
/Users/ramsesrosales/Documents/bin/AndrewScripts/RemDupRemAmb.py -f {}_CDS.fasta -o {}"
```

then move it back to palmetto with a filezilla software.

# Post Annotation Filters

## Concatenate and Clean

``` bash
#this is the previous corrections on the data before start with Chimera Killer

#to revise the number of toxins 
for i in `cat list1`
do echo $i
grep -c '>' $i/12venomancer/${i}_Toxins_cds_RedundancyFiltered.fasta
grep -c '>' $i/12venomancer/${i}_Toxins_cds_SPfiltered_RedundancyFiltered.fasta
grep -c '>' $i/12venomancer/${i}_PutativeToxins_cds_SPfiltered.fasta
grep -c '>' $i/12venomancer/${i}_Toxins.fasta
done

#seen the number of contigs, the {}_Toxins.fasta is not concatenate correctly, I have to concatenate again with a for loop

interactive
Cerrophidion

for i in `cat list1`
do echo $i 
cd $i/12venomancer
cat ${i}_PutativeToxins_cds_SPfiltered.fasta ${i}_Toxins_cds_RedundancyFiltered.fasta > ${i}_Toxins.fasta
perl -p -i -e 's/>/>TOXIN_/g' ${i}_Toxins.fasta
cd ../..
done


#convine toxins and nontoxins of ToxCodAn, Venomancer

for i in `cat list1`
do echo $i 
mkdir $i/14completeannotation
cd $i/14completeannotation
cat ../12venomancer/${i}_Toxins.fasta ../13nontoxins/${i}_NonToxins_annotated.fasta > ${i}_annotated.fasta
cd ../..
done
```

``` bash
#use parallel and CDHIT to concatenate manual and venomancer annotation and remove redundancy
#PBS -N merge_annotation
#PBS -l select=13:ncpus=20:interconnect=fdr:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

source .bash_profile

module load anaconda3/5.1.0-gcc/8.3.1
source activate bio

cd $PBS_O_WORKDIR

parallel -a list1 --sshloginfile $PBS_NODEFILE -j1 " source .bash_profile
module load anaconda3/5.1.0-gcc/8.3.1
source activate bio
cd /zfs/venom/Ramses/Cerrophidion
cd {}/14completeannotation
cat {}_annotated.fasta ../11finalannotation/{}_clean.fasta > {}_ma_toxcodan.fasta
cd-hit-est -i {}_ma_toxcodan.fasta -o {}_annotated_reduced.fasta -d 0 -c 1.00 -M 95000
cd ../.." 
```

## Chimera Killer

modification to chimera killer code

``` bash
##########picard modified with new syntaxis
https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition)

command = bwa + " mem -M -t " + str(procs) + " -R \'@RG\\tID:" + input.name + "\\tSM:" + reads.name + "' " + input.name + " " + reads.name + " | " + grepNM  + " > tmp1.sam"
subprocess.call(command,shell=True)
# Create a sorted bam file
command = picard + " SortSam -INPUT tmp1.sam -OUTPUT tmp2.bam -SORT_ORDER coordinate -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
subprocess.call(command,shell=True)
command = picard + " BuildBamIndex -INPUT tmp2.bam -USE_JDK_DEFLATER true -USE_JDK_INFLATER=true"
subprocess.call(command,shell=True)
# Remove overclipped reads
command = picard + " CreateSequenceDictionary -REFERENCE " + input.name + " -OUTPUT " + input.name.split(".")[0] + ".dict -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
subprocess.call(command,shell=True)
command = samtools + " faidx " + input.name
subprocess.call(command,shell=True)
```

``` bash
#PBS -N ChimeraKiller
#PBS -l select=13:ncpus=20:interconnect=fdr:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

source .bash_profile

module load anaconda3/5.1.0-gcc/8.3.1
source activate bio

cd $PBS_O_WORKDIR

parallel -a list1 --sshloginfile $PBS_NODEFILE -j1 " source .bash_profile
module load anaconda3/5.1.0-gcc/8.3.1
source activate chimerakiller_env 
cd /zfs/venom/Ramses/Cerrophidion
mkdir {}/15chimerakiller
cd {}/15chimerakiller
cp ../14completeannotation/{}_annotated_reduced.fasta .
ChimeraKiller_v0.7.3.py -i {}_annotated_reduced.fasta  -r ../03merged/{}_pear.assembled.fastq -d 0.75 -p 18"
```

Notes: use the full path for ChimeraKiller input

some characters in the names can produce errors, as “=”

move the 15chimerakiller to Desktop, check the coverage plots of the
toxins and move any souspicious chimera from good to bad folder and any
contig that look god from bad to good

# Consense Transcriptomes

## Cluster contigs

``` bash
interactive #alias to start a interactive job
Cerrophidion #alias to go to directory with the samples

for i in `cat list1`
do echo $i
mkdir $i/16transcriptome
cd $i/15chimerakiller/fastas/good
cat *.fasta >> /Users/ramsesrosales/Desktop/Cerrophidion/${i}/16transcriptome/${i}_transcriptome_v0.fasta
cd /Users/ramsesrosales/Desktop/Cerrophidion/${i}/16transcriptome/
cd-hit-est -i ${i}_transcriptome_v0.fasta -o ${i}_transcriptome_v1.fasta -d 0 -c 0.99 -M 15000
done
```

## Filtering Transcriptomes (CDS\_filter and RSEM\_combiner)

``` bash
ls -d *Cgod* > ID_Cgodm
ls -d *Ctzot* > ID_Ctzot
ls -d *Cpetl* > ID_Cpetl
#ls -d *Csasa* > ID_sasa do this at the start is more efficient, I didn't

ls ID_* >species
perl -p -i -e 's/ID_//g' Species

for i in `cat Species`
do mkdir -p ${i}_transcriptome/01transcriptome
cat $i-*/16transcriptome1/*transcriptome_v1.fasta >> ${i}_transcriptome/01transcriptome/${i}_consensus_v0.fasta
cd ${i}_transcriptome/01transcriptome
cd-hit-est -i ${i}_consensus_v0.fasta -o ${i}_consensus_v1.fasta -d 0 -c 0.98 -M 15000
mkdir ../02rsem
cd ../02rsem
rsem-prepare-reference --bowtie2 ../01transcriptome/${i}_consensus_v1.fasta ${i}_Reference
cd ../..
done

for i in `cat Species`
do echo $i
cd ${i}_transcriptome/02rsem
for j in `cat ../../ID_${i}`
do echo $j
rsem-calculate-expression --bowtie2 ../../$j/03merged/${j}_pear.assembled.fastq ${i}_Reference $j
rsem-plot-model $j ${j}_diagnostic.pdf
done
cd ../..
done 

#for Csasai

mkdir Csasai_transcriptome
cd Csasai_transcriptome
mkdir 01transcriptome
cp ../Csasa-LIAP244/16transcriptome1/*transcriptome_v1.fasta  /01transcriptome/
mkdir 02rsem
cd 02rsem
rsem-prepare-reference --bowtie2 ../16transcriptome/*_consensus_v1.fasta Csasa-LIAP244_Reference
rsem-calculate-expression --bowtie2 ../../Csasa-LIAP244/03merged/Csasa-LIAP244_pear.assembled.fastq Csasa_Reference Csasa-LIAP244
rsem-plot-model Csasa-LIAP244 Csasa-LIAP244_diagnostic.pdf
```

moved consensus transcriptomes to a directory called
final\_transcriptomes then used script CDS\_filter\_v0.1.py (ask
mason.501@osu.edu) to remove CDS with internal stop codons without stop
codon or sequences of incorrect length(not divisible by 3)

``` bash
echo Csasa >> species
echo Csasa-LIAP244 >> ID_Csasa
cd final_transcriptomes
for i in `cat ../Species`
mkdir other_filter
cd other_filter
mkdir $i
cd $i
cp ../../${i}*_v1.fasta .
CDS_filter_v0.1.py -f ${i}*_v1.fasta -o $i -hs T
rm ${i}**_consensus*_v1.fasta
cd ../..
```

Then used the script RSEM\_combiner.R
(<https://github.com/RhettRautsaw/Bioinformatics/blob/master/scripts/RSEM_combiner.R>)
to filter the Toxins with expression below the 95% percentile

``` bash
for i in `cat Species` 
cd ${i}_transcriptome/02rsem
RSEM_combiner.R -g TOXIN -c 0.95 -o $i
cd ../../
```

## Clean consense

We then made a second consensus trancriptome using the output of both
filters above

``` bash
qsub Rsem_cleaned_2
```

``` bash
#PBS -N Rsem_cleaned_2
#PBS -l select=13:ncpus=20:interconnect=fdr:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

source .bashrc
module load parallel/20190222-gcc/8.3.1
cd $PBS_O_WORKDIR

parallel -a Species --sshloginfile $PBS_NODEFILE -j1 "echo {}
source .bashrc
module load anaconda3/2021.05-gcc/8.3.1
source activate bio
cd /zfs/venom/Ramses/Cerrophidion
cd {}_transcriptome
pwd
mkdir 02rsem_clean2
cp /zfs/venom/Ramses/Cerrophidion/final_transcriptomes/other_filter/{}/{}_cds_filtered.fasta 02rsem_clean2
cd 02rsem_clean2
cat ../02rsem/*.list > remove.list
RemoveSeqs.py -f {}_cds_filtered.fasta -l remove.list -o {}_consensus_V2.fasta
rsem-prepare-reference --bowtie2 {}_consensus_V2.fasta {}_Reference
cd ../.."
```

## RSEM, Plots, and Differential expression

We run RSEM with the cleaned consensus transcriptome as referece, then
we used the results to make plots and run DESeq2 and edgeR programs.

``` bash
qsub 
```

``` bash
#PBS -N Dif_Exp_2
#PBS -l select=4:ncpus=20:interconnect=fdr:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

source .bashrc
module load parallel/20190222-gcc/8.3.1
cd $PBS_O_WORKDIR

parallel -a Species  --sshloginfile $PBS_NODEFILE -j1 " source .bashrc
module load r/4.0.2-gcc/8.3.1
cd /zfs/venom/Ramses/Cerrophidion
cd {}_transcriptome/02rsem_clean2
cp ../../ID_{} .
data_frame_consense.R -l ID_{}
perl -p -i -e 's/Cgodm-CLP2903/Ctzot-CLP2903/g' ID_{}_rsem_Cerrophidion.csv
perl -p -i -e 's/Cgodm-CLP2903/Ctzot-CLP2903/g' ID_{}
/zfs/venom/Ramses/bin/RamsesScripts/Rscripts/consense.R -i ID_{}_rsem_Cerrophidion.csv -s {}
/zfs/venom/Ramses/bin/RamsesScripts/Rscripts/consense_expected_count.R -i ID_{}_rsem_Cerrophidion.csv -s {}
/zfs/venom/Ramses/bin/RamsesScripts/Rscripts/ToxNames.R -i {}_consense_df.csv -m {}_consense_ExpCount_df.csv -s {}
/zfs/venom/Ramses/bin/RamsesScripts/Rscripts/FancyPlots3.R -i ID_{}_rsem_Cerrophidion.csv -n ID_{}
/zfs/venom/Ramses/bin/RamsesScripts/Rscripts/FancyPlotsAverage.R -i {}_consense_df.csv -n {}_Average
mkdir Plots
mv *transcriptome.pdf Plots
cd ../../"


parallel -a pop_s_list  --sshloginfile $PBS_NODEFILE -j1 " source .bashrc
module load r/4.0.2-gcc/8.3.1
source activate Dif_Exp_env
cd /zfs/venom/Ramses/Cerrophidion
cd {}_transcriptome/02rsem_clean2
cp ../02rsem/NS_list .
/zfs/venom/Ramses/bin/RamsesScripts/Rscripts/FancyPlotAverageNS.R -i {}_consense_df.csv -l NS_list -n {}
cd ../../"


parallel -a pop_s_list  --sshloginfile $PBS_NODEFILE -j1 " source .bashrc
module load anaconda3/2021.05-gcc/8.3.1
source activate Dif_Exp_env
cd /zfs/venom/Ramses/Cerrophidion
cd {}_transcriptome/02rsem_clean2
/zfs/venom/Ramses/bin/RamsesScripts/Rscripts/Dif_Exp_Cerrophidion.R -i {}_consense_ExpCount_df_original.csv -m  ../../Cerrophidion_specimens1.csv -s {} 
conda deactivate Dif_Exp_env
module purge
module load r/4.0.2-gcc/8.3.1
cd Dif_Exp
/zfs/venom/Ramses/bin/RamsesScripts/Rscripts/Dif_Exp_Tab.R -i ../{}_consense_ExpCount_df_original.csv -s {} -p T
/zfs/venom/Ramses/bin/RamsesScripts/Rscripts/Cerrophidion_HeatMap_AverageExpressionOrder.R -i ../{}_consense_df.csv -m ../../../Cerrophidion_specimens1.csv -d {}_DifExp_tab.csv -s {} -p T  -t {}_DifExp_TF_tab.csv -c T -w 1500 -e 1800
cd ../../.."


module load anaconda3/2021.05-gcc/8.3.1
source activate Dif_Exp_env
cd /zfs/venom/Ramses/Cerrophidion
cd Cpetl_transcriptome/02rsem_clean2
Dif_Exp_Cerrophidion.R -i Cpetl_consense_ExpCount_df_original.csv -m  ../../Cerrophidion_specimens1.csv -s Cpetl -t SVL
conda deactivate Dif_Exp_env
module purge
module load r/4.0.2-gcc/8.3.1
cd Dif_Exp
Dif_Exp_Tab.R -i ../Cpetl_consense_ExpCount_df_original.csv -s Cpetl -p F
Cerrophidion_HeatMap_AverageExpressionOrder.R -i ../Cpetl_consense_df.csv -m ../../../Cerrophidion_specimens1.csv -d Cpetl_DifExp_tab.csv -s Cgodm -p F  -t Cpetl_DifExp_TF_tab.csv -c T -w 1200 -e 1500
cd ../../..
```

heatmap for figure in the paper

``` bash
#scritp that works for Cgodm only, it can be modified in the script changing rn<-read_csv("../Cgodm_ToxNames.csv") by the path to the file with the ordered toxin names.
cd Cgodm_transcriptome/02rsem_clean2
cd Dif_Exp
Cerrophidion_HeatMap_order_difExp.R -i ../Cgodm_consense_df.csv -m ../../Cerrophidion_specimens1.csv -d Cgodm_DifExp_tab.csv -s Cgodm -p T  -t Cgodm_DifExp_TF_tab.csv -c T -w 1500 -e 1800
cd ../../..
```

stacked plots for the figures, download the <species>\_consense\_df.csv
files to local computer in the same folder

``` bash
stacked_toxin_plot.R -f F -i Cgodm_consense_df.csv -o Cgodm -n 6 -t tf_list
# tf_list is a list with the names of the genes that were identify as crotoxin subunit homologs and the new name of the toxin family 
stacked_toxin_plot.R -f F -i Ctzot>_consense_df.csv -o Ctzot -n 4
stacked_toxin_plot.R -f F -i Cpetl_consense_df.csv -o Cpetl -n 3
```

content of tf\_list

``` bash
TOXIN_extenderContig120||Cgodm-CLP2360trinity_PLA2-3a_trinityContig14327_PLA2,PLA2_neuro
TOXIN_extenderContig449||Cgodm-CLP2360trinity_PLA2-2a_trinityContig14324_PLA2,PLA2_neuro
```

# PLA2s Phylogeny

## Get PLA2s secuences

We extracted all the PLA2s sequences from the cleaned consensus
transcriptomes the shell script we used also rename the fasta using the
list we got from the script ToxNames.R which makes a list of the toxins
ordered by expression.

``` bash
mkdir PLA2s
./extractTF.sh PLA2 /zfs/venom/Ramses/Cerrophidion/PLA2s
```

code in ./extractTF.sh this script is specific for my folder
organization

``` bash
#!/bin/bash

###this program expect two arguments
###the first one is the name of the toxin family
### the second one is the output path respect to the Cerrophidion folder
### it uses a list with the abbreviation of the species

A=$(pwd)

for i in `cat Species`
do echo $i
cd ${i}_transcriptome/02rsem_clean2
grep -A1 -e $1 ${i}_consensus_V2.fasta | grep -A1 'TOXIN' > ${i}_${1}.fasta
grep -e $1 ${i}_ToxNames.csv > ${i}_${1}_ref.list
perl -p -i -e "s/TOXIN/${i}_TOXINS/g" ${i}_${1}.fasta
perl -p -i -e "s/TOXIN/${i}_TOXINS/g" ${i}_${1}_ref.list
cd ../../
mv ${i}_transcriptome/02rsem_clean2/${i}_${1}.fasta $2/${i}_${1}.fasta
mv ${i}_transcriptome/02rsem_clean2/${i}_${1}_ref.list $2/${i}_${1}_ref.list
module load r/4.0.2-gcc/8.3.1
cd $A/$2
echo old_name,new_name > tmp
cat ${i}_${1}_ref.list >> tmp
mv tmp ${i}_${1}_ref.list 
Rename_Fasta_Headers.R -i ${i}_${1}.fasta -t ${i}_${1}_ref.list
cat ${i}_${1}_renamed.fasta >> Cerrophidion_${1}.fasta
cd $A
done
```

We used a script to download some samples from genebank and extract the
CDS, we had a list of the accession numbers. Then we concatenated the
fasta files and used a script to clean the names and make a new fasta
with the simplyfied names.

``` bash
#need bio environment
for i in `cat list`
do echo $i 
import_gb.py -a $i -d nucleotide -e <your email addres> -f gb
extract_CDS.py -i ${i}.gb -f gb -o $i
done

cat *.fasta >> other_PLA2.fasta
./clean_names.sh other_PLA.fata
Rename_Fasta_Headers.R -i other_PLA2.fasta -t reference_names.csv -o T
```

We got Bothriechis nigroviridis and B. nubestris and other Agkistrodon,
Crotalus and Sistrurus samples from mason.501@osu.edu We got Mixcoatlus
melanurum sequence from

then we moved all the sequences in one folder and concatenated them

``` bash
cat Bothriechis.fasta Cerrophidion_PLA2.fasta melanurotoxin.fasta Rooted_Crotalus_PLA2_Alignment_2017.fasta >> all.fasta
./change_names.sh all.fasta
```

align the sequences with mafft, then clean and trim alignment with
CIAlign and trimal finally add the outgroup again and aling again with
mafft

``` bash
./run_tree.sh
```

make sure the name of the outgroup is “P.bivitatus”.

``` bash
#!/bin/bash

module load anaconda3/2021.05-gcc/8.3.1
source activate bio

mafft --auto --adjustdirectionaccurately --thread 2 all.fasta > all.aln.fasta
CIAlign --infile all.aln.fasta --outfile_stem all --remove_divergent --remove_divergent_minperc 0.80 --remove_insertions --crop_ends --remove_short
trimal -in all_cleaned.fasta -out all.trimal.fasta -strictplus
grep -A1 "P._b" all.fasta >> all.trimal.fasta
mafft --auto --adjustdirectionaccurately --thread 2 all.trimal.fasta > all.trim.aln.fasta


iqtree -s all.trim.aln.fasta -B 1000 -seed 12345 -T 8
```

``` bash
```

Edit names on the tree

``` bash
cp all.trim.aln.fasta.contree Figure.tree
Change_tree_names.R -i Figure.tree -g T
#edit the names in excel <tree_names.csv>
Change_tree_names.R -i Figure.tree -n tree_names_new.csv -o Figure -r T 
```

the final tree was these

``` bash
(B._nigroviridis_2a:2.1331e-06,((((((((((((((((((((B._nigroviridis_3a:2.1331e-06,B._nigroviridis_17:2.1331e-06)100:0.0028336081,B._nubestris_3a:2.1331e-06)100:0.0397104179,(((C._scutulatus_2:0.0028361065,C._tigris_2:2.1331e-06)50:2.1331e-06,C._durissus_CrotA:0.0085633854)99:0.0121405855,(C._horridus_II_2:0.0057673886,S._tergeminus_A:0.0115421807)97:0.0022821263)100:0.0530829202)95:0.0175033583,C._godmani_11:0.0346068164)85:0.0199898999,(M._melanurum_A:0.0833360079,G._intermedius_A:0.0467984905)87:0.0158679557)85:0.0134285049,((C._godmani_4:0.0257378673,(B._pictus:0.0413843607,((B._erythromelas:2.1331e-06,(B._diporus_1:0.0028591714,B._diporus_2:0.0143865518)71:0.0028354662)49:2.1367e-06,B._neuwiedi:0.0294824764)97:0.0100349163)76:0.0182346512)61:0.0155862628,(G._halys:0.0286405145,((G._shedaoensis:0.0158168937,G._intermedius_1:0.0116856732)63:2.1331e-06,G._strauchi:0.0054914921)53:0.0031325371)97:0.0104269535)59:0.014703096)100:0.0411309865,(((((((C._adamanteus_1:0.0058348734,(C._atrox_1:0.0027698473,C._scutulatus_1:0.0057846184)100:0.0174924151)97:0.0029639528,C._horridus_II_3:0.0057829774)100:0.0390898209,C._mitchellii_4:0.0285198891)100:0.0117559886,S._miliarius_2:0.0515850585)92:0.0154121905,((C._cerastes_1:0.033071398,C._horridus_I_1:0.0269609382)73:0.0024164802,(C._lepidusi_1:2.1331e-06,C._lepidusi_3:0.023381102)98:0.0086618174)74:0.0054147314)99:0.0112748903,S._catenatus_1:0.0415389169)95:0.0158263426,((A._contortrix_3:0.0322779646,A._piscivorus_2:0.0034829143)100:0.0102368098,(A._contortrix_1:0.0372020297,A._piscivorus_3:0.0333423573)100:0.0425037587)99:0.0170556215)95:0.0164187968)66:0.0203272312,C._petlalcalensis_50:0.0792680038)62:0.0073659197,(C._petlalcalensis_1:2.1331e-06,C._tzotzilorum_1:0.0115062216)53:2.0941e-06)60:0.0159710796,C._sasai_3:0.0292490941)60:0.0109803158,C._godmani_1:2.4337e-06)100:0.092397015,C._godmani_43:2.1435e-06)99:0.084744414,(((C._godmani_35:0.0013568707,C._sasai_4:0.0130613107)100:0.0728656023,C._godmani_110:0.000945821)83:0.010554649,P._bivitatus:1.364220179)53:0.0299781524)51:0.0225711065,(((C._godmani_3:2.1331e-06,C._sasai_2:0.0028451758)56:2.1331e-06,(((((((C._godmani_41:0.0025483299,C._godmani_26:0.0730297871)81:0.0060869284,C._godmani_42:0.0397303178)79:0.0222570506,(C._tzotzilorum_3:0.0073131334,C._tzotzilorum_40:0.003252526)100:0.0194447685)100:0.0430175706,(A._contortrix_2:0.0171190593,A._piscivorus_1:0.0222262294)100:0.0715169582)78:0.0193077249,(C._atrox_3:0.0282813154,C._mitchellii_1:0.0072109254)100:0.066123409)97:0.008690011,((C._godmani_8:0.0186852628,C._tzotzilorum_35:0.0044561468)98:0.0085284479,C._tzotzilorum_39:0.0100419913)98:0.0091038658)98:0.014909511,C._tzotzilorum_73:0.0071204548)75:0.0014215294)100:0.0858799409,C._godmani_57:0.0012470277)98:0.0894698323)100:0.0516495152,C._godmani_18:0.0176376126)66:0.0027228403,(((C._godmani_27:0.0484833511,C._tzotzilorum_58:0.0305695055)96:0.0241121701,C._petlalcalensis_49:0.0651210956)69:0.0162691181,(C._petlalcalensis_9:0.0028761464,C._tzotzilorum_19:0.0027980651)99:0.0059979352)66:0.0027304518)56:0.0057659649,((((C._atrox_2:0.0057420239,C._mitchellii_3:0.005722935)97:0.0111592531,(S._miliarius_1:0.0148178743,S._catenatus_2:0.032927937)94:0.0057708752)86:0.0029668156,C._mitchellii_2:0.0262004993)88:0.0025945327,C._lepidusi_2:0.0152081766)91:0.0268045205)100:0.0720355745,((((C._tigris_1:0.0057788265,(((C._scutulatus_3:0.00571258,C._durissus_CrotB:0.0264083051)99:0.0028676106,C._basiliscus:2.1331e-06)84:2.1331e-06,C._oreganus_B:0.0085763874)62:0.0027807687)95:0.0137857377,C._horridus_II_1:0.009624571)54:0.0103009722,S._tergeminus_B:0.017184214)94:0.0218541616,G._intermedius_B:0.0435920278)100:0.0364616416)98:0.0142307934,C._godmani_16:0.0299343427)100:0.0225865958,B._nubestris_2a:0.0013128332)100:0.0072349257,B._nigroviridis:2.1331e-06);
```

Get Isoelectric Point

``` bash
Translate_RARG.py -f all.fasta -n all.prot
# I edited the file manually changing the CDS of the sequences from genbank files
IsoelectricPoint_df.py -f all.prot.fasta -o all_IP
```

Further edition of the tree and the final figure was done in fig.tree
and Inkscape

# Sequence Diversity 1

## Busco and Phylogenetics

### Busco

busco 5 installed in a conda environment and run with a pbs file

``` bash
qsub Busco2.pbs
```

``` bash
#PBS -N Busco_locus
#PBS -l select=1:ncpus=20:interconnect=fdr:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

source .bash_profile

module load anaconda3/2020.07-gcc/8.3.1
source activate Busco_env

cd $PBS_O_WORKDIR

parallel -a ../list1 -j13 "cd /zfs/venom/Ramses/Cerrophidion/Busco/
busco -i /zfs/venom/Ramses/Cerrophidion/{}/07assembly/{}_assembly_reduced.fasta -o {}_Busco -m transcriptome -l /zfs/venom/Ramses/bin/Venomancer/non_toxin_models/tetrapoda_odb10"
```

run BuscoCleaner.py available at
<https://github.com/RhettRautsaw/Bioinformatics/tree/master/scripts>

``` bash
#pbs script for this
#PBS -N Busco_cleaner
#PBS -l select=1:ncpus=20:interconnect=fdr:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

module load anaconda3/2020.07-gcc/8.3.1
source activate Busco_cleaner_env

cd $PBS_O_WORKDIR
cd phylogenomics
for i in `cat ../../list1`
do BuscoCleaner.py -f $i -n $I
done 
```

### Prepare sequences

merge busco sequences for all individuals

``` bash
./tmp.list_phylogenomics.sh
```

code of tmp.list\_phylogenomics.sh

``` bash
#!/bin/bash

for i in $(ls | grep 'run_C')
do echo $i
cd $i/BuscoCleaner
ls >> /zfs/venom/Ramses/Cerrophidion/Busco/phylogenomics/tmp.list
cd /zfs/venom/Ramses/Cerrophidion/Busco
done


cd phylogenomics
perl -p -i -e "s/.fasta//g" tmp.list
sort tmp.list | uniq > tmp.list1
sort tmp.list | uniq -c > tmp.list_counts

mkdir BuscoSeq

for i in $(cat tmp.list1)
do echo $i
cat ../run*/BuscoCleaner/${i}.fasta >> BuscoSeq/${i}_Cerrophidion.fasta
done
```

``` bash
mkdir all_genes
mv BuscoSeq all_genes
```

### All genes

#### Alignment, Cleaning and Triminig

``` bash
#PBS -N GeneTrees
#PBS -l select=1:ncpus=20:interconnect=fdr:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

cd /zfs/venom/Ramses/Cerrophidion/Busco/phylogenomics/all_genes

./clean.sh
```

``` bash
#!/bin/bash

mkdir 01_aln
mkdir 02_CIAlign
mkdir 03_treeshrink
mkdir 04_CIAlign2
mkdir 05_trimal
mkdir 06_genetrees
mkdir 07_speciestrees

module load anaconda3/2020.07-gcc/8.3.1
source activate bio
parallel -a ../tmp.list1 -j 20 "echo {}
       mafft --auto --adjustdirectionaccurately --thread 2 BuscoSeq/{}_Cerrophidion.fasta > 01_aln/{}.aln.fasta
       CIAlign --infile 01_aln/{}.aln.fasta --outfile_stem 02_CIAlign/{} --remove_divergent --remove_divergent_minperc 0.80 --remove_insertions --crop_ends --remove_short
       cd 02_CIAlign
       iqtree -s {}_cleaned.fasta -bb 1000 -seed 12345 -T 2
       cd ../03_treeshrink
       mkdir {}
       cp ../02_CIAlign/{}_cleaned.fasta {}/input.fasta
       cp ../02_CIAlign/{}_cleaned.fasta.contree {}/input.tree"

conda deactivate


source activate treeshrink_env
run_treeshrink.py -i 03_treeshrink -t input.tree -a input.fasta
conda deactivate

source activate bio
parallel -a ../tmp.list1 -j 20 --verbose "echo {}
       CIAlign --infile 03_treeshrink/{}/output.fasta --outfile_stem 04_CIAlign2/{} --remove_divergent --remove_insertions --crop_ends --remove_short
       trimal -in 04_CIAlign2/{}_cleaned.fasta -out 05_trimal/{}.fasta -strictplus"
```

#### Gene trees and Species trees

``` bash
#PBS -N Species_tree
#PBS -l select=1:ncpus=20:interconnect=fdr:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

cd /zfs/venom/Ramses/Cerrophidion/Busco/phylogenomics/all_genes

./gene_tree.sh
```

``` bash
#!/bin/bash

module load anaconda3/2020.07-gcc/8.3.1
source activate bio

cd 05_trimal
ls *.fasta > ../gene_list
cd ..

perl -p -i -e "s/\.fasta//g" gene_list

cd 06_genetrees

parallel -a ../gene_list -j 20 --verbose "echo {}
cp ../05_trimal/{}.fasta .
iqtree -s {}.fasta -bb 1000 -seed 12345 -T 2
rm {}.fasta"

cd ../07_speciestrees
cat ../06_genetrees/*.treefile >> Cerrophidion_Busco_genes.tre

for i in `cat ../gene_list`
do echo $i
perl -pi -e "s/${i}\|//g" Cerrophidion_Busco_genes.tre
done

java -jar /zfs/venom/Ramses/bin/ASTRAL/Astral/astral.5.7.7.jar -i Cerrophidion_Busco_genes.tre -o Cerrophidion_astral.tree
```

### Core-Genes

#### Subset core-genes

``` bash
grep "14 " tmp.list_counts | cut -d " " -f 7 > core_gene_list
```

``` bash
./subset_trimal.sh core_genes core_gene_list
```

``` bash
#!/bin/bash
#first argument name of the new directory
#second argument list of genes to subset from all_genes/trimal

mkdir $1

cd $1
mkdir 05_trimal
mkdir 06_genetrees
mkdir 07_speciestrees

for i in $(cat ../$2)
do echo $i
cp ../all_genes/05_trimal/${i}.fasta 05_trimal
done
```

#### Concatenated sequences

``` bash
./concatenate.sh
```

the Concatenate4PartitionFinder.py was obtained from mason.501@osu.edu

``` bash
#/bin/bash

cd 08_concatenatedtree
cp ../05_trimal/*  .

for i in `cat ../gene_list`
do echo $i
perl -pi -e "s/${i}\|//g" ${i}.fasta
done

cd ..

Concatenate4PartitionFinder.py -d 08_concatenatedtree -o Cerrophidion_core_genes_concat.phy
```

### Scaling tree with IQtree

``` bash
cd /08_concatenatedtree/Concatenated/
mkdir scaled_tree
```

edit file partitions.txt for IQtree, to to this add the next lines in
the begining:

``` bash
#nexus
begin sets;'
```

then add the next before each gene name

``` bash
charset <genename> = <start>-<finish>;
...
...
charset <genename> = <start>-<finish>;
```

you can add these part with one script

``` bash
perl -pi -e 's/^(.*)/^charset $1/g' partitions.txt
```

at the end add:

``` bash
end;
```

Use IQtree to make a tree, for these use the species tree to fix the
topology and the concatenated sequences and the partitions to scale the
lenght of the branches

``` bash
cd scaled_tree
cp ../../../../all_genes/07_speciestrees/Cerrophidion_astral.tree .
cp ../Cerro_core_alignment.phy .

iqtree -s Cerro_core_alignment.phy -g Cerrophidion_astral.tree -p ../partitions.txt -m MFP -B 1000
```

# Sequence Diversity 2

## Variant Calling

``` bash
Cerrophidion ## alias to go to the starting directory
mkdir Variants
cd Variants
```

``` bash
qsub Cgodm_VariantCalling2.pbs
```

``` bash
#PBS -N Variant_Cgodm_2
#PBS -l select=1:ncpus=20:interconnect=any:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

source .bash_profile

module load anaconda3/5.1.0-gcc/8.3.1
source activate VariantsCalling_env

cd $PBS_O_WORKDIR

mkdir Cgodm
cp /zfs/venom/Ramses/Cerrophidion/Cgodm_transcriptome/02rsem_clean2/Cgodm_consensus_V2.fasta Cgodm/ref.fa
cd Cgodm
mkdir RESULTS

bwa index ref.fa 

picard CreateSequenceDictionary -REFERENCE ref.fa -OUTPUT ref.dict
samtools faidx ref.fa
#sudo pip install pyfaidx
faidx --transform bed ref.fa > ref.bed


parallel -a ../ID_Cgodm -j20 "source .bash_profile
    module load anaconda3/5.1.0-gcc/8.3.1
    source activate VariantsCalling_env
    cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm
    mkdir {}
    cd {}
    bwa mem -M -t 2 -R '@RG\tID:Cgodm\tSM:{}' ../ref.fa /zfs/venom/Ramses/Cerrophidion/{}/03merged/{}_pear.assembled.fastq > {}_aln.sam
    cat {}_aln.sam | grep -E 'NM:i:[0-2][[:space:]]|^@' > {}_filtered.sam
    picard SortSam -INPUT {}_filtered.sam -OUTPUT {}.bam -SORT_ORDER coordinate
    picard BuildBamIndex -INPUT {}.bam
    rm {}*.sam
    java -jar /zfs/venom/Ramses/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../ref.fa -I {}.bam -o realigner.intervals
    java -jar /zfs/venom/Ramses/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T IndelRealigner -R ../ref.fa -I {}.bam -targetIntervals realigner.intervals  -o {}_realigned.bam
    java -jar /zfs/venom/Ramses/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T PrintReads -R ../ref.fa -I {}_realigned.bam  -rf OverclippedRead --filter_is_too_short_value 120 --do_not_require_softclips_both_ends -rf MappingQuality -mmq 40 -rf ReadLength -minRead 120 -maxRead 500 -o {}_realigned_filter.bam
    gatk HaplotypeCaller -R ../ref.fa -I {}_realigned_filter.bam -O {}.g.vcf -ERC GVCF
    gatk GenotypeGVCFs -R ../ref.fa -V {}.g.vcf -O output.vcf
    gatk SelectVariants -R ../ref.fa -V output.vcf -select-type SNP -O raw_snps.vcf
    gatk VariantFiltration -R ../ref.fa -V raw_snps.vcf --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name 'FILTER' -O filtered_snps.vcf
    gatk SelectVariants -R ../ref.fa -V output.vcf -select-type INDEL -O raw_indels.vcf
    gatk VariantFiltration -R ../ref.fa -V raw_indels.vcf --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' --filter-name 'FILTER' -O filtered_indels.vcf
    grep -E '^#|PASS' filtered_snps.vcf > {}_SNPs.vcf
    grep -E '^#|PASS' filtered_indels.vcf > {}_INDELs.vcf
    bedtools coverage -a ../ref.bed -b {}.bam -d > {}_coverage.txt
    bedtools genomecov -bga -ibam {}.bam -g ../ref.bed > tmp.bed
    grep -w 0$ tmp.bed > 0cov.bed
    rm tmp.bed
    cp {}_SNPs.vcf ../RESULTS/
    cp {}_INDELs.vcf ../RESULTS/
    cp {}_coverage.txt ../RESULTS/
    cp 0cov.bed ../RESULTS/{}_0cov.bed
    whatshap phase --reference ../ref.fa -o {}_phased.vcf {}_SNPs.vcf {}.bam
    cp {}_phased.vcf ../RESULTS/
    bgzip {}_phased.vcf
    tabix {}_phased.vcf.gz
    bcftools consensus -H 1 -f ../ref.fa {}_phased.vcf.gz > haplotype1.fasta
    bcftools consensus -H 2 -f ../ref.fa {}_phased.vcf.gz > haplotype2.fasta
    bedtools maskfasta -fi haplotype1.fasta -bed 0cov.bed -fo {}_allele1.fasta -mc -
    bedtools maskfasta -fi haplotype2.fasta -bed 0cov.bed -fo {}_allele2.fasta -mc -
    rm haplotype1.fasta haplotype2.fasta
    cp {}_allele*.fasta ../RESULTS/
    cp {}_phased.vcf.gz.tbi ../RESULTS/
    cd .."

cd RESULTS/
mkdir SNPs INDELs PhasedAlleles Coverage PhasedSNPs
mv *_SNPs.vcf SNPs/
mv *_INDELs.vcf INDELs/
mv *_phased.vcf PhasedSNPs/
mv *_allele*.fasta PhasedAlleles/
mv *_coverage.txt Coverage/
mv *_0cov.bed Coverage/
```

## Selection Analysis

### filter coverage

``` bash
mv Cgodm Cgodm_goodone #changed the name of the folder
cd Cgodm_goodone
```

``` bash
qsub Cgodm_filter.pbs
```

filter\_cov.py is a python script that filter the contigs by its
coverage, the ouput is a csv file with the name of the gene and a column
that tells if a gene is present or absent based on the percentage of
locus with a threshold coverage, in this case 5 % or more with 0
coverage.

``` bash
#PBS -N filter_Cgodm
#PBS -l select=6:ncpus=20:interconnect=any:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

source .bash_profile

module load gnuparallel/20210222

cd $PBS_O_WORKDIR

parallel -a ../ID_Cgodm --sshloginfile $PBS_NODEFILE -j1 "echo {}
    source /home/ramsesr/.bashrc
    module load anaconda3/5.1.0-gcc/8.3.1
    source activate VariantsCalling_env
    cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm
    cd {}
    pwd
    filter_cov.py -i {}_coverage.txt -o {}"
```

use the .csv to get a list of the present genes

``` bash
cd /RESULTS/Coverage
grep -h ",present" Cgodm*.csv | sort | uniq -c >present_genes.csv
#grep -E "[3456] " present_genes.csv | cut -d " " -f 8 > tmp
grep -E "[3456] " present_genes.csv> tmp

#### this section is not necesarry as I used the cleaned consensus as reference.
cp ../../../../Cgodm_transcriptome/02rsem_clean/remove.list .

#### this is necesary as we use the file tmp_final later
cat tmp > tmp_final

#### this section is not necesarry as I used the cleaned consensus as reference.
for i in `cat remove.list`
do grep -c "TOXIN" tmp_final
grep "${i}" tmp_final
grep -v "${i}" tmp_final > tmp1
cp tmp1 tmp_final
grep -c "TOXIN" tmp_final
done
rm tmp1
#################
```

### Run Analysis(Tajimas’ D, *F*<sub>*S**T*</sub>, *π* and snpEff)

#### snpEff reference

Previous to run snpEff do this

``` bash
###move consense transcriptome to Geneious, chose option of making a list
###select one then cmnd + a to select all, then annotate all, use the name CDS as name 
### save as genes.gff
####remove the comment lines with this code 
perl -pi -e 's/^#.*\n//g' genes.gff
#### remove the "created by" with this code 
perl -pi -e 's/created\ by\=ramsesrosales//g' genes.gff
##Add some information to the end of each line 
perl -pi -e 's/^(.*?)(\t.*)/$1$2gene_id "$1"; transcript_id "1"/gm' genes.gff
#### send the file to palmetto
scp genes.gff ramsesr@xfer01-ext.palmetto.clemson.edu:/zfs/venom/Ramses/bin/snpEff/data/Cgodm

###also copy ref.fa into the snpEff/data/Cgodm folder and name it sequences.fa

#####go to the snpEff package
cd /home/ramsesr/.conda/pkgs/snpeff-5.0-hdfd78af_1/share/snpeff-5.0-1

#change the path of the databases for the path you are going to use to save your files in this case (/zfs/veom/Ramses/bin/snpEff/data/Cgodm)
nano snpEff.config 
#change data.dir = ./data for 
data.dir = /zfs/veom/Ramses/bin/snpEff/data/
echo "Cgodm.genome : Cgodm" >> snpEff.config

##now we have to build the data base
snpEff build -gff3 -v Cgodm

#####now this should work
snpEff Cgodm Combined.vcf > Annotated.vcf
```

#### Run Analysis

``` bash
qsub Cgodm_vcf3.pbs
```

``` bash
#PBS -N Cgodm_vcf_tools3.0
#PBS -l select=20:ncpus=20:interconnect=fdr:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

module load anaconda3/2021.05-gcc/8.3.1
source activate VariantsCalling_env

cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone

##first a test with Cgodm
#in Variant/Cgodm
mkdir test
cp -r RESULTS/PhasedSNPs test/

cd test/PhasedSNPs
####tabix need a compresed file to work, first I need to run bgzip
for i in *_phased.vcf
do echo $i
bgzip $i
done
##now tabix should work
for i in *
do echo $i
tabix -fp vcf $i
done

# Merge SNPs into a single file
vcf-merge -R 0/0 *_phased.vcf.gz  > Combined.vcf
bgzip -c Combined.vcf > Combined.vcf.gz
tabix -fp vcf Combined.vcf.gz
mkdir Analysis_filtered
mv Combined.* Analysis_filtered/
cd Analysis_filtered

## Tajima's D and Nucleotide Diversity per Gene
mkdir Genes NDiversity TajimasD_perGene RESULTS
cp /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/RESULTS/Coverage/tmp_final ./Genes.txt
perl -p -i -e "s/^.*[\d]\ //g" Genes.txt
perl -p -i -e 's/(.*),present[\s]*/$1\n/g' Genes.txt


###code to look to do the Genes.txt file
##contig=<ID=
##trinityContig9992||sp|Q8IXM2.1|,assembly=ref.fa,length=531.vcf

####run parallel
conda deactivate 
module purge 
module load parallel/20190222-gcc/8.3.1

parallel -a Genes.txt --sshloginfile $PBS_NODEFILE -j 1 --verbose "echo {}
source .bashrc
module load anaconda3/2021.05-gcc/8.3.1
source activate VariantsCalling_env
cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/test/PhasedSNPs/Analysis_filtered
tabix -h Combined.vcf.gz "{}" > Genes/"{}".vcf
vcftools --vcf Genes/"{}".vcf --TajimaD 8000 --out TajimasD_perGene/"{}"
vcftools --vcf Genes/"{}".vcf --site-pi --out NDiversity/"{}"
tail -n+2 TajimasD_perGene/"{}".Tajima.D >> RESULTS/TajimaD_perGene.txt
tail -n+2 NDiversity/"{}".sites.pi >> RESULTS/Pi_perSite.txt"
 
module purge
module load anaconda3/2021.05-gcc/8.3.1
source activate VariantsCalling_env
cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/test/PhasedSNPs/Analysis_filtered

# TajimaD per site
mkdir TajimaD_perSite
vcftools --vcf Combined.vcf --TajimaD 1 --out TajimaD_perSite/TajimaD_perSite
cp TajimaD_perSite/TajimaD_perSite.Tajima.D RESULTS/TajimaD_perSite.txt

## Fst Between Populations
### get Toxin and Nontoxins.txt
grep 'TOXIN' Genes.txt > Toxins.txt
cp Genes.txt  Nontoxins.txt
perl -p -i -e 's/TOXIN.*\n//g' Nontoxins.txt

# Fst per site between populations (mean per gene in R)
tabix -h Combined.vcf.gz `cat Toxins.txt` > CombinedToxins.vcf
tabix -h Combined.vcf.gz `cat Nontoxins.txt` > CombinedNontoxins.vcf

#### make the list of individuals for each population

printf 'Cgodm-CLP2359\nCgodm-CLP2360\nCgodm-CLP2362\nCgodm-CLP2388\n' > Cgodm-South
printf "Cgodm-CLP2377\nCgodm-CLP2378\n" > Cgodm-North

####calculate fst
mkdir Fst
vcftools --vcf CombinedToxins.vcf --weir-fst-pop Cgodm-North --weir-fst-pop Cgodm-South --out Fst/Toxin_SPvNP &> Fst/Toxin_SPvNP.log
vcftools --vcf CombinedNontoxins.vcf --weir-fst-pop Cgodm-North --weir-fst-pop Cgodm-South --out Fst/Nontoxin_SPvNP &> Fst/Nontoxin_SPvNP.log

### add a tab and the name of the comparation, in this case is just one
perl -p -i -e "s/\n/\t SPvNP\n/g" Fst/*_SPvNP.weir.fst
#perl -p -i -e "s/SPvNP\t SPvNP\n/SPvNP\n/g" Fst/*_SPvNP.weir.fst
tail -n+2 Fst/*_SPvNP.weir.fst >> Fst.txt

###this was in a loop in the original code
A="$(grep "mean Fst" Fst/Nontoxin_SPvNP.log)"
B="$(grep "mean Fst" Fst/Toxin_SPvNP.log)"
C="$(grep "weighted Fst" Fst/Nontoxin_SPvNP.log)"
D="$(grep "weighted Fst" Fst/Toxin_SPvNP.log)"
printf "SPvNP \t $A \t $B \t $C \t $D \n" >> RESULTS/Fst_PopMeans.txt


grep -v "==" Fst.txt | grep -v '^$' > Fst_perSite.txt

#####now this should work
snpEff Cgodm Combined.vcf > Annotated.vcf

####now following the html of Ccera
mkdir PopulationAnalyses
cd PopulationAnalyses
###make a directory for each population
mkdir Cgodm-North
mkdir Cgodm-South

###repeat the previous steps to get all the data for each population
for i in `cat ../Cgodm-North`
do echo $i
cp ../../${i}* Cgodm-North
done

for i in `cat ../Cgodm-South`
do echo $i
cp ../../${i}* Cgodm-South
done

echo Cgodm-North >>pop_list
echo Cgodm-South >>pop_list

for i in `cat pop_list`
do echo $i
cd $i
vcf-merge -R 0/0 *_phased.vcf.gz > Combined.vcf
bgzip -c Combined.vcf > Combined.vcf.gz
tabix -fp vcf Combined.vcf.gz
mkdir Analysis
cp ../../Genes.txt Analysis/
cd Analysis/
mkdir Genes NDiversity TajimasD
    parallel -a Genes.txt -j 20 --verbose "echo {}
    tabix -h ../Combined.vcf.gz {} > Genes/"{}".vcf
    vcftools --vcf Genes/"{}".vcf --TajimaD 8000 --out TajimasD/"{}"
    vcftools --vcf Genes/"{}".vcf --site-pi --out NDiversity/"{}"
    tail -n+2 TajimasD/"{}".Tajima.D >> CombinedTajimaD.txt
    tail -n+2 NDiversity/"{}".sites.pi >> CombinedPi.txt"
cd ../../
done

cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/test/PhasedSNPs/Analysis_filtered

##### Tajima's D Synonymous vs. Nonsynonymous Mutations
grep "synonymous_variant" Annotated.vcf > Synonymous.vcf
grep -v "synonymous_variant" Annotated.vcf > Nonsynonymous.vcf
# Manually edit the Synonymous.vcf to add header back into document

# Headers are in line  7276
grep -n '#CHROM' Nonsynonymous.vcf > tmp
perl -p -i -e 's/([0-9]*):.*/$1/g' tmp
grep -m `cat tmp` -e '.*' Nonsynonymous.vcf > tmp.vcf
cat Synonymous.vcf >> tmp.vcf
rm Synonymous.vcf
mv tmp.vcf Synonymous.vcf


bgzip -c Synonymous.vcf > Synonymous.vcf.gz
tabix -fp vcf Synonymous.vcf.gz
bgzip -c Nonsynonymous.vcf > Nonsynonymous.vcf.gz
tabix -fp vcf Nonsynonymous.vcf.gz
mkdir TajimasD_Synonymous
mkdir TajimasD_Nonsynonymous
mkdir Genes_Synonymous
mkdir Genes_Nonsynonymous

conda deactivate 
module purge 
module load parallel/20190222-gcc/8.3.1

parallel -a Genes.txt --sshloginfile $PBS_NODEFILE -j 1 --verbose "echo {}
    source .bashrc
    module load anaconda3/2021.05-gcc/8.3.1
    source activate VariantsCalling_env
    cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/test/PhasedSNPs/Analysis_filtered
    tabix -h Synonymous.vcf.gz {} > Genes_Synonymous/"{}".vcf
    vcftools --vcf Genes_Synonymous/"{}".vcf --TajimaD 8000 --out TajimasD_Synonymous/"{}"
    tail -n+2 TajimasD_Synonymous/"{}".Tajima.D >> RESULTS/TajimasD_Synonymous.txt
    tabix -h Nonsynonymous.vcf.gz {} > Genes_Nonsynonymous/"{}".vcf
    vcftools --vcf Genes_Nonsynonymous/"{}".vcf --TajimaD 8000 --out TajimasD_Nonsynonymous/"{}"
    tail -n+2 TajimasD_Nonsynonymous/"{}".Tajima.D >> RESULTS/TajimasD_Nonsynonymous.txt"
```

to get the lenght of the transcripts in the same order than in the
Genes.txt file runt this

``` bash
cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/test/PhasedSNPs/Analysis_filtered
./lenght.sh
```

content of the script

``` bash
#!/bin/bash

##this code gets the values of lenght in the same order then Genes.txt

grep  '##contig=<ID=' Combined.vcf > tmp1
for i in `cat Genes.txt`
do echo $i
grep -e ${i} tmp1 >> length.txt
done

perl -p -i -e 's/##contig.*length=([\d]*)>/$1/g' length.txt
```

the tables in the results have random blank lines, to correct that run
the next script

``` bash
cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/test/PhasedSNPs/Analysis_filtered/RESULTS
./correct_tables.sh 
```

content of the script

``` bash
#correct TXT tables, for some reason have randon extra spaces and I dont know why

tail -n +2 -q ../TajimasD_perGene/* > TajimaD_perGene.txt
tail -n +2 -q ../TajimasD_Nonsynonymous/* > TajimasD_Nonsynonymous.txt
tail -n +2 -q ../TajimasD_Synonymous/* > TajimasD_Synonymous.txt
tail -n +2 -q ../TajimaD_perSite/* > TajimaD_perSite.txt
tail -n +2 -q ../NDiversity/* > Pi_perSite.txt
```

### Run Analysis(HyPhy BUSTED)

#### Clip tree

``` bash
# renames the tree scaled in Best to test.tree
# -i input, -l list with the tip labels to keep, -m mode(keep, or clip), -o output
Clip_trees.R -i Cerro_core_alignment.phy.contree -l ../../../../../../ID_Cgodm -m keep -o Cgodm_IQtree_scaled.tree
```

``` bash
(Cgodm-CLP2359:0.0010376778,((Cgodm-CLP2360:0.0010322517,(Cgodm-CLP2377:0.0014207813,Cgodm-CLP2378:0.0013584519)0.46/100:0.0020232332)100:0.000263894,Cgodm-CLP2362:0.0010620274)0.99/100:0.0001413664,Cgodm-CLP2388:0.0009634051)0.99;
```

#### Run Analysis

``` bash
qsub Cgodm_HyPhy_busted_4.pbs
```

The scripts Linearize.py and RemoveStop.py are available in
<https://github.com/RhettRautsaw/Bioinformatics/tree/master/scripts>

``` bash
#PBS -N Cgodm_HyPhy_busted4.0
#PBS -l select=20:ncpus=20:interconnect=fdr:mem=100gb,walltime=72:00:00
#PBS -M ramsesr@g.clemson.edu
#PBS -m abe
#PBS -j oe

#######
cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/
cd test
mkdir HyPhy
cd HyPhy
mkdir genes

cp /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/RESULTS/Coverage/tmp_final ./Genes.txt
perl -p -i -e "s/^.*[\d]\ //g" Genes.txt
perl -p -i -e 's/(.*),present[\s]*/$1\n/g' Genes.txt


cp /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/RESULTS/PhasedAlleles/*_allele1.fasta .

#download Linearize to my path, it is in the bioinformatics module in Rhett's github
#add to the path

#use Linearize.py to make a file with only one line for each sequence
module load anaconda3/2021.05-gcc/8.3.1
source activate bio

for i in *.fasta
do echo $i
Linearize.py -f $i
done


#add the name of the sample at the start of each gene
for i in `cat /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/ID_Cgodm`
do echo $i
perl -p -i -e "s/>/>${i}_/g" ${i}_allele1.fasta
done



## use the Genes list and the list of the samples to extract the sequence of each sample of a specific gene and save it in a fasta file
for i in `cat Genes.txt`
do echo $i
for j in `cat ../../../ID_Cgodm`
do echo $j
## in grep -h to make sure it does not include the name of the file
# -A1 means that it includes the line that match the pattern and a line after it
# -e means that the pattern is a regexp, it makes sure grep recognise ${i} as a varoable
grep -h -A1 -e ${i} ${j}_allele1.fasta >> genes/Cgodm_${i}.fasta
done
done

##this needs bio environment
#remove the stop codons

ls genes >> tmp.list
for i in   `cat tmp.list`
do echo $i
cd genes
RemoveStop.py -f $i
cd ..
done

#remove the anotation of each genes and only let the name of the individuals
for i in `cat ../../../ID_Cgodm`
do echo $i
perl -p -i -e "s/>${i}_.*/>${i}/g" genes/*noStop.fasta
done

module load parallel/20190222-gcc/8.3.1
cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/test/HyPhy
mkdir Busted_Busco_Scaled_tree

echo '(Cgodm-CLP2359:0.0010376778,((Cgodm-CLP2360:0.0010322517,(Cgodm-CLP2377:0.0014207813,Cgodm-CLP2378:0.0013584519)0.46/100:0.0020232332)100:0.000263894,Cgodm-CLP2362:0.0010620274)0.99/100:0.0001413664,Cgodm-CLP2388:0.0009634051)0.99;' > busted_busco_scaled.tree
mv busted_busco_scaled.tree Busted_Busco_Scaled_tree

parallel -a tmp1.list --sshloginfile $PBS_NODEFILE -j1 "echo {}
source .bashrc
module load anaconda3/2021.05-gcc/8.3.1
source activate HyPhy_env
cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/test/HyPhy
cd Busted_Busco_Scaled_tree
pwd
cp ../genes/"{}" /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/test/HyPhy/Busted_Busco_Scaled_tree
hyphy busted --alignment "{}" --tree busted_busco_scaled.tree
rm "{}"
echo busten model finalized"

cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/test/HyPhy/Busted_Busco_Scaled_tree
/zfs/venom/Ramses/bin/RamsesScripts/busted_maketable.sh
```

code of busted\_maketable.sh

this script make a table from the resulting .json files

``` bash
#!/bin/bash

echo -e gen,LRT,p-value,CW1,CW2,CW3,UCW1,UCW2,UCW3 > Busted_Results
for i in *.BUSTED.json
do echo $i
grep -A29 -e 'Constrained model' ${i} > tmp.file
grep -A29 -e 'Unconstrained model' ${i} > tmp1.file
grep -n -e ".*" tmp.file > tmp1.1.file
grep -n -e ".*" tmp1.file > tmp1.2.file
echo 0`grep -e '^21' tmp1.1.file`,0`grep -e '^25' tmp1.1.file`,0`grep -e '^29' tmp1.1.file`,0`grep -e '^21' tmp1.2.file`,0`grep -e '^25' tmp1.2.file`,0`grep -e '^29' tmp1.2.file` > tmp2.file
perl -p -i -e "s/21\: //g" tmp2.file
perl -p -i -e "s/25\: //g" tmp2.file
perl -p -i -e "s/29\: //g" tmp2.file
perl -p -i -e "s/\"omega\"\://g" tmp2.file
perl -p -i -e "s/,\n/\n/g" tmp2.file
perl -p -i -e "s/,00\./,0\./g" tmp2.file
perl -p -i -e "s/^00\./0\./g" tmp2.file
perl -p -i -e 's/,0([0-9]+)\./,$1\./g' tmp2.file
perl -p -i -e 's/0([0-9]+),,/$1,,/g' tmp2.file
echo $i,`grep -e '\"LRT\"\:' $i`,`grep -e '\"p-value\"\:' $i`, `cat tmp2.file` >> Busted_Results
done
perl -p -i -e "s/^Cgodm_//g" Busted_Results
perl -p -i -e "s/_noStop\.fasta\.BUSTED\.json//g" Busted_Results
perl -p -i -e "s/\"LRT\"\://g" Busted_Results
perl -p -i -e "s/\"p-value\"\://g" Busted_Results
perl -p -i -e "s/,,/,/g" Busted_Results
```

The program RemoveStops cut the names of some of my files, this will
make and error later when we try to relate the results from BUSTED with
the rest of the results, so to avoid that I ran these script to correct
the changed names withing the directory with the Busted\_Results file.

``` bash
./correct_table.sh
```

The directory assumes you have a file called Genes.txt with the right
names from the genes in the parent directory.

``` bash
mv Busted_Results Busted_Results_original
head -n 1 Busted_Results_original > Busted_Results
for i in $(tail -n +2 Busted_Results_original | cut -d "," -f 1)
do echo $i
grep $i ../Genes.txt 
echo $(grep ${i} ../Genes.txt),$(grep ${i} Busted_Results_original | cut -d "," -f 2-9) >> Busted_Results
done
```

## Prepare Final Data

### Correct tables

The tables of Tajima’s D and Pi have some random changes of line, so I
did a small script to generate the corrected tables

``` bash
cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/test/PhasedSNPs/Analysis_filtered/RESULTS
./correct_tables.sh
```

``` bash
#correct TXT tables, for some reason have randon extra spaces and I dont know why

tail -n +2 -q ../TajimasD_perGene/* > TajimaD_perGene.txt
tail -n +2 -q ../TajimasD_Nonsynonymous/* > TajimasD_Nonsynonymous.txt
tail -n +2 -q ../TajimasD_Synonymous/* > TajimasD_Synonymous.txt
tail -n +2 -q ../TajimaD_perSite/* > TajimaD_perSite.txt
tail -n +2 -q ../NDiversity/* > Pi_perSite.txt
```

The Busted\_Results table had a problem with the names as the program
RemoveStop.py change the ouput names if there is a “.” somewere in the
name of the fasta file, run the next script to fix these, it uses the
Gene.txt list that has to be in the parent directory.

``` bash
cd /zfs/venom/Ramses/Cerrophidion/Variants/Cgodm_goodone/test/HyPhy/Busted_Busco_Scaled_tree
./correct_table.sh
```

``` bash
mv Busted_Results Busted_Results_original
head -n 1 Busted_Results_original > Busted_Results
for i in $(tail -n +2 Busted_Results_original | cut -d "," -f 1)
do echo $i
grep $i ../Genes.txt 
echo $(grep ${i} ../Genes.txt),$(grep ${i} Busted_Results_original | cut -d "," -f 2-9) >> Busted_Results
done
```

### Merge tables

In the same file put the next files and a list with the names of the
files (TDlist) The snpEff where put in Excel to be saved as csv In Excel
eliminate the first line and the \# at the start of the headers in
snpEff\_genes.csv use a perl pie to eliminate the first CDS\_

``` bash
perl -pi -e "s/^CDS_//g" snpEff_genes.csv
```

for the snpEff\_Summary.csv we copyed the table from the html, excel
added some weird characters at the end of each cell I use find and
replace in VSC to eliminate that character. The file Busted\_Results is
a csv, I just add the extention

``` bash
TajimaD_perGene.txt
TajimasD_Nonsynonymous.txt
TajimasD_Synonymous.txt
snpEff_Summary.csv
snpEff_genes.csv
Pi_perSite.txt
Fst_perSite.txt
Busted_Results.csv
```

We also need to copy the Genes.txt and the length.txt file for the first
script that concatenate the files in TDlist in one data base

The second script adds the Rsem results and the differential expression
result, for these we need to have both files in the directory

> -   Cgodm\_TajimasD.csv
> -   Cgodm\_consense\_df.csv

``` bash
r_4
make_FD_1.R -g Genes.txt -e length.txt -l TDlist -o Cgodm -s linux

make_FD_2_final.R -t Cgodm_TajimasD.csv -r Cgodm_consense_df.csv -d Cgodm_DifExp_tab.csv -o Cgodm
```

# Final

Here we end the scripts in linux, we add a second rmarkdown that include
the Ranalysis using the Cgodm\_finaldata.csv file.
