---
layout: page
title: Theoretical and Practical HiC Workshop
subtitle: Quality control of NGS data
minutes: 5
---

> ## Learning objectives
>
> *   Understand the FastQ format.
> *   Understand the concept of sequence quality. 
> *   Learn how to use software to measure and improve the quality of massive sequencing data.

## Downloading the data

The first thing we do when receiving our data is to download them. Our collaborator has provided the [link](https://drive.google.com/open?id=0B9ZVSRlHL8cIbm5EUXczdzM4a2M). Once downloaded we will copy them into our docker container using the following command. Remember to replace the id by your docker id:

~~~ {.bash}
$ docker cp FastQC_Short.tar.gz a646764e06a0:/usr/local/data/FastQC_Short.tar.gz
~~~

They are usually in a compressed format. We can decompress the data using the `tar` command:

~~~ {.bash}
$ tar -xvf FastQC_Short.tar.gz
~~~

This command decompresses a directory called `FastQC_Short`.
We enter that directory:

~~~ {.bash}
$ cd /usr/local/data/FastQC_Short
~~~
 
Revise its contents:

~~~ {.bash}
$ ls
~~~

~~~ {.output}
Partial_SRR2467141.fastq 
Partial_SRR2467142.fastq 
Partial_SRR2467143.fastq 
Partial_SRR2467144.fastq 
Partial_SRR2467145.fastq 
Partial_SRR2467146.fastq
Partial_SRR2467147.fastq
Partial_SRR2467148.fastq
Partial_SRR2467149.fastq
Partial_SRR2467150.fastq
Partial_SRR2467151.fastq
~~~

> ## The name of fastq files {.callout}
>
> In this case, the name of the files has been assigned in a
> arbitrary but often contain important information.
> For example, the Illumina files generally have the following format:
>
> ~~~ {.output}
> <sample name>_<barcode>_L<lane(3 digits)>_R<read>_<set number (3 digits)>.fastq.gz
> ~~~
> Example: NA10831_ATCACG_L002_R1_001.fastq.gz
>
> It is ideal to respect these names to avoid losing information.
 
Let's review one of the files using the head command. It usually shows us the first
10 lines of a file but with the `-n` flag will show us the first 12:

~~~ {.bash}
$ head -n 12 Partial_SRR2467141.fastq 
~~~

~~~ {.output}
@SRR2467141.1 SALLY:472:C6NCDACXX:8:1101:1392:1873 length=101
NTTATTTTTTTCGTTCTTCTTGAAGAAGACGTTACCTACGGCGTATTTGCCCATCTCAGGTATGTCTAGATCAAGATCTAACTTGAATTCTCTTTTCATAA
+SRR2467141.1 SALLY:472:C6NCDACXX:8:1101:1392:1873 length=101
#1:=BDDDHDHB?<?CFGGGC9@FF@GGGG>EEEDGDHGFHGE;AEFH>AC@D;@B>C>CCC@C>>DDCC3:>AA5>CC>>CCD>@CCDCDCCCCC@C@>C
@SRR2467141.2 SALLY:472:C6NCDACXX:8:1101:1326:1950 length=101
CACCCATTGACTGGCCAAATGCCCCATTATTTTGAGTATTGTTATTTCCAAATAAACTGTTACTATTACTGCCAGCGGCAGAAGTGAATCCACAGATCGGA
+SRR2467141.2 SALLY:472:C6NCDACXX:8:1101:1326:1950 length=101
??@DFFFFHHFH?GEHEHIDEHCCFEHIE@GIGCHG<DGGIHIGIGGIGHIIFIHFFGDHGGIIHIGIIIIGGEH@EADB>=AA>3@>CCCCCCBC@C###
@SRR2467141.3 SALLY:472:C6NCDACXX:8:1101:1477:1959 length=101
CATCTTTTCTTTAGGCACATCATCCTGATAAGTGTACTTACCAGGATATATACCATCGGTATTGATGTTATCGGCATCACATAAAACTAATTCACCAGAAA
+SRR2467141.3 SALLY:472:C6NCDACXX:8:1101:1477:1959 length=101
@CCFFFFFHHHHHJJJIIJJIJJJJJJEIJJJIIJJJIJIJJIJJIIIJJJIIIHIIIJDHIJJIGJJIJJJJJIHHHFFFFFEEEEEEDDEDDDDDDDDD
~~~

This file is in `fastq` format. This is the format that most modern sequencing platforms generate. This format is a modification of the sequence format
[Fasta](https://en.wikipedia.org/wiki/FASTA_format). You can recognize them because they end
with the extension `.fq` or `.fastq`.

FastQ files are plain ASCII text files that contain both the sequences
detected by the sequencer as well as information about the quality of each one
of the sequenced nucleotides.

> ## Do not change fastq files! {.callout}
>
> The fastq files that you receive after a sequencing experiment will be
> requested when you want to publish an article. Ideally, make a
> couple of copies on different devices as soon as you receive them and *do not modify them manually*.
> 
> It is also recommended that, if you have the information at hand, you create a small
> README file (plain text) that is in the same directory as your fastq files
> and describe how the experiment was performed, how many replicas there are of each condition, which
> indexes and controls were used, etc. This will help you greatly when trying to
> interpret your results.

We will take a single sequence as an example:

~~~ {.output}
@SRR2467141.1 SALLY:472:C6NCDACXX:8:1101:1392:1873 length=101
NTTATTTTTTTCGTTCTTCTTGAAGAAGACGTTACCTACGGCGTATTTGCCCATCTCAGGTATGTCTAGATCAAGATCTAACTTGAATTCTCTTTTCATAA
+SRR2467141.1 SALLY:472:C6NCDACXX:8:1101:1392:1873 length=101
#1:=BDDDHDHB?<?CFGGGC9@FF@GGGG>EEEDGDHGFHGE;AEFH>AC@D;@B>C>CCC@C>>DDCC3:>AA5>CC>>CCD>@CCDCDCCCCC@C@>C
~~~

We see that each sequence is represented by 4 lines.

1. **Sequence name** - This line begins with the symbol `@`, followed by the read name. 
2. The nucleotide **sequence**. 
3. **Second title** - This line starts with the symbol `+`. Generally the information is the same as in the first line but it can also be blank as long as it begins with the symbol `+`.
4. **Quality information** - It contains an ASCII character string. Every
one of these characters corresponds to a nucleotide of the sequence and represents the quality of the same. The score of each nucleotide indicates the
level of confidence that you have in that base. High levels indicate that the base has been reported correctly while low levels suggest
that there is uncertainty about the actual sequence in this position.

The name of the sequence also contains important information. In
In this case, the Illumina data provide the following information:

~~~ {.output}
@<instrument>:<run>:<flow cell identifier>:<lane>:<square>:<x-coord>:<y-coord> <read>:<if the ead was filtered>:<control number>:<sequence index> length=<sequence size>
~~~

The level of quality represented in the fourth line is called Phred score. In its
original format a Phred score is simply the logarithm of the probabilities of
error:

~~~ {.output}
Phred score = - 10 * log10(error probability)
~~~

| Quality score | Error probability |
|:----------------------|:----------------------|
| Q40 0.0001 | (1 en 10,000) |
| Q30 0.001 | (1 en 1,000) |
| Q20 0.01 | (1 en 100) |
| Q10 0.1 | (1 en 10) |

Wait, if the Phred scores are numbers, why do the quality scores shown
in line four are a mixture of alphanumeric symbols?

~~~ {.output}
#1:=BDDDHDHB?<?CFGGGC9@FF@GGGG>EEEDGDHGFHGE;AEFH>AC@D;@B>C>CCC@C>>DDCC3:>AA5>CC>>CCD>@CCDCDCCCCC@C@>C
~~~

This is because these scores are encoded in ASCII (American Standard Code for Informational Interchange).
In short, ASCII is a coding system that equates a number to a character
alphanumeric. For example, the character 'A' is represented by the number 65 in the table
ASCII code, while '%' is represented by the number 37. This system allows us to
represent a total of 256 different characters.

Why use ASCII? Given the huge amount of data produced during massive parallel sequencing, it's all about reducing the data to the maximum. The ASCII system allows us to represent
two-digit numbers in a single bit (8 bits), which reduces the space they occupy
these files. If you are interested in the subject, you can read more about binary coding.

Finally, since this quality coding system was invented, which is already
used with Sanger-type sequencing, different companies have "slipped" the scales
of ASCII conversion because what is important to verify with the company or laboratory
in which version of this scale have been coded to use the correct conversion. Some
of the quality control tools that we will use infer the type of coding
used directly from the data provided.

## Verifying the quality of the sequences

We will verify the quality of the sequences using the program [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
This is one of the most used programs for this type of analysis,
since it provides a global vision of how the data is structured, as well as
being pretty fast.

FastQC reads a series of sequencing files and produces a report of
quality control for each one, which consists of a number of different modules,
each of which helps us identify different problems in our data.

Let's analyze our first sequencing file,
first we created a directory to store our results:

~~~ {.bash}
$ mkdir QUAL
~~~

And we perform our quality analysis using FastQC:

~~~ {.bash}
$ fastqc -O ./QUAL/ Partial_SRR2467141.fastq 
~~~

~~~ {.output}
Started analysis of Partial_SRR2467141.fastq
Approx 20% complete for Partial_SRR2467141.fastq
Approx 40% complete for Partial_SRR2467141.fastq
Approx 60% complete for Partial_SRR2467141.fastq
Approx 80% complete for Partial_SRR2467141.fastq
Analysis complete for Partial_SRR2467141.fastq
~~~ 

Once finished, let's review the result:

~~~ {.bash}
$ cd QUAL
$ ls
~~~

~~~ {.output}
Partial_SRR2467141_fastqc.html 
Partial_SRR2467141_fastqc.zip
~~~

The easiest way to explore these results is by opening the html file in
your browser. You can download them to a working directory in your computer using the command:

~~~ {.bash}
$ docker cp a646764e06a0:/usr/local/data/FastQC_Short/QUAL/Partial_SRR2467141_fastqc.html .
~~~
You can open it by double clicking on the file.

The results should be similar to [this](http://liz-fernandez.github.io/transcriptome_analysis/Partial_SRR2467141_fastqc.html).

This page contains a lot of information broken down into the following sections:

1. **Basic Statistics** -  The basic statistics of the experiment.
2. **Per base sequence quality** - Box diagrams showing the quality of each base.
3. **Per tile sequence quality** - It contains the box (tile) from which each sequence comes. It only appears if the analyzes are made in an Illumina library that retains its original identifiers.
4. **Per sequence quality scores** - It allows to see if a subgroup of sequences has bad quality.
5. **Per base sequence content** - Shows the proportion of each base in each position. 
6. **Per sequence GC content** - Compares the observed GC distribution with a modeled GC distribution.
7. **Per base N content** - Indicates how many nucleotides could not be interpreted and are indicated with an N.
8. **Sequence Length Distribution** - Shows the sequence length distribution.
9. **Sequence Duplication Levels** - Indicate how many times each sequence is repeated. It is done only in the first 100,000 sequences.
10. **Overrepresented sequences** - Shows over represented sequences.
11. **Adapter Content** - Shows the content of adapters in the library.

Let's review the results once more by clicking on the [link](http://liz-fernandez.github.io/transcriptome_analysis/Partial_SRR2467141_fastqc.html). 

In general we see that most of the criteria have a green tick of quality control pass, with the exception of 
per base sequence content, GC content per base and overrepresented sequences.

Although we did not provide this information, the program has inferred the coding
ASCII (Sanger / Illumina 1.9) as well as the number of sequences, their size and
its GC content.

When we observe the **quality per base**, we see that most of the bases are in the area
green or have good quality. Each base is represented by a box diagram that
provides a good idea about the distribution of quality scores in all
the sequences. It is evident that the quality decreases more markedly
at the end of the sequence. This is known as small errors accumulate while
the longer the sequence, this is one of the reasons why sequencers based
in nucleotide incorporation they have a limit on the maximum length of their readings.

The graph showing the **quality per tile** is almost completely blue showing that the average distribution of errors in the tile is very similar.

The graph of **quality score per sequence** shows that most of the sequences have high scores.

However, the **distribution of nucleotides per base** shows that, at the beginning of the sequence, the distribution is very disparate. 

We also observe that:

* The **GC distribution per sequence** is somewhat different to the hypothetical distribution.
* There are few **ambiguous bases (Ns)**.
* All sequences have the same **size**.
* Levels of **duplication of sequences** are very low.
* There are a few **over represented sequences**. Mostly Illumina PCR primers.
* There are no **adapters**.

When running FastQC, a compressed file called
`Partial_SRR2467141_fastqc.zip`. 

We unpack this file using the following command:

~~~ {.bash}
$ unzip Partial_SRR2467141_fastqc.zip
~~~

~~~ {.output}
Archive:  Partial_SRR2467141_fastqc.zip
   creating: Partial_SRR2467141_fastqc/
   creating: Partial_SRR2467141_fastqc/Icons/
   creating: Partial_SRR2467141_fastqc/Images/
  inflating: Partial_SRR2467141_fastqc/Icons/fastqc_icon.png
  inflating: Partial_SRR2467141_fastqc/Icons/warning.png
  inflating: Partial_SRR2467141_fastqc/Icons/error.png
  inflating: Partial_SRR2467141_fastqc/Icons/tick.png
  inflating: Partial_SRR2467141_fastqc/summary.txt
  inflating: Partial_SRR2467141_fastqc/Images/per_base_quality.png
  inflating: Partial_SRR2467141_fastqc/Images/per_tile_quality.png
  inflating: Partial_SRR2467141_fastqc/Images/per_sequence_quality.png
  inflating: Partial_SRR2467141_fastqc/Images/per_base_sequence_content.png
  inflating: Partial_SRR2467141_fastqc/Images/per_sequence_gc_content.png
  inflating: Partial_SRR2467141_fastqc/Images/per_base_n_content.png
  inflating: Partial_SRR2467141_fastqc/Images/sequence_length_distribution.png
  inflating: Partial_SRR2467141_fastqc/Images/duplication_levels.png
  inflating: Partial_SRR2467141_fastqc/Images/adapter_content.png
  inflating: Partial_SRR2467141_fastqc/fastqc_report.html
  inflating: Partial_SRR2467141_fastqc/fastqc_data.txt
  inflating: Partial_SRR2467141_fastqc/fastqc.fo
~~~

If we enter the newly created directory, we can see the report in plain text format,
both the summary:

~~~ {.bash}
$ cd Partial_SRR2467141_fastqc
$ more summary.txt
~~~

~~~ {.output}
PASS    Basic Statistics        Partial_SRR2467141.fastq
PASS    Per base sequence quality       Partial_SRR2467141.fastq
PASS    Per tile sequence quality       Partial_SRR2467141.fastq
PASS    Per sequence quality scores     Partial_SRR2467141.fastq
FAIL    Per base sequence content       Partial_SRR2467141.fastq
PASS    Per sequence GC content Partial_SRR2467141.fastq
PASS    Per base N content      Partial_SRR2467141.fastq
PASS    Sequence Length Distribution    Partial_SRR2467141.fastq
PASS    Sequence Duplication Levels     Partial_SRR2467141.fastq
PASS    Overrepresented sequences       Partial_SRR2467141.fastq
PASS    Adapter Content Partial_SRR2467141.fastq
PASS    Kmer Content    Partial_SRR2467141.fastq
~~~

and the full report file:

~~~ {.bash}
$ more fastqc_data.txt
~~~

We have omitted the content since it is very long.
Being able to review the results in flat files is extremely useful when
results are in a server that does not have a graphical interface.

FastQC programmers have compiled examples of quality sequencing data
diverse Let's analyze them:

* [Good data - Illumina](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)
* [Bad data - Illumina](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)
* [Dimer contamination](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/RNA-Seq_fastqc.html)
* [Small RNAs with read through the adaptor ](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/small_rna_fastqc.html)

> ## How valid is our analysis? {.challenge}
>
> Since the data we used only less than 200,000 sequences, can we
> trust that the complete data will have similar quality? Why or why not?

> ## How do the rest of the files look? {.challenge}
>
> Obtain the fastQC report html files for all other three fastq files.
>

## Cleaning the sequences

There are a number of tools to clean sequences and different researchers have different preferences. Similar to other sequencing problems there is an active development of programs for this purpose.


We will use a program called [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). 
This is a Java based program which, like FastQC, is very fast.
It also has filtering modes for single or pair-end reads.

The program performs various functions, including: 

* Removal of adapters
* Removal of low quality bases at 5 'or 3' of the read
* Scanning of the sequence through a window of 4 nucleotides, cutting when the average quality of the base drops of 15.
* Delete sequences shorter than a particular size
* Cut a specific number of nucleotides at the 5 'or 3' end
* Conversion of quality scores between different versions of Phred.

The specific functions according to the manual are:

* ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
* SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
* LEADING: Cut bases off the start of a read, if below a threshold quality
* TRAILING: Cut bases off the end of a read, if below a threshold quality
* CROP: Cut the read to a specified length
* HEADCROP: Cut the specified number of bases from the start of the read
* MINLEN: Drop the read if it is below a specified length
* TOPHRED33: Convert quality scores to Phred-33
* TOPHRED64: Convert quality scores to Phred-64

The configuration changes depending on whether the sequences are single end or paired-end.

You can find more detailed descriptions of their functions in the [manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).

Let's see how it works. The software is already installed on the docker container.

Let's go back to the folder where our fastq files are.

~~~ {.bash}
$ cd ../..
~~~ 
 
The quality of our test sequences is quite good, but we saw that the quality low at the 3 'end, as well as the disproportionate presence of bases at the end 5'. Let's clean our sequences using Trimmomatic.

We will ask Trimmomatic to cut:

* The Illumina TruSeq single end (SE) adapters from the 3 'end
* The first 10 nucleotides of the 5' end
* Bases of poor quality using a window of 4 nucleotides with a minimum average of Phred of 15
* Sequences less than 60 nucleotides.

~~~ {.bash}
$ trimmomatic SE Partial_SRR2467141.fastq Trimmed_Partial_SRR2467141.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 HEADCROP:10 SLIDINGWINDOW:4:15 MINLEN:60
~~~ 

~~~ {.output}
Automatically using 1 threads
java.io.FileNotFoundException: /home/sfernandez/TRANSCRIPTOMICA/FastQC_Short/TruSeq3-SE.fa (No such file or directory)
	at java.io.FileInputStream.open0(Native Method)
	at java.io.FileInputStream.open(FileInputStream.java:195)
	at java.io.FileInputStream.<init>(FileInputStream.java:138)
	at org.usadellab.trimmomatic.fasta.FastaParser.parse(FastaParser.java:54)
	at org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer.loadSequences(IlluminaClippingTrimmer.java:110)
	at org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer.makeIlluminaClippingTrimmer(IlluminaClippingTrimmer.java:71)
	at org.usadellab.trimmomatic.trim.TrimmerFactory.makeTrimmer(TrimmerFactory.java:32)
	at org.usadellab.trimmomatic.Trimmomatic.createTrimmers(Trimmomatic.java:59)
	at org.usadellab.trimmomatic.TrimmomaticSE.run(TrimmomaticSE.java:303)
	at org.usadellab.trimmomatic.Trimmomatic.main(Trimmomatic.java:85)
Quality encoding detected as phred33
Input Reads: 5000 Surviving: 4856 (97.12%) Dropped: 144 (2.88%)
TrimmomaticSE: Completed successfully
~~~ 

We can see that there were very few sequences that were discarded but still are things that could have affected our analysis and/or assembly. We also see that there is an error:

~~~ {.output}
java.io.FileNotFoundException: /home/sfernandez/CLASS_TEST/FastQC_Short/TruSeq3-SE.fa (No such file or directory)
~~~ 

This is because the sequences of the adapters are not in the same directory as our data. One way to compensate for this is to download the necessary files (these files are provided with trimmomatic but we will download them as we are working on the server)

Download this file and run Trimmomatic again:

~~~ {.bash}
$ wget https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-SE.fa
~~~

~~~ {.bash}
$ trimmomatic SE Partial_SRR2467141.fastq Trimmed_Partial_SRR2467141.fastq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 HEADCROP:10 SLIDINGWINDOW:4:15 MINLEN:60
~~~ 
~~~ {.output}
Automatically using 1 threads
Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
Using Long Clipping Sequence: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
ILLUMINACLIP: Using 0 prefix pairs, 2 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Quality encoding detected as phred33
Input Reads: 5000 Surviving: 4847 (96.94%) Dropped: 153 (3.06%)
TrimmomaticSE: Completed successfully
~~~ 

This result shows us that the program was executed properly and eliminated sequences
which contain the adapters that we indicate. Comparing with our previous results
we see that it found 9 sequences with these adapters.

If we use the command `ls` we can see that a result file called `Trimmed_Partial_SRR2467141.fastq`.

~~~ {.bash}
$ ls
~~~ 

~~~ {.output}
Partial_SRR2467141.fastq
Partial_SRR2467142.fastq
Partial_SRR2467143.fastq
Partial_SRR2467144.fastq
Partial_SRR2467145.fastq
Partial_SRR2467146.fastq
Partial_SRR2467147.fastq
Partial_SRR2467148.fastq
Partial_SRR2467149.fastq
Partial_SRR2467150.fastq
Partial_SRR2467151.fastq
QUAL
**Trimmed_Partial_SRR2467141.fastq**
TruSeq3-SE.fa
~~~ 

Let's now review the results of our sequence cleaning. First we will see the first lines of the clean file.

~~~ {.bash}
$ head -n 12 Trimmed_Partial_SRR2467141.fastq
~~~ 

~~~ {.output}
@SRR2467141.1 SALLY:472:C6NCDACXX:8:1101:1392:1873 length=101
TCGTTCTTCTTGAAGAAGACGTTACCTACGGCGTATTTGCCCATCTCAGGTATGTCTAGATCAAGATCTAACTTGAATTCTCTTTTCATAA
+SRR2467141.1 SALLY:472:C6NCDACXX:8:1101:1392:1873 length=101
HB?<?CFGGGC9@FF@GGGG>EEEDGDHGFHGE;AEFH>AC@D;@B>C>CCC@C>>DDCC3:>AA5>CC>>CCD>@CCDCDCCCCC@C@>C
@SRR2467141.2 SALLY:472:C6NCDACXX:8:1101:1326:1950 length=101
CTGGCCAAATGCCCCATTATTTTGAGTATTGTTATTTCCAAATAAACTGTTACTATTACTGCCAGCGGCAGAAGTGAATCCACAGATC
+SRR2467141.2 SALLY:472:C6NCDACXX:8:1101:1326:1950 length=101
FH?GEHEHIDEHCCFEHIE@GIGCHG<DGGIHIGIGGIGHIIFIHFFGDHGGIIHIGIIIIGGEH@EADB>=AA>3@>CCCCCCBC@C
@SRR2467141.3 SALLY:472:C6NCDACXX:8:1101:1477:1959 length=101
TTAGGCACATCATCCTGATAAGTGTACTTACCAGGATATATACCATCGGTATTGATGTTATCGGCATCACATAAAACTAATTCACCAGAAA
+SRR2467141.3 SALLY:472:C6NCDACXX:8:1101:1477:1959 length=101
HHHJJJIIJJIJJJJJJEIJJJIIJJJIJIJJIJJIIIJJJIIIHIIIJDHIJJIGJJIJJJJJIHHHFFFFFEEEEEEDDEDDDDDDDDD
~~~

We can see that even in the first lines cuts were made, mostly in the 3' end.

Let's now review the global results using FastQC:

~~~ {.bash}
$ fastqc -O ./QUAL/ Trimmed_Partial_SRR2467141.fastq 
~~~

~~~ {.output}
Started analysis of Trimmed_Partial_SRR2467141.fastq
Approx 20% complete for Trimmed_Partial_SRR2467141.fastq
Approx 40% complete for Trimmed_Partial_SRR2467141.fastq
Approx 60% complete for Trimmed_Partial_SRR2467141.fastq
Approx 80% complete for Trimmed_Partial_SRR2467141.fastq
Analysis complete for Trimmed_Partial_SRR2467141.fastq
~~~ 

Open the html result [Trimmed_Partial_SRR2467141_fastqc.html](http://liz-fernandez.github.io/transcriptome_analysis/Trimmed_Partial_SRR2467141_fastqc.html).

> ## What are the differences between raw and clean sequences? {.challenge}
>
> Compare html files between raw and clean sequences. Summarizes:
> What differences are there?
> What do you think is the cause of these differences?

We save our data using:

~~~ {.bash}
docker commit a646764e06a0 lizfernandez/hic-langebio:day1
~~~





