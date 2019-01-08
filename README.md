# Telomeres Shorten with Age
Human telomeres are estimated to be 5,000 - 15,000 base pairs at birth (Sanders 2013). The end replication problem shortens telomeres by approximately 50 bp with each round of cell division (Proctor 2002, Suda 2002). 

#### Simple 1 Telomere Shortening Model
Here's a simple mathematical model for one telomere:

```{r}
# These variables are the starting point. The starting telomere length is 5000. 50 bp are lost / division. 
starting_length = 10000
current_length <- starting_length
bp_loss_per_division = 50
current_division <- 0
number_of_cells <- 1
# This simplified model divides until there are no telomeres left. We will shortly see why this is not the case. 
while(current_length > 0) {
  # There will be 200 divisions for this case. I don't want to overwhelm the screen, so I'm limiting to #1 and #s divisible by 20.
  if(current_division == 0) {
    print(paste("There are ", number_of_cells, " cells after ", current_division, " doublings with a telomere length of ", current_length, " bp."))
  }
  else if (current_division %% 20 == 0) {
    print(paste("There are ", number_of_cells, " cells after ", current_division, " doublings with a telomere length of ", current_length, " bp."))
  }
  # The cells double, the telomere length is shortened by 50 bp, and the new division is initiated.
  number_of_cells <- number_of_cells * 2
  current_length <- current_length - 50
  current_division <- current_division + 1
}
# Here's the final situation.
print(paste("There are ", number_of_cells, " cells after ", current_division, " doublings with a telomere length of ", current_length, " bp."))
```

Here are the last couple of lines of output:

```sh
[1] "There are  1  cells after  0  doublings with a telomere length of  10000  bp."
[1] "There are  1048576  cells after  20  doublings with a telomere length of  9000  bp."
[1] "There are  1099511627776  cells after  40  doublings with a telomere length of  8000  bp."
[1] "There are  1152921504606846976  cells after  60  doublings with a telomere length of  7000  bp."
[1] "There are  1.20892581961463e+24  cells after  80  doublings with a telomere length of  6000  bp."
[1] "There are  1.26765060022823e+30  cells after  100  doublings with a telomere length of  5000  bp."
[1] "There are  1.32922799578492e+36  cells after  120  doublings with a telomere length of  4000  bp."
[1] "There are  1.39379657490816e+42  cells after  140  doublings with a telomere length of  3000  bp."
[1] "There are  1.4615016373309e+48  cells after  160  doublings with a telomere length of  2000  bp."
[1] "There are  1.53249554086589e+54  cells after  180  doublings with a telomere length of  1000  bp."
[1] "There are  1.60693804425899e+60  cells after  200  doublings with a telomere length of  0  bp."
```

#### Simple 1 Telomere Model & Damage Checkpoint
The situation is more complicated than that model. A DNA damage checkpoint will be triggered at around 5k bp. If cell cycle checkpoints are intact, the cell will senesce at around 5 kb. If oncogenic changes have occurred, the cell will divide until around 3k bp. There will be genomic instability, which will lead to more mutations and eventual cell death UNLESS a telomere maintenance mechanism stabilizes the telomeres (Harley 2008, Shay 2012). 

![Harley_2008_Box1a](/Assets/Harley_2008_Box1a.jpg "Harley_2008_Box1a")


Here's that slightly more complicated model:

```r
# These variables are the starting point. The starting telomere length is 5000. 50 bp are lost / division. 
starting_length = 10000
current_length <- starting_length
bp_loss_per_division = 50
current_division <- 0
number_of_cells <- 1
mutations <- 0
dividing <- TRUE
# This simple model divides until senescence is reached. 
while(dividing) {
  # There will be 200 divisions for this case. I don't want to overwhelm the screen, so I'm limiting to #1 and #s divisible by 20.
  if(current_division == 0) {
    print(paste("There is one ", number_of_cells, " cell with a telomere length of ", current_length, " bp and ", mutations, " mutations."))
  }
  else if (current_division %% 20 == 0) {
    # Shay 2012 model uses 1 mutation / 20 divisions 
    mutations <- mutations + 1
    print(paste("There are ", number_of_cells, " cells after ", current_division, " doublings with a telomere length of ", current_length, " bp and ", mutations, " mutations."))
  }
  # cells senesce around a telomere length of 5000 bp Harley 2008 
  else if (current_length < 5000) {
    dividing <- FALSE
    print("This cell has senesced at a telomere length of 5,000 bp")
  }
  # The cells double, the telomere length is shortened by 50 bp, and the new division is initiated.
  number_of_cells <- number_of_cells * 2
  current_length <- current_length - 50
  current_division <- current_division + 1
}
```


Here are the last couple of lines of output:

```sh
[1] "There is one  1  cell with a telomere length of  10000  bp and  0  mutations."
[1] "There are  1048576  cells after  20  doublings with a telomere length of  9000  bp and  1  mutations."
[1] "There are  1099511627776  cells after  40  doublings with a telomere length of  8000  bp and  2  mutations."
[1] "There are  1152921504606846976  cells after  60  doublings with a telomere length of  7000  bp and  3  mutations."
[1] "There are  1.20892581961463e+24  cells after  80  doublings with a telomere length of  6000  bp and  4  mutations."
[1] "There are  1.26765060022823e+30  cells after  100  doublings with a telomere length of  5000  bp and  5  mutations."
[1] "This cell has senesced at a telomere length of 5,000 bp"
```

In the previous model we saw that 5 mutations occurred over 100 population doublings. I based this on an assumption from one of Jerry Shay's oncogenic models. Cellular senescence stopped the cells in my current model from developing more mutations. It is thought that 8-15 mutations are enough for a cell to become cancerous (Shay 2012).

![Shay_2012_Mutations_Cancer](/Assets/Shay_2012_Mutations_Cancer.jpg "Shay_2012_Mutations_Cancer")


#### 23 Telomere Model
But there isn't just one telomere! Human cells have 23 pairs of chromosomes, so lets limit the view to the total telomere length for 23 human chromosomes. We will dive into p and q arms in a later section. It's been reported that the DNA damage checkpoint will be triggered by the shortest telomere (Harley 2008) ... so you should be questioning the validity of this model! Don't worry, it's about to get a lot more complicated :) 

```r
# picking a length of 6000 bp for each chromosome
starting_lengths = as.list(rep(6000, 23))
chromosome_names <- c((1:22), "X")
barplot(as.numeric(starting_lengths), names.arg=chromosome_names, las=2, main ="Telomere Lengths", xlab="Chromosome", ylab="Telomere Lenghts (bp)", ylim=c(0,7000))
abline(h=5000, col="red", lwd=3)
abline(h=3000, col="black", lwd=3)
```

![23_Chromosome_Simple_Lengths](/Assets/23_Chromosome_Simple_Lengths.jpg "23_Chromosome_Simple_Lengths")

#### Telomere Length is Linearly Related to Chromosome Length
A DNA damage checkpoint will get triggered by the shortest telomere (Harley 2008). So a better mathematical model would take into account the distribution of telomere lengths and each individual telomere's shortening. See table 1 of Suda 2002 for the estimated telomere lengths of chromosomes 1-22 for the GM130B cell line. It is a telomerase positive (TEL+) spontaneously immortalized lymphoblast line from a male human. 

```r
# Suda 2002 Table 1 chromosome size (megabases) telomere length (bp).
chromosome_lengths <- c(246, 202, 193, 184, 173, 160, 146, 125, 110, 103, 100, 93, 85, 81, 62, 67, 48, 52)
telomere_lengths <- c(5681, 4987, 5018, 4589, 4302, 4127, 3922, 3708, 4045, 3624, 3460, 3109, 3077, 3007, 2750, 2913, 2735, 2806)
# these should be 18 long cause that's how many rows Suda 2002 has
length(chromosome_lengths)
length(telomere_lengths)
# plot the relationship
plot(chromosome_lengths, telomere_lengths, main ="Telomere and Chromosome Lengths are Linearly Related", xlab="Chromosome Lengths (Mb)", ylab="Telomere Lenghts (bp)")
```

![telomere_chromosome_linear](/Assets/telomere_chromosome_linear.jpg "telomere_chromosome_linear")

Suda 2002 didn't determine the lengths of the sex chromosomes. I want to create a formula to estimate the telomere lengths of the sex chromosomes. Here's how I did a linear regression in R:

```r
# this is a linear regression of the variables
lm(telomere_lengths ~ chromosome_lengths)
```

```sh
Call:
lm(formula = telomere_lengths ~ chromosome_lengths)

Coefficients:
       (Intercept)  chromosome_lengths  
           1920.25               14.93  
```

Here's that same initial plot, but I have included the line of fit:

```{r}
# here I'm using abline to add the line of regression
plot(chromosome_lengths, telomere_lengths, main ="Telomere Lengths are Linearly Related to Chromosome Lengths", xlab="Chromosome Lengths (Mb)", ylab="Telomere Lenghts (bp)")
abline(1920.25, 14.93)
```

![line_of_fit_chromosome_telomere](/Assets/line_of_fit_chromosome_telomere.jpg "line_of_fit_chromosome_telomere")

```r
# I'm assuming that this relationship will hold true for chromosomes X + Y
# AND that p/q is 1 in telomerase (Perrem 2001) 
# Morton 1991 Parameters has lengths of each chromosome, BUT, keep in mind that the chromosome Mb are slightly different from Suda 2002
chromosome_X_length <- 14.93*164 + 1920.25
chromosome_Y_length <- 14.93*59 + 1920.25
# Suda 2002 used a FACS technique that couldn't differentiate these chromosomes 1=2, 9=10=11=12
# I have assumed that their lengths are close enough (cause they're chromosomes are close in length)
chromosome_1_length <- telomere_lengths[1]
chromosome_2_length <- telomere_lengths[1]
chromosome_9_length <- telomere_lengths[8]
chromosome_10_length <- telomere_lengths[8]
chromosome_11_length <-telomere_lengths[8]
chromosome_12_length <- telomere_lengths[8]
# this could've been done more delicately, lol
# here's the complete list of lengths for autosomal chromosomes 1-22 and sex chromosomes x and y
chrom_1 <- chromosome_1_length
chrom_2 <- chromosome_2_length
chrom_3 <- telomere_lengths[2]
chrom_4 <- telomere_lengths[3]
chrom_5 <- telomere_lengths[4]
chrom_6 <- telomere_lengths[5]
chrom_7 <- telomere_lengths[6]
chrom_8 <- telomere_lengths[7]
chrom_9 <- chromosome_9_length
chrom_10 <- chromosome_10_length
chrom_11 <- chromosome_11_length
chrom_12 <- chromosome_12_length
chrom_13 <- telomere_lengths[9]
chrom_14 <- telomere_lengths[10]
chrom_15 <- telomere_lengths[11]
chrom_16 <- telomere_lengths[12]
chrom_17 <- telomere_lengths[13]
chrom_18 <- telomere_lengths[14]
chrom_19 <- telomere_lengths[15]
chrom_20 <- telomere_lengths[16]
chrom_21 <- telomere_lengths[17]
chrom_22 <- telomere_lengths[18]
chrom_X <- chromosome_X_length
chrom_Y <- chromosome_Y_length
# I should've just gone with chromosome X OR Y cause I'm only using one representative of the chromosomes 1-22
# BUT, I like how the bar graph looks this way :)
telomerase_chromosome_lengths <- c(chrom_1, chrom_2, chrom_3, chrom_4, chrom_5, chrom_6, chrom_7, chrom_8, chrom_9, chrom_10, chrom_11, chrom_12,chrom_13,chrom_14,chrom_15,chrom_16,chrom_17,chrom_18,chrom_19,chrom_20, chrom_21, chrom_22, chrom_X, chrom_Y)
# here are the labels for the chromosomes that I'm using 
telomerase_chromosome_names <- c((1:22), "X", "Y")
# and finally the barplot
barplot(telomerase_chromosome_lengths, names.arg=telomerase_chromosome_names, las=2, main ="GM130B Telomere Lengths", xlab="Chromosome", ylab="Telomere Lenghts (bp)", ylim=c(0,6000))
```

![GM130B_Telomere_Lengths](/Assets/GM130B_Telomere_Lengths.jpg "GM130B_Telomere_Lengths")

Suda 2002 estimated that the GM130B cell line could divide somewhere between 15 and 42 times. Here's how they arrived at that conclusion:

1) Identify the shortest and average telomere lengths
	* chromosome 21 is 2735 bp
2) Identify the average telomere lengths
	* They claim it's 3998 bp (It's actually 3770 bp ... check their math)
3) Subtract the 2000 bp subtelomeric region
	* chromosome 21 is down to 735 bp and the average is down to 1998 bp
4) Divide those bp numbers by the bp lost/division
	* chromosome 21 735/48 = 15.3 PD & 1998/48=41.6 PD

It's well-established nowadays that the shortest telomere is what determines when cellular senescence will start. Here's a model of the short telomere cellular senescence model from Suda 2002:

```r
# this is the condition that, while TRUE, runs the model
above_zero <- TRUE
division_number <- 0
# assuming subtelomere is 2000 bp
telomerase_chromosome_lengths_play <- telomerase_chromosome_lengths - 2000
while(above_zero) {
  barplot(telomerase_chromosome_lengths_play, names.arg=telomerase_chromosome_names, las=2, main = paste("GM130B Senescence After ", division_number, " Divisions"), xlab="Chromosome", ylab="Telomere Lenghts (bp)", ylim=c(0,3000))
  telomerase_chromosome_lengths_play <- telomerase_chromosome_lengths_play - 48
  division_number <- division_number + 1
  if(sum((telomerase_chromosome_lengths_play) < 0) > 0) {
    above_zero <- FALSE
  }
}
```

Here's that final senescent state:

![GM130B_Senesces_After_15PD](/Assets/GM130B_Senesces_After_15PD.jpg "GM130B_Senesces_After_15PD")

#### Modeling WT Telomere Shortening
The Suda 2002 is a very interesting case with an immortalized cell. This is my attempt to model the initial telomere length state for a human at birth and follow it through until senescence. Note that I'm using the Harley 2008 5 kb red line of senescence and the black 3kb line of crisis.

```r
# I wanna make sure that I don't change the length list that I created earlier
telomerase_chromosome_lengths_play <- telomerase_chromosome_lengths
# the range of telomere lengths from GM130B is 2946
max(telomerase_chromosome_lengths_play) - min(telomerase_chromosome_lengths_play)
# scaling that range up to 10000 bp with multiplication yields a range of 5185
telomerase_chromosome_lengths_START <- telomerase_chromosome_lengths * 10000/max(telomerase_chromosome_lengths)
max(telomerase_chromosome_lengths_START) - min(telomerase_chromosome_lengths_START)
# 5185 bp seemed a bit high for the range (just intuition), so I used an additive instead
telomerase_chromosome_lengths_START <- telomerase_chromosome_lengths + 10000 - max(telomerase_chromosome_lengths)
max(telomerase_chromosome_lengths_START) - min(telomerase_chromosome_lengths_START)

# How long till senescence if starting max is 10K AND maintaining telomere length range of GM130B?
above_five_thousand <- TRUE
division_number <- 0
telomerase_chromosome_lengths_play <- telomerase_chromosome_lengths_START
while(above_five_thousand) {
  #print(telomerase_chromosome_lengths_play)
  #print(division_number)
  barplot(telomerase_chromosome_lengths_play, names.arg=telomerase_chromosome_names, las=2, main = paste("hCell Telomere Shortening After ", division_number, " Divisions"), xlab="Chromosome", ylab="Telomere Lenghts (bp)", ylim=c(0,10000))
  abline(h=5000, col="red", lwd=3)
  abline(h=3000, col="black", lwd=3)
  division_number <- division_number + 1
  telomerase_chromosome_lengths_play <- telomerase_chromosome_lengths_play - 48
  if(sum(telomerase_chromosome_lengths_play < 5000) > 0) {
    above_five_thousand <- FALSE
  }
}
```

Here's the initial state:
![hCell_Shortening_After_0_Divisions](/Assets/hCell_Shortening_After_0_Divisions.jpg "hCell_Shortening_After_0_Divisions")

This model divided 42 times before senescing.

![hCell_Shortening_After_42_Divisions](/Assets/hCell_Shortening_After_42_Divisions.jpg "hCell_Shortening_After_42_Divisions")



# Telomere Maintenance Mechanism (TMM) Selection
The selection of TMM seems to mainly depend upon telomerase chromatin compaction and mutations in ATRX & p53 (Gocha 2013). 
 
![ALT_TEL_Permissive_Mutations](/Assets/ALT_TEL_Permissive_Mutations.jpg "ALT_TEL_Permissive_Mutations")

(Gocha 2013)

# Telomerase Extends Telomeres
Telomerase adds 5'-GGTTAG-3' (Harley 2008). Telomerase extends the shortest telomeres first (Harley 2008, Cristofari 2006). The G-rich strand is 5'-GGTTAG-3' and the C-rich strand is 3'-CCAATC-5'. Telomerase adds 5'-GGTTAG-3' (Harley 2008). The telomerase enzyme TERT uses the telomerase RNA 3'-CAAUCCCAAUC-5' as a template for the extension (Gavory 2002). There is a G-rich single stranded telomeric overhang of 130-210 nucleotides (Cesare 2010). telomerase, which is a reverse transcriptase that adds repetitive telomeric DNA (TTAGGG)n to the ends of the chromosomes (Allsopp 2001). 

#### Model of Different Telomerase Levels
Telomerase overexpression might have an upper limit of 0.8 kb/division ... I'm not certain, but that was reported in Cristofari 2006. I might want to create a future update that limits to 800 bases/division just to keep cellular resources in mind. I don't know how model the telomerase activity in Python with a great deal of biological accuracy ... 

1. How much RNA template is available?
2. How long will a telomerase enzyme be active?
3. How do you model 3-D interactions?
etc.,

so I've decided to simplify the telomerase model to adding back the 48 bp lost HALF of the time. Harley 2008 box 1b (below) was the inspiration of this idea, but keep in mind that things are WAY more complicated than this. Yes, it's a spherical cow kind of model, but what could I do better? Seriously, message me if you've got an idea, cause I'll try it out :) 

![Harley_2008_box1b](/Assets/Harley_2008_box1b.jpg "Harley_2008_box1b")

```{r}
# How long till senescence if starting max is 10K AND maintaining telomere lengths of GM130B?
# AND telomerase is periodically turned on? AND senescence checkpoints work.
# Telomerase extends the shortest telomeres first (Harley 2008, Cristofari 2006)
# I'm not worried about that when telomere shortening overpowers telomere lengthening
above_five_thousand <- TRUE
division_number <- 0
telomerase_chromosome_lengths_play <- telomerase_chromosome_lengths_START
while(above_five_thousand) {
  #print(telomerase_chromosome_lengths_play)
  #print(division_number)
  barplot(telomerase_chromosome_lengths_play, names.arg=telomerase_chromosome_names, las=2, main = paste("Trasient Telomerase Telomere Shortening After ", division_number, " Divisions"), xlab="Chromosome", ylab="Telomere Lenghts (bp)", ylim=c(0,10000))
  abline(h=5000, col="red", lwd=3)
  abline(h=3000, col="black", lwd=3)
  division_number <- division_number + 1
  telomerase_chromosome_lengths_play <- telomerase_chromosome_lengths_play - 48
  if(division_number %% 2){
    telomerase_chromosome_lengths_play <- telomerase_chromosome_lengths_play + 48
  }
  if(sum(telomerase_chromosome_lengths_play < 5000) > 0) {
    above_five_thousand <- FALSE
  }
  # this indicates immortality and I don't want a runaway while loop
  if(division_number > 100) {
    break
  }
}
```

![Transient_Telomerase_HSC](/Assets/Transient_Telomerase_HSC.jpg "Transient_Telomerase_HSC")


This is when telomerase is equal to the telomere shortening.

```{r}
# what about a cancer with telomerase always turned on?
# telomerase + cells still have the trimming mechanism turned on
# Pickett 2009  HT1090 +TEL stuck at 7.5 kb, HeLa stuck at 4.5 kb
# overexpression of hTER HeLa 8 kb, HT1080 hTR increasingly heteogenous 9.4 kb
# I could make a complicated function to make sure things stay matched, but you get the idea ...
# I'll just make shortening equal to lengthening
above_five_thousand <- TRUE
division_number <- 0
telomerase_chromosome_lengths_play <- telomerase_chromosome_lengths_START
while(above_five_thousand) {
  #print(telomerase_chromosome_lengths_play)
  #print(division_number)
  barplot(telomerase_chromosome_lengths_play, names.arg=telomerase_chromosome_names, las=2, main = paste("Immortalized Telomerase Telomere Shortening After ", division_number, " Divisions"), xlab="Chromosome", ylab="Telomere Lenghts (bp)", ylim=c(0,10000))
  abline(h=5000, col="red", lwd=3)
  abline(h=3000, col="black", lwd=3)
  division_number <- division_number + 1
  telomerase_chromosome_lengths_play <- telomerase_chromosome_lengths_play - 48
  telomerase_chromosome_lengths_play <- telomerase_chromosome_lengths_play + 48
  if(sum(telomerase_chromosome_lengths_play < 5000) > 0) {
    above_five_thousand <- FALSE
  }
  # this indicates immortality and I don't want a runaway while loop
  if(division_number > 300) {
    break
  }
}
```

Note that the code breaks at 300 population doublings because it'd go on forever otherwise. Telomerase does seem to have an upper limit for the telomere lengths that are possible. Picket 2009 found that HeLa and HT1080 stayed at 4.5 and 7 kbp under normal telomerase conditions, BUT that they went up to 8 and 9.4 kbp with hTR overexpression. There appears to be a natural telomere trimming mechanism under normal conditions (Picket 2009).

![Immortalized_Telomerase](/Assets/Immortalized_Telomerase.jpg "Immortalized_Telomerase")

#### Modeling the Single Stranded G-rich Tail
Another thing that I left out of earlier models was the G-rich and C-rich strands. I've just been describing the telomeres as having a length. They have sequence too! The G-rich strand is 5'-GGTTAG-3' and the C-rich strand is 3'-CCAATC-5'. Telomerase adds 5'-GGTTAG-3' (Harley 2008). The telomerase enzyme TERT uses the telomerase RNA 3'-CAAUCCCAAUC-5' as a template for the extension (Gavory 2002). There is a G-rich single stranded telomeric overhang of 130-210 nucleotides (Cesare 2010).

I made a simple model of the G'rich overhang in R.

```r
setwd("/media/david/Linux/ALT_Introns_Exons_and_Promoters/Telomere_Math_Models")
library(msa)
library(seqinr)
G_overhang_sequence <- "GGTTAG"
C_overhang_sequence <- "CCAATC"
reverse_complement_C_overhang_sequence <- "CTAACC"
# 22*6 = 132 (min single stranded overhang is 130)
# assuming 10 kbp telomere lengths, so 1666 * 6 9996 is for C strand
# then 1666+22 is for G strand
full_G_overhang_sequence <- paste(replicate(1688, G_overhang_sequence), collapse = "")
full_C_overhang_sequence <- paste(replicate(1666, reverse_complement_C_overhang_sequence), collapse = "")
write.fasta(sequences = full_G_overhang_sequence, names = "full_G_overhang_sequence", file.out = "full_G_overhang_sequence.fasta", open = "w", nbchar = 70, as.string = FALSE)
write.fasta(sequences = full_C_overhang_sequence, names = "full_C_overhang_sequence", file.out = "full_C_overhang_sequence.fasta", open = "w", nbchar = 70, as.string = FALSE)
```

I can't figure out how to make a pretty DNA complementary base pairing figure ... any ideas? Please let me know :) I made a simple, shortened, text version:

![G_overhang](/Assets/G_overhang.jpg "G_overhang")

NOTE: I NEED TO MAKE SURE I'VE GOT THE RIGHT NUMBER FOR TELOMERE VS. SUBTELOMERE ... CHECK RESTRICTION ENZYME DATA AND READ THESE TWO PAPERS: 
Riethman 2004 Mapping and Initial Analysis of Human Subtelomeric Sequence Assemblies
Riethman 2009 Human subtelomeric copy number variations
Mefford 2002 The complex structure and dynamic evolution of human subtelomeres
Suda 2002 Interchromosomal Telomere Length Variation
Graakjaer 2003 The pattern of chromosome-specific variations in telomere length in humans is determined by inherited, telomere-near factors and is maintained throughout life

#### Model of Inhibiting Telomerase 
When DNA is getting copied, small end pieces aren't able to be fully copied. This is known as the End Replication Problem and is a result of Okazaki fragments incompletely covering the end of the chromosome. This is an issue for cell types that need to divide a lot, like Hematopoietic stem cells (HSC). HSCs have plenty of division to do so they can keep up with all the differentiation to make blood cells. This is why they express telomerase, which is a reverse transcriptase that adds repetitive telomeric DNA (TTAGGG)n to the ends of the chromosomes (Allsopp 2001). Did you notice that I'm citing a paper from the Weissman lab? ;)  
 
![telomeres-shorten-with-age](/Assets/telomeres-shorten-with-age.jpg "telomeres-shorten-with-age")

(Finkel 2007)

Cells that divide a lot will senesce in the absence of an active telomere maintenance mechanism (Shay 2012). That's why Telomerase is great for stem cells and other cell types that need to divide a lot. HOWEVER, it's also part of how the majority of cancers become immortal (Cesare 2010). This is why several companies are working on telomere-based anti-cancer therapies:

* Telomerase enzyme inhibitor: 
	* GRN163L (Geron)
* Telomerase active immunotherapy:	
	* GRNVAC1 (Geron)
	* GV1001 (Pharmexa)
	* P540-548 (Gemvax)
	* Vx01 (Vaxon Biotech)
	* TLI (Cosmo Bioscience) 

(Shay 2012, Harley 2008)

This model is the same as the telomere shortening + telomerase model AND you can choose a cell division number for telomerase to be inhibited at. There are two scenarios that I want to model:

1) cancers with relatively short telomeres (these will divide out without telomerase)
2) cancers with very long telomeres (these might be able to keep dividing in the absence of telomerase)
3) cancers with telomeres that are equilavent to stem cell telomeres (anti-telomerase therapy will exhaust the stem cell pools before the cancer)

I bring up point 3 because one problem with telomerase inhibitors is that they have been reported to cause problems with blood stem cells (Hu 2017). This update only needed three lines of code. I added a varible for the division number for telomerase to stop at and an if condition that sets telomerase activity to 0 AFTER that division number has been reached. 

```r
I NEED TO FINISH MAKING THIS
```

# Alternative Lengthening of Telomeres (ALT) Extends Telomeres
Stem cell telomerase isn't the only problem with the telomerase inhibition approach. Approximately 10-15% of cancers use the Alterantive Lengthening of Telomeres (ALT) to extend telomeres, some cancers won't be treated with these anti-telomerase therapies. But, it's way worse of a problem than that! Tumors have been reported to use both ALT and TEL simultaneously (Gocha 2013) AND in vitro inhibition of telomerase selects for ALT activity (Sahin 2012). 

![TEL_ALT_Reversible](/Assets/TEL_ALT_Reversible.jpg "TEL_ALT_Reversible")

(Shay 2012)

On an interesting note, there is a high frequency of cancers with a Mesenchymal Stem Cell origin using the ALT pathway (77% malignant fibrous histocytomas, 47-66% of osteosarcomas (Lafferty-Whyte 2009), 21.4% of liposarcomas (Venturini APB ALT LIPOSARCOMA 2008). This may be a result of the telomerase promoter being repressed with chromatin compaction (Atkinson 2005). Telomerase inhibition MIGHT not be alone in potentially causing problems for stem cells. There is some evidence to suggest that stem cells use ALT (Kalmbach 2014, Huang 2014). To make things worse, some cancers appear to extensively divide without a telomere maintenance mechanism (Dagg 2017).  

![stem_cell_ALT.jpg](/Assets/stem_cell_ALT.jpg "stem_cell_ALT.jpg")

(Kalmbach 2014)

#### ALT Telomeres are Long and Heterogenous
ALT telomeres have long, heterogenous telomeres. They can range from less than 1 kbp (Rogan 1995) all the way up to 50 kbp (Bryan 1995). There is a rapid addition of telomeric sequences to short telomeres AND a rapid deletion of telomeric sequences from the long telomeres. This suggests recombinogenic behavior (Murnane 1994). ALT cells seem to be extending their telomeres through this process, but the mechanism isn't completely understood yet (Cesare 2010). There are common losses and duplications of chromosomes (Sakellariou 2013). The ratio of p/q telomeric arms range from 10-0.1.

I think there are still more details that I should take into account for this idea, but I'm not sure when I'll have the time for that. The simplified model I have in mind will:

1) range in initialized telomere length from 500 bp to 50,000 bp
2) have a p/q telomeric arm ratio betwen 10-0.1
3) take the 92 different telomere endings into account

Forgive me for not representing actual ALT activity. It was pretty tough just getting the code to inititate ALT telomeres!

```r
# human cells have 92 telomeres at the G0 phase of the cell cycle.
# 2 pairs of chromosomes * 23 chromosomes * 2 (p and q arms) = 92 telomeres
# here are all of the names that I wrote out by hand ... boring! I should've written a function in Python to do this, lol!
autosome_pairs_sex_chromosomes_and_arms <- c("c1a1p", "c1a1q","c2a1p","c2a1q","c3a1p","c3a1q","c4a1p","c4a1q","c5a1p","c5a1q","c6a1p","c6a1q","c7a1p","c7a1q","c8a1p","c8a1q","c9a1p","c9a1q","c10a1p","c10a1q","c11a1p","c11a1q","c12a1p","c12a1q","c13a1p","c13a1q","c14a1p","c14a1q","c15a1p","c15a1q","c16a1p","c16a1q","c17a1p","c17a1q","c18a1p","c18a1q","c19a1p","c19a1q","c20a1p","c20a1q","c21a1p","c21a1q","c22a1p","c22a1q","cXp","cXq","c1a2p","c1a2q","c2a2p","c2a2q","c3a2p","c3a2q","c4a2p","c4a2q","c5a2p","c5a2q","c6a2p","c6a2q","c7a2q", "c8a2p","c8a2q","c7a2p","c9a2q","c9a2p","c102q","c102p","c11a2q","c11a2p","c12a2p","c12a2q", "c13a2p","c13a2q", "c142p","c142q", "c15a2p", "c15a2q","c16a2p","c16a2q","c172p","c172q","c18a2p","c18a2q","c19a2p", "c19a2q","c20a2p","c20a2q", "c21a2p","c21a2q","c22a2p","c22a2q","cYp","cYq")
length(autosome_pairs_sex_chromosomes_and_arms)
# it creates 92 telomeres
ALT_telomere_lengths <- vector("list", 92)
# these are placeholders for values that will be used throughout the program
current_telomere_length_total <- 0
current_telomere_length_Q <- 0
current_telomere_length_P <- 0
# i keeps track of the current interation. Note that each run has two iterations of adding to the tel length list
i <- 1
generate_ALT_telomeres <- function() {
  # runs for 46 items (cause 92 total and creating p and q w/ each loop)
  for(length in 1:46){
    # determine the total telomere length between 500 and 50,000 bp
    current_telomere_length <- runif(1, min=500, max=50000)
    # determine the p/q ratio between 0.1 and 10
    P_to_Q_ratio <- runif(1, 0.1, 10)
    # get the q arm length. This is some cool math! BUT, it might not be obvious
    # there are two formulas here
    # 1) p arm + q arm = total length
    # 2) p arm / q arm = ratio => p arm = ratio * q arm
    # substituting 2 into 1: ratio * q arm + q arm = total length
    current_telomere_length_Q <- current_telomere_length/(P_to_Q_ratio+1)
    # get the p arm length
    current_telomere_length_P <- current_telomere_length - current_telomere_length_Q
    # get the sum (this number is already known. I'm using it for error checking)
    current_sum_p_q <- current_telomere_length_Q+current_telomere_length_P
    # get the p/q ratio (this number is already known. I'm using it for error checking)
    determine_P_to_Q_ratio <- current_telomere_length_P / current_telomere_length_Q
    # print out the math to check by hand
    print(paste("current total length is ", current_telomere_length, "Q is ", current_telomere_length_Q, "P is ", current_telomere_length_P, " sum is ", current_sum_p_q))
    print(paste("current P_to_Q_ratio is ", P_to_Q_ratio, "Q is ", current_telomere_length_Q, "P is ", current_telomere_length_P, " P_to_Q_ratio is ", determine_P_to_Q_ratio))
    print("")
    # store current p and q arms
    ALT_telomere_lengths[i] <- current_telomere_length_P 
    i <- i + 1
    ALT_telomere_lengths[i] <- current_telomere_length_Q
    i <- i + 1
  }
  return(ALT_telomere_lengths)
}
ALT_telomere_lengths <- generate_ALT_telomeres()
```

Here's the printed output. Note that I made it caluclate and return all of the values, so I could check them by eye. The math in the code isn't immediately obvious and it's always important to check your code and your math for errors! Here are the first couple of lines of code. Give them a math check!

```sh
[1] "current total length is  29483.7727698032 Q is  5957.23538429581 P is  23526.5373855074  sum is  29483.7727698032"
[1] "current P_to_Q_ratio is  3.94923750159796 Q is  5957.23538429581 P is  23526.5373855074  P_to_Q_ratio is  3.94923750159796"
[1] ""
[1] "current total length is  5732.48566687107 Q is  1917.02364226407 P is  3815.462024607  sum is  5732.48566687107"
[1] "current P_to_Q_ratio is  1.99030514829792 Q is  1917.02364226407 P is  3815.462024607  P_to_Q_ratio is  1.99030514829792"
[1] ""
[1] "current total length is  40303.4766089404 Q is  6558.95629193797 P is  33744.5203170024  sum is  40303.4766089404"
[1] "current P_to_Q_ratio is  5.14480030282866 Q is  6558.95629193797 P is  33744.5203170024  P_to_Q_ratio is  5.14480030282866"
[1] ""
[1] "current total length is  32052.6402043179 Q is  9641.5347792759 P is  22411.105425042  sum is  32052.6402043179"
[1] "current P_to_Q_ratio is  2.32443339552265 Q is  9641.5347792759 P is  22411.105425042  P_to_Q_ratio is  2.32443339552265"
[1] ""
[1] "current total length is  26283.9130933862 Q is  4428.17189769478 P is  21855.7411956914  sum is  26283.9130933862"
[1] "current P_to_Q_ratio is  4.93561264120508 Q is  4428.17189769478 P is  21855.7411956914  P_to_Q_ratio is  4.93561264120508"
```

Here are those telomeres plotted as a bar graph:

```r
# This is a disgusting mess of overhwelming data
barplot(as.numeric(ALT_telomere_lengths), names.arg=autosome_pairs_sex_chromosomes_and_arms, las=2, main="ALT+ Telomere Lengths for Chromosomes", xlab="Chromosome", ylab="Telomere Lengths (bp)", ylim=c(0,50000))
abline(h=5000, col="red", lwd=3)
abline(h=3000, col="black", lwd=3)
```

Gosh! This looks terrible! The axes labels are covering over the telomere lengths AND the chromosome names. I'll need to tidy up ... actually, I'm just going to look at telomeres 1-46

![ALT_ugly_telomeres](/Assets/ALT_ugly_telomeres.jpg "ALT_ugly_telomeres")

```r
barplot(as.numeric(ALT_telomere_lengths[1:46]), names.arg=autosome_pairs_sex_chromosomes_and_arms[1:46], las=2, main="ALT+ Telomeres 1-46", ylab="Telomere Lengths (bp)", ylim=c(0,50000))
abline(h=5000, col="red", lwd=3)
abline(h=3000, col="black", lwd=3)
```

![ALT_telomeres_1_46](/Assets/ALT_telomeres_1_46.jpg "ALT_telomeres_1_46")

```r
barplot(as.numeric(ALT_telomere_lengths[47:92]), names.arg=autosome_pairs_sex_chromosomes_and_arms[47:92], las=2, main="ALT+ Telomeres 47-92", ylab="Telomere Lengths (bp)", ylim=c(0,50000))
abline(h=5000, col="red", lwd=3)
abline(h=3000, col="black", lwd=3)
```

![ALT_telomeres_47_92](/Assets/ALT_telomeres_47_92.jpg "ALT_telomeres_47_92")


#### ALT Telomeres Have C-rich Overhangs
There is some degree of 5' C-rich overhang in ALT cells. 

![C_overhang](/Assets/C_overhang.jpg "C_overhang")

^I need to do a bit more reading to make an accurate model of this. Here are some good papers to check out:

* Oganesian 2011 Mammalian 5 0 C-Rich Telomeric Overhangs Are a Mark of Recombination-Dependent Telomere Maintenance
* Min 2017 Alternative lengthening of telomeres can be maintained by preferential elongation of lagging strands
* Mao 2016 Homologous recombination-dependent repair of telomeric DSBs in proliferating human cells

The altered C-rich overhangs involved in ALT will likely prevent POT-1 from binding effectively to the telomeres. This is a great segway into the ALT literature! For example, POT-1 deficiency creates ALT in C. elegans!

# POT-1 Deficiency Creates ALT+ C. Elegans Strains
Telomeres cap linear chromosomes because DNA polymerase can't completely copy chromosomes. Telomerase adds telomeric repeats to the ends of linear chromosomes with reverse transcription. Most human cancers have long, heterogenous telomeres. Telomere shortening leads to senescence and potentially crisis. Cancer emerges as part of massive cell death and genomic rearrangements after crisis. 10-15% of cancers are estimated to use ALT (Cheng 2012). 

ALT can happen in Caenorhabditis elegans! Mammalian POT1 has homologs in C. elegans as pot-1 (CeOB2) and pot-2 (CeOB1). What's the deal with the reversing of 1 and 2? That's how it's reported in the paper ... it's odd. pot-1 mutant C. elegans have HUGE telomere lengths while pot-2 mutants have normal telomere lengths. The authors of Cheng 2012 created a variety of mutants in C. elegans.  The trt-1 C. elegans mutant has a deletion in telomerase reverse transcriptase. trt-1 & pot-2 absence led to ALT+ Caenorhabditis elegans with normal telomere lengths. trt-1 and pot-1 mutants were found to have long, heterogenous telomere lengths like those seen in human ALT. Here is the survival figure showing that C. elegans can survive in the absence of telomerase reverse transcriptase.

![Celegans_ALT_Generation_Survival](/Assets/Celegans_ALT_Generation_Survival.jpg "Celegans_ALT_Generation_Survival")

(Cheng 2012)

#### Multiple Sequence Alignment of pot-1 Genes
YES, pot-2 was the central point of the paper, but it won't be as fun to play with because it only has one isoform. I picked pot-1 cause there is a lot of cool stuff to play with. There were a lot of workup steps to get all of the sequences ... It would take a long while to review them. Essentially, I looked up the proteins on UniProt and then grabbed the DNA files from NCBI GenBank and WormBase. Check out the Celegans_POT1_ALT folder for the file names of everything. The file containing all the C. elegans genes is Celegans_POT1_genes.fasta. I used the R package "msa" for multiple sequence alignment with this code:

```r
library(msa)
Celegans_POT1_genes <- "/media/david/Linux/Introns_Exons_and_Promoters/Celegans_POT1_ALT/DNA/Celegans_POT1_genes.fasta"
Celegans_POT1_genes_DNA <- readDNAStringSet(Celegans_POT1_genes)
Celegans_POT1_gene_alignment <- msa(Celegans_POT1_genes_DNA)
msaPrettyPrint(Celegans_POT1_gene_alignment, output="pdf", showNames="none",
showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
```

The aligned sequences aren't very pretty ... I decided not to include sequence labels cause it shortened the available nucleotide space for each new line. Here's part of the output for you to get the idea of the work:

![Celegans_POT1_gene_alignment](/Assets/Celegans_POT1_gene_alignment.jpg "Celegans_POT1_gene_alignment")

#### Multiple Sequence Alignment of pot-1 Proteins
I grabbed all the C elegans pot-1 isoform sequences from UniProt. You can check them out in Celegans_POT1_ALT/Protein. The file containing all of the sequences is Celegans_POT1_Proteins.fasta. I aligned all of the proteins with code that is similar to the DNA alignment code:

```r
Celegans_POT1_proteins <- "/media/david/Linux/Introns_Exons_and_Promoters/Celegans_POT1_ALT/Protein/Celegans_POT1_Proteins.fasta"
Celegans_POT1_proteins_AA <- readAAStringSet(Celegans_POT1_proteins)
Celegans_POT1_protein_alignment <- msa(Celegans_POT1_proteins_AA)
Celegans_POT1_protein_alignment
msaPrettyPrint(Celegans_POT1_protein_alignment, output="pdf", showNames="none",
showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
```

This alignment looks great! You can see all of the alignments between the different pot-1 isoforms :)

![Celegans_POT1_protein_alignment](/Assets/Celegans_POT1_protein_alignment.jpg "Celegans_POT1_protein_alignment")

#### Displaying pot-1 Open Reading Frames
I used code from the BioPython Tutorial to identify the C. elegans pot-1 gene Open Reading Frames AND to report the translated proteins! Compare this to the last section to see that the the DNA -> Protein translations all have the correct lengths! The code is from http://biopython.org/DIST/docs/tutorial/Tutorial.html in the section titled "20.1.13. Identifying open reading frames". You should check this code out! It can identify Open Reading Frames in the +/- strand AND in three different reading frames, SO IT DOES ALL 6 FRAMES!!! I re-used the code four times (instead of making a function, haha). Here's part of the code that I used:

```python
#!/usr/bin/env python

from Bio import SeqIO
# record = SeqIO.read("NC_005816.fna", "fasta")
file_0 = "NM_001361730.1_Caenorhabditis_elegans_pot-1_gene_homolog.fasta"
file_1 = "NM_001361731.1_Caenorhabditis_elegans_pot-1_gene_homolog.fasta"
file_2 = "NM_001361732.1_Caenorhabditis_elegans_pot-1_gene_homolog.fasta"
file_3 = "NM_066157.3_pot-1_NCBI_DNA_matchesP42001_Celegans.fasta"

print("")
print("NM_001361730.1_Caenorhabditis_elegans_pot-1_gene_homolog.fasta")
record = SeqIO.read(file_0, "fasta")
table = 11
min_pro_len = 150

for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
    for frame in range(3):
        length = 3 * ((len(record)-frame) // 3) #Multiple of three
        for pro in nuc[frame:frame+length].translate(table).split("*"):
            if len(pro) >= min_pro_len:
                print("%s...%s - length %i, strand %i, frame %i" \
                % (pro[:30], pro[-3:], len(pro), strand, frame))
```

![Displaying_pot-1_Open_Reading_Frames](/Assets/Displaying_pot-1_Open_Reading_Frames.jpg "Displaying_pot-1_Open_Reading_Frames")

#### Discussing C. elegans pot-1 Alternative Splicing
WormBase has three isoforms for pot-1 in C. elegans https://wormbase.org/species/c_elegans/gene/WBGene00015105#0-9g-3

![WormBase_pot-1_Celegans_Isoforms](/Assets/WormBase_pot-1_Celegans_Isoforms.jpg "WormBase_pot-1_Celegans_Isoforms")

Transcript B0280.10a.1 is 1216 nucleotides in length and codes for a 400 amino acid protien. The WormBase curators used RNA-seq data from Boeck 2016 to alter the original WormBase entry to include the published alternate intron splicing. You can see from the table that Exons 1-10 are part of this pot-1 isoform.

![B0280.10a.1_Celegans_pot-1_isoformA_1203NT_400AA](/Assets/B0280.10a.1_Celegans_pot-1_isoformA_1203NT_400AA.jpg "B0280.10a.1_Celegans_pot-1_isoformA_1203NT_400AA")

Transcript B0280.10b.1 is 462 nucleotides in length and it codes for a 153 amino acid protein. You can see from the table that Exons 1-4 are part of this pot-1 isoform. Boeck 2016 goes into more detail about the alternative intron that causes alternative splicing here. 

![B0280.10b.1_Celegans_pot-1_isoformB_462NT_153AA](/Assets/B0280.10b.1_Celegans_pot-1_isoformB_462NT_153AA.jpg "B0280.10b.1_Celegans_pot-1_isoformB_462NT_153AA")

Transcript B0280.10c.1 is 1140 nucleotides long and it codes for a 379 amino acid protein. You can see from the table that Exons 1-10 make it into this protein isoform.

![B0280.10c.1_Celegans_pot-1_isoformC_1140NT_379AA](/Assets/B0280.10c.1_Celegans_pot-1_isoformC_1140NT_379AA.jpg "B0280.10c.1_Celegans_pot-1_isoformC_1140NT_379AA")

# ATRX Exon Deletion is Common in ALT
This project can be found in the Human_ATRX_ALT folder. ATRX gene mutations are found in a range of cancers. 10-15% of cancers are estimated to use ALT. ALT involves homologous recombination-based telomere elongation. Inactivating mutations in either ATRX or DAXX are found in many cancers. Depletion of ATRX seems insufficient to trigger ALT, but it does seem to play a key role in the ALT pathway. The absence of ATRX might lead to the failure of stalled replication forks to get resolved. The required fork restart would require homologus recombination and could jumpstart the ALT pathway (Clynes 2013). ALT involves a template-based lengthening of telomeres with homologous recombination. The genetic and epigenetic changes are not full understood. Lovejoy 2012 reported that ATRX gene mutations are a common feature of ALT. Specifically 19/22 ALT+ cell lines had an issue with the expression of ATRX or DAXX (Lovejoy 2012). See the Lovejoy 2012 supplementary information for the Excel table of Exon deletions in ALT cell lines. 
![ATRX_Prevents_Fork_Collapse](/Assets/ATRX_Prevents_Fork_Collapse.jpg "ATRX_Prevents_Fork_Collapse")

(Clynes 2013)

#### Getting ATRX DNA
Searching Ensembl for human ATRX yielded ATRX-201 and ATRX-202. I picked ATRX-201 cause it has 35 exons (which matches the Lovejoy 2012 paper). It was Ensembly ENST00000373344.10. Ensembl refseq switch to NCBI Reference Sequence yielded NM_000489.5 for the gene. I saved it as NM_000489.5_homo_sapiens_ATRX_Gene.fasta.

#### Removing ATRX Exons 2-29 
See the Lovejoy 2012 supplementary Excel table for a list of commonly missing ATRX Exons. I decided to play with the U2OS variant because that is a cell line that I used to grow :) U2OS is missing ATRX exons 2-29. NCBI says exon 2 is [236:348] and exon 29 is 6542..6719 https://www.ncbi.nlm.nih.gov/nuccore/NM_000489. I used R to remove those exons.

```r
library(seqinr)
WT_hATRX_Gene <- read.fasta("NM_000489.5_homo_sapiens_ATRX_Gene.fasta")
WT_hATRX_Gene_Nucleotides <- WT_hATRX_Gene[[1]]
length(WT_hATRX_Gene_Nucleotides)
typeof(WT_hATRX_Gene_Nucleotides)
WT_hATRX_Gene_Nucleotides
U2OS_hATRX_Gene_Nucleotide_FIRST <- WT_hATRX_Gene_Nucleotides[1:235]
U2OS_hATRX_Gene_Nucleotide_SECOND <- WT_hATRX_Gene_Nucleotides[6720:length(WT_hATRX_Gene_Nucleotides)]
U2OS_ATRX_Characters <- c(U2OS_hATRX_Gene_Nucleotide_FIRST, U2OS_hATRX_Gene_Nucleotide_SECOND)
U2OS_ATRX_DNAstring <- DNAString(paste(toupper(U2OS_ATRX_Characters), collapse = ""))
```

#### Sequence Alignment of WT ATRX to Mutant ATRX
The R msa package can't handle the full length of the ATRX gene, so I shortened it down to 400 nucleotides.
```{r}
# limit it to 400 ... that's more than enough to see exon absence
U2OS_ATRX_DNA_Short <- U2OS_ATRX_DNAstring[1:400]
WT_hATRX_DNA_Short <- toupper(WT_hATRX_Gene_Nucleotides[1:400])
write.fasta(sequences = U2OS_ATRX_DNA_Short, names = "U2OS_ATRX_DNA_Short", file.out = "U2OS_ATRX_DNA_Short.fasta", open = "w", nbchar = 70, as.string = FALSE)
write.fasta(sequences = WT_hATRX_DNA_Short, names = "WT_hATRX_DNA_Short", file.out = "WT_hATRX_DNA_Short.fasta", open = "w", nbchar = 70, as.string = FALSE)
library(msa)
# both_ATRX_Sequences <- read.fasta("WT_and_U2OS_hATRX.fasta")
both_ATRX_Sequences_SHORT <- "both_ATRX_Sequences_SHORT.fasta"
# typeof(both_ATRX_Sequences)
#both_ATRX_DNAStringSet <- readDNAStringSet(both_ATRX_Sequences)
#both_ATRX_Sequences_Alignment <- msa(both_ATRX_DNAStringSet)
both_ATRX_Sequences_SHORT_StringSet <- readDNAStringSet(both_ATRX_Sequences_SHORT)
both_ATRX_Sequences_Alignment_SHORT <- msa(both_ATRX_Sequences_SHORT_StringSet)
both_ATRX_Sequences_Alignment_SHORT
msaPrettyPrint(both_ATRX_Sequences_Alignment_SHORT, output="pdf", showNames="none",
showLogo="none", askForOverwrite=FALSE, verbose=TRUE)
#texi2pdf("both_ATRX_Sequences_Alignment.tex", clean=TRUE)
```
You can see that the sequences are the same until postion 236. That is where the Exon deletion for U2OS starts! 

![ATRX_Exon_Deletion_Alignment](/Assets/ATRX_Exon_Deletion_Alignment.jpg "ATRX_Exon_Deletion_Alignment")

# Variants Repeats are Found in ALT Telomeres
POT1 doesn't appear to play a major role in ALT telomeres. Here's a paper that showed relatively unchanged POT1 protein levels in a pre vs. post-ALT cancer (CITATION NEEDED). What's much more interested is the distribution of variant telomeric repeats found in ALT telomeres (Conomos 2012). Telomerase cells appear to have mostly TTAGGG repeats in their telomeres, but ALT telomeres also seem to have TCAGGG repeats. There is a very interesting model of ALT that invovles these TCAGGG repeats getting bound by nuclear receptors (Conomos 2012).

Here are common telomeric sequences found in TEL+:

```r
# Conomos 2012 Figure S1 C HeLa
# Pie Chart with Percentages
slices <- c(889, 74, (0 + 8 + 94 + 8 + 6 + 7 + 4 + 2 + 471)) 
lbls <- c("TTAGGG", "TCAGGG", "Other")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=rainbow(length(lbls)),
   main="HeLa Average Number of Telomeric Repeats")
```

![HeLa_Variant_Repeats](/Assets/HeLa_Variant_Repeats.jpg "HeLa_Variant_Repeats")

Here are common telomeric sequences found in ALT+:

```r
# Conomos 2012 Figure S1 C 
# Pie Chart with Percentages
slices <- c(10165, 4466, (528 + 332 + 241 + 108 + 76 + 73 + 63 + 31 + 3231)) 
lbls <- c("TTAGGG", "TCAGGG", "Other")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=rainbow(length(lbls)),
   main="WI38-VA13/2RA Average Number of Telomeric Repeats")
```
![WI38_VA13_2RA_Variant_Repeat](/Assets/WI38_VA13_2RA_Variant_Repeat.jpg "WI38_VA13_2RA_Variant_Repeat")

This is the nuclear receptor recruitment model first presented in Conomos 2012:

![Conomos_2012_Variant_Repeat_ALT_Model](/Assets/Conomos_2012_Variant_Repeat_ALT_Model.jpg "Conomos_2012_Variant_Repeat_ALT_Model")

# TERT Promoter Compaction is Found in ALT
ALT cells commonly have long, heterogenous telomere lengths (Kumakura 2005). Mouse embryonic stem cells deficient for DNMT have HUGE telomeres. Under normal conditions, mouse subtelomeres are heavily methylated, BUT that is not the case in mESC deficient for DNMT. The lack of DNMT increased the rate of telomeric sister chromatid exchanges (T-SCE), and ALT-associated Promyelocytic Nuclear Bodies (APBs). T-SCE and APBs are both common features of ALT activity. The authors concluded that the increased telomeric recombination MIGHT lead to telomere length changes, BUT they do not exclude the involvement of telomerase in the weirdly long telomeres that were seen (Gonzalo 2006). Luckily, I found these two other papers that go into more detail about TERT chromatin compaction in ALT!

Atkinson 2005 found that chromatin modifications of hTR and hTERT promoters were commonly found in ALT activity. Treatment of ALT+ cells with 5-AZC or Trichostatin A lead to chromatin remodeling of hTR and hTERT. This induced telomerase expression. Interestingly enough they found that mehtylated Lys20 Histon H4 was not associated with gene expression, BUT does seem to be ALT specific (Atkinson 2005). This might be a new marker of ALT activity! Acetylation of H3K9 and methylation of H3K4 is known to be associated with an open chromatin conformation. In Kumakura 2005, the authors found that ALT+ cells had H3K9 methylation and low levels of H3K4 methylation and H3K9+H3K14 acetylations. The ratio of H3K9 methylation / H3K4 methylation was different across ALT+ and TEL+ cell lines. They found that treating ALT+ cells with TSC or 5-AZC caused a reversion from complete to partial methylation of the CpG islands on the hTERT promoter. They switched an E6CL TEL+ line to TEL- and it was able to grow for well over 240 population doublings (Kumakura 2005). That's some ALT activity right there!
![H3K9_H3K4_Methylation_Ratio](/Assets/H3K9_H3K4_Methylation_Ratio.jpg "H3K9_H3K4_Methylation_Ratio")

(Kumakura 2005)

#### Getting The hTERT Sequence
The UniProt entry O14746 is for hTERT https://www.uniprot.org/uniprot/O14746. Following GeneID 7015 gets TERT telomerase reverse transcriptase for humans https://www.ncbi.nlm.nih.gov/gene/7015. Note that the reverse arrow on TERT indicates that the sequence is on the reverse strand (this will become important later). I downloaded the FASTA as "NC_000005.10_hChrom5_TERT_CpG_Start.fasta". A quick text search shows that the start codon is at position 59. Take care with this sequence cause it's the reverse complement of the actual sequence!

![hTERT_NCBI_Reverse_Strand](/Assets/hTERT_NCBI_Reverse_Strand.jpg "hTERT_NCBI_Reverse_Strand")

#### Obtaining hTERT WITH the CpG Island Region
Stay with me ... we're about to dig a bit into the literature! The hTERT sequence that I grabbed from NCBI DOES NOT contain the CpG island for hTERT. It doesn't even contain the normal promoter region for hTERT! Cong 1999 reports that the core hTERT promoter region is from -330 to +361 bp of the ATG start codon. HOWEVER, Kumakura 2005 found the hTERT CpG island to be from 654 bp upstream to 510 bp downstream of the ATG start codon, so this is the actual region that I need to grab! 

Grabbing the hTERT FASTA sequence from https://www.ncbi.nlm.nih.gov/gene/7015 INITIALLY is from: 1253167 to: 1295047. Checking Cong 1999 and The FASTA file, I can see that the hTERT start codon, AND a bit more of that region, of "ATGCCGCGCGCT" is at the end of the first FASTA line, which is 59 in from the left (59 is A of ATG). CpG is 654 bp upstream of the transcriptoin start site, SO going 595 back from current start site 1253167-595 = 1252572. 1252572 should be the start of the CpG island, RIGHT?!? WRONG!!! ... Why aren't I getting any more nucleotides before the current start read?!?!? OH!!! I'm looking at the reverse complement, lol! ;) 

It should be  1295047 + 595 = 1295642 YESSSSSSS, that's right :) Now I have more at the beginniig, so "atgccgcgcgctccccgct" is fully searchable! This is the new range:
https://www.ncbi.nlm.nih.gov/nuccore/NC_000005.10?report=fasta&from=1253167&to=1295642&strand=true. I saved the FASTA as NC_000005.10_hChrom5_TERT_CpG_Start.fasta. 

#### Analyzing the Alleged CpG Promoter Region
Dessain 2000 reported that the hTERT CpG island has a GC content of 74% and a CG:GC ratio of 0.87. Is that what I get for the same region?!? I wrote code in R to get the GC content and CG:GC ratio of the hTERT CpG promoter region that I identified. I didn't comment my code ... I am sorry. Note that I picked i=1164 cause Kumakura 2005 says that the hTERT CpG island to be from 654 bp upstream to 510 bp downstream of the ATG start codon, which is 654 + 510 = 1164 :) Here's R code that I didn't bother commenting :( 

```r
dna_file <- read.fasta("/media/david/Linux/Introns_Exons_and_Promoters/hTERT_CpG_Island/NC_000005.10_hChrom5_TERT_CpG_Start.fasta")
individual_characters <- dna_file[[1]]
i <- 1
CG_count <- 0
GC_count <- 0
g_count <- 0
c_count <- 0
for (letter in individual_characters) {
  #print(letter)
  i <- i + 1
  if (letter == "g") {
    g_count <- g_count + 1
  }
  if (letter == "c") {
    c_count <- c_count + 1
  }
  if ((letter == "c") & (last_letter == "g")) {
    GC_count <- GC_count + 1
  }
  if ((letter == "g") & (last_letter == "c")) {
    CG_count <- CG_count + 1
  }
  last_letter <- letter
  if (i >= 1164) {
    break
  }
}
print(100*(c_count+g_count)/1164)
print(GC_count)
print(CG_count)
print(CG_count/GC_count)
```
The lazily unlabeled output is:

[1] 76.03093
[1] 167
[1] 141
[1] 0.8443114

Cool-ness! My GC content is 76 % and the ratio of CpG/GpC is 0.84. Recall that Dessain 2000 reported that the hTERT CpG island has a GC content of 74% and a CG:GC ratio of 0.87. I'm 2 % off of the GC content and 0.03 off of the CpG/GpC. THAT'S PRETTY GOOD FOR REPLICATING DATA FROM A PAPER THAT IS ALMOST TWO DECADES OLD :D But, what if that was just random luck? I re-ran that same R code on the region AFTER the CpG island and I got wildly different data. Here's that code (yes, I should've used a function; I am lazy, lol):

```r
dna_file <- read.fasta("/media/david/Linux/Introns_Exons_and_Promoters/hTERT_CpG_Island/NC_000005.10_hChrom5_TERT_CpG_Start.fasta")
individual_characters <- dna_file[[1]]
#individual_characters[5]
i <- 1164
CG_count <- 0
GC_count <- 0
g_count <- 0
c_count <- 0
for (letter in individual_characters) {
  if (i <1164) {
    next
  }
  #print(letter)
  i <- i + 1
  if (letter == "g") {
    g_count <- g_count + 1
  }
  if (letter == "c") {
    c_count <- c_count + 1
  }
  if ((letter == "c") & (last_letter == "g")) {
    GC_count <- GC_count + 1
  }
  if ((letter == "g") & (last_letter == "c")) {
    CG_count <- CG_count + 1
  }
  last_letter <- letter
  
  
  if (i >= 42476) {
    break
  }
}
print(100*(c_count+g_count)/41312)
print(GC_count)
print(CG_count)
print(CG_count/GC_count)
```

The lazily unlabeled output is:

[1] 58.55926
[1] 3047
[1] 1654
[1] 0.542829

My GC content is 58.5 % and the ratio of CpG/GpC is 0.54. Recall that Dessain 2000 reported that the hTERT CpG island has a GC content of 74% and a CG:GC ratio of 0.87. THE REGION THAT IS NOT A CpG ISLAND IS 15.5% off of the GC content and 0.33 off of the CpG/GpC. I could dig into this with more statistical rigor, but I think you get the idea. I'M EXCITED!!! This was a really cool biological programming exercise!!! :D

# STN1 Mutation Triggers ALT in Yeast
ADDED STUFF
Counter 1996 The roles of telomeres and telomerase in cell life span
Yeast also can undergo a cellular catastrophe
analogous to the crisis induced in transformed human cells by the loss of telomeric DNA. Disruption
of any one of the S. ceret'isiae genes EST1, KEMI,
or TCL1 (TER1 in K. lactus) results in a loss of
telomeric DNA at the rate of approximately ~ 3-5
bp/generation. This shortening is accompanied by a
gradual increase in cell and chromosome loss. By
approximately generation 100, most cells perish
[7,57,58,162].
Counter 1996 The roles of telomeres and telomerase in cell life span
ADDED STUFF


ALT is a recombination-based telomere maintenace mechanism used by some human cancers to maintain cellular immortality. ALT cells typically have widly long, heterogenous telomeres. The exact molecular mechanism involved in this pathway are unknown. Iyer 2005 found that a STN1 gene mutation can initiate ALT in the yeast Kluyveromyces lactis. These ALT-like yeast cells experience a rapid telomere shortening when WT STN1 is re-introduced (Iyer 2005). There aren't any neat figures in this paper to talk about, so I'm going to try to replicate the protein multiple sequence alignment from Iyer 2005 figure 4. STN1 is aligned from K. lactis S. cerevisiae (Sc) and Candida glabrata (Cgl).

![Yeast_Protein_Alignment](/Assets/Yeast_Protein_Alignment.jpg "Yeast_Protein_Alignment")

#### Obtaining STN1 From 3 Yeast Organisms
The paper states that the S. cerevisiae (Sc) and Candida glabrata (Cgl) GenBank accession numbers are P_38960 and XP_448655. HOWEVER, it wasn't that easy to find the sequences cause those accession numbers are from back in 2005. NCBI says "The following term was not found in Nucleotide: P_38960." The XP_448655 is here: https://www.ncbi.nlm.nih.gov/protein/XP_448655. I had trouble finding the sequence for K. lactis they were talking about. Here's the crazy search term that I used on NCBI:

and it's the only result (other than whole chromosome chunks) w/ this crazy search term I made:
(((stn1) NOT "Pyrenophora tritici-repentis"[porgn:__txid45151] NOT "Fusarium fujikuroi"[porgn:__txid5127] NOT "[Candida] glabrata"[porgn:__txid5478] NOT "Hortaea werneckii"[porgn:__txid91943] NOT "Saccharomyces cerevisiae"[porgn:__txid4932]) NOT "Metarhizium robertsii"[porgn:__txid568076] NOT "Fusarium sp. FIESC_5 CS3069"[porgn:__txid1318460] NOT "Fusarium pseudograminearum CS3487"[porgn:__txid1318458] NOT "Fusarium pseudograminearum CS3427"[porgn:__txid1318457] NOT "Fusarium pseudograminearum CS3220"[porgn:__txid1318456] NOT "Fonsecaea multimorphosa"[porgn:__txid979981] NOT "Cladophialophora immunda"[porgn:__txid569365] NOT "Aspergillus nidulans FGSC A4"[porgn:__txid227321] NOT "Candida viswanathii"[porgn:__txid5486] NOT "Zygosaccharomyces bailii"[porgn:__txid4954] NOT "Metarhizium anisopliae"[porgn:__txid5530] NOT "Aspergillus flavus"[porgn:__txid5059] NOT "Talaromyces atroroseus"[porgn:__txid1441469] NOT "[Candida] auris"[porgn:__txid498019] NOT "Zygosaccharomyces rouxii"[porgn:__txid4956] NOT "[Candida] boidinii"[porgn:__txid5477] NOT "Komagataella phaffii"[porgn:__txid460519] NOT "Aspergillus fumigatus"[porgn:__txid746128] NOT "Candida albicans SC5314"[porgn:__txid237561] NOT "Yarrowia lipolytica"[porgn:__txid4952]) AND "Kluyveromyces lactis"[porgn:__txid28985] 

... I'm pretty sure that it's NCBI Reference Sequence: XM_452728.1, titled "Kluyveromyces lactis uncharacterized protein (KLLA0_C11825g), partial mRNA" BECAUSE that entry has this note "cerevisiae YDR082W STN1 Protein involved in telomere length regulation functions in telomere metabolism during late S phase." S. ceriviasa yeast was easier to find cause it's got stn1p in the title :) https://www.ncbi.nlm.nih.gov/nuccore/NM_001180390.1
                     
#### Attempt to Replicate Lyer 2005 Figure 4
I used the R msa library to align the three yeast organism STN1 proteins that were obtained above. I'm not too happy with this alignment :( The authors state that "The protein
sequences were aligned in ClustalW using default values." HOWEVER, the R msa package that I used was set to use ClustalW default values and it didn't replicate Figure 4 from Lyer 2005. The protein sequences seem to match (at least by eye). I'm not really sure where I went wrong here. I've only played with sequence alignment a little bit, so I'm probably doing something wrong. Anyway, here's the code:

```r
setwd("/media/david/Linux/Introns_Exons_and_Promoters/Yeast_STN1_ALT/Proteins")

library(msa)
All_Yeast_STN1_Protein <- "All_Yeast_STN1_Protein.fasta"
All_Yeast_STN1_Protein <- readDNAStringSet(All_Yeast_STN1_Protein)
All_Yeast_STN1_Protein_alignment <- msa(All_Yeast_STN1_Protein)
All_Yeast_STN1_Protein_alignment
msaPrettyPrint(All_Yeast_STN1_Protein_alignment, output="pdf", showNames="none",
showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
```

![Aligning_Yeast_STN1_Proteins](/Assets/Aligning_Yeast_STN1_Proteins.jpg "Aligning_Yeast_STN1_Proteins")


#### Aligning Human, Yeast, and Frog STN1
Cohen 2002 reported that Xenopus laevis form extrachromosomal circular telomeric DNA. This is commonly associated with ALT! I briefly looked around for reptile ALT and this is the closest thing I could find. I'm not convinced that the activity reported by Cohen 2002 was actually ALT. This alignment isn't really related to anything. I just made it mostly for fun ... well, it's kinda connected to ALT and I know a herpetologist that might enjoy seeing frogs getting included in this repo ;) 

Cohen 2002 Formation of extrachromosomal circles from telomeric DNA in Xenopus laevis 
```r
Klactis_Xenopus_Human_STN1_Proteins <- "Klactis_Xenopus_Human_STN1_Proteins.fasta"
Klactis_Xenopus_Human_STN1_Proteins_AA <- readAAStringSet(Klactis_Xenopus_Human_STN1_Proteins)
Klactis_Xenopus_Human_STN1_Proteins_AA_alignment <- msa(Klactis_Xenopus_Human_STN1_Proteins_AA)
Klactis_Xenopus_Human_STN1_Proteins_AA_alignment
msaPrettyPrint(Klactis_Xenopus_Human_STN1_Proteins_AA_alignment, output="pdf", showNames="none",
showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
```
![Klactis_Xenopus_Human_STN1_Proteins_AA_alignment](/Assets/Klactis_Xenopus_Human_STN1_Proteins_AA_alignment.jpg "Klactis_Xenopus_Human_STN1_Proteins_AA_alignment")

# Citations
NOTE: I AM WAYYYYYY BEHIND ON UPDATING THE CITATIONS LIST ... SORRY! Please message me if you have any questions about the papers that I mentioned throughout this review :)
* Cheng 2012 Caenorhabditis elegans POT-2 telomere protein represses a mode of alternative lengthening of telomeres with normal telomere lengths
* Cong 1999 The human telomerase catalytic subunit hTERT: organization of the gene and characterization of the promoter
* Lovejoy 2012 PLoS Genet Loss of ATRX, genome instability, altered DNA damage response hallmarks of ALT pathway
* Clynes 2013 Curr Opin Genet Dev ATRX and the replication of structured DNA
* Gonzalo 2006 DNA methyltransferases control telomere length and telomere recombination in mammalian cells.pdf
* Atkinson 2005 Lack of Telomerase Gene Expression in Alternative Lengthening of Telomere Cells Is Associated with Chromatin Remodeling of the hTR and hTERT Gene Promoters
* Kumakura 2005 Reversible Conversion of Immortal Human Cells from Telomerase-Positive to Telomerase-Negative Cells
* Cong 1999 The human telomerase catalytic subunit hTERT: organization of the gene and characterization of the promoter
* Dessain 2000 Methylation of the Human Telomerase Gene CpG Island
* Cheng 2012 Caenorhabditis elegans POT-2 telomere protein represses a mode of alternative lengthening of telomeres with normal telomere lengths
* Boeck 2016 The time resolved transcriptome of C. elegans
* Iyer 2005 A Mutation in the STN1 Gene Triggers an Alternative Lengthening of Telomere-Like Runaway Recombinational Telomere Elongation and Rapid Deletion in Yeast
* Cohen 2002 Formation of extrachromosomal circles from telomeric DNA in Xenopus laevis 
* Hu 2017 Imetelstat, a Telomerase Inhibitor, Is Capable of Depleting Myelofibrosis Hematopoietic Stem Cells and Progenitor Cells
* Weissman 2001 Telomere Shortening Accompanies Increased Cell Cycle Activity during Serial Transplantation of Hematopoietic Stem Cells
* Cesare 2010 Alternative lengthening of telomeres: models, mechanisms and implications
* Shay 2012 Cancer and Telomeres −− An ALTernative to Telomerase
* Harley 2008 Telomerase and cancer therapeutics
* Sahin 2012 Antitelomerase therapy provokes ALT and mitochondrial adaptive mechanisms in cancer
* Gocha 2013 Human Sarcomas Are Mosaic for Telomerase-Dependent and Telomerase-Independent Telomere Maintenance Mechanisms
* Kalmbach 2014 Telomere Length Reprogramming in Embryos and Stem Cells
* Huang 2014 Telomere regulation in pluripotent stem cells
* Dagg 2017 Extensive Proliferation of Human Cancer Cells with Ever-Shorter Telomeres
* Suda 2002 Interchromosomal Telomere Length Variation
* Cristofari 2006 Telomere length homeostasis requires that telomerase levels are limiting
* Gavory 2002 Minimum length requirement of the alignment domain of human telomerase RNA to sustain catalytic activity in vitro
* Murnane 1994 Telomere dynamics in an immortal human cell line
* Telomere elongation in immortal human cells without detectable telomerase activity



