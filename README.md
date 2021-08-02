# Model-Contextualization-Pipeline
A tutorial for contextualizing metabolic models, using Pseudomonas Aeruginosa as an example


## Table of Contents
```
1.) Introduction 

2.) Find a transcriptomic data set which provides specific context for a given model

3.) Download and clean data 

4.) Use RIPTiDe to contextualize the model using the data 

5.) Troubleshooting Potential Issues 
```



1.) Introduction
```
Welcome to the Model Contextualization Pipeline! THe purpose for this repository is to provide a tutorial for contextualizing metabolic models with RIPTiDe. I walk through the process step by step, using a P. Aeruginosa model and relevant transcriptomes as an example. Essentially, the project can be split into 4 major phases:

Data Collection: If one desires to contextualize a model for a specific environment, it is important to find an appropriate, updated model first. Then, he or she must find transcriptomic data sets that appropriately reflect the transcriptional activity of an organism in a relevant environment. The model, or GENRE, and the transcriptomes are not necessarily easy to find. I provide the resources I used to find the materials for my project, but the resources are not limited to those alone. 

Data Proccessing: Transcriptomic data sets are extremely large, and usually a high power computer is needed. I am using the University of Virginia's HPC, Rivanna. The data must be made congruent enough to fit with the GENRE before it can be passed through RIPTiDe. 

RIPTiDe Processing: The GENRE and the transcriptome are passed through RIPTiDe to output the contextualized model. 


The aim of this tutorial is to guide those through the process of using and understanding RIPTiDe, so more analysis of contextualized models can be done. 


Download Anaconda: https://www.anaconda.com/products/individual#Downloads
Download jupyter notebook in Anaconda prompt: pip install jupyterlab
Download Ubuntu: https://ubuntu.com/download/desktop
```


2.) Data Collection
```
For computational modeling, I used CobraPy. For this reason, my tutorial is limited to CobraPy syntax, which can be found in section 2.1

Table Of Contents:
Genre
Transcriptome
Genome




GENRE
The first step is to find an up-to-date metabolic model of a given organism (GENRE).
GENRE's are mathematical representations of all of the metabolic processes included in a given organism.
The model will be used in RIPTiDe as the foundational basis of the more contextualized model. 

There are many models for use on BioModels, a public repository for mathematical models. There are also easy to access models on GitHub. 

BioModels website: https://www.ebi.ac.uk/biomodels/
GitHub ModelSEED: https://github.com/ModelSEED/ModelSEEDDatabase
Ex. of GitHub as a modeling resource: https://github.com/cdanielmachado/carveme



For my project, I accessed my model via GitHub. The model is of P. Aeruginosa UCBPP-PA14 and is named iPau21. 

Models can be in many formats, the most common of which are .json and .xml. Each of these require different syntax in order to import into Jupyter Notebook.
See "CobraPy Modeling" section for more detail
Keep track of where exactly you save the model on your computer, as you will need to move your transcriptomic data there later. 

Transcriptome
Next, find a transcriptome of the organism in a paritular environment that would provide the desired context for the model. 
Appropriate transcriptomes can demand extensive searching. However, they can be found in many places.
Using Google Scholar to find relevant research articles is an effective strategy, as the transcriptomic data sets are sometimes attached. NCBI SRA(Sequence Read Archives) 
is the definitive resource for transcriptomic data sets, but it can be challanging to find one that is ideal for the desired context. For that reason, it might be more effective 
to search for articles instead of navigating the SRA website directly.

NCBI SRA website:https://www.ncbi.nlm.nih.gov/sra
Google Scholar:https://scholar.google.com/

For my project, I wished to compare the contextualized models of P. Aeruginosa in different infection environments in the human body. Using the article 
"Genetically diverse Pseudomonas aeruginosa populations display similar transcriptomic profiles in a cystic fibrosis explanted lung," I found a link to the transcriptome in the 
SRA database. THe transcriptome satisfied the projects needs. It provided data that could contextualize the iPau21 model to simulate P. Aeruginosa metabolic activity
in the human lung. 

I also used "Major Transcriptome Changes Accompany the Growth of Pseudomonas aeruginosa in Blood from Patients with Severe Thermal Injuries" to find a transcriptome that
could provide context for P. Aeruginosa infection in a burn wound. The actual dataset was also stored in the NCBI SRA database. 

These two transcriptomes provided context for infections in different environments, which allowed me to compare the differences in the metabolic activity in these environments.


Genome
A complete reference genome of the organism in question is needed to processes the transcriptomic data. NCBI genome is a great place to find a genome that would serve this purpose.

NCBI Genome: https://www.ncbi.nlm.nih.gov/genome
```
2.1) CobraPy Modeling

# General Documentation


## Stuff to import at beginning

Load packages & model
```
import cobra
import xlrd
import pandas as pd
import pickle
from libsbml import *
from pprint import pprint
from cobra.flux_analysis import gapfill
from copy import *
from cobra import *
from cobra.flux_analysis import *
from LJD_Functions_DDP import *
from DDP_Functions  import *
from time import time
```                                  

### To import a model
```
model =cobra.io.load_json_model("file_name")    #imports the specific model if .json
model =cobra.io.read_sbml_model("file_name")   #If .xml file
```
### To optimmize the model
```
model.optimize()
```
### check the objective function via
```
model.summary()
```
### To change the reaction fluxes
 This is a standard way to ensure nutrition is entering the cell. If it isn't, the cell won't grow
```
for i in model.exchanges:
         i.lower_bound = -1000  
```
 You can also use model.reactions or any other group of reactions you would like to change the fluxes of


### To change the media
```model = changeMedia_PA_LJD(model,1)```
#### IMPORTANT- if you have trouble using this function like I did, open the LJD file and just copy the code for the 
specific media specifically. Sometimes you must change the exchange fluxes manually. 



## retrieving reactions

### returns all reactions using loop
```
for x in model.reactions:
    print(x.id+"\t"+x.name)  
```
### or more simply
```
model.reactions
```
### returns specific reactions using ID

```
model.reactions.get_by_id("") 
```
### or 

```
model.reactions.ID
```
Use ModelSeed database to look up reactions- can find name, ID, metabolites, all info involved






## retreving genes

### returns all genes using loop
```
for x in model.genes:
    print(x.id+"\t"+x.name)  
```
### or more simply
```
model.genes
```
### returns specific genes using ID

```
model.genes.get_by_id("") 
```
### or 

```
model.genes.ID

```
Use Kegg to look up gene information





## retreving metabolites

### returns all metabolites using loop
```
for x in model.metabolites:
    print(x.id+"\t"+x.name)  
```
### or more simply
```
model.metabolites
```
### returns specific metabolites using ID

```
model.metabolites.get_by_id("") 
```
### or 

```
model.metabolites.ID
```

Use ModelSeed to look up metabolites



    
## Methods and their physical correlative:
```
Exchange reactions- introducing metabolites to the model- tracks what the organism is producing/consuming
Transport reactions- metabolite movement from extracellular space to cytosol and vice versa
Metabolic reactions- Chemical reactions in the cytosolic compartment that represent the organisms metabolism

Exchange rxns- in or out of extracellular space
Transport- in or out of cytosolic space
metabolic rxns- conversions inside the cytosol

lower bound of Exchange reactions(negative)- nutrition entering model- (-1000 = 1000 into the cell)
upper bound of Exchange reactions(positive) - nutrition exiting model- (1000 = 1000 out of the cell)        
```         
         
         
### How to check flux of certain reactions after optimizing-

```
model.optimize().fluzes.reactionID #individual reaction and flux

model.optimize().fluxes  #fluxes for every reaction 
```
### How to export data to an excel file-
```
soln = model.optimize().fluxes
soln.to_excel(“solution.xlsx”) #Save pandas dataframe to excel
```


### How to change objective funciton-
```
model.objective = model.reactions.rxnID
```
### How to add rxns to model-
```
reaction1 = Reaction(‘Desired Rxn name’)
reaction1.name = ‘Desired Rxn name’
reaction1.subsystem = ‘subsystem of desire’
reaction1.lower_bound = 0
reaction1.upper_bound = 1000
reaction1.add_metabolites({
   model.metabolites.ID: -1.0- represents stoichiometric coefficient})
model.add_reactions([reaction1])
```


### How to knock out components of metabolism
```
with model:
    model.(can be genes or reactions here).ID.knock_out()
    print(model.optimize())
    ```
    
```
 This is used to compare how essential a gene or reaction is to the objective function and other metabolic processes
 
 3.) Data Processing
 Because of the large computational demands of processing transcriptomic data, a HPC system is often used. Remote HPC servers can be accesssed remotely off of a base platform. For my project I used Ubuntu, so this section will include only Ubuntu syntax. See section 3.1 for basic Ubuntu commands. 

### Data Processing Table of Contents:
1.)Index Genome
2.) Downloading Transcriptomic Data sets
3.) Downloading sickle
4.)Processing reads-single and paired 
5.)Moving Output .tsv 











### Index Genome
This first step is to find and download the reference genome. This should be done locally. When downloading the genome, the file should be in FASTA format, having a .fasta or .fa extention. When downloading off of NCBI, the file is often stored in a .fna.gz format. To convert the file into the appropriate format, the file must first be unzipped. In Ubuntu, locate the location of the downloaded file and move it to the folder you are working in: I simply used my base directory. 


My code, which moves the downloaded genome from NCBI into my workspace:
```
mv /mnt/c/Users/rglad/Downloads/GCF_000006765.1_ASM676v1_genomic.fna.gz /home/rgladwell2
```

Now, when ls is called in the home directory, the file name should appear. We can now unzip the file uzing the command gunzip filename.fna.gz. When ls is called now, the file name should appear without the .gz format. Next the file should be converted from .fna to .fa. Use the command cp filename.fna filename.fa to accomplish this. 

My example code:
```
gunzip GCF_000006765.1_ASM676v1_genomic.fna.gz
 cp GCF_000006765.1_ASM676v1_genomic.fna GCF_000006765.1_ASM676v1_genomic.fa
 ```
  Use the rm command to get rid of any files you no longer need, such as the remaining .fna file.
 
 It is now time to upload the .fa file to an HPC. For my HPC, I signed into sftp(secure file transfer protocol) on the rivanna server then used the put command to put the file in my scratch directory.

```
put GCF_000006765.1_ASM676v1_genomic.fa /scratch/username/
```
 When executing ls in the scratch directory, the file should now appear
 The index genome can now be constructed. Run the following script in order to generate the index:
 
 ```
 module load gcc
module load bowtie2
bowtie2-build filename.fa Desired_filename_genes
```
My example code:
```
module load gcc
module load bowtie2
bowtie2-build bowtie2-build GCF_000006765.1_ASM676v1_genomic.fa Index_genome
```
When running ls, 6 new files ending in .bt2 should have appeared bearing the name that you chose in the last step:
```
Index_genome.1.bt2                  
Index_genome.2.bt2                  
Index_genome.3.bt2                   
Index_genome.4.bt2                  
Index_genome.rev.1.bt2              
Index_genome.rev.2.bt2
```
This is the index Genome.
Some troubleshooting might be required here. For my project, the reference genome I chose created an invalid reference genome, so I had to download gene data from another webstie and format the file to work with my model. I used the same process as above, just with a different file. 

### Downloading Transcriptomic data sets
After finding the appropriate transcriptomic data sets, they will need to be downloaded onto the hpc. 

Every transcriptome on sra NCBI has a corresponding run key. These numbers strings typically start with "SRR." Make note of the run key(including the characters at the beginning), as this is what is used to download the transcriptome. 

Use the following code to create a script to downlaod any transcriptome to a desired hpc:
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 120:00:00#
SBATCH --mem=150000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard

module load sratoolkit/2.10.5
fastq-dump --split-files ${1}
```
I named my script sra.slurm. After writing the script use the following code to execut it for the desired data sets:
```
for x in SRR#; do sbatch scriptname $x; done
```
My example code:
``` 
for x in SRR2080030 SRR7777123; do sbatch sra.slurm $x; done
```
This will download and save the transcriptomes as their run key.fastq. Some download as single end reads and others download as paired end reads. Single ends will download in one file while paired ends will download in two. 

These are the transcriptomic data sets. 


### Downloading sickle
Next, on your local Ubuntu, download sickle. The tutorial link is below:
```
https://github.com/najoshi/sickle
```
After downloading, move the file to your hpc server, to the same directory as your transcriptomic data sets. Sickle is used soon. Note that you also have to download zlib for using sickle. Links to zlib here:
```
http://www.zlib.net/
```

### Processing reads
Now, we want to construct a .csv or .tsv file that contains the gene ID, length, and # of hits. To do this , create the following script:
so first use nano "script name here" to create the script location-for uniformity, I have given all files the name of their purpose- name them whatever you like as long as you maintain consistency.

First is the script for single end reads. 
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 50:00:00
#SBATCH --mem=150000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard

module load gcc/9.2.0 bowtie2/2.2.9
module load samtools/1.10

#Quality trim samples based on quality scores
#first is the location of your sickle file.
/location/of/sickle se -t sanger \
        -f input_SRR_file.fastq \
        -o intermediate_trim_file.fastq \
       
       
#Single Reads
bowtie2 -q --fr -x Index_genome \
        -U intermediate_trim_file.fastq \
        -S single.reads.sam
samtools view -h single.reads.sam -o single.reads.bam
samtools sort single.reads.bam -o single.reads.sort.bam
rm single.reads.bam single.reads.sam


#Sort again
samtools sort single.reads.sort.bam -o single.reads.sort2.bam

#readable format
samtools idxstats single.reads.sort2.bam > Output.tsv
```
Here is my sample code for an example single read:
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 50:00:00
#SBATCH --mem=150000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard


module load gcc/9.2.0 bowtie2/2.2.9
module load samtools/1.10

#Quality trim samples based on quality scores
/scratch/rag6yy/sickle se -t sanger \
        -f SRR7777123_1.fastq \
        -o CFLung.trim.fastq \

#Single Reads
bowtie2 -q --fr -x Pseudomonas_genes \
        -U CFLung.trim.fastq \
        -S single.reads.sam
samtools view -h single.reads.sam -o single.reads.bam
samtools sort single.reads.bam -o single.reads.sort.bam
rm single.reads.bam single.reads.sam

#Sort again
samtools sort single.reads.sort.bam -o single.reads.sort2.bam

#readable format
samtools idxstats single.reads.sort2.bam > CF.lung.tsv
```
Here is my generalized code for a list of single reads:
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 50:00:00
#SBATCH --mem=150000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard


module load gcc/9.2.0 bowtie2/2.2.9
module load samtools/1.10

#Quality trim samples based on quality scores
/scratch/rag6yy/sickle se -t sanger \
        -f ${1}_1.fastq \
        -o ${1}.trim.fastq \

#Single Reads
bowtie2 -q --fr -x Pseudomonas_genes \
        -U ${1}.trim.fastq \
        -S ${1}single.reads.sam
samtools view -h ${1}single.reads.sam -o ${1}single.reads.bam
samtools sort ${1}single.reads.bam -o ${1}single.reads.sort.bam
rm ${1}single.reads.bam

#Sort again
samtools sort ${1}single.reads.sort.bam -o ${1}single.reads.sort2.bam


#readable format
samtools idxstats ${1}single.reads.sort2.bam > ${1}.tsv
rm ${1}single.reads.sort2.bam ${1}single.reads.sort.bam
rm ${1}single.reads.sam
```
Then run:
```
for x in SRR7777123(or more, you could do all of your single read files here); do sbatch processing_single_reads $x; done
```
Next is the script for paired end reads- it is not all that different. 
```
#!/bin/bash
#SBATCH--nodes=1
#SBATCH--ntasks-per-node=16
#SBATCH -t 50:00:00
#SBATCH --mem=150000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard

module load gcc/9.2.0 bowtie2/2.2.9
module load samtools/1.10

#Quality trim samples based on quality scores
/location/of/sickle pe -t sanger -q 25 -l 21 \
        -f  Input_transcriptome_1.fastq \
        -r  Input_transcriptome_2.fastq \
        -o Input_transcriptome_1.trim.fastq \
        -p Input_transcriptome_2.trim.fastq \
        -s Orphan_reads.fastq
        
        
 #paired reads
bowtie2 -q --fr -x Index_genome
        -1 Input_transcriptome_1.trim.fastq \
        -2 Input_transcriptome_2.trim.fastq \
        -S pair.reads.sam
samtools view -bS pair.reads.sam > pair.reads.bam
samtools sort -o pair.reads.sort.bam pair.reads.bam
rm pair.reads.sam pair.reads.bam

#Orphaned reads
bowtie2 -q -x Index_genome
        -U Orphan_reads.fastq \
        -S orphanReads.sam
samtools view -bS orphanReads.sam > orphanReads.bam
samtools sort -o orphanReads.sort.bam orphanReads.bam
rm orphanReads.bam orphanReads.sam

#merge Alignments
samtools merge reads.merge.bam pair.reads.sort.bam orphanReads.sort.bam
samtools sort -o reads.merge.sort.bam reads.merge.bam
rm reads.merge.bam orphanReads.sort.bam pair.reads.sort.bam

#Readable format
samtools idxstats reads.merge.sort.bam > Output.tsv
rm reads.merge.sort.bam
```
Here is my sample code for a paired ends:
```
#!/bin/bash
#SBATCH--nodes=1
#SBATCH--ntasks-per-node=16
#SBATCH -t 50:00:00
#SBATCH --mem=150000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard

module load gcc/9.2.0 bowtie2/2.2.9
module load samtools/1.10

#Quality trim samples based on quality scores
/scratch/rag6yy/sickle pe -t sanger -q 25 -l 21 \
        -f  SRR2080030_1.fastq \
        -r  SRR2080030_2.fastq \
        -o BurnWound1.trim.fastq \
        -p BurnWound2.trim.fastq \
        -s BurnOrphan.fastq
        
#paired reads
bowtie2 -q --fr -x Pseudomonas_genes \
        -1 BurnWound1.trim.fastq \
        -2 BurnWound2.trim.fastq \
        -S pair.reads.sam
samtools view -bS pair.reads.sam > pair.reads.bam
samtools sort -o pair.reads.sort.bam pair.reads.bam
rm pair.reads.sam pair.reads.bam

#Orphaned reads
bowtie2 -q -x Pseudomonas_genes \
        -U BurnOrphan.fastq \
        -S orphanReads.sam
samtools view -bS orphanReads.sam > orphanReads.bam
samtools sort -o orphanReads.sort.bam orphanReads.bam
rm orphanReads.bam orphanReads.sam

#merge Alignments
samtools merge reads.merge.bam pair.reads.sort.bam orphanReads.sort.bam
samtools sort -o reads.merge.sort.bam reads.merge.bam
rm reads.merge.bam orphanReads.sort.bam pair.reads.sort.bam

#Readable format
samtools idxstats reads.merge.sort.bam > burn_pseudomonas.tsv
rm reads.merge.sort.bam
```
Here is my sample code for a list of paired reads:
```
#!/bin/bash
#SBATCH--nodes=1
#SBATCH--ntasks-per-node=16
#SBATCH -t 50:00:00
#SBATCH --mem=150000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard

module load gcc/9.2.0 bowtie2/2.2.9
module load samtools/1.10

#Quality trim samples based on quality scores
/scratch/rag6yy/sickle pe -t sanger -q 25 -l 21 \
        -f  ${1}_1.fastq \
        -r  ${1}_2.fastq \
        -o ${1}1.trim.fastq \
        -p ${1}2.trim.fastq \
        -s ${1}Orphan.fastq

#paired reads
bowtie2 -q --fr -x Pseudomonas_genes \
        -1 ${1}1.trim.fastq \
        -2 ${1}2.trim.fastq \
        -S ${1}pair.reads.sam
samtools view -bS ${1}pair.reads.sam > ${1}pair.reads.bam
samtools sort -o ${1}pair.reads.sort.bam ${1}pair.reads.bam
rm ${1}pair.reads.sam ${1}pair.reads.bam

#Orphaned reads
bowtie2 -q -x Pseudomonas_genes \
        -U ${1}Orphan.fastq \
        -S ${1}orphanReads.sam
samtools view -bS ${1}orphanReads.sam > ${1}orphanReads.bam
samtools sort -o ${1}orphanReads.sort.bam ${1}orphanReads.bam
rm ${1}orphanReads.bam ${1}orphanReads.sam

#merge Alignments
samtools merge ${1}reads.merge.bam ${1}pair.reads.sort.bam ${1}orphanReads.sort.bam
samtools sort -o ${1}reads.merge.sort.bam ${1}reads.merge.bam
rm ${1}reads.merge.bam ${1}orphanReads.sort.bam ${1}pair.reads.sort.bam

#Readable format
samtools idxstats ${1}reads.merge.sort.bam > ${1}.tsv
rm ${1}reads.merge.sort.bam
```
### Moving output
It's now time to move output files to jupyter notebook for analysis with RIPTiDe. Use the same sftp procedure as before but instead of put, use get. 
get /current/location/Output.tsv /desired/location/path
Here is my sample code:
```
 get /scratch/rag6yy/CF.lung.tsv /mnt/c/Users/rglad/Documents/REU Modeling
 ```
 Here's how to move all .tsv at once:
 ```
 cd /location/of/files
get *.tsv /local/location
```
My sample code:
```
cd /scratch/rag6yy/
get *.tsv /mnt/c/Users/rglad/Documents
```
Once you get your output .tsv files into the same location as your GENRE, you are ready to create a contextualized model!


3.1) Ubuntu Basics
Links for Ubuntu tutorials
``` 
https://ubuntu.com/tutorials/install-ubuntu-desktop#1-overview

https://ubuntu.com/tutorials/command-line-for-beginners#1-overview
```
Common commands or processses I used:
```
cd filelocation
```
moves into a certain location
```
cd /scratch/username 
```
moves me into my scratch or cd /home/username moves me into my home

```
pwd 
```
outputs the current location
```
pwd
```
might output /home/usernmae if that is the current location

```
nano filename
```
creates a new file with the given name or moves into existing file with that name
```
nano workspace 
```
will either create a file named workspace or, if workspace already exists, move into it. 
```
rm filename
```
 deletes files with the file name
```
rm workspace
```
would delete the workspace file

```
mv filename location
```
moves a file to a desired location
```
mv workspace /home/username
```
moves the workspace file into my home directory
```
sudo apt-get 
```
or 
```
pip.exe install 
```
are ways to install packages. Each package is different and it is best to look up how to install them on ubuntu specifically.

```
sqeue -u username
```
checks how long a process has been going for.

4.) RIPTiDe
## Creating a contextualized model
A similar walkthrough of RIPTiDe can be found at the original GitHub repository:  
```
https://github.com/mjenior/riptide
```

By now, you should have a location where your processed transcriptomic data and GENRE are stored. 

Start by downloading RIPTiDe in the anaconda prompt using: pip install riptide==3.3.5

RIPTiDe works with .json models- if you have a .xml model simply use the following code to convert the model:
```
model =cobra.io.read_sbml_model("Model.xml")
cobra.io.save_json_model(model, "Model.json")
```
Here's my code:
```
model =cobra.io.read_sbml_model("iPau21.xml")
cobra.io.save_json_model(model, "iPau21.json")
```
Next we can enter jupyter notebook in the location where our materials are stored. Import everything for RIPTiDe and modeling:
```
import xlrd
import pandas as pd
import pickle
from libsbml import *
from pprint import pprint
from cobra.flux_analysis import gapfill
from copy import *
from cobra import *
from cobra.flux_analysis import *
from LJD_Functions_DDP import *
from DDP_Functions  import *
from time import time
from riptide import *

import cobra
import csv
import numpy as np
```
Next we need to read in our tsv with the following command: 
```
transcriptomic_data=riptide.read_transcription_file("Output.tsv")
```
Now to create the contextualized model:
```
model =cobra.io.load_json_model("Model.json")
for i in model.exchanges:
         i.lower_bound = -1000
         i.upper_bound = 1000
Contextualized_Model = riptide.contextualize(model=model,transcriptome=transcriptomic_data)
```
Here's my sample code:
```
lung2=riptide.read_transcription_file("CF.lung.tsv")

model =cobra.io.load_json_model("iPau21.json")
#I set the lower bound to -1000 and upper to 1000 to allow nutrition into my model because by default it is 0, and it is recommended to do this when working with RIPTiDe.
for i in model.exchanges:
         i.lower_bound = -1000
         i_upper_bound = 1000
PA_lung1 = riptide.contextualize(model=model,transcriptome=lung2)
```
CONGRATULATIONS! You have created a contextualized metabolic model!

5.) Troubleshooting
### ISSUES WITH GETTING A OCCUPIED .tsv FILE:

The first thing to do if your .tsv is empty, or if the output looks unexpected, is the comment out any rm commands in your processing scripts. Then, go through all intermediate files and make sure they are occupied with expected output. If you find one of the imtermediate files looks strange, focus your attention on the code right before you output that file. If intermediate files look fine, see below for more assisstance.

### Invalid reference genome:
If your .tsv isn't full and your intermediate files look fine, the problem might lie in your reference genome. Go through all your referene genome files and make sure they are occupied and look normal. They should all contain a significant amount of data- having a small amount of data in a .bt2 file might indicate the problem is with the reference genome. 

If you downloaded your file off of NCBI but you suspect it is not working, you might try downloading a raw genes read to use as your genome. Go to NCBI Gene and find the exact organism you want. Click the small link under the large link (Pardon the vague description). This should take you to a page with a table named Entrez records. Under direct links on the Gene section, click the link provided. A list of organisms should come up. Pick the one you need and scroll down until you find a reference assembly seciton. There should be a 'Genomic' section with a FASFTA file you can download. Attempt using this as your reference genome. (Make sure to recreate the .bt2 files and rewrite your code!)

If that doesn't work, consider getting in touch with a researcher who has interacted with the organism of your choice extensively. They might have a reference genome that they use that you can use.



### Improper Code:
Make sure all of your code is exact. Small details can throw the code off. For example, in my project I had one extra dollar sign in my generalized code that completely disabled my code. Go over your code finely if you think it might be a potential cause of issue.

### .fastq files
The transcriptomic files can download improperly on rare occasion, or the availible information might be subpar. Consider deleteing and redownloading the data. If this doesn't work, try replicating the process with different datasets. If it works with different data, consider contacting the researchers responsible for collecting the data for your dataset to ask for assistance. 

### Overloading an HPC
If you submit too many jobs at once, they may all ubruptly fail if the strain on the HPC is too much. If you think the processing is failing due to the HPC, try submitting less jobs at once. If that doesn't work but you still suspect the HPC, get in touch with the HPC's support. 
