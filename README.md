# Epigenomics Project - Analysis of ENCODE ChIP-Seq Data for Transcription Factors

After porting to Unix run the following commands **once**:
```bash
dos2unix pipeline.sh
chmod +x pipeline.sh
```

## To run the script:
```bash
./pipeline.sh encode.bed blacklist.bed chromHMM.bed rep1.bam rep2.bam control.bam
```

Make sure that every file is already unzipped:
```bash
gunzip encode.bed.gz
gunzip blacklist.bed.gz
...
```

Information on the available options and generic tool documentation is available under the ```-h``` flag:
```bash
./pipeline.sh -h
```

If you want the pipeline to generate more beautiful boxplots than the R standard, simply include your own R script in the working directory, and make sure to name it ```boxplots.R```.