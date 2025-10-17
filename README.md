After porting to Unix run the following commands **once**:
```bash
dos2unix pipeline.sh
chmod +x pipeline.sh
```

To run the script:
```bash
./pipeline.sh encode.bed rep1.bam rep2.bam control.bam
```

Make sure the ENCODE results are already unzipped:
```bash
gunzip encode.bed.gz
```