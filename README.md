After porting to Unix run the following commands **once**:
```bash
dos2unix pipeline.sh
chmod +x pipeline.sh
```

To run the script:
```bash
./pipeline.sh encode.bed blacklist.bed rep1.bam rep2.bam control.bam
```

Make sure the ENCODE results and the blacklist are already unzipped:
```bash
gunzip encode.bed.gz
gunzip blacklist.bed.gz
```

Information on the available options and generic tool documentation is available under the ```-h``` flag:
```bash
./pipeline.sh -h
```