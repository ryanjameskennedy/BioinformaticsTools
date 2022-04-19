# Lineage_retrieval v0.0.1
## Dependencies (latest)
* Python 3
* ete3

### Create lineage environment
```
conda create --name lineage
conda activate lineage
conda install -c bioconda ete3
```

### Add this to your PATH
```
export PATH="/gpfs1/scratch/ryan/tools/MEGAN:$PATH"
```

### Execute
```
lineage_retrieval.py -l Lineage_SpeciesAll.txt -i Species_SpeciesAll.txt -o .txt
```