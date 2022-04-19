# Lineage v0.0.1
## Dependencies (latest)
* Python 3
* ete3

### Create lineage environment
```
conda create --name lineage python
conda activate lineage
pip install ete3
```

### Add this to your PATH
```
export PATH="/gpfs1/scratch/ryan/tools/lineage:$PATH"
```

### Execute
```
lineage.py -i input.csv -o output.csv
```