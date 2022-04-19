
# Taxonomic identification pipeline
## Dependencies (latest)
* biopython
* blast
* pandas
* requests
* R

### Create taxonomy environment
```
conda create --name taxonomy biopython blast blast-legacy pandas requests r seaborn-base
```

### Change directory to OrthoANI environment
```
cd path/to/miniconda/envs/taxonomy
```

### Install OrthoANI
```
pip install orthoani
```