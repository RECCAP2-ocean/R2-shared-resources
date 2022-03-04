# RECCAP2-ocean shared resources

This is a code repository where shared code can be uploaded. Contact Luke if you've got code that you think would be useful for the entire RECCAP2 ocean group. 

## Installation and activation

```bash
git clone https://github.com/RECCAP2-ocean/RECCAP2-shared-resources.git 

cd RECCAP2-shared-resources

conda env create -f environment.yml
conda activate reccap2
```

## Contents

- `data/regions`: masks for the RECCAP study
- `data/reccap2_data.yml`: a file containing download info for RECCAP2 ocean data. Used by `scripts.data`
  
- `scripts/regions`: scripts to make the masks
- `script/data`: For an example see `notebooks/reccap2-load-ocean-data.ipynb`
  1. Downloads data from the server (may take some time) using the configuration in the `data/reccap2_data.yml` file
  2. Open and transform the data. Gives summary of changes made to homogenise data
  3. Merge data based on variable

- `scripts/plot`: sets the configuration for plotting to have a nice simple style. Also supports mapping with `cartopy`
