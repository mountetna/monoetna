# Purpose  
The package `geoline` is aimed to automate submissions of biological data to GEO repository.  
## Install  
```bash
pip install .
```
Also requires `magby` to be installed (https://github.com/mountetna/monoetna/tree/master/etna/packages/magby)


## Usage
In a `jupyter` notebook

```python
from geoline import Geoline

url = 'https://magma.ucsf.edu'
token = 'YOUR_JANUS_TOKEN_HERE'
project = 'PROJECT_NAME_HERE'  # name of a project you want to submit data from. It should exist in magma
gl = Geoline.Geoline(url, token, project)

assay = 'ASSAY_NAME_HERE'  # accepted assays 'rna_seq', 'dna_seq', 'sc_seq'
primary_model = 'MODEL_NAME_HERE'  # name of a model that has most of the data (e.g. rna_seq) 
geo_meta = gl.seq_workflow(assay, primary_model)
```

You will be prompted with a series of interactive questions. Please follow the instructions in these questions.  
Output - `pandas.DataFrame` which can be further parsed and/or exported. Use this DataFrame to populate GEO metadata table
