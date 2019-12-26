---
title: Testing recommendation engine on a UC signature(PZ349) from creeds and recommending
  top 3 datasets.
authors:
tags:
- Recommendation
- UC

tldr: Ran the recommendation engine on UC signature from creeds. Manually curated
  dataset and metadata were used which were prepared earlier. Recommended top 3 datasets.
---
```python
import requests
from bs4 import BeautifulSoup
import pandas as pd
from itertools import chain
GEO_URL = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
```
### Function to get CREEDS response for a signature id 


```python
CREEDS_URL = 'http://amp.pharm.mssm.edu/CREEDS/'
def get_creeds_response(sig_id):
    parameter = {"id" : sig_id}
    response = requests.get(CREEDS_URL + 'api', params=parameter)
    if response.status_code == 200:
        output = response.json()
    else:
        output = 'CREEDS server could not return output'
    return output
```
### Function to run creeds for up and down genes


```python
def run_creeds(up_genes, down_genes):
    payload = {
        'up_genes': up_genes,
        'dn_genes': down_genes,
        'direction': 'similar',
        'db_version': 'v1.0'
    }

    r = requests.post(CREEDS_URL + 'search', json=payload)
    print(r.status_code)
    response = r.json()

    creeds_response_df = pd.DataFrame.from_records(response).set_index('id')
    return creeds_response_df
```
### Function to get up, down and query organism from creeds response returned


```python
def get_creeds_sig_info(creeds_response_dict):
    desired_data = ['do_id', 'geo_id', 'cell_type', 'platform', 'disease_name', 'organism']
    interested_dict = dict([(i, creeds_response_dict.get(i)) for i in desired_data])
    down_genes = [i[0] for i in creeds_response_dict['down_genes']]
    up_genes = [j[0] for j in creeds_response_dict['up_genes']]
    organism_dict = {'human': 'hs', 'mouse': 'mm', 'rat': "rt"}
    query_org = organism_dict[interested_dict['organism']]
    return up_genes, down_genes, query_org
```
### Get creeds response for UC signature ID P439 


```python
creeds_data = get_creeds_response("dz:P439")
creeds_data
up_genes, down_genes, query_org = get_creeds_sig_info(creeds_data)
creeds_response_df = run_creeds(up_genes, down_genes)
creeds_response_df.head()
```
    200






<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>geo_id</th>
      <th>name</th>
      <th>signed_jaccard</th>
    </tr>
    <tr>
      <th>id</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>dz:475</th>
      <td>GSE16515</td>
      <td>[pancreatic cancer, http://disease-ontology.or...</td>
      <td>0.08240</td>
    </tr>
    <tr>
      <th>dz:610</th>
      <td>GSE19650</td>
      <td>[pancreatic invasive intraductal papillary-muc...</td>
      <td>0.07633</td>
    </tr>
    <tr>
      <th>dz:813</th>
      <td>GSE53431</td>
      <td>[psoriasis, http://disease-ontology.org/term/D...</td>
      <td>0.07299</td>
    </tr>
    <tr>
      <th>dz:587</th>
      <td>GSE11223</td>
      <td>[ulcerative colitis, http://disease-ontology.o...</td>
      <td>0.07106</td>
    </tr>
    <tr>
      <th>gene:2360</th>
      <td>GSE11759</td>
      <td>[Hnf4a, http://www.ncbi.nlm.nih.gov/gene/15378]</td>
      <td>0.06766</td>
    </tr>
  </tbody>
</table>
</div>



#### take top signature from the responses returned by CREEEDS


```python
bestsig_id = creeds_response_df.index.values[0]
bestsig_id
```




    'dz:475'




```python
bestsig_id2 = creeds_response_df.index.values[1]
bestsig_id2
```




    'dz:610'



### Function to get summary and title from GSE ID


```python
def get_summary_and_title(gse_id):
    url = GEO_URL + gse_id
   
    page = requests.get(url)
    soup = BeautifulSoup(page.text, "html.parser")
    title = soup.find("td", text="Title").find_next_sibling("td").text
    # pltfm_org = soup.find("td ", text="Platform organism").find_next_sibling("td").text
    # sample_org = soup.find("td", text="Sample organism").find_next_sibling("td").text
    exp_type = soup.find("td", text="Experiment type").find_next_sibling("td").text
    abstract = soup.find("td", text="Summary").find_next_sibling("td").text

    dataset_meta = {'gse_id': gse_id, 'title': title, 'exp_type': exp_type, 'abstract': abstract}
    return(dataset_meta)
```

```python
def query_jensen_api(input_string):
   query_string = input_string.replace(" ", "+")
   url = 'http://tagger.jensenlab.org/GetEntities?document=' + \
       query_string + '&entity_types=-2+-25+-26+-27+-21+-22+-23+0+-1+-3+ \
       -11+-24+-28+-29+-30+-31+-36&format=tsv'
   response = requests.get(url)
   response_jensen = pd.DataFrame([x.split('\t') for x in str(
       response.text).split("\n")], columns=["Name", "Annotation","Identifier"])
   response_jensen_wo_duplicates = response_jensen.drop_duplicates(["Name"])
   return response_jensen_wo_duplicates
```

```python
def get_metadata(sig_id):
   creeds_response_dict = get_creeds_response(sig_id)
   gse_id = creeds_response_dict['geo_id']
   data = get_summary_and_title(gse_id)
   jensen_output = query_jensen_api(data['abstract'])
   annotated_data = (annotate_biomedical_entities(jensen_output))
   return annotated_data

```

```python
biomedicalTermsJensenAnnotated = {'APO_phenotypes': -28,
                                  'BTO_Tissues': -25,
                                  'DOID_Diseases': -26,
                                  'ENVO_environments': -27,
                                  'FYPO_phenotypes': -29,
                                  'GOBiologicalProcess': -21,
                                  'GOCellularComponent': -22,
                                  'GOMolecularFunction': -23,
                                  'GOOther': -24,
                                  'MPheno_phenotypes': -30,
                                  'NBO_behaviors': -31,
                                  'NCBI_Chemicals': -1,
                                  'NCBI_Species': -2,
                                  'NCBI_Species_Proteins': -3,
                                  'Wikipedia': -11,
                                  'mammalian_phenotypes': -36}
```

```python
def annotate_biomedical_entities(response_jensen_wo_duplicates):
    """
    This function annotate words along with biomedical entities

    Parameters
    ----------
    response_jensen_wo_duplicates : pandas dataframe
         pandas dataframe having word along with annotation and identifier

    Returns
    -------
    dict
        dictionary where each entity contains list of
        words from input string as values
    """

    annotated_dict = dict()
    for index, row in response_jensen_wo_duplicates.iterrows():
        k = [key for (key, value) in biomedicalTermsJensenAnnotated.
             items() if int(row[1]) == value]
        if int(row[1]) > 1:
            k = ["Genes"]
        m = "Not_Known" if len(k) == 0 else k[0]
        annotated_dict.setdefault(m, [])
        annotated_dict[m].append(row[0])
    return annotated_dict

```
### get metadata for top response returned by CREEDS


```python
creeds_response_df.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>geo_id</th>
      <th>name</th>
      <th>signed_jaccard</th>
    </tr>
    <tr>
      <th>id</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>dz:475</th>
      <td>GSE16515</td>
      <td>[pancreatic cancer, http://disease-ontology.or...</td>
      <td>0.08240</td>
    </tr>
    <tr>
      <th>dz:610</th>
      <td>GSE19650</td>
      <td>[pancreatic invasive intraductal papillary-muc...</td>
      <td>0.07633</td>
    </tr>
    <tr>
      <th>dz:813</th>
      <td>GSE53431</td>
      <td>[psoriasis, http://disease-ontology.org/term/D...</td>
      <td>0.07299</td>
    </tr>
    <tr>
      <th>dz:587</th>
      <td>GSE11223</td>
      <td>[ulcerative colitis, http://disease-ontology.o...</td>
      <td>0.07106</td>
    </tr>
    <tr>
      <th>gene:2360</th>
      <td>GSE11759</td>
      <td>[Hnf4a, http://www.ncbi.nlm.nih.gov/gene/15378]</td>
      <td>0.06766</td>
    </tr>
  </tbody>
</table>
</div>




```python
top5_sigs = {}
```

```python
best_sig_metadata = get_metadata(bestsig_id2)
best_sig_metadata
```




    {'BTO_Tissues': ['pancreas'],
     'DOID_Diseases': ['adenoma', 'pancreatic cancer', 'carcinoma', 'neoplasm'],
     'GOBiologicalProcess': ['immune response']}




```python
for sig_id in creeds_response_df.index.values[0:8]:
    print(sig_id)
    if('gene' not in sig_id):
        top5_sigs[sig_id] = get_metadata(sig_id)
```
    dz:475
    dz:610
    dz:813
    dz:587
    gene:2360
    dz:609
    dz:604
    drug:3415



```python
for sig_id in creeds_response_df.index.values[0:8]:
    print(sig_id)
    if('gene' not in sig_id):
        metadata_dict = get_metadata(sig_id)
        disease_list = metadata_dict['DOID_Diseases']
        if('Disease' in disease_list):
            disease_list.remove('Disease')
        elif('disease' in disease_list):
            disease_list.remove('disease')
        tissue_list = metadata_dict.get('BTO_Tissues','')
        disease_list = [x.lower() for x in disease_list]
        tissue_list = [x.lower() for x in tissue_list]
        top5_sigs[sig_id] = (disease_list, tissue_list)
```
    dz:475
    dz:610
    dz:813
    dz:587
    gene:2360
    dz:609
    dz:604
    drug:3415



```python
top5_sigs
```




    {'dz:610': (['adenoma', 'pancreatic cancer', 'carcinoma', 'neoplasm'],
      ['pancreas']),
     'dz:813': (['psoriasis', 'psoriasis', 'skin disorder'], ['epidermis']),
     'dz:587': (['ulcerative colitis'], ['colon epithelial']),
     'dz:609': (['adenoma', 'pancreatic cancer', 'carcinoma', 'neoplasm'],
      ['pancreas']),
     'dz:604': (['pancreatic cancer', 'pancreatic ductal adenocarcinoma'], []),
     'drug:3415': (['cancer', 'colon cancer'],
      ['colorectal', 'caco-2 cells', 'caco-2']),
     'dz:475': (['pancreatic tumor'], [])}




```python
metadata_disease_columns = set(chain.from_iterable([top5_sigs[sig_id][0] for sig_id in top5_sigs.keys()]))
metadata_tissue_columns = set(chain.from_iterable([top5_sigs[sig_id][1] for sig_id in top5_sigs.keys()]))
```

```python
metadata_disease_columns
```




    {'adenoma',
     'cancer',
     'carcinoma',
     'colon cancer',
     'neoplasm',
     'pancreatic cancer',
     'pancreatic ductal adenocarcinoma',
     'pancreatic tumor',
     'psoriasis',
     'skin disorder',
     'ulcerative colitis'}




```python
metadata_tissue_columns
```




    {'caco-2',
     'caco-2 cells',
     'colon epithelial',
     'colorectal',
     'epidermis',
     'pancreas'}




```python
metadata_disease_columns_JENSEN = ['JENSEN_' + x for x in set.union(metadata_disease_columns, metadata_tissue_columns)]
metadata_disease_columns_JENSEN
```




    ['JENSEN_neoplasm',
     'JENSEN_psoriasis',
     'JENSEN_colon cancer',
     'JENSEN_colorectal',
     'JENSEN_colon epithelial',
     'JENSEN_caco-2',
     'JENSEN_skin disorder',
     'JENSEN_pancreatic tumor',
     'JENSEN_pancreas',
     'JENSEN_ulcerative colitis',
     'JENSEN_pancreatic ductal adenocarcinoma',
     'JENSEN_caco-2 cells',
     'JENSEN_cancer',
     'JENSEN_carcinoma',
     'JENSEN_adenoma',
     'JENSEN_pancreatic cancer',
     'JENSEN_epidermis']



### Get metadata from CREEDS


```python
creeds_metadata_dict_best = {}
```

```python
for sig_id in creeds_response_df.index.values[0:8]:
    if('gene' not in sig_id):
        creeds_response = get_creeds_response(bestsig_id)['cell_type']
        cell_type_list = [creeds_response]
        if isinstance(creeds_response, list):
            cell_type_list = [x['name'] for x in creeds_response]
        cell_type_list = [x.lower() for x in cell_type_list]
        creeds_metadata_dict_best[sig_id] = cell_type_list
```

```python
creeds_metadata_dict_best
```




    {'dz:475': ['pancreatic tissue'],
     'dz:610': ['pancreatic tissue'],
     'dz:813': ['pancreatic tissue'],
     'dz:587': ['pancreatic tissue'],
     'dz:609': ['pancreatic tissue'],
     'dz:604': ['pancreatic tissue'],
     'drug:3415': ['pancreatic tissue']}




```python
creeds_metadata_celltype_columns = set(chain.from_iterable([creeds_metadata_dict_best[sig_id] for sig_id in creeds_metadata_dict_best.keys()]))
creeds_metadata_celltype_columns_CREEDS = list(['CREEDS_' + x for x in creeds_metadata_celltype_columns])
creeds_metadata_celltype_columns_CREEDS
```




    ['CREEDS_pancreatic tissue']



### READ metadata file made earlier to get similarity scores


```python
import numpy as np
import pandas as pd
import requests
```

```python
meta_data = pd.read_csv('UC_signatures_metadata.csv')
```

```python
#### Check how many columns are common
np.intersect1d(list(meta_data.columns), list())
```

```python
np.intersect1d(list(meta_data.columns), metadata_disease_columns_JENSEN)
```




    array(['JENSEN_cancer', 'JENSEN_colon epithelial',
           'JENSEN_ulcerative colitis'],
          dtype='<U50')




```python
vec_user_sig = pd.Series(0, index=meta_data.columns[1:])
vec_user_sig[np.intersect1d(list(meta_data.columns), metadata_disease_columns_JENSEN)] = 1
vec_user_sig
```




    JENSEN_mucosal                                        0
    JENSEN_ileum                                          0
    JENSEN_colorectal cancer                              0
    JENSEN_epithelium                                     0
    JENSEN_cancer                                         1
    JENSEN_colon epithelial                               1
    JENSEN_inflammatory bowel disease                     0
    JENSEN_colonic epithelium                             0
    JENSEN_colitis                                        0
    JENSEN_crohn's disease                                0
    JENSEN_mucosa                                         0
    JENSEN_colonic mucosal                                0
    JENSEN_sigmoid colon                                  0
    JENSEN_ulcerative colitis                             1
    JENSEN_colonic mucosa                                 0
    JENSEN_pbmc                                           0
    JENSEN_bowel                                          0
    CREEDS_intestinal mucosa                              0
    CREEDS_colon                                          0
    CREEDS_colonic mucosa - non-inflamed                  0
    CREEDS_pbmc (peripheral blood mononuclear cells)      0
    CREEDS_colon structure (body structure)               0
    CREEDS_ascending colon                                0
    CREEDS_intestinal biopsies                            0
    CREEDS_intestines                                     0
    CREEDS_sigmoid colon                                  0
    CREEDS_descending colon                               0
    CREEDS_intestine - large intestine - colon (mmhcc)    0
    CREEDS_colonic mucosa - inflamed                      0
    CREEDS_interleukin-27 complex                         0
    CREEDS_sigmoid colons (mucosal biopsy)                0
    CREEDS_ileum                                          0
    CREEDS_interleukin-12 complex                         0
    CREEDS_body tissue                                    0
    CREEDS_mucous membrane                                0
    CREEDS_interleukin-23 complex                         0
    CREEDS_rectum                                         0
    CREEDS_terminal ileum                                 0
    CREEDS_distal part of ileum                           0
    CREEDS_colon mucosa                                   0
    CREEDS_peripheral blood mononuclear cell              0
    dtype: int64




```python
meta_data.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Unnamed: 0</th>
      <th>JENSEN_mucosal</th>
      <th>JENSEN_ileum</th>
      <th>JENSEN_colorectal cancer</th>
      <th>JENSEN_epithelium</th>
      <th>JENSEN_cancer</th>
      <th>JENSEN_colon epithelial</th>
      <th>JENSEN_inflammatory bowel disease</th>
      <th>JENSEN_colonic epithelium</th>
      <th>JENSEN_colitis</th>
      <th>...</th>
      <th>CREEDS_ileum</th>
      <th>CREEDS_interleukin-12 complex</th>
      <th>CREEDS_body tissue</th>
      <th>CREEDS_mucous membrane</th>
      <th>CREEDS_interleukin-23 complex</th>
      <th>CREEDS_rectum</th>
      <th>CREEDS_terminal ileum</th>
      <th>CREEDS_distal part of ileum</th>
      <th>CREEDS_colon mucosa</th>
      <th>CREEDS_peripheral blood mononuclear cell</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>dz:188</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>dz:249</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>dz:264</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>dz:454</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>dz:585</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 42 columns</p>
</div>




```python
meta_data_df = meta_data.set_index('Unnamed: 0')
```

```python
meta_data_df = meta_data_df.reset_index()
```

```python
user_ratings = pd.read_csv('manual_rating_creeds.csv')
```

```python
user_ratings.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Unnamed: 0</th>
      <th>GSE1009</th>
      <th>GSE11194</th>
      <th>GSE11223</th>
      <th>GSE11237</th>
      <th>GSE11618</th>
      <th>GSE12452</th>
      <th>GSE13760</th>
      <th>GSE1420</th>
      <th>GSE15102</th>
      <th>...</th>
      <th>GSE6475</th>
      <th>GSE65144</th>
      <th>GSE6688</th>
      <th>GSE6731</th>
      <th>GSE67561</th>
      <th>GSE6787</th>
      <th>GSE7657</th>
      <th>GSE7664</th>
      <th>GSE9128</th>
      <th>GSE9452</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>dz:454</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>2.0</td>
      <td>NaN</td>
      <td>5.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>dz:P2532</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>dz:P2531</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>5.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>dz:586</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>5.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2.0</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>5.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>5.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>dz:585</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>5.0</td>
      <td>2.0</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>5.0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 102 columns</p>
</div>




```python
user_ratings = user_ratings.set_index('Unnamed: 0')
```

```python
user_ratings.shape
```




    (9, 101)




```python
sig_dataset = user_ratings.index
```

```python
s= sig_dataset.tolist()
s
```




    ['dz:454',
     'dz:P2532',
     'dz:P2531',
     'dz:586',
     'dz:585',
     'dz:188',
     'dz:P442',
     'dz:264',
     'dz:249']




```python
metadata_df_new = meta_data_df[meta_data_df['Unnamed: 0'].isin(s)]
metadata_df_new
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Unnamed: 0</th>
      <th>JENSEN_mucosal</th>
      <th>JENSEN_ileum</th>
      <th>JENSEN_colorectal cancer</th>
      <th>JENSEN_epithelium</th>
      <th>JENSEN_cancer</th>
      <th>JENSEN_colon epithelial</th>
      <th>JENSEN_inflammatory bowel disease</th>
      <th>JENSEN_colonic epithelium</th>
      <th>JENSEN_colitis</th>
      <th>...</th>
      <th>CREEDS_ileum</th>
      <th>CREEDS_interleukin-12 complex</th>
      <th>CREEDS_body tissue</th>
      <th>CREEDS_mucous membrane</th>
      <th>CREEDS_interleukin-23 complex</th>
      <th>CREEDS_rectum</th>
      <th>CREEDS_terminal ileum</th>
      <th>CREEDS_distal part of ileum</th>
      <th>CREEDS_colon mucosa</th>
      <th>CREEDS_peripheral blood mononuclear cell</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>dz:188</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>dz:249</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>dz:264</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>dz:454</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>dz:585</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>5</th>
      <td>dz:586</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>30</th>
      <td>dz:P442</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>31</th>
      <td>dz:P2531</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>32</th>
      <td>dz:P2532</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>9 rows × 42 columns</p>
</div>




```python
metadata_df_new = metadata_df_new.set_index('Unnamed: 0')
```

```python
metadata_df_new
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>JENSEN_mucosal</th>
      <th>JENSEN_ileum</th>
      <th>JENSEN_colorectal cancer</th>
      <th>JENSEN_epithelium</th>
      <th>JENSEN_cancer</th>
      <th>JENSEN_colon epithelial</th>
      <th>JENSEN_inflammatory bowel disease</th>
      <th>JENSEN_colonic epithelium</th>
      <th>JENSEN_colitis</th>
      <th>JENSEN_crohn's disease</th>
      <th>...</th>
      <th>CREEDS_ileum</th>
      <th>CREEDS_interleukin-12 complex</th>
      <th>CREEDS_body tissue</th>
      <th>CREEDS_mucous membrane</th>
      <th>CREEDS_interleukin-23 complex</th>
      <th>CREEDS_rectum</th>
      <th>CREEDS_terminal ileum</th>
      <th>CREEDS_distal part of ileum</th>
      <th>CREEDS_colon mucosa</th>
      <th>CREEDS_peripheral blood mononuclear cell</th>
    </tr>
    <tr>
      <th>Unnamed: 0</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>dz:188</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>dz:249</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>dz:264</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>dz:454</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>dz:585</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>dz:586</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>dz:P442</th>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>dz:P2531</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>dz:P2532</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>9 rows × 41 columns</p>
</div>




```python
meta_data_df = metadata_df_new
```
### Cosine Similarity


```python
from scipy import spatial
```

```python
meta_data_df.head()
```




<div>
<style>
    .dataframe thead tr:only-child th {
        text-align: right;
    }

    .dataframe thead th {
        text-align: left;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>JENSEN_mucosal</th>
      <th>JENSEN_ileum</th>
      <th>JENSEN_colorectal cancer</th>
      <th>JENSEN_epithelium</th>
      <th>JENSEN_cancer</th>
      <th>JENSEN_colon epithelial</th>
      <th>JENSEN_inflammatory bowel disease</th>
      <th>JENSEN_colonic epithelium</th>
      <th>JENSEN_colitis</th>
      <th>JENSEN_crohn's disease</th>
      <th>...</th>
      <th>CREEDS_ileum</th>
      <th>CREEDS_interleukin-12 complex</th>
      <th>CREEDS_body tissue</th>
      <th>CREEDS_mucous membrane</th>
      <th>CREEDS_interleukin-23 complex</th>
      <th>CREEDS_rectum</th>
      <th>CREEDS_terminal ileum</th>
      <th>CREEDS_distal part of ileum</th>
      <th>CREEDS_colon mucosa</th>
      <th>CREEDS_peripheral blood mononuclear cell</th>
    </tr>
    <tr>
      <th>Unnamed: 0</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>dz:188</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>dz:249</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>dz:264</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>dz:454</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>dz:585</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>...</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 41 columns</p>
</div>



### Calculating similarity of all signatures with user signature


```python
list_of_signs = list(meta_data_df.index)

sim_score = []

user_sign = np.array(vec_user_sig)
for i in range(len(list_of_signs)):
    print(i)     
    similarity = 1 - spatial.distance.cosine(meta_data_df.iloc[i,:],np.array(vec_user_sig))
    sim_score.append(similarity)
    
```
    0
    1
    2
    3
    4
    5
    6
    7
    8



```python
sim_score
```




    [0.0,
     0.28867513459481287,
     0.2581988897471611,
     0.28867513459481287,
     0.66666666666666674,
     0.66666666666666674,
     0.19245008972987532,
     0.28867513459481287,
     0.28867513459481287]



### Creating a dict for signature and their similarity scores


```python
sim_sign_dict = {}

for i in range(len(sim_score)):
    sim_sign_dict[list_of_signs[i]] = sim_score[i]
```

```python
sim_sign_dict
```




    {'dz:188': 0.0,
     'dz:249': 0.28867513459481287,
     'dz:264': 0.2581988897471611,
     'dz:454': 0.28867513459481287,
     'dz:585': 0.66666666666666674,
     'dz:586': 0.66666666666666674,
     'dz:P442': 0.19245008972987532,
     'dz:P2531': 0.28867513459481287,
     'dz:P2532': 0.28867513459481287}



### Sorted Similarity


```python
def sorted_sim_score_dict(sim_scores):
    
    sorted_score = sorted(sim_sign_dict.items(), key=lambda value: value[1], reverse=True)
    return sorted_score

```

```python
sort_sim_score = sorted_sim_score_dict(sim_score)
sort_sim_score
```




    [('dz:585', 0.66666666666666674),
     ('dz:586', 0.66666666666666674),
     ('dz:249', 0.28867513459481287),
     ('dz:454', 0.28867513459481287),
     ('dz:P2531', 0.28867513459481287),
     ('dz:P2532', 0.28867513459481287),
     ('dz:264', 0.2581988897471611),
     ('dz:P442', 0.19245008972987532),
     ('dz:188', 0.0)]



### Top-K neighbors of similar signature


```python
def top_k_neighbors(sorted_score, k):
    
    top_k = sorted_score[0:k]
    
    return top_k
```

```python
top_k = top_k_neighbors(sort_sim_score, 3)
top_k
```




    [('dz:585', 0.66666666666666674),
     ('dz:586', 0.66666666666666674),
     ('dz:249', 0.28867513459481287)]



### Filtering from Top-k


```python
def filter_top_k(top_k):
    
    filtered_top_k =[]
    
    for i in range(len(top_k)):
        
        if top_k[i][1] > 0.1:
            
            filtered_top_k.append(top_k[i])
            
    return filtered_top_k
```

```python
filtered_top_k = filter_top_k(top_k)
filtered_top_k
```




    [('dz:585', 0.66666666666666674),
     ('dz:586', 0.66666666666666674),
     ('dz:249', 0.28867513459481287)]



### Extracting the signature and the similarity scores


```python
sim_sig = [x[0] for x in filtered_top_k]
sim_scr = [x[1] for x in filtered_top_k]
```

```python
sim_sig
```




    ['dz:585', 'dz:586', 'dz:249']




```python
sim_scr
```




    [0.66666666666666674, 0.66666666666666674, 0.28867513459481287]



### Dataset scores


```python
dataset_names = user_ratings.columns
```

```python
dataset_names
```




    Index(['GSE1009', 'GSE11194', 'GSE11223', 'GSE11237', 'GSE11618', 'GSE12452',
           'GSE13760', 'GSE1420', 'GSE15102', 'GSE15471',
           ...
           'GSE6475', 'GSE65144', 'GSE6688', 'GSE6731', 'GSE67561', 'GSE6787',
           'GSE7657', 'GSE7664', 'GSE9128', 'GSE9452'],
          dtype='object', length=101)




```python
len(dataset_names)
```




    101




```python
len(dataset_names.unique())
```




    101




```python
total_score_list = []

for i in range(0, len(user_ratings.columns)):
    score = []
    
    for j in range(len(sim_sig)):
        rating = user_ratings.loc[sim_sig[j]][dataset_names[i]] 
        score.append(rating*sim_scr[j])
        
    total_score = np.nansum(score)/sum(sim_scr)
    total_score_list.append(total_score)
    print('{0}:{1}'.format(dataset_names[i], total_score))
    
```
    GSE1009:0.0
    GSE11194:0.0
    GSE11223:4.110130617987622
    GSE11237:0.8220261235975244
    GSE11618:0.0
    GSE12452:0.0
    GSE13760:0.0
    GSE1420:0.8220261235975244
    GSE15102:0.3559477528049514
    GSE15471:0.8220261235975244
    GSE15481:1.2330391853962863
    GSE15568:0.0
    GSE15811:0.0
    GSE1629:2.055065308993811
    GSE16464:0.0
    GSE16515:0.8220261235975244
    GSE1676:0.0
    GSE16962:0.0
    GSE17025:0.0
    GSE18560:0.0
    GSE18567:0.0
    GSE1919:0.8220261235975244
    GSE19650:0.0
    GSE19675:0.0
    GSE1987:0.0
    GSE2077:0.0
    GSE21329:0.0
    GSE21422:0.0
    GSE22366:0.0
    GSE22606:0.0
    GSE22619:2.055065308993811
    GSE23031:0.0
    GSE24514:0.0
    GSE2503:0.0
    GSE26104:1.6440522471950487
    GSE26299:0.0
    GSE26817:1.6440522471950487
    GSE2684:0.0
    GSE27318:0.0
    GSE28185:0.0
    GSE29110:0.3559477528049514
    GSE29145:0.0
    GSE29691:0.0
    GSE3112:1.6440522471950487
    GSE3140:0.0
    GSE31432:0.533921629207427
    GSE32323:0.533921629207427
    GSE32924:0.533921629207427
    GSE3365:0.0
    GSE33709:0.0
    GSE34305:0.0
    GSE34526:0.0
    GSE34619:0.0
    GSE3467:0.533921629207427
    GSE34860:0.0
    GSE35543:0.7118955056099028
    GSE3624:1.2330391853962863
    GSE36700:0.0
    GSE3744:0.0
    GSE3868:0.0
    GSE38713:0.8898693820123784
    GSE40207:0.0
    GSE4183:2.1229085674086647
    GSE43696:1.6440522471950487
    GSE43741:0.0
    GSE44590:0.0
    GSE45452:0.0
    GSE4587:0.0
    GSE46448:0.0
    GSE46449:0.0
    GSE47685:0.8898693820123784
    GSE47751:0.8220261235975244
    GSE48466:1.6440522471950487
    GSE49153:0.0
    GSE5007:0.0
    GSE51808:0.0
    GSE52471:0.0
    GSE5339:0.0
    GSE53431:0.0
    GSE54645:0.0
    GSE5504:0.0
    GSE5550:0.0
    GSE5681:0.0
    GSE5788:0.0
    GSE6012:0.0
    GSE6088:0.0
    GSE61304:0.0
    GSE62632:0.0
    GSE6281:0.3559477528049514
    GSE6400:0.0
    GSE6443:0.7118955056099028
    GSE6475:0.533921629207427
    GSE65144:0.0
    GSE6688:0.0
    GSE6731:2.9449346910061895
    GSE67561:0.0
    GSE6787:0.0
    GSE7657:0.0
    GSE7664:0.0
    GSE9128:0.0
    GSE9452:5.0


### Creating a dict for dataset and there scores


```python
dataset_score_dict = {}
for i in range(len(user_ratings.columns)):
    #print(str(i+1))
    dataset_score_dict[dataset_names[i]] = total_score_list[i]
```

```python
dataset_score_dict
```




    {'GSE1009': 0.0,
     'GSE11194': 0.0,
     'GSE11223': 4.1101306179876218,
     'GSE11237': 0.82202612359752436,
     'GSE11618': 0.0,
     'GSE12452': 0.0,
     'GSE13760': 0.0,
     'GSE1420': 0.82202612359752436,
     'GSE15102': 0.35594775280495139,
     'GSE15471': 0.82202612359752436,
     'GSE15481': 1.2330391853962863,
     'GSE15568': 0.0,
     'GSE15811': 0.0,
     'GSE1629': 2.0550653089938109,
     'GSE16464': 0.0,
     'GSE16515': 0.82202612359752436,
     'GSE1676': 0.0,
     'GSE16962': 0.0,
     'GSE17025': 0.0,
     'GSE18560': 0.0,
     'GSE18567': 0.0,
     'GSE1919': 0.82202612359752436,
     'GSE19650': 0.0,
     'GSE19675': 0.0,
     'GSE1987': 0.0,
     'GSE2077': 0.0,
     'GSE21329': 0.0,
     'GSE21422': 0.0,
     'GSE22366': 0.0,
     'GSE22606': 0.0,
     'GSE22619': 2.0550653089938109,
     'GSE23031': 0.0,
     'GSE24514': 0.0,
     'GSE2503': 0.0,
     'GSE26104': 1.6440522471950487,
     'GSE26299': 0.0,
     'GSE26817': 1.6440522471950487,
     'GSE2684': 0.0,
     'GSE27318': 0.0,
     'GSE28185': 0.0,
     'GSE29110': 0.35594775280495139,
     'GSE29145': 0.0,
     'GSE29691': 0.0,
     'GSE3112': 1.6440522471950487,
     'GSE3140': 0.0,
     'GSE31432': 0.53392162920742703,
     'GSE32323': 0.53392162920742703,
     'GSE32924': 0.53392162920742703,
     'GSE3365': 0.0,
     'GSE33709': 0.0,
     'GSE34305': 0.0,
     'GSE34526': 0.0,
     'GSE34619': 0.0,
     'GSE3467': 0.53392162920742703,
     'GSE34860': 0.0,
     'GSE35543': 0.71189550560990278,
     'GSE3624': 1.2330391853962863,
     'GSE36700': 0.0,
     'GSE3744': 0.0,
     'GSE3868': 0.0,
     'GSE38713': 0.88986938201237842,
     'GSE40207': 0.0,
     'GSE4183': 2.1229085674086647,
     'GSE43696': 1.6440522471950487,
     'GSE43741': 0.0,
     'GSE44590': 0.0,
     'GSE45452': 0.0,
     'GSE4587': 0.0,
     'GSE46448': 0.0,
     'GSE46449': 0.0,
     'GSE47685': 0.88986938201237842,
     'GSE47751': 0.82202612359752436,
     'GSE48466': 1.6440522471950487,
     'GSE49153': 0.0,
     'GSE5007': 0.0,
     'GSE51808': 0.0,
     'GSE52471': 0.0,
     'GSE5339': 0.0,
     'GSE53431': 0.0,
     'GSE54645': 0.0,
     'GSE5504': 0.0,
     'GSE5550': 0.0,
     'GSE5681': 0.0,
     'GSE5788': 0.0,
     'GSE6012': 0.0,
     'GSE6088': 0.0,
     'GSE61304': 0.0,
     'GSE62632': 0.0,
     'GSE6281': 0.35594775280495139,
     'GSE6400': 0.0,
     'GSE6443': 0.71189550560990278,
     'GSE6475': 0.53392162920742703,
     'GSE65144': 0.0,
     'GSE6688': 0.0,
     'GSE6731': 2.9449346910061895,
     'GSE67561': 0.0,
     'GSE6787': 0.0,
     'GSE7657': 0.0,
     'GSE7664': 0.0,
     'GSE9128': 0.0,
     'GSE9452': 5.0}




```python
total_score_list
```




    [0.0,
     0.0,
     4.1101306179876218,
     0.82202612359752436,
     0.0,
     0.0,
     0.0,
     0.82202612359752436,
     0.35594775280495139,
     0.82202612359752436,
     1.2330391853962863,
     0.0,
     0.0,
     2.0550653089938109,
     0.0,
     0.82202612359752436,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.82202612359752436,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     2.0550653089938109,
     0.0,
     0.0,
     0.0,
     1.6440522471950487,
     0.0,
     1.6440522471950487,
     0.0,
     0.0,
     0.0,
     0.35594775280495139,
     0.0,
     0.0,
     1.6440522471950487,
     0.0,
     0.53392162920742703,
     0.53392162920742703,
     0.53392162920742703,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.53392162920742703,
     0.0,
     0.71189550560990278,
     1.2330391853962863,
     0.0,
     0.0,
     0.0,
     0.88986938201237842,
     0.0,
     2.1229085674086647,
     1.6440522471950487,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.88986938201237842,
     0.82202612359752436,
     1.6440522471950487,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     0.35594775280495139,
     0.0,
     0.71189550560990278,
     0.53392162920742703,
     0.0,
     0.0,
     2.9449346910061895,
     0.0,
     0.0,
     0.0,
     0.0,
     0.0,
     5.0]



### Sorting the dataset with their scores


```python
def sorted_dataset_score_dict(dataset_scores):
    
    sorted_dataset_score = sorted(dataset_scores.items(), key=lambda value: value[1], reverse=True)
    return sorted_dataset_score
```

```python
sorted_score = sorted_dataset_score_dict(dataset_score_dict)
```

```python
sorted_score
```




    [('GSE9452', 5.0),
     ('GSE11223', 4.1101306179876218),
     ('GSE6731', 2.9449346910061895),
     ('GSE4183', 2.1229085674086647),
     ('GSE1629', 2.0550653089938109),
     ('GSE22619', 2.0550653089938109),
     ('GSE26104', 1.6440522471950487),
     ('GSE26817', 1.6440522471950487),
     ('GSE3112', 1.6440522471950487),
     ('GSE43696', 1.6440522471950487),
     ('GSE48466', 1.6440522471950487),
     ('GSE15481', 1.2330391853962863),
     ('GSE3624', 1.2330391853962863),
     ('GSE38713', 0.88986938201237842),
     ('GSE47685', 0.88986938201237842),
     ('GSE11237', 0.82202612359752436),
     ('GSE1420', 0.82202612359752436),
     ('GSE15471', 0.82202612359752436),
     ('GSE16515', 0.82202612359752436),
     ('GSE1919', 0.82202612359752436),
     ('GSE47751', 0.82202612359752436),
     ('GSE35543', 0.71189550560990278),
     ('GSE6443', 0.71189550560990278),
     ('GSE31432', 0.53392162920742703),
     ('GSE32323', 0.53392162920742703),
     ('GSE32924', 0.53392162920742703),
     ('GSE3467', 0.53392162920742703),
     ('GSE6475', 0.53392162920742703),
     ('GSE15102', 0.35594775280495139),
     ('GSE29110', 0.35594775280495139),
     ('GSE6281', 0.35594775280495139),
     ('GSE1009', 0.0),
     ('GSE11194', 0.0),
     ('GSE11618', 0.0),
     ('GSE12452', 0.0),
     ('GSE13760', 0.0),
     ('GSE15568', 0.0),
     ('GSE15811', 0.0),
     ('GSE16464', 0.0),
     ('GSE1676', 0.0),
     ('GSE16962', 0.0),
     ('GSE17025', 0.0),
     ('GSE18560', 0.0),
     ('GSE18567', 0.0),
     ('GSE19650', 0.0),
     ('GSE19675', 0.0),
     ('GSE1987', 0.0),
     ('GSE2077', 0.0),
     ('GSE21329', 0.0),
     ('GSE21422', 0.0),
     ('GSE22366', 0.0),
     ('GSE22606', 0.0),
     ('GSE23031', 0.0),
     ('GSE24514', 0.0),
     ('GSE2503', 0.0),
     ('GSE26299', 0.0),
     ('GSE2684', 0.0),
     ('GSE27318', 0.0),
     ('GSE28185', 0.0),
     ('GSE29145', 0.0),
     ('GSE29691', 0.0),
     ('GSE3140', 0.0),
     ('GSE3365', 0.0),
     ('GSE33709', 0.0),
     ('GSE34305', 0.0),
     ('GSE34526', 0.0),
     ('GSE34619', 0.0),
     ('GSE34860', 0.0),
     ('GSE36700', 0.0),
     ('GSE3744', 0.0),
     ('GSE3868', 0.0),
     ('GSE40207', 0.0),
     ('GSE43741', 0.0),
     ('GSE44590', 0.0),
     ('GSE45452', 0.0),
     ('GSE4587', 0.0),
     ('GSE46448', 0.0),
     ('GSE46449', 0.0),
     ('GSE49153', 0.0),
     ('GSE5007', 0.0),
     ('GSE51808', 0.0),
     ('GSE52471', 0.0),
     ('GSE5339', 0.0),
     ('GSE53431', 0.0),
     ('GSE54645', 0.0),
     ('GSE5504', 0.0),
     ('GSE5550', 0.0),
     ('GSE5681', 0.0),
     ('GSE5788', 0.0),
     ('GSE6012', 0.0),
     ('GSE6088', 0.0),
     ('GSE61304', 0.0),
     ('GSE62632', 0.0),
     ('GSE6400', 0.0),
     ('GSE65144', 0.0),
     ('GSE6688', 0.0),
     ('GSE67561', 0.0),
     ('GSE6787', 0.0),
     ('GSE7657', 0.0),
     ('GSE7664', 0.0),
     ('GSE9128', 0.0)]



### Selecting the top n datasets


```python
def top_n_dataset(sorted_score, k):
    top_k = sorted_score[0:k]
    
    return top_k
```
### Top n datasets recommended 


```python
top_n_dataset = top_n_dataset(sorted_score, 3)
```

```python
top_n_dataset
```




    [('GSE9452', 5.0),
     ('GSE11223', 4.1101306179876218),
     ('GSE6731', 2.9449346910061895)]




```python
get_summary_and_title('GSE9452')
```




    {'gse_id': 'GSE9452',
     'title': 'Definition of an ulcerative colitis preinflammatory state',
     'exp_type': 'Expression profiling by array',
     'abstract': 'The samples are a part of a study aiming at diagnosing ulcerative colitis from genome-wide gene expression analysis of the colonic mucosa. Colonic mucosal samples were collected as endoscopic pinch biopsies from ulcerative colitis patients and from control subjects. Samples with and without macroscopic signs of inflammation were collected from the patients.Keywords: Disease state analysis'}




```python
get_summary_and_title('GSE11223')
```




    {'gse_id': 'GSE11223',
     'title': 'Colon biopsies from UC patients and healthy controls',
     'exp_type': 'Expression profiling by array',
     'abstract': 'Transcriptional profiling of colon epithelial biopsies from ulcerative colitis patients and healthy control donors.Study aims to survey and analyze variation from disease in different GI regions.Keywords: disease state analysis'}




```python
get_summary_and_title('GSE6731')
```




    {'gse_id': 'GSE6731',
     'title': 'Genome-wide gene expression differences between Crohnâ\x80\x99s and ulcerative colitis from endoscopic pinch biopsies: ',
     'exp_type': 'Expression profiling by array',
     'abstract': 'Ulcerative colitis (UC) and Crohnâ\x80\x99s disease (CD) are inflammatory bowel diseases (IBD) with variable, overlapping clinical features and complex pathophysiologies. To identify pathogenic processes underlying these disease subtypes, using single endoscopic pinch biopsies to estabolish 36 expression profiles, we elucidated gene expression patterns of active and inactive areas of UC and CD, and compared these to infectious colitis and healthy controls. Keywords: RNA'}
