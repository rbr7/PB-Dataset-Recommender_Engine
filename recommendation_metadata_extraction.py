# import requests
# from html.parser import HTMLParser
# from bs4 import BeautifulSoup
# import pandas as pd
# import numpy as np
# import becas
# import pickle
# import json
# import copy
# import pymongo
from sys import argv
import logging

log_file = 'Overall_gse_ids_log_file.log'
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - level=logging.INFO %(message)s')

script, gse_id = argv

gse_id = argv[1]

count = 0

try:
    # getting the summary/abstract

    import get_summary_or_abstract

    dataset_meta = get_summary_or_abstract.get_summary_and_title(gse_id)

    count+=1

    # getting the metadata from becas

    import calling_becas_api

    metadata_becas = calling_becas_api.call_becas_api(dataset_meta)

    count+=1

    # processing the data
    import processing_metadata

    metadata_becas_processed = processing_metadata.metadata_processing(metadata_becas)

    ### Feeding the data in mongodb
    import feed_mongo

    doc_feeding = feed_mongo.insert_mongo_doc(metadata_becas_processed)

    print('congrats! all done...\n')
    logging.info(gse_id + ' ' + 'SUCCESS')
except:
    logging.info(gse_id + ' ' + 'FAILURE')
    
