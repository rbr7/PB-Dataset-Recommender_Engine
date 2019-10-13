def insert_mongo_doc(metadata_becas_processed):
    
    import pymongo
    import logging
    
#     log_file = 'Overall_gse_ids_log_file.log'
#     logging.basicConfig(filename=log_file, level=logging.DEBUG, format='%(asctime)s - level=logging.INFO %(message)s')
    
    myclient = pymongo.MongoClient("mongodb://localhost:27017/")
    
    # creating the database
    gsedb = myclient["Metadata_becas_gse"]
    
    # creating the collection
    gsecol = gsedb['gse_ids']

    count= 0
    
    metadata_further_processed = {}
    
    gse = list(metadata_becas_processed.keys())[0]
    
    metadata_further_processed['gse_id'] = gse
    metadata_further_processed['entities'] = metadata_becas_processed[gse]['entities']
    metadata_further_processed['ids'] = metadata_becas_processed[gse]['ids']

    x = gsecol.insert_one(metadata_further_processed)
    
    count += 1
    
    print('Fed {} document for {}'.format(count, gse))
