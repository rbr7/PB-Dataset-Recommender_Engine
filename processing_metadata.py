count = 0

def metadata_processing(metadata_becas):
    
    import copy
    import logging
    
    metadata_becas_processed = copy.deepcopy(metadata_becas)
    gse = list(metadata_becas_processed.keys())[0]
    
    ent_processed = {}
    
    for ent in metadata_becas_processed[gse]['entities']:
        
        name_source_meta = ent.split('|')[0:2]
        
        ent_processed[name_source_meta[0]] = name_source_meta[1].split(';')
   
    ent_further_processed = {}
    
    for s in list(ent_processed.keys()):
        
        if '.'.lower() in s.lower():
            
            new_name = s.replace('.','')
            ent_further_processed[new_name] = ent_processed[s]
            
        else:
            ent_further_processed[s] = ent_processed[s]
    
    if len(ent_processed) != len(ent_further_processed):
        print('something fishy with this {}'.format(gse))
    
    metadata_becas_processed[gse]['entities'] = ent_further_processed
    
    ids_processed = {}
    
    for ident in metadata_becas_processed[gse]['ids'].keys():
        name_ref = metadata_becas_processed[gse]['ids'][ident]
        
        ref_id = {}

        try: 
            ref_id['reference'] = metadata_becas_processed[gse]['ids'][ident]['refs']
            ref_id['identity'] = ident
        
        
            ids_processed[metadata_becas_processed[gse]['ids'][ident]['name']] = ref_id
        except TypeError:
            print('NoneType object is not subscriptable')
            #logging.exception("Exception occurred")
            
    ids_further_processed = {}
    
    for s in list(ids_processed.keys()):
        
        if '.'.lower() in s.lower():
            
            new_name = s.replace('.','')
            
            ids_further_processed[new_name] = ids_processed[s]
            
        else:
            ids_further_processed[s] = ids_processed[s]
           
    metadata_becas_processed[gse]['ids'] = ids_further_processed
    
    print('metadata for {} processed'.format(gse))
    #logging.info('metadata for' ' ' + gse + ' ' + 'processed')
    
    return metadata_becas_processed

