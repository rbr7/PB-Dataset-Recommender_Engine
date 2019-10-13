def call_becas_api(geo_info):
    
    import becas
    import logging
    import requests
    
    becas.email = 'rb7.itr7@gmail.com'
    
    metadata_becas ={}
    invalid_geo_ids = []
    
    geo_id = geo_info['gse_id']
    
    #logging.info(geo_id)
    
    count = 0
    summary_count = 0
    abstract_count = 0
    
    summary = geo_info['summary']
    pmid = geo_info['PMID']
    
    if pmid == None:
        results = becas.annotate_text(summary)
        summary_count += 1
        print('No PMID for {}'.format(geo_id))
        print('Extracting summary')
        
        #logging.info('No PMID for' + ' ' + geo_id + ' so extracting summary')
        
    else:
        url_test = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id="+geo_info['PMID']+"&retmode=json&rettype=abstract"
                                         
        Response_test = requests.get(url_test)
        abstract = Response_test.text.replace("\n", "")
        
        try:
            results = becas.annotate_text(abstract)
            print('Extracting abstract')
            #logging.info('PMID available for' ' ' + geo_id + ' ' + 'so extracting abstract')
        except becas.BecasException:
            invalid_geo_ids.append(geo_id)
            summary_count-= 1
            print('500 Server Error: Internal Server Error for BecasException for {}'.format(geo_ids))
            #logging.exception("Exception occurred")
        
        abstract_count += 1
    
    count+= 1
    print('Fetched info of {} {} id'.format(count, geo_id))
    
    results = becas.annotate_text(summary)
    
    metadata_becas[geo_id] = results
     
    return metadata_becas
