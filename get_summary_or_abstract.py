def get_summary_and_title(gse_id):
    
    import requests
    from html.parser import HTMLParser
    from bs4 import BeautifulSoup
    import logging
    
    GEO_URL = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
    
#     Log_File = 'GSE_log_file.log'
    
#     logging.basicConfig(filename=Log_File, level=logging.INFO, format='%(asctime)s - %(levelname)s %(message)s')
    
    url = GEO_URL + gse_id

    page = requests.get(url)
    soup = BeautifulSoup(page.text, "html.parser")
    title = soup.find("td", text="Title").find_next_sibling("td").text
    # pltfm_org = soup.find("td ", text="Platform organism").find_next_sibling("td").text
    # sample_org = soup.find("td", text="Sample organism").find_next_sibling("td").text
    exp_type = soup.find("td", text="Experiment type").find_next_sibling("td").text
    summary = soup.find("td", text="Summary").find_next_sibling("td").text
    PMID_Probe = soup.find("td", text="Citation(s)")
    if PMID_Probe!=None:
        PMID = soup.find("td", text="Citation(s)").find_next_sibling("td").text
    else:
        PMID = None

    dataset_meta = {'gse_id': gse_id, 'title': title, 'exp_type': exp_type, 'summary': summary, 'PMID':PMID}
    print('got the summary and title for {}'.format(gse_id))
    #logging.info('got the summary and title for' + ' ' + gse_id)
    
    return(dataset_meta)