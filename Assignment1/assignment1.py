import multiprocessing as mp
from Bio import Entrez
import argparse as ap
import time
from pathlib import Path

def get_references(pmid):
    '''function to get referenced documents 
    input: pmid
    output: list of referenced documents ids'''
    Entrez.email = 'badmusanisat@gmail.com'  
    results = Entrez.read(Entrez.elink(dbfrom="pubmed",
                                    db="pmc",
                                    LinkName="pubmed_pmc_refs",
                                    id=pmid,
                                    api_key='8fa896ca3cd1a5e694493b053a03429e4d08'))
    references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]
    return references

def download_article(pmid):
    '''function to download an article
    input: pmid
    output: article in xml format'''
    Entrez.email = 'badmusanisat@gmail.com'  
    handle = Entrez.efetch(db="pmc", id=pmid, rettype="XML", retmode="text",
                            api_key='8fa896ca3cd1a5e694493b053a03429e4d08')
    path = Path(__file__).parent.absolute()
    output_path = path/'output'
    with open(f'{output_path}/{pmid}.xml', 'wb') as file:
        file.write(handle.read())



if __name__ == '__main__':
  
    argparser = ap.ArgumentParser(description="Script that downloads (default) 10 articles referenced by the given PubMed ID concurrently.")
    argparser.add_argument("-n", action="store",
                           dest="n", required=False, type=int,
                           help="Number of references to download concurrently.")
    argparser.add_argument("pubmed_id", action="store", type=str, nargs=1, help="Pubmed ID of the article to harvest for references to download.")
    args = argparser.parse_args()
    print("Getting: ", args.pubmed_id)
    
    #print ("Start : %s" % time.ctime())
    
    #get pmid from command
    pmid =  str(args.pubmed_id)
    #print(pmid)
    
    # get articles referenced in the article with the given pmid
    references = get_references(pmid)
    #print(references)
    
    #download articles using multiprocessing
    cpu = int(mp.cpu_count() / 2)
    with mp.Pool(cpu) as pool:
        results = pool.map(download_article, references[0:10])
    
    #print ("End : %s" % time.ctime() )



