import multiprocessing as mp
import Bio
from Bio import Entrez
import argparse as ap
from pathlib import Path
import os

def get_references(pmid):
    Entrez.email = 'badmusanisat@gmail.com'  
    results = Entrez.read(Entrez.elink(dbfrom="pubmed",
                                    db="pmc",
                                    LinkName="pubmed_pmc_refs",
                                    id=pmid,
                                    api_key='8fa896ca3cd1a5e694493b053a03429e4d08'))
    references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]
    return references
    #print(references)
    


def download_article(pmid):
    Entrez.email = 'badmusanisat@gmail.com'  
    handle = Entrez.efetch(db="pmc", id=pmid, rettype="XML", retmode="text",
                            api_key='8fa896ca3cd1a5e694493b053a03429e4d08')
    
    with open(f'/Users/anisah/Programming3/Programming3/Assignment1/output/{pmid}.xml', 'wb') as file:
        file.write(handle.read())



if __name__ == '__main__':
    argparser = ap.ArgumentParser(description="Script that downloads (default) 10 articles referenced by the given PubMed ID concurrently.")
    argparser.add_argument("-n", action="store",
                           dest="n", required=False, type=int,
                           help="Number of references to download concurrently.")
    argparser.add_argument("pubmed_id", action="store", type=str, nargs=1, help="Pubmed ID of the article to harvest for references to download.")
    args = argparser.parse_args()
    print("Getting: ", args.pubmed_id)
    pmid =  str(args.pubmed_id)
    print(pmid)
    
    references = get_references(pmid)
    print(references)
    
    
    with mp.Pool(4) as pool:
        results = pool.map(download_article, references[0:10])
    



