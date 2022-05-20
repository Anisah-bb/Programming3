import multiprocessing as mp
from Bio import Entrez
import argparse as ap
import time
from pathlib import Path
import pickle

#pmid = 30049270
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
def get_authors(pmid):
    '''funtion to get authors of each xml file
    '''
    Entrez.email = 'badmusanisat@gmail.com' 
    results = Entrez.read(Entrez.esummary(dbfrom="pubmed",
                                        db="pmc",
                                        LinkName="pubmed_pmc_refs",
                                        id=pmid,
                                        api_key='8fa896ca3cd1a5e694493b053a03429e4d08'))
    author_list = []
    author_list = [author for author  in results[0]["AuthorList"]]
    author_tup = tuple(author_list)
    path = Path(__file__).parent.absolute()
    output_path = path/'output'
    with open(f'{output_path}/{pmid}.authors.pickle', 'wb') as file:
        pickle.dump(author_tup, file)


if __name__ == '__main__':
  
    # argparser = ap.ArgumentParser(description="Script that downloads (default) 10 articles referenced by the given PubMed ID concurrently.")
    # argparser.add_argument("-n", action="store",
    #                        dest="n", required=False, type=int,
    #                        help="Number of references to download concurrently.")
    # argparser.add_argument("pubmed_id", action="store", type=str, nargs=1, help="Pubmed ID of the article to harvest for references to download.")
    # args = argparser.parse_args()
    # print("Getting: ", args.pubmed_id)
    
    #print ("Start : %s" % time.ctime())
    
    #get pmid from command
    pmid =  8767730
    #print(pmid)
    
    # get articles referenced in the article with the given pmid
    references = get_references(pmid)
    #print(references)
    get_authors(pmid)
    


