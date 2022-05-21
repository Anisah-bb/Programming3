import multiprocessing as mp
from Bio import Entrez
import argparse as ap
import sys
import time
from pathlib import Path
import pickle

#pmid = 30049270
def make_directory():
    '''function to create out put directory
    output: a output folder
    '''
    path = Path(__file__).parent.absolute()
    output_path = path/'output'
    try:
        output_path.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        print("Ouput Folder already exists")
    return output_path

def get_references(pmid):
    '''function to get referenced documents 
    input: pmid
    output: list of referenced documents ids
    '''
    Entrez.email = 'badmusanisat@gmail.com'  
    results = Entrez.read(Entrez.elink(dbfrom="pubmed",
                                    db="pmc",
                                    LinkName="pubmed_pmc_refs",
                                    id=pmid,
                                    api_key='8fa896ca3cd1a5e694493b053a03429e4d08'))
    references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]
    return references

def get_authors(pmid):
    '''funtion to get authors of the referenced articles
    '''
    Entrez.email = 'badmusanisat@gmail.com' 
    results = Entrez.read(Entrez.esummary(db="pubmed",
                                        id=pmid,
                                        api_key='8fa896ca3cd1a5e694493b053a03429e4d08'))
    author_list = []
    author_list = [author for author  in results[0]["AuthorList"]]
    #print(f"Author list: {author_list}")
    author_tup = tuple(author_list)
    return author_tup

def write_pickle(pmid):
    '''function to write out the authors list as a tuple into a pickle 
    file
    input: pmid
    output: file(s) containing authors as a tuple
    '''
    output_path = make_directory()
    author_tup = get_authors(pmid)
    with open(f'{output_path}/{pmid}.authors.pickle', 'wb') as file:
        pickle.dump(author_tup, file)


if __name__ == '__main__':
    #pmid = 30049270 # 8767730
    pmid = sys.argv[-1]
    references = get_references(pmid)
    for id in references:
        write_pickle(id)
    
   
    

    


