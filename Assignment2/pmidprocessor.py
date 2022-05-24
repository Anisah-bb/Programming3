import argparse as ap
from pathlib import Path
import pickle
from Bio import Entrez


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

def get_authors(pmid):
    '''funtion to get authors of the referenced articles
    '''
    try:
        Entrez.email = 'badmusanisat@gmail.com'
        results = Entrez.read(Entrez.esummary(db="pubmed",id=pmid,
                                            api_key='8fa896ca3cd1a5e694493b053a03429e4d08'))
        author_list = []
        author_list = [author for author  in results[0]["AuthorList"]]
        #print(f"Author list: {author_list}")
        author_tup = tuple(author_list)
    except RuntimeError:
        author_tup = (None)
    return author_tup

def write_pickle(pmid):
    '''function to write out the authors list as a tuple into a pickle
    file
    input: pmid
    output: file(s) containing authors as a tuple
    '''
    output_path = make_directory()
    author_tup = get_authors(pmid)
    download_article(pmid)
    with open(f'{output_path}/{pmid}.authors.pickle', 'wb') as file:
        pickle.dump(author_tup, file)

if __name__ == '__main__':
    #pmid = 30049270 # 8767730
    argparser = ap.ArgumentParser(
                                description="Script that downloads (default) 10 articles referenced by the given PubMed ID concurrently.")
    argparser.add_argument("-n", action="store",dest="n", required=False, type=int,
                           help="Number of references to download concurrently.")
    argparser.add_argument("pubmed_id", action="store", type=str, nargs=1,
                             help="Pubmed ID of the article to harvest for references to download.")
    args = argparser.parse_args()
    print("Getting: ", args.pubmed_id)
    pmid =  str(args.pubmed_id)
    references = get_references(pmid)
    for pmid in references:
        write_pickle(pmid)
