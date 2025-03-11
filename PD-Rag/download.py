import requests
import os
import pandas as pd
import time
import re
from bs4 import BeautifulSoup
import random
from Bio import Entrez

# Configuration
OUTPUT_FOLDER = "downloaded_papers"
ARTICLES_CSV = "articles.csv"
SEARCH_TERM = "parkinson's disease"
MAX_RESULTS = 300  # Number of papers to search for
SEARCH_SORT = "date"  # Sort by publication date (most recent first)

# Set your email for NCBI/Entrez API (required by NCBI)
Entrez.email = "your_email@example.com"  # CHANGE THIS TO YOUR EMAIL

# User agents to make requests appear more natural
USER_AGENTS = [
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
    'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/14.1.1 Safari/605.1.15',
    'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/92.0.4515.107 Safari/537.36'
]

def get_random_user_agent():
    """Return a random user agent."""
    return random.choice(USER_AGENTS)

def clean_filename(filename):
    """Remove invalid characters from filename."""
    return re.sub(r'[\\/*?:"<>|]', "", filename)

def search_pubmed_central():
    """
    Search PubMed for papers on Parkinson's disease, check which ones are in PMC,
    and create a CSV file with the available papers.
    """
    print(f"Searching PubMed for papers on '{SEARCH_TERM}'...")
    
    # First, search PubMed for papers matching our criteria
    handle = Entrez.esearch(
        db="pubmed",
        term=SEARCH_TERM,
        retmax=MAX_RESULTS,
        sort=SEARCH_SORT
    )
    record = Entrez.read(handle)
    handle.close()
    
    pmid_list = record["IdList"]
    print(f"Found {len(pmid_list)} papers in PubMed.")
    
    # Check which papers are available in PubMed Central
    pmc_papers = []
    
    for i, pmid in enumerate(pmid_list):
        print(f"Checking paper {i+1}/{len(pmid_list)} (PMID: {pmid})...")
        
        try:
            # Get detailed information about the paper
            handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
            record = Entrez.read(handle)
            handle.close()
            
            if 'PubmedArticle' not in record:
                print(f"  No PubMed article data for PMID {pmid}")
                continue
            
            article_data = record['PubmedArticle'][0]
            
            # Extract article details
            article_meta = article_data.get('MedlineCitation', {}).get('Article', {})
            
            # Get title
            title = article_meta.get('ArticleTitle', 'Unknown Title')
            
            # Get author information
            author_list = article_meta.get('AuthorList', [])
            
            if author_list and 'LastName' in author_list[0]:
                first_author = author_list[0]['LastName']
                if len(author_list) > 1:
                    author_display = f"{first_author}_et_al"
                else:
                    author_display = first_author
            else:
                author_display = "Unknown"
            
            # Get publication year
            pub_date = article_meta.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
            pub_year = pub_date.get('Year', '0000')
            
            # Get DOI
            doi = None
            if 'ELocationID' in article_meta:
                for location in article_meta['ELocationID']:
                    if location.attributes.get('EIdType') == 'doi':
                        doi = str(location)
                        break
            
            if not doi:
                # Try to get DOI from ArticleIdList
                article_ids = article_data.get('PubmedData', {}).get('ArticleIdList', [])
                for article_id in article_ids:
                    if article_id.attributes.get('IdType') == 'doi':
                        doi = str(article_id)
                        break
            
            # Check if in PMC by searching for links to PMC
            pmc_id = None
            
            # Check if it's in PMC using ELink
            handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
            link_results = Entrez.read(handle)
            handle.close()
            
            # Check if we have PMC links
            if link_results and link_results[0].get('LinkSetDb'):
                links = link_results[0]['LinkSetDb']
                for link_db in links:
                    if link_db.get('DbTo') == 'pmc':
                        if link_db.get('Link') and link_db['Link'][0].get('Id'):
                            pmc_id = link_db['Link'][0]['Id']
                            break
            
            if not pmc_id:
                # Alternative method: check ArticleIdList for PMC IDs
                article_ids = article_data.get('PubmedData', {}).get('ArticleIdList', [])
                for article_id in article_ids:
                    if article_id.attributes.get('IdType') == 'pmc':
                        pmc_id = str(article_id).replace('PMC', '')
                        break
            
            if pmc_id:
                print(f"  ✓ Found in PMC! (PMC{pmc_id})")
                
                # Add to our list of available papers
                pmc_papers.append({
                    'PMID': pmid,
                    'PMC': pmc_id,
                    'DOI': doi,
                    'Title': title,
                    'Author': author_display,
                    'Year': pub_year
                })
            else:
                print(f"  ✗ Not available in PMC")
            
            # Pause briefly to avoid overloading NCBI's servers
            time.sleep(0.5)
            
        except Exception as e:
            print(f"  Error processing PMID {pmid}: {e}")
    
    # Create a DataFrame and save to CSV
    if pmc_papers:
        df = pd.DataFrame(pmc_papers)
        df.to_csv(ARTICLES_CSV, index=False)
        print(f"\nSaved {len(pmc_papers)} PMC-available papers to {ARTICLES_CSV}")
        return df
    else:
        print("No papers available in PMC were found.")
        return None

def download_from_pmc(pmc_id, author, year, title, pmid):
    """Download a paper from PubMed Central using its PMC ID."""
    try:
        # Create filename based on author and year
        filename = f"{clean_filename(author)}_{year}_PMC{pmc_id}.pdf"
        filepath = os.path.join(OUTPUT_FOLDER, filename)
        
        # PMC download URL format
        pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/pdf"
        
        headers = {
            'User-Agent': get_random_user_agent(),
            'Accept': 'text/html,application/xhtml+xml,application/xml,application/pdf',
            'Accept-Language': 'en-US,en;q=0.9'
        }
        
        print(f"Attempting to download: {author}_{year} - {title}")
        print(f"URL: {pmc_url}")
        
        # Make the request to PMC
        response = requests.get(pmc_url, headers=headers, timeout=30)
        
        # If we get a redirect, follow it
        if response.history:
            print(f"Followed {len(response.history)} redirects to {response.url}")
        
        # Check if we got a valid response
        if response.status_code != 200:
            print(f"Error: Received status code {response.status_code}")
            return False
        
        # Verify we got a PDF
        content_type = response.headers.get('Content-Type', '')
        if 'application/pdf' not in content_type and not response.url.endswith('.pdf'):
            print(f"Warning: Response doesn't appear to be a PDF. Content type: {content_type}")
            
            # Check if we got an HTML page that might contain a link to the PDF
            if 'text/html' in content_type:
                soup = BeautifulSoup(response.text, 'html.parser')
                
                # Look for PDF download link
                pdf_links = []
                for a in soup.find_all('a', href=True):
                    href = a['href']
                    if href.endswith('.pdf') or '/pdf/' in href:
                        pdf_links.append(href)
                
                if pdf_links:
                    # Try the first PDF link
                    pdf_url = pdf_links[0]
                    if not pdf_url.startswith('http'):
                        # Handle relative URLs
                        base_url = response.url.split('/pmc/')[0] if '/pmc/' in response.url else 'https://www.ncbi.nlm.nih.gov'
                        pdf_url = base_url + pdf_url if pdf_url.startswith('/') else base_url + '/' + pdf_url
                    
                    print(f"Found PDF link: {pdf_url}")
                    pdf_response = requests.get(pdf_url, headers=headers, timeout=30)
                    
                    if pdf_response.status_code == 200 and ('application/pdf' in pdf_response.headers.get('Content-Type', '') or pdf_url.endswith('.pdf')):
                        # Save the PDF
                        with open(filepath, 'wb') as f:
                            f.write(pdf_response.content)
                        
                        print(f"Successfully downloaded: {filename}")
                        return True
            
            # Try an alternative URL format
            alt_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/pdf/nihms-{pmid}.pdf"
            print(f"Trying alternative URL: {alt_url}")
            alt_response = requests.get(alt_url, headers=headers, timeout=30)
            
            if alt_response.status_code == 200 and 'application/pdf' in alt_response.headers.get('Content-Type', ''):
                # Save the PDF
                with open(filepath, 'wb') as f:
                    f.write(alt_response.content)
                
                print(f"Successfully downloaded from alternative URL: {filename}")
                return True
            
            return False
        
        # Save the PDF
        with open(filepath, 'wb') as f:
            f.write(response.content)
        
        print(f"Successfully downloaded: {filename}")
        return True
    
    except Exception as e:
        print(f"Error downloading from PMC: {e}")
        return False

def download_papers():
    """Download papers from the CSV file."""
    # Create output directory if it doesn't exist
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
    
    # Load CSV file
    try:
        df = pd.read_csv(ARTICLES_CSV)
        print(f"Loaded {len(df)} articles from {ARTICLES_CSV}")
    except Exception as e:
        print(f"Error loading CSV file: {e}")
        return
    
    # Track statistics
    successful = 0
    failed = 0
    
    # Process each article
    for idx, row in df.iterrows():
        try:
            pmc_id = row['PMC']
            author = row['Author']
            year = row['Year']
            title = row['Title']
            pmid = row['PMID']
            
            print(f"\nProcessing {idx+1}/{len(df)}: {author}_{year}")
            
            # Try to download from PMC
            if download_from_pmc(pmc_id, author, year, title, pmid):
                successful += 1
            else:
                failed += 1
                # Save failed downloads for later
                with open("failed_downloads.txt", "a") as f:
                    f.write(f"PMC{pmc_id},{author},{year},{pmid}\n")
            
            # Delay between requests
            wait_time = random.uniform(2, 4)  
            print(f"Waiting {wait_time:.1f} seconds before next download...")
            time.sleep(wait_time)
            
        except Exception as e:
            print(f"Unexpected error processing row {idx+1}: {e}")
            failed += 1
            time.sleep(2)
    
    # Print summary
    print("\n" + "="*50)
    print(f"Download Summary:")
    print(f"- Attempted: {len(df)} papers")
    print(f"- Successfully downloaded: {successful} papers")
    print(f"- Failed: {failed} papers")
    print(f"- Success rate: {successful/len(df)*100:.1f}%")
    if failed > 0:
        print(f"- Failed downloads saved to failed_downloads.txt")
    print("="*50)

if __name__ == "__main__":
    print("PubMed Central Paper Downloader")
    print("="*50)
    print("This script searches for and downloads open access papers from PubMed Central.")
    print("="*50 + "\n")
    
    # Step 1: Search PubMed and create CSV file
    papers_df = search_pubmed_central()
    
    if papers_df is not None and len(papers_df) > 0:
        # Step 2: Download papers
        download_papers()
    else:
        print("No papers available to download.")