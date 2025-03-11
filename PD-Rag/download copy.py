import requests
import os
import pandas as pd
import time
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from Bio import Entrez
import re

# Set up PubMed search parameters
Entrez.email = "your_email@example.com"  # Replace with your email

# Phase 1: Generate the articles.csv file with URLs and citation info
def generate_articles_csv():
    print("Phase 1: Generating articles.csv with Parkinson's disease paper URLs...")
    
    # Search for Parkinson's disease papers, sorted by date (most recent first)
    print("Searching PubMed for the 300 most recent papers on Parkinson's disease...")
    handle = Entrez.esearch(db="pubmed", term="parkinson's disease", retmax=300, sort="date")
    record = Entrez.read(handle)
    handle.close()
    
    # Get the list of PubMed IDs
    pmid_list = record['IdList']
    print(f"Found {len(pmid_list)} recent articles on Parkinson's disease.")
    
    # Fetch article details including URLs
    urls = []
    authors = []
    years = []
    titles = []
    
    for i, pmid in enumerate(pmid_list):
        print(f"Fetching details for article {i+1}/{len(pmid_list)}")
        try:
            # Get article full record to find DOI and citation info
            handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
            record = Entrez.read(handle)
            handle.close()
            
            # Set default values
            article_url = None
            author_string = "Unknown"
            pub_year = "0000"
            article_title = "Unknown_Title"
            
            if 'PubmedArticle' in record and 'MedlineCitation' in record['PubmedArticle'][0]:
                article = record['PubmedArticle'][0]['MedlineCitation']
                
                # Get article title
                if 'Article' in article and 'ArticleTitle' in article['Article']:
                    article_title = article['Article']['ArticleTitle']
                    titles.append(article_title)
                else:
                    titles.append("Unknown_Title")
                
                # Get publication year
                if 'Article' in article and 'Journal' in article['Article'] and 'JournalIssue' in article['Article']['Journal']:
                    journal_issue = article['Article']['Journal']['JournalIssue']
                    if 'PubDate' in journal_issue:
                        pub_date = journal_issue['PubDate']
                        if 'Year' in pub_date:
                            pub_year = pub_date['Year']
                            years.append(pub_year)
                        else:
                            years.append("0000")
                    else:
                        years.append("0000")
                else:
                    years.append("0000")
                
                # Get authors
                if 'Article' in article and 'AuthorList' in article['Article']:
                    author_list = article['Article']['AuthorList']
                    if len(author_list) > 0:
                        first_author = author_list[0]
                        if 'LastName' in first_author:
                            if len(author_list) > 1:
                                author_string = f"{first_author['LastName']}_et_al"
                            else:
                                author_string = first_author['LastName']
                            authors.append(author_string)
                        else:
                            authors.append("Unknown")
                    else:
                        authors.append("Unknown")
                else:
                    authors.append("Unknown")
                
                # Get DOI URL
                if 'Article' in article and 'ELocationID' in article['Article']:
                    for location in article['Article']['ELocationID']:
                        if location.attributes['EIdType'] == 'doi':
                            article_url = f"https://doi.org/{location}"
                            break
            
            # If we found a URL, add it to our list
            if article_url:
                urls.append(article_url)
                print(f"Added: {author_string}_{pub_year} - {article_title}")
            else:
                print(f"No URL found for PMID: {pmid}")
                if i < len(authors):  # Remove the author entry if we didn't add a URL
                    authors.pop()
                if i < len(years):
                    years.pop()
                if i < len(titles):
                    titles.pop()
            
            # Sleep to avoid overloading PubMed's servers
            time.sleep(0.5)
            
        except Exception as e:
            print(f"Error fetching data for PMID {pmid}: {e}")
    
    # Create a dataframe and save URLs to CSV
    if urls:
        # Ensure all lists are the same length
        min_length = min(len(urls), len(authors), len(years), len(titles))
        df = pd.DataFrame({
            'URL': urls[:min_length],
            'Author': authors[:min_length],
            'Year': years[:min_length],
            'Title': titles[:min_length]
        })
        df.to_csv("articles.csv", index=False)
        print(f"Saved {len(urls)} URLs to articles.csv")
        return True
    else:
        print("No URLs found to save.")
        return False

# Phase 2: Download PDFs from the CSV file
# Phase 2: Download PDFs from the CSV file
def download_from_csv():
    print("\nPhase 2: Downloading papers from articles.csv")
    
    # Folder to save the PDFs
    save_folder = "downloaded_articles"
    
    # Create the folder if it doesn't exist
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    
    # Load the CSV file
    try:
        df = pd.read_csv("articles.csv")
        print(f"Loaded {len(df)} articles from articles.csv")
    except Exception as e:
        print(f"Error loading CSV file: {e}")
        return
    
    # Define headers to mimic a real browser request
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
    }
    
    # Create a session with retry logic
    session = requests.Session()
    retries = Retry(total=3, backoff_factor=1, status_forcelist=[502, 503, 504])
    session.mount('https://', HTTPAdapter(max_retries=retries))
    
    # Count existing PDFs before we start
    initial_count = len([f for f in os.listdir(save_folder) if f.lower().endswith('.pdf')])
    
    def download_pdf(url, author, year, index):
        try:
            # Create citation-style filename
            file_name = f"{author}_{year}.pdf"
            
            # Clean filename of any invalid characters
            file_name = re.sub(r'[\\/*?:"<>|]', "", file_name)
            
            # Add index if there are multiple papers by the same author in the same year
            base_path = os.path.join(save_folder, file_name)
            file_path = base_path
            
            # If file exists with same name, add a letter suffix (a, b, c, etc.)
            if os.path.exists(file_path):
                name, ext = os.path.splitext(file_name)
                suffix = 'a'
                while os.path.exists(os.path.join(save_folder, f"{name}_{suffix}{ext}")):
                    # Move to next letter: a -> b -> c, etc.
                    suffix = chr(ord(suffix) + 1)
                file_name = f"{name}_{suffix}{ext}"
                file_path = os.path.join(save_folder, file_name)
            
            # Get the content from the URL with headers and timeout
            print(f"Downloading: {author}_{year} from {url}")
            response = session.get(url, headers=headers, timeout=20)
            response.raise_for_status()  # Check if the request was successful
            
            # Save the PDF file
            with open(file_path, "wb") as file:
                file.write(response.content)
            
            print(f"Downloaded: {file_name}")
            return True
        except requests.exceptions.Timeout:
            print(f"Timeout occurred while downloading {url}")
            return False
        except requests.exceptions.RequestException as e:
            print(f"Failed to download {url}: {e}")
            return False
    
    # Loop through all URLs and download PDFs with a 2-second delay
    successful_downloads = 0
    for idx, row in df.iterrows():
        if download_pdf(row['URL'], row['Author'], row['Year'], idx + 1):
            successful_downloads += 1
        time.sleep(2)
    
    # Count total PDFs after downloads
    final_count = len([f for f in os.listdir(save_folder) if f.lower().endswith('.pdf')])
    new_downloads = final_count - initial_count
    
    print(f"Download summary:")
    print(f"- Attempted to download: {len(df)} papers")
    print(f"- Successfully downloaded: {new_downloads} new PDFs")
    print(f"- Total PDFs in '{save_folder}': {final_count}")

# Main execution
if __name__ == "__main__":
    if generate_articles_csv():  # Generate the CSV first
        download_from_csv()      # Then download from it
    else:
        print("Failed to generate articles.csv, cannot proceed with downloads.")