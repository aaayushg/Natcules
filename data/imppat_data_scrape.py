import csv
import time

import requests
from bs4 import BeautifulSoup

base_url = "https://cb.imsc.res.in/imppat/phytochemical-detailedpage/IMPHY{:06d}"
output_file = "imppat_smiles.csv"
headers = {
    "User-Agent": "Mozilla/5.0"
}

def extract_smiles(soup):
    for strong in soup.find_all("strong"):
        if "SMILES:" in strong.text:
            text_tag = strong.find_next("text")
            if text_tag:
                return text_tag.text.strip()
    return None

with open(output_file, mode='w', newline='', encoding='utf-8') as file:
    writer = csv.writer(file)
    writer.writerow(["IMPHY_ID", "SMILES"])

    for i in range(18001):  # 0 to 18000 inclusive
        url = base_url.format(i)
        try:
            response = requests.get(url, headers=headers, timeout=10)
            if response.status_code == 200:
                soup = BeautifulSoup(response.text, "html.parser")
                smiles = extract_smiles(soup)
                if smiles:
                    writer.writerow([f"IMPHY{i:06d}", smiles])
                    print(f"[{i}] Found: {smiles}")
                else:
                    print(f"[{i}] SMILES not found")
            else:
                print(f"[{i}] Page not found (Status {response.status_code})")
        except Exception as e:
            print(f"[{i}] Error: {e}")
        time.sleep(0.3)  # Throttle requests to avoid server overload

