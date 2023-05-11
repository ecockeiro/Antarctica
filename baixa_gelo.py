import urllib.request
from bs4 import BeautifulSoup

url = 'https://thredds.met.no/thredds/catalog/osisaf/met.no/ice/conc/catalog.html'

# Faz o download do catálogo do OSISAF
with urllib.request.urlopen(url) as response:
    html = response.read()

# Analisa o HTML do catálogo usando a biblioteca BeautifulSoup
soup = BeautifulSoup(html, 'html.parser')

# Encontra todos os links para arquivos NetCDF no catálogo
links = soup.find_all('a', href=True)
for link in links:
    # Verifica se o link aponta para um arquivo NetCDF e se o ano está entre 2010 e 2019
    if link['href'].endswith('.nc') and ('2010' in link['href'] or '2011' in link['href'] or '2012' in link['href'] or '2013' in link['href'] or '2014' in link['href'] or '2015' in link['href'] or '2016' in link['href'] or '2017' in link['href'] or '2018' in link['href'] or '2019' in link['href']):
        # Faz o download do arquivo NetCDF
        file_url = url[:-12] + link['href']
        urllib.request.urlretrieve(file_url, link['href'])
        print('Arquivo ' + link['href'] + ' baixado com sucesso!')