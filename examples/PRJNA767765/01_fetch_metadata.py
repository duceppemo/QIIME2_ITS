#!/usr/bin/env python3
"""Step 1: build the QIIME2 metadata file and the read-download manifest for
BioProject PRJNA767765, straight from the public archives.

Two public, unauthenticated web APIs, standard library only:
  * ENA's file report  -> one row per sequencing run: fastq URLs + MD5s
  * NCBI BioSample     -> each sample's name and attributes (country,
                          environment, collection date, elevation, replicate)

Writes, next to this script (or into the folder given as first argument):
  metadata.tsv           the QIIME2 sample metadata `qiime2-its -m` takes
  download_manifest.tsv  sample-id / run accession / fastq URLs / MD5s, read
                         by 02_download_reads.sh

Both files are also committed to the repository, so this step is optional --
it is here so nothing about the example is hand-made or unverifiable.
"""
import csv
import io
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

BIOPROJECT = 'PRJNA767765'

# The BioProject also holds a laboratory incubation experiment (microcosms of
# one soil spiked with BTEX or metals, plus its controls). This example keeps
# the 53 field samples only -- the natural-soil survey.
EXCLUDED_SAMPLE_PREFIXES = ('BTEX_soil', 'control_soil', 'metal_soil')

ENA_FILEREPORT = ('https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run'
                  '&fields=run_accession,sample_accession,library_layout,fastq_ftp,fastq_md5&format=tsv')
NCBI_EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

METADATA_COLUMNS = {  # column -> QIIME2 type
    'country': 'categorical',
    'env-medium': 'categorical',
    'env-local-scale': 'categorical',
    'collection-date': 'categorical',
    'elevation-m': 'numeric',
    'replicate': 'categorical',
}


def fetch(url, data=None):
    with urllib.request.urlopen(url, data=data, timeout=120) as response:
        return response.read().decode('utf-8')


def fetch_runs():
    text = fetch(ENA_FILEREPORT.format(accession=BIOPROJECT))
    return list(csv.DictReader(io.StringIO(text), delimiter='\t'))


def fetch_biosamples(accessions):
    """{BioSample accession: (sample name, {attribute: value})}"""
    body = urllib.parse.urlencode({'db': 'biosample', 'id': ','.join(accessions),
                                   'rettype': 'full', 'retmode': 'xml'}).encode()
    root = ET.fromstring(fetch(NCBI_EFETCH, data=body))
    biosamples = {}
    for biosample in root.findall('BioSample'):
        name = next((el.text for el in biosample.findall('.//Ids/Id') if el.get('db_label') == 'Sample name'), None)
        attributes = {(a.get('harmonized_name') or a.get('attribute_name')): a.text
                      for a in biosample.findall('.//Attributes/Attribute')}
        biosamples[biosample.get('accession')] = (name, attributes)
    return biosamples


def strip_envo_term(value):
    """ "contaminated soil [ENVO:00002116]" -> "contaminated soil" """
    return (value or '').split(' [ENVO')[0].strip()


def main():
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    out_dir.mkdir(parents=True, exist_ok=True)

    runs = fetch_runs()
    biosamples = fetch_biosamples(sorted({run['sample_accession'] for run in runs}))

    rows = []
    for run in runs:
        name, attributes = biosamples[run['sample_accession']]
        if name.startswith(EXCLUDED_SAMPLE_PREFIXES):
            continue
        if run['library_layout'] != 'PAIRED':
            raise SystemExit(f'{run["run_accession"]} is not paired-end -- the archive record changed?')
        geo = attributes.get('geo_loc_name') or ''
        rows.append({
            'sample-id': name,
            'run-accession': run['run_accession'],
            'fastq_ftp': run['fastq_ftp'],
            'fastq_md5': run['fastq_md5'],
            'country': geo.split(':')[0].strip(),
            'env-medium': strip_envo_term(attributes.get('env_medium')),
            'env-local-scale': strip_envo_term(attributes.get('env_local_scale')),
            'collection-date': attributes.get('collection_date') or '',
            'elevation-m': attributes.get('elev') or '',
            'replicate': attributes.get('replicate') or '',
        })
    rows.sort(key=lambda row: row['sample-id'])

    sample_ids = [row['sample-id'] for row in rows]
    if len(set(sample_ids)) != len(sample_ids):
        raise SystemExit('A sample has more than one run -- the archive record changed?')

    with open(out_dir / 'metadata.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['sample-id', *METADATA_COLUMNS])
        writer.writerow(['#q2:types', *METADATA_COLUMNS.values()])
        for row in rows:
            writer.writerow([row['sample-id'], *(row[column] for column in METADATA_COLUMNS)])

    with open(out_dir / 'download_manifest.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['sample-id', 'run-accession', 'fastq_ftp', 'fastq_md5'])
        for row in rows:
            writer.writerow([row[key] for key in ('sample-id', 'run-accession', 'fastq_ftp', 'fastq_md5')])

    print(f'{len(rows)} field samples -> {out_dir / "metadata.tsv"}, {out_dir / "download_manifest.tsv"}')


if __name__ == '__main__':
    main()
