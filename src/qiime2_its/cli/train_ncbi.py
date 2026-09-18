"""Build a QIIME2 classifier from an NCBI nucleotide query (or accession list) plus taxonomy."""
import argparse
import http.client
import os
from pathlib import Path
from time import sleep, time

from Bio import Entrez

from qiime2_its import downloader, env_checks, qiime_wrapper, taxonomy, timing
from qiime2_its._version import __version__

TAXDUMP_URL = 'https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'
ACC2TAXID_URL = 'https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
DEAD_ACC2TAXID_URL = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz'

BATCH_SIZE = 300


def _efetch_fasta(**kwargs):
    """efetch with one retry on a transient network error."""
    try:
        return Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', **kwargs)
    except (http.client.IncompleteRead, ValueError) as e:
        print(f'Network error ({e}). Retrying once...')
        return Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', **kwargs)


def download_sequences(query, seq_file, email, api_key):
    """`query` is either an NCBI search term or a path to a file of one accession per line."""
    Entrez.email = email
    Entrez.api_key = api_key
    Entrez.sleep_between_tries = 15
    Entrez.max_tries = 3

    with open(seq_file, 'w') as out_handle:
        if os.path.isfile(query):
            acc_list = [line.strip() for line in open(query) if line.strip()]
            count = len(acc_list)
            print(f'Downloading {count} records...')
            for start in range(0, count, BATCH_SIZE):
                end = min(count, start + BATCH_SIZE)
                print(f'\tDownloading record {start + 1} to {end} ({count})')
                fetch_handle = _efetch_fasta(id=','.join(acc_list[start:end]))
                out_handle.write(fetch_handle.read())
                fetch_handle.close()
                sleep(0.5)
        else:
            print(f"Searching NCBI for '{query}'.")
            search_handle = Entrez.esearch(db='nucleotide', term=query, idtype='acc', usehistory='y')
            search_results = Entrez.read(search_handle)
            search_handle.close()

            count = int(search_results['Count'])
            webenv = search_results['WebEnv']
            query_key = search_results['QueryKey']

            print(f'Downloading {count} records...')
            for start in range(0, count, BATCH_SIZE):
                end = min(count, start + BATCH_SIZE)
                print(f'\tDownloading record {start + 1} to {end} ({count})')
                fetch_handle = _efetch_fasta(retstart=start, retmax=BATCH_SIZE, webenv=webenv,
                                              query_key=query_key, idtype='acc')
                out_handle.write(fetch_handle.read())
                fetch_handle.close()
                sleep(0.5)


def run(query, output_folder, threads, email, api_key, taxdump, acc2taxid, dead_acc2taxid):
    env_checks.check_qiime2_env_active()

    output_folder = Path(output_folder)
    if not query:
        raise ValueError('Your query is empty.')
    output_folder.mkdir(parents=True, exist_ok=True)
    threads = env_checks.clamp_cpu(threads)

    t_zero = time()
    seq_file = output_folder / 'seq.fasta'
    if seq_file.exists():
        # Skipping a re-download is a deliberate resume optimization (an
        # NCBI query can be slow/rate-limited), but silently reusing this
        # file is dangerous: it doesn't know whether it's actually a
        # complete download for *this* query, a partial one left behind by
        # a crashed prior run, or a stale one from a different query that
        # happened to reuse this output folder -- so make that risk explicit
        # rather than silently training on whatever's there.
        print(f'{seq_file} already exists, reusing it instead of re-downloading -- delete it '
              f'(or use a different -o) to force a fresh download for this query.')
    else:
        download_sequences(query, seq_file, email, api_key)

    start = time()
    if taxdump:
        print('Extracting taxdump.tar.gz...', end='', flush=True)
        downloader.extract_targz(taxdump, output_folder)
    else:
        print('Downloading taxdump.tar.gz...', end='', flush=True)
        downloader.download(TAXDUMP_URL, output_folder / 'taxdump.tar.gz')
        downloader.extract_targz(output_folder / 'taxdump.tar.gz', output_folder)
    print(f' took {timing.format_elapsed(time() - start)}')

    if not acc2taxid:
        start = time()
        print('Downloading nucl_gb.accession2taxid.gz...', end='', flush=True)
        downloader.download(ACC2TAXID_URL, output_folder / 'nucl_gb.accession2taxid.gz')
        acc2taxid = output_folder / 'nucl_gb.accession2taxid.gz'
        print(f' took {timing.format_elapsed(time() - start)}')

    if not dead_acc2taxid:
        start = time()
        print('Downloading dead_nucl.accession2taxid.gz...', end='', flush=True)
        downloader.download(DEAD_ACC2TAXID_URL, output_folder / 'dead_nucl.accession2taxid.gz')
        dead_acc2taxid = output_folder / 'dead_nucl.accession2taxid.gz'
        print(f' took {timing.format_elapsed(time() - start)}')

    start = time()
    print('Extracting accession numbers from fasta file...', end='', flush=True)
    acc_file = output_folder / 'acc.list'
    acc_dict = taxonomy.extract_accessions_from_fasta(seq_file, acc_file)
    print(f' took {timing.format_elapsed(time() - start)}')

    start = time()
    print('Parsing nucl_gb.accession2taxid.gz...', end='', flush=True)
    acc2taxid_dict = taxonomy.parse_accession2taxid(acc2taxid, acc_dict)
    print(f' took {timing.format_elapsed(time() - start)}')

    start = time()
    print('Parsing dead_nucl.accession2taxid.gz...', end='', flush=True)
    acc2taxid_dict.update(taxonomy.parse_accession2taxid(dead_acc2taxid, acc_dict))
    print(f' took {timing.format_elapsed(time() - start)}')

    start = time()
    print('Finding taxID for accessions...', end='', flush=True)
    taxid_file = output_folder / 'taxid.list'
    missing = taxonomy.accessions_to_taxids(acc_dict, acc2taxid_dict, taxid_file)
    print(f' took {timing.format_elapsed(time() - start)}')
    if missing:
        print(f'\nThe following accessions were not found in the taxdump files: {", ".join(missing)}')

    start = time()
    print('Writing taxonomy...', end='', flush=True)
    taxonomy_file = output_folder / 'taxonomy.txt'
    id_dict = taxonomy.parse_id_table(taxid_file)  # {taxid: accession}
    taxonomy.write_taxonomy_file(id_dict, taxonomy_file, output_folder / 'nodes.dmp',
                                  output_folder / 'names.dmp', output_folder / 'merged.dmp')
    print(f' took {timing.format_elapsed(time() - start)}')

    qiime2_seq = seq_file.with_suffix(seq_file.suffix + '.qza')
    qiime2_taxo = taxonomy_file.with_suffix(taxonomy_file.suffix + '.qza')
    classifier_file = output_folder / 'seq_ncbi.qza'

    print('Importing sequences into QIIME2...')
    qiime_wrapper.import_sequences(seq_file, qiime2_seq)
    print('Importing taxonomy into QIIME2...')
    qiime_wrapper.import_taxonomy(taxonomy_file, qiime2_taxo)
    print('Training classifier...')
    qiime_wrapper.train_naive_bayes_classifier(qiime2_seq, qiime2_taxo, classifier_file)

    print(f'Done (total time: {timing.format_elapsed(time() - t_zero)}).')


def build_parser():
    parser = argparse.ArgumentParser(description="Download DNA sequences from NCBI and add taxonomy for QIIME2.")
    parser.add_argument('-q', '--query',
                         metavar='"txid4762[Organism:exp] AND (\\"internal transcribed spacer\\"[Title])"',
                         required=True, type=str,
                         help='NCBI query string, OR a text file with one accession number per line. Mandatory.')
    parser.add_argument('-o', '--output', metavar='/output_folder/', required=True, type=str,
                         help='Output folder. Mandatory.')
    parser.add_argument('-t', '--threads', metavar='4', default=4, type=int,
                         help='Number of CPUs. Default is 4. Optional.')
    parser.add_argument('-e', '--email', metavar='your.email@example.org', default='your.email@example.org',
                         type=str, help='Your email address. Optional.')
    parser.add_argument('-a', '--api-key', metavar='API_KEY', default=None, type=str,
                         help='Your NCBI API key. Allows up to 10 requests per second instead of 3. Optional.')
    parser.add_argument('--taxdump', metavar='/path/to/taxdump.tar.gz', default=None, type=str,
                         help='Path to a downloaded taxdump.tar.gz. Downloaded otherwise. Optional.')
    parser.add_argument('--acc2taxid', metavar='/path/to/nucl_gb.accession2taxid.gz', default=None, type=str,
                         help='Path to a downloaded nucl_gb.accession2taxid.gz. Downloaded otherwise. Optional.')
    parser.add_argument('--dead-acc2taxid', metavar='/path/to/dead_nucl.accession2taxid.gz', default=None,
                         type=str,
                         help='Path to a downloaded dead_nucl.accession2taxid.gz. Downloaded otherwise. Optional.')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    return parser


def main():
    args = build_parser().parse_args()
    run(args.query, args.output, args.threads, args.email, args.api_key,
        args.taxdump, args.acc2taxid, args.dead_acc2taxid)


if __name__ == '__main__':
    main()
