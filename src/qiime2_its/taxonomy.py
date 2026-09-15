"""NCBI taxonomy parsing and expansion.

Shared by the ncbi- and fasta-based classifier trainers (previously duplicated
almost verbatim between the two scripts).
"""
import gzip

UNCLASSIFIED_TAXID = '12908'  # NCBI taxid for "unclassified sequences"

_RANK_TO_CODE = {
    'species': 's',
    'genus': 'g',
    'family': 'f',
    'order': 'o',
    'clade': 'c',
    'phylum': 'p',
    'superkingdom': 'k',
}


def _opener(path):
    return gzip.open if str(path).endswith('.gz') else open


def parse_id_table(id_table_path):
    """Parse a 2-column (accession, taxid) TSV into {taxid: accession}."""
    id_dict = {}
    with _opener(id_table_path)(id_table_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) != 2:
                raise ValueError('"id_table" must have exactly two tab-separated columns (accession, taxid)')
            acc, taxid = fields
            id_dict[taxid] = acc
    return id_dict


def extract_accessions_from_fasta(input_fasta, acc_file):
    """Write one accession per line (from fasta headers) to `acc_file`; return {accession: ''}."""
    acc_dict = {}
    with open(input_fasta) as in_fh, open(acc_file, 'w') as out_fh:
        for line in in_fh:
            line = line.rstrip()
            if not line.startswith('>'):
                continue
            acc = line[1:].split()[0]
            out_fh.write(f'{acc}\n')
            acc_dict[acc] = ''
    return acc_dict


def parse_accession2taxid(acc2taxid_path, accessions):
    """Parse an NCBI accession2taxid[.gz] file, keeping only entries in `accessions`.

    `accessions` may include versioned (e.g. "AB123456.1") accessions; matching
    is done on the unversioned accession. Returns {unversioned_accession: taxid}.
    """
    wanted = {acc.split('.')[0] for acc in accessions}
    result = {}
    with _opener(acc2taxid_path)(acc2taxid_path, 'rt') as f:
        next(f, None)  # header
        for line in f:
            fields = line.split('\t')
            acc, taxid = fields[0], fields[2]
            if acc in wanted:
                result[acc] = taxid
    return result


def accessions_to_taxids(acc_dict, acc2taxid_dict, output_file, unclassified_taxid=UNCLASSIFIED_TAXID):
    """Write "accession\\ttaxid" for every key of `acc_dict`, falling back to
    `unclassified_taxid` for accessions missing from `acc2taxid_dict`.

    Returns the list of accessions that had no taxid match.
    """
    missing = []
    with open(output_file, 'w') as f:
        for acc in acc_dict:
            taxid = acc2taxid_dict.get(acc.split('.')[0])
            if taxid is None:
                missing.append(acc)
                taxid = unclassified_taxid
            f.write(f'{acc}\t{taxid}\n')
    return missing


def parse_nodes_dmp(nodes_file):
    """Parse NCBI taxdump nodes.dmp into {taxid: (parent_taxid, rank)}."""
    node_dict = {}
    with open(nodes_file) as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            taxid, parent, rank = fields[0], fields[2], fields[4]
            node_dict[taxid] = (parent, rank)
    return node_dict


def parse_names_dmp(names_file):
    """Parse NCBI taxdump names.dmp into {taxid: scientific_name}."""
    names_dict = {}
    with open(names_file) as f:
        for line in f:
            if 'scientific name' not in line:
                continue
            fields = line.split('\t')
            taxid, name = fields[0], fields[2]
            names_dict.setdefault(taxid, name)
    return names_dict


def apply_merged_taxids(id_dict, merged_file):
    """Remap taxids in `id_dict` per NCBI taxdump merged.dmp. Returns a new dict."""
    id_dict = dict(id_dict)
    with open(merged_file) as f:
        for line in f:
            fields = line.split('\t')
            old_taxid, new_taxid = fields[0], fields[2]
            if old_taxid in id_dict:
                id_dict[new_taxid] = id_dict.pop(old_taxid)
    return id_dict


def lineage_string(taxid, node_dict, names_dict):
    """Build a "k__x;p__x;...;s__x" QIIME2-style lineage string for `taxid`."""
    taxo = {code: 'unidentified' for code in 'kpcofgs'}
    current = taxid
    while current != '1' and current in node_dict:
        parent, rank = node_dict[current]
        code = _RANK_TO_CODE.get(rank)
        if code:
            name = names_dict.get(current, 'unidentified')
            taxo[code] = name.replace(' ', '_') if code == 's' else name
        current = parent
    return 'k__{k};p__{p};c__{c};o__{o};f__{f};g__{g};s__{s}'.format(**taxo)


def write_taxonomy_file(id_dict, taxonomy_file, nodes_file, names_file, merged_file):
    """Write a QIIME2 HeaderlessTSVTaxonomyFormat file from {taxid: accession} + taxdump files."""
    node_dict = parse_nodes_dmp(nodes_file)
    names_dict = parse_names_dmp(names_file)
    id_dict = apply_merged_taxids(id_dict, merged_file)

    with open(taxonomy_file, 'w') as f:
        for taxid, acc in id_dict.items():
            f.write(f'{acc}\t{lineage_string(taxid, node_dict, names_dict)}\n')
