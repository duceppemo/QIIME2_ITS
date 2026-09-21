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
    'class': 'c',
    'phylum': 'p',
    'kingdom': 'k',
}


def _opener(path):
    return gzip.open if str(path).endswith('.gz') else open


def parse_id_table(id_table_path):
    """Parse a 2-column (accession, taxid) TSV into {accession: taxid}.

    Keyed by accession, not taxid: many reference sequences share one taxid
    (every ITS record of the same species), so a {taxid: accession} mapping
    silently keeps only the last accession per taxon -- confirmed against the
    bundled validation data, where 20 reference sequences produced only 4
    taxonomy lines, leaving the classifier trained on one sequence per taxon.
    """
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
            id_dict[acc] = taxid
    return id_dict


def read_fasta_ids(input_fasta):
    """Sequence IDs (header text up to the first whitespace) of a fasta[.gz], in file order."""
    with _opener(input_fasta)(input_fasta, 'rt') as f:
        return [line[1:].split()[0] for line in f if line.startswith('>') and line[1:].strip()]


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
    """Remap the taxids in an {accession: taxid} dict per NCBI taxdump
    merged.dmp (old taxid -> the taxid it was merged into). Returns a new dict."""
    wanted = set(id_dict.values())
    merged = {}
    with open(merged_file) as f:
        for line in f:
            fields = line.split('\t')
            if len(fields) < 3:
                continue
            old_taxid, new_taxid = fields[0], fields[2]
            if old_taxid in wanted:
                merged[old_taxid] = new_taxid
    return {acc: merged.get(taxid, taxid) for acc, taxid in id_dict.items()}


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
    """Write a QIIME2 HeaderlessTSVTaxonomyFormat file (one line per
    accession) from an {accession: taxid} dict + taxdump files."""
    node_dict = parse_nodes_dmp(nodes_file)
    names_dict = parse_names_dmp(names_file)
    id_dict = apply_merged_taxids(id_dict, merged_file)

    lineages = {}  # many accessions share a taxid -- walk each lineage once
    with open(taxonomy_file, 'w') as f:
        for acc, taxid in id_dict.items():
            if taxid not in lineages:
                lineages[taxid] = lineage_string(taxid, node_dict, names_dict)
            f.write(f'{acc}\t{lineages[taxid]}\n')


def max_lineage_depth(taxonomy_tsv_path):
    """Greatest number of ';'-separated rank fields across every lineage
    string in a QIIME2 taxonomy.tsv export (header + "id\\tlineage\\t..."
    rows, as rewritten by biom_utils.rewrite_taxonomy_header).

    `qiime taxa collapse --p-level N` fails outright if N exceeds this for
    every feature -- a classifier's assignments commonly don't reach genus/
    species for reads unrelated to its training set, so a fixed level should
    be capped to what's actually present rather than assumed.
    """
    max_depth = 0
    with open(taxonomy_tsv_path) as f:
        next(f, None)  # header
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 2:
                continue
            max_depth = max(max_depth, len(fields[1].split(';')))
    return max_depth
