"""WBGene identity and explicit protein associations, never accession nodes."""
import csv
import gzip
import re
from collections import Counter, defaultdict

WB = re.compile(r'WBGene\d{8}\Z')
UP = re.compile(r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?\Z')


def read_mapping(path, header=False):
    with open(path, newline='') as f:
        rows = list(csv.reader(f, delimiter='\t'))
    if header:
        if not rows or rows.pop(0) != ['WBGene', 'UniProt']:
            raise ValueError('Expected exact header WBGene<TAB>UniProt')
    if not rows:
        raise ValueError('Empty input mapping')
    valid, errors = [], []
    for i, row in enumerate(rows, start=2 if header else 1):
        if len(row) != 2 or not WB.fullmatch(row[0]) or not UP.fullmatch(row[1]):
            errors.append({'source_row': i, 'raw': row, 'reason': 'invalid_columns_or_identifier_syntax'})
        else:
            valid.append({'source_row': i, 'WBGene': row[0], 'UniProt': row[1]})
    return valid, errors


def mapping_qc(rows):
    pairs = Counter((r['WBGene'], r['UniProt']) for r in rows)
    genes, proteins = defaultdict(set), defaultdict(set)
    for g, p in pairs:
        genes[g].add(p)
        proteins[p].add(g)
    return {'valid_rows': len(rows), 'unique_WBGene': len(genes), 'unique_UniProt': len(proteins),
            'unique_pairs': len(pairs), 'duplicate_pair_extra_rows': sum(n-1 for n in pairs.values()),
            'one_gene_multiple_accessions': {k: sorted(v) for k,v in genes.items() if len(v)>1},
            'one_accession_multiple_genes': {k: sorted(v) for k,v in proteins.items() if len(v)>1},
            'orthogroups': None, 'orthogroup_status': 'not_provided'}


def read_gene_info(path, taxon):
    """A GeneID is resolvable only if all matching current xrefs agree."""
    mapping = defaultdict(set)
    with gzip.open(path, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not {'#tax_id', 'GeneID', 'dbXrefs'} <= set(reader.fieldnames or []):
            raise ValueError('NCBI gene_info schema missing required columns')
        for row in reader:
            if row['#tax_id'] == str(taxon):
                mapping[row['GeneID']].update(re.findall(r'WBGene\d{8}(?!\d)', row['dbXrefs']))
    return dict(mapping)


def resolve(gene_id, mapping):
    values = mapping.get(gene_id, set())
    if len(values) == 1:
        return next(iter(values)), 'resolved'
    return '', 'ambiguous' if values else 'unresolved'
