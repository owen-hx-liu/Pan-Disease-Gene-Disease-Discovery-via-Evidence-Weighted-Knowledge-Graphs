#!/usr/bin/env python3
"""Biomedical Identifier Auto-Detection"""

import re
import json
import gzip
import csv
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict

COLUMN_HINTS = {
    'gene': ['gene', 'hgnc', 'symbol', 'ensembl', 'gene_id', 'entrez', 'geneid', 'ncbi_gene', 'refseq', 'mgi', 'rgd', 'transcript', 'ensg', 'ensembl_gene', 'ensembl_transcript'],
    'disease': ['disease', 'disorder', 'condition', 'mondo', 'doid', 'omim', 'mesh', 'icd', 'icd10', 'icd11', 'umls', 'orphanet', 'ordo', 'snomed', 'snomedct', 'orpha', 'orpha_id', 'orpha_code', 'orphan_id', 'orphan_code', 'orphanet_id', 'orphadata', 'orphanumber', 'orphacode'],
    'phenotype': ['phenotype', 'hpo', 'hp', 'human_phenotype', 'mp', 'mammalian_phenotype', 'trait', 'zp', 'zebrafish_phenotype', 'pat', 'phenotypic_abnormality'],
    'therapy': ['therapy', 'treatment', 'medical_action', 'maxo', 'procedure', 'intervention', 'therapeutic_procedure', 'preventative_therapy', 'palliative_care', 'complementary_therapy', 'alternative_therapy', 'medical_procedure'],
    'protein': ['protein', 'uniprot', 'uniprotkb', 'ensembl', 'refseq', 'pdb', 'interpro', 'pfam', 'accession', 'uniprot_ac', 'uniprot_id', 'ensembl_protein'],
    'compound': ['drug', 'compound', 'chemical', 'chembl', 'drugbank', 'pubchem', 'pubchem_cid', 'kegg', 'kegg_compound', 'chebi', 'inchi', 'inchikey', 'smiles', 'rxcui'],
    'variant': ['variant', 'snp', 'rsid', 'dbsnp', 'clinvar', 'cosmic', 'mutation', 'hgvs', 'clinvar_id', 'cosmic_id', 'gnomad', 'af', 'allele_frequency'],
    'pathway': ['pathway', 'kegg', 'reactome', 'go', 'gene_ontology', 'biological_process', 'wiki_pathway', 'wp', 'pathway_id'],
    'anatomy': ['anatomy', 'tissue', 'organ', 'uberon', 'cell', 'cl', 'cell_type', 'fma', 'cell_ontology', 'efo', 'uberon_id'],
    'taxonomy': ['taxonomy', 'taxon', 'species', 'organism', 'ncbi', 'ncbitaxon', 'taxid', 'taxonomy_id']
}

IDENTIFIER_PATTERNS = {
    'gene': {
        'HGNC': r'^HGNC:\d+$',
        'HGNC_Symbol': r'^[A-Z][A-Z0-9\-]{1,9}$',
        'Ensembl_Gene': r'^ENSG\d{11}$',
        'Ensembl_Transcript': r'^ENST\d{11}$',
        'NCBIGene': r'^(NCBIGene:|ENTREZGENE:)?\d{3,9}$',
        'RefSeq_Gene': r'^(NM_|NR_|XM_|XR_)\d+(\.\d+)?$',
        'MGI': r'^MGI:\d+$',
        'RGD': r'^RGD:\d+$',
        'ZFIN': r'^ZFIN:ZDB-GENE-\d+-\d+$',
        'FlyBase': r'^FBgn\d{7}$',
        'VGNC': r'^VGNC:\d+$'
    },
    'disease': {
        'MONDO': r'^MONDO:\d{7}$',
        'DOID': r'^DOID:\d+$',
        'OMIM': r'^(OMIM:)?\d{6,7}$',
        'Orphanet': r'^(ORPHA:)?\d{1,7}$',
        'SNOMED': r'^(SNOMEDCT:)?\d{6,18}$',
        'UMLS': r'^(UMLS:)?C\d{7}$',
        'MeSH': r'^(MESH:)?[CD]\d{6,9}$',
        'EFO': r'^EFO:\d{7}$'
    },
    'phenotype': {
        'HPO': r'^HP:\d{7}$',
        'MP': r'^MP:\d{7}$',
        'ZP': r'^ZP:\d{7}$',
        'UBERON': r'^UBERON:\d{7}$'
    },
    'therapy': {
        'MAXO': r'^MAXO:\d{7}$'
    },
    'variant': {
        'dbSNP': r'^rs\d+$',
        'ClinVar_RCV': r'^(ClinVar:)?RCV\d+(\.\d+)?$',
        'ClinVar_VCV': r'^(ClinVar:)?VCV\d+(\.\d+)?$',
        'COSMIC': r'^(COSMIC:)?COSM\d+$',
        'HGVS_genomic': r'^g\.\d+[A-Z>del]+',
        'HGVS_coding': r'^c\.\d+[A-Z>del]+',
        'HGVS_protein': r'^p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}',
        'gnomAD': r'^\d+-\d+-[ACGT]-[ACGT]$'
    },
    'protein': {
        'UniProtKB': r'^(UniProtKB:)?([A-NR-Z][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z]\d[A-Z]\d[A-Z]\d)$',
        'PDB': r'^(PDB:)?[0-9][A-Z0-9]{3}$',
        'Ensembl_Protein': r'^ENSP\d{11}$',
        'RefSeq_Protein': r'^(NP_|XP_|YP_|WP_)\d+(\.\d+)?$',
        'InterPro': r'^IPR\d{6}$',
        'Pfam': r'^PF\d{5}$'
    },
    'compound': {
        'PubChem_CID': r'^(CID:)?\d{4,9}$',
        'ChEBI': r'^CHEBI:\d+$',
        'DrugBank': r'^DB\d{5}$',
        'ChEMBL': r'^CHEMBL\d+$',
        'KEGG_Compound': r'^C\d{5}$',
        'CAS': r'^(CAS:)?\d{2,7}-\d{2}-\d$',
        'InChI': r'^InChI=.+',
        'InChIKey': r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$',
        'SMILES': r'^[CNOPSFIBrClcnops0-9\[\]\(\)=#@\+\-\\\/\.]{10,}$',
        'RxCUI': r'^\d{6,7}$'
    },
    'pathway': {
        'Reactome': r'^R-[A-Z]{3}-\d+(-\d+)?$',
        'KEGG_Pathway': r'^(KEGG:)?(hsa|mmu|rno|dme)\d{5}$',
        'WikiPathways': r'^WP:\d+$',
        'GO': r'^GO:\d{7}$'
    },
    'anatomy': {
        'UBERON': r'^UBERON:\d{7}$',
        'FMA': r'^FMA:\d+$',
        'CL': r'^CL:\d{7}$'
    },
    'taxonomy': {
        'NCBITaxon': r'^(NCBITaxon:)?\d{4,7}$'
    }
}

# Compile all regex patterns
for category in IDENTIFIER_PATTERNS:
    IDENTIFIER_PATTERNS[category] = {
        name: re.compile(pattern, re.IGNORECASE) 
        for name, pattern in IDENTIFIER_PATTERNS[category].items()
    }


def detect_identifiers(raw_data_path='data/raw', output_csv='identifier_detection_results.csv'):
    """Main function to detect biomedical identifiers in datasets"""
    raw_path = Path(raw_data_path)
    
    if not raw_path.exists():
        print(f"Error: Directory '{raw_data_path}' does not exist", file=sys.stderr)
        print(f"Please create the directory or specify a valid path.", file=sys.stderr)
        return 1
    
    datasets = [d for d in raw_path.iterdir() if d.is_dir() and not d.name.startswith('.')]
    
    if not datasets:
        print(f"Warning: No dataset directories found in '{raw_data_path}'", file=sys.stderr)
        print(f"Expected structure: {raw_data_path}/dataset_name/data_files", file=sys.stderr)
        return 1
    
    all_results = []
    
    print("=" * 80)
    print("BIOMEDICAL IDENTIFIER DETECTION")
    print("=" * 80 + "\n")
    
    for dataset_path in sorted(datasets):
        dataset_id = dataset_path.name
        print(f"📁 {dataset_id}")
        print("-" * 80)
        
        all_detections = defaultdict(set)
        column_info = defaultdict(lambda: defaultdict(set))
        file_count = 0
        
        for file_path in dataset_path.rglob('*'):
            if not file_path.is_file():
                continue
            
            ext = file_path.suffix.lower()
            if ext == '.gz':
                ext = '.' + file_path.stem.split('.')[-1] if '.' in file_path.stem else ''
            
            if ext not in ['.csv', '.tsv', '.json', '.jsonl', '.ndjson', '.xml']:
                continue
            
            file_count += 1
            try:
                detections, col_mapping = _scan_file(file_path, ext)
                for cat, ids in detections.items():
                    all_detections[cat].update(ids)
                for col, col_dets in col_mapping.items():
                    for cat, ids in col_dets.items():
                        column_info[col][cat].update(ids)
            except Exception as e:
                print(f"  ⚠️  {file_path.name}: {e}")
        
        print(f"  Files scanned: {file_count}\n")
        
        # Print detected identifiers
        for cat in COLUMN_HINTS.keys():
            if cat in all_detections and all_detections[cat]:
                ids = ', '.join(sorted(all_detections[cat]))
                print(f"  🔍 {cat.upper():<12} {ids}")
                all_results.append({
                    'dataset': dataset_id,
                    'category': cat,
                    'identifier_types': ids
                })
        
        # Print column information
        if column_info:
            print("\n  Fields/Columns:")
            for col, cats in sorted(column_info.items())[:10]:
                col_ids = [id_type for cat in cats.values() for id_type in cat]
                print(f"    • {col:<30} → {', '.join(col_ids)}")
            if len(column_info) > 10:
                print(f"    ... +{len(column_info) - 10} more")
        
        if not all_detections:
            print("  ⚠️  No identifiers detected")
        
        print()
    
    # Write results to CSV
    if all_results:
        try:
            with open(output_csv, 'w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=['dataset', 'category', 'identifier_types'])
                writer.writeheader()
                writer.writerows(all_results)
            print(f"✅ Results saved to: {output_csv}")
            print(f"   Total records: {len(all_results)}")
        except Exception as e:
            print(f"❌ Error writing CSV: {e}", file=sys.stderr)
            return 1
    else:
        print("⚠️  No identifiers detected in any dataset. CSV not created.")
    
    print("\n" + "=" * 80)
    print("✓ Detection complete!")
    print("=" * 80)
    
    return 0


def _scan_file(file_path, ext):
    """Scan a file and return detected identifiers"""
    is_gz = file_path.suffix.lower() == '.gz'
    
    if ext in ['.csv', '.tsv']:
        return _scan_tabular(file_path, is_gz, '\t' if ext == '.tsv' else ',')
    elif ext in ['.json', '.jsonl', '.ndjson']:
        return _scan_json(file_path, is_gz, ext in ['.jsonl', '.ndjson'])
    elif ext == '.xml':
        return _scan_xml(file_path, is_gz)
    
    return {}, {}


def _scan_xml(file_path, is_gz):
    """Scan XML data files"""
    opener = gzip.open if is_gz else open
    detections = defaultdict(set)
    column_mapping = defaultdict(lambda: defaultdict(set))
    
    try:
        with opener(file_path, 'rb') as f:
            tree = ET.parse(f)
            root = tree.getroot()
        
        # Collect all text values from XML elements and attributes
        element_data = defaultdict(list)
        
        def extract_text(element, path=''):
            """Recursively extract text from XML elements"""
            # Get element tag (without namespace)
            tag = element.tag.split('}')[-1] if '}' in element.tag else element.tag
            current_path = f"{path}/{tag}" if path else tag
            
            # Extract text content
            if element.text and element.text.strip():
                element_data[tag].append(element.text.strip())
            
            # Extract attributes
            for attr_name, attr_value in element.attrib.items():
                if attr_value and attr_value.strip():
                    element_data[f"{tag}@{attr_name}"].append(attr_value.strip())
            
            # Recursively process children
            for child in element:
                extract_text(child, current_path)
        
        extract_text(root)
        
        # Scan collected data
        for element_name, values in element_data.items():
            if not values:
                continue
            
            # Sample values for detection
            sample_values = values[:100]
            
            # Determine which categories to check based on element name
            element_lower = element_name.lower()
            matching_categories = []
            
            for category, hints in COLUMN_HINTS.items():
                if any(hint in element_lower for hint in hints):
                    matching_categories.insert(0, category)
                else:
                    matching_categories.append(category)
            
            # Detect identifiers
            for cat in matching_categories:
                detected = _detect_type(sample_values, cat)
                if detected:
                    detections[cat].update(detected)
                    column_mapping[element_name][cat].update(detected)
        
    except ET.ParseError as e:
        # If XML parsing fails, try to extract text values directly
        try:
            with opener(file_path, 'rt', encoding='utf-8', errors='ignore') as f:
                content = f.read(100000)  # Read first 100KB
                
                # Extract values from XML tags using regex
                tag_pattern = r'<([^/>]+)>([^<]+)</\1>'
                matches = re.findall(tag_pattern, content)
                
                for tag, value in matches[:200]:
                    if value and value.strip():
                        tag_lower = tag.lower()
                        
                        # Match to categories
                        for category, hints in COLUMN_HINTS.items():
                            if any(hint in tag_lower for hint in hints):
                                detected = _detect_type([value.strip()], category)
                                if detected:
                                    detections[category].update(detected)
                                    column_mapping[tag][category].update(detected)
        except:
            pass
    
    return detections, column_mapping


def _scan_tabular(file_path, is_gz, delimiter):
    """Scan tabular data files (CSV/TSV)"""
    opener = gzip.open if is_gz else open
    detections = defaultdict(set)
    column_mapping = defaultdict(lambda: defaultdict(set))
    
    with opener(file_path, 'rt', encoding='utf-8', errors='ignore') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        
        try:
            rows = [next(reader)]
            for _ in range(199):
                try:
                    rows.append(next(reader))
                except StopIteration:
                    break
        except:
            return detections, column_mapping
        
        if not rows:
            return detections, column_mapping
        
        fieldnames = list(rows[0].keys())
        
        # Match columns to categories based on hints
        hint_matched_cols = defaultdict(list)
        for category, hints in COLUMN_HINTS.items():
            for col in fieldnames:
                if any(hint in col.lower() for hint in hints):
                    hint_matched_cols[category].append(col)
        
        # Scan each column
        for col in fieldnames:
            values = [row.get(col, '') for row in rows if row.get(col)][:100]
            
            # Prioritize categories that match column hints
            categories = []
            for category in COLUMN_HINTS.keys():
                if col in hint_matched_cols.get(category, []):
                    categories.insert(0, category)
                else:
                    categories.append(category)
            
            for cat in categories:
                detected = _detect_type(values, cat)
                if detected:
                    detections[cat].update(detected)
                    column_mapping[col][cat].update(detected)
    
    return detections, column_mapping


def _scan_json(file_path, is_gz, is_lines):
    """Scan JSON data files"""
    opener = gzip.open if is_gz else open
    detections = defaultdict(set)
    column_mapping = defaultdict(lambda: defaultdict(set))
    
    with opener(file_path, 'rt', encoding='utf-8', errors='ignore') as f:
        objects = []
        
        if is_lines:
            for i, line in enumerate(f):
                if i >= 200:
                    break
                try:
                    objects.append(json.loads(line))
                except:
                    pass
        else:
            try:
                data = json.load(f)
                if isinstance(data, list):
                    objects = data[:200]
                elif isinstance(data, dict):
                    objects = [data]
            except:
                pass
        
        for obj in objects:
            if not isinstance(obj, dict):
                continue
            
            for key, value in obj.items():
                if not value:
                    continue
                
                if isinstance(value, str):
                    values = [value]
                elif isinstance(value, list):
                    values = value
                else:
                    values = []
                
                values = [str(v) for v in values if v][:100]
                key_lower = key.lower()
                
                # Prioritize categories that match key hints
                matching_categories = []
                for category, hints in COLUMN_HINTS.items():
                    if any(hint in key_lower for hint in hints):
                        matching_categories.insert(0, category)
                    else:
                        matching_categories.append(category)
                
                for cat in matching_categories:
                    detected = _detect_type(values, cat)
                    if detected:
                        detections[cat].update(detected)
                        column_mapping[key][cat].update(detected)
    
    return detections, column_mapping


def _detect_type(values, category):
    """Detect identifier types in a list of values"""
    if not values or category not in IDENTIFIER_PATTERNS:
        return set()
    
    detected = set()
    patterns = IDENTIFIER_PATTERNS[category]
    
    for name, pattern in patterns.items():
        matches = sum(1 for v in values[:50] if v and pattern.match(str(v).strip()))
        if matches >= 3:
            detected.add(name)
    
    return detected


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Detect biomedical identifiers in datasets',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s
  %(prog)s --input data/raw --output results.csv
  %(prog)s -i /path/to/data -o identifiers.csv
        '''
    )
    
    parser.add_argument(
        '-i', '--input',
        default='data/raw',
        help='Path to raw data directory (default: data/raw)'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='identifier_detection_results.csv',
        help='Output CSV file path (default: identifier_detection_results.csv)'
    )
    
    args = parser.parse_args()
    
    exit_code = detect_identifiers(args.input, args.output)
    sys.exit(exit_code)



