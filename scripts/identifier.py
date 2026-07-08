"""
STEP 4: Biomedical Identifier Auto-Detection
Comprehensive detection of all major biomedical identifier systems
Uses verified biological entity keywords
Outputs to CSV file
"""

import re
import json
import gzip
import csv
from pathlib import Path
from collections import defaultdict
import sys

# Verified biological entity keywords
COLUMN_HINTS = {
    'gene': [
        'gene', 'hgnc', 'symbol', 'ensembl', 'gene_id', 'entrez', 'geneid', 'ncbi_gene', 
        'refseq', 'mgi', 'rgd', 'transcript', 'ensg', 'ensembl_gene', 'ensembl_transcript'
    ],
    'disease': [
        'disease', 'disorder', 'condition', 'mondo', 'doid', 'omim', 'mesh', 
        'icd', 'icd10', 'icd11', 'umls', 'orphanet', 'ordo', 'snomed', 'snomedct'
    ],
    'phenotype': [
        'phenotype', 'hpo', 'hp', 'human_phenotype', 'mp', 'mammalian_phenotype',
        'trait', 'zp', 'zebrafish_phenotype', 'pat', 'phenotypic_abnormality'
    ],
    'therapy': [
        'therapy', 'treatment', 'medical_action', 'maxo', 'procedure', 'intervention',
        'therapeutic_procedure', 'preventative_therapy', 'palliative_care', 
        'complementary_therapy', 'alternative_therapy', 'medical_procedure'
    ],
    'protein': [
        'protein', 'uniprot', 'uniprotkb', 'ensembl', 'refseq', 'pdb', 'interpro', 
        'pfam', 'accession', 'uniprot_ac', 'uniprot_id', 'ensembl_protein'
    ],
    'compound': [
        'drug', 'compound', 'chemical', 'chembl', 'drugbank', 'pubchem', 'pubchem_cid', 
        'kegg', 'kegg_compound', 'chebi', 'inchi', 'inchikey', 'smiles', 'rxcui'
    ],
    'variant': [
        'variant', 'snp', 'rsid', 'dbsnp', 'clinvar', 'cosmic', 'mutation', 'hgvs', 
        'clinvar_id', 'cosmic_id', 'gnomad', 'af', 'allele_frequency'
    ],
    'pathway': [
        'pathway', 'kegg', 'reactome', 'go', 'gene_ontology', 'biological_process', 
        'wiki_pathway', 'wp', 'pathway_id'
    ],
    'anatomy': [
        'anatomy', 'tissue', 'organ', 'uberon', 'cell', 'cl', 'cell_type', 'fma', 
        'cell_ontology', 'efo', 'uberon_id'
    ],
    'taxonomy': [
        'taxonomy', 'taxon', 'species', 'organism', 'ncbi', 'ncbitaxon', 'taxid', 
        'taxonomy_id'
    ]
}

# VALIDATED IDENTIFIER PATTERNS
IDENTIFIER_PATTERNS = {
    'gene': {
        'HGNC_Symbol': re.compile(r'^[A-Z][A-Z0-9\-]{1,9}$'),
        'HGNC': re.compile(r'^HGNC:\d+$'),
        'Ensembl_Gene': re.compile(r'^ENSG\d{11}$'),
        'Ensembl_Transcript': re.compile(r'^ENST\d{11}$'),
        'Entrez_Gene': re.compile(r'^\d{3,9}$'),
        'RefSeq_Gene': re.compile(r'^(NM_|NR_|XM_|XR_)\d+(\.\d+)?$'),
        'MGI': re.compile(r'^MGI:\d+$'),
        'RGD': re.compile(r'^RGD:\d+$'),
    },
    'disease': {
        'MONDO': re.compile(r'^MONDO:\d{7}$'),
        'DOID': re.compile(r'^DOID:\d+$'),
        'OMIM': re.compile(r'^\d{6}$'),
        'MeSH': re.compile(r'^[CD]\d{6,7}$'),
        'ICD10': re.compile(r'^[A-TV-Z]\d{2}(\.\d{1,4})?$'),
        'ICD11': re.compile(r'^[0-9A-Z]{2,10}$'),
        'UMLS': re.compile(r'^C\d{7}$'),
        'Orphanet': re.compile(r'^ORPHA:\d+$'),
        'SNOMED': re.compile(r'^\d{6,18}$'),
    },
    'phenotype': {
        'HPO': re.compile(r'^HP:\d{7}$'),
        'MP': re.compile(r'^MP:\d{7}$'),
        'ZP': re.compile(r'^ZP:\d{7}$'),
    },
    'therapy': {
        'MAXO': re.compile(r'^MAXO:\d{7}$'),
    },
    'protein': {
        'UniProt': re.compile(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$'),
        'Ensembl_Protein': re.compile(r'^ENSP\d{11}$'),
        'RefSeq_Protein': re.compile(r'^(NP_|XP_|YP_|WP_)\d+(\.\d+)?$'),
        'PDB': re.compile(r'^[0-9][A-Z0-9]{3}$'),
        'InterPro': re.compile(r'^IPR\d{6}$'),
        'Pfam': re.compile(r'^PF\d{5}$'),
    },
    'compound': {
        'DrugBank': re.compile(r'^DB\d{5}$'),
        'ChEMBL': re.compile(r'^CHEMBL\d+$'),
        'PubChem_CID': re.compile(r'^\d{4,9}$'),
        'KEGG_Compound': re.compile(r'^C\d{5}$'),
        'ChEBI': re.compile(r'^CHEBI:\d+$'),
        'InChI': re.compile(r'^InChI=1S?/[A-Za-z0-9\.]+'),
        'InChIKey': re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$'),
        'SMILES': re.compile(r'^[CNOPSFIBrClcnops0-9\[\]\(\)=#@\+\-\\\/\.]{10,}$'),
        'RxCUI': re.compile(r'^[1-9]\d{4,6}$'),
    },
    'variant': {
        'dbSNP': re.compile(r'^rs\d+$'),
        'ClinVar_RCV': re.compile(r'^RCV\d{9}(\.\d+)?$'),
        'ClinVar_VCV': re.compile(r'^VCV\d{9}(\.\d+)?$'),
        'COSMIC': re.compile(r'^COSM\d+$'),
        'HGVS': re.compile(r'^(c|g|m|n|p|r)\.[A-Z0-9*\-_>]+'),
    },
    'pathway': {
        'KEGG_Pathway': re.compile(r'^(hsa|mmu|rno|dme|cel|sce)\d{5}$'),
        'Reactome': re.compile(r'^R-[A-Z]{3}-\d+(-\d+)?$'),
        'GO': re.compile(r'^GO:\d{7}$'),
        'WikiPathways': re.compile(r'^WP\d+$'),
    },
    'anatomy': {
        'UBERON': re.compile(r'^UBERON:\d{7}$'),
        'CL': re.compile(r'^CL:\d{7}$'),
        'FMA': re.compile(r'^FMA:\d+$'),
        'EFO': re.compile(r'^EFO:\d{7}$'),
    },
    'taxonomy': {
        'NCBITaxon': re.compile(r'^(NCBITaxon:)?\d{1,7}$'),
    }
}

def detect_identifiers(raw_data_path='data/raw', output_path='data/registry/identifier_candidates.csv'):
    """Main detection pipeline - outputs to CSV"""
    
    raw_path = Path(raw_data_path)
    
    # Verify input directory exists
    if not raw_path.exists():
        print(f"❌ Error: Directory not found: {raw_data_path}")
        print(f"   Current directory: {Path.cwd()}")
        sys.exit(1)
    
    datasets = [d for d in raw_path.iterdir() if d.is_dir() and not d.name.startswith('.')]
    
    if not datasets:
        print(f"❌ Error: No datasets found in {raw_data_path}")
        sys.exit(1)
    
    print("="*80)
    print("COMPREHENSIVE BIOMEDICAL IDENTIFIER DETECTION")
    print("="*80)
    print()
    
    results = []
    
    for dataset_path in sorted(datasets):
        dataset_id = dataset_path.name
        print(f"📁 Dataset: {dataset_id}")
        
        all_detections = defaultdict(set)
        column_info = defaultdict(lambda: defaultdict(set))
        
        file_count = 0
        for file_path in dataset_path.rglob('*'):
            if not file_path.is_file():
                continue
            
            ext = file_path.suffix.lower()
            if ext == '.gz':
                ext = '.' + file_path.stem.split('.')[-1] if '.' in file_path.stem else ''
            
            if ext not in ['.csv', '.tsv', '.json', '.jsonl', '.ndjson']:
                continue
            
            file_count += 1
            
            try:
                detections, col_mapping = _scan_file(file_path)
                for category, ids in detections.items():
                    all_detections[category].update(ids)
                
                for col_name, col_detections in col_mapping.items():
                    for category, ids in col_detections.items():
                        for id_type in ids:
                            column_info[col_name][category].add(id_type)
                            
            except Exception as e:
                print(f"  ⚠️  Error scanning {file_path.name}: {e}")
        
        print(f"  Files scanned: {file_count}")
        
        # Build result row
        result_row = {'dataset_id': dataset_id}
        
        for category in ['gene', 'disease', 'phenotype', 'therapy', 'protein', 'compound', 'variant', 'pathway', 'anatomy', 'taxonomy']:
            column_name = f'possible_{category}_ids'
            if category in all_detections and all_detections[category]:
                result_row[column_name] = '|'.join(sorted(all_detections[category]))
                print(f"  ✓ {category.upper()}: {result_row[column_name]}")
            else:
                result_row[column_name] = 'NA'
        
        # Add column-level summary
        col_summary = []
        for col_name, categories in list(column_info.items())[:5]:
            for category, ids in categories.items():
                col_summary.append(f"{col_name}({','.join(ids)})")
        result_row['example_columns'] = '; '.join(col_summary) if col_summary else 'NA'
        
        results.append(result_row)
        print()
    
    # Create output directory
    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Write CSV
    fieldnames = ['dataset_id', 'possible_gene_ids', 'possible_disease_ids', 'possible_phenotype_ids', 
                  'possible_therapy_ids', 'possible_protein_ids', 'possible_compound_ids', 
                  'possible_variant_ids', 'possible_pathway_ids', 'possible_anatomy_ids', 
                  'possible_taxonomy_ids', 'example_columns']
    
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    print("="*80)
    print(f"✓ Detection complete!")
    print(f"✓ Results saved to: {output_file}")
    print(f"✓ Total datasets: {len(results)}")
    print("="*80)

def _scan_file(file_path):
    """Scan single file and return detected identifier types + column mapping"""
    
    ext = file_path.suffix.lower()
    is_gz = ext == '.gz'
    
    if is_gz:
        ext = '.' + file_path.stem.split('.')[-1] if '.' in file_path.stem else ''
    
    if ext in ['.csv', '.tsv']:
        return _scan_tabular(file_path, is_gz, '\t' if ext == '.tsv' else ',')
    elif ext in ['.json', '.jsonl', '.ndjson']:
        return _scan_json(file_path, is_gz, ext in ['.jsonl', '.ndjson'])
    
    return {}, {}

def _scan_tabular(file_path, is_gz, delimiter):
    """Fast tabular file scanner - TWO-PHASE detection"""
    
    opener = gzip.open if is_gz else open
    detections = defaultdict(set)
    column_mapping = defaultdict(lambda: defaultdict(set))
    
    try:
        with opener(file_path, 'rt', encoding='utf-8', errors='ignore') as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            
            try:
                first_row = next(reader)
                fieldnames = list(first_row.keys())
            except:
                return detections, column_mapping
            
            # PHASE 1: Check column names for hints
            hint_matched_cols = defaultdict(list)
            for category, hints in COLUMN_HINTS.items():
                for col in fieldnames:
                    col_lower = col.lower()
                    if any(hint in col_lower for hint in hints):
                        hint_matched_cols[category].append(col)
            
            # Sample up to 100 rows
            sample_rows = [first_row]
            for i, row in enumerate(reader):
                if i >= 99:
                    break
                sample_rows.append(row)
            
            # PHASE 2: Pattern matching
            all_categories = list(COLUMN_HINTS.keys())
            
            for col in fieldnames:
                values = [row.get(col, '') for row in sample_rows if row.get(col)][:50]
                
                categories_to_check = []
                for category in all_categories:
                    if col in hint_matched_cols.get(category, []):
                        categories_to_check.insert(0, category)
                    else:
                        categories_to_check.append(category)
                
                for category in categories_to_check:
                    detected = _detect_type(values, category)
                    if detected:
                        detections[category].update(detected)
                        column_mapping[col][category].update(detected)
    except Exception as e:
        pass
    
    return detections, column_mapping

def _scan_json(file_path, is_gz, is_lines):
    """Fast JSON scanner - TWO-PHASE detection"""
    
    opener = gzip.open if is_gz else open
    detections = defaultdict(set)
    column_mapping = defaultdict(lambda: defaultdict(set))
    
    try:
        with opener(file_path, 'rt', encoding='utf-8', errors='ignore') as f:
            if is_lines:
                for i, line in enumerate(f):
                    if i >= 100:
                        break
                    try:
                        obj = json.loads(line)
                        obj_detections, obj_cols = _check_json_obj(obj)
                        for category, ids in obj_detections.items():
                            detections[category].update(ids)
                        for col_name, col_dets in obj_cols.items():
                            for category, ids in col_dets.items():
                                column_mapping[col_name][category].update(ids)
                    except:
                        pass
            else:
                try:
                    data = json.load(f)
                    if isinstance(data, list):
                        for obj in data[:100]:
                            obj_detections, obj_cols = _check_json_obj(obj)
                            for category, ids in obj_detections.items():
                                detections[category].update(ids)
                            for col_name, col_dets in obj_cols.items():
                                for category, ids in col_dets.items():
                                    column_mapping[col_name][category].update(ids)
                    elif isinstance(data, dict):
                        obj_detections, obj_cols = _check_json_obj(data)
                        for category, ids in obj_detections.items():
                            detections[category].update(ids)
                        for col_name, col_dets in obj_cols.items():
                            for category, ids in col_dets.items():
                                column_mapping[col_name][category].update(ids)
                except:
                    pass
    except Exception as e:
        pass
    
    return detections, column_mapping

def _check_json_obj(obj):
    """Check single JSON object - TWO-PHASE"""
    
    if not isinstance(obj, dict):
        return {}, {}
    
    detections = defaultdict(set)
    column_mapping = defaultdict(lambda: defaultdict(set))
    
    for key, value in obj.items():
        if not value:
            continue
        
        key_lower = key.lower()
        values = [value] if isinstance(value, str) else (value if isinstance(value, list) else [])
        values = [str(v) for v in values if v][:50]
        
        matching_categories = []
        for category, hints in COLUMN_HINTS.items():
            if any(hint in key_lower for hint in hints):
                matching_categories.insert(0, category)
            else:
                matching_categories.append(category)
        
        for category in matching_categories:
            detected = _detect_type(values, category)
            if detected:
                detections[category].update(detected)
                column_mapping[key][category].update(detected)
    
    return detections, column_mapping

def _detect_type(values, category):
    """Pattern-based type detection - requires 3+ matches"""
    
    if not values or category not in IDENTIFIER_PATTERNS:
        return set()
    
    detected = set()
    patterns = IDENTIFIER_PATTERNS[category]
    
    for name, pattern in patterns.items():
        matches = 0
        for v in values[:30]:
            if v and pattern.match(str(v).strip()):
                matches += 1
                if matches >= 3:
                    detected.add(name)
                    break
    
    return detected

if __name__ == '__main__':
    detect_identifiers()

