import os
import gzip
import json
import csv
from pathlib import Path
from collections import defaultdict
import pandas as pd
from typing import Set, Dict, Any, Tuple, List
import xml.etree.ElementTree as ET

class DatasetScanner:
    """Memory-efficient dataset scanner for large-scale biomedical data."""
    
    def __init__(self, raw_data_path: Path):
        self.raw_data_path = raw_data_path
        self.registry_path = raw_data_path.parent / "registry"
        self.registry_path.mkdir(exist_ok=True)
    
    def scan_all_datasets(self) -> pd.DataFrame:
        """Main entry point: scan all datasets and build registry."""
        
        dataset_folders = self._discover_datasets()
        total_datasets = len(dataset_folders)
        
        print(f"Auto-detected {total_datasets} datasets in {self.raw_data_path}")
        print("="*80)
        
        datasets = []
        for idx, dataset_folder in enumerate(dataset_folders, 1):
            print(f"[{idx}/{total_datasets}] Scanning {dataset_folder.name}...")
            metadata = self.scan_dataset(dataset_folder)
            datasets.append(metadata)
        
        df = pd.DataFrame(datasets)
        output_path = self.registry_path / "dataset_registry.csv"
        df.to_csv(output_path, index=False)
        
        print("="*80)
        print(f"✓ Registry created: {output_path}")
        print(f"✓ Total datasets processed: {len(df)}")
        
        return df
    
    def _discover_datasets(self) -> List[Path]:
        """Auto-detect all dataset folders."""
        dataset_folders = [
            folder for folder in self.raw_data_path.iterdir()
            if folder.is_dir() and not folder.name.startswith('.')
        ]
        return sorted(dataset_folders)
    
    def scan_dataset(self, dataset_path: Path) -> Dict[str, Any]:
        """Scan a single dataset folder and extract metadata."""
        dataset_id = dataset_path.name
        
        files = self._discover_files(dataset_path)
        
        if not files:
            print(f"  ⚠ Warning: No valid files found in {dataset_id}")
            return self._empty_metadata(dataset_id)
        
        print(f"  Found {len(files)} file(s)")
        
        file_formats = set()
        total_size_bytes = 0
        total_rows = 0
        all_columns = set()
        
        for file_path in files:
            total_size_bytes += file_path.stat().st_size
            format_info = self._get_file_format(file_path)
            file_formats.add(format_info['extension'])
            
            try:
                rows, columns = self._parse_file(file_path, format_info)
                total_rows += rows
                all_columns.update(columns)
            except Exception as e:
                print(f"  ⚠ Warning: Could not parse {file_path.name}: {e}")
        
        return {
            'dataset_id': dataset_id,
            'file_formats': ','.join(sorted(file_formats)),
            'file_size_mb': round(total_size_bytes / (1024 * 1024), 2),
            'row_count': total_rows if total_rows > 0 else 'NA',
            'num_columns': len(all_columns),
            'column_names': json.dumps(sorted(list(all_columns)))
        }
    
    def _discover_files(self, dataset_path: Path) -> List[Path]:
        """Auto-detect all valid files in a dataset folder."""
        supported_extensions = {
            'csv', 'tsv', 'json', 'jsonl', 'ndjson', 'xml', 'obo', 'gz'
        }
        
        files = []
        for file_path in dataset_path.rglob("*"):
            if not file_path.is_file():
                continue
            
            suffixes = [s[1:] for s in file_path.suffixes]
            if any(ext in supported_extensions for ext in suffixes):
                files.append(file_path)
        
        return sorted(files)
    
    def _empty_metadata(self, dataset_id: str) -> Dict[str, Any]:
        """Return empty metadata for datasets with no valid files."""
        return {
            'dataset_id': dataset_id,
            'file_formats': 'none',
            'file_size_mb': 0.0,
            'row_count': 'NA',
            'num_columns': 0,
            'column_names': '[]'
        }
    
    def _get_file_format(self, file_path: Path) -> Dict[str, Any]:
        """Determine file format, handling compression."""
        suffixes = file_path.suffixes
        is_compressed = suffixes[-1] == '.gz' if suffixes else False
        
        if is_compressed and len(suffixes) >= 2:
            extension = suffixes[-2][1:]
        elif suffixes:
            extension = suffixes[-1][1:]
        else:
            extension = 'unknown'
        
        return {'extension': extension, 'is_compressed': is_compressed, 'path': file_path}
    
    def _parse_file(self, file_path: Path, format_info: Dict) -> Tuple[int, Set[str]]:
        """Parse a file and return (row_count, column_names)."""
        ext = format_info['extension']
        is_compressed = format_info['is_compressed']
        
        if ext in ['csv', 'tsv']:
            return self._parse_csv_tsv(file_path, ext, is_compressed)
        elif ext == 'json':
            return self._parse_json(file_path, is_compressed)
        elif ext in ['jsonl', 'ndjson']:
            return self._parse_jsonl(file_path, is_compressed)
        elif ext == 'xml':
            return self._parse_xml(file_path, is_compressed)
        elif ext == 'obo':
            return self._parse_obo(file_path, is_compressed)
        else:
            return 0, set()
    
    def _open_file(self, file_path: Path, is_compressed: bool, mode='rt'):
        """Open file, handling compression transparently."""
        if is_compressed:
            return gzip.open(file_path, mode, encoding='utf-8', errors='ignore')
        else:
            return open(file_path, mode, encoding='utf-8', errors='ignore')
    
    def _parse_csv_tsv(self, file_path: Path, ext: str, is_compressed: bool):
        delimiter = '\t' if ext == 'tsv' else ','
        
        with self._open_file(file_path, is_compressed) as f:
            reader = csv.reader(f, delimiter=delimiter)
            try:
                header = next(reader)
                columns = set(header)
            except StopIteration:
                return 0, set()
            
            row_count = sum(1 for _ in reader)
        
        return row_count, columns
    
    def _parse_json(self, file_path: Path, is_compressed: bool):
        with self._open_file(file_path, is_compressed) as f:
            try:
                data = json.load(f)
            except json.JSONDecodeError:
                return 0, set()
        
        columns = set()
        row_count = 0
        
        if isinstance(data, list):
            row_count = len(data)
            for item in data[:1000]:
                if isinstance(item, dict):
                    columns.update(item.keys())
        elif isinstance(data, dict):
            columns.update(data.keys())
            row_count = 1
        
        return row_count, columns
    
    def _parse_jsonl(self, file_path: Path, is_compressed: bool):
        columns = set()
        row_count = 0
        max_sample = 10000
        
        with self._open_file(file_path, is_compressed) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                try:
                    obj = json.loads(line)
                    if isinstance(obj, dict) and row_count < max_sample:
                        columns.update(obj.keys())
                    row_count += 1
                except json.JSONDecodeError:
                    continue
        
        return row_count, columns
    
    def _parse_xml(self, file_path: Path, is_compressed: bool):
        columns = set()
        max_tags = 1000
        
        try:
            with self._open_file(file_path, is_compressed) as f:
                for event, elem in ET.iterparse(f, events=('start',)):
                    if len(columns) < max_tags:
                        columns.add(elem.tag)
                    elem.clear()
        except ET.ParseError:
            pass
        
        return 0, columns
    
    def _parse_obo(self, file_path: Path, is_compressed: bool):
        columns = set()
        row_count = 0
        
        with self._open_file(file_path, is_compressed) as f:
            for line in f:
                line = line.strip()
                
                if line.startswith('['):
                    row_count += 1
                
                if ':' in line and not line.startswith('['):
                    key = line.split(':', 1)[0].strip()
                    if key:
                        columns.add(key)
        
        return row_count, columns


if __name__ == "__main__":
    script_dir = Path.cwd()
    raw_data_path = script_dir / "data" / "raw"
    
    if not raw_data_path.exists():
        raise FileNotFoundError(
            f"Expected data/raw/ folder not found at: {raw_data_path}\n"
            f"Current directory: {script_dir}"
        )
    
    subdirs = [d for d in raw_data_path.iterdir() if d.is_dir()]
    if len(subdirs) == 0:
        print(f"⚠ WARNING: {raw_data_path} exists but is empty!")
        print(f"\nRun 'python sample_data_generator.py' to create test data.")
        exit(1)
    
    scanner = DatasetScanner(raw_data_path)
    registry_df = scanner.scan_all_datasets()
    
    print("\n" + "="*80)
    print("REGISTRY SUMMARY")
    print("="*80)
    
    if len(registry_df) > 0:
        print(registry_df.to_string(index=False))
        print("\n" + "="*80)
        print(f"Total datasets: {len(registry_df)}")
        
        if 'file_size_mb' in registry_df.columns:
            print(f"Total size: {registry_df['file_size_mb'].sum():.2f} MB")
        
        if 'row_count' in registry_df.columns:
            numeric_rows = registry_df[registry_df['row_count'] != 'NA']['row_count']
            if len(numeric_rows) > 0:
                print(f"Total rows: {numeric_rows.sum():,}")
        
        print(f"Output: data/registry/dataset_registry.csv")

