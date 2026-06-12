from abc import abstractmethod, ABC
from .taxonomy import Taxonomy, TaxonNode
import sys
import time
import os
import io
import numpy as np
import pandas as pd


class DataLoader(ABC):
    """
    Abstract base class for data loaders. It defines a standard interface
    for loading taxonomic data from different file formats.
    """
    def __init__(self, file_path_or_data):
        self.source = file_path_or_data
        self._initialize_taxonomies()

    @abstractmethod
    def _initialize_taxonomies(self) :
        """
        Loads and parses all data upon instantiation. Must be implemented by subclasses.
        """
        pass

    @abstractmethod
    def get_table(self, rank, metric_type='proportion', filter_taxon=None, lineage_type='name', include_taxid=True):
        """
        Returns a formatted DataFrame for a given rank from the pre-loaded data.

        Args:
            rank (str): The taxonomic rank to load.
            metric_type (str): The metric to use ('proportion' or 'count').
            filter_taxon (str, optional): A higher-level taxon to filter by. Defaults to None.
            lineage_type (str, optional): The type of lineage to retrieve ('name' or 'full_lineage'). Defaults to 'name'.
            include_taxid (bool, optional): Whether to include the 'taxid' column in the output. Defaults to True.
        """
        return pd.DataFrame()
        
class SingleSampleLoader(DataLoader):
    """단일 샘플/데이터셋을 로드하여 단일 Taxonomy 객체를 생성합니다."""
    def _initialize_taxonomies(self):
        self.taxonomy = self._parse_source(self.source)

    @abstractmethod
    def _parse_source(self, source_item):
        """소스를 파싱하여 하나의 Taxonomy 객체를 반환합니다."""
        pass

class MultiSampleLoader(DataLoader):
    
    @abstractmethod
    def _parse_source(self, source_item):
        """Parses a single source item and returns a Taxonomy object."""
        pass
    
    def _initialize_taxonomies(self):    
        """Loads data from all sources and creates a Taxonomy object for each."""
        self.taxonomies = {}
        sources = self.source if isinstance(self.source, list) else [self.source]

        for source_item in sources:
            sample_name = "string_input"
            if not (isinstance(source_item, str) and '\n' in source_item):
                sample_name = os.path.splitext(os.path.basename(source_item))[0].split('.')[0]
            
            tax_obj = self._parse_source(source_item)
            if tax_obj:
                self.taxonomies[sample_name] = tax_obj
                
    def get_table(self, rank, metric_type='proportion', filter_taxon=None, lineage_type='name', include_taxid=True): 
        """
        Filters the pre-loaded taxonomy objects and pivots the data to wide format.
        """
        all_samples_data = []
        
        for sample_name, tax_obj in self.taxonomies.items():
            start_node = tax_obj.root
            if filter_taxon:
                start_node = tax_obj.find_node(filter_taxon)
                if not start_node:
                    print(f"Warning: filter_taxon '{filter_taxon}' not found in sample '{sample_name}'.")
                    continue
            
            target_nodes = tax_obj.get_descendant_nodes(rank, start_node=start_node)
            
            for node in target_nodes:
                taxon_info = node.name
                if lineage_type == 'full_lineage':
                    lineage = tax_obj.get_full_lineage(node.name)
                    taxon_info = '|'.join(lineage) if lineage else node.name
                elif lineage_type == 'full_lineage_with_taxid':
                    lineage_with_taxids = tax_obj.get_full_lineage_with_taxids(node.name)
                    if lineage_with_taxids:
                        taxon_info = '|'.join([f"{name}({taxid})" for name, taxid in lineage_with_taxids])
                    else:
                        # Fallback to node name if lineage is not found
                        taxon_info = node.name

                all_samples_data.append({
                    'taxon': taxon_info,
                    'taxid': node.attributes.get('taxid', 'N/A'),
                    'sample': sample_name,
                    'proportion': node.attributes.get('proportion', 0.0),
                    'count': node.attributes.get('count', 0)
                })
                
        if not all_samples_data:
            return pd.DataFrame()
        
        long_df = pd.DataFrame(all_samples_data)
        
        if include_taxid:
            wide_df = long_df.pivot_table(
                index=['taxon', 'taxid'], 
                columns='sample', 
                values=metric_type, 
                fill_value=0
            )
        else:
            wide_df = long_df.pivot_table(index='taxon', columns='sample', values=metric_type, fill_value=0)

        return wide_df
        

class KrakenReportLoader(MultiSampleLoader):
    """Loads data specifically from Kraken-style reports, with hierarchical filtering."""
    def _parse_source(self, source_item):
        try:
            content = ""
            if isinstance(source_item, str) and '\n' in source_item:
                content = source_item
            else:
                content = open(source_item, 'r', encoding='utf-8').read()
            
            tax_obj = Taxonomy()
            node_stack = [(tax_obj.root, -1)]

            for line in content.strip().split('\n'):
                if line.strip().startswith('#') or not line.strip():
                    continue
                
                parts = line.split('\t')
                if len(parts) != 6: continue
                
                name_field = parts[5]
                indentation = len(name_field) - len(name_field.lstrip(' '))
                name = name_field.strip()
                rank = parts[3].strip()
                taxid = parts[4].strip()
                
                new_node = TaxonNode(name, rank)
                new_node.attributes['proportion'] = float(parts[0])
                new_node.attributes['count'] = int(parts[1])
                new_node.attributes['taxid'] = taxid

                while indentation <= node_stack[-1][1]:
                    node_stack.pop()
                
                parent_node = node_stack[-1][0]
                parent_node.children[name] = new_node
                node_stack.append((new_node, indentation))

            return tax_obj
        except Exception as e:
            print(f"An error occurred while parsing Kraken source '{source_item}': {e}")
            return None


class PhantaLoader(MultiSampleLoader):
    """
    Loads data from Phanta-style reports (wide format) into Taxonomy objects.
    """
    
    def _parse_source(self, source_item):
        pass
    
    def _initialize_taxonomies(self):
        """Overridden _load_all to handle wide format."""
        self.taxonomies = {}
        source_file = self.source[0] if isinstance(self.source, list) else self.source
        
        try:
            content = ""
            if isinstance(source_file, str) and '\n' in source_file:
                content = source_file
            else:
                with open(source_file, 'r', encoding='utf-8') as f:
                    content = f.read()

            df = pd.read_csv(io.StringIO(content), sep='\t')
            sample_columns = df.columns[2:]

            for sample in sample_columns:
                tax_obj = Taxonomy()
                for _, row in df.iterrows():
                    lineage_str = row.iloc[0]
                    lineage_id_str = row.iloc[1]
                    attributes = {
                        'proportion': pd.to_numeric(row[sample], errors='coerce'),
                        'count': pd.to_numeric(row[sample], errors='coerce')
                    }
                    self._add_lineage_to_taxonomy(tax_obj, lineage_str, lineage_id_str, attributes)
                self.taxonomies[sample] = tax_obj
        except Exception as e:
            print(f"An error occurred while parsing Phanta source '{source_file}': {e}")
    
    def _add_lineage_to_taxonomy(self, tax_obj, lineage_str, lineage_id_str, attributes):
        """Helper to add a single lineage to a Taxonomy object."""
        current_node = tax_obj.root
        lineage_ids = lineage_id_str.split('|')
        
        for i, segment in enumerate(lineage_str.split('|')):
            try:
                rank, name = segment.split('_', 1)
                taxid = lineage_ids[i] if i < len(lineage_ids) else 'N/A'
                if name not in current_node.children:
                    current_node.children[name] = TaxonNode(name, rank)
                    current_node.children[name].attributes['taxid'] = taxid
                current_node = current_node.children[name]
            except ValueError:
                continue
        if attributes:
            current_node.attributes.update(attributes)
    


# class VitaLoader(DataLoader):
#     """
#     Loads data from VITA-style reports (wide format, semicolon-separated, with name and ID lineages).
#     """
#     def _load_all(self):
#         self.taxonomies = {}
#         source_file = self.source[0] if isinstance(self.source, list) else self.source
        
#         try:
#             content = ""
#             if isinstance(source_file, str) and '\n' in source_file:
#                 content = source_file
#             else:
#                 with open(source_file, 'r', encoding='utf-8') as f:
#                     content = f.read()
            
#             df = pd.read_csv(io.StringIO(content), sep='\t')
#             # The third column onwards are samples
#             sample_columns = df.columns[2:]

#             for sample in sample_columns:
#                 tax_obj = Taxonomy()
#                 for _, row in df.iterrows():
#                     lineage_name_str = str(row.iloc[0])
#                     lineage_id_str = str(row.iloc[1])
                    
#                     count_val = pd.to_numeric(row[sample], errors='coerce')
#                     count = 0 if pd.isna(count_val) else int(count_val)
                    
#                     attributes = {
#                         'proportion': 1.0 if count > 0 else 0.0,
#                         'count': count,
#                     }
#                     self._add_lineage_to_taxonomy(tax_obj, lineage_name_str, lineage_id_str, attributes)
#                 self.taxonomies[sample] = tax_obj
#         except Exception as e:
#             print(f"An error occurred while parsing VITA source '{source_file}': {e}")

#     def _add_lineage_to_taxonomy(self, tax_obj, lineage_name_str, lineage_id_str, attributes):
#         rank_map = {
#             'su': 'superkingdom', 'ki': 'kingdom', 'ph': 'phylum',
#             'cl': 'class', 'or': 'order', 'fa': 'family',
#             'ge': 'genus', 'sp': 'species', 'no': 'no rank'
#         }
#         current_node = tax_obj.root
        
#         lineage_id_parts = [p.strip().split('_')[-1] for p in lineage_id_str.split(';')]
#         lineage_name_parts = lineage_name_str.split(';')

#         for i, segment in enumerate(lineage_name_parts):
#             segment = segment.strip()
#             if not segment: continue
            
#             parts = segment.split('_', 1)
#             if len(parts) == 2:
#                 rank_abbr, name = parts
#                 rank = rank_map.get(rank_abbr, 'no rank')
#             else:
#                 rank = 'no rank'
#                 name = segment
            
#             taxid = lineage_id_parts[i] if i < len(lineage_id_parts) else name
#             if name not in current_node.children:
#                 current_node.children[name] = TaxonNode(name, rank)
            
#             current_node = current_node.children[name]
#             current_node.attributes['taxid'] = taxid

#         current_node.attributes.update(attributes)
    
#     def get_table(self, rank, metric_type='proportion', filter_taxon=None, lineage_type='name', include_taxid_column=True):

#         all_samples_data = []
#         for sample_name, tax_obj in self.taxonomies.items():
#             start_node = tax_obj.root
#             if filter_taxon:
#                 start_node = tax_obj.find_node(filter_taxon)
#                 if not start_node:
#                     print(f"Warning: filter_taxon '{filter_taxon}' not found in sample '{sample_name}'.")
#                     continue
            
#             target_nodes = tax_obj.get_descendant_nodes(rank, start_node=start_node)

#             for node in target_nodes:
#                 if lineage_type == 'name':
#                     taxon_info = node.name
#                 elif lineage_type == 'full_lineage':
#                     taxon_info = '|'.join(tax_obj.get_full_lineage(node.name))
#                 elif lineage_type == 'full_lineage_with_taxid':
#                     lineage_with_taxids = tax_obj.get_full_lineage_with_taxids(node.name)
#                     taxon_info = '|'.join([f"{name}({taxid})" for name, taxid in lineage_with_taxids])
#                 else:
#                     taxon_info = node.name  # Default to name

#                 all_samples_data.append({
#                     'taxon': taxon_info,
#                     'taxid': node.attributes.get('taxid', 'N/A'),
#                     'sample': sample_name,

#                     'proportion': node.attributes.get('proportion', 0.0),
#                     'count': node.attributes.get('count', 0)
#                 })
        
#         if not all_samples_data:
#             return pd.DataFrame()
        

#         long_df = pd.DataFrame(all_samples_data)
        
#         wide_df = long_df.pivot_table(
#             index='taxon',
#             columns='sample', 
#             values=metric_type, 
#             fill_value=0,
#             aggfunc=np.sum
#         )

#         id_map = long_df.drop_duplicates(subset=['taxon'])[['taxon', 'taxid']].set_index('taxon')
#         final_df = wide_df.join(id_map)
#         cols = ['taxid'] + [col for col in final_df if col != 'taxid']
#         return final_df[cols]
    
# class GanonLoader(MultiSourceLoader):
#     """Loads data from Ganon-style reports, with hierarchical filtering."""
#     def _parse_source(self, source_item):
#         try:
#             content = ""
#             if isinstance(source_item, str) and '\n' in source_item:
#                 content = source_item
#             else:
#                 with open(source_item, 'r', encoding='utf-8') as f:
#                     content = f.read()
            
#             tax_obj = Taxonomy()
#             id_to_node = {'1': tax_obj.root}

#             lines = content.strip().split('\n')
#             if 'unclassified' in lines[0].lower():
#                 lines = lines[1:]

#             for line in lines:
#                 parts = line.split('\t')
#                 if len(parts) < 9: continue
#                 rank, target_id, _, name = parts[0], parts[1], parts[2], parts[3]
                
#                 if target_id in id_to_node:
#                      node_to_update = id_to_node[target_id]
#                 else:
#                      node_to_update = TaxonNode(name.strip(), rank.strip())
#                      id_to_node[target_id] = node_to_update
                
#                 node_to_update.attributes['proportion'] = float(parts[8])
#                 node_to_update.attributes['count'] = int(parts[7])
#                 node_to_update.attributes['taxid'] = target_id

#             for line in lines:
#                 parts = line.split('\t')
#                 if len(parts) < 9: continue
#                 target_id, lineage = parts[1], parts[2]
                
#                 if target_id == '1': continue
                
#                 lineage_ids = [lid for lid in str(lineage).split('|') if lid]
#                 parent_id = lineage_ids[-2] if len(lineage_ids) > 1 else '1'

#                 if parent_id in id_to_node and target_id in id_to_node:
#                     parent_node = id_to_node[parent_id]
#                     child_node = id_to_node[target_id]
#                     parent_node.children[child_node.name] = child_node
#                 # else:
#                 #     print(f"Warning: Could not link node {target_id} ({name}) to parent {parent_id}. Parent not found.")

#             return tax_obj
#         except Exception as e:
#             print(f"An error occurred while parsing Ganon source '{source_item}': {e}")
#             return None
