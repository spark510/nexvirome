
class TaxonNode:
    """Represents a single node in a taxonomy tree with flexible attributes."""
    def __init__(self, name, rank):
        self.name = name
        self.rank = rank
        self.children = {}  # {child_name: TaxonNode}
        self.attributes = {} # Flexible dictionary for additional data
        
    def __str__(self) -> str:
        return f"TaxonNode(name={self.name}, rank={self.rank}, attributes={self.attributes}), children={len(self.children)})"


class Taxonomy:
    """Manages the hierarchical structure of a taxonomy."""
    def __init__(self):
        self.root = TaxonNode("root", "no rank")
        self.root.attributes['taxid'] = '1'
        
    def find_node(self, taxon_name, start_node=None):
        """Recursively finds a node by its name, starting from a given node or the root."""
        if start_node is None:
            start_node = self.root
        
        if start_node.name.lower() == taxon_name.lower():
            return start_node
        
        for child in start_node.children.values():
            found_node = self.find_node(taxon_name, child)
            if found_node:
                return found_node
        return None
    
    def get_descendant_nodes(self, rank, start_node=None):
        """Recursively gets all descendant nodes of a specific rank."""
        if start_node is None:
            start_node = self.root
        
        descendant_nodes = []
        if start_node.rank.lower() == rank.lower():
            descendant_nodes.append(start_node)
        
        for child in start_node.children.values():
            descendant_nodes.extend(self.get_descendant_nodes(rank, child))
        
        return descendant_nodes
    
    def print_tree(self, start_node=None, prefix=""):
        """Prints the taxonomy tree structure for debugging purposes."""
        if start_node is None:
            start_node = self.root
            print(f"{start_node.rank}: {start_node.name} {start_node.attributes}")

        for i, child_node in enumerate(start_node.children.values()):
            connector = "└── " if i == len(start_node.children) - 1 else "├── "
            print(f"{prefix}{connector}{child_node.rank}: {child_node.name} {child_node.attributes}")
            if child_node.children:
                new_prefix = prefix + ("    " if i == len(start_node.children) - 1 else "│   ")
                self.print_tree(child_node, new_prefix)

    def find_path_to_node(self, taxon_name, start_node=None):
        """Recursively finds a path of nodes to a taxon by its name."""
        if start_node is None:
            start_node = self.root

        if start_node.name.lower() == taxon_name.lower():
            return [start_node]

        for child in start_node.children.values():
            path = self.find_path_to_node(taxon_name, child)
            if path:
                return [start_node] + path
        return None

    def get_full_lineage(self, taxon_name):
        """Returns the full lineage of a taxon as a list of names."""
        path = self.find_path_to_node(taxon_name)
        if not path:
            return []
        # Exclude root from lineage string unless it's the only node
        return [node.name for node in (path[1:] if len(path) > 1 else path)]

    def get_full_lineage_with_taxids(self, taxon_name):
        """Returns the full lineage of a taxon with taxids as a list of tuples (name, taxid)."""
        path = self.find_path_to_node(taxon_name)
        if not path:
            return []
        # Exclude root from lineage string unless it's the only node
        path_to_return = path[1:] if len(path) > 1 else path
        return [(node.name, node.attributes.get('taxid', 'N/A')) for node in path_to_return]