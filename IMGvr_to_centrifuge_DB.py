import pandas as pd

# Configuration
input_header = "/media/eternus1/nfs/projects/databases/IMG_VR/All_nucleotides.headers"     # Header mapping file
metadata_file = "/media/eternus1/nfs/projects/databases/IMG_VR/All_Sequence_information.tsv"  # Metadata file
conversion_table = "conversion_table.txt"   # Conversion table for centrifuge
name_table = "name_table.txt"               # Name table for centrifuge
taxonomy_tree = "taxonomy_tree.txt"         # Taxonomy tree for centrifuge

def load_headers(header_file):
    """Load the header file mapping sequence IDs."""
    print("Loading headers...")
    headers = pd.read_csv(header_file, sep="\t", header=None, names=["header", "offset"])
    headers["sequence_id"] = headers["header"].str.split("|").str[0]
    return headers

def load_metadata(metadata_file):
    """Load and parse the metadata file."""
    print("Loading metadata...")
    metadata = pd.read_csv(metadata_file, sep="\t", low_memory=False)
    required_cols = [
            "## UViG",
            "Taxon_oid",
            "Taxonomic classification",
        ]
    metadata = metadata[required_cols]
    metadata['## UViG']=metadata[['## UViG', 'Taxon_oid']].agg('|'.join, axis=1)
    metadata = metadata[['## UViG', 'Taxonomic classification']].fillna('Unclassified')
    metadata = metadata.replace(r'\|Reference', '', regex=True)
    metadata.columns = ["sequence_id", "taxonomy"]
    metadata['sequence_id'] = metadata['sequence_id'].apply(lambda x: '|'.join(x.split('|')[:2]))
    metadata['taxonomy'] = metadata['taxonomy'].str.rstrip(';')
    return metadata

def parse_taxonomy(metadata):
    """Parse taxonomy information to create taxonomy IDs and hierarchy, handling skipped levels."""
    print("Parsing taxonomy...")
    taxonomy_hierarchy = {}
    taxonomy_id_map = {}
    current_id = 2  # Start taxonomy IDs from 2
    
    # Create taxonomy ID mapping and hierarchy
    for taxonomy in metadata["taxonomy"].dropna().unique():
        parts = taxonomy.split(";")
        parent_id = 1  # Root
        
        # Process each level in the taxonomy
        for i in range(len(parts)):
            sub_taxonomy = ";".join(parts[:i + 1])
            
            # Skip empty levels
            if not parts[i].strip():
                continue
            
            # Add unique taxonomy levels
            if sub_taxonomy not in taxonomy_id_map:
                taxonomy_id_map[sub_taxonomy] = current_id
                taxonomy_hierarchy[current_id] = [parent_id] if parent_id else [1]  # Ensure list
                parent_id = current_id  # Update parent ID for next level
                current_id += 1
            else:
                parent_id = taxonomy_id_map[sub_taxonomy]  # Use existing ID for parent
                #taxonomy_hierarchy[parent_id].append(current_id)  # Add current to the parent's list
    
    return taxonomy_id_map, taxonomy_hierarchy

def write_output(metadata, taxonomy_id_map, taxonomy_hierarchy, conversion_table, name_table, taxonomy_tree):
    """Write the conversion table, name table, and taxonomy tree."""
    print("Writing output files...")
    
    # Add taxonomy IDs to metadata
    metadata["taxonomy_id"] = metadata["taxonomy"].map(taxonomy_id_map)

    # Write conversion table
    with open(conversion_table, "w") as conv_file:
        for _, row in metadata.iterrows():
            conv_file.write(f"{row['sequence_id']}\t{row['taxonomy_id']}\n")
    print(f"Conversion table written to: {conversion_table}")

    # Write name table
    with open(name_table, "w") as name_file:
        name_file.write("1\t|\troot\t|\t\t|\tscientific name\t|\n")
        for taxonomy, tax_id in taxonomy_id_map.items():
            name_file.write(f"{tax_id}\t|\t{taxonomy}\t|\t\t|\tscientific name\t|\n")
    print(f"Name table written to: {name_table}")

    # Write taxonomy tree
    with open(taxonomy_tree, "w") as tree_file:
        tree_file.write('1\t|\t1\t|\troot\n')
        for tax_id, parent_ids in taxonomy_hierarchy.items():
            parent_id = parent_ids[-1] if parent_ids else 1  # Default to root ID (1)
            tree_file.write(f"{tax_id}\t|\t{parent_id}\t|\tno rank\n")
    print(f"Taxonomy tree written to: {taxonomy_tree}")

def main():
    # Step 1: Load headers and metadata
    headers = load_headers(input_header)
    metadata = load_metadata(metadata_file)

    # Step 2: Parse taxonomy information
    taxonomy_id_map, taxonomy_hierarchy = parse_taxonomy(metadata)

    # Step 3: Write output files
    write_output(metadata, taxonomy_id_map, taxonomy_hierarchy, conversion_table, name_table, taxonomy_tree)

if __name__ == "__main__":
    main()