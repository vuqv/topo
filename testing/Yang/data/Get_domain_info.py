import requests

def get_cath_domains(pdb_id):
    """Retrieve CATH domains for a given PDB ID."""
    url = f"http://www.cathdb.info/api/rest/domain_list/pdb/{pdb_id.lower()}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()  # Returns a list of CATH domains for the PDB
    else:
        print(f"Failed to retrieve CATH domains for {pdb_id}")
        return None

def get_uniprot_sequence(uniprot_id):
    """Retrieve the canonical sequence from UniProt given a UniProt ID."""
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        fasta_data = response.text
        sequence = ''.join(fasta_data.splitlines()[1:])  # Remove the FASTA header
        return sequence
    else:
        print(f"Failed to retrieve UniProt sequence for {uniprot_id}")
        return None

def get_sifts_mapping(pdb_id, chain_id):
    """Retrieve PDB to UniProt residue mappings for a given PDB ID and chain."""
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        # Extract the mapping for the specified chain
        mappings = data[pdb_id.lower()]['UniProt']
        for uniprot_id, details in mappings.items():
            for chain_mapping in details['mappings']:
                if chain_mapping['chain_id'] == chain_id:
                    return chain_mapping  # Returns the UniProt <-> PDB residue mappings
    else:
        print(f"Failed to retrieve SIFTS mapping for {pdb_id} chain {chain_id}")
        return None

def map_cath_to_uniprot(cath_domains, sifts_mapping, canonical_sequence):
    """Map CATH domains to the UniProt canonical sequence."""
    uniprot_domains = []
    for domain in cath_domains:
        pdb_start = domain['residues']['start']['res_num']
        pdb_end = domain['residues']['end']['res_num']
        
        # Map PDB residues to UniProt residues using the SIFTS mapping
        uniprot_start = sifts_mapping['unp_start']
        uniprot_end = sifts_mapping['unp_end']
        
        # Extract the corresponding UniProt sequence
        domain_sequence = canonical_sequence[uniprot_start-1:uniprot_end]
        uniprot_domains.append({
            'cath_id': domain['domain'],
            'uniprot_start': uniprot_start,
            'uniprot_end': uniprot_end,
            'sequence': domain_sequence
        })
    return uniprot_domains

def get_cath_domains_mapped_to_uniprot(uniprot_id, pdb_id, chain_id):
    """Retrieve CATH domains and map them to the UniProt canonical sequence."""
    # Step 1: Retrieve CATH domains for the PDB ID
    cath_domains = get_cath_domains(pdb_id)
    if not cath_domains:
        return None

    # Step 2: Retrieve the canonical sequence from UniProt
    canonical_sequence = get_uniprot_sequence(uniprot_id)
    if not canonical_sequence:
        return None

    # Step 3: Retrieve SIFTS mapping between PDB and UniProt residues
    sifts_mapping = get_sifts_mapping(pdb_id, chain_id)
    if not sifts_mapping:
        return None

    # Step 4: Map the CATH domains to the UniProt canonical sequence
    uniprot_domains = map_cath_to_uniprot(cath_domains, sifts_mapping, canonical_sequence)

    return uniprot_domains

# Example usage
uniprot_id = "P77214"  # Replace with your UniProt ID
pdb_id = "1ZEQ"  # Replace with your PDB ID
chain_id = "X"  # Replace with your chain ID

mapped_domains = get_cath_domains_mapped_to_uniprot(uniprot_id, pdb_id, chain_id)
if mapped_domains:
    for domain in mapped_domains:
        print(f"CATH Domain: {domain['cath_id']}, UniProt Start: {domain['uniprot_start']}, UniProt End: {domain['uniprot_end']}")
        print(f"Domain Sequence: {domain['sequence']}")
else:
    print("Failed to retrieve or map CATH domains.")

