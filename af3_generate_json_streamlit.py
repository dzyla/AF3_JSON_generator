from ast import mod
import streamlit as st
import requests
import json
import datetime

ptm_mapping = {
        'Phosphoserine': ('CCD_SEP', 'S'),
        'Phosphothreonine': ('CCD_TPO', 'T'),
        'Phosphotyrosine': ('CCD_PTR', 'Y'),
        'N6-acetyllysine': ('CCD_ALY', 'K'),
        '3-hydroxyproline': ('CCD_HYP', 'P'),
        'N6-methyllysine': ('CCD_MLY', 'K'),
        'N6,N6-dimethyllysine': ('CCD_M2L', 'K'),
        'N6,N6,N6-trimethyllysine': ('CCD_M3L', 'K'),
        'Omega-N-methylarginine': ('CCD_AGM', 'R'),
        'Phosphohistidine': ('CCD_HIP', 'H'),
        'Hydroxyproline': ('CCD_HYP', 'P'),
        'Pyrrolidone carboxylic acid': ('CCD_PCA', 'Q'),
        'Asymmetric dimethylarginine': ('CCD_ADM', 'R'),
        'Symmetric dimethylarginine': ('CCD_SDM', 'R'),
        '5-hydroxylysine': ('CCD_HYL', 'K'),
        'Hydroxylysine': ('CCD_HYL', 'K'),
        'N6-succinyllysine': ('CCD_SNN', 'K'),
        'S-palmitoyl cysteine': ('CCD_PTC', 'C'),
        'Phospholysine': ('CCD_P1L', 'K'),
    }

known_ions = ['MG', 'ZN', 'CL', 'CA', 'NA', 'MN', 'K', 'FE', 'CU', 'CO']
af3_ligands = ['ADP', 'ATP', 'AMP', 'GTP', 'GDP', 'FAD', 'NAD', 'NAP', 'NDP', 'HEM', 'HEC', 'PLM', 'OLA', 'MYR', 'CIT', 'CLA', 'CHL', 'BCL', 'BCB']
valid_dna_bases = {'A', 'T', 'C', 'G'}
valid_rna_bases = {'A', 'U', 'C', 'G'}

# Function to truncate long text
def truncate_text(text, max_length=50):
    return text if len(text) <= max_length else text[:max_length] + '...'

def search_uniprot(query):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&fields=accession,id,protein_name,gene_primary,organism_name,length&format=tsv"
    response = requests.get(url)
    lines = response.text.strip().split('\n')
    headers = lines[0].split('\t')
    results = [dict(zip(headers, line.split('\t'))) for line in lines[1:]]
    return results

def get_protein_data(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    response = requests.get(url)
    return response.json()

def get_regions(protein_data):
    regions = []
    for feature in protein_data.get('features', []):
        if feature['type'] in ['Chain', 'Signal peptide', 'Propeptide', 'Transit peptide', 'Initiator methionine', 'Peptide']:
            region = {}
            region['type'] = feature['type']
            region['description'] = feature.get('description', '')
            position_info = feature['location']
            if 'start' in position_info and 'end' in position_info:
                start_pos = int(position_info['start']['value'])
                end_pos = int(position_info['end']['value'])
                region['start'] = start_pos
                region['end'] = end_pos
                regions.append(region)
    return regions

def get_protein_features(protein_data, sequence_start, selected_sequence):
    glycans = []
    modifications = []
    glycan_positions = set()
    modification_positions = set()
    for feature in protein_data.get('features', []):
        position_info = feature['location']
        if 'start' in position_info and 'end' in position_info:
            if position_info['start']['value'] == position_info['end']['value']:
                position = int(position_info['start']['value'])
            else:
                continue  # Skip ranges for simplicity
        elif 'position' in position_info:
            position = int(position_info['position']['value'])
        else:
            continue

        if not (sequence_start <= position <= sequence_start + len(selected_sequence) - 1):
            continue  # Skip if position not in selected sequence

        adjusted_position = position - sequence_start + 1  # Adjust position for fragment
        amino_acid = selected_sequence[adjusted_position - 1]  # Adjust for 0-based indexing

        if feature['type'] == 'Glycosylation':
            description = feature.get('description', '')
            expected_residues = []

            if 'N-linked' in description:
                expected_residues.append('N')
            elif 'O-linked' in description:
                expected_residues.extend(['S', 'T'])
            else:
                expected_residues.append('N')

            if amino_acid not in expected_residues:
                st.warning(f"Expected residue {expected_residues} at position {adjusted_position}, but found '{amino_acid}'. Skipping glycosylation at this position.")
                continue

            glycan = {
                "residues": 'NAG(NAG(MAN(MAN(MAN)(MAN(MAN)(MAN)))))',
                "position": adjusted_position
            }
            if glycan not in glycans:
                glycans.append(glycan)
                glycan_positions.add(adjusted_position)
                    
        elif feature['type'] == 'Modified residue':
            description = feature.get('description', '')
            if description in ptm_mapping:
                ptm_type, expected_residue = ptm_mapping[description]
                if adjusted_position in glycan_positions or adjusted_position in modification_positions:
                    continue
                if amino_acid != expected_residue:
                    st.warning(f"Expected residue '{expected_residue}' at position {adjusted_position} for modification '{description}', but found '{amino_acid}'. Skipping modification at this position.")
                    continue
                modification = {
                    "ptmType": ptm_type,
                    "ptmPosition": adjusted_position
                }
                modifications.append(modification)
                modification_positions.add(adjusted_position)

    glycans.sort(key=lambda x: x['position'])
    modifications.sort(key=lambda x: x['ptmPosition'])

    return glycans, modifications

def get_ions_and_ligands(protein_data):
    ions = {}
    ligands = {}

    for feature in protein_data.get('features', []):
        if feature['type'] == 'Binding site':
            ligand_info = feature.get('ligand', {})
            ligand_name = ligand_info.get('name', '')
            if ligand_name:
                ligand_name = ligand_name.upper()
                ligand_symbol = ligand_name.split()[0].split('(')[0]

                if ligand_symbol in known_ions:
                    if ligand_symbol in ions:
                        ions[ligand_symbol] += 1
                    else:
                        ions[ligand_symbol] = 1
                else:
                    if ligand_symbol in af3_ligands:
                        if ligand_symbol in ligands:
                            ligands['CCD_'+ligand_symbol] += 1
                        else:
                            ligands['CCD_'+ligand_symbol] = 1
    return ions, ligands

def main():
    st.set_page_config(page_title='Uniprot to AlphaFold3 JSON Generator', page_icon=":comet:", layout="centered")
    
    st.title("Uniprot to AlphaFold3 JSON Generator :comet: (Beta)")

    model_seeds = []
    sequences = []
    ions_in_sequences = {}
    ligands_in_sequences = {}
    total_tokens = 0
    
    multi_protein = st.radio("Prediction of single protein or protein complex?", ('single', 'complex'), index=0)

    if multi_protein == 'single':
        num_proteins = 1
    else:
        num_proteins = st.number_input("Enter the number of different proteins:", min_value=1, value=2, step=1)

    for i in range(num_proteins):
        st.subheader(f"Protein {i+1}")

        # Query input
        query = st.text_input(f"Enter protein name, ID, or description for Protein {i+1}:")

        if query:
            results = search_uniprot(query)

            if not results:
                st.error("No results found.")
                continue

            # Simplify options for select box
            options = []
            option_map = {}
            for idx, entry in enumerate(results[:10]):
                entry_id = entry['Entry']
                gene_name = entry.get('Gene Names', 'N/A').split(' ')[0]
                organism = entry['Organism']
                protein_name = truncate_text(entry['Protein names'].split(' (')[0], max_length=40)
                option_text = f"{entry_id}: {protein_name} ({organism})"
                options.append(option_text)
                option_map[option_text] = idx

            selected_option = st.selectbox(f"Select a protein for Protein {i+1}:", options)
            selected_index = option_map[selected_option]
            selected_entry = results[selected_index]
            
            accession = selected_entry['Entry']
            
            # Display detailed information about the selected protein
            protein_details_markdown = f"""
            **Selected Protein Details:**

            - **Entry ID**: {selected_entry['Entry']}
            - **uniProtkbId**: {selected_entry['Entry Name']}
            - **Protein Name**: {selected_entry['Protein names']}
            - **Gene Names**: {selected_entry.get('Gene Names', 'N/A')}
            - **Organism**: {selected_entry['Organism']}
            - **Length**: {selected_entry['Length']}
            """

            st.markdown(protein_details_markdown)

            # Get protein data
            protein_data = get_protein_data(accession)
            sequence = protein_data['sequence']['value']
            length = len(sequence)

            # Whole sequence checkbox
            st.markdown("""---""")
            whole_sequence = st.checkbox("Use whole sequence?", value=True, key=f"whole_sequence_{i}")

            if whole_sequence:
                selected_sequence = sequence
                sequence_start = 1
                sequence_end = length
            else:
                selection_method = st.radio("Select sequence by:", ('Start and end positions', 'UniProt-defined regions'), key=f"selection_method_{i}")

                if selection_method == 'Start and end positions':
                    sequence_start = st.number_input("Start position:", min_value=1, max_value=length, value=1)
                    sequence_end = st.number_input("End position:", min_value=1, max_value=length, value=length)
                    selected_sequence = sequence[sequence_start - 1:sequence_end]
                else:
                    regions = get_regions(protein_data)
                    if not regions:
                        st.warning("No defined regions available. Please enter start and end positions manually.")
                        sequence_start = st.number_input("Start position:", min_value=1, max_value=length, value=1)
                        sequence_end = st.number_input("End position:", min_value=1, max_value=length, value=length)
                        selected_sequence = sequence[sequence_start - 1:sequence_end]
                    else:
                        region_options = [f"{region['type']} ({truncate_text(region['description'], 30)}): {region['start']}-{region['end']}" for region in regions]
                        selected_region_option = st.selectbox("Select a region:", region_options, key=f"selected_region_{i}")
                        selected_region_index = region_options.index(selected_region_option)
                        selected_region = regions[selected_region_index]
                        sequence_start = selected_region['start']
                        sequence_end = selected_region['end']
                        selected_sequence = sequence[sequence_start - 1:sequence_end]

            
            glycans, modifications = get_protein_features(protein_data, sequence_start, selected_sequence)

            ions, ligands = get_ions_and_ligands(protein_data)

            for ion_name, count in ions.items():
                if ion_name in ions_in_sequences:
                    ions_in_sequences[ion_name] += count
                else:
                    ions_in_sequences[ion_name] = count

            count = st.number_input(f"Number of copies for Protein {i+1}:", min_value=1, value=1, step=1)
            total_tokens += len(selected_sequence)*count

            protein_chain = {
                "sequence": selected_sequence,
                "count": count
            }
            c1,c2 = st.columns(2)
            if c1.checkbox("Add glycans to the protein chain?", key=f"glycans_{i}", value=True):
                if glycans:
                    # add here glycan modification option if needed
                    
                    protein_chain["glycans"] = glycans
                    total_tokens += sum([85 for glycan in glycans])*count # rough estimate of glycan tokens
                    
            if c2.checkbox("Add modifications to the protein chain?", key=f"modifications_{i}", value=True):
                if modifications:
                    protein_chain["modifications"] = modifications
                    total_tokens += sum([20 for modification in modifications])*count # rough estimate of modification tokens

            sequences.append({"proteinChain": protein_chain})

            for ion_name, count in ions.items():
                if ion_name in ions_in_sequences:
                    ions_in_sequences[ion_name] += count
                else:
                    ions_in_sequences[ion_name] = count
            for ligand_name, count in ligands.items():
                if ligand_name in ligands_in_sequences:
                    ligands_in_sequences[ligand_name] += count
                else:
                    ligands_in_sequences[ligand_name] = count
    
    st.markdown("""---""")
    
    c1,c2 = st.columns(2)
    
    if c1.checkbox("Add ions and ligands to the sequences?", value=True):
        for ion_name, count in ions_in_sequences.items():
            sequences.append({"ion": {"ion": ion_name, "count": count}})
            total_tokens += count

    if c2.checkbox("Add ligands to the sequences?", value=True):
        for ligand_name, count in ligands_in_sequences.items():
            sequences.append({"ligand": {"ligand": ligand_name, "count": count}})
            total_tokens += count*100 # rough estimate of ligand tokens

    # ** New Section: Add Custom Sequences, Ions, Ligands, or Nucleic Acid (DNA, RNA) **
    st.markdown("## Add Custom Input")
    
    c1, c2 = st.columns(2)
    if c1.checkbox("Add custom input?", key="add_custom_input"):
        number_of_own_sequences = c2.number_input("Enter the number of custom sequences, ions, ligands, or nucleic acids to add:", min_value=1, value=1, step=1, key="num_custom_sequences")
        
        for n in range(number_of_own_sequences):
            custom_type = st.selectbox(f"What would you like to add? (Custom Input {n+1})", 
                                    ["Protein Sequence", "DNA Sequence", "RNA Sequence", "Ion", "Ligand"], 
                                    key=f"custom_type_{n}")

            if custom_type == "Protein Sequence":
                custom_sequence = st.text_area(f"Enter protein sequence for Custom Input {n+1}:", key=f"protein_sequence_{n}")
                custom_count = st.number_input(f"Number of copies for Custom Input {n+1}:", min_value=1, value=1, step=1, key=f"protein_count_{n}")
                if custom_sequence:
                    sequences.append({"proteinChain": {"sequence": custom_sequence, "count": custom_count}})
                    total_tokens += len(custom_sequence) * custom_count

            elif custom_type == "DNA Sequence":
                custom_sequence = st.text_area(f"Enter DNA sequence for Custom Input {n+1}:", key=f"dna_sequence_{n}")
                custom_count = st.number_input(f"Number of copies for Custom Input {n+1}:", min_value=1, value=1, step=1, key=f"dna_count_{n}")
                
                # Validate DNA sequence
                if custom_sequence:
                    if set(custom_sequence.upper()).issubset(valid_dna_bases):
                        sequences.append({"dnaSequence": {"sequence": custom_sequence, "count": custom_count}})
                        total_tokens += len(custom_sequence) * custom_count
                    else:
                        st.error(f"Invalid DNA sequence: '{custom_sequence}'. DNA can only contain the bases: A, T, C, G.")

            elif custom_type == "RNA Sequence":
                custom_sequence = st.text_area(f"Enter RNA sequence for Custom Input {n+1}:", key=f"rna_sequence_{n}")
                custom_count = st.number_input(f"Number of copies for Custom Input {n+1}:", min_value=1, value=1, step=1, key=f"rna_count_{n}")
                
                # Validate RNA sequence
                if custom_sequence:
                    if set(custom_sequence.upper()).issubset(valid_rna_bases):
                        sequences.append({"rnaSequence": {"sequence": custom_sequence, "count": custom_count}})
                        total_tokens += len(custom_sequence) * custom_count
                    else:
                        st.error(f"Invalid RNA sequence: '{custom_sequence}'. RNA can only contain the bases: A, U, C, G.")

            elif custom_type == "Ion":
                custom_ion = st.text_input(f"Enter ion symbol (e.g., CA, MG) for Custom Input {n+1}:", key=f"ion_symbol_{n}").upper()
                custom_count = st.number_input(f"Number of ions for Custom Input {n+1}:", min_value=1, value=1, step=1, key=f"ion_count_{n}")
                if custom_ion:
                    if custom_ion in known_ions:
                        sequences.append({"ion": {"ion": custom_ion, "count": custom_count}})
                        total_tokens += custom_count
                    else:
                        st.error(f"Unknown ion symbol: {custom_ion}. Please use one of the following: {', '.join(known_ions)}")

            elif custom_type == "Ligand":
                custom_ligand = st.text_input(f"Enter ligand symbol (e.g., ATP, HEM) for Custom Input {n+1}:", key=f"ligand_symbol_{n}").upper()
                custom_count = st.number_input(f"Number of ligands for Custom Input {n+1}:", min_value=1, value=1, step=1, key=f"ligand_count_{n}")
                if custom_ligand:
                    if custom_ligand in af3_ligands:
                        sequences.append({"ligand": {"ligand": "CCD_" + custom_ligand, "count": custom_count}})
                        total_tokens += custom_count
                    else:
                        st.error(f"Unknown ligand symbol: {custom_ligand}. Please use one of the following: {', '.join(af3_ligands)}")

    
    if sequences:
        # Generate output
        timestamp = datetime.datetime.now().strftime('%y%m%d_%H%M%S')
        try:
            job_name = f"{selected_entry['Entry Name']}_{timestamp}"
        except:
            job_name = f"custom_job_{timestamp}"

        output = [
            {
                "name": job_name,
                "modelSeeds": model_seeds,
                "sequences": sequences
            }
        ]

        if total_tokens > 5000:
            st.error(f"Total tokens ({total_tokens}) exceed the maximum limit of 5000. Please reduce the number of copies or remove some chains.")
        else:
            st.progress(total_tokens/5000, f"Total tokens: {total_tokens} (max 5000)")

        st.subheader("Generated JSON:")
        st.json(output)

        filename_job_name = job_name.replace(':', '').replace(' ', '_')
        filename = f"{filename_job_name}_job_request.json"
        
        st.download_button(
            label="Download JSON",
            data=json.dumps(output, indent=4),
            file_name=filename,
            mime="application/json"
        )

        st.markdown("""
        ---
        **[Open the AlphaFold3](https://golgi.sandbox.google.com/) website and upload the generated JSON file to run the prediction.**"""
        )

    st.markdown("""
    ---
    <div style="text-align: center;">
        <strong><a href="https://github.com/dzyla/AF3_JSON_generator">Source code</a></strong> | Developed by <a href="https://dzyla.com">Dawid Zyla</a>
    </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
