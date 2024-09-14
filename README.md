# UniProt to AlphaFold3 JSON Generator

This application allows users to generate JSON files for input into the AlphaFold3 protein structure prediction platform. It extracts protein sequence data from UniProt and formats it to be compatible with AlphaFold3, supporting features like glycosylation, post-translational modifications (PTMs), and ion/ligand binding information.

## Features
- Query protein data from UniProt by name, ID, or description.
- Select entire or specific regions of a protein sequence.
- Add glycosylation and post-translational modifications.
- Include known ions and ligands in the sequence.
- Generate custom protein, DNA, RNA sequences, ions, and ligands.
- Download the generated JSON for use in AlphaFold3.

## Usage
1. Install required dependencies:
    ```bash
    pip install streamlit requests
    ```
2. Run the application locally:
    ```bash
    streamlit run app.py
    ```
3. Use the interface to select proteins, customize sequences, and generate the JSON file.

## Example

1. Query a protein from UniProt.
2. Select regions, modifications, or ligands.
3. Generate and download the JSON file.
4. Upload the JSON to [AlphaFold3](https://golgi.sandbox.google.com/) for structure prediction.

## Development
Suggestions and bug reports are always welcome! use the Issues to report bug, functionality request or just say hi.
