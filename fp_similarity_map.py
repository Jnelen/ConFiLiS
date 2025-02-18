#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Similarity Map Generator using RDKit

This script computes similarity maps for a given query molecule against a reference molecule
from a provided ConFiLiS molecular library.

More information regarding Similarity Maps: https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-5-43

Usage:
    python fp_similarity_map.py --query <query_smiles> --id <compound_id> --library <library_path> [options]

Arguments:
    --query         SMILES string of the query molecule.
    --id            Compound ID to find in the library (should be between quotes).
    --library       Path to the library directory containing *_smiles.pkl.gz files.
    --output        Output directory to save similarity maps (default: output/SimilarityMap).
    --types         Fingerprint types for similarity mapping (default: Morgan_1).
                    Available types: morgan_1, morgan_2, morgan_3, ap, tt, or 'all' for all.

Author:
    Jochem Nelen (jnelen@ucam.edu)
"""

import argparse
import gzip
import glob
import pickle
import os
import re
import sys
from typing import Dict, Set

from rdkit import Chem
from rdkit.Chem.Draw import SimilarityMaps


def sanitize_filename(smiles: str) -> str:
    """
    Converts a SMILES string into a safe filename format.

    Args:
        smiles (str): Input SMILES string.

    Returns:
        str: Sanitized filename-safe SMILES string.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        smiles = Chem.MolToSmiles(mol, canonical=True)  # Canonicalize SMILES
        return re.sub(r"[^\w\-_.]", "_", smiles)  # Replace invalid filename characters
    except Exception as e:
        print(f"Error processing SMILES: {e}")
        sys.exit(1)


def load_library(library_path: str) -> Dict[str, str]:
    """
    Loads a molecular library from a compressed pickle file.

    Args:
        library_path (str): Path to the library directory containing *_smiles.pkl.gz files.

    Returns:
        Dict[str, str]: Dictionary mapping compound IDs to SMILES strings.

    Raises:
        FileNotFoundError: If no valid pickle file is found.
    """
    try:
        smiles_dict_files = glob.glob(f"{library_path}/*_smiles.pkl.gz")
        if not smiles_dict_files:
            raise FileNotFoundError(f"No *_smiles.pkl.gz file found in {library_path}")

        smiles_dict_path = smiles_dict_files[0]  # Use the first matching file
        with gzip.open(smiles_dict_path, "rb") as zip_dict:
            lookup_dict = pickle.load(zip_dict)

        if not isinstance(lookup_dict, dict):
            raise ValueError("Loaded library file does not contain a valid dictionary.")

        return lookup_dict

    except Exception as e:
        print(f"Error loading library: {e}")
        sys.exit(1)


def generate_similarity_map(
    query_mol: Chem.Mol, ref_mol: Chem.Mol, method, output_path: str
) -> None:
    """
    Generates a similarity map image for two molecules.

    Args:
        query_mol (Chem.Mol): Query molecule.
        ref_mol (Chem.Mol): Reference molecule.
        method (function): Fingerprint method for similarity mapping.
        output_path (str): Path to save the generated image.

    Returns:
        None
    """
    try:
        fig, _ = SimilarityMaps.GetSimilarityMapForFingerprint(
            query_mol, ref_mol, method
        )
        fig.savefig(output_path, bbox_inches="tight")
        print(f"Saved: {output_path}")
    except Exception as e:
        print(f"Error generating similarity map ({output_path}): {e}")


def main():
    """
    Main function that processes input arguments, loads the molecular library,
    generates similarity maps, and saves them to the specified directory.

    Args:
        None (arguments are parsed from command line)

    Returns:
        None
    """
    parser = argparse.ArgumentParser(
        description="Generate similarity maps using RDKit."
    )
    parser.add_argument(
        "--query", required=True, help="SMILES string of the query molecule (should be between quotes)"
    )
    parser.add_argument(
        "--id", required=True, help="Compound ID to find in the library"
    )
    parser.add_argument(
        "--library", required=True, help="Path to the library directory"
    )
    parser.add_argument(
        "--output",
        default="output/SimilarityMap",
        help="Directory to save output images",
    )
    parser.add_argument(
        "--types",
        nargs="+",
        default=["Morgan_1"],
        help="Fingerprint types to generate similarity maps for (case insensitive). Use 'all' to generate all types.",
    )

    args = parser.parse_args()

    # Define valid fingerprint types
    valid_types: Set[str] = {"morgan_1", "morgan_2", "morgan_3", "ap", "tt"}

    # Normalize fingerprint types (case-insensitive)
    selected_types = {t.lower() for t in args.types}

    # If "all" is selected, use all valid types
    if "all" in selected_types:
        selected_types = valid_types
    else:
        selected_types &= valid_types  # Keep only valid types

    if not selected_types:
        print(
            f"Error: No valid fingerprint types selected. Choose from {valid_types} or use 'all'."
        )
        sys.exit(1)

    # Load the compound library
    lookup_dict = load_library(args.library)

    # Validate compound ID
    if args.id not in lookup_dict:
        print(f"Error: Compound ID '{args.id}' not found in the library.")
        sys.exit(1)

    # Convert SMILES to RDKit molecules
    query_mol = Chem.MolFromSmiles(args.query)
    ref_mol = Chem.MolFromSmiles(lookup_dict[args.id])

    if query_mol is None or ref_mol is None:
        print("Error: Invalid SMILES input or reference molecule.")
        sys.exit(1)

    # Ensure output directory exists
    os.makedirs(args.output, exist_ok=True)

    # Sanitize SMILES for filename usage
    smiles_filename = sanitize_filename(args.query)

    # Generate similarity maps
    if "morgan_1" in selected_types:
        generate_similarity_map(
            query_mol,
            ref_mol,
            lambda m, idx: SimilarityMaps.GetMorganFingerprint(
                m, atomId=idx, radius=1, fpType="count"
            ),
            os.path.join(args.output, f"Morgan_1_{args.id}_{smiles_filename}.png"),
        )

    if "morgan_2" in selected_types:
        generate_similarity_map(
            query_mol,
            ref_mol,
            lambda m, idx: SimilarityMaps.GetMorganFingerprint(
                m, atomId=idx, radius=2, fpType="count"
            ),
            os.path.join(args.output, f"Morgan_2_{args.id}_{smiles_filename}.png"),
        )

    if "morgan_3" in selected_types:
        generate_similarity_map(
            query_mol,
            ref_mol,
            lambda m, idx: SimilarityMaps.GetMorganFingerprint(
                m, atomId=idx, radius=3, fpType="count"
            ),
            os.path.join(args.output, f"Morgan_3_{args.id}_{smiles_filename}.png"),
        )

    if "ap" in selected_types:
        generate_similarity_map(
            query_mol,
            ref_mol,
            SimilarityMaps.GetAPFingerprint,
            os.path.join(args.output, f"AP_{args.id}_{smiles_filename}.png"),
        )

    if "tt" in selected_types:
        generate_similarity_map(
            query_mol,
            ref_mol,
            SimilarityMaps.GetTTFingerprint,
            os.path.join(args.output, f"TT_{args.id}_{smiles_filename}.png"),
        )


if __name__ == "__main__":
    main()
