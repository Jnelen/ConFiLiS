#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ConFiLiS Library Creator

This script reads an input file (SDF or CSV/TSV) containing compounds,
computes various fingerprints in parallel, and saves the resulting fingerprint dictionaries
as compressed pickle files.

Usage:
    python library_creator.py input_file output_directory [options]

Options:
    --delimiter       Delimiter for CSV/TSV files (auto-detected if not specified)
    --smiles_index    Column index for SMILES (default: 0)
    --name_index      Column index for compound name (default: None, auto-determined)
    --header          Skip header row in CSV/TSV files
    --num_workers     Number of parallel workers (default: 1)
    --fingerprints    List of fingerprint types to compute.
                      Supported fingerprints: avalon, atom-pair, cdk-substructure, circular, ECFP4, ECFP6, klekota-roth, mol2vec, pubchem, rdk-maccs, rdkit, shortestpath, topological-torsion
    --silent          Disable progress bars

Author: Jochem Nelen (jnelen@ucam.edu)
"""

import argparse
import csv
import gzip
import logging
import os
import pickle
import time
from typing import Dict, List, Optional, Tuple

import multiprocessing as mp
import concurrent.futures

import numpy as np

from PyFingerprint.fingerprint import get_fingerprint
from rdkit import Chem
from tqdm import tqdm


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="ConFiLiS Fingerprint Library Creator")
    parser.add_argument("input_file", help="Input file (SDF or CSV/TSV)")
    parser.add_argument("output_path", help="Output directory")
    parser.add_argument(
        "-f",
        "--fingerprints",
        "--fps",
        nargs="+",
        required=True,
        type=str,
        help="List of fingerprint types to compute. Supported fingerprints: avalon, atom-pair, cdk-substructure, circular, ECFP4, ECFP6, klekota-roth, mol2vec, pubchem, rdk-maccs, rdkit, shortestpath, topological-torsion",
    )
    parser.add_argument(
        "--delimiter",
        default=None,
        help="Delimiter for CSV/TSV files (auto-detected if not specified)",
    )
    parser.add_argument(
        "--smiles_index",
        type=int,
        default=0,
        help="Column index for SMILES (default: 0)",
    )
    parser.add_argument(
        "--name_index",
        type=int,
        default=None,
        help="Column index for compound name (default: None, auto-determined)",
    )
    parser.add_argument(
        "--header",
        action="store_true",
        default=False,
        help="Skip header row in CSV/TSV files",
    )
    parser.add_argument(
        "--silent", action="store_true", default=False, help="Disable progress bars"
    )
    parser.add_argument(
        "--num_workers",
        "--cores",
        "-c",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1)",
    )

    return parser.parse_args()


def detect_delimiter(file_path: str) -> str:
    """
    Detect the delimiter of a CSV/TSV file by reading a sample.

    Args:
        file_path (str): Path to the input file.

    Returns:
        str: Detected delimiter character.
    """
    with open(file_path, "r", newline="") as f:
        sample = f.read(2048)
        sniffer = csv.Sniffer()
        detected_delimiter = sniffer.sniff(sample).delimiter
    return detected_delimiter


def read_library(
    input_file: str,
    delimiter: Optional[str] = None,
    smiles_index: int = 0,
    name_index: Optional[int] = None,
    header: bool = False,
) -> Dict[str, str]:
    """
    Read the input library file and return a dictionary mapping compound names to SMILES strings.

    Args:
        input_file (str): Path to the input file (SDF or CSV/TSV).
        delimiter (Optional[str]): Delimiter for CSV/TSV files. Auto-detected if None.
        smiles_index (int): Column index for SMILES.
        name_index (Optional[int]): Column index for compound name. If None, uses SMILES as name.
        header (bool): If True, skip the header row in CSV/TSV files.

    Returns:
        Dict[str, str]: Mapping from compound name to SMILES string.
    """
    smiles_dict: Dict[str, str] = {}

    if input_file.lower().endswith(".sdf"):
        # Read SDF file using RDKit
        suppl = Chem.SDMolSupplier(input_file)
        for mol in suppl:
            if mol is None:
                continue
            mol_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
            mol_name = (
                mol.GetProp("_Name").strip()
                if mol.HasProp("_Name") and mol.GetProp("_Name").strip()
                else mol_smiles
            )
            smiles_dict[mol_name] = mol_smiles
    else:
        # Auto-detect delimiter if not provided
        if delimiter is None:
            delimiter = detect_delimiter(input_file)
            logging.info(f"Auto-detected delimiter: '{delimiter}'")

        with open(input_file, newline="") as f:
            reader = csv.reader(f, delimiter=delimiter)
            if header:
                next(reader, None)  # Skip header row
            for row_idx, row in enumerate(reader, start=1):
                try:
                    smiles = row[smiles_index].strip()
                    compound_name = (
                        row[name_index].strip()
                        if name_index is not None and len(row) > name_index
                        else smiles
                    )
                    smiles_dict[compound_name] = smiles
                except IndexError:
                    logging.error(f"Error processing row {row_idx}: {row}")
                    continue

    if not smiles_dict:
        logging.warning("No compounds were detected in the input file.")
    return smiles_dict


def compute_fingerprint_for_type(
    args_tuple: Tuple[str, Dict[str, str], str, str, bool],
) -> Tuple[str, Dict[str, np.ndarray]]:
    """
    Compute fingerprints for all compounds using the specified fingerprint type.

    Args:
        args_tuple: A tuple containing:
            - fingerprint (str): Fingerprint type.
            - smiles_dict (Dict[str, str]): Mapping of compound names to SMILES.
            - output_path (str): Output directory (unused here).
            - base_name (str): Base name for output files (unused here).
            - silent (bool): Whether to suppress progress bars.

    Returns:
        Tuple[str, Dict[str, np.ndarray]]: The fingerprint type and a mapping from compound names to fingerprint arrays.
    """
    fingerprint, smiles_dict, output_path, base_name, silent = args_tuple
    fp_dict: Dict[str, np.ndarray] = {}
    fingerprint_lower = fingerprint.lower()

    if fingerprint_lower == "ecfp4":
        fp_func = lambda s: get_fingerprint(s, "morgan", depth=2).to_numpy()
    elif fingerprint_lower == "ecfp6":
        fp_func = lambda s: get_fingerprint(s, "morgan", depth=3).to_numpy()
    else:
        fp_func = lambda s: get_fingerprint(s, fingerprint).to_numpy()

    iterator = (
        tqdm(
            smiles_dict.items(),
            desc=f"Processing {fingerprint}",
            leave=True,
            dynamic_ncols=True,
            disable=silent,
        )
        if not silent
        else smiles_dict.items()
    )

    for compound_name, mol_smiles in iterator:
        try:
            fp_dict[compound_name] = fp_func(mol_smiles)
        except Exception as e:
            msg = f"Error computing {fingerprint} for {compound_name}: {e}"
            if not silent:
                tqdm.write(msg)
            else:
                logging.error(msg)

    return fingerprint, fp_dict


def main() -> None:
    """
    Main function to parse arguments, read the compound library,
    compute fingerprints in parallel, and save the fingerprint dictionaries as compressed pickle files.
    """
    start_time = time.time()

    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    args = parse_arguments()

    # Auto-determine name_index if not provided.
    if args.name_index is None:
        args.name_index = 1 if args.smiles_index == 0 else 0
    logging.info(
        f"Using column {args.smiles_index} for SMILES and column {args.name_index} for names."
    )

    # Read the input library.
    smiles_dict = read_library(
        args.input_file, args.delimiter, args.smiles_index, args.name_index, args.header
    )
    num_compounds = len(smiles_dict)

    base_name = os.path.splitext(os.path.basename(args.input_file))[0]

    # Ensure output directory exists.
    os.makedirs(args.output_path, exist_ok=True)

    logging.info(f"Detected {num_compounds} compounds.")

    smiles_dict_path = os.path.join(args.output_path, f"{base_name}_smiles.pkl.gz")
    with gzip.open(smiles_dict_path, "wb") as f:
        pickle.dump(smiles_dict, f)
        logging.info(f"Saved the smiles dictionary to {smiles_dict_path}")

    logging.info(
        f"Computing {len(args.fingerprints)} fingerprints: {args.fingerprints}"
    )

    # Prepare tasks: each task is a tuple for a fingerprint type. If there are multiple workers always set silent to True.
    tasks: List[Tuple[str, Dict[str, str], str, str, bool]] = [
        (
            fp,
            smiles_dict,
            args.output_path,
            base_name,
            args.silent or args.num_workers > 1,
        )
        for fp in args.fingerprints
    ]

    if args.num_workers > 1:
        mp.set_start_method("spawn", force=True)
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=min(args.num_workers, len(args.fingerprints))
        ) as executor:
            futures = {
                executor.submit(compute_fingerprint_for_type, task): task[0]
                for task in tasks
            }
            with tqdm(
                total=len(futures),
                desc="Overall fingerprint progress",
                disable=args.silent,
            ) as pbar:
                for future in concurrent.futures.as_completed(futures):
                    fingerprint = futures[future]
                    try:
                        fp_type, fp_dict = future.result()
                        output_file = os.path.join(
                            args.output_path, f"{base_name}_{fp_type}.pkl.gz"
                        )
                        with gzip.open(output_file, "wb") as f:
                            pickle.dump(fp_dict, f)
                        logging.info(
                            f"Saved {len(fp_dict)} {fp_type} fingerprints to {output_file}"
                        )
                    except Exception as exc:
                        logging.error(
                            f"Fingerprint {fingerprint} generated an exception: {exc}"
                        )
                    pbar.update(1)
    else:
        # Sequential processing: compute and save each fingerprint immediately.
        for task in tasks:
            fp_type, fp_dict = compute_fingerprint_for_type(task)
            output_file = os.path.join(
                args.output_path, f"{base_name}_{fp_type}.pkl.gz"
            )
            try:
                with gzip.open(output_file, "wb") as f:
                    pickle.dump(fp_dict, f)
                logging.info(
                    f"Saved {len(fp_dict)} {fp_type} fingerprints to {output_file}"
                )
            except Exception as e:
                logging.error(
                    f"Failed to save {fp_type} fingerprints to {output_file}: {e}"
                )

    total_time = time.time() - start_time
    logging.info(
        f"Finished processing {num_compounds} compounds in {total_time:.2f} seconds."
    )


if __name__ == "__main__":
    main()
