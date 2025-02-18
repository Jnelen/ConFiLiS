#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ConFiLiS: Consensus Fingerprints for Ligand-based Screening

This command-line script compares a query compound (provided as a SMILES string)
against compounds in one or more fingerprint libraries. For each provided library directory,
the script automatically identifies fingerprint files using glob, loads the appropriate libraries,
computes Euclidean distances between the query fingerprint and each library compound, and aggregates
the results with Z-score normalization.

Usage:
    python screen.py -i "<SMILES>" -l <library_dir1> [<library_dir2> ...] -f <fingerprint1> [<fingerprint2> ...] [--cores <num>]

Arguments:
    -i, --input, --smiles
         The SMILES string for the query compound.

    -n, --molname, --name
         Optional compound name. If not provided, the SMILES string is used as the name.

    -l, --libraries, --libs
         One or more library directories to search. Each directory is expected to contain fingerprint files
         (e.g. lib1_circular.pkl.gz) where the fingerprint type is identified by splitting the filename on '_' and
         extracting the last part (before the extension). The keys are stored in lowercase.

    -f, --fingerprints, --fps
         Some examples of supported fingerprints:
         avalon, atom-pair, cdk-substructure, circular, ECFP4, ECFP6, klekota-roth, mol2vec, pubchem, rdk-maccs, rdkit, shortestpath, topological-torsion

    -o, --output, --output_path
         Output directory or file path prefix. Defaults to 'output' if not specified.

    -c, --cores, --num_workers
         Number of worker processes to use (default: 1).

    --shortened_rows
         Maximum number of rows for the shortened output file (default: 1000).

Author: Jochem Nelen (jnelen@ucam.edu)
"""

import argparse
import glob
import gzip
import logging
import os
import pickle
import sys
import time
from typing import Dict, List

import numpy as np
import pandas as pd
from PyFingerprint.fingerprint import get_fingerprint

import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed

from tqdm import tqdm


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="ConFiLiS: Consensus Fingerprints for Ligand-based Screening",
        epilog="""Usage example:
    python screen.py -i "C1=CC=CC=C1" -l ./lib1 ./lib2 -f circular rdk-maccs --cores 4
    This compares the benzene SMILES against the libraries located in ./lib1 and ./lib2 using circular and rdk-maccs fingerprints.
    """,
    )
    parser.add_argument(
        "-i",
        "--input",
        "--smiles",
        required=True,
        type=str,
        help="Input SMILES string for the query compound.",
    )
    parser.add_argument(
        "-n",
        "--molname",
        "--name",
        required=False,
        default=None,
        type=str,
        help="Optional compound name. If not provided, the SMILES string is used as the name.",
    )
    parser.add_argument(
        "-l",
        "--libraries",
        "--libs",
        nargs="+",
        required=True,
        type=str,
        help="List of library directories to compare against. Each directory should contain fingerprint files (e.g. test_circular.pkl.gz).",
    )
    parser.add_argument(
        "-f",
        "--fingerprints",
        "--fps",
        nargs="+",
        required=True,
        type=str,
        help="List of fingerprint types to use. Supported fingerprints: avalon, atom-pair, cdk-substructure, circular, ECFP4, ECFP6, klekota-roth, mol2vec, pubchem, rdk-maccs, rdkit, shortestpath, topological-torsion",
    )
    parser.add_argument(
        "-o",
        "--output",
        "--output_path",
        required=False,
        type=str,
        help="Output directory or file path prefix. Defaults to 'output' if not specified.",
    )
    parser.add_argument(
        "-c",
        "--cores",
        "--num_workers",
        type=int,
        default=1,
        help="Number of worker processes to use (default: 1).",
    )
    parser.add_argument(
        "--shortened_rows",
        type=int,
        default=1000,
        help="Maximum number of rows for the shortened output file (default: 1000).",
    )

    return parser.parse_args()


def compute_fingerprint(smiles: str, fp_type: str) -> np.ndarray:
    """Compute the fingerprint for a given SMILES string and fingerprint type."""
    try:
        if fp_type.lower() == "ecfp6":
            return get_fingerprint(smiles, "morgan", depth=3).to_numpy()
        elif fp_type.lower() == "ecfp4":
            return get_fingerprint(smiles, "morgan", depth=2).to_numpy()
        else:
            return get_fingerprint(smiles, fp_type).to_numpy()
    except Exception as e:
        raise ValueError(
            f"Failed to compute {fp_type} fingerprint for SMILES '{smiles}': {e}"
        )


def load_library_file(library_file: str) -> Dict[str, np.ndarray]:
    """Load the fingerprint library from a gzip-compressed pickle file."""
    if not os.path.isfile(library_file):
        raise FileNotFoundError(f"Library file not found: {library_file}")
    try:
        with gzip.open(library_file, "rb") as f:
            library_data = pickle.load(f)
        return library_data
    except Exception as e:
        raise Exception(f"Error loading library file {library_file}: {e}")


def process_fingerprint(
    query_fp: np.ndarray, library_data: Dict[str, np.ndarray]
) -> pd.DataFrame:
    """Compute Euclidean distances and perform z-score normalization."""
    if not library_data:
        return pd.DataFrame(columns=["raw_distance", "normalized_distance"])
    compounds = list(library_data.keys())
    fps = np.array([library_data[c] for c in compounds])
    diff = fps - query_fp
    dists = np.linalg.norm(diff, axis=1)
    mean_val = np.mean(dists)
    std_val = np.std(dists)
    normed = np.zeros_like(dists) if std_val == 0 else (dists - mean_val) / std_val
    result_df = pd.DataFrame(
        {"raw_distance": dists, "normalized_distance": normed}, index=compounds
    )
    return result_df


def aggregate_results(
    all_results: Dict[str, pd.DataFrame], fp_types: List[str]
) -> pd.DataFrame:
    """Aggregate individual fingerprint results into one DataFrame with averages."""
    df_list: List[pd.DataFrame] = []
    raw_columns: List[str] = []
    norm_columns: List[str] = []
    for fp in fp_types:
        df = all_results.get(fp)
        if df is None or df.empty:
            continue
        df = df.rename(
            columns={"raw_distance": fp, "normalized_distance": f"{fp}_normalized"}
        )
        df_list.append(df)
        raw_columns.append(fp)
        norm_columns.append(f"{fp}_normalized")
    aggregated_df = pd.concat(df_list, axis=1)
    aggregated_df["Total_Average"] = aggregated_df[raw_columns].mean(
        axis=1, skipna=True
    )
    aggregated_df["Consensus_Average"] = aggregated_df[norm_columns].mean(
        axis=1, skipna=True
    )
    aggregated_df = aggregated_df.sort_values("Consensus_Average", na_position="last")
    aggregated_df = aggregated_df.reset_index().rename(columns={"index": "Name"})
    return aggregated_df


def write_results_to_csv(
    df: pd.DataFrame, fp_types: List[str], output_prefix: str, shortened_rows: int
) -> None:
    """Write aggregated results to CSV files."""
    output_file = f"{output_prefix}.csv"
    try:
        df.to_csv(output_file, index=False, float_format="%.3f")
        logging.info(f"Results written to {output_file}")
    except Exception as e:
        logging.error(f"Error writing to {output_file}: {e}")
    if shortened_rows > 0 and len(df) > shortened_rows:
        short_file = f"{output_prefix}_short.csv"
        try:
            df.head(shortened_rows).to_csv(short_file, index=False, float_format="%.3f")
            logging.info(f"Shortened results written to {short_file}")
        except Exception as e:
            logging.error(f"Error writing to {short_file}: {e}")


def process_single_fingerprint(
    library_dir: str, compound_name: str, query_smiles: str, fp_type: str
):
    """
    Process a single fingerprint type for a given library directory.
    Returns a tuple (library_name, fp_type, result DataFrame) or None.
    """
    library_name = os.path.basename(os.path.normpath(library_dir))
    lib_file = None
    pattern = os.path.join(library_dir, "*.pkl.gz")
    for file_path in glob.glob(pattern):
        fname = os.path.basename(file_path)
        parts = fname.split("_")
        if len(parts) < 2:
            continue
        file_fp = parts[-1].split(".")[0].lower()
        if file_fp == fp_type.lower():
            lib_file = file_path
            break

    if lib_file is None:
        logging.warning(
            f"Fingerprint '{fp_type}' not found in library directory '{library_dir}'."
        )
        return None

    try:
        library_data = load_library_file(lib_file)
    except Exception as e:
        logging.error(f"Error loading library file {lib_file}: {e}")
        return None

    try:
        query_fp = compute_fingerprint(query_smiles, fp_type)
    except ValueError as e:
        logging.error(e)
        return None

    logging.info(
        f"Processing {library_name} for fingerprint '{fp_type}' using file: {lib_file}"
    )
    fp_results = process_fingerprint(query_fp, library_data)
    if fp_results.empty:
        logging.warning(
            f"No valid distances computed for fingerprint '{fp_type}' in library '{library_name}'."
        )
        return None

    return (library_name, fp_type, fp_results)


def process_library_results(
    lib_name: str,
    results: Dict[str, pd.DataFrame],
    output_dir: str,
    compound_name: str,
    shortened_rows: int,
) -> None:
    """Aggregate results and write CSV output for a given library."""
    if not results:
        logging.error(
            f"No valid results obtained for library '{lib_name}'. Skipping calculations."
        )
        return
    aggregated_results = aggregate_results(results, list(results.keys()))
    if os.path.isdir(output_dir):
        output_prefix = os.path.join(output_dir, f"ConFiLiS_{compound_name}_{lib_name}")
    else:
        base, _ = os.path.splitext(output_dir)
        output_prefix = f"{base}_{lib_name}"
    write_results_to_csv(
        aggregated_results, list(results.keys()), output_prefix, shortened_rows
    )
    logging.info(
        f"Finished processing library '{lib_name}'. Total compounds processed: {len(aggregated_results)}\n"
    )


def main() -> None:
    start_time = time.time()
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    args = parse_arguments()
    query_smiles = args.input
    compound_name = args.molname if args.molname else query_smiles

    # Sort fingerprint types.
    fp_types = sorted([fp.lower() for fp in args.fingerprints])

    logging.info(f"Processing compound '{compound_name}' with fingerprints: {fp_types}")
    logging.info(f"Comparing against library directories: {args.libraries}\n")

    # Validate library directories.
    for lib in args.libraries:
        if not os.path.isdir(lib):
            logging.error(f"Library directory does not exist: {lib}")
            sys.exit(1)

    # Set up output directory.
    output_dir = args.output if args.output else "output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    total_tasks = len(args.libraries) * len(fp_types)
    num_workers = args.cores
    if num_workers > total_tasks:
        logging.info(
            f"Requested number of cores ({num_workers}) exceeds maximum required ({total_tasks}). Using {total_tasks} cores instead."
        )
        num_workers = total_tasks
    else:
        logging.info(
            f"Using {num_workers} worker process{'es' if num_workers != 1 else ''} for processing.{' (sequential mode)' if num_workers == 1 else ''}\n"
        )

    if num_workers == 1:
        # Sequential mode.
        for lib in args.libraries:
            library_name = os.path.basename(os.path.normpath(lib))
            logging.info(f"Starting processing for library '{library_name}'")
            library_results: Dict[str, pd.DataFrame] = {}
            tasks_for_lib = [(lib, compound_name, query_smiles, fp) for fp in fp_types]
            for task in tqdm(tasks_for_lib, desc=f"Processing {library_name}"):
                res = process_single_fingerprint(*task)
                if res is None:
                    continue
                _, fp, df = res
                library_results[fp] = df

            if not library_results:
                logging.error(
                    f"No results obtained for library '{library_name}'. Skipping..."
                )
                continue

            aggregated_results = aggregate_results(
                library_results, list(library_results.keys())
            )
            if os.path.isdir(output_dir):
                output_prefix = os.path.join(
                    output_dir, f"ConFiLiS_{compound_name}_{library_name}"
                )
            else:
                base, _ = os.path.splitext(output_dir)
                output_prefix = f"{base}_{library_name}"
            write_results_to_csv(
                aggregated_results,
                list(library_results.keys()),
                output_prefix,
                args.shortened_rows,
            )
            logging.info(
                f"Finished processing library '{library_name}'. Total compounds processed: {len(aggregated_results)}\n"
            )
    else:
        # Multicore mode: process and write each library as soon as its tasks finish.
        mp.set_start_method("spawn", force=True)
        library_results_done: Dict[str, Dict[str, pd.DataFrame]] = {}
        # Initialize task counters per library.
        library_tasks_completed: Dict[str, int] = {
            os.path.basename(os.path.normpath(lib)): 0 for lib in args.libraries
        }
        processed_libraries = set()
        futures = []
        future_to_lib = {}

        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            for lib in args.libraries:
                lib_name = os.path.basename(os.path.normpath(lib))
                for fp in fp_types:
                    fut = executor.submit(
                        process_single_fingerprint, lib, compound_name, query_smiles, fp
                    )
                    futures.append(fut)
                    future_to_lib[fut] = lib_name

            # Unified progress bar across all tasks.
            for future in tqdm(
                as_completed(futures), total=total_tasks, desc="Processing Fingerprints"
            ):
                lib_name = future_to_lib[future]
                try:
                    res = future.result()
                except Exception as e:
                    logging.error(f"Error in task for library '{lib_name}': {e}")
                    res = None

                library_tasks_completed[lib_name] += 1
                if res is not None:
                    _, fp, df = res
                    library_results_done.setdefault(lib_name, {})[fp] = df

                # If all tasks for a library have finished, process it immediately.
                if (library_tasks_completed[lib_name] == len(fp_types)) and (
                    lib_name not in processed_libraries
                ):
                    process_library_results(
                        lib_name,
                        library_results_done.get(lib_name, {}),
                        output_dir,
                        compound_name,
                        args.shortened_rows,
                    )
                    processed_libraries.add(lib_name)

    total_runtime = time.time() - start_time
    logging.info(f"Total runtime: {total_runtime:.2f} seconds")


if __name__ == "__main__":
    main()
