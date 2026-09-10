#!/usr/bin/env python
"""
Parallel wrapper for ternary_model_prediction.py

Splits the decoy list file into N parts and runs N instances of
ternary_model_prediction.py simultaneously, then merges all outputs
back to the original expected files (rmsd.txt, ternary*.pdb).

Usage example (replace the original ternary_model_prediction.py call):
    python ../run_parallel_ternary.py \
        -la ./atomlist/auto_linker_atom_list.txt \
        -da ./atomlist/auto_decoy_atom_list.txt \
        -dl decoy_list.txt \
        -ll linker_list.txt \
        -wd ./atomlist/auto_decoy_atom_list_delete.txt \
        -c 0.4 \
        -t default \
        -r rmsd.txt \
        -n 8
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from math import ceil
from concurrent.futures import ThreadPoolExecutor, as_completed


def split_file(file_path, n, output_dir, prefix="decoy_list"):
    """Split a text file into n roughly equal parts, return list of output paths.
    
    Relative paths in the file are converted to absolute paths based on the
    directory containing file_path, so subprocesses work correctly from any cwd.
    """
    # Determine the base directory for resolving relative paths
    base_dir = os.path.dirname(os.path.abspath(file_path))

    with open(file_path, "r") as f:
        raw_lines = f.readlines()

    # Convert relative paths to absolute paths
    lines = []
    for line in raw_lines:
        stripped = line.strip()
        if stripped:
            if not os.path.isabs(stripped):
                stripped = os.path.normpath(os.path.join(base_dir, stripped))
            lines.append(stripped + "\n")
        else:
            lines.append(line)

    total = len(lines)
    chunk_size = max(1, ceil(total / n))
    out_files = []

    for i in range(n):
        start = i * chunk_size
        end = min(start + chunk_size, total)
        if start >= total:
            break
        out_name = os.path.join(output_dir, "{}_{}.txt".format(prefix, i + 1))
        with open(out_name, "w") as fout:
            fout.writelines(lines[start:end])
        out_files.append(out_name)

    return out_files


def merge_rmsd_files(rmsd_files, output_path):
    """Concatenate multiple temporary RMSD files into one."""
    with open(output_path, "w") as fout:
        for f in rmsd_files:
            if os.path.exists(f) and os.path.getsize(f) > 0:
                with open(f, "r") as fin:
                    fout.write(fin.read())


def collect_pdb_files(pdb_dirs, output_dir, prefix="ternary"):
    """Collect all .pdb files from temp dirs, rename sequentially in output_dir."""
    pdb_files = []
    for d in pdb_dirs:
        if not os.path.isdir(d):
            continue
        for fname in sorted(os.listdir(d)):
            if fname.endswith(".pdb"):
                pdb_files.append(os.path.join(d, fname))

    for idx, src in enumerate(pdb_files):
        dst = os.path.join(output_dir, "{}{}.pdb".format(prefix, idx))
        shutil.copy2(src, dst)


def run_instance(script_path, args_dict, work_dir, instance_id):
    """Run one instance of ternary_model_prediction.py and return its results."""
    decoy_list = args_dict["decoy_list_parts"][instance_id]
    rmsd_out = os.path.join(work_dir, "rmsd_{}.txt".format(instance_id))
    pdb_out_dir = os.path.join(work_dir, "pdb_out_{}".format(instance_id))
    os.makedirs(pdb_out_dir, exist_ok=True)

    cmd = [
        sys.executable, script_path,
        "-la", args_dict["linker_alignment"],
        "-da", args_dict["decoy_alignment"],
        "-dl", decoy_list,
        "-ll", args_dict["linker_list"],
        "-c", str(args_dict["cutoff"]),
        "-t", args_dict["ternary"],
        "-r", rmsd_out,
    ]

    if args_dict.get("warheads_delete"):
        cmd += ["-wd", args_dict["warheads_delete"]]
    if args_dict.get("linker_delete"):
        cmd += ["-ld", args_dict["linker_delete"]]
    if args_dict.get("alignment_iterations"):
        cmd += ["-ai", str(args_dict["alignment_iterations"])]

    result = subprocess.run(
        cmd,
        cwd=pdb_out_dir,
        capture_output=True,
        text=True,
    )

    return {
        "instance_id": instance_id,
        "rmsd_file": rmsd_out,
        "pdb_dir": pdb_out_dir,
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Run ternary_model_prediction.py in parallel "
                    "(split decoy list by N)."
    )

    # Original ternary_model_prediction.py arguments
    parser.add_argument("-la", "--linker_alignment", required=True,
                        help="Text file of linker atoms for the alignment")
    parser.add_argument("-da", "--decoy_alignment", required=True,
                        help="Text file of decoy atoms for the alignment")
    parser.add_argument("-dl", "--decoy_list", required=True,
                        help="Decoy PDB files listed in a single txt file")
    parser.add_argument("-ll", "--linker_list", required=True,
                        help="Linker PDB files listed in a single txt file")
    parser.add_argument("-c", "--cutoff", type=float, default=0.4,
                        help="Cut off value of RMSD (default: 0.4)")
    parser.add_argument("-r", "--rmsd", default="rmsd.txt",
                        help="Output RMSD file (default: rmsd.txt)")
    parser.add_argument("-wd", "--warheads_delete", default="",
                        help="Warhead atoms to be deleted")
    parser.add_argument("-ld", "--linker_delete", default="",
                        help="Linker atoms to be deleted")
    parser.add_argument("-t", "--ternary",
                        default="default",
                        choices=["none", "default", "specify"],
                        help="Output ternary mode (default: default)")
    parser.add_argument("-ai", "--alignment_iterations", type=int, default=5,
                        help="Number of iterations used for alignment (default: 5)")
    # Parallelization argument
    parser.add_argument("-n", "--num_parallel", type=int, default=4,
                        help="Number of parallel instances (default: 4)")

    args = parser.parse_args()
    n = args.num_parallel

    print("[run_parallel_ternary] Splitting {} into {} parts ...".format(
        args.decoy_list, n
    ))

    # Create a temporary working directory for split files and intermediate outputs
    with tempfile.TemporaryDirectory(prefix="ternary_parallel_") as tmp_dir:

        # Resolve relative paths in linker_list.txt -> temp file with absolute paths
        linker_list_resolved = os.path.join(tmp_dir, "linker_list_resolved.txt")
        linkers_base_dir = os.path.dirname(os.path.abspath(args.linker_list))
        with open(args.linker_list, "r") as fin:
            with open(linker_list_resolved, "w") as fout:
                for line in fin:
                    stripped = line.strip()
                    if stripped:
                        if not os.path.isabs(stripped):
                            stripped = os.path.normpath(os.path.join(
                                linkers_base_dir, stripped))
                        fout.write(stripped + "\n")
        print("[run_parallel_ternary] Resolved linker_list -> {} (abs paths)".format(
            linker_list_resolved))
        # 1. Split the decoy list
        decoy_parts = split_file(args.decoy_list, n, tmp_dir, prefix="decoy_list")

        actual_n = len(decoy_parts)
        if actual_n < n:
            print("[run_parallel_ternary] Warning: only need {} part(s) "
                  "for {} lines".format(actual_n, args.decoy_list))

        # 2. Prepare arguments dict for worker processes
        script_dir = os.path.dirname(os.path.abspath(__file__))
        script_path = os.path.join(script_dir, "ternary_model_prediction.py")

        # Convert relative paths to absolute paths so worker subprocesses
        # work correctly regardless of cwd settings
        def resolve_path(p):
            if p and not os.path.isabs(p):
                return os.path.abspath(p)
            return p

        args_dict = {
            "decoy_list_parts": decoy_parts,
            "decoy_alignment": resolve_path(args.decoy_alignment),
            "linker_alignment": resolve_path(args.linker_alignment),
            "linker_list": linker_list_resolved,  # already absolute paths
            "cutoff": args.cutoff,
            "ternary": args.ternary,
            "warheads_delete": (
                resolve_path(args.warheads_delete) if args.warheads_delete else ""
            ),
            "linker_delete": (
                resolve_path(args.linker_delete) if args.linker_delete else ""
            ),
            "alignment_iterations": args.alignment_iterations,
        }

        # 3. Run all instances in parallel
        print("[run_parallel_ternary] Launching {} parallel instances ...".format(
            actual_n
        ))
        futures = []
        with ThreadPoolExecutor(max_workers=actual_n) as executor:
            for i in range(actual_n):
                future = executor.submit(
                    run_instance, script_path, args_dict, tmp_dir, i
                )
                futures.append(future)

            # Collect results
            results = []
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                status = "OK" if result["returncode"] == 0 else "FAILED"
                print("[run_parallel_ternary] Instance {}/{} finished ({})".format(
                    result["instance_id"] + 1, actual_n, status
                ))

        # Check for failures
        failed = [r for r in results if r["returncode"] != 0]
        if failed:
            print("[run_parallel_ternary] Warning: {} instance(s) failed:".format(
                len(failed)
            ))
            for f in failed:
                # Show full error (last 1500 chars) for debugging
                stderr_short = f["stderr"].strip()
                if len(stderr_short) > 1500:
                    stderr_short = "...[truncated]...\n" + stderr_short[-1500:]
                print("  Instance {}:\n{}".format(
                    f["instance_id"] + 1, stderr_short
                ))

        # 4. Merge RMSD files
        rmsd_output = resolve_path(args.rmsd)
        rmsd_files = sorted([r["rmsd_file"] for r in results])
        merge_rmsd_files(rmsd_files, rmsd_output)
        print("[run_parallel_ternary] Merged RMSD results -> {}".format(
            rmsd_output
        ))

        # 5. Collect and rename all PDB files
        pdb_dirs = sorted([r["pdb_dir"] for r in results])
        output_cwd = os.getcwd()
        collect_pdb_files(pdb_dirs, output_cwd, prefix="ternary")
        pdb_count = len([f for f in os.listdir(output_cwd)
                         if f.startswith("ternary") and f.endswith(".pdb")])
        print("[run_parallel_ternary] Collected {} ternary PDB file(s) "
              "-> {}".format(pdb_count, output_cwd))

    # temp dir cleaned up automatically
    print("[run_parallel_ternary] All done.")


if __name__ == "__main__":
    main()
