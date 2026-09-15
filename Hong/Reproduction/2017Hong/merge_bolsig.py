"""
Merge three BOLSIG+ cross-section files into one combined database.
Modifications:
  - N2: merge N2(A3,v0-4), N2(A3,v5-9), N2(A3,v10-) into N2(A3) by summing
        cross-sections at each energy point; remove N2(v1res)
  - H2: rename specific process names to match ZDPlasKin kinet.inp BOLSIG refs
  - N:  keep as-is
"""
import os
import re
from pathlib import Path

# File paths
base_dir = Path(r"G:\Kimi_project\ZDPlaskin模拟\Reproduction\2017Hong\data")
files = {
    "N2": base_dir / "N2_for_bolsig.txt",
    "H2": base_dir / "H2_for_bolsig.txt",
    "N":  base_dir / "N_for_bolsig.txt",
}
out_file = base_dir / "bolsigdb.dat"


def parse_bolsig_file(filepath):
    """
    Parse a BOLSIG+ cross-section file into a list of process dicts.
    Each dict has keys: 'type', 'name', 'threshold', 'data' (list of (energy, cs) tuples).
    """
    processes = []
    with open(filepath, "r", encoding="utf-8") as f:
        lines = f.readlines()

    i = 0
    n = len(lines)

    while i < n:
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        # Skip separator lines (dashes, equals, etc.)
        if re.match(r"^[\-=\*]+$", line):
            i += 1
            continue

        # If the line looks like a process type (ELASTIC, EXCITATION, IONIZATION, etc.)
        # BOLSIG+ process types are typically uppercase keywords
        if re.match(r"^[A-Z][A-Z\s]+$", line):
            proc_type = line
            i += 1
            if i >= n:
                break
            proc_name = lines[i].strip()
            i += 1
            if i >= n:
                break
            # threshold line
            threshold_line = lines[i].strip()
            try:
                threshold = float(threshold_line)
            except ValueError:
                # sometimes extra text after threshold, keep as string or parse first token
                m = re.match(r"([+-]?\d+\.?\d*[eE][+-]?\d+|\d+\.?\d*)", threshold_line)
                threshold = float(m.group(1)) if m else threshold_line
            i += 1

            # skip separator line(s)
            while i < n and re.match(r"^[\-=\*]+$", lines[i].strip()):
                i += 1
            if i >= n:
                break

            # read data until next process type or separator or EOF
            data = []
            while i < n:
                raw = lines[i].strip()
                if not raw:
                    i += 1
                    continue
                # If this line is a separator, it might be followed by next process
                if re.match(r"^[\-=\*]+$", raw):
                    # peek ahead: if next non-empty line is uppercase process type, break
                    j = i + 1
                    while j < n and not lines[j].strip():
                        j += 1
                    if j < n and re.match(r"^[A-Z][A-Z\s]+$", lines[j].strip()):
                        break
                    # otherwise it's a mid-process separator? skip it
                    i += 1
                    continue
                # If this line is a process type, break
                if re.match(r"^[A-Z][A-Z\s]+$", raw):
                    break

                # parse two columns: energy and cross-section
                parts = raw.split()
                if len(parts) >= 2:
                    try:
                        e = float(parts[0])
                        cs = float(parts[1])
                        data.append((e, cs))
                    except ValueError:
                        pass
                i += 1

            processes.append({
                "type": proc_type,
                "name": proc_name,
                "threshold": threshold,
                "data": data,
            })
            continue
        else:
            i += 1

    return processes


def merge_a3_processes(processes):
    """
    Find and merge the three N2(A3,v0-4/v5-9/v10-) processes into a single N2(A3) process.
    Remove the N2(v1res) process.
    """
    a3_names = ["N2 -> N2(A3,v0-4)", "N2 -> N2(A3,v5-9)", "N2 -> N2(A3,v10-)"]
    a3_procs = []
    new_procs = []
    removed_names = {"N2 -> N2(v1res)"}

    for p in processes:
        if p["name"] in a3_names:
            a3_procs.append(p)
        elif p["name"] in removed_names:
            continue
        else:
            new_procs.append(p)

    if len(a3_procs) != 3:
        raise RuntimeError(
            f"Expected 3 A3 sub-processes, found {len(a3_procs)}: {[p['name'] for p in a3_procs]}"
        )

    # Build energy -> list of cross-sections dict
    # Union of all energy points, sorted
    energy_map = {}
    for p in a3_procs:
        for e, cs in p["data"]:
            if e not in energy_map:
                energy_map[e] = [0.0, 0.0, 0.0]

    # Sum cross-sections per sub-process at each energy point
    for idx, p in enumerate(a3_procs):
        for e, cs in p["data"]:
            energy_map[e][idx] += cs

    merged_data = []
    for e in sorted(energy_map.keys()):
        total_cs = sum(energy_map[e])
        merged_data.append((e, total_cs))

    # Threshold: use the minimum of the three thresholds (or the first process threshold)
    # BOLSIG+ typically uses the lowest threshold for merged excitations
    merged_threshold = min(float(p["threshold"]) for p in a3_procs)

    merged_proc = {
        "type": "EXCITATION",
        "name": "N2 -> N2(A3)",
        "threshold": merged_threshold,
        "data": merged_data,
    }
    new_procs.append(merged_proc)
    return new_procs


def rename_h2_processes(processes):
    """
    Rename H2 processes to match ZDPlasKin kinet.inp BOLSIG references.
    """
    rename_map = {
        "H2 -> H2(B3SIG)(8.9eV)": "H2 -> H2(b3)",
        "H2 -> H2(B1SIG)(11.3eV)": "H2 -> H2(B1)",
        "H2 -> H2(C3PI)(11.75eV)": "H2 -> H2(c3)",
        "H2 -> H2(A3SIG)(11.8eV)": "H2 -> H2(a3)",
        "H2 -> H2(C1PI)(12.4eV)": "H2 -> H2(C1)",
        "H2 -> H2(V1)(0.516eV)": "H2 -> H2(v1)",
        "H2 -> H2(V2)(1eV)": "H2 -> H2(v2)",
        "H2 -> H2(V3)(1.5eV)": "H2 -> H2(v3)",
    }
    for p in processes:
        if p["name"] in rename_map:
            p["name"] = rename_map[p["name"]]
    return processes


def write_bolsig_file(processes, filepath):
    """
    Write merged processes in BOLSIG+ format.
    """
    with open(filepath, "w", encoding="utf-8") as f:
        first = True
        for p in processes:
            if not first:
                f.write("------------------------------\n")
            first = False

            f.write(f"{p['type']}\n")
            f.write(f"{p['name']}\n")
            # Format threshold in scientific notation similar to input
            thr = float(p["threshold"])
            f.write(f"{thr:.3e}\n")
            f.write("-----------------------------\n")
            for e, cs in p["data"]:
                # Use consistent formatting: energy in 4-decimal scientific, CS in 4-decimal scientific
                f.write(f"{e:.4e}\t{cs:.4e}\n")


def main():
    # Parse all three files
    all_procs = []
    for label, path in files.items():
        procs = parse_bolsig_file(path)
        print(f"Parsed {label}: {len(procs)} processes")
        all_procs.extend(procs)

    # Separate by species
    n2_procs = [p for p in all_procs if p["name"] == "N2" or p["name"].startswith("N2 ->")]
    h2_procs = [p for p in all_procs if p["name"] == "H2" or p["name"].startswith("H2 ->")]
    n_procs = [p for p in all_procs if p["name"] == "N" or p["name"].startswith("N ->")]

    print(f"  N2 processes: {len(n2_procs)}")
    print(f"  H2 processes: {len(h2_procs)}")
    print(f"  N processes:  {len(n_procs)}")

    # Apply modifications
    n2_procs = merge_a3_processes(n2_procs)
    h2_procs = rename_h2_processes(h2_procs)
    # N processes stay as-is

    # Combine in order: N2, H2, N
    final_procs = n2_procs + h2_procs + n_procs

    # Write output
    write_bolsig_file(final_procs, out_file)

    # Report
    size_bytes = os.path.getsize(out_file)
    print(f"\nFinal merged file: {out_file}")
    print(f"  Total processes: {len(final_procs)}")
    print(f"  File size: {size_bytes:,} bytes ({size_bytes/1024:.2f} KB)")

    # Print summary of all processes
    print(f"\nProcess list ({len(final_procs)} total):")
    for i, p in enumerate(final_procs, 1):
        print(f"  {i:2d}. {p['type']:12s} | {p['name']}")


if __name__ == "__main__":
    main()
