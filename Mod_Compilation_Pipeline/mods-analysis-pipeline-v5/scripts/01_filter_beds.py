# scripts/01_filter_beds.py

import os
import argparse
import subprocess
from glob import glob

def main(sample_dirs, mod_codes, output_dir, filter_script):
    os.makedirs(output_dir, exist_ok=True)

    # Validate filter script exists
    if not os.path.isfile(filter_script):
        raise FileNotFoundError(f"Missing filter script: {filter_script}")

    for sample_dir in sample_dirs:
        sample_name = os.path.basename(sample_dir.rstrip("/"))
        print(f"\n📦 Processing: {sample_name}")

        for mod_code in mod_codes:
            in_bed = os.path.join(sample_dir, f"{sample_name}_{mod_code}_filtered_mod.bed")
            if not os.path.isfile(in_bed):
                print(f"  ⚠️ Missing BED: {in_bed}")
                continue

            out_bed = os.path.join(output_dir, f"{sample_name}_{mod_code}_filtered_passonly.bed")

            cmd = [
                "python3", filter_script,
                "--in_bed", in_bed,
                "--mod_code", mod_code,
                "--out_bed", out_bed
            ]

            try:
                subprocess.run(cmd, check=True)
                print(f"  ✅ Wrote: {out_bed}")
            except subprocess.CalledProcessError as e:
                print(f"  ❌ Error filtering {in_bed}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter modkit BED files across samples")
    parser.add_argument("--sample_dirs", nargs='+', required=True, help="List of sample directories")
    parser.add_argument("--mod_codes", nargs='+', required=True, help="List of mod codes to filter")
    parser.add_argument("--output_dir", required=True, help="Directory to store filtered BEDs")
    parser.add_argument("--filter_script", required=True, help="Path to the filter_bed_sites.py script")

    args = parser.parse_args()
    main(args.sample_dirs, args.mod_codes, args.output_dir, args.filter_script)

