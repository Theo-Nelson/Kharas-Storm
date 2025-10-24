import os
import argparse

def passes_filter(fields):
    try:
        valid_coverage = int(fields[9])
        count_modified = int(fields[11])
        count_fail = int(fields[15])
        count_diff = int(fields[16])
    except (ValueError, IndexError):
        return False

    if count_diff > 3 * valid_coverage:
        return False
    if count_modified <= count_fail + 1:
        return False
    return True

def filter_bed(in_bed, out_bed, mod_code):
    with open(in_bed, 'r') as infile, open(out_bed, 'w') as outfile:
        for line in infile:
            if line.startswith("#") or line.startswith("chrom"):
                outfile.write(line)
                continue

            fields = line.strip().split("\t")
            if len(fields) < 18 or fields[3] != mod_code:
                continue

            if passes_filter(fields):
                outfile.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter modkit BED file by custom criteria.")
    parser.add_argument("--in_bed", required=True, help="Input modkit BED file")
    parser.add_argument("--mod_code", required=True, help="Mod code to filter (e.g. a, m, 17596, etc.)")
    parser.add_argument("--out_bed", required=True, help="Output BED file with only passing rows")

    args = parser.parse_args()
    filter_bed(args.in_bed, args.out_bed, args.mod_code)

