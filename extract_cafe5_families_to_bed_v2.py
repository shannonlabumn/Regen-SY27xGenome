#!/usr/bin/env python3
import argparse, csv, re, sys
from pathlib import Path
from collections import defaultdict

ID_IN_ANGLE = re.compile(r'^(.*?)(?:<\s*(\d+)\s*>)?$')
COORD_ID = re.compile(r'^(?P<chr>[^_\s]+)_(?P<start>\d+)_(?P<end>\d+)_(?P<strand>[+-])$')

def parse_header_name(col):
    col = col.strip()
    m = ID_IN_ANGLE.match(col)
    if m:
        return ((m.group(1) or "").strip(), m.group(2))
    return (col, None)

def load_change(change_path):
    rows = []
    with open(change_path, newline='') as fh:
        rdr = csv.reader(fh, delimiter='\t')
        header = next(rdr)
        for r in rdr:
            if not r: continue
            rows.append(r)
    return header, rows

def find_matching_columns(header, node_query):
    node_query = str(node_query).strip()
    hits = []
    for i, col in enumerate(header):
        if i == 0:  # FamilyID
            continue
        name, nid = parse_header_name(col)
        if (node_query == col) or (node_query == name) or (node_query == (nid or "")):
            hits.append(i)
    return hits

def collect_families_for_node(change_header, change_rows, col_idxs):
    expanded, contracted = set(), set()
    for row in change_rows:
        fam = row[0]
        for c in col_idxs:
            try:
                v = float(row[c])
            except Exception:
                v = 0.0
            if v > 0:
                expanded.add(fam)
            elif v < 0:
                contracted.add(fam)
    return expanded, contracted

# ---------- significance filter ----------
def load_family_pvalues(family_results_path, p_col_name=None):
    if family_results_path is None:
        return {}
    with open(family_results_path, newline='') as fh:
        rdr = csv.reader(fh, delimiter='\t')
        header = next(rdr)
        cols_low = [c.lstrip('#').lower() for c in header]  # strip leading '#'
        # find FamilyID
        fam_idx = None
        for i,c in enumerate(cols_low):
            if c in ('familyid','family','id'):
                fam_idx = i; break
        if fam_idx is None:
            sys.exit("ERROR: Could not find FamilyID column in Base_family_results.txt")
        # find p-value
        p_idx = None
        if p_col_name:
            p_col_low = p_col_name.lower()
            for i,c in enumerate(cols_low):
                if c == p_col_low:
                    p_idx = i; break
            if p_idx is None:
                sys.exit(f"ERROR: p-value column '{p_col_name}' not found.")
        else:
            for i,c in enumerate(cols_low):
                if ('p' in c) and ('value' in c):
                    p_idx = i; break
            if p_idx is None:
                for key in ('pvalue','p.val','p'):
                    if key in cols_low:
                        p_idx = cols_low.index(key); break
        if p_idx is None:
            sys.exit("ERROR: Could not infer a p-value column. Use --p_col.")
        pvals = {}
        for row in rdr:
            if not row or len(row) <= max(fam_idx, p_idx): continue
            fam = row[fam_idx]
            try:
                p = float(row[p_idx])
            except Exception:
                continue
            pvals[fam] = p
        return pvals

def filter_by_pvalue(families, pval_map, alpha):
    if not pval_map:
        return set(families)
    return {f for f in families if (f in pval_map and pval_map[f] <= alpha)}

# ---------- Orthogroups ----------
def load_orthogroups(orthogroups_tsv):
    og = {}
    with open(orthogroups_tsv, newline='') as fh:
        rdr = csv.reader(fh, delimiter='\t')
        hdr = next(rdr)
        for row in rdr:
            if not row: continue
            fam = row[0]
            per = {}
            for i in range(1, len(hdr)):
                per[hdr[i]] = row[i] if i < len(row) else ""
            og[fam] = per
    return og

def genes_for_families(og_map, families, species_cols):
    species_to_genes = {sp: {} for sp in species_cols}
    for fam in families:
        per = og_map.get(fam, {})
        for sp in species_cols:
            raw = per.get(sp, "")
            if not raw: continue
            for tok in re.split(r'[,\s;]+', raw.strip()):
                if not tok: continue
                species_to_genes[sp][tok] = fam
    return species_to_genes

# ---------- coords from GFF (+chrom sizes) ----------
def parse_gff3_to_gene_coords(gff3_path):
    gene_coords = {}
    mrna_coords = {}
    chrom_max = defaultdict(int)

    def extract_id(attr):
        m = re.search(r'ID=([^;]+)', attr);
        if m: return m.group(1)
        m = re.search(r'gene_id=([^;]+)', attr);
        if m: return m.group(1)
        m = re.search(r'Name=([^;]+)', attr);
        if m: return m.group(1)
        return None

    with open(gff3_path) as fh:
        for line in fh:
            if not line or line.startswith('#'): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9: continue
            chrom, src, ftype, start, end, score, strand, phase, attrs = parts
            try:
                s = int(start); e = int(end)
            except:
                continue
            if e > chrom_max[chrom]:
                chrom_max[chrom] = e
            if ftype == 'gene':
                gid = extract_id(attrs)
                if gid:
                    gene_coords[gid] = (chrom, s, e, strand)
            elif ftype in ('mRNA','transcript'):
                tid = extract_id(attrs)
                if tid:
                    mrna_coords[tid] = (chrom, s, e, strand)
    merged = dict(mrna_coords)
    merged.update(gene_coords)  # prefer 'gene'
    return merged, chrom_max

def load_fai_sizes(fai_path):
    chrom_max = {}
    with open(fai_path) as fh:
        for line in fh:
            if not line.strip(): continue
            f = line.split('\t')
            chrom = f[0]; size = int(f[1])
            chrom_max[chrom] = size
    return chrom_max

# ---------- coordinate-ID fallback ----------
def coords_from_coord_id(gene_id):
    m = COORD_ID.match(gene_id)
    if not m:
        return None
    chrom = m.group('chr')
    s = int(m.group('start'))
    e = int(m.group('end'))
    strand = m.group('strand')
    if e < s:  # just in case
        s, e = e, s
    return (chrom, s, e, strand)

# ---------- write BED ----------
def write_bed(out_path, gene_to_family, coord_map, allow_coord_ids=True):
    n_found = n_missing = n_coord = 0
    with open(out_path, 'w') as out:
        for gene, fam in sorted(gene_to_family.items()):
            info = coord_map.get(gene)
            if not info:
                # try stripping transcript suffixes
                base = re.sub(r'(\.t\d+|\.\d+)$', '', gene)
                info = coord_map.get(base)
            if not info and allow_coord_ids:
                info = coords_from_coord_id(gene)
                if info: n_coord += 1
            if not info:
                n_missing += 1
                continue
            chrom, s1, e1, strand = info
            start0 = max(0, s1-1)
            out.write(f"{chrom}\t{start0}\t{e1}\t{gene}|{fam}\t0\t{strand if strand in ('+','-') else '.'}\n")
            n_found += 1
    return n_found, n_missing, n_coord

def write_bed_slopped(out_path, gene_to_family, coord_map, chrom_max, window_bp, allow_coord_ids=True):
    n_found = n_missing = n_coord = 0
    with open(out_path, 'w') as out:
        for gene, fam in sorted(gene_to_family.items()):
            info = coord_map.get(gene)
            if not info:
                base = re.sub(r'(\.t\d+|\.\d+)$', '', gene)
                info = coord_map.get(base)
            if not info and allow_coord_ids:
                info = coords_from_coord_id(gene)
                if info: n_coord += 1
            if not info:
                n_missing += 1
                continue
            chrom, s1, e1, strand = info
            start0 = max(0, (s1-1) - window_bp)
            end1   = e1 + window_bp
            if chrom in chrom_max:
                end1 = min(end1, chrom_max[chrom])
                start0 = max(0, start0)
            out.write(f"{chrom}\t{start0}\t{end1}\t{gene}|{fam}\t0\t{strand if strand in ('+','-') else '.'}\n")
            n_found += 1
    return n_found, n_missing, n_coord

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(
        description="Extract expanded/contracted families for a node; optional p-value filter; output gene lists and BED (with ±window slop). Supports coordinate-style IDs."
    )
    ap.add_argument("--change", required=True, help="Base_change.tab")
    ap.add_argument("--node",   required=True, help="Node to extract (exact header, name before <ID>, or numeric ID, e.g. 6)")
    ap.add_argument("--out_prefix", default="node")

    # significance options
    ap.add_argument("--family_results", help="Base_family_results.txt (optional)")
    ap.add_argument("--alpha", type=float, default=0.05, help="p-value cutoff (default 0.05)")
    ap.add_argument("--p_col", help="Column name to use as p-value in Base_family_results.txt (optional)")

    # mapping families -> genes
    ap.add_argument("--orthogroups", help="OrthoFinder Orthogroups.tsv (optional, needed for gene/BED outputs)")
    ap.add_argument("--species", nargs='+', help="One or more species columns from Orthogroups.tsv to export genes/BED for")

    # gene -> coordinates
    ap.add_argument("--gff", help="Gene annotation GFF3 for the species (optional; improves mapping & gives chrom sizes)")
    ap.add_argument("--fai", help="FASTA .fai index for chrom sizes (optional; overrides/in addition to GFF sizes)")

    # pre-slop window for BEDs
    ap.add_argument("--window", type=int, default=0,
                    help="If >0, also write slopped BEDs expanded by ±WINDOW bp and clamped by chrom sizes.")

    args = ap.parse_args()

    # 1) families from Base_change.tab
    header, rows = load_change(args.change)
    col_idxs = find_matching_columns(header, args.node)
    if not col_idxs:
        example_cols = ', '.join(header[:8])
        sys.exit(f"ERROR: No column matched '{args.node}'. Example headers: {example_cols} ...")
    expanded, contracted = collect_families_for_node(header, rows, col_idxs)

    # 2) p-value filter (optional)
    if args.family_results:
        pvals = load_family_pvalues(args.family_results, args.p_col)
        expanded = filter_by_pvalue(expanded, pvals, args.alpha)
        contracted = filter_by_pvalue(contracted, pvals, args.alpha)

    # 3) write family lists
    Path(f"{args.out_prefix}_expanded_families.txt").write_text("\n".join(sorted(expanded)) + ("\n" if expanded else ""))
    Path(f"{args.out_prefix}_contracted_families.txt").write_text("\n".join(sorted(contracted)) + ("\n" if contracted else ""))
    print(f"[INFO] Families written:")
    print(f"  {args.out_prefix}_expanded_families.txt : {len(expanded)} families")
    print(f"  {args.out_prefix}_contracted_families.txt : {len(contracted)} families")

    # 4) families -> genes -> BED
    if args.orthogroups and args.species:
        og = load_orthogroups(args.orthogroups)

        # coord sources
        coord_map = {}
        chrom_max = {}
        if args.gff:
            coord_map, chrom_max = parse_gff3_to_gene_coords(args.gff)
        if args.fai:
            chrom_max = load_fai_sizes(args.fai)  # override/augment sizes with .fai

        for sp in args.species:
            exp_map = genes_for_families(og, expanded, [sp])[sp]
            con_map = genes_for_families(og, contracted, [sp])[sp]

            # TSV gene lists
            with open(f"{args.out_prefix}_expanded_genes_{sp}.tsv", "w") as out:
                out.write("GeneID\tFamilyID\n")
                for g, fam in sorted(exp_map.items()):
                    out.write(f"{g}\t{fam}\n")
            with open(f"{args.out_prefix}_contracted_genes_{sp}.tsv", "w") as out:
                out.write("GeneID\tFamilyID\n")
                for g, fam in sorted(con_map.items()):
                    out.write(f"{g}\t{fam}\n")
            print(f"[INFO] Genes written for {sp}: "
                  f"{args.out_prefix}_expanded_genes_{sp}.tsv ({len(exp_map)} genes), "
                  f"{args.out_prefix}_contracted_genes_{sp}.tsv ({len(con_map)} genes)")

            # BEDs (with coordinate-ID fallback)
            base_exp_bed = f"{args.out_prefix}_expanded_{sp}.bed"
            base_con_bed = f"{args.out_prefix}_contracted_{sp}.bed"
            n1, m1, c1 = write_bed(base_exp_bed, exp_map, coord_map, allow_coord_ids=True)
            n2, m2, c2 = write_bed(base_con_bed, con_map, coord_map, allow_coord_ids=True)
            print(f"[INFO] BED for {sp}: {base_exp_bed} (found {n1}, missing {m1}, coordID {c1}); "
                  f"{base_con_bed} (found {n2}, missing {m2}, coordID {c2})")

            if args.window and args.window > 0:
                exp_slop = f"{args.out_prefix}_expanded_{sp}_w{args.window}bp.bed"
                con_slop = f"{args.out_prefix}_contracted_{sp}_w{args.window}bp.bed"
                n3, m3, c3 = write_bed_slopped(exp_slop, exp_map, coord_map, chrom_max, args.window, allow_coord_ids=True)
                n4, m4, c4 = write_bed_slopped(con_slop, con_map, coord_map, chrom_max, args.window, allow_coord_ids=True)
                print(f"[INFO] Slopped BED (±{args.window} bp) for {sp}: "
                      f"{exp_slop} (found {n3}, missing {m3}, coordID {c3}); "
                      f"{con_slop} (found {n4}, missing {m4}, coordID {c4})")
    else:
        if args.gff or args.fai:
            print("[WARN] --gff/--fai provided but --orthogroups/--species missing; skipping BED output.")

if __name__ == "__main__":
    main()
