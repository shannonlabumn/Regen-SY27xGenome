#!/usr/bin/env python3
"""
CAFE5 tree figure (v15.1): manual branch drawing so scaling affects branches + labels.

Adds:
  --branch_lw   (thicker branches)
  --x_scale     (stretch horizontal axis / branch lengths)

Keeps everything from v14:
  families vs genes, label shifts, species alignment, fig size, y_scale, etc.
"""
import argparse, io, re, sys
from pathlib import Path
import pandas as pd
from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

UNICODE_MINUS = "−"
LEAF_WITH_ID  = re.compile(r"^(.*?)(?:<\s*(\d+)\s*>).*?$")
INT_IN_BRACKS = re.compile(r"<\s*(\d+)\s*>")

# ---------------------- helpers ----------------------
def extract_id_from_token(token: str):
    m = INT_IN_BRACKS.search(str(token))
    return m.group(1) if m else None

def split_header(colname: str):
    m = LEAF_WITH_ID.match(colname.strip())
    if m:
        sp = (m.group(1) or "").strip()
        return extract_id_from_token(colname), (sp if sp else None)
    return extract_id_from_token(colname), None

# ---------------------- load data --------------------
def load_families_present(count_path: str):
    df = pd.read_csv(count_path, sep="\t", engine="python")
    if "FamilyID" not in df.columns:
        sys.exit("ERROR: 'FamilyID' not found in Base_count.tab")
    out = {}
    for col in df.columns:
        if col == "FamilyID": continue
        node_id, _ = split_header(col)
        if not node_id: continue
        s = df[col]
        out[node_id] = int((s.fillna(0) > 0).sum())
    return out

def load_clade_results(clade_path: str):
    df = pd.read_csv(clade_path, sep="\t", engine="python")
    df.rename(columns={c: c.lstrip("#") for c in df.columns}, inplace=True)
    need = {"Taxon_ID","Increase","Decrease"}
    if not need.issubset(df.columns):
        sys.exit("ERROR: Expected columns 'Taxon_ID, Increase, Decrease' in Base_clade_results.txt")
    out = {}
    for _, row in df.iterrows():
        nid = extract_id_from_token(row["Taxon_ID"])
        if not nid: continue
        out[nid] = (int(row["Increase"]), int(row["Decrease"]))
    return out

def load_gene_gains_losses(change_path: str):
    df = pd.read_csv(change_path, sep="\t", engine="python")
    if "FamilyID" not in df.columns:
        sys.exit("ERROR: 'FamilyID' not found in Base_change.tab")
    out = {}
    for col in df.columns:
        if col == "FamilyID": continue
        node_id, _ = split_header(col)
        if not node_id: continue
        s = df[col].fillna(0)
        gains  = int(s.clip(lower=0).sum())
        losses = int((-s.clip(upper=0)).sum())
        out[node_id] = (gains, losses)
    return out

# ---------------------- tree math --------------------
def assign_unit_branch_lengths(tree):
    for cl in tree.find_clades(order="preorder"):
        if cl.branch_length is None:
            cl.branch_length = 1.0

def get_x_positions(tree):
    assign_unit_branch_lengths(tree)
    return tree.depths()  # dict: clade -> cumulative branch length

def get_y_positions(tree):
    y = {}
    for i, tip in enumerate(tree.get_terminals()):
        y[tip] = float(i)
    def calc_y(clade):
        if clade in y: return y[clade]
        ys = [calc_y(c) for c in clade.clades]
        y[clade] = sum(ys)/len(ys)
        return y[clade]
    calc_y(tree.root)
    return y

def preorder(clade):
    out = [clade]
    for c in getattr(clade, "clades", []):
        out.extend(preorder(c))
    return out

def extract_node_id_from_clade(clade):
    for token in (clade.name, getattr(clade, "confidence", None), str(clade)):
        nid = extract_id_from_token(token) if token is not None else None
        if nid: return nid
    return ""

# ---------------------- drawing primitives ------------
def measure_line(ax, x_left, y, chunks, fontsize):
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    inv = ax.transData.inverted()
    cursor_x = x_left; total_w = 0.0; max_h = 0.0; ghosts = []
    for txt, color in chunks:
        t = ax.text(cursor_x, y, txt, color=color, fontsize=fontsize,
                    va="baseline", ha="left", zorder=-1, alpha=0.0)
        bb = t.get_window_extent(renderer=renderer)
        (x0,y0),(x1,y1) = inv.transform([(bb.x0,bb.y0),(bb.x1,bb.y1)])
        w,h = x1-x0, y1-y0
        total_w += w; max_h = max(max_h,h); cursor_x += w; ghosts.append(t)
    for g in ghosts: g.remove()
    return total_w, max_h

def draw_line(ax, x_left, y, chunks, fontsize, z=8):
    cursor_x = x_left
    for txt, color in chunks:
        ax.text(cursor_x, y, txt, color=color, fontsize=fontsize,
                va="baseline", ha="left", zorder=z)
        w,_ = measure_line(ax, cursor_x, y, [(txt,color)], fontsize)
        cursor_x += w

def draw_box(ax, x_left, y_center, width, height, pad_x, pad_y, alpha=0.25, z=6):
    rect = FancyBboxPatch((x_left - pad_x, y_center - (height/2.0 + pad_y)),
                          width + 2*pad_x, height + 2*pad_y,
                          boxstyle="round,pad=0.0,rounding_size=0.15",
                          linewidth=0.0, facecolor="lightgrey", alpha=alpha, zorder=z)
    ax.add_patch(rect)

def chunks_line(n_fams:int, plus:int, minus:int, unit_label:str, show_unit:bool):
    suffix = f" {unit_label})" if show_unit else ")"
    return [
        (f"{n_fams} families (", "black"),
        (f"+{plus}", "forestgreen"),
        (f", {UNICODE_MINUS}{minus}{suffix}", "firebrick"),
    ]

# ---------------------- main draw ---------------------
def draw_tree(report_cafe, count_tab, clade_results, out_png, out_svg, out_pdf,
              fontsize, paren_mode, show_paren_unit, change_tab,
              species_align, species_dx_frac, species_xshift_frac, species_yshift,
              label_yshift, label_xshift_frac,
              fig_width, fig_height, y_scale, x_scale, branch_lw):

    fam_present   = load_families_present(count_tab)
    incdec_fams   = load_clade_results(clade_results)
    gainsloss_gen = load_gene_gains_losses(change_tab) if (paren_mode=="genes" and change_tab) else {}

    # parse trees from report
    text = Path(report_cafe).read_text()
    ids_nw = tree_nw = None
    for line in text.splitlines():
        if line.startswith("# IDs of nodes:"): ids_nw  = line.replace("# IDs of nodes:","").strip()
        elif line.startswith("Tree:"):          tree_nw = line.replace("Tree:","").strip()
    if not ids_nw or not tree_nw:
        sys.exit("ERROR: Could not parse '# IDs of nodes:' and 'Tree:' from Base_report.cafe")

    ids_tree  = Phylo.read(io.StringIO(ids_nw), "newick")
    disp_tree = Phylo.read(io.StringIO(tree_nw), "newick")
    disp_tree.ladderize()

    # positions
    x = get_x_positions(disp_tree)
    y = get_y_positions(disp_tree)

    # apply scaling
    if x_scale != 1.0:
        x = {k: v*x_scale for k,v in x.items()}
    if y_scale != 1.0:
        y = {k: v*y_scale for k,v in y.items()}

    # map clade -> parent (Biopython Phylo does not store .up)
    parent_map = {}
    for parent in disp_tree.find_clades(order="preorder"):
        for child in parent.clades:
            parent_map[child] = parent

    # id mapping
    disp_list = preorder(disp_tree.root)
    ids_list  = preorder(ids_tree.root)
    if len(disp_list) != len(ids_list):
        sys.exit("ERROR: display vs id tree sizes differ; cannot align.")
    node_ids = [extract_node_id_from_clade(c) for c in ids_list]

    total_depth = max(x.values()) if x else 1.0

    # figure
    fig = plt.figure(figsize=(fig_width, fig_height), frameon=False)
    ax  = fig.add_subplot(111)

    # draw branches manually so scaling applies
    for cl in disp_tree.find_clades(order="preorder"):
        parent = parent_map.get(cl, None)
        if parent is None:
            continue  # root
        px = x[parent]
        cx = x[cl]
        py = y[parent]
        cy = y[cl]
        # horizontal segment at child y
        ax.plot([px, cx], [cy, cy], color="k", lw=branch_lw, solid_capstyle="round")
        # vertical segment at parent x
        ax.plot([px, px], [py, cy], color="k", lw=branch_lw, solid_capstyle="round")

    # frame cosmetics
    ax.set_ylabel("")
    ax.tick_params(axis="y", which="both", left=False, right=False, labelleft=False)
    for b in ["left","right","top"]:
        ax.spines[b].set_visible(False)
    ax.set_title("Gene family expansions and contractions per lineage", fontsize=fontsize+2)

    # species names
    tip_max_x = max(x.get(t,0.0) for t in disp_tree.get_terminals())
    dx_species = (species_dx_frac + species_xshift_frac) * total_depth
    name_end = {}

    if species_align.lower()=="column":
        sx = tip_max_x + dx_species
        for tip in disp_tree.get_terminals():
            raw = tip.name or ""
            m = LEAF_WITH_ID.match(raw)
            pretty = (m.group(1).strip() if m else raw).replace("_"," ")
            yt = y[tip] - species_yshift
            ax.text(sx, yt, pretty, fontsize=fontsize, va="baseline", ha="left", zorder=9)
            w,_ = measure_line(ax, sx, yt, [(pretty,"black")], fontsize)
            name_end[tip] = sx + w
    else:
        for tip in disp_tree.get_terminals():
            raw = tip.name or ""
            m = LEAF_WITH_ID.match(raw)
            pretty = (m.group(1).strip() if m else raw).replace("_"," ")
            sx = x[tip] + dx_species
            yt = y[tip] - species_yshift
            ax.text(sx, yt, pretty, fontsize=fontsize, va="baseline", ha="left", zorder=9)
            w,_ = measure_line(ax, sx, yt, [(pretty,"black")], fontsize)
            name_end[tip] = sx + w

    pad_x = 0.02*total_depth
    pad_y = 0.14

    # callouts
    for d_clade, nid in zip(disp_list, node_ids):
        if not nid: continue
        n_fams = fam_present.get(nid, 0)
        if paren_mode == "families":
            plus, minus = incdec_fams.get(nid, (0,0))
            unit_label = " families"
        else:
            plus, minus = gainsloss_gen.get(nid, (0,0))
            unit_label = " genes"
        parts = chunks_line(n_fams, plus, minus, unit_label, show_paren_unit)

        if d_clade.is_terminal():
            xa = name_end.get(d_clade, (max(name_end.values()) if name_end else 0)) + dx_species
            ya = y[d_clade] - label_yshift
            w,h = measure_line(ax, xa, ya, parts, fontsize)
            draw_box(ax, xa, ya, w, h, pad_x, pad_y, alpha=0.25, z=6)
            draw_line(ax, xa, ya, parts, fontsize, z=8)
        else:
            node_x = x[d_clade] + (label_xshift_frac * total_depth)
            ya = y[d_clade] - label_yshift
            w,h = measure_line(ax, 0.0, ya, parts, fontsize)
            xl = node_x - w/2.0
            draw_box(ax, xl, ya, w, h, pad_x, pad_y, alpha=0.25, z=6)
            draw_line(ax, xl, ya, parts, fontsize, z=8)

    # limits
    ymin = min(y.values()) if y else 0.0
    ymax = max(y.values()) if y else 1.0
    ax.set_ylim(ymin - 1, ymax + 1)
    xmax = max(x.values()) if x else 1.0
    ax.set_xlim(0, xmax*1.1)
    ax.axis("off")

    # save
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    if out_svg: fig.savefig(out_svg, bbox_inches="tight")
    if out_pdf: fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

# ---------------------- CLI ---------------------------
def str2bool(s): return str(s).lower() in ("1","true","yes","y")

def main():
    ap = argparse.ArgumentParser(description="Draw CAFE5 tree with full layout controls (v15.1).")
    ap.add_argument("--report", required=True, help="Base_report.cafe")
    ap.add_argument("--count",  required=True, help="Base_count.tab")
    ap.add_argument("--clade",  required=True, help="Base_clade_results.txt")
    ap.add_argument("--change", required=False, default=None,
                    help="Base_change.tab (only used if --paren_mode genes)")

    ap.add_argument("-o","--outpng", default="cafe5_expansion_contraction.png")
    ap.add_argument("--svg", dest="outsvg", default=None)
    ap.add_argument("--pdf", dest="outpdf", default=None)
    ap.add_argument("--fontsize", type=int, default=12)

    ap.add_argument("--paren_mode", choices=["families","genes"], default="families")
    ap.add_argument("--show_paren_unit", type=str2bool, default=False)

    ap.add_argument("--species_align", choices=["per-tip","column"], default="per-tip")
    ap.add_argument("--species_dx_frac",    type=float, default=0.012)
    ap.add_argument("--species_xshift_frac",type=float, default=0.0)
    ap.add_argument("--species_yshift",     type=float, default=0.18)

    ap.add_argument("--label_yshift",       type=float, default=0.18)
    ap.add_argument("--label_xshift_frac",  type=float, default=0.0)

    ap.add_argument("--fig_width",  type=float, default=26.0)
    ap.add_argument("--fig_height", type=float, default=16.0)
    ap.add_argument("--y_scale",    type=float, default=1.0,
                    help="multiply vertical spacing (nodes & labels)")
    ap.add_argument("--x_scale",    type=float, default=1.0,
                    help="multiply horizontal branch lengths")
    ap.add_argument("--branch_lw",  type=float, default=2.0,
                    help="branch line width")

    args = ap.parse_args()
    if args.paren_mode=="genes" and not args.change:
        sys.exit("ERROR: --paren_mode genes requires --change Base_change.tab")

    draw_tree(args.report, args.count, args.clade, args.outpng, args.outsvg, args.outpdf,
              args.fontsize, args.paren_mode, args.show_paren_unit, args.change,
              args.species_align, args.species_dx_frac, args.species_xshift_frac, args.species_yshift,
              args.label_yshift, args.label_xshift_frac,
              args.fig_width, args.fig_height, args.y_scale, args.x_scale, args.branch_lw)

if __name__ == "__main__":
    main()
