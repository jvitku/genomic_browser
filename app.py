"""
gnomAD Mitochondrial Variant Browser
Flask backend: parses VCF at startup, serves REST API + frontend.
"""

import gzip
import math
import os
import re
from flask import Flask, jsonify, render_template, request
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

VCF_PATH = os.path.join(os.path.dirname(__file__), "gnomad.genomes.v3.1.sites.chrM.vcf.bgz")

# Population and haplogroup order from VCF header
POPULATIONS = ["afr", "ami", "amr", "asj", "eas", "fin", "nfe", "oth", "sas", "mid"]
POP_NAMES = {
    "afr": "African/African American",
    "ami": "Amish",
    "amr": "Latino/Admixed American",
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "fin": "Finnish",
    "nfe": "Non-Finnish European",
    "oth": "Other",
    "sas": "South Asian",
    "mid": "Middle Eastern",
}
HAPLOGROUPS = [
    "A", "B", "C", "D", "E", "F", "G", "H", "HV", "I", "J", "K",
    "L0", "L1", "L2", "L3", "L4", "L5", "M", "N", "P", "R", "T", "U", "V", "W", "X", "Y", "Z",
]

# Known mitochondrial gene coordinates (GRCh38 chrM)
MT_GENES = {
    "MT-TF":   (577, 647),
    "MT-RNR1": (648, 1601),
    "MT-TV":   (1602, 1670),
    "MT-RNR2": (1671, 3229),
    "MT-TL1":  (3230, 3304),
    "MT-ND1":  (3307, 4262),
    "MT-TI":   (4263, 4331),
    "MT-TQ":   (4329, 4400),
    "MT-TM":   (4402, 4469),
    "MT-ND2":  (4470, 5511),
    "MT-TW":   (5512, 5579),
    "MT-TA":   (5587, 5655),
    "MT-TN":   (5657, 5729),
    "MT-TC":   (5761, 5826),
    "MT-TY":   (5826, 5891),
    "MT-CO1":  (5904, 7445),
    "MT-TS1":  (7446, 7514),
    "MT-TD":   (7518, 7585),
    "MT-CO2":  (7586, 8269),
    "MT-TK":   (8295, 8364),
    "MT-ATP8": (8366, 8572),
    "MT-ATP6": (8527, 9207),
    "MT-CO3":  (9207, 9990),
    "MT-TG":   (9991, 10058),
    "MT-ND3":  (10059, 10404),
    "MT-TR":   (10405, 10469),
    "MT-ND4L": (10470, 10766),
    "MT-ND4":  (10760, 12137),
    "MT-TH":   (12138, 12206),
    "MT-TS2":  (12207, 12265),
    "MT-TL2":  (12266, 12336),
    "MT-ND5":  (12337, 14148),
    "MT-ND6":  (14149, 14673),
    "MT-TE":   (14674, 14742),
    "MT-CYB":  (14747, 15887),
    "MT-TT":   (15888, 15953),
    "MT-TP":   (15956, 16023),
}

# VEP field order from VCF header
VEP_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature",
    "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position",
    "Protein_position", "Amino_acids", "Codons", "ALLELE_NUM", "DISTANCE", "STRAND",
    "VARIANT_CLASS", "MINIMISED", "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL", "TSL",
    "APPRIS", "CCDS", "ENSP", "SWISSPROT", "TREMBL", "UNIPARC", "GENE_PHENO",
    "SIFT", "PolyPhen", "DOMAINS", "HGVS_OFFSET", "MOTIF_NAME", "MOTIF_POS",
    "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "LoF", "LoF_filter", "LoF_flags", "LoF_info",
]


def parse_info(info_str):
    """Parse VCF INFO field into a dict."""
    d = {}
    for item in info_str.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            d[key] = val
        else:
            d[item] = True
    return d


def parse_vep(vep_str):
    """Parse first canonical VEP annotation, return dict of key fields."""
    if not vep_str:
        return {}
    annotations = vep_str.split(",")
    # Prefer canonical annotation
    best = None
    for ann_str in annotations:
        parts = ann_str.split("|")
        ann = {}
        for i, field_name in enumerate(VEP_FIELDS):
            ann[field_name] = parts[i] if i < len(parts) else ""
        if ann.get("CANONICAL") == "YES":
            best = ann
            break
        if best is None and ann.get("SYMBOL"):
            best = ann
    return best or {}


def safe_float(val, default=None):
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def safe_int(val, default=None):
    try:
        return int(val)
    except (ValueError, TypeError):
        return default


def parse_pipe_list(val):
    """Parse pipe-separated list of numbers."""
    if not val:
        return []
    return [safe_float(x) for x in val.split("|")]


def load_vcf(path):
    """Parse bgzipped VCF into list of variant dicts."""
    variants = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, vid, ref, alt, qual, filt, info_str = fields[:8]
            pos = int(pos)
            info = parse_info(info_str)
            vep = parse_vep(info.get("vep", ""))

            # Parse population frequencies
            pop_an = parse_pipe_list(info.get("pop_AN", ""))
            pop_ac_hom = parse_pipe_list(info.get("pop_AC_hom", ""))
            pop_ac_het = parse_pipe_list(info.get("pop_AC_het", ""))
            pop_af_hom = parse_pipe_list(info.get("pop_AF_hom", ""))
            pop_af_het = parse_pipe_list(info.get("pop_AF_het", ""))

            pop_data = {}
            for i, p in enumerate(POPULATIONS):
                pop_data[p] = {
                    "AN": pop_an[i] if i < len(pop_an) else None,
                    "AC_hom": pop_ac_hom[i] if i < len(pop_ac_hom) else None,
                    "AC_het": pop_ac_het[i] if i < len(pop_ac_het) else None,
                    "AF_hom": pop_af_hom[i] if i < len(pop_af_hom) else None,
                    "AF_het": pop_af_het[i] if i < len(pop_af_het) else None,
                }

            # Parse haplogroup data
            hap_an = parse_pipe_list(info.get("hap_AN", ""))
            hap_ac_hom = parse_pipe_list(info.get("hap_AC_hom", ""))
            hap_af_het = parse_pipe_list(info.get("hap_AF_het", ""))
            hap_af_hom = parse_pipe_list(info.get("hap_AF_hom", ""))

            hap_data = {}
            for i, h in enumerate(HAPLOGROUPS):
                hap_data[h] = {
                    "AN": hap_an[i] if i < len(hap_an) else None,
                    "AC_hom": hap_ac_hom[i] if i < len(hap_ac_hom) else None,
                    "AF_het": hap_af_het[i] if i < len(hap_af_het) else None,
                    "AF_hom": hap_af_hom[i] if i < len(hap_af_hom) else None,
                }

            # Heteroplasmy histogram
            hl_hist = parse_pipe_list(info.get("hl_hist", ""))

            variant = {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "filter": filt,
                "variant_collapsed": info.get("variant_collapsed", ""),
                # Allele counts
                "AN": safe_int(info.get("AN")),
                "AC_hom": safe_int(info.get("AC_hom")),
                "AC_het": safe_int(info.get("AC_het")),
                "AF_hom": safe_float(info.get("AF_hom")),
                "AF_het": safe_float(info.get("AF_het")),
                # Quality
                "dp_mean": safe_float(info.get("dp_mean")),
                "mq_mean": safe_float(info.get("mq_mean")),
                "tlod_mean": safe_float(info.get("tlod_mean")),
                # VEP annotation
                "consequence": vep.get("Consequence", ""),
                "impact": vep.get("IMPACT", ""),
                "symbol": vep.get("SYMBOL", ""),
                "gene_id": vep.get("Gene", ""),
                "biotype": vep.get("BIOTYPE", ""),
                "hgvsc": vep.get("HGVSc", ""),
                "hgvsp": vep.get("HGVSp", ""),
                "protein_position": vep.get("Protein_position", ""),
                "amino_acids": vep.get("Amino_acids", ""),
                "codons": vep.get("Codons", ""),
                "variant_class": vep.get("VARIANT_CLASS", ""),
                "sift": vep.get("SIFT", ""),
                "polyphen": vep.get("PolyPhen", ""),
                "lof": vep.get("LoF", ""),
                "lof_filter": vep.get("LoF_filter", ""),
                "lof_flags": vep.get("LoF_flags", ""),
                # Pathogenicity
                "mitotip_score": safe_float(info.get("mitotip_score")),
                "mitotip_trna_prediction": info.get("mitotip_trna_prediction", ""),
                "pon_mt_trna_prediction": info.get("pon_mt_trna_prediction", ""),
                "pon_ml_probability_of_pathogenicity": safe_float(
                    info.get("pon_ml_probability_of_pathogenicity")
                ),
                # Flags
                "hap_defining_variant": info.get("hap_defining_variant", False) is True,
                "common_low_heteroplasmy": info.get("common_low_heteroplasmy", False) is True,
                # Population & haplogroup (stored separately, served in detail view)
                "pop_data": pop_data,
                "hap_data": hap_data,
                "hl_hist": hl_hist,
                # Filters info
                "filters_info": info.get("filters", ""),
                "excluded_AC": safe_int(info.get("excluded_AC")),
            }
            variants.append(variant)
    return variants


# ---------- Load data at startup ----------
print("Loading VCF...")
VARIANTS = load_vcf(VCF_PATH)
print(f"Loaded {len(VARIANTS)} variants")

# Build indexes for fast lookup
VARIANT_INDEX = {}  # (pos, ref, alt) -> variant
GENE_VARIANTS = {}  # gene -> [indices]
for i, v in enumerate(VARIANTS):
    VARIANT_INDEX[(v["pos"], v["ref"], v["alt"])] = v
    gene = v["symbol"]
    if gene:
        GENE_VARIANTS.setdefault(gene, []).append(i)

# Collect unique values for filter dropdowns
ALL_CONSEQUENCES = sorted({v["consequence"] for v in VARIANTS if v["consequence"]})
ALL_IMPACTS = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
ALL_FILTERS = sorted({v["filter"] for v in VARIANTS})
ALL_GENES = sorted(GENE_VARIANTS.keys())


# ---------- API endpoints ----------

@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/summary")
def api_summary():
    total = len(VARIANTS)
    by_impact = {}
    by_consequence = {}
    by_filter = {}
    pass_count = 0
    for v in VARIANTS:
        imp = v["impact"] or "Unknown"
        by_impact[imp] = by_impact.get(imp, 0) + 1
        con = v["consequence"] or "Unknown"
        by_consequence[con] = by_consequence.get(con, 0) + 1
        f = v["filter"]
        by_filter[f] = by_filter.get(f, 0) + 1
        if f == "PASS":
            pass_count += 1

    return jsonify({
        "total": total,
        "pass_count": pass_count,
        "by_impact": by_impact,
        "by_consequence": dict(sorted(by_consequence.items(), key=lambda x: -x[1])),
        "by_filter": by_filter,
    })


@app.route("/api/density")
def api_density():
    """Return variant positions for density track (lightweight)."""
    bins = 200
    chrom_length = 16569
    bin_size = chrom_length / bins
    counts = [0] * bins
    for v in VARIANTS:
        b = min(int(v["pos"] / bin_size), bins - 1)
        counts[b] += 1
    return jsonify({"counts": counts, "bins": bins, "chrom_length": chrom_length})


@app.route("/api/genes")
def api_genes():
    genes = []
    for gene_name, coords in MT_GENES.items():
        variant_count = len(GENE_VARIANTS.get(gene_name, []))
        genes.append({
            "name": gene_name,
            "start": coords[0],
            "end": coords[1],
            "variant_count": variant_count,
            "type": "tRNA" if gene_name.startswith("MT-T") else (
                "rRNA" if gene_name.startswith("MT-RNR") else "protein_coding"
            ),
        })
    genes.sort(key=lambda g: g["start"])
    return jsonify({"genes": genes, "chrom_length": 16569})


@app.route("/api/variants")
def api_variants():
    # Parse filter parameters
    pos_start = safe_int(request.args.get("pos_start"))
    pos_end = safe_int(request.args.get("pos_end"))
    gene = request.args.get("gene", "").strip()
    consequence = request.args.get("consequence", "").strip()
    impact = request.args.get("impact", "").strip()
    af_min = safe_float(request.args.get("af_min"))
    af_max = safe_float(request.args.get("af_max"))
    filt = request.args.get("filter", "").strip()
    search = request.args.get("search", "").strip()
    sort_by = request.args.get("sort_by", "pos")
    sort_order = request.args.get("sort_order", "asc")
    page = max(1, safe_int(request.args.get("page"), 1))
    per_page = min(500, max(1, safe_int(request.args.get("per_page"), 50)))

    # Filter
    results = VARIANTS
    if pos_start is not None:
        results = [v for v in results if v["pos"] >= pos_start]
    if pos_end is not None:
        results = [v for v in results if v["pos"] <= pos_end]
    if gene:
        genes_filter = {g.strip() for g in gene.split(",")}
        results = [v for v in results if v["symbol"] in genes_filter]
    if consequence:
        results = [v for v in results if v["consequence"] == consequence]
    if impact:
        impacts_filter = {i.strip() for i in impact.split(",")}
        results = [v for v in results if v["impact"] in impacts_filter]
    if af_min is not None:
        results = [v for v in results if (v["AF_hom"] or 0) + (v["AF_het"] or 0) >= af_min]
    if af_max is not None:
        results = [v for v in results if (v["AF_hom"] or 0) + (v["AF_het"] or 0) <= af_max]
    if filt:
        if filt == "PASS":
            results = [v for v in results if v["filter"] == "PASS"]
        else:
            results = [v for v in results if filt in v["filter"]]
    if search:
        s = search.upper()
        # Support chrM:start-end format
        coord_match = re.match(r"(?:CHRM:)?(\d+)-(\d+)", s)
        if coord_match:
            s_start, s_end = int(coord_match.group(1)), int(coord_match.group(2))
            results = [v for v in results if s_start <= v["pos"] <= s_end]
        else:
            results = [
                v for v in results
                if s in v["symbol"].upper()
                or s in str(v["pos"])
                or s in v["ref"].upper()
                or s in v["alt"].upper()
                or s in v["consequence"].upper()
                or s in v.get("variant_collapsed", "").upper()
            ]

    # Sort
    sort_key = sort_by if sort_by in (
        "pos", "ref", "alt", "symbol", "consequence", "impact", "AF_hom", "AF_het", "AN", "filter"
    ) else "pos"
    reverse = sort_order == "desc"
    results.sort(key=lambda v: (v.get(sort_key) is None, v.get(sort_key, "")), reverse=reverse)

    # Paginate
    total = len(results)
    total_pages = max(1, math.ceil(total / per_page))
    start = (page - 1) * per_page
    page_results = results[start : start + per_page]

    # Strip heavy fields for list view
    slim = []
    for v in page_results:
        slim.append({
            "pos": v["pos"],
            "ref": v["ref"],
            "alt": v["alt"],
            "filter": v["filter"],
            "symbol": v["symbol"],
            "consequence": v["consequence"],
            "impact": v["impact"],
            "AF_hom": v["AF_hom"],
            "AF_het": v["AF_het"],
            "AN": v["AN"],
            "AC_hom": v["AC_hom"],
            "AC_het": v["AC_het"],
            "variant_collapsed": v["variant_collapsed"],
            "hap_defining_variant": v["hap_defining_variant"],
            "mitotip_trna_prediction": v["mitotip_trna_prediction"],
        })

    return jsonify({
        "variants": slim,
        "total": total,
        "page": page,
        "pages": total_pages,
        "per_page": per_page,
    })


@app.route("/api/variant/<int:pos>/<ref>/<alt>")
def api_variant_detail(pos, ref, alt):
    v = VARIANT_INDEX.get((pos, ref, alt))
    if not v:
        return jsonify({"error": "Variant not found"}), 404
    # Return everything
    return jsonify(v)


@app.route("/api/filters")
def api_filters():
    """Return available filter values for dropdowns."""
    return jsonify({
        "consequences": ALL_CONSEQUENCES,
        "impacts": ALL_IMPACTS,
        "filters": ALL_FILTERS,
        "genes": ALL_GENES,
        "populations": POPULATIONS,
        "pop_names": POP_NAMES,
    })


if __name__ == "__main__":
    app.run(debug=True, port=5001, threaded=True)
