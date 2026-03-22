#!/usr/bin/env python3
"""
alphafold_surprise_me.py
A pipeline that fetches a random AlphaFold homodimer prediction,
visualizes confidence, and retrieves related PubMed literature.

Requirements:
    uv pip install requests biopython matplotlib numpy
"""

import random
import time
import sys
import json
import textwrap

import requests
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from Bio import Entrez

# ── Configuration ──────────────────────────────────────────────────────────────
Entrez.email = "your.email@example.com"   # NCBI requires a valid email
MAX_UNIPROT_RESULTS = 200                  # pool to sample from
MAX_PUBMED_RESULTS  = 5


# ── Step 1: Sample a random homodimeric protein from UniProt ───────────────────

def fetch_homodimer_pool(size: int = MAX_UNIPROT_RESULTS) -> list[dict]:
    """
    Query UniProt for reviewed human proteins annotated as homodimers.
    Returns a list of dicts with 'accession', 'geneName', 'proteinName'.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": 'reviewed:true organism_id:9606 cc_subunit:"homodimer"',
        "fields": "accession,gene_names,protein_name",
        "format": "json",
        "size": size,
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    results = resp.json().get("results", [])

    pool = []
    for entry in results:
        acc = entry.get("primaryAccession", "")
        gene = ""
        genes = entry.get("genes", [])
        if genes:
            gene_names = genes[0].get("geneName", {})
            gene = gene_names.get("value", "")
        protein_name = (
            entry.get("proteinDescription", {})
                 .get("recommendedName", {})
                 .get("fullName", {})
                 .get("value", "Unknown protein")
        )
        if acc:
            pool.append({"accession": acc, "gene": gene, "protein": protein_name})

    return pool


def pick_random_entry(pool: list[dict]) -> dict:
    return random.choice(pool)


# ── Step 2: Fetch AlphaFold complex prediction ─────────────────────────────────

def fetch_alphafold_entry(uniprot_id: str) -> dict | None:
    """
    Query the AlphaFold DB API for a given UniProt accession.
    Returns the first entry dict if found, else None.
    Handles both monomer (AF-P00520-F1) and complex (AF-00000...) entries.
    """
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    resp = requests.get(url, timeout=30)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    entries = resp.json()
    if not entries:
        return None

    # Prefer a complex entry (isComplex=True) when available
    complex_entries = [e for e in entries if e.get("isComplex")]
    return complex_entries[0] if complex_entries else entries[0]


def download_confidence(entry: dict) -> dict:
    """Download per-residue pLDDT confidence JSON."""
    url = entry["plddtDocUrl"]
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    return resp.json()


def download_pae(entry: dict) -> dict | None:
    """
    Download PAE JSON.
    For complexes, the top-level payload is a list; unwrap it.
    Returns None if no PAE is available.
    """
    pae_url = entry.get("paeDocUrl")
    if not pae_url:
        return None
    resp = requests.get(pae_url, timeout=60)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    data = resp.json()
    # Complexes wrap the PAE dict in a single-element list
    if isinstance(data, list):
        data = data[0]
    return data


# ── Step 3: Visualize confidence ───────────────────────────────────────────────

def plot_confidence(entry: dict, conf: dict, pae_data: dict | None,
                    protein_name: str, gene: str, out_file: str = "surprise_me.png") -> None:
    """
    Generate a figure with:
      - Top: per-residue pLDDT bar chart (colour-coded by confidence tier)
      - Bottom: inter-chain PAE heatmap (if complex PAE is available)
    """
    is_complex = entry.get("isComplex", False)
    plddt_scores  = conf["confidenceScore"]
    residue_nums  = conf["residueNumber"]
    n             = len(plddt_scores)

    # Colour each residue by confidence category
    colours = []
    for s in plddt_scores:
        if s > 90:
            colours.append("#0053D6")   # very high — dark blue
        elif s > 70:
            colours.append("#65CBF3")   # high — light blue
        elif s > 50:
            colours.append("#FFDB13")   # low — yellow
        else:
            colours.append("#FF7D45")   # very low — orange

    fig = plt.figure(figsize=(14, 8 if pae_data else 5))
    title = f"{protein_name}"
    if gene:
        title += f"  [{gene}]"
    if is_complex:
        title += "  — homodimer"
    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)

    if pae_data:
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1.6], hspace=0.45)
    else:
        gs = gridspec.GridSpec(1, 1)

    # pLDDT bar chart
    ax1 = fig.add_subplot(gs[0])
    ax1.bar(residue_nums, plddt_scores, color=colours, width=1.0, linewidth=0)
    ax1.axhline(90, color="#0053D6", linewidth=0.7, linestyle="--", alpha=0.5)
    ax1.axhline(70, color="#65CBF3", linewidth=0.7, linestyle="--", alpha=0.5)
    ax1.axhline(50, color="#FFDB13", linewidth=0.7, linestyle="--", alpha=0.5)
    ax1.set_ylim(0, 100)
    ax1.set_xlabel("Residue position", fontsize=9)
    ax1.set_ylabel("pLDDT", fontsize=9)
    ax1.set_title("Per-residue confidence (pLDDT)", fontsize=10)

    # Chain boundary markers for complexes
    if is_complex and "chains" in conf:
        for chain in conf["chains"][:-1]:   # vertical line after each chain except the last
            ax1.axvline(chain["sequenceEnd"] + 0.5, color="black", linewidth=1.2, linestyle="-")

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#0053D6", label=">90 Very high"),
        Patch(facecolor="#65CBF3", label="70–90 High"),
        Patch(facecolor="#FFDB13", label="50–70 Low"),
        Patch(facecolor="#FF7D45", label="<50 Very low"),
    ]
    ax1.legend(handles=legend_elements, loc="lower right", fontsize=7, ncol=2)

    # PAE heatmap
    if pae_data:
        pae_matrix  = np.array(pae_data["predicted_aligned_error"])
        max_pae     = pae_data.get("max_predicted_aligned_error", 30)
        ax2         = fig.add_subplot(gs[1])
        im          = ax2.imshow(pae_matrix, cmap="viridis_r", vmin=0, vmax=max_pae,
                                 origin="upper", aspect="auto")
        plt.colorbar(im, ax=ax2, label="PAE (Å)", fraction=0.025, pad=0.02)
        ax2.set_xlabel("Scored residue", fontsize=9)
        ax2.set_ylabel("Aligned residue", fontsize=9)
        ax2.set_title("Predicted Aligned Error (PAE)", fontsize=10)

        # Draw chain boundary lines on PAE plot
        if "chains" in pae_data:
            for chain in pae_data["chains"][:-1]:
                boundary = chain["sequenceEnd"]
                ax2.axhline(boundary - 0.5, color="white", linewidth=1, linestyle="--", alpha=0.6)
                ax2.axvline(boundary - 0.5, color="white", linewidth=1, linestyle="--", alpha=0.6)

        # Summarize inter-chain confidence if complex
        if is_complex and "chains" in pae_data and len(pae_data["chains"]) >= 2:
            split = pae_data["chains"][0]["sequenceEnd"]
            inter = [pae_matrix[i][j]
                     for i in range(split)
                     for j in range(split, len(pae_matrix))]
            mean_inter = np.mean(inter)
            ax2.set_xlabel(
                f"Scored residue  (inter-chain mean PAE: {mean_inter:.1f} Å)", fontsize=9
            )

    plt.savefig(out_file, dpi=150, bbox_inches="tight")
    print(f"  Saved figure → {out_file}")


# ── Step 4: Search PubMed via NCBI Entrez ─────────────────────────────────────

def search_pubmed(query: str, max_results: int = MAX_PUBMED_RESULTS) -> list[dict]:
    """
    Search PubMed with Biopython Entrez.
    Returns a list of dicts with 'pmid', 'title', 'authors', 'journal', 'year', 'abstract'.
    """
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
    record = Entrez.read(handle)
    handle.close()

    pmids = record.get("IdList", [])
    if not pmids:
        return []

    time.sleep(0.4)   # NCBI rate-limit courtesy pause

    handle = Entrez.efetch(db="pubmed", id=",".join(pmids), rettype="xml", retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    articles = []
    for article in records.get("PubmedArticle", []):
        med  = article["MedlineCitation"]
        info = med["Article"]

        # PMID
        pmid = str(med.get("PMID", ""))

        # Title
        title = str(info.get("ArticleTitle", "No title"))

        # Authors
        author_list = info.get("AuthorList", [])
        authors = []
        for a in author_list[:3]:
            last  = a.get("LastName", "")
            init  = a.get("Initials", "")
            if last:
                authors.append(f"{last} {init}".strip())
        author_str = ", ".join(authors)
        if len(author_list) > 3:
            author_str += " et al."

        # Journal + year
        journal = str(info.get("Journal", {}).get("Title", ""))
        pub_date = info.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
        year = str(pub_date.get("Year", pub_date.get("MedlineDate", "")[:4]))

        # Abstract
        abstract_obj = info.get("Abstract", {}).get("AbstractText", [""])
        if isinstance(abstract_obj, list):
            abstract = " ".join(str(s) for s in abstract_obj)
        else:
            abstract = str(abstract_obj)

        articles.append({
            "pmid": pmid,
            "title": title,
            "authors": author_str,
            "journal": journal,
            "year": year,
            "abstract": abstract[:400] + "…" if len(abstract) > 400 else abstract,
        })

    return articles


def print_literature(articles: list[dict], protein_name: str, gene: str) -> None:
    query_label = gene if gene else protein_name
    print(f"\n{'─'*70}")
    print(f"  PubMed results for: {query_label}")
    print(f"{'─'*70}")
    if not articles:
        print("  No results found.")
        return
    for i, a in enumerate(articles, 1):
        print(f"\n  [{i}] {a['title']}")
        print(f"      {a['authors']}  |  {a['journal']} ({a['year']})")
        print(f"      PMID: {a['pmid']}  →  https://pubmed.ncbi.nlm.nih.gov/{a['pmid']}/")
        if a["abstract"]:
            wrapped = textwrap.fill(a["abstract"], width=66, initial_indent="      ",
                                    subsequent_indent="      ")
            print(wrapped)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print("AlphaFold Surprise Me Pipeline")
    print("=" * 70)

    # ── 1. Sample a random homodimer ──────────────────────────────────────────
    print("\n[1/4] Fetching homodimer pool from UniProt…")
    pool = fetch_homodimer_pool()
    if not pool:
        print("  ERROR: UniProt returned no results.")
        sys.exit(1)
    print(f"  Pool size: {len(pool)} reviewed human homodimers")

    entry_meta = pick_random_entry(pool)
    uniprot_id   = entry_meta["accession"]
    gene         = entry_meta["gene"]
    protein_name = entry_meta["protein"]
    print(f"\n  Selected: {protein_name}")
    print(f"  UniProt: {uniprot_id}  |  Gene: {gene or 'N/A'}")
    print(f"  UniProt page: https://www.uniprot.org/uniprotkb/{uniprot_id}")

    # ── 2. Fetch AlphaFold prediction ─────────────────────────────────────────
    print("\n[2/4] Fetching AlphaFold prediction…")
    af_entry = fetch_alphafold_entry(uniprot_id)

    if af_entry is None:
        print(f"  No AlphaFold entry found for {uniprot_id}. Re-run to try another protein.")
        sys.exit(0)

    is_complex = af_entry.get("isComplex", False)
    af_id      = af_entry.get("modelEntityId") or af_entry.get("entryId", "")
    tool_used  = af_entry.get("toolUsed", "AlphaFold")

    print(f"  AlphaFold ID:  {af_id}")
    print(f"  Is complex:    {is_complex}"
          + ("  (monomer entry — complex prediction not yet in DB)" if not is_complex else ""))
    print(f"  Tool used:     {tool_used}")
    print(f"  Structure URL: {af_entry.get('cifUrl', 'N/A')}")
    print(f"  AFDB page:     https://alphafold.ebi.ac.uk/entry/{uniprot_id}")

    # ── 3. Visualize structure ─────────────────────────────────────────────────
    print("\n[3/4] Downloading confidence data and generating figure…")
    conf = download_confidence(af_entry)
    pae_data = download_pae(af_entry)

    out_file = f"surprise_{uniprot_id}.png"
    plot_confidence(af_entry, conf, pae_data, protein_name, gene, out_file=out_file)

    # Print summary statistics
    scores = conf["confidenceScore"]
    avg_plddt  = np.mean(scores)
    frac_vh    = sum(1 for s in scores if s > 90) / len(scores)
    print(f"  Residues: {len(scores)}  |  Mean pLDDT: {avg_plddt:.1f}  |  >90: {frac_vh:.0%}")

    if is_complex and pae_data and "chains" in pae_data:
        pae_matrix = np.array(pae_data["predicted_aligned_error"])
        split      = pae_data["chains"][0]["sequenceEnd"]
        inter      = [pae_matrix[i][j]
                      for i in range(split)
                      for j in range(split, len(pae_matrix))]
        print(f"  Inter-chain mean PAE: {np.mean(inter):.1f} Å  "
              f"({'confident' if np.mean(inter) < 10 else 'moderate/uncertain'} interface)")

    # ── 4. Search PubMed ───────────────────────────────────────────────────────
    print("\n[4/4] Searching PubMed…")
    pubmed_query = gene if gene else protein_name
    # Try increasingly broad queries until we get results
    for query in [
        f"{pubmed_query} protein structure homodimer",
        f"{pubmed_query} protein structure",
        pubmed_query,
    ]:
        articles = search_pubmed(query)
        if articles:
            break
        time.sleep(0.3)
    print_literature(articles, protein_name, gene)

    print(f"\n{'='*70}")
    print("Done! Open the saved figure to explore the structure confidence.")
    print(f"  Figure: {out_file}")


if __name__ == "__main__":
    main()
