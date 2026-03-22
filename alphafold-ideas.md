## Entry-Level Project Ideas

These are approachable for students or early-career researchers with basic Python/bioinformatics skills:

### Structural Biology / Visualization
1. **"Before and after" gallery** — pick 20–50 proteins where the monomer looked disordered but the homodimer is well-folded. Visualize with PyMOL or Mol* and write up the biology. Great for a science communication piece.
2. **Interface residue predictor** — using the new structures, extract which residues form the dimer interface (e.g. via PISA or a simple distance cutoff). Train a simple ML classifier on sequence features to predict interface residues.

### Drug Discovery / Disease
3. **WHO pathogen interface mining** — for each WHO priority pathogen, identify the homodimer interfaces in the new data and cross-reference with known inhibitor binding sites (ChEMBL). Flag novel pocketable interfaces.
4. **Human disease gene complexes** — intersect the 1.7M homodimers with OMIM or ClinVar disease-associated genes. Build a table of disease-linked homodimers with no previously known structure.

### Comparative / Evolutionary
5. **Conservation of dimer interfaces** — pick a protein family (e.g. a kinase family) and compare the interface geometry across species in the 20 covered organisms. Are the interfaces more conserved than the rest of the surface?
6. **Taxonomic interface diversity** — cluster the 1.7M structures by interface shape (e.g. using contact maps) and ask: do closely related species have more similar interfaces?

### Data Science / ML
7. **Confidence score distribution analysis** — the 18M lower-confidence structures are now downloadable. What sequence/structural features correlate with low confidence? Build a predictor of AlphaFold complex confidence from sequence alone.
8. **Interface geometry clustering** — use dimensionality reduction (PCA/UMAP) on interface contact maps to find recurring "archetypes" of homodimer geometry. Does this recapitulate known structural families?

### Beginner-Friendly (no wet lab, just Python + public data)
9. **Validation project** — for proteins where an experimental homodimer structure exists in the PDB, compare the AlphaFold prediction vs. the crystal structure. Calculate RMSD, plot accuracy vs. pLDDT score.
