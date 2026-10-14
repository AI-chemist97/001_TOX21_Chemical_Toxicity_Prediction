# Python
"""Check whether chemistry-based substructure groups have a higher active rate in Tox21.

The four hypotheses below were written before looking at the numbers:
    H1  quinone / hydroquinone / catechol -> SR-ARE   (redox cycling gives oxidative stress)
    H2  organotin                         -> SR-MMP   (organotins damage the mitochondrial membrane)
    H3  3+ fused aromatic rings           -> NR-AhR   (flat aromatics bind the AhR receptor)
    H4  nitroaromatic                     -> SR-ARE   (electrophilic stress)

Usage:
    python domain_analysis.py        # writes domain_enrichment.csv
"""
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
HERE = Path(__file__).parent
TARGETS = ["NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase", "NR-ER", "NR-ER-LBD",
           "NR-PPAR-gamma", "SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53"]

# group name -> SMARTS patterns (a molecule is in the group if any pattern matches)
SMARTS = {
    "hydroquinone_or_catechol": ["[OX2H]c1ccc([OX2H])cc1", "[OX2H]c:c[OX2H]"],
    "quinone": ["[#6]1(=O)[#6]=,:[#6][#6](=O)[#6]=,:[#6]1", "[#6]1(=O)[#6](=O)[#6]=,:[#6][#6]=,:[#6]1"],
    "nitroaromatic": ["c[N+](=O)[O-]", "cN(=O)=O"],
    "phenol": ["c[OX2H]"],
    "halogenated_aromatic": ["c[Cl,Br,I]"],
    "organotin": ["[Sn;$([Sn]-[#6])]"],
    "heavy_metal": ["[Hg,Cd,Pb,As,Cu,Ni,Co,Cr,Zn,Sb,Tl,Se]"],
    "quaternary_ammonium": ["[NX4+]"],
    "long_alkyl_chain": ["[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]"],
    "carboxylic_acid": ["C(=O)[OX2H,OX1-]"],
}
PATTERNS = {k: [Chem.MolFromSmarts(s) for s in v] for k, v in SMARTS.items()}


def max_fused_aromatic_rings(mol):
    """Number of rings in the biggest group of aromatic rings that share 2+ atoms (fused rings)."""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()
             if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    best, seen = 0, set()
    for start in range(len(rings)):
        if start in seen:
            continue
        group, todo = {start}, [start]
        while todo:
            cur = todo.pop()
            for j in range(len(rings)):
                if j not in group and len(rings[cur] & rings[j]) >= 2:
                    group.add(j)
                    todo.append(j)
        seen |= group
        best = max(best, len(group))
    return best


def main():
    if not (HERE / "tox21.csv").exists():
        import tox21_baseline
        tox21_baseline.load_data(None)      # downloads tox21.csv next to this file
    df = pd.read_csv(HERE / "tox21.csv")
    rows = []
    for smi in df["smiles"]:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rows.append(None)
            continue
        flag = {k: any(mol.HasSubstructMatch(p) for p in v if p is not None) for k, v in PATTERNS.items()}
        flag["polycyclic_aromatic_3plus"] = max_fused_aromatic_rings(mol) >= 3
        rows.append(flag)

    ok = [i for i, r in enumerate(rows) if r is not None]
    flags = pd.DataFrame([rows[i] for i in ok], index=df.index[ok])
    labels = df.loc[ok, TARGETS]
    print("parsed", len(ok), "of", len(df))

    out = []
    for group in flags.columns:
        in_group = flags[group]
        for t in TARGETS:
            y = labels[t]
            tested = y.notna()                      # missing labels are left out
            grp = tested & in_group
            rest = tested & ~in_group
            if grp.sum() == 0:
                continue
            rate_grp = y[grp].mean()
            rate_rest = y[rest].mean()
            out.append({"group": group, "target": t, "n_group": int(grp.sum()),
                        "active_in_group": int(y[grp].sum()), "rate_group": round(rate_grp, 3),
                        "rate_rest": round(rate_rest, 3), "fold": round(rate_grp / rate_rest, 2)})
    out = pd.DataFrame(out)
    out.to_csv(HERE / "domain_enrichment.csv", index=False)

    print("group sizes:", {g: int(flags[g].sum()) for g in flags.columns})
    print("\nThe four hypotheses")
    for g, t in [("quinone", "SR-ARE"), ("hydroquinone_or_catechol", "SR-ARE"), ("organotin", "SR-MMP"),
                 ("polycyclic_aromatic_3plus", "NR-AhR"), ("nitroaromatic", "SR-ARE")]:
        print(out[(out.group == g) & (out.target == t)].to_string(index=False, header=False))
    print("\nLargest differences (groups with 15+ compounds)")
    big = out[out.n_group >= 15].sort_values("fold", ascending=False).head(15)
    print(big.to_string(index=False))


if __name__ == "__main__":
    main()
