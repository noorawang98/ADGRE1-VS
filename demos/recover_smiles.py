#!/usr/bin/env python3
"""Get trustworthy SMILES for the DrugBank/FDA docking ligands by matching against
/home/workspace/MD/MD/drugbank.sdf (which HAS bond orders).

Why this works where earlier attempts failed: PDBQT loses bond order, so obabel emits
radical aromatic rings ([C]1[C]=[C]...) and RDKit's 3D perception (rdDetermineBonds) does not
resolve valences for these OEChem files. Matching on a *neutralised single-bond skeleton*
canonical SMILES ignores bond order and charge entirely, so the radical-annotated PDBQT graph
still matches its proper SDF counterpart -- and the SDF then supplies a valid SMILES.

  ~/anaconda3/envs/AmberTools24/bin/python demo/ligand_smiles_match.py --limit 0
Outputs demo/ligand_smiles.tsv (ligand_id, smiles, hac, source, match)
"""
import argparse, os, subprocess, sys, tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SDF = "/home/workspace/MD/MD/drugbank.sdf"
FDA_DIR = os.path.expanduser("~/scripts/vinalc_test/data/fda")
DB_DIR = "/home/workspace/MD/MD/drugbank_mk_pdbqt"


def skeleton(mol):
    """canonical SMILES of the neutralised, all-single-bond graph -> bond-order free key"""
    from rdkit import Chem
    rw = Chem.RWMol(mol)
    for a in rw.GetAtoms():
        a.SetFormalCharge(0); a.SetIsAromatic(False)
        a.SetNumExplicitHs(0); a.SetNoImplicit(True)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
    try:
        return Chem.MolToSmiles(rw.GetMol())
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--out", default=os.path.join(ROOT, "demo", "ligand_smiles.tsv"))
    args = ap.parse_args()

    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    # ---- index the SDF by skeleton ----
    idx = {}
    supp = Chem.SDMolSupplier(SDF, sanitize=True, removeHs=True)
    n = 0
    for mol in supp:
        n += 1
        if mol is None:
            continue
        sk = skeleton(mol)
        if not sk:
            continue
        try:
            smi = Chem.MolToSmiles(mol)
        except Exception:
            continue
        idx.setdefault(sk, (smi, mol.GetNumHeavyAtoms()))
    print("SDF records: %d -> %d distinct skeletons" % (n, len(idx)))

    # ---- convert our PDBQT ligands and match ----
    jobs = []
    for d in (FDA_DIR, DB_DIR):
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            if f.endswith(".pdbqt"):
                jobs.append((f[:-6], os.path.join(d, f)))
    if args.limit:
        jobs = jobs[:args.limit]
    print("ligands to resolve: %d" % len(jobs))

    def work(job):
        lid, path = job
        tmp = tempfile.mktemp(suffix=".smi")
        try:
            subprocess.run(["obabel", path, "-osmi", "-O", tmp], capture_output=True, timeout=90)
            if not os.path.exists(tmp):
                return lid, "", 0, "obabel_fail", ""
            toks = open(tmp).readline().split()
            if not toks:
                return lid, "", 0, "obabel_empty", ""
            m = Chem.MolFromSmiles(toks[0], sanitize=False)
            if m is None:
                return lid, "", 0, "rdkit_fail", ""
            hac = m.GetNumAtoms()
            sk = skeleton(m)
            hit = idx.get(sk) if sk else None
            if hit:
                return lid, hit[0], hit[1], "sdf_match", sk
            return lid, "", hac, "no_match", sk or ""
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)

    with ThreadPoolExecutor(args.workers) as ex:
        rows = list(ex.map(work, jobs))
    c = Counter(r[3] for r in rows)
    print("match stats: %s" % dict(c))
    with open(args.out, "w") as fh:
        fh.write("ligand_id\tsmiles\thac\tsource\tmatch\n")
        for lid, smi, hac, src, _ in rows:
            fh.write("%s\t%s\t%d\t%s\t\n" % (lid, smi, hac, src))
    print("out -> %s" % os.path.relpath(args.out, ROOT))


if __name__ == "__main__":
    main()
