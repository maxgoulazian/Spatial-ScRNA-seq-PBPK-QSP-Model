#!/usr/bin/env python3
"""
Renumber an antibody PDB to Chothia using precomputed ANARCI/ANARCII CSV output.

Assumptions
-----------
- Input PDB has chain A = original light chain, chain B = original heavy chain.
  (Override with --light-chain-id / --heavy-chain-id.)
- CSV columns after metadata (Name, Chain, Score, Query start, Query end, ...) are
  positional labels: numeric or numeric with insertion letters (e.g. "52A").
- Cell values are one-letter amino acids; "-" = gap / no residue.
- We map non-gap positions in CSV, in column order, to residues in the
  corresponding PDB chain encountered in file order.

Outputs
-------
- Renumbered PDB with chain IDs L (light) and H (heavy).
- Mapping report TSV (old_chain old_resseq old_icode old_resname -> new_chain new_resseq new_icode chothia_label aa_csv aa_pdb match_bool).

Author: ChatGPT (generated for user workflow)
"""

import argparse
import csv
import sys
from collections import namedtuple
from typing import List, Tuple, Dict, Optional
import pandas as pd
import re
import pathlib

# --- AA maps ---
AA3_TO_1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    "SEC":"U","PYL":"O",  # rarely used
}
AA1_TO_3 = {v:k for k,v in AA3_TO_1.items()}
AA1_TO_3.setdefault("X","UNK")
AA3_TO_1.setdefault("MSE","M")  # selenomethionine -> M


# ----------------- CSV PARSING -----------------
META_COLS = ["Name","Chain","Score","Query start","Query end"]

def _col_sort_key(label: str) -> Tuple[int,str]:
    """
    Sort positional column labels like '52','52A','82B','113','113AA'.
    Returns (int_part, insertion_str) where insertion_str sorts lexicographically.
    """
    m = re.match(r"^(\d+)([A-Za-z]*)$", label)
    if not m:
        # Unexpected; push to far right
        return (10**9, label)
    num = int(m.group(1))
    ins = m.group(2) or ""
    return (num, ins)

def load_anarci_csv(path: str,
                    light_row_idx: Optional[int]=None,
                    heavy_row_idx: Optional[int]=None) -> Dict[str, List[Tuple[str,str,str]]]:
    """
    Load ANARCI/ANARCII CSV and return dict {'L':[(label,num,aa1),...], 'H':[...]}.
    label is e.g. '52A'; num is numeric part as str; aa1 single letter.
    """
    df = pd.read_csv(path, dtype=str)  # keep as strings
    # Normalize column names (strip spaces)
    df.columns = [c.strip() for c in df.columns]

    # Determine row indices
    if light_row_idx is None or heavy_row_idx is None:
        # attempt by Chain column
        lcandidates = df.index[df["Chain"].str.upper().str.startswith(("K","L"))].tolist()
        hcandidates = df.index[df["Chain"].str.upper().str.startswith("H")].tolist()
        if light_row_idx is None:
            light_row_idx = lcandidates[0] if lcandidates else 0
        if heavy_row_idx is None:
            heavy_row_idx = hcandidates[0] if hcandidates else (1 if df.shape[0]>1 else 0)

    light_row = df.iloc[light_row_idx]
    heavy_row = df.iloc[heavy_row_idx]

    # Positional columns
    pos_cols = [c for c in df.columns if c not in META_COLS]

    # Sort them by numeric/insertion order
    pos_cols.sort(key=_col_sort_key)

    def build_list(row) -> List[Tuple[str,str,str]]:
        out = []
        for c in pos_cols:
            aa = str(row[c]).strip() if c in row else "-"
            if aa in ("", "nan", "NaN"):
                aa = "-"
            if aa == "-":
                continue
            # accept only first character
            aa1 = aa[0].upper()
            # parse number + insertion
            m = re.match(r"^(\d+)([A-Za-z]*)$", c)
            if not m:
                continue
            num = m.group(1)
            ins = m.group(2)  # may be ''
            out.append((c, num, ins, aa1))
        return out

    L_list = build_list(light_row)
    H_list = build_list(heavy_row)

    return {"L": L_list, "H": H_list}


# ----------------- PDB PARSING / WRITING -----------------
PDBAtom = namedtuple("PDBAtom",
                     "raw line_idx record serial name alt loc chain resname resseq icode x y z occ bfac element charge")

def parse_pdb_atoms(path: str) -> List[PDBAtom]:
    atoms = []
    with open(path) as fh:
        for i,line in enumerate(fh):
            rec = line[0:6].strip()
            if rec not in ("ATOM","HETATM"):
                continue
            serial = line[6:11]
            name   = line[12:16]
            alt    = line[16:17]
            resname= line[17:20]
            loc    = line[20:21]
            chain  = line[21:22]
            resseq = line[22:26]
            icode  = line[26:27]
            x      = line[30:38]
            y      = line[38:46]
            z      = line[46:54]
            occ    = line[54:60]
            bfac   = line[60:66]
            element= line[76:78] if len(line)>=78 else "  "
            charge = line[78:80] if len(line)>=80 else "  "
            atoms.append(PDBAtom(line, i, rec, serial, name, alt, loc, chain, resname,
                                 resseq, icode, x,y,z, occ, bfac, element, charge))
    return atoms

def group_residues(atoms: List[PDBAtom], chain_id: str) -> List[Tuple[str,str,str,List[PDBAtom]]]:
    """
    Return list of residues for given chain in file order.
    Each entry: (resname, resseq_str, icode_char, atom_list_for_residue)
    """
    reslist = []
    curkey = None
    curatoms = []
    for a in atoms:
        if a.chain != chain_id:
            continue
        key = (a.resname, a.resseq.strip(), a.icode.strip())
        if curkey is None:
            curkey = key
        if key != curkey:
            reslist.append((*curkey, curatoms))
            curatoms = []
            curkey = key
        curatoms.append(a)
    if curkey is not None:
        reslist.append((*curkey, curatoms))
    return reslist

def aa3_to_1(resname: str) -> str:
    return AA3_TO_1.get(resname.upper().strip(), "X")

def format_pdb_atom(a: PDBAtom,
                    new_chain: str,
                    new_resseq: int,
                    new_icode: str) -> str:
    """
    Return reformatted PDB atom line with updated chain, resseq, insertion code.
    """
    # PDB formatting fields
    # resseq width 4 right-aligned; icode single char
    resseq_str = f"{new_resseq:>4d}"
    icode_char = (new_icode or " ")[0]
    newline = (
        f"{a.record:<6s}"
        f"{int(a.serial):5d} "
        f"{a.name:<4s}"
        f"{a.alt:<1s}"
        f"{a.resname:>3s} "
        f"{new_chain:1s}"
        f"{resseq_str}{icode_char}"
        f"   "  # gap for segid in col 28-30
        f"{float(a.x):8.3f}{float(a.y):8.3f}{float(a.z):8.3f}"
        f"{float(a.occ):6.2f}{float(a.bfac):6.2f}          "
        f"{a.element:>2s}{a.charge:>2s}\n"
    )
    # raw numeric parse robust fallback: if x etc not numeric, fall back to original coords
    try:
        float(a.x); float(a.y); float(a.z)
    except ValueError:
        newline = (
            a.raw[:21] + new_chain +
            resseq_str + icode_char +
            a.raw[27:]  # remainder
        )
    return newline


# ----------------- MAPPING / RENUMB -----------------
def map_chain(pdb_reslist, csv_list, chain_label, ignore_aa_mismatch=False):
    """
    Map residues: pdb_reslist list of (resname, old_resseq, old_icode, atoms)
    csv_list list of tuples (label, num, ins, aa1)
    Returns list of mapping records for the min(len(pdb_reslist), len(csv_list)).
    """
    n = min(len(pdb_reslist), len(csv_list))
    mapped = []
    warnings = []
    for i in range(n):
        resname, old_resseq, old_icode, atoms = pdb_reslist[i]
        label, num, ins, aa_csv = csv_list[i]
        aa_pdb = aa3_to_1(resname)
        match = (aa_pdb == aa_csv)
        if not match and not ignore_aa_mismatch:
            warnings.append(
                f"[{chain_label}] AA mismatch at {i}: PDB {aa_pdb} vs CSV {aa_csv} (label {label})"
            )
        # Choose insertion char
        icode = ins[-1] if ins else ""
        mapped.append({
            "old_chain": atoms[0].chain,
            "old_resseq": old_resseq,
            "old_icode": old_icode,
            "old_resname": resname,
            "new_chain": chain_label,
            "new_resseq": int(num),
            "new_icode": icode,
            "chothia_label": label,
            "aa_csv": aa_csv,
            "aa_pdb": aa_pdb,
            "match": match
        })
    # trailing residues in pdb_reslist beyond csv_list
    extra = pdb_reslist[n:]
    if extra:
        warnings.append(
            f"[{chain_label}] {len(extra)} PDB residues have no Chothia mapping (trailing)."
        )
    # trailing positions in csv_list beyond pdb_reslist
    extra_csv = csv_list[n:]
    if extra_csv:
        warnings.append(
            f"[{chain_label}] CSV has {len(extra_csv)} positions beyond PDB length."
        )
    return mapped, extra, warnings

def renumber_pdb(pdb_path, csv_path, out_path,
                 report_path=None,
                 light_chain_id="A",
                 heavy_chain_id="B",
                 ignore_aa_mismatch=False,
                 crop_to_mapped=False,
                 light_row_idx=None,
                 heavy_row_idx=None):
    atoms = parse_pdb_atoms(pdb_path)
    pdb_L = group_residues(atoms, light_chain_id)
    pdb_H = group_residues(atoms, heavy_chain_id)
    csv_data = load_anarci_csv(csv_path,
                               light_row_idx=light_row_idx,
                               heavy_row_idx=heavy_row_idx)
    map_L, extra_L, warn_L = map_chain(pdb_L, csv_data["L"], "L",
                                       ignore_aa_mismatch=ignore_aa_mismatch)
    map_H, extra_H, warn_H = map_chain(pdb_H, csv_data["H"], "H",
                                       ignore_aa_mismatch=ignore_aa_mismatch)

    warnings = warn_L + warn_H

    # For extras: assign sequential numbers continuing from last mapped if not cropping
    def continue_seq(start_resseq, extra_reslist, chain_label):
        maps = []
        num = start_resseq
        for resname, old_resseq, old_icode, atoms in extra_reslist:
            num += 1
            maps.append({
                "old_chain": atoms[0].chain,
                "old_resseq": old_resseq,
                "old_icode": old_icode,
                "old_resname": resname,
                "new_chain": chain_label,
                "new_resseq": num,
                "new_icode": "",
                "chothia_label": f"{num}",
                "aa_csv": "",
                "aa_pdb": aa3_to_1(resname),
                "match": True
            })
        return maps

    if not crop_to_mapped:
        if extra_L:
            last = map_L[-1]["new_resseq"] if map_L else 0
            map_L += continue_seq(last, extra_L, "L")
        if extra_H:
            last = map_H[-1]["new_resseq"] if map_H else 0
            map_H += continue_seq(last, extra_H, "H")

    full_map = map_L + map_H
    # Build lookup: (old_chain, old_resseq, old_icode) -> (new_chain,new_resseq,new_icode)
    lookup = {}
    for m in full_map:
        key = (m["old_chain"],
               str(m["old_resseq"]).strip(),
               str(m["old_icode"]).strip())
        lookup[key] = (m["new_chain"], m["new_resseq"], m["new_icode"])

    # Rewrite PDB
    out_lines = []
    with open(pdb_path) as fh:
        for line in fh:
            rec = line[0:6].strip()
            if rec in ("ATOM","HETATM"):
                chain = line[21:22]
                resseq = line[22:26].strip()
                icode = line[26:27].strip()
                key = (chain, resseq, icode)
                if key in lookup:
                    new_chain, new_resseq, new_icode = lookup[key]
                    # parse numeric safe
                    try:
                        serial = int(line[6:11])
                    except ValueError:
                        serial = 0
                    # reformat via our PDBAtom parse (we stored earlier but simpler here)
                    a = PDBAtom(line, -1, rec, line[6:11], line[12:16], line[16:17],
                                line[20:21], chain, line[17:20],
                                resseq, icode,
                                line[30:38], line[38:46], line[46:54],
                                line[54:60], line[60:66],
                                line[76:78] if len(line)>=78 else "  ",
                                line[78:80] if len(line)>=80 else "  ")
                    out_lines.append(format_pdb_atom(a, new_chain, new_resseq, new_icode))
                else:
                    if crop_to_mapped:
                        # skip this atom
                        continue
                    # leave line unchanged but optionally rename chain?
                    if chain == light_chain_id:
                        new_chain = "L"
                    elif chain == heavy_chain_id:
                        new_chain = "H"
                    else:
                        new_chain = chain
                    # don't change residue numbering
                    a = PDBAtom(line, -1, rec, line[6:11], line[12:16], line[16:17],
                                line[20:21], chain, line[17:20],
                                resseq, icode,
                                line[30:38], line[38:46], line[46:54],
                                line[54:60], line[60:66],
                                line[76:78] if len(line)>=78 else "  ",
                                line[78:80] if len(line)>=80 else "  ")
                    out_lines.append(format_pdb_atom(a, new_chain, int(resseq or 0), icode))
            else:
                out_lines.append(line)

    with open(out_path, "w") as out_fh:
        out_fh.writelines(out_lines)

    if report_path:
        cols = ["old_chain","old_resseq","old_icode","old_resname",
                "new_chain","new_resseq","new_icode","chothia_label",
                "aa_csv","aa_pdb","match"]
        pd.DataFrame(full_map)[cols].to_csv(report_path, sep="\t", index=False)

    return warnings


# ----------------- CLI -----------------
def main():
    ap = argparse.ArgumentParser(description="Renumber antibody PDB to Chothia using ANARCI CSV.")
    ap.add_argument("--pdb", required=True, help="Input PDB file (chain A=LC, B=HC unless overridden).")
    ap.add_argument("--csv", required=True, help="ANARCI CSV file.")
    ap.add_argument("--out", required=True, help="Output renumbered PDB.")
    ap.add_argument("--report", help="Write mapping report TSV.")
    ap.add_argument("--light-chain-id", default="A", help="Original PDB light-chain ID.")
    ap.add_argument("--heavy-chain-id", default="B", help="Original PDB heavy-chain ID.")
    ap.add_argument("--light-row-idx", type=int, help="Row index in CSV to use as light chain.")
    ap.add_argument("--heavy-row-idx", type=int, help="Row index in CSV to use as heavy chain.")
    ap.add_argument("--ignore-aa-mismatch", action="store_true", help="Do not warn/fail on AA mismatches.")
    ap.add_argument("--crop-to-mapped", action="store_true", help="Drop atoms not mapped to Chothia positions.")
    args = ap.parse_args()

    warns = renumber_pdb(
        pdb_path=args.pdb,
        csv_path=args.csv,
        out_path=args.out,
        report_path=args.report,
        light_chain_id=args.light_chain_id,
        heavy_chain_id=args.heavy_chain_id,
        ignore_aa_mismatch=args.ignore_aa_mismatch,
        crop_to_mapped=args.crop_to_mapped,
        light_row_idx=args.light_row_idx,
        heavy_row_idx=args.heavy_row_idx,
    )

    if warns:
        sys.stderr.write("\n".join(["WARN: "+w for w in warns]) + "\n")

if __name__ == "__main__":
    main()

