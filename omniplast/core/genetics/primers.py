"""Diseno y control de calidad de cebadores (primers) para qPCR/PCR.

Catalogo de genes funcionales de degradacion de plasticos en config.GENE_CATALOG.
REGLA DE HONESTIDAD: no se inventan secuencias de cebadores para genes funcionales
(son taxon-especificos y degenerados); el usuario los disena desde secuencias de
referencia con este modulo o consulta la fuente primaria.
"""
import numpy as np

from ...config import GENE_CATALOG, PRIMER_DISCLAIMER, UNIVERSAL_PRIMERS


def tm_wallace(seq):
    """Tm aproximada Wallace: 2(A+T) + 4(G+C) en °C (oligos cortos <14 nt)."""
    seq = seq.upper()
    return 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))


def tm_gc_salt(seq, na_conc_mM=50):
    """Tm por %GC ajustada por sal (formula simple para oligos 14-70 nt)."""
    seq = seq.upper()
    if len(seq) == 0:
        return np.nan
    gc = gc_content(seq) / 100
    tm = 64.9 + 41 * (len(seq) * gc - 16.4) / len(seq)
    return tm + 16.6 * np.log10(max(na_conc_mM, 1) / 1000)


def gc_content(seq):
    seq = seq.upper()
    if not seq:
        return np.nan
    return 100 * (seq.count("G") + seq.count("C")) / len(seq)


def reverse_complement(seq):
    return seq.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _max_complement_run(a, b):
    """Maximo tramo de complementariedad contigua entre dos secuencias."""
    rb = reverse_complement(b)
    best = 0
    for shift in range(-len(a) + 1, len(rb)):
        matched = 0
        for i in range(len(a)):
            j = i + shift
            if 0 <= j < len(rb) and a[i] == rb[j] and a[i] in "ACGT":
                matched += 1
                best = max(best, matched)
            else:
                matched = 0
    return best


def dimer_check(seq1, seq2=None):
    """Deteccion de dimeros de cebador: tramos complementarios >= 4 nt."""
    seq2 = seq2 or seq1
    runs = {"3p_end_1": _max_complement_run(seq1[-5:], seq2),
            "full_1": _max_complement_run(seq1, seq2)}
    risk = "alto" if runs["3p_end_1"] >= 5 else "medio" if runs["3p_end_1"] >= 4 else "bajo"
    return {"max_tramo_extremo_3p": runs["3p_end_1"], "max_tramo_completo": runs["full_1"],
            "riesgo_dimero": risk,
            "criterio": "Tramo complementario ≥5 nt en extremo 3' = riesgo alto"}


def hairpin_check(seq, min_stem=4):
    """Deteccion naive de horquillas: palindromo invertido con stem >= min_stem."""
    seq = seq.upper()
    found = []
    for loop in range(3, 9):
        for start in range(0, len(seq) - 2 * min_stem - loop + 1):
            stem1 = seq[start:start + min_stem]
            stem2_rc = reverse_complement(seq[start + min_stem + loop:start + 2 * min_stem + loop])
            if stem1 == stem2_rc:
                found.append({"posicion": start, "stem": stem1, "largo_loop": loop})
                break
        if found:
            break
    return {"horquilla_probable": bool(found), "detalle": found[:3],
            "nota": "Detección simplificada; confirmar con OligoAnalyzer/herramienta dedicada"}


def primer_penalty(seq, ideal_tm=60, ideal_gc=50):
    gc = gc_content(seq)
    tm = tm_gc_salt(seq)
    tm_pen = abs(tm - ideal_tm) if tm == tm else 99
    gc_pen = abs(gc - ideal_gc)
    end_gc_pen = 0 if seq.upper().endswith(("G", "C")) else 0.5
    return {"secuencia": seq, "longitud": len(seq), "tm": round(float(tm), 1),
            "gc_pct": round(gc, 1), "tm_wallace": tm_wallace(seq),
            "penalizacion": round(float(tm_pen * 0.6 + gc_pen * 0.4 + end_gc_pen), 2)}


def design_primers_from_sequence(sequence, product_min=80, product_max=250, tm_min=57,
                                 tm_max=63, max_candidates=10):
    """Propuesta simple de parejas forward/reverse desde una secuencia objetivo."""
    sequence = "".join(str(sequence).upper().split())
    results = []
    for f_start in range(0, len(sequence) - product_min):
        for f_len in range(18, 25):
            fwd = sequence[f_start:f_start + f_len]
            if len(fwd) < f_len:
                continue
            tm_f = tm_gc_salt(fwd)
            if not (tm_min <= tm_f <= tm_max) or not (45 <= gc_content(fwd) <= 60):
                continue
            r_end = f_start + product_max
            for r_len in range(18, 25):
                rev = reverse_complement(sequence[r_end - r_len:r_end])
                tm_r = tm_gc_salt(rev)
                if not (tm_min <= tm_r <= tm_max) or not (45 <= gc_content(rev) <= 60):
                    continue
                product = f_start + product_max - f_start
                actual_product = (r_end - r_len) - (f_start + f_len) + f_len + r_len
                dimer = dimer_check(fwd, rev)
                if dimer["riesgo_dimero"] == "alto":
                    continue
                results.append({
                    "forward": fwd, "reverse": rev,
                    "tm_forward": round(float(tm_f), 1), "tm_reverse": round(float(tm_r), 1),
                    "gc_f": round(gc_content(fwd), 1), "gc_r": round(gc_content(rev), 1),
                    "producto_estimado": int(r_end - f_start),
                    "penalizacion_total": round(primer_penalty(fwd)["penalizacion"] +
                                                primer_penalty(rev)["penalizacion"], 2)})
                break
            if len(results) >= max_candidates:
                break
        if len(results) >= max_candidates:
            break
    results.sort(key=lambda x: x["penalizacion_total"])
    return {"test": "Diseño de cebadores", "candidatos": results,
            "nota": PRIMER_DISCLAIMER,
            "interpretacion": f"{len(results)} parejas candidatas; validar experimentalmente "
                              "(gradiente de Tm, electroforesis, Sanger del amplicón)"}


def gene_panel_info():
    """Catalogo del panel de genes de plastisfera con referencias."""
    rows = []
    for key, g in GENE_CATALOG.items():
        rows.append({"clave": key, "gen": g["gen"], "organismo_modelo": g["organismo"],
                     "funcion": g["funcion"], "referencia": g["ref"]})
    return rows, PRIMER_DISCLAIMER


def universal_primers_table():
    return UNIVERSAL_PRIMERS
