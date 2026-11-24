"""
Perfil longitudinal (rasante): PVI, tangentes y curvas verticales parabólicas.

Criterios (diapositivas "Diseño vertical"):
  * Pendiente máxima de tangente: tabla 4.1 (por categoría y VTR).  Pendiente mínima 0,3-0,5 %.
  * Longitud mínima de tangente vertical: tabla 4.3 (según velocidad específica).
  * Curva vertical simétrica parabólica; L = máx(criterio de seguridad por Dp, criterio de operación 0,6·V),
    redondeada hacia arriba a múltiplos de 10 m.
  * VECV = VECH donde la curva vertical coincide con una curva horizontal; en tangente horizontal, VETH
    (aquí: la mayor VECH de las curvas adyacentes)  [VERIFICAR contra la tabla 2.x del manual].

Método:
  1. Rasante inicial punto a punto por programación lineal (mínimo corte+relleno, |pendiente| <= máx).
  2. PVI iniciales por simplificación (Douglas-Peucker) de esa rasante.
  3. Ajuste de las cotas de los PVI por LP (pendientes entre mín y máx, con el signo ya definido),
     descontando el efecto de las curvas verticales de la iteración anterior.
  4. Longitud de cada curva por normas; si una curva no cabe (tangente mínima, bordes) se elimina el PVI
     más débil (menor A) y se repite hasta converger.
"""
from __future__ import annotations
import math
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack

from . import normas
from .alineamiento import douglas_peucker


# ----------------------------------------------------------------------------- terreno en el eje
def terreno_en_eje(terr, al, ds: float = 5.0):
    """Estaciones cada `ds` m (más la final), cota del TIN y coordenadas del eje."""
    S = al.muestrear(1.0)
    s_dense = S[:, 2] * (al.longitud / S[-1, 2])           # ajusta cuerdas a longitud de elementos
    s = np.unique(np.append(np.arange(0.0, al.longitud, ds), al.longitud))
    x = np.interp(s, s_dense, S[:, 0])
    y = np.interp(s, s_dense, S[:, 1])
    return s, terr.z(x, y), x, y


# ----------------------------------------------------------------------------- geometría vertical
def cota_diseno(s, s_pvi, z_pvi, curvas):
    """Cota de la rasante en las abscisas `s`: tangentes (interpolación entre PVI) + corrección parabólica.
    Dentro de una curva de longitud L y diferencia algebraica A (m/m): z = z_tangente + A·d²/(2L),
    con d = distancia al extremo más cercano de la curva (PCV o PTV)."""
    s = np.asarray(s, float)
    z = np.interp(s, s_pvi, z_pvi)
    for c in curvas:
        if c["L"] <= 0:
            continue
        x = s - (c["s_pvi"] - c["L"] / 2.0)
        m = (x >= 0) & (x <= c["L"])
        d = np.minimum(x[m], c["L"] - x[m])
        z[m] += c["A"] * d * d / (2.0 * c["L"])
    return z


def _velocidad_vertical(al, s_a, s_b, vtr):
    """VECV: si el tramo [s_a, s_b] coincide con curvas horizontales, la mayor VECH; si cae en tangente
    horizontal, VETH = mayor VECH de las curvas adyacentes; sin curvas, la VTR."""
    sobre = [c for c in al.curvas if not (c["s_fin"] < s_a or c["s_ini"] > s_b)]
    if sobre:
        return max(c["V"] for c in sobre)
    ant = [c for c in al.curvas if c["s_fin"] < s_a]
    sig = [c for c in al.curvas if c["s_ini"] > s_b]
    vs = ([ant[-1]["V"]] if ant else []) + ([sig[0]["V"]] if sig else [])
    return max(vs) if vs else float(vtr)


def _curvas(s_pvi, z_pvi, al, vtr, L_prev, paso):
    out = []
    for i in range(1, len(s_pvi) - 1):
        g1 = (z_pvi[i] - z_pvi[i - 1]) / (s_pvi[i] - s_pvi[i - 1])
        g2 = (z_pvi[i + 1] - z_pvi[i]) / (s_pvi[i + 1] - s_pvi[i])
        A = g2 - g1
        tipo = "convexa" if A < 0 else "concava"
        V = _velocidad_vertical(al, s_pvi[i] - L_prev[i] / 2, s_pvi[i] + L_prev[i] / 2, vtr)
        if abs(A) < 1e-9:
            out.append(dict(i=i, s_pvi=float(s_pvi[i]), z_pvi=float(z_pvi[i]), g1=g1, g2=g2, A=A, tipo=tipo,
                            V=normas.v_snap(V), L=0.0, Lmin=0.0, K=float("inf"), Kmin=0.0, Dp=0.0))
            continue
        r = normas.l_curva_vertical(abs(A) * 100.0, tipo, V, paso)
        out.append(dict(i=i, s_pvi=float(s_pvi[i]), z_pvi=float(z_pvi[i]), g1=g1, g2=g2, A=A, tipo=tipo,
                        V=r["V"], L=r["L"], Lmin=r["Lmin"], K=r["K"], Kmin=r["Kmin"], Dp=r["Dp"]))
    return out


def _conflictos(s_pvi, curvas, A_min=0.003):
    """Lista de (índice de PVI a eliminar, |A|, motivo)."""
    cand = []
    for c in curvas:
        if abs(c["A"]) < A_min:
            cand.append((c["i"], abs(c["A"]), f"deflexión vertical despreciable ({abs(c['A']) * 100:.2f} %)"))
        elif (c["s_pvi"] - c["L"] / 2 < s_pvi[0] - 1e-6) or (c["s_pvi"] + c["L"] / 2 > s_pvi[-1] + 1e-6):
            cand.append((c["i"], abs(c["A"]), "la curva se sale del perfil"))
    for a, b in zip(curvas[:-1], curvas[1:]):
        gap = (b["s_pvi"] - b["L"] / 2) - (a["s_pvi"] + a["L"] / 2)
        req = normas.ltv_min(max(a["V"], b["V"]))
        opuestas = a["tipo"] != b["tipo"]
        if not (gap >= req or (opuestas and gap >= 0)):     # convexa-cóncava contiguas: se permite sin tangente [DIAP]
            débil = a["i"] if abs(a["A"]) <= abs(b["A"]) else b["i"]
            cand.append((débil, min(abs(a["A"]), abs(b["A"])), f"tangente de {gap:.0f} m < {req:.0f} m (tabla 4.3)"))
    return cand


# ----------------------------------------------------------------------------- optimización (LP)
def _rasante_lp(s, zt, gmax, c_cut, c_fill):
    """Rasante punto a punto de mínimo costo de corte/relleno con |pendiente| <= gmax."""
    K = len(s)
    ds = np.diff(s)
    I = coo_matrix((np.ones(K), (np.arange(K), np.arange(K))), shape=(K, 3 * K)).tocsr()
    Iu = coo_matrix((np.ones(K), (np.arange(K), K + np.arange(K))), shape=(K, 3 * K)).tocsr()
    Iv = coo_matrix((-np.ones(K), (np.arange(K), 2 * K + np.arange(K))), shape=(K, 3 * K)).tocsr()
    A_eq = I + Iu + Iv
    r = np.arange(K - 1)
    D = coo_matrix((np.r_[-np.ones(K - 1), np.ones(K - 1)], (np.r_[r, r], np.r_[r, r + 1])),
                   shape=(K - 1, 3 * K)).tocsr()
    A_ub = vstack([D, -D]).tocsr()
    b_ub = np.r_[gmax * 0.999 * ds, gmax * 0.999 * ds]
    c = np.r_[np.zeros(K), np.full(K, c_cut), np.full(K, c_fill)]
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=zt,
                  bounds=[(None, None)] * K + [(0, None)] * (2 * K), method="highs")
    if not res.success:
        raise RuntimeError("Rasante inicial sin solución: " + res.message)
    return res.x[:K]


def _lp_pvis(s_pvi, s, zt, off, gmax, gmin, signs, c_cut, c_fill):
    """Cotas de los PVI (tangentes) que minimizan corte+relleno. `off` = efecto de curvas verticales previas."""
    m, K = len(s_pvi), len(s)
    seg = np.clip(np.searchsorted(s_pvi, s, side="right") - 1, 0, m - 2)
    t = (s - s_pvi[seg]) / (s_pvi[seg + 1] - s_pvi[seg])
    rows = np.repeat(np.arange(K), 4)
    cols = np.column_stack([seg, seg + 1, m + np.arange(K), m + K + np.arange(K)]).ravel()
    vals = np.column_stack([1 - t, t, np.ones(K), -np.ones(K)]).ravel()
    A_eq = coo_matrix((vals, (rows, cols)), shape=(K, m + 2 * K)).tocsr()
    b_eq = zt - off
    L = np.diff(s_pvi)
    lo = np.where(signs > 0, gmin, -gmax)
    hi = np.where(signs < 0, -gmin, gmax)
    r = np.arange(m - 1)
    D = coo_matrix((np.r_[-np.ones(m - 1), np.ones(m - 1)], (np.r_[r, r], np.r_[r, r + 1])),
                   shape=(m - 1, m + 2 * K)).tocsr()
    A_ub = vstack([D, -D]).tocsr()
    b_ub = np.r_[hi * L, -lo * L]
    c = np.r_[np.zeros(m), np.full(K, c_cut), np.full(K, c_fill)]
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=[(None, None)] * m + [(0, None)] * (2 * K), method="highs")
    if not res.success:
        raise RuntimeError("Ajuste de PVI sin solución: " + res.message)
    return res.x[:m]


# ----------------------------------------------------------------------------- resultado
@dataclass
class Perfil:
    s: np.ndarray
    zt: np.ndarray
    zd: np.ndarray
    s_pvi: np.ndarray
    z_pvi: np.ndarray
    curvas: list
    gmax: float
    gmin: float
    vtr: float
    violaciones: list = field(default_factory=list)
    avisos: list = field(default_factory=list)
    log: list = field(default_factory=list)
    stats: dict = field(default_factory=dict)
    modo: str = "eje"          # "eje": mínimo corte/relleno en el eje; "objetivo": ajustado a una rasante objetivo

    def pvis(self):
        """Lista [(abscisa, cota, L_curva)] para LandXML (extremos con L=0)."""
        L = {c["i"]: c["L"] for c in self.curvas}
        return [(float(a), float(b), float(L.get(i, 0.0))) for i, (a, b) in enumerate(zip(self.s_pvi, self.z_pvi))]

    def tabla(self):
        filas = []
        for i, (a, b) in enumerate(zip(self.s_pvi, self.z_pvi)):
            c = next((c for c in self.curvas if c["i"] == i), None)
            filas.append(dict(PVI=i, abscisa=a, cota=b,
                              pend_entrada_pct=100 * c["g1"] if c else float("nan"),
                              pend_salida_pct=100 * c["g2"] if c else 100 * (self.z_pvi[i + 1] - b) / (self.s_pvi[i + 1] - a) if i < len(self.s_pvi) - 1 else float("nan"),
                              A_pct=100 * abs(c["A"]) if c else 0.0, tipo=c["tipo"] if c else "-", V_kmh=c["V"] if c else float("nan"),
                              L_m=c["L"] if c else 0.0, Lmin_m=c["Lmin"] if c else 0.0,
                              K=c["K"] if c and c["L"] > 0 else float("nan")))
        return filas


def disenar_perfil(terr, al, categoria, vtr, ds=5.0, tol_dp=1.0, c_cut=1.0, c_fill=1.4,
                   gmax=None, gmin=normas.PEND_MIN, paso_L=10.0, max_iter=80, A_min=0.003, z_ref=None) -> Perfil:
    """z_ref=(abscisas, cotas): rasante OBJETIVO (p. ej. la del optimizador por volúmenes). Si se da, los PVI se
    ajustan a ella (L1, costos c_cut/c_fill) en lugar de seguir el terreno del eje; el perfil sigue cumpliendo
    pendientes, tangentes mínimas y curvas verticales, y las estadísticas se miden contra el terreno real."""
    if not getattr(al, "valido", True):
        raise ValueError("El alineamiento horizontal es inválido (curvas solapadas); no se puede diseñar el perfil: "
                         + "; ".join(al.solapes[:3]))
    s, zt, _, _ = terreno_en_eje(terr, al, ds)
    gmax = gmax if gmax is not None else normas.pend_media_max(categoria, vtr)
    zfit = zt if z_ref is None else np.interp(s, np.asarray(z_ref[0], float), np.asarray(z_ref[1], float))
    z_lp = _rasante_lp(s, zfit, gmax, c_cut, c_fill)
    idx = douglas_peucker(np.column_stack([s, z_lp]), tol_dp)
    s_pvi = s[idx].astype(float)
    signs = np.zeros(len(s_pvi) - 1)
    off = np.zeros(len(s))
    L_prev = np.full(len(s_pvi), 100.0)
    log = [f"PVI iniciales por simplificación (tol {tol_dp} m): {len(s_pvi)}"]
    z_pvi, curvas = None, []
    for it in range(max_iter):
        z_pvi = _lp_pvis(s_pvi, s, zfit, off, gmax, gmin, signs, c_cut, c_fill)
        curvas = _curvas(s_pvi, z_pvi, al, vtr, L_prev, paso_L)
        conf = _conflictos(s_pvi, curvas, A_min)
        if conf:
            k, _, motivo = min(conf, key=lambda t: t[1])
            log.append(f"iter {it}: se elimina PVI en K{s_pvi[k] / 1000:.3f} ({motivo})")
            s_pvi = np.delete(s_pvi, k)
            signs = np.zeros(len(s_pvi) - 1)
            off = np.zeros(len(s))
            L_prev = np.full(len(s_pvi), 100.0)
            continue
        g = np.diff(z_pvi) / np.diff(s_pvi)
        tend = np.diff(np.interp(s_pvi, s, zfit))
        nuevos = np.where(np.abs(g) >= gmin * 0.999, np.sign(g), np.sign(tend))
        nuevos[nuevos == 0] = 1.0
        off_new = cota_diseno(s, s_pvi, z_pvi, curvas) - np.interp(s, s_pvi, z_pvi)
        conv = (np.max(np.abs(off_new - off)) < 0.02) and np.array_equal(nuevos, signs)
        off, signs = off_new, nuevos
        L_prev = np.array([0.0] + [c["L"] for c in curvas] + [0.0])
        if conv:
            log.append(f"convergió en la iteración {it}")
            break
    else:
        log.append(f"AVISO: no convergió en {max_iter} iteraciones")

    zd = cota_diseno(s, s_pvi, z_pvi, curvas)
    pf = Perfil(s=s, zt=zt, zd=zd, s_pvi=s_pvi, z_pvi=z_pvi, curvas=curvas, gmax=gmax, gmin=gmin, vtr=vtr, log=log)
    pf.modo = "eje" if z_ref is None else "objetivo"
    if z_ref is None:
        pf._area_lp = float(np.trapezoid(np.abs(zt - z_lp), s))   # cota inferior: rasante libre punto a punto
    _verificar(pf, A_min)
    return pf


# ----------------------------------------------------------------------------- verificación y estadísticas
def _verificar(pf: Perfil, A_min):
    tol = 1e-6
    g = np.diff(pf.z_pvi) / np.diff(pf.s_pvi)
    for i, gi in enumerate(g):
        if abs(gi) > pf.gmax + tol:
            pf.violaciones.append(f"Tangente {i}-{i + 1}: pendiente {gi * 100:.2f} % > máx {pf.gmax * 100:.1f} %")
        if abs(gi) < pf.gmin - tol:
            pf.violaciones.append(f"Tangente {i}-{i + 1}: pendiente {abs(gi) * 100:.2f} % < mín {pf.gmin * 100:.1f} %")
    for _, _, motivo in _conflictos(pf.s_pvi, pf.curvas, A_min):
        pf.violaciones.append("Curvas verticales: " + motivo)
    # avisos: longitud crítica por tramos continuos de igual signo con pendiente >= 3 %
    i = 0
    while i < len(g):
        if abs(g[i]) >= 0.03:
            j = i
            while j + 1 < len(g) and np.sign(g[j + 1]) == np.sign(g[i]) and abs(g[j + 1]) >= 0.03:
                j += 1
            largo = pf.s_pvi[j + 1] - pf.s_pvi[i]
            media = (pf.z_pvi[j + 1] - pf.z_pvi[i]) / largo
            lc = normas.l_critica(media)
            if largo > lc:
                sentido = "ascenso" if media > 0 else "ascenso en sentido contrario"
                pf.avisos.append(f"K{pf.s_pvi[i] / 1000:.3f}-K{pf.s_pvi[j + 1] / 1000:.3f}: {largo:.0f} m a {abs(media) * 100:.1f} % supera la "
                                 f"longitud crítica aprox. ({lc:.0f} m, {sentido}); evaluar carril de ascenso")
            i = j + 1
        else:
            i += 1
    neta = abs(pf.z_pvi[-1] - pf.z_pvi[0]) / (pf.s_pvi[-1] - pf.s_pvi[0])
    if neta > pf.gmax:
        pf.avisos.append(f"Pendiente media del corredor {neta * 100:.1f} % > {pf.gmax * 100:.1f} %")
    dz = pf.zt - pf.zd
    corte = np.where(dz > 0, dz, 0.0)
    relleno = np.where(dz < 0, -dz, 0.0)
    pf.stats = dict(area_total_cota_inferior_m2=getattr(pf, "_area_lp", float("nan")), n_pvi=len(pf.s_pvi), n_curvas=sum(1 for c in pf.curvas if c["L"] > 0),
                    pend_media_corredor_pct=100 * neta,
                    pend_max_usada_pct=100 * float(np.max(np.abs(g))) if len(g) else 0.0,
                    area_corte_m2=float(np.trapezoid(corte, pf.s)), area_relleno_m2=float(np.trapezoid(relleno, pf.s)),
                    corte_max_m=float(corte.max()), relleno_max_m=float(relleno.max()),
                    corte_medio_m=float(corte.mean()), relleno_medio_m=float(relleno.mean()))


def graficar(pf: Perfil, ruta: str):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(13, 8), gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
    ax.plot(pf.s, pf.zt, color="green", lw=1.5, label="Terreno en el eje")
    ax.plot(pf.s, pf.zd, color="red", lw=2.2, label="Rasante (con curvas verticales)")
    ax.fill_between(pf.s, pf.zt, pf.zd, where=pf.zt > pf.zd, color="red", alpha=0.2, label="Corte")
    ax.fill_between(pf.s, pf.zt, pf.zd, where=pf.zt <= pf.zd, color="blue", alpha=0.2, label="Relleno")
    ax.plot(pf.s_pvi, pf.z_pvi, "k^", ms=6, label="PVI")
    for c in pf.curvas:
        if c["L"] > 0:
            ax.axvspan(c["s_pvi"] - c["L"] / 2, c["s_pvi"] + c["L"] / 2, color="orange", alpha=0.15)
    for i in range(len(pf.s_pvi) - 1):
        gi = (pf.z_pvi[i + 1] - pf.z_pvi[i]) / (pf.s_pvi[i + 1] - pf.s_pvi[i])
        xm = (pf.s_pvi[i] + pf.s_pvi[i + 1]) / 2
        ax.annotate(f"{gi * 100:+.2f} %", (xm, np.interp(xm, pf.s_pvi, pf.z_pvi)), textcoords="offset points",
                    xytext=(0, 12), ha="center", fontsize=8, color="darkred")
    ax.set_ylabel("Cota (m)"); ax.legend(loc="upper left", fontsize=8); ax.grid(alpha=0.3)
    ax.set_title(f"Perfil longitudinal - VTR {pf.vtr:.0f} km/h - pendiente máx. {pf.gmax * 100:.0f} %")
    ax2.fill_between(pf.s, 0, pf.zt - pf.zd, where=pf.zt >= pf.zd, color="red", alpha=0.5)
    ax2.fill_between(pf.s, 0, pf.zt - pf.zd, where=pf.zt < pf.zd, color="blue", alpha=0.5)
    ax2.axhline(0, color="k", lw=0.6); ax2.set_ylabel("Corte (+) / Relleno (-) en eje (m)"); ax2.set_xlabel("Abscisa (m)")
    ax2.grid(alpha=0.3)
    plt.tight_layout(); fig.savefig(ruta, dpi=100); plt.close(fig)
