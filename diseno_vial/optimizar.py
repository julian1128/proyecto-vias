"""
Rasante optimizada por volúmenes (balance de masas).

La rasante base (perfil.py) minimiza corte y relleno solo EN EL EJE. En ladera eso deja secciones con mucho corte
y poco terraplén (excedente de material). Aquí se optimiza con el costo de la SECCIÓN COMPLETA:

  costo(rasante) = c_corte·Vc + c_relleno·Vf + acarreo + botadero(excedente) + préstamo(faltante)

Método:
  1. Para cada estación y cada cota candidata de rasante (terreno del eje ± rango) se calculan las áreas reales de
     corte y relleno con la plantilla completa (peralte, sobreancho, cuneta, taludes).
  2. Programación dinámica (exacta sobre la malla): ruta de costo mínimo sujeta a |pendiente| ≤ máxima.
     La malla de cotas es ABSOLUTA y fina (misma para todas las estaciones): con una malla relativa al terreno,
     el incremento máximo alcanzable por estación queda por debajo del límite y la rasante no puede seguir
     una subida sostenida a la pendiente máxima. Se da un nivel de holgura a la pendiente (el ajuste posterior
     de PVI vuelve a imponer el límite exacto).
     El balance se maneja con un multiplicador λ: costo_nodo = c_c·Ac + c_f·Af + λ·(Ac − f·Af).
     λ solo tiene sentido en (−c_corte, c_relleno/f): fuera de ese rango excavar o rellenar saldría 'gratis'.
  3. Para varios λ se obtiene una rasante objetivo; los PVI y curvas verticales se ajustan a ella cumpliendo las normas
     (perfil.disenar_perfil con z_ref) y se evalúa el costo REAL con acarreo. Se elige la mejor, incluida la base,
     así que el resultado nunca es peor que la rasante base.
"""
from __future__ import annotations
import dataclasses
import types
import numpy as np

from . import normas, perfil as PF, secciones as SC


def dp_rasante(costo, z, dz_max):
    """Ruta de costo mínimo.  costo: (n, m) costo del nodo (estación k, nivel j); z: (n, m) cotas absolutas;
    dz_max: (n-1,) máximo cambio de cota entre estaciones consecutivas.  Devuelve (índices j_k, costo total)."""
    costo, z, dz_max = np.asarray(costo, float), np.asarray(z, float), np.asarray(dz_max, float)
    n, m = costo.shape
    acum = costo[0].copy()
    prev = np.zeros((n, m), dtype=int)
    for k in range(1, n):
        ok = np.abs(z[k][None, :] - z[k - 1][:, None]) <= dz_max[k - 1] + 1e-9      # (anterior, actual)
        cand = np.where(ok, acum[:, None], np.inf)
        prev[k] = np.argmin(cand, axis=0)
        acum = cand[prev[k], np.arange(m)] + costo[k]
    j = int(np.argmin(acum))
    if not np.isfinite(acum[j]):
        raise RuntimeError("Sin rasante factible dentro del rango de cotas con esa pendiente máxima.")
    total = float(acum[j])
    ruta = [j]
    for k in range(n - 1, 0, -1):
        j = int(prev[k][j]); ruta.append(j)
    return np.array(ruta[::-1]), total


def dp_malla(costo, W):
    """Ruta de costo mínimo sobre una malla ABSOLUTA de cotas común a todas las estaciones.
    costo: (n, M) con inf donde el nivel no está permitido; W: (n-1,) niveles máximos entre estaciones consecutivas.
    Usa mínimo deslizante (O(n·M)). Devuelve (índices de nivel por estación, costo total)."""
    from scipy.ndimage import minimum_filter1d
    costo = np.asarray(costo, float)
    n, M = costo.shape
    acum = np.empty((n, M))
    acum[0] = costo[0]
    for k in range(1, n):
        ventana = minimum_filter1d(acum[k - 1], size=2 * int(W[k - 1]) + 1, mode="constant", cval=np.inf)
        acum[k] = ventana + costo[k]
    j = int(np.argmin(acum[-1]))
    if not np.isfinite(acum[-1][j]):
        raise RuntimeError("Sin rasante factible dentro del rango de cotas con esa pendiente máxima.")
    total = float(acum[-1][j])
    ruta = [j]
    for k in range(n - 1, 0, -1):
        lo, hi = max(j - int(W[k - 1]), 0), min(j + int(W[k - 1]) + 1, M)
        j = lo + int(np.argmin(acum[k - 1][lo:hi]))
        ruta.append(j)
    return np.array(ruta[::-1]), total


def _resumen_candidato(sec):
    r = sec.resumen
    return dict(costo=r["costo"]["total"], corte_m3=r["corte_m3"], relleno_m3=r["relleno_compacto_m3"],
                excedente_m3=r["excedente_m3"], faltante_m3=r["faltante_m3"], acarreo_long_m3=r["acarreo_long_m3"],
                dma_m=r["dma_m"], costo_detalle=r["costo"])


def _factible(pf, s_pvi, z_pvi, al, vtr, L_prev, paso_L=10.0, A_min=0.003):
    """Curvas verticales si la configuración cumple pendiente máx./mín., tangentes mínimas y curvas dentro del
    perfil; None si no."""
    g = np.diff(z_pvi) / np.diff(s_pvi)
    if np.any(np.abs(g) > pf.gmax + 1e-9) or np.any(np.abs(g) < pf.gmin - 1e-9):
        return None
    curvas = PF._curvas(s_pvi, z_pvi, al, vtr, L_prev, paso_L)
    if PF._conflictos(s_pvi, curvas, A_min):
        return None
    return curvas


def pulir_pvis(terr, al, pf, st, prep, vtr, pasos_cota=(2.0, 1.0, 0.5, 0.25, 0.1), pasos_sta=(40.0, 20.0, 10.0),
               max_barridos=8, tol_rel=5e-4, verbose=False):
    """Búsqueda local (coordenadas) sobre las cotas y abscisas de los PVI minimizando el costo REAL de la sección
    completa (con acarreo, excedente y préstamo). Cada paso se acepta solo si baja el costo y se mantiene factible
    (pendientes, tangentes mínimas, curvas verticales). Devuelve (Perfil pulido, Secciones, nº de evaluaciones)."""
    s_pvi, z_pvi = pf.s_pvi.copy(), pf.z_pvi.copy()
    L_prev = np.array([0.0] + [c["L"] for c in pf.curvas] + [0.0])

    def costo(sp, zp):
        curvas = _factible(pf, sp, zp, al, vtr, L_prev)
        if curvas is None:
            return np.inf, None, None
        zd = PF.cota_diseno(pf.s, sp, zp, curvas)
        tmp = types.SimpleNamespace(s=pf.s, zd=zd, zt=pf.zt)
        sec = SC.calcular(terr, al, tmp, st, prep=prep)
        return sec.resumen["costo"]["total"], curvas, sec

    mejor, curvas, sec = costo(s_pvi, z_pvi)
    if not np.isfinite(mejor):
        return pf, None, 0
    n_eval = 1
    inicial = mejor
    for barrido in range(max_barridos):
        previo = mejor
        for i in range(len(s_pvi)):
            for paso in pasos_cota:
                for signo in (+1.0, -1.0):
                    zp = z_pvi.copy(); zp[i] += signo * paso
                    c, cu, se = costo(s_pvi, zp); n_eval += 1
                    if c < mejor - 1e-9:
                        mejor, z_pvi, curvas, sec = c, zp, cu, se
            if 0 < i < len(s_pvi) - 1:                            # los extremos conservan su abscisa
                for paso in pasos_sta:
                    for signo in (+1.0, -1.0):
                        sp = s_pvi.copy(); sp[i] += signo * paso
                        if not (sp[i - 1] < sp[i] < sp[i + 1]):
                            continue
                        c, cu, se = costo(sp, z_pvi); n_eval += 1
                        if c < mejor - 1e-9:
                            mejor, s_pvi, curvas, sec = c, sp, cu, se
        if verbose:
            print(f"   barrido {barrido + 1}: costo {mejor:,.0f} ({100 * (inicial - mejor) / inicial:.2f} % menos que el inicial)")
        if (previo - mejor) / previo < tol_rel:
            break
    zd = PF.cota_diseno(pf.s, s_pvi, z_pvi, curvas)
    nuevo = dataclasses.replace(pf, s_pvi=s_pvi, z_pvi=z_pvi, curvas=curvas, zd=zd, violaciones=[], avisos=[],
                                stats={}, modo=pf.modo + "+pulido",
                                log=pf.log + [f"pulido por costo real de sección: {n_eval} evaluaciones"])
    PF._verificar(nuevo, 0.003)
    sec = SC.calcular(terr, al, nuevo, st, prep=prep)
    return nuevo, sec, n_eval


def optimizar_rasante(terr, al, categoria, vtr, st: SC.SeccionTipo, base=None, rango=10.0, paso_offset=0.25,
                      paso_malla=0.05, lambdas=None, tol_dp_opciones=(1.0,), ds=10.0, pulir=True, verbose=False):
    """Devuelve dict: pf, sec (la mejor), base (pf, sec), candidatos (tabla), mejora_pct, aviso_malla, dp (λ, costo, N)."""
    prep = SC.preparar(terr, al, st, ds)
    s = prep.s
    zt_eje = terr.z(prep.x, prep.y)
    offs = np.arange(-rango, rango + paso_offset / 2, paso_offset)
    Ac, Af = SC.tabla_areas(prep, zt_eje, offs)
    gmax = normas.pend_media_max(categoria, vtr)
    dk = np.diff(s)
    w = np.r_[dk[0] / 2, (dk[:-1] + dk[1:]) / 2, dk[-1] / 2]           # pesos de áreas extremas
    # malla absoluta de cotas y áreas interpoladas en ella (las áreas varían suavemente con la cota)
    niveles = np.arange(np.floor(zt_eje.min() - rango), np.ceil(zt_eje.max() + rango) + paso_malla / 2, paso_malla)
    n, M = len(s), len(niveles)
    AcF = np.zeros((n, M)); AfF = np.zeros((n, M)); valido = np.zeros((n, M), dtype=bool)
    for k in range(n):
        h = niveles - zt_eje[k]
        m = np.abs(h) <= rango
        valido[k] = m
        AcF[k, m] = np.interp(h[m], offs, Ac[k]); AfF[k, m] = np.interp(h[m], offs, Af[k])
    W = np.ceil(gmax * dk / paso_malla - 1e-9).astype(int) + 1         # +1 nivel de holgura (ver docstring)
    if base is None:
        pf0 = PF.disenar_perfil(terr, al, categoria, vtr)
        base = (pf0, SC.calcular(terr, al, pf0, st, prep=prep))
    pf0, sec0 = base
    f = st.factor_terraplen
    if lambdas is None:                                                # rango donde los costos netos siguen siendo > 0
        lambdas = list(np.linspace(-0.9 * st.c_corte, 0.95 * st.c_relleno / f, 7))
    else:
        lambdas = list(lambdas)
    dp_tabla = []

    def objetivo(lam):
        nodo = np.where(valido, w[:, None] * (st.c_corte * AcF + st.c_relleno * AfF + lam * (AcF - f * AfF)), np.inf)
        ruta, total = dp_malla(nodo, W)
        kk = np.arange(n)
        N = float((w * (AcF[kk, ruta] - f * AfF[kk, ruta])).sum())      # balance neto (corte - f·relleno)
        dp_tabla.append(dict(lam=float(lam), costo_nodo=total, N=N,
                             Vc=float((w * AcF[kk, ruta]).sum()), Vf=float((w * AfF[kk, ruta]).sum())))
        return ruta, niveles[ruta]

    candidatos = [dict(lam=None, tol_dp=None, modo="base (mínimo corte/relleno en el eje)", **_resumen_candidato(sec0),
                       violaciones=len(pf0.violaciones), _pf=pf0, _sec=sec0)]
    tocan_borde = False

    def evaluar(lam, tol):
        nonlocal tocan_borde
        try:
            ruta, zobj = objetivo(lam)           # puede no haber rasante factible en la banda de cotas (ladera muy empinada)
        except (RuntimeError, ValueError) as e:
            if verbose:
                print(f"   λ={lam:+.2f}: sin rasante objetivo ({e})")
            return None
        try:
            pf = PF.disenar_perfil(terr, al, categoria, vtr, z_ref=(s, zobj), c_cut=1.0, c_fill=1.0, tol_dp=tol)
            sec = SC.calcular(terr, al, pf, st, prep=prep)
        except Exception as e:                                          # candidato inviable: se descarta
            if verbose:
                print(f"   λ={lam:+.2f} tol={tol}: descartado ({e})")
            return None
        if np.any(np.abs(niveles[ruta] - zt_eje) >= rango - paso_malla):
            tocan_borde = True
        c = dict(lam=float(lam), tol_dp=float(tol), modo="volúmenes", **_resumen_candidato(sec),
                 violaciones=len(pf.violaciones), _pf=pf, _sec=sec)
        candidatos.append(c)
        return c

    tol0 = tol_dp_opciones[0]
    for lam in lambdas:
        evaluar(lam, tol0)
    # refinamiento: λ intermedios alrededor del mejor y otras tolerancias de simplificación de PVI
    validos = [c for c in candidatos if c["lam"] is not None and c["violaciones"] == 0]
    if validos:
        mejor_lam = min(validos, key=lambda c: c["costo"])["lam"]
        paso_l = (lambdas[1] - lambdas[0]) if len(lambdas) > 1 else 0.5
        for lam in (mejor_lam - paso_l / 2, mejor_lam + paso_l / 2):
            evaluar(lam, tol0)
        for tol in tol_dp_opciones[1:]:
            evaluar(mejor_lam, tol)
    aptos = [c for c in candidatos if c["violaciones"] == 0 or c["lam"] is None]
    if pulir:                                   # búsqueda local sobre el costo real, desde la base y desde el mejor candidato DP
        puntos = [candidatos[0]] + [min([c for c in aptos if c["lam"] is not None], key=lambda c: c["costo"])] \
            if any(c["lam"] is not None for c in aptos) else [candidatos[0]]
        for c0 in puntos:
            pf_p, sec_p, n_ev = pulir_pvis(terr, al, c0["_pf"], st, prep, vtr, verbose=verbose)
            if sec_p is None:
                continue
            c = dict(lam=c0["lam"], tol_dp=c0["tol_dp"], modo=("base" if c0["lam"] is None else "volúmenes") + " + pulido",
                     **_resumen_candidato(sec_p), violaciones=len(pf_p.violaciones), _pf=pf_p, _sec=sec_p, evaluaciones=n_ev)
            candidatos.append(c)
        aptos = [c for c in candidatos if c["violaciones"] == 0 or c["lam"] is None]
    mejor = min(aptos, key=lambda c: c["costo"])
    base_costo = candidatos[0]["costo"]
    return dict(pf=mejor["_pf"], sec=mejor["_sec"], base=base, candidatos=candidatos, mejor=mejor,
                mejora_pct=100.0 * (base_costo - mejor["costo"]) / base_costo,
                aviso_malla=("Alguna rasante candidata toca el límite de la malla de cotas; ampliar `rango`." if tocan_borde else None),
                dp=dp_tabla, prep=prep)


def graficar_comparacion(opt, ruta):
    """Perfil y diagrama de masas: rasante base frente a la optimizada."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pf0, sec0 = opt["base"]
    pf1, sec1 = opt["pf"], opt["sec"]
    fig, axs = plt.subplots(3, 1, figsize=(13, 11), sharex=True, gridspec_kw={"height_ratios": [2.4, 1, 1.6]})
    ax = axs[0]
    ax.plot(pf0.s, pf0.zt, color="green", lw=1.3, label="Terreno en el eje")
    ax.plot(pf0.s, pf0.zd, color="0.45", lw=2, ls="--", label=f"Rasante base (costo {opt['candidatos'][0]['costo']:,.0f})")
    ax.plot(pf1.s, pf1.zd, color="red", lw=2.2, label=f"Rasante optimizada (costo {opt['mejor']['costo']:,.0f}; {opt['mejora_pct']:.1f} % menos)")
    ax.plot(pf1.s_pvi, pf1.z_pvi, "k^", ms=6, label="PVI optimizados")
    ax.set_ylabel("Cota (m)"); ax.grid(alpha=0.3); ax.legend(fontsize=8, loc="upper left")
    ax.set_title("Rasante base (mínimo corte/relleno en el eje) frente a rasante optimizada por costo real de la sección")
    ax = axs[1]
    ax.plot(pf0.s, pf0.zd - pf0.zt, color="0.45", ls="--", label="base")
    ax.plot(pf1.s, pf1.zd - pf1.zt, color="red", label="optimizada")
    ax.axhline(0, color="k", lw=0.6); ax.set_ylabel("Rasante - terreno (m)\n(+ relleno / - corte)"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax = axs[2]
    ax.plot(sec0.masas["s"], sec0.masas["ordenada"], color="0.45", ls="--", lw=1.8,
            label=f"base: excedente {sec0.resumen['excedente_m3']:.0f} m³, acarreo {sec0.resumen['acarreo_long_m3']:.0f} m³")
    ax.plot(sec1.masas["s"], sec1.masas["ordenada"], color="red", lw=2,
            label=f"optimizada: excedente {sec1.resumen['excedente_m3']:.0f} m³, faltante {sec1.resumen['faltante_m3']:.0f} m³, acarreo {sec1.resumen['acarreo_long_m3']:.0f} m³")
    ax.axhline(0, color="k", lw=0.6); ax.set_ylabel("Ordenada de masas (m³)"); ax.set_xlabel("Abscisa (m)")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    plt.tight_layout(); fig.savefig(ruta, dpi=95); plt.close(fig)
