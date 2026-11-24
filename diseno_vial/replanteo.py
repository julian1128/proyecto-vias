"""
Replanteo del trazado en planta según las condiciones de la vía.

La línea de ceros por compás solo depende de la pendiente máxima. Aquí se generan RUTAS CANDIDATAS:
  * compás (línea de ceros clásica de pendiente máxima),
  * rutas de costo mínimo con el costo transversal de la sección real (calzada, bermas, cuneta, taludes,
    espesor de pavimento) y distinto peso α del movimiento de tierras frente al costo por metro de vía
    (α = 0: la más corta factible; α alto: prefiere laderas suaves aunque alargue el trazado).
Cada candidata se diseña COMPLETA (alineamiento con espirales, rasante, secciones, volúmenes y masas) y se compara por
costo total real. Criterio de elección: primero menos observaciones normativas (horizontales y verticales), luego
menor costo.  Al cambiar el ancho de calzada o el tipo de corte, la ruta elegida puede cambiar.
"""
from __future__ import annotations
import numpy as np
from scipy.spatial import cKDTree

from . import normas, linea_ceros as LC, alineamiento as A, perfil as PF, secciones as SC, optimizar as OPT


def _evaluar(terr, nombre, P, categoria, vtr, st, e_max, ancho_carril, plano, tol_pi):
    """Diseña completas todas las variantes de ajuste de PI de la ruta y devuelve la mejor (menos observaciones,
    luego menor costo). Guarda también el costo de la variante que elegiría el criterio de fidelidad geométrica."""
    fila = dict(nombre=nombre, P=P, valida=False, error=None, variantes=0)
    try:
        vars_ = A.variantes(P, vtr, e_max=e_max, ancho_carril=ancho_carril, terreno_plano=plano, tol_pi=tol_pi)
    except (RuntimeError, ValueError) as e:            # inviabilidad esperada; cualquier otro error es un defecto y debe verse
        fila["error"] = str(e)[:120]
        return fila
    if not vars_:
        fila["error"] = "ninguna variante de alineamiento es válida (curvas solapadas)"
        return fila
    resultados = []
    for al in vars_:
        try:
            pf = PF.disenar_perfil(terr, al, categoria, vtr)
            sec = SC.calcular(terr, al, pf, st)
        except (RuntimeError, ValueError):
            continue
        obs = len(al.violaciones) + len(pf.violaciones) + (1 if sec.resumen.get("eje_en_vacio_m", 0.0) > 0 else 0)
        resultados.append((obs, sec.resumen["costo"]["total"], al, pf, sec))
    if not resultados:
        fila["error"] = "ninguna variante admite perfil y secciones"
        return fila
    obs, costo, al, pf, sec = min(resultados, key=lambda t: (t[0], t[1]))
    r = sec.resumen
    fila.update(valida=True, al=al, pf=pf, sec=sec, costo=costo, longitud=al.longitud, obs=obs,
                obs_h=len(al.violaciones), obs_v=len(pf.violaciones), curvas=len(al.curvas), corte=r["corte_m3"],
                relleno=r["relleno_compacto_m3"], excedente=r["excedente_m3"], faltante=r["faltante_m3"],
                desv=al.desv_media, variantes=len(resultados), tol_pi=al.tol_pi,
                costo_peor_variante=max(t[1] for t in resultados))
    return fila


def _separacion(Pa, Pb):
    """Distancia media y máxima entre dos trazados (m)."""
    d = cKDTree(Pb[:, :2]).query(A._densificar(Pa[:, :2], 10.0))[0]
    return float(d.mean()), float(d.max())


def replantear(terr, origen, destino, categoria, vtr, st, equidistancia=5.0, e_max=0.08, ancho_carril=3.65,
               terreno_plano=False, tol_pi=None, alfas=(0.0, 0.5, 1.0, 2.0, 4.0), verbose=False):
    """Devuelve dict: candidatas (lista de filas), mejor (fila), referencia (la del compás, si es válida)."""
    pend = normas.pend_media_max(categoria, vtr)
    ext = max(terr.xmax - terr.xmin, terr.ymax - terr.ymin)
    tol_pi = tol_pi or max(ext / 80.0, 5.0)
    rutas = []
    try:
        P0, _ = LC.linea_de_ceros(terr, origen, destino, pend, dh=equidistancia)
        rutas.append(("compás (pendiente máx.)", P0))
    except RuntimeError as e:
        if verbose:
            print("   compás sin solución:", e)
    tabla = LC.tabla_costo_transversal(st)
    for a in alfas:
        try:
            P, _ = LC.ruta_costo_minimo(terr, origen, destino, pend, st, alfa=a, tabla=tabla)
            rutas.append((f"costo mínimo α={a:g}" + (" (más corta)" if a == 0 else ""), P))
        except RuntimeError as e:
            if verbose:
                print(f"   costo mínimo α={a}: sin solución ({e})")
    if not rutas:
        raise RuntimeError("Ninguna ruta candidata es factible con esta pendiente máxima y esta sección.")
    cands = [_evaluar(terr, n, P, categoria, vtr, st, e_max, ancho_carril, terreno_plano, tol_pi) for n, P in rutas]
    validas = [c for c in cands if c["valida"]]
    if not validas:
        raise RuntimeError("Ninguna ruta candidata produjo un alineamiento válido: " +
                           "; ".join(f"{c['nombre']}: {c['error']}" for c in cands))
    mejor = min(validas, key=lambda c: (c["obs"], c["costo"]))
    ref = next((c for c in cands if c["nombre"].startswith("compás") and c["valida"]), None)
    for c in validas:
        c["sep_mejor_media"], c["sep_mejor_max"] = _separacion(c["P"], mejor["P"])
    return dict(candidatas=cands, mejor=mejor, referencia=ref)
