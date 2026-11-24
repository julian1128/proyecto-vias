"""
Optimización local del trazado en planta: posición de los PI y radios de las curvas, por costo real.

El replanteo elige entre rutas ya existentes; pero el costo depende sobre todo de cómo se cierran las curvas
(la misma ruta varía hasta 2× entre ajustes de PI). Aquí, partiendo de un alineamiento válido, se mueven los PI
(± paso en x e y, con pasos decrecientes) y se cambia el radio de cada curva (×1.25 y ÷1.25), evaluando cada
configuración con el diseño completo: curvas con espirales, rasante, secciones, volúmenes y masas.
Criterio de aceptación (el de siempre): menos observaciones normativas y, a igualdad, menor costo total.
Se rechazan configuraciones con curvas solapadas, PI fuera del levantamiento o chaflanes con cota extrapolada.
"""
from __future__ import annotations
import numpy as np

from . import alineamiento as A, perfil as PF, secciones as SC


def evaluar_configuracion(terr, pis, radios, P_ref, categoria, vtr, st, e_max, ancho_carril, terreno_plano, tol_pi=0.0):
    """Diseño completo de una configuración de PI/radios. None si es inviable."""
    try:
        al = A.disenar_desde_pis(pis, P_ref, vtr, e_max=e_max, ancho_carril=ancho_carril,
                                 terreno_plano=terreno_plano, radios=radios, tol_pi=tol_pi)
        if not al.valido:
            return None
        if np.any(terr.tri.find_simplex(al.pis) < 0):
            return None
        pf = PF.disenar_perfil(terr, al, categoria, vtr)
        sec = SC.calcular(terr, al, pf, st)
    except (RuntimeError, ValueError):
        return None
    if any("fuera del área levantada" in a for a in sec.avisos):
        return None
    obs = len(al.violaciones) + len(pf.violaciones) + (1 if sec.resumen.get("eje_en_vacio_m", 0.0) > 0 else 0)
    return dict(key=(obs, sec.resumen["costo"]["total"]), al=al, pf=pf, sec=sec)


def _mejor(nuevo, actual, tol=1e-6):
    if nuevo is None:
        return False
    if nuevo["key"][0] != actual["key"][0]:
        return nuevo["key"][0] < actual["key"][0]
    return nuevo["key"][1] < actual["key"][1] - tol


def optimizar_trazado(terr, al0, P_ref, categoria, vtr, st, e_max=0.08, ancho_carril=3.65, terreno_plano=False,
                      pasos=None, max_evals=120, verbose=False, mover_pi=True, mover_radios=True):
    """Devuelve dict: al, pf, sec (el mejor), key0/key1 (obs, costo), mejora_pct, n_eval, historial.
    mover_pi / mover_radios permiten buscar solo en una de las dos familias de variables."""
    ext = max(terr.xmax - terr.xmin, terr.ymax - terr.ymin)
    pasos = list(pasos) if pasos is not None else [ext * 0.02, ext * 0.01, ext * 0.005, ext * 0.0025]
    args = (P_ref, categoria, vtr, st, e_max, ancho_carril, terreno_plano, al0.tol_pi)
    pis = al0.pis.copy()
    radios = [c["R"] for c in al0.curvas]
    actual = evaluar_configuracion(terr, pis, radios, *args)
    n_eval = 1
    if actual is None:                                    # no se pudo reevaluar: se devuelve el original sin tocar
        return dict(al=al0, pf=None, sec=None, key0=None, key1=None, mejora_pct=0.0, n_eval=1, historial=[])
    key0 = actual["key"]
    historial = [(n_eval, actual["key"])]

    def aceptar(nuevo, pis_n):
        nonlocal actual, pis, radios
        actual, pis = nuevo, pis_n
        radios = [c["R"] for c in nuevo["al"].curvas]
        historial.append((n_eval, nuevo["key"]))
        if verbose:
            print(f"   eval {n_eval:3d}: obs {nuevo['key'][0]}, costo {nuevo['key'][1]:,.0f}")

    for paso in pasos:
        mejoro = True
        while mejoro and n_eval < max_evals:
            mejoro = False
            for i in (range(1, len(pis) - 1) if mover_pi else []):  # mover cada PI interior en los 4 sentidos
                for d in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    if n_eval >= max_evals:
                        break
                    cand = pis.copy()
                    cand[i] = cand[i] + paso * np.array(d, float)
                    r = evaluar_configuracion(terr, cand, radios, *args); n_eval += 1
                    if _mejor(r, actual):
                        aceptar(r, cand); mejoro = True
            for i in (range(len(radios)) if mover_radios else []):  # cambiar el radio de cada curva
                for f in (1.25, 1 / 1.25):
                    if n_eval >= max_evals:
                        break
                    rad = list(radios); rad[i] = radios[i] * f
                    r = evaluar_configuracion(terr, pis, rad, *args); n_eval += 1
                    if _mejor(r, actual):
                        aceptar(r, pis); mejoro = True
    return dict(al=actual["al"], pf=actual["pf"], sec=actual["sec"], key0=key0, key1=actual["key"],
                mejora_pct=100.0 * (key0[1] - actual["key"][1]) / key0[1] if actual["key"][0] <= key0[0] else 0.0,
                n_eval=n_eval, historial=historial)
