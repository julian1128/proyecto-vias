"""
Alineamiento horizontal: de la línea de ceros a tangentes + PI + curvas espiralizadas.

Pasos (diapositivas "Trazado del alineamiento inicial"):
  1. Tangentes tan largas como sea posible: simplificación de la línea de ceros (Douglas-Peucker)
     y ajuste de una recta por tramo. Los PI son la intersección de rectas consecutivas, por eso
     NO coinciden con la línea de ceros: la curva es la que se ajusta a ella.
  2. Velocidad específica de cada curva (VECH) según el segmento recto anterior.
  3. Radio: de Rmín (AASHTO, e_max) hacia arriba, el que mejor se ajusta a la línea de ceros
     respetando espacio de tangentes, entretangencia mínima y segmento circular mínimo.
  4. Espiral: longitud mínima por rampa de peralte / comodidad / 2 s; sin espiral si R > 1000 m.
"""
from __future__ import annotations
import math
from dataclasses import dataclass, field
import numpy as np
from scipy.spatial import cKDTree

from . import normas

LC_MIN = 20.0          # longitud mínima del segmento circular [DIAP]
PI_DUPLICADO = 0.5     # PI a menos de esta distancia (m) se consideran el mismo punto


# ------------------------------------------------------------------ utilidades geométricas
def azimut(v) -> float:
    """Azimut en grados desde el norte, horario. v = (dx_este, dy_norte)."""
    return math.degrees(math.atan2(v[0], v[1])) % 360.0


def _unit(v):
    n = max(float(np.hypot(*v)), 1e-12)                # un vector nulo da (0,0) en lugar de NaN
    return np.asarray(v, float) / n


def _normal_izq(u):
    return np.array([-u[1], u[0]])


def douglas_peucker(P: np.ndarray, tol: float) -> list[int]:
    keep = {0, len(P) - 1}
    stack = [(0, len(P) - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1:
            continue
        ab = P[b] - P[a]
        L = np.hypot(*ab)
        seg = P[a + 1:b] - P[a]
        dist = np.abs(ab[0] * seg[:, 1] - ab[1] * seg[:, 0]) / L if L > 0 else np.hypot(*seg.T)
        k = int(np.argmax(dist))
        if dist[k] > tol:
            m = a + 1 + k
            keep.add(m); stack += [(a, m), (m, b)]
    return sorted(keep)


def _ajuste_recta(P):
    c = P.mean(axis=0)
    _, _, vt = np.linalg.svd(P - c)
    u = vt[0]
    if np.dot(u, P[-1] - P[0]) < 0:
        u = -u
    return c, u


def _interseccion(c1, u1, c2, u2):
    A = np.array([u1, -u2]).T
    if abs(np.linalg.det(A)) < 1e-9:
        return None
    t = np.linalg.solve(A, c2 - c1)
    return c1 + t[0] * u1


def _pis_desde_indices(P2, idx, tol):
    rectas = [_ajuste_recta(P2[a:b + 1]) for a, b in zip(idx[:-1], idx[1:])]
    pis = []
    for k in range(len(rectas) - 1):
        pi = _interseccion(*rectas[k], *rectas[k + 1])
        v = P2[idx[k + 1]]
        if pi is None or np.hypot(*(pi - v)) > max(8 * tol, 80.0):
            pi = v
        pis.append(pi)
    return np.vstack([P2[0]] + pis + [P2[-1]]) if pis else np.vstack([P2[0], P2[-1]])


def _deflexiones(pts):
    u = np.diff(pts, axis=0); u /= np.maximum(np.linalg.norm(u, axis=1), 1e-12)[:, None]
    return np.degrees(np.arccos(np.clip((u[:-1] * u[1:]).sum(1), -1, 1)))


def ajustar_pis(P2: np.ndarray, tol: float = 6.0, defl_min_deg: float = 2.0, defl_max_deg: float = 70.0):
    """PI a partir de la línea de ceros (Nx2). Devuelve [BP, PI1..PIn, EP].
    - quiebres casi colineales (< defl_min) se eliminan
    - deflexiones grandes (> defl_max) se dividen en varios PI (curvas más suaves y viables)"""
    idx = douglas_peucker(P2, tol)
    for _ in range(200):
        pts = _pis_desde_indices(P2, idx, tol)
        if len(pts) <= 2:
            return pts
        # PI coincidentes (la división de deflexiones grandes puede reutilizar un mismo vértice): se elimina el vértice
        # sobrante del bucle (no el punto de `pts`, para no desalinear los índices con `idx`)
        dup = np.where(np.hypot(*np.diff(pts, axis=0).T) < PI_DUPLICADO)[0]
        if len(dup) and len(idx) > 2:
            quitar = {int(j) + 1 if int(j) + 1 < len(pts) - 1 else int(j) for j in dup}
            for k in sorted((k for k in quitar if 0 < k < len(idx) - 1), reverse=True):
                del idx[k]
            continue
        ang = _deflexiones(pts)
        peq = np.where(ang < defl_min_deg)[0]
        if len(peq):
            # se eliminan de una vez todos los quiebres casi colineales (la ruta raster puede traer decenas)
            for k in sorted({int(j) + 1 for j in peq}, reverse=True):
                del idx[k]
            continue
        grandes = np.where(ang > defl_max_deg)[0]
        if len(grandes):
            k = grandes[0] + 1                       # vértice idx[k] concentra el cambio de rumbo
            a, b = idx[k - 1], idx[k + 1]
            nuevos = [a + j for j in douglas_peucker(P2[a:b + 1], tol * 0.3)]
            if len(nuevos) <= 3 and set(nuevos) <= set(idx):
                nuevos = sorted(set(nuevos) | {(a + idx[k]) // 2, (idx[k] + b) // 2})
            nidx = sorted(set(idx) | set(nuevos))
            if nidx == idx:
                return pts
            idx = nidx
            continue
        return pts
    return pts


# ------------------------------------------------------------------ espiral (clotoide)
def espiral_xy(s, Ls, R):
    th = s * s / (2.0 * R * Ls)
    x = s * (1 - th**2 / 10 + th**4 / 216 - th**6 / 9360)
    y = s * (th / 3 - th**3 / 42 + th**5 / 1320 - th**7 / 75600)
    return x, y


def espiral_params(Ls, R):
    th = Ls / (2 * R)
    Xs, Ys = espiral_xy(Ls, Ls, R)
    p = Ys - R * (1 - math.cos(th))
    k = Xs - R * math.sin(th)
    return dict(theta=th, Xs=Xs, Ys=Ys, p=p, k=k,
                tan_larga=Xs - Ys / math.tan(th) if th > 1e-9 else Ls,
                tan_corta=Ys / math.sin(th) if th > 1e-9 else 0.0)


# ------------------------------------------------------------------ diseño de una curva
def curva_en_pi(pi, u1, u2, R, Le):
    """Geometría de una curva espiral-circular-espiral en un PI. u1/u2 tangentes unitarias de entrada/salida."""
    cross = u1[0] * u2[1] - u1[1] * u2[0]
    sgn = 1 if cross > 0 else -1                      # +1 izquierda (antihorario)
    delta = math.acos(np.clip(np.dot(u1, u2), -1, 1))
    if Le > 0:
        sp = espiral_params(Le, R)
        th, p, k = sp["theta"], sp["p"], sp["k"]
    else:
        sp, th, p, k = None, 0.0, 0.0, 0.0
    Ts = (R + p) * math.tan(delta / 2) + k
    n1, n2 = _normal_izq(u1), _normal_izq(u2)
    TE = pi - Ts * u1
    ST = pi + Ts * u2
    C = TE + k * u1 + (R + p) * sgn * n1
    if Le > 0:
        SC = TE + sp["Xs"] * u1 + sp["Ys"] * sgn * n1
        CS = ST - sp["Xs"] * u2 + sp["Ys"] * sgn * n2
    else:
        SC, CS = TE, ST
    Lc = R * (delta - 2 * th)
    return dict(sgn=sgn, delta=delta, theta=th, Ts=Ts, TE=TE, SC=SC, CS=CS, ST=ST, C=C,
                Lc=Lc, Le=Le, R=R, sp=sp, PI=pi, u1=u1, u2=u2)


def _muestras_curva(c, ds=2.0):
    """Puntos densos TE -> ST de una curva (para medir ajuste y dibujar)."""
    pts = []
    sgn, R, Le = c["sgn"], c["R"], c["Le"]
    n1, n2 = _normal_izq(c["u1"]), _normal_izq(c["u2"])
    if Le > 0:
        for s in np.linspace(0, Le, max(int(Le / ds), 3)):
            x, y = espiral_xy(s, Le, R)
            pts.append(c["TE"] + x * c["u1"] + y * sgn * n1)
    a0 = math.atan2(*(c["SC"] - c["C"])[::-1]); da = sgn * (c["delta"] - 2 * c["theta"])
    for t in np.linspace(0, 1, max(int(abs(c["Lc"]) / ds), 3)):
        a = a0 + da * t
        pts.append(c["C"] + R * np.array([math.cos(a), math.sin(a)]))
    if Le > 0:
        for s in np.linspace(Le, 0, max(int(Le / ds), 3)):
            x, y = espiral_xy(s, Le, R)
            pts.append(c["ST"] - x * c["u2"] + y * sgn * n2)
    return np.array(pts)


# ------------------------------------------------------------------ diseño del alineamiento
@dataclass
class Alineamiento:
    pis: np.ndarray
    curvas: list
    elementos: list            # dicts: linea / espiral / curva, con estacionamiento
    longitud: float
    violaciones: list = field(default_factory=list)
    vtr: float = 0.0
    solapes: list = field(default_factory=list)   # curvas que se pisan (geometría inválida)
    valido: bool = True
    tol_pi: float = 0.0                            # tolerancia de simplificación finalmente usada
    desv_media: float = 0.0                        # distancia media eje <-> línea de ceros (m)
    desv_max: float = 0.0
    avisos: list = field(default_factory=list)

    def muestrear(self, ds=2.0):
        """Polilínea densa [x, y, abscisa] del eje."""
        out, s0 = [], 0.0
        if not self.elementos:
            return np.empty((0, 3))
        for e in self.elementos:
            if e["tipo"] == "linea":
                a, b = np.array(e["ini"]), np.array(e["fin"])
                if not np.isfinite(np.hypot(*(b - a))) or np.hypot(*(b - a)) > 1e6:
                    raise ValueError("Geometría absurda en el alineamiento (elemento recto de más de 1.000 km): "
                                     "curvas mal condicionadas")
                n = max(int(np.hypot(*(b - a)) / ds), 1)
                pts = [a + (b - a) * t for t in np.linspace(0, 1, n + 1)]
            else:
                pts = e["_pts"]
            for p in pts:
                out.append(p)
        P = np.array(out)
        dd = np.r_[0, np.hypot(*np.diff(P, axis=0).T)]
        keep = dd > 1e-9; keep[0] = True
        P = P[keep]
        s = np.r_[0, np.cumsum(np.hypot(*np.diff(P, axis=0).T))]
        return np.column_stack([P, s])


def _quitar_pi_colineales(pis, tol_deg=0.05):
    keep = [0]                                         # PI colineales (deflexión ~0) no generan curva
    for i in range(1, len(pis) - 1):
        if np.hypot(*(pis[i] - pis[keep[-1]])) < PI_DUPLICADO:     # PI repetido
            continue
        u1, u2 = _unit(pis[i] - pis[keep[-1]]), _unit(pis[i + 1] - pis[i])
        if np.degrees(np.arccos(np.clip(np.dot(u1, u2), -1, 1))) > tol_deg:
            keep.append(i)
    if len(keep) > 1 and np.hypot(*(pis[-1] - pis[keep[-1]])) < PI_DUPLICADO:
        keep.pop()                                     # el último PI interior coincide con el punto final
    keep.append(len(pis) - 1)
    return pis[keep]


def _disenar_una(P_ceros: np.ndarray, vtr: float, e_max: float = 0.08, ancho_carril: float = 3.65,
                 bombeo: float = 0.02, terreno_plano: bool = False, tol_pi: float = 6.0,
                 n_radios: int = 14, factor_radio_max: float = 8.0) -> Alineamiento:
    pis = _quitar_pi_colineales(ajustar_pis(P_ceros[:, :2], tol=tol_pi))
    return disenar_desde_pis(pis, P_ceros, vtr, e_max, ancho_carril, bombeo, terreno_plano, n_radios,
                             factor_radio_max, tol_pi=tol_pi)


def _alineamiento_invalido(pis, vtr, tol_pi, msg):
    return Alineamiento(pis=np.asarray(pis, float), curvas=[], elementos=[], longitud=float("nan"), violaciones=[msg],
                        vtr=vtr, solapes=[msg], valido=False, tol_pi=tol_pi,
                        desv_media=float("inf"), desv_max=float("inf"))


def disenar_desde_pis(pis, P_ref: np.ndarray, vtr: float, e_max: float = 0.08, ancho_carril: float = 3.65,
                      bombeo: float = 0.02, terreno_plano: bool = False, n_radios: int = 14,
                      factor_radio_max: float = 8.0, radios=None, tol_pi: float = 0.0) -> Alineamiento:
    """Diseña las curvas (espiral-circular-espiral) para PI dados. `radios`: lista (un valor por PI interior, o None)
    para fijar el radio de cada curva (se limita a Rmín de su velocidad específica); sin fijar, se elige el radio
    que mejor se ajusta a la ruta de referencia `P_ref`."""
    pis = np.asarray(pis, float)
    P2 = P_ref[:, :2]
    if len(pis) >= 2 and np.any(np.hypot(*np.diff(pis, axis=0).T) < PI_DUPLICADO):
        return _alineamiento_invalido(pis, vtr, tol_pi, "PI coincidentes (a menos de 0,5 m): el trazado tiene un tramo "
                                                       "recto de longitud nula")
    dens = _densificar(P2, 2.0)
    arbol = cKDTree(dens)
    curvas, viol = [], []
    prev_end = pis[0].copy()           # fin del elemento anterior (BP o ST)
    vech_ant, ndelta = None, len(pis) - 2
    for i in range(1, len(pis) - 1):
        pi = pis[i]
        v1, v2 = pi - pis[i - 1], pis[i + 1] - pi
        u1, u2 = _unit(v1), _unit(v2)
        delta = math.acos(np.clip(np.dot(u1, u2), -1, 1))
        d_libre_ant = float(np.hypot(*(pi - prev_end)))
        L_sig = float(np.hypot(*v2))
        cap_sig = 0.5 * L_sig if i < len(pis) - 2 else 0.95 * L_sig
        V = normas.vech(vtr, vech_ant, d_libre_ant, math.degrees(delta))
        V = float(min(max(round(V / 10) * 10, 20), 130))
        Rmin = normas.r_min(V, e_max)
        mejor = None
        fijo = radios[i - 1] if (radios is not None and i - 1 < len(radios)) else None
        if fijo is not None:                                   # radio impuesto (nunca menor que el mínimo de la velocidad)
            Rs = [max(float(fijo), Rmin)]
        elif math.degrees(delta) < 6.0:                        # deflexión pequeña: circular simple, R grande
            Lmin = 30 * (10 - math.degrees(delta))             # [VERIFICAR] tabla 3.10 del manual
            Rs = [max(Lmin / delta, 1000.0)]
        else:
            Rs = list(Rmin * np.geomspace(1.0, factor_radio_max, n_radios))
        for R in Rs:
            e = normas.peralte(R, V, e_max, bombeo)
            Le = 0.0
            if normas.requiere_espiral(R):
                Le = math.ceil(normas.longitud_espiral_min(V, R, e, ancho_carril, bombeo)["Le_min"] / 5) * 5
            c = curva_en_pi(pi, u1, u2, R, Le)
            ok_circ = c["Lc"] >= (LC_MIN if Le > 0 else 0.0)
            ok_prev = d_libre_ant >= c["Ts"] + normas.entretangencia_min(V, Le > 0, terreno_plano)
            ok_sig = c["Ts"] <= cap_sig
            if not (ok_circ and ok_prev and ok_sig):
                continue
            pts = _muestras_curva(c)
            dev = float(arbol.query(pts)[0].mean())
            if mejor is None or dev < mejor[0] - 0.05:
                mejor = (dev, c, e, V)
        if mejor is None:                                       # sin solución normativa: Rmín y se reporta
            R = Rs[0]
            e = normas.peralte(R, V, e_max, bombeo)
            Le = (math.ceil(normas.longitud_espiral_min(V, R, e, ancho_carril, bombeo)["Le_min"] / 5) * 5
                  if normas.requiere_espiral(R) else 0.0)
            Le = min(Le, R * delta * 0.98) if Le > 0 else 0.0
            c = curva_en_pi(pi, u1, u2, R, Le)
            if (not np.isfinite(c["Ts"])) or c["Ts"] > 20.0 * max(d_libre_ant, L_sig, 1000.0):
                # curva de retorno (deflexión cercana a 180°): una curva simple espiral-circular-espiral no la resuelve
                # y su tangente Ts = R·tan(Δ/2) se dispara; se declara inviable sin construir esa geometría
                msg = (f"PI-{i}: deflexión de {math.degrees(delta):.0f}° (curva de retorno): no se puede resolver con una "
                       f"curva simple (tangente {c['Ts']:.3g} m); requiere una curva de volteo")
                al_inv = _alineamiento_invalido(pis, vtr, tol_pi, msg)
                al_inv.violaciones = viol + [msg]
                return al_inv
            motivo = []
            if Le > 0 and c["Lc"] < LC_MIN: motivo.append(f"segmento circular {c['Lc']:.0f} m < {LC_MIN:.0f} m")
            if d_libre_ant < c["Ts"] + normas.entretangencia_min(V, Le > 0, terreno_plano): motivo.append("entretangencia anterior insuficiente")
            if c["Ts"] > cap_sig: motivo.append("tangente siguiente insuficiente")
            viol.append(f"PI-{i}: {'; '.join(motivo) or 'sin ajuste'} con R={R:.0f} m (V={V:.0f} km/h)")
            mejor = (float(arbol.query(_muestras_curva(c))[0].mean()), c, e, V)
        dev, c, e, V = mejor
        c.update(e=e, V=V, dev=dev, idx=i)
        curvas.append(c)
        prev_end = c["ST"]; vech_ant = V
    # elementos con estacionamiento
    elems, s = [], 0.0
    cursor = pis[0].copy()
    for c in curvas:
        L = float(np.hypot(*(c["TE"] - cursor)))
        elems.append(dict(tipo="linea", ini=tuple(cursor), fin=tuple(c["TE"]), s0=s, longitud=L)); s += L
        c["s_ini"] = s                                   # abscisa del TE (inicio de la curva espiralizada)
        if c["Le"] > 0:
            sp = c["sp"]
            elems.append(dict(tipo="espiral", entrada=True, ini=tuple(c["TE"]), fin=tuple(c["SC"]), pi=tuple(c["TE"] + sp["tan_larga"] * c["u1"]),
                              longitud=c["Le"], R=c["R"], sgn=c["sgn"], sp=sp, s0=s,
                              _pts=_muestras_espiral(c, True)))
            s += c["Le"]
        elems.append(dict(tipo="curva", ini=tuple(c["SC"]), fin=tuple(c["CS"]), centro=tuple(c["C"]), R=c["R"], sgn=c["sgn"],
                          delta=c["delta"] - 2 * c["theta"], longitud=c["Lc"], s0=s, curva=c, _pts=_muestras_arco(c)))
        s += c["Lc"]
        if c["Le"] > 0:
            elems.append(dict(tipo="espiral", entrada=False, ini=tuple(c["CS"]), fin=tuple(c["ST"]), pi=tuple(c["ST"] - c["sp"]["tan_larga"] * c["u2"]),
                              longitud=c["Le"], R=c["R"], sgn=c["sgn"], sp=c["sp"], s0=s,
                              _pts=_muestras_espiral(c, False)))
            s += c["Le"]
        c["s_fin"] = s                                   # abscisa del ST (fin de la curva espiralizada)
        cursor = c["ST"]
    L = float(np.hypot(*(pis[-1] - cursor)))
    elems.append(dict(tipo="linea", ini=tuple(cursor), fin=tuple(pis[-1]), s0=s, longitud=L)); s += L
    # Un elemento recto cuyo avance sobre su tangente es negativo indica curvas que se pisan.
    solapes = []
    rectas = [e for e in elems if e["tipo"] == "linea"]
    defl = _deflexiones(pis) if len(pis) > 2 else np.array([])
    for k, e in enumerate(rectas):
        u = _unit(pis[k + 1] - pis[k])
        avance = float(np.dot(np.array(e["fin"]) - np.array(e["ini"]), u))
        if avance < -1e-6:
            adyacente = float(defl[max(k - 1, 0):k + 1].max()) if len(defl) else 0.0
            pista = (f"; deflexión de {adyacente:.0f}°: curva de retorno, requiere una curva de volteo"
                     if adyacente > 150.0 else "")
            solapes.append(f"tangente {k} (PI-{k} a PI-{k + 1}): las curvas se pisan {abs(avance):.1f} m{pista}")
    al = Alineamiento(pis=pis, curvas=curvas, elementos=elems, longitud=s, violaciones=viol, vtr=vtr,
                      solapes=solapes, valido=not solapes, tol_pi=tol_pi)
    d = cKDTree(al.muestrear(2.0)[:, :2]).query(_densificar(P2, 10.0))[0]      # fidelidad a la línea de ceros
    al.desv_media, al.desv_max = float(d.mean()), float(d.max())
    return al


def _intentos(P_ceros, vtr, e_max, ancho_carril, bombeo, terreno_plano, tol_pi, n_radios, factor_radio_max, factores_tol):
    return [_disenar_una(P_ceros, vtr, e_max, ancho_carril, bombeo, terreno_plano, tol_pi * f, n_radios, factor_radio_max)
            for f in factores_tol]


FACTORES_TOL = (1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0)


def variantes(P_ceros: np.ndarray, vtr: float, e_max: float = 0.08, ancho_carril: float = 3.65,
              bombeo: float = 0.02, terreno_plano: bool = False, tol_pi: float = 6.0,
              n_radios: int = 14, factor_radio_max: float = 8.0,
              factores_tol=(0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0)):
    """Alineamientos VÁLIDOS distintos que se obtienen al variar la tolerancia de simplificación de PI (sin duplicados).
    Permite elegir por costo real en vez de por fidelidad geométrica."""
    vistos, out = set(), []
    for al in _intentos(P_ceros, vtr, e_max, ancho_carril, bombeo, terreno_plano, tol_pi, n_radios, factor_radio_max, factores_tol):
        if not al.valido:
            continue
        clave = tuple(np.round(al.pis, 0).ravel()) + tuple(round(c["R"]) for c in al.curvas)
        if clave in vistos:
            continue
        vistos.add(clave); out.append(al)
    return out


def disenar(P_ceros: np.ndarray, vtr: float, e_max: float = 0.08, ancho_carril: float = 3.65,
            bombeo: float = 0.02, terreno_plano: bool = False, tol_pi: float = 6.0,
            n_radios: int = 14, factor_radio_max: float = 8.0,
            factores_tol=FACTORES_TOL) -> Alineamiento:
    """Diseña el alineamiento. Si con la tolerancia pedida las curvas no caben (solapes u observaciones
    normativas), reintenta con tolerancias mayores (menos PI, tangentes más largas) y se queda con la
    primera solución sin observaciones; si ninguna lo logra, con la de geometría válida y menos
    observaciones; si todas se pisan, devuelve la menos mala marcada como inválida."""
    intentos = _intentos(P_ceros, vtr, e_max, ancho_carril, bombeo, terreno_plano, tol_pi, n_radios, factor_radio_max, factores_tol)
    extension = float(max(np.ptp(P_ceros[:, 0]), np.ptp(P_ceros[:, 1])))
    umbral = max(50.0, 0.03 * extension)                      # desviación máxima aceptable respecto a la línea de ceros
    validos = [a for a in intentos if a.valido]
    if not validos:
        mejor = min(intentos, key=lambda a: (len(a.solapes), len(a.violaciones)))
        mejor.avisos.append("Todas las tolerancias producen curvas solapadas: el trazado no es viable con estas especificaciones.")
        return mejor
    # Se prefiere un eje fiel al terreno: solo cuentan los diseños que no se alejan más de `umbral` de la línea
    # de ceros. Una recta sin curvas "cumple" normas por vacío pero ignora el terreno: no se elige por eso.
    cerca = [a for a in validos if a.desv_max <= umbral]
    pool = cerca or validos
    mejor = min(pool, key=lambda a: (len(a.violaciones), a.desv_media, a.tol_pi))
    if not cerca:
        mejor.avisos.append(f"Ningún diseño válido se mantiene a menos de {umbral:.0f} m de la línea de ceros "
                            f"(mejor: {mejor.desv_max:.0f} m): movimiento de tierras alto, revisar especificaciones.")
    if mejor.violaciones:
        mejor.avisos.append(f"El mejor diseño tiene {len(mejor.violaciones)} observación(es) normativa(s).")
    return mejor
    # Se prefiere un eje fiel al terreno: solo cuentan los diseños que no se alejan más de `umbral` de la línea
    # de ceros. Una recta sin curvas "cumple" normas por vacío pero ignora el terreno: no se elige por eso.
    cerca = [a for a in validos if a.desv_max <= umbral]
    pool = cerca or validos
    mejor = min(pool, key=lambda a: (len(a.violaciones), a.desv_media, a.tol_pi))
    if not cerca:
        mejor.avisos.append(f"Ningún diseño válido se mantiene a menos de {umbral:.0f} m de la línea de ceros "
                            f"(mejor: {mejor.desv_max:.0f} m): movimiento de tierras alto, revisar especificaciones.")
    if mejor.violaciones:
        mejor.avisos.append(f"El mejor diseño tiene {len(mejor.violaciones)} observación(es) normativa(s).")
    return mejor


def _densificar(P, ds):
    out = []
    for a, b in zip(P[:-1], P[1:]):
        n = max(int(np.hypot(*(b - a)) / ds), 1)
        out += [a + (b - a) * t for t in np.linspace(0, 1, n, endpoint=False)]
    out.append(P[-1])
    return np.array(out)


def _muestras_espiral(c, entrada, ds=2.0):
    Le, R, sgn = c["Le"], c["R"], c["sgn"]
    if entrada:
        n1 = _normal_izq(c["u1"]); ss = np.linspace(0, Le, max(int(Le / ds), 3))
        return [c["TE"] + espiral_xy(s, Le, R)[0] * c["u1"] + espiral_xy(s, Le, R)[1] * sgn * n1 for s in ss]
    n2 = _normal_izq(c["u2"]); ss = np.linspace(Le, 0, max(int(Le / ds), 3))
    return [c["ST"] - espiral_xy(s, Le, R)[0] * c["u2"] + espiral_xy(s, Le, R)[1] * sgn * n2 for s in ss]


def _muestras_arco(c, ds=2.0):
    a0 = math.atan2(*(c["SC"] - c["C"])[::-1]); da = c["sgn"] * (c["delta"] - 2 * c["theta"])
    return [c["C"] + c["R"] * np.array([math.cos(a0 + da * t), math.sin(a0 + da * t)])
            for t in np.linspace(0, 1, max(int(abs(c["Lc"]) / ds), 3))]


def tabla_curvas(al: Alineamiento):
    """Filas de resumen por curva (para el reporte)."""
    rows = []
    for c in al.curvas:
        rows.append(dict(PI=c["idx"], deflexion_deg=math.degrees(c["delta"]), V_kmh=c["V"], R_m=c["R"],
                         peralte_pct=100 * c["e"], Le_m=c["Le"], Lc_m=c["Lc"], Ts_m=c["Ts"],
                         desviacion_media_m=c["dev"]))
    return rows
