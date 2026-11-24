"""
Línea de ceros (línea de pendiente) por el método del compás.

Procedimiento de las diapositivas: desde el punto inicial se abre el compás a
    d = equidistancia / pendiente
y se marca la siguiente curva de nivel (cota ± equidistancia); en cada paso hay dos
intersecciones posibles. Criterios de selección: menor longitud, menor variación de
pendiente y menores ángulos de deflexión.

Aquí la ramificación se resuelve con una búsqueda A*: costo = distancia + penalización por
deflexión + penalización por contrapendiente. Cuando el terreno es más plano que la pendiente
máxima, la línea sigue recta hacia el destino (pendiente < máxima, sobre el terreno).
"""
from __future__ import annotations
import heapq, math
import numpy as np

from .normas import longitud_virtual


def _ang_diff(a, b):
    d = (a - b + math.pi) % (2 * math.pi) - math.pi
    return abs(d)


def linea_de_ceros(terr, origen, destino, pend: float, dh: float = 1.0,
                   w_defl: float = 1.5, w_contra: float = 1.0, paso_ang_deg: float = 3.0,
                   max_nodos: int = 150_000):
    """Devuelve (P Nx3 [x,y,z sobre el terreno], info)."""
    d = dh / pend
    z0 = float(terr.z([origen[0]], [origen[1]])[0])
    zD = float(terr.z([destino[0]], [destino[1]])[0])
    ang = np.radians(np.arange(0, 360 + paso_ang_deg, paso_ang_deg))
    ca, sa = np.cos(ang), np.sin(ang)

    X, Y, Z, G, PAR, AZ = [origen[0]], [origen[1]], [z0], [0.0], [-1], [None]
    h = lambda x, y: math.hypot(destino[0] - x, destino[1] - y)
    heap = [(h(*origen), 0)]
    mejor = {}
    contador = 0
    meta = None
    while heap and contador < max_nodos:
        _, i = heapq.heappop(heap)
        x, y, z, g, az0 = X[i], Y[i], Z[i], G[i], AZ[i]
        dist = h(x, y)
        if dist <= d:                                   # intento de cierre directo al destino
            if abs(zD - z) / max(dist, 1e-9) <= pend * 1.02:
                meta = i
                break
        cand = []
        px, py = x + d * ca, y + d * sa
        zg = terr.z(px, py)
        for sg in (+1, -1):
            f = zg - (z + sg * dh)
            idx = np.where(f[:-1] * f[1:] <= 0)[0]
            for k in idx:
                t = f[k] / (f[k] - f[k + 1]) if f[k] != f[k + 1] else 0.0
                a = ang[k] + t * (ang[k + 1] - ang[k])
                nx, ny = x + d * math.cos(a), y + d * math.sin(a)
                if terr.dentro(nx, ny) and not terr.en_vacio([nx], [ny])[0][0]:
                    contra = (sg * (zD - z)) < 0
                    cand.append((nx, ny, z + sg * dh, a, w_contra * d if contra else 0.0))
        # paso directo hacia el destino cuando el terreno es más plano que la pendiente máxima
        a_d = math.atan2(destino[1] - y, destino[0] - x)
        paso = min(d, dist)
        nx, ny = x + paso * math.cos(a_d), y + paso * math.sin(a_d)
        nz = float(terr.z([nx], [ny])[0])
        if paso > 0 and abs(nz - z) / paso <= pend and terr.dentro(nx, ny) and not terr.en_vacio([nx], [ny])[0][0]:
            cand.append((nx, ny, nz, a_d, 0.0))
        for nx, ny, nz, a, pen in cand:
            paso_real = math.hypot(nx - x, ny - y)
            dd = 0.0 if az0 is None else _ang_diff(a, az0)
            ng = g + paso_real * (1 + w_defl * dd / math.pi) + pen
            key = (round(nx / (d * 0.5)), round(ny / (d * 0.5)), round(nz / dh))
            if key in mejor and mejor[key] <= ng:
                continue
            mejor[key] = ng
            X.append(nx); Y.append(ny); Z.append(nz); G.append(ng); PAR.append(i); AZ.append(a)
            heapq.heappush(heap, (ng + h(nx, ny), len(X) - 1))
            contador += 1
    if meta is None:
        raise RuntimeError("La línea de ceros no alcanza el destino con esa pendiente. "
                           "Prueba otra pendiente, una equidistancia menor o revisa origen/destino.")
    idxs = []
    j = meta
    while j != -1:
        idxs.append(j); j = PAR[j]
    idxs = idxs[::-1]
    P = np.column_stack([[X[k] for k in idxs], [Y[k] for k in idxs], [Z[k] for k in idxs]])
    P = np.vstack([P, [destino[0], destino[1], zD]])        # cierre al destino
    seg = np.hypot(*np.diff(P[:, :2], axis=0).T)
    dz = np.diff(P[:, 2])
    info = {"longitud": float(seg.sum()), "pend_media": float(np.abs(dz).sum() / seg.sum()),
            "ascensos": float(dz[dz > 0].sum()), "d_compas": d, "nodos": contador,
            "bruce": longitud_virtual(list(zip(seg, dz)))}
    return P, info


# =============================================================================================
# Ruta de costo mínimo con el costo transversal de la sección real (replanteo según la vía)
# =============================================================================================
_MOV16 = [(1, 0), (0, 1), (-1, 0), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1),
          (1, 2), (2, 1), (1, -2), (2, -1), (-1, 2), (-2, 1), (-1, -2), (-2, -1)]


def tabla_costo_transversal(st, s_max: float = 1.2, paso: float = 0.01, xmax: float = 60.0, dd: float = 0.25):
    """Costo por metro de vía (corte + relleno de la sección completa, eje a nivel del terreno) en función de la
    pendiente transversal natural s. Usa la MISMA plantilla que el cálculo de secciones (calzada, bermas,
    bombeo, cuneta, taludes, espesor de pavimento), sobre un terreno plano inclinado. inf si el talud no cierra."""
    from . import secciones as SC
    d = np.arange(0.0, xmax + dd / 2, dd)
    sg = np.arange(0.0, s_max + paso / 2, paso)
    costo = np.empty(len(sg))
    for i, s in enumerate(sg):
        up = SC._lado(d, s * d, 0.0, -st.bombeo, st.carril, st)
        dn = SC._lado(d, -s * d, 0.0, -st.bombeo, st.carril, st)
        if up["sin_chaflan"] or dn["sin_chaflan"]:
            costo[i] = np.inf
        else:
            costo[i] = st.c_corte * (up["Ac"] + dn["Ac"]) + st.c_relleno * (up["Af"] + dn["Af"])
    return sg, costo


def _interp_costo(sg, costo, s):
    fin = np.where(np.isfinite(costo), costo, 1e12)
    v = np.interp(s, sg, fin)
    v = np.where((s > sg[-1]) | (v > 1e11), np.inf, v)
    return v


def _malla_ruta(terr, h, margen_celdas: int = 3, evitar_vacios: bool = True):
    """Malla de cotas suavizada y máscara de nodos VÁLIDOS: dentro del casco del TIN y alejados `margen_celdas`
    de su borde (fuera del levantamiento la cota es una extrapolación y el suavizado distorsiona las pendientes)."""
    from scipy.ndimage import gaussian_filter, binary_erosion
    xs, ys, X, Y, Z = terr.malla(h)
    Zs = gaussian_filter(Z, 1.0, mode="nearest")
    gy, gx = np.gradient(Zs, h)
    Q = np.column_stack([X.ravel(), Y.ravel()])
    dentro = (terr.tri.find_simplex(Q) >= 0).reshape(X.shape)
    valido = binary_erosion(dentro, iterations=margen_celdas, border_value=0)
    if evitar_vacios:                                  # nunca se traza por zonas sin datos (la cota allí sería inventada)
        valido &= ~terr.en_vacio(Q[:, 0], Q[:, 1])[0].reshape(X.shape)
    return xs, ys, Zs, gx, gy, valido


HOLGURA_PENDIENTE = 1.15       # factor sobre la pendiente máxima en las aristas de la malla (ver ruta_costo_minimo)


def _grafo_seccion(Zs, gx, gy, h, pend, sg, ct, c_long, alfa, valido=None):
    from scipy.sparse import coo_matrix
    ny, nx = Zs.shape
    idx = np.arange(nx * ny).reshape(ny, nx)
    R, C, W = [], [], []
    for di, dj in _MOV16:
        rs = slice(max(0, -dj), ny - max(0, dj)); cs = slice(max(0, -di), nx - max(0, di))
        rd = slice(max(0, dj), ny - max(0, -dj)); cd = slice(max(0, di), nx - max(0, -di))
        L = h * math.hypot(di, dj)
        tx, ty = di / math.hypot(di, dj), dj / math.hypot(di, dj)
        dz = Zs[rd, cd] - Zs[rs, cs]
        ok = np.abs(dz) / L <= pend
        if valido is not None:
            ok &= valido[rs, cs] & valido[rd, cd]
        gxm = 0.5 * (gx[rs, cs] + gx[rd, cd]); gym = 0.5 * (gy[rs, cs] + gy[rd, cd])
        s_t = np.abs(gxm * (-ty) + gym * tx)                     # pendiente transversal al eje de la arista
        ct_t = _interp_costo(sg, ct, s_t)
        ok &= np.isfinite(ct_t)                                  # el talud no cierra con el terreno: no es construible (con cualquier α)
        w = L * (c_long + np.where(np.isfinite(ct_t), alfa * np.where(np.isfinite(ct_t), ct_t, 0.0), 0.0))
        R.append(idx[rs, cs][ok]); C.append(idx[rd, cd][ok]); W.append(w[ok])
    R, C, W = np.concatenate(R), np.concatenate(C), np.concatenate(W)
    return coo_matrix((W, (R, C)), shape=(nx * ny, nx * ny)).tocsr()


def ruta_costo_minimo(terr, origen, destino, pend: float, st, alfa: float = 1.0, h: float | None = None,
                      tabla=None, devolver_grafo: bool = False, holgura: float = HOLGURA_PENDIENTE):
    """Ruta de costo mínimo (Dijkstra, 16 direcciones) con |pendiente longitudinal natural| <= pend y costo por arista
    = L·(c_long + alfa·costo transversal de la sección). alfa=0: ruta más corta factible; alfa grande: prioriza
    laderas suaves. Devuelve (P Nx3, info)."""
    from scipy.sparse.csgraph import dijkstra
    ext = max(terr.xmax - terr.xmin, terr.ymax - terr.ymin)
    h = h or max(ext / 250.0, 4.0)
    sg, ct = tabla if tabla is not None else tabla_costo_transversal(st)
    xs, ys, Zs, gx, gy, valido = _malla_ruta(terr, h)
    # holgura: la malla discreta y el ruido de pendiente del TIN (sobre todo si viene de curvas de nivel) hacen inalcanzable el
    # límite exacto; la rasante posterior sigue limitada a `pend` exacto (los excesos locales se resuelven con corte/relleno)
    G = _grafo_seccion(Zs, gx, gy, h, pend * holgura, sg, ct, st.c_long, alfa, valido)

    def nodo(p, nombre):
        j = int(np.argmin(np.abs(ys - p[1]))); i = int(np.argmin(np.abs(xs - p[0])))
        if not valido[j, i]:
            raise RuntimeError(f"El {nombre} queda fuera del área levantada, en una zona sin datos o demasiado cerca de su borde "
                               f"(a menos de {3 * h:.0f} m): no hay datos de terreno confiables ahí.")
        return j * len(xs) + i
    o, d_ = nodo(origen, "origen"), nodo(destino, "destino")
    dist, pred = dijkstra(G, directed=True, indices=o, return_predecessors=True)
    if not np.isfinite(dist[d_]):
        raise RuntimeError("Sin ruta factible con esa pendiente máxima y esta sección (¿talud que no cierra con el terreno?)")
    ruta = [d_]
    while ruta[-1] != o:
        ruta.append(pred[ruta[-1]])
    ruta = ruta[::-1]
    nx = len(xs)
    px = np.array([xs[n % nx] for n in ruta]); py = np.array([ys[n // nx] for n in ruta])
    pz = terr.z(px, py)
    # extremos exactos: origen y destino reales
    px[0], py[0], px[-1], py[-1] = origen[0], origen[1], destino[0], destino[1]
    P = np.column_stack([px, py, terr.z(px, py)])
    seg = np.hypot(*np.diff(P[:, :2], axis=0).T)
    info = dict(longitud=float(seg.sum()), costo_modelo=float(dist[d_]), alfa=alfa, h=h, nodos=ruta)
    if devolver_grafo:
        info["grafo"] = G
    return P, info
