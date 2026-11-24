"""Nube de puntos -> TIN -> curvas de nivel. Clasificación del terreno (INVÍAS)."""
from __future__ import annotations
import math
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

from . import normas

_ALIAS = {"x": {"x", "este", "east", "easting", "e"},
          "y": {"y", "norte", "north", "northing", "n"},
          "z": {"z", "cota", "elevacion", "elevación", "altura", "elev", "h"}}


def cargar_puntos(ruta: str, hoja=0) -> np.ndarray:
    """XLSX / CSV / TXT (X Y Z o Este Norte Cota) o LAS/LAZ. Devuelve array (n,3)."""
    ext = ruta.lower().rsplit(".", 1)[-1]
    if ext in ("xlsx", "xlsm", "xls"):
        hoja = int(hoja) if str(hoja).isdigit() else hoja
        df = pd.read_excel(ruta, sheet_name=hoja)
        cols = {k: next((c for c in df.columns if str(c).strip().lower() in v), None) for k, v in _ALIAS.items()}
        sub = df[[cols["x"], cols["y"], cols["z"]]] if None not in cols.values() else df.iloc[:, :3]
        return sub.apply(pd.to_numeric, errors="coerce").dropna().to_numpy(float)
    if ext in ("las", "laz"):
        import laspy
        las = laspy.read(ruta)
        return np.column_stack([las.x, las.y, las.z])
    # CSV / TXT: separador detectado, con o sin encabezado (columnas por nombre si coinciden con los alias)
    df = pd.read_csv(ruta, sep=None, engine="python", header=0)
    cols = {k: next((c for c in df.columns if str(c).strip().lower() in v), None) for k, v in _ALIAS.items()}
    if None in cols.values():                          # sin encabezado reconocible: las tres primeras columnas
        df = pd.read_csv(ruta, sep=None, engine="python", header=None)
        sub = df.iloc[:, :3]
    else:
        sub = df[[cols["x"], cols["y"], cols["z"]]]
    return sub.apply(pd.to_numeric, errors="coerce").dropna().to_numpy(float)


def adelgazar(pts: np.ndarray, celda: float) -> np.ndarray:
    """Promedia puntos por celda (reduce ruido y tamaño del TIN)."""
    x0, y0 = pts[:, 0].min(), pts[:, 1].min()
    k = (np.floor((pts[:, 0] - x0) / celda).astype(np.int64) * 1_000_003
         + np.floor((pts[:, 1] - y0) / celda).astype(np.int64))
    df = pd.DataFrame(pts, columns=["x", "y", "z"]); df["k"] = k
    return df.groupby("k").mean().to_numpy()


class Terreno:
    """TIN de Delaunay con interpolación lineal (como la superficie de Civil 3D)."""

    def __init__(self, pts: np.ndarray, celda_tin: float | None = 2.0):
        pts = np.asarray(pts, float)
        self.d_vacio = self._umbral_vacio(pts)               # calibrado con los datos ORIGINALES, antes de adelgazar
        self._arbol = None
        self.pts = adelgazar(pts, celda_tin) if celda_tin else pts
        self.tri = Delaunay(self.pts[:, :2])
        self._lin = LinearNDInterpolator(self.tri, self.pts[:, 2])
        self._near = NearestNDInterpolator(self.pts[:, :2], self.pts[:, 2])
        self.xmin, self.ymin = self.pts[:, 0].min(), self.pts[:, 1].min()
        self.xmax, self.ymax = self.pts[:, 0].max(), self.pts[:, 1].max()

    def z(self, x, y):
        x, y = np.asarray(x, float), np.asarray(y, float)
        v = np.atleast_1d(self._lin(x, y)).astype(float)
        m = np.isnan(v)
        if m.any():
            v[m] = self._near(np.atleast_1d(x)[m], np.atleast_1d(y)[m])
        return v

    @staticmethod
    def _umbral_vacio(pts: np.ndarray, minimo: float = 50.0, factor: float = 5.0) -> float:
        """Distancia (m) a partir de la cual un punto se considera SIN DATOS: el mayor entre 50 m y 5 veces el espaciado
        mediano al vecino más cercano de los datos originales (muestra de hasta 20.000 puntos)."""
        from scipy.spatial import cKDTree
        xy = pts[:, :2]
        if len(xy) > 20000:
            xy = xy[np.random.default_rng(0).choice(len(xy), 20000, replace=False)]
        d, _ = cKDTree(pts[:, :2]).query(xy, k=2)
        return max(minimo, factor * float(np.median(d[:, 1])))

    def en_vacio(self, x, y, d_max: float | None = None):
        """True donde (x, y) está a más de `d_max` m (por defecto `self.d_vacio`) de cualquier dato medido: ahí la cota es
        solo una interpolación (o una extrapolación). Sirve igual para nube de puntos y para curvas de nivel, donde la
        longitud de las aristas del TIN no distingue un hueco de una zona llana entre curvas. Devuelve (máscara, umbral)."""
        from scipy.spatial import cKDTree
        if self._arbol is None:
            self._arbol = cKDTree(self.pts[:, :2])
        umbral = float(d_max) if d_max is not None else self.d_vacio
        d, _ = self._arbol.query(np.column_stack([np.atleast_1d(x), np.atleast_1d(y)]))
        return d > umbral, umbral

    def dentro(self, x, y, margen=0.0):
        return (self.xmin + margen <= x <= self.xmax - margen) and (self.ymin + margen <= y <= self.ymax - margen)

    def malla(self, h: float):
        xs = np.arange(self.xmin, self.xmax + h / 2, h); ys = np.arange(self.ymin, self.ymax + h / 2, h)
        X, Y = np.meshgrid(xs, ys)
        return xs, ys, X, Y, self.z(X.ravel(), Y.ravel()).reshape(X.shape)


def pendiente_transversal_predominante(terr: Terreno, h: float = 5.0) -> float:
    """Pendiente natural predominante (%) = mediana del módulo del gradiente sobre una malla, SOLO en las celdas con datos
    (a menos de `terr.d_vacio` de una medición): las zonas sin datos se interpolan con triángulos planos enormes y falsean
    la mediana (con curvas de nivel reales daba 11,5 % cuando en la zona levantada es 23 %)."""
    from scipy.ndimage import gaussian_filter
    _, _, X, Y, Z = terr.malla(h)
    gy, gx = np.gradient(gaussian_filter(Z, 1.0, mode="nearest"), h)
    con_datos = ~terr.en_vacio(X.ravel(), Y.ravel())[0].reshape(X.shape)
    g = np.hypot(gx, gy)
    if con_datos.sum() < 50:                           # casi nada cubierto: se usa toda la malla
        return float(np.median(g) * 100.0)
    return float(np.median(g[con_datos]) * 100.0)


def clasificar(terr: Terreno) -> tuple[str, float]:
    p = pendiente_transversal_predominante(terr)
    return normas.clasificar_terreno(p), p


# ------------------------------------------------------------------ curvas de nivel
def chaikin(P: np.ndarray, it: int = 2) -> np.ndarray:
    for _ in range(it):
        if len(P) < 3:
            break
        Q, R = 0.75 * P[:-1] + 0.25 * P[1:], 0.25 * P[:-1] + 0.75 * P[1:]
        P = np.vstack([P[:1], np.column_stack([Q, R]).reshape(-1, 2), P[-1:]])
    return P


def curvas_de_nivel(terr: Terreno, intervalo: float, max_arista: float = 40.0):
    """Curvas exactas sobre el TIN. Devuelve ([(cota, polilínea Nx2)], triangulación mpl)."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import tri
    s, xy = terr.tri.simplices, terr.pts[:, :2]
    a = xy[s]
    arista = np.max(np.linalg.norm(a - np.roll(a, 1, axis=1), axis=2), axis=1)
    tr = tri.Triangulation(xy[:, 0], xy[:, 1], s, mask=arista > max_arista)
    z = terr.pts[:, 2]
    niveles = np.arange(math.ceil(z.min() / intervalo) * intervalo, z.max(), intervalo)
    fig, ax = plt.subplots()
    cs = ax.tricontour(tr, z, niveles)
    out = [(float(lv), chaikin(np.asarray(sg))) for lv, segs in zip(cs.levels, cs.allsegs)
           for sg in segs if len(sg) >= 4]
    plt.close(fig)
    return out, tr
