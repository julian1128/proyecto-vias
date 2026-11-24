"""
Sección transversal, peralte, sobreancho, chaflanes, volúmenes y diagrama de masas.

Método (diapositivas "Diseño de sección transversal" y "Análisis de volúmenes y movimiento de tierras"):
  * Banca: calzada (n carriles) + bermas; en curva se añade el SOBREANCHO por el lado interno, creciendo
    linealmente a lo largo de la espiral.
  * Peralte: en curvas espiralizadas el borde externo gira respecto al eje: aplanamiento (-bombeo -> 0), luego
    bombeo contrario (0 -> +bombeo) y finalmente giro conjunto hasta e; la transición ocupa la longitud de la espiral.
  * Cuneta triangular (talud interior 3:1) en corte; talud de corte desde el fondo de cuneta; talud de
    terraplén desde el borde de la berma.
  * Chaflán: intersección del talud con el terreno natural (TIN) a lo largo de la normal al eje.
  * Cota de trabajo = cota roja (rasante) - cota negra (terreno): + relleno / - corte.
  * Áreas por integración numérica entre terreno y plantilla; volúmenes por áreas extremas.
  * Diagrama de masas con factor de terraplén compacto (la ordenada acumula corte - f·relleno).
Datos marcados [VERIFICAR] son valores típicos que no vienen en las diapositivas.
"""
from __future__ import annotations
import math
from dataclasses import dataclass, field, replace

import numpy as np

from . import normas

# --------------------------------------------------------------------------- vehículos de diseño [DIAP]
VEHICULOS = {
    # tractocamión 3S2 de la tabla de las diapositivas
    "3S2": dict(A=1.22, L1=5.95, L2=0.0, L3=12.97, u=2.59),
    # vehículo del ejemplo resuelto de las diapositivas (se usa para verificar las fórmulas)
    "ejemplo_diapositiva": dict(A=2.62, L1=6.71, L2=1.89, L3=4.02, u=2.59),
}

# Pendiente de talud de corte (H por 1V) según material  [VERIFICAR: valores típicos; el manual da una tabla por
# tipo de material y altura de corte]. Terraplén 1.5H:1V por defecto (las diapositivas usan 1:1 a 2:1).
MATERIALES_CORTE = {"roca_sana": 0.25, "roca_fracturada": 0.5, "suelo_comun": 1.0, "suelo_blando": 1.5}


RELACION_TIERRAS_ALERTA = 5.0     # aviso si el costo de corte + relleno supera este múltiplo del costo de la vía


@dataclass
class SeccionTipo:
    calzada: float = 7.30           # ancho total de calzada (tabla 5.2)
    n_carriles: int = 2
    berma: float = 1.8              # por lado (tabla 5.4)
    bombeo: float = 0.02            # tabla 5.3 (superficie asfáltica o concreto)
    berma_pend: float = 0.04        # pendiente transversal de la berma en tangente  [VERIFICAR]
    e_pav: float = 0.0              # espesor de la estructura de pavimento (el movimiento de tierras se mide a subrasante)
    cuneta_ancho: float = 0.80      # borde de berma -> fondo de cuneta (m), como en el taller N°5
    cuneta_talud: float = 3.0       # talud interior de la cuneta (H por 1V) = 3:1
    cuneta_prof_sub: float | None = None   # profundidad bajo subrasante; si se da, ancho = talud·(e_pav + prof)
    tc_h: float = 1.0               # talud de corte (H por 1V)
    tf_h: float = 1.5               # talud de terraplén (H por 1V)
    vehiculo: str = "3S2"
    factor_terraplen: float = 1.25  # factor de terraplén compacto (taller N°5: 1.25-1.30)
    c_long: float = 40.0            # costo por m de vía (estructura, drenaje...) [unidades relativas]
    c_corte: float = 1.0            # costo por m³ de corte
    c_relleno: float = 1.4          # costo por m³ de relleno compactado
    c_acarreo: float = 0.002        # costo por m³·m de acarreo longitudinal  [VERIFICAR]
    c_prestamo: float = 2.0         # costo por m³ de préstamo (material faltante)
    c_botadero: float = 0.5         # costo por m³ de excedente a botadero

    @property
    def carril(self) -> float:
        return self.calzada / self.n_carriles


# --------------------------------------------------------------------------- sobreancho y visibilidad
def holgura_C(At: float) -> float:
    """Holgura lateral C (m) según el ancho de calzada en tangente (interpola 6.00/6.60/7.20) [DIAP]."""
    return float(np.interp(At, [6.0, 6.6, 7.2], [0.60, 0.75, 0.90]))


def sobreancho(R: float, V: float, At: float, n: int = 2, vehiculo="3S2") -> dict:
    """S = Ac - At ;  Ac = n·(U + C) + (n-1)·FA + Z ;  U = u + Rc - √(Rc² - (L1+L2+L3)²) ;
    FA = √(Rc² + A·(2·L1 + A)) - Rc ;  Z = 0.1·√(V/Rc)   [DIAP]. Rc = R (radio del eje), como en el ejemplo resuelto."""
    v = VEHICULOS[vehiculo] if isinstance(vehiculo, str) else vehiculo
    Rc = float(R)
    suma = v["L1"] + v["L2"] + v["L3"]
    U = v["u"] + Rc - math.sqrt(max(Rc * Rc - suma * suma, 0.0))
    FA = math.sqrt(Rc * Rc + v["A"] * (2 * v["L1"] + v["A"])) - Rc
    Z = 0.1 * math.sqrt(V / Rc)
    C = holgura_C(At)
    Ac = n * (U + C) + (n - 1) * FA + Z
    return dict(U=U, FA=FA, Z=Z, C=C, Ac=Ac, S=max(Ac - At, 0.0))


def radio_carril_interior(R: float, c: float, S: float) -> float:
    """R'c = R - S/2 - (c + S/2)/2 (radio de la trayectoria del carril interior; ejemplo de las diapositivas)."""
    return R - S / 2.0 - (c + S / 2.0) / 2.0


def flecha_M(Rc: float, Dp: float) -> float:
    """Ordenada media M = Rc·(1 - cos(28.65·Dp/Rc)) entre la trayectoria y el obstáculo para ver Dp  [DIAP]."""
    return Rc * (1.0 - math.cos(math.radians(28.65 * Dp / Rc)))


# --------------------------------------------------------------------------- peralte
def pendientes_peralte(f: float, e: float, m: float):
    """Pendientes transversales (positivo = el borde queda más alto que el eje) del carril EXTERNO y del INTERNO
    en una fracción f∈[0,1] de la transición (0 = bombeo normal, 1 = peralte pleno e):
      externo: -m -> e linealmente ;  interno: -m mientras el externo no pase de +m, luego -externo."""
    so = -m + f * (e + m)
    si = -m if so <= m else -so
    return so, si


# --------------------------------------------------------------------------- diagrama de masas
def diagrama_masas(s, Vc, Vf, f: float):
    """s: abscisas de las secciones (n+1); Vc, Vf: volúmenes de corte y relleno (compacto) por tramo (n).
    Ordenada de masa = acumulado de (corte - f·relleno). Devuelve dict con ordenadas, lóbulos cerrados,
    volumen de acarreo longitudinal (suma de máximos), distancia media de acarreo, compensación transversal,
    excedente y faltante."""
    s = np.asarray(s, float)
    Vc, Vf = np.asarray(Vc, float), np.asarray(Vf, float) * f
    neto = Vc - Vf
    M = np.concatenate([[0.0], np.cumsum(neto)])
    transversal = float(np.minimum(Vc, Vf).sum())
    lobulos, i, n = [], 0, len(M)
    # cruces por cero: se recorre cada tramo y se corta donde cambia el signo
    pts_s, pts_m = [s[0]], [M[0]]
    for k in range(n - 1):
        a, b = M[k], M[k + 1]
        if a * b < 0:
            t = a / (a - b)
            pts_s.append(s[k] + t * (s[k + 1] - s[k])); pts_m.append(0.0)
        pts_s.append(s[k + 1]); pts_m.append(b)
    pts_s, pts_m = np.array(pts_s), np.array(pts_m)
    ceros = np.where(np.isclose(pts_m, 0.0, atol=1e-9))[0]
    cerrados = 0.0
    area_cerrada = 0.0
    areas, vmax = [], []
    for a, b in zip(ceros[:-1], ceros[1:]):
        seg_s, seg_m = pts_s[a:b + 1], pts_m[a:b + 1]
        if np.all(np.abs(seg_m) < 1e-9):
            continue
        area = float(np.trapezoid(np.abs(seg_m), seg_s))
        vm = float(np.abs(seg_m).max())
        lobulos.append(dict(s0=float(seg_s[0]), s1=float(seg_s[-1]), vmax=vm, area=area,
                            dma=area / vm if vm > 0 else 0.0, signo=float(np.sign(seg_m[np.argmax(np.abs(seg_m))]))))
        cerrados += vm
        area_cerrada += area
    fin = float(M[-1])
    # la parte abierta del final (si la ordenada no vuelve a cero) es excedente (+) o faltante (-)
    if len(ceros) and ceros[-1] < len(pts_m) - 1:
        seg_s, seg_m = pts_s[ceros[-1]:], pts_m[ceros[-1]:]
        area_abierta = float(np.trapezoid(np.abs(seg_m), seg_s))
    else:
        area_abierta = 0.0
    return dict(ordenada=M, s=s, lobulos=lobulos, acarreo_long_m3=cerrados, acarreo_m3m=area_cerrada,
                dma_m=(area_cerrada / cerrados) if cerrados > 0 else 0.0, transversal_m3=transversal,
                excedente_m3=max(fin, 0.0), faltante_m3=max(-fin, 0.0), ordenada_final=fin,
                area_abierta_m3m=area_abierta)


def volumenes_areas_extremas(s, A):
    """Volumen por tramo = (A1 + A2)/2 · distancia  (áreas extremas)."""
    s, A = np.asarray(s, float), np.asarray(A, float)
    return (A[:-1] + A[1:]) / 2.0 * np.diff(s)


# --------------------------------------------------------------------------- geometría de un lado de la sección
def _lado(d, zt, zr, s_lane, cw, st: SeccionTipo):
    """Plantilla y chaflán de un lado. d: distancias al eje (0..xmax); zt: terreno en esas distancias;
    zr: rasante en el eje; s_lane: pendiente transversal del carril; cw: ancho del carril de ese lado."""
    b = st.berma
    dE = cw + b
    z_edge = zr + s_lane * cw
    sb = s_lane if abs(s_lane) >= st.berma_pend else -st.berma_pend
    zE = z_edge + sb * b                                   # borde exterior de la berma (superficie)
    zsubE = zE - st.e_pav                                  # nivel de subrasante en el borde
    # --- plataforma a nivel de subrasante
    T = np.empty_like(d)
    dentro = d <= dE
    zsup = np.where(d <= cw, zr + s_lane * d, z_edge + sb * (d - cw))
    T[dentro] = (zsup - st.e_pav)[dentro]
    ztE = float(np.interp(dE, d, zt))
    modo = "corte" if ztE >= zsubE else "relleno"
    puntos = [(0.0, zr), (cw, z_edge), (dE, zE)]
    dB = zB = None
    fuera = ~dentro
    if modo == "corte":
        if st.cuneta_prof_sub is not None:
            prof = st.e_pav + st.cuneta_prof_sub
        else:
            prof = max(st.cuneta_ancho / st.cuneta_talud, st.e_pav) if st.cuneta_ancho > 0 else st.e_pav
        zB = zE - prof
        dB = dE + st.cuneta_talud * prof
        inner = zE - (d - dE) / st.cuneta_talud
        tramo_c = fuera & (d <= dB)
        T[tramo_c] = np.minimum(zsubE, inner[tramo_c]) if prof > 0 else zsubE
        tramo_t = fuera & (d > dB)
        T[tramo_t] = zB + (d[tramo_t] - dB) / st.tc_h
        if prof > 0:
            puntos.append((dB, zB))
    else:
        T[fuera] = np.minimum(zsubE, zE - (d[fuera] - dE) / st.tf_h)
    diff = zt - T
    # --- chaflán: primer cruce (después del borde) donde el terreno alcanza el talud
    idx = np.where(fuera)[0]
    if modo == "corte":
        cruz = idx[np.where(diff[idx] <= 0)[0]]
    else:
        cruz = idx[np.where(diff[idx] >= 0)[0]]
    sin_chaflan = len(cruz) == 0
    if sin_chaflan:
        k = len(d) - 1
        dday = float(d[-1]); zday = float(zt[-1])
        dd_, df_ = d, diff
    else:
        k = int(cruz[0])
        if k > 0 and diff[k] != diff[k - 1]:
            t = diff[k - 1] / (diff[k - 1] - diff[k])
        else:
            t = 1.0
        dday = float(d[k - 1] + t * (d[k] - d[k - 1]))
        zday = float(zt[k - 1] + t * (zt[k] - zt[k - 1]))
        dd_ = np.append(d[:k], dday)
        df_ = np.append(diff[:k], 0.0)
    Ac = float(np.trapezoid(np.maximum(df_, 0.0), dd_))
    Af = float(np.trapezoid(np.maximum(-df_, 0.0), dd_))
    puntos.append((dday, zday))
    return dict(modo=modo, dE=dE, zE=zE, d_chaflan=dday, z_chaflan=zday, Ac=Ac, Af=Af,
                sin_chaflan=sin_chaflan, puntos=puntos, dB=dB, zB=zB)


# --------------------------------------------------------------------------- cálculo a lo largo del eje
@dataclass
class Secciones:
    s: np.ndarray
    x: np.ndarray
    y: np.ndarray
    zr: np.ndarray
    zt_eje: np.ndarray
    cota_trabajo: np.ndarray          # cota roja - cota negra (+ relleno / - corte)
    Ac: np.ndarray
    Af: np.ndarray
    izq: list
    der: list
    nrm: np.ndarray
    peralte_ext: np.ndarray           # pendiente del carril externo (+ sube hacia afuera)
    sobreancho: np.ndarray            # sobreancho aplicado (m) por estación
    lado_sobreancho: np.ndarray       # +1 izquierdo, -1 derecho, 0 ninguno
    tipo: list
    vol_corte: np.ndarray
    vol_relleno: np.ndarray
    masas: dict
    resumen: dict
    avisos: list
    visibilidad: list
    st: SeccionTipo
    perfiles: dict = field(default_factory=dict)       # terreno por estación (para dibujar)

    def puntos_seccion(self, k):
        """Polilínea de la plantilla (izq→der) de la estación k: [(offset, cota)], offset < 0 a la izquierda."""
        izq = [(-d, z) for d, z in self.izq[k]["puntos"]][::-1]
        der = [(d, z) for d, z in self.der[k]["puntos"]]
        return izq + der[1:]            # el punto del eje (0, zr) aparece una sola vez


def _estado_transversal(curvas_info, s, m):
    """Pendiente del carril izquierdo/derecho y sobreancho (lado, valor) en la abscisa s."""
    sL = sR = -m
    sob, lado = 0.0, 0
    for c in curvas_info:
        a0, Le, Lc = c["a0"], c["Le"], c["Lc"]
        a1, a2, a3 = a0 + Le, a0 + Le + Lc, a0 + 2 * Le + Lc
        if not (a0 <= s <= a3) or Le <= 0:
            continue
        f = (s - a0) / Le if s < a1 else (1.0 if s <= a2 else (a3 - s) / Le)
        so, si = pendientes_peralte(f, c["e"], m)
        if c["sgn"] > 0:                # curva a la izquierda: centro a la izquierda -> interno = izquierdo, externo = derecho
            sL, sR = si, so; lado = +1
        else:
            sL, sR = so, si; lado = -1
        sob = c["S"] * f
    return sL, sR, sob, lado


@dataclass
class Preparacion:
    """Todo lo que NO depende de la rasante: transectos de terreno, estado de peralte y sobreancho por estación.
    Permite evaluar muchas cotas de rasante sin volver a muestrear el TIN (lo usa el optimizador)."""
    s: np.ndarray
    x: np.ndarray
    y: np.ndarray
    nrm: np.ndarray
    off: np.ndarray
    ZL: np.ndarray
    ZR: np.ndarray
    sL: np.ndarray
    sR: np.ndarray
    cL: np.ndarray
    cR: np.ndarray
    sob: np.ndarray
    lado: np.ndarray
    curvas_info: list
    st: SeccionTipo
    dd: float


def preparar(terr, al, st: SeccionTipo, ds: float = 10.0, xmax: float = 60.0, dd: float = 0.25) -> Preparacion:
    if not getattr(al, "valido", True):
        raise ValueError("El alineamiento es inválido; no se pueden calcular secciones.")
    s = np.unique(np.append(np.arange(0.0, al.longitud, ds), al.longitud))
    D = al.muestrear(1.0)
    sd = D[:, 2] * (al.longitud / D[-1, 2])
    x = np.interp(s, sd, D[:, 0]); y = np.interp(s, sd, D[:, 1])
    tx = np.interp(s, sd, np.gradient(D[:, 0], sd)); ty = np.interp(s, sd, np.gradient(D[:, 1], sd))
    nt = np.hypot(tx, ty); tx, ty = tx / nt, ty / nt
    nrm = np.column_stack([-ty, tx])                                    # normal izquierda
    off = np.arange(0.0, xmax + dd / 2, dd)
    XL = x[:, None] + nrm[:, 0, None] * off[None, :]; YL = y[:, None] + nrm[:, 1, None] * off[None, :]
    XR = x[:, None] - nrm[:, 0, None] * off[None, :]; YR = y[:, None] - nrm[:, 1, None] * off[None, :]
    ZL = terr.z(XL.ravel(), YL.ravel()).reshape(XL.shape)
    ZR = terr.z(XR.ravel(), YR.ravel()).reshape(XR.shape)
    curvas_info = []
    for c in al.curvas:
        sob = sobreancho(c["R"], c["V"], st.calzada, st.n_carriles, st.vehiculo)
        curvas_info.append(dict(a0=c["s_ini"], Le=c["Le"], Lc=c["Lc"], e=c["e"], sgn=c["sgn"], S=sob["S"],
                                R=c["R"], V=c["V"], idx=c["idx"], sob=sob))
    n = len(s)
    sL = np.zeros(n); sR = np.zeros(n); sobv = np.zeros(n); lado = np.zeros(n)
    for k in range(n):
        sL[k], sR[k], sobv[k], lado[k] = _estado_transversal(curvas_info, s[k], st.bombeo)
    cL = st.carril + np.where(lado > 0, sobv, 0.0)
    cR = st.carril + np.where(lado < 0, sobv, 0.0)
    return Preparacion(s=s, x=x, y=y, nrm=nrm, off=off, ZL=ZL, ZR=ZR, sL=sL, sR=sR, cL=cL, cR=cR, sob=sobv,
                       lado=lado, curvas_info=curvas_info, st=st, dd=dd)


def seccion_en(prep: Preparacion, k: int, zr: float):
    """Plantilla, chaflanes y áreas de la estación k para una rasante `zr` en el eje (izquierda, derecha)."""
    li = _lado(prep.off, prep.ZL[k], zr, prep.sL[k], prep.cL[k], prep.st)
    ld = _lado(prep.off, prep.ZR[k], zr, prep.sR[k], prep.cR[k], prep.st)
    return li, ld


def tabla_areas(prep: Preparacion, zt_eje, offsets):
    """Áreas de corte y relleno (n estaciones × m cotas) para rasantes = terreno del eje + offsets."""
    n, m = len(prep.s), len(offsets)
    Ac = np.zeros((n, m)); Af = np.zeros((n, m))
    for k in range(n):
        for j, h in enumerate(offsets):
            li, ld = seccion_en(prep, k, zt_eje[k] + h)
            Ac[k, j] = li["Ac"] + ld["Ac"]; Af[k, j] = li["Af"] + ld["Af"]
    return Ac, Af


def calcular(terr, al, pf, st: SeccionTipo, ds: float = 10.0, xmax: float = 60.0, dd: float = 0.25,
             prep: Preparacion | None = None) -> Secciones:
    """Secciones transversales cada `ds` m a lo largo del eje con la rasante `pf` (Perfil)."""
    prep = prep or preparar(terr, al, st, ds, xmax, dd)
    s, x, y, nrm, off, ZL, ZR = prep.s, prep.x, prep.y, prep.nrm, prep.off, prep.ZL, prep.ZR
    dd = prep.dd
    zr = np.interp(s, pf.s, pf.zd)
    n = len(s)
    Ac = np.zeros(n); Af = np.zeros(n)
    izq, der, tipo = [], [], []
    per_ext = np.zeros(n)
    avisos = []
    n_sin = 0
    n_ladera = 0
    for k in range(n):
        li, ld = seccion_en(prep, k, zr[k])
        izq.append(li); der.append(ld)
        Ac[k] = li["Ac"] + ld["Ac"]; Af[k] = li["Af"] + ld["Af"]
        tipo.append("corte" if (li["modo"] == "corte" and ld["modo"] == "corte")
                    else "terraplén" if (li["modo"] == "relleno" and ld["modo"] == "relleno") else "mixta")
        per_ext[k] = max(prep.sL[k], prep.sR[k])
        n_sin += int(li["sin_chaflan"]) + int(ld["sin_chaflan"])
        # terraplén sobre ladera > 20 %: requiere escalones (diapositivas)
        pend_t = abs(ZL[k][int(5 / dd)] - ZR[k][int(5 / dd)]) / 10.0
        if (li["Af"] + ld["Af"]) > 0.5 and pend_t > 0.20:
            n_ladera += 1
    sob_arr, lado_arr = prep.sob, prep.lado
    curvas_info = prep.curvas_info
    cota_trab = zr - np.interp(s, pf.s, pf.zt)
    vc = volumenes_areas_extremas(s, Ac)
    vf = volumenes_areas_extremas(s, Af)
    masas = diagrama_masas(s, vc, vf, st.factor_terraplen)
    L = float(s[-1])
    costo = dict(via=st.c_long * L, corte=st.c_corte * vc.sum(), relleno=st.c_relleno * vf.sum(),
                 acarreo=st.c_acarreo * masas["acarreo_m3m"], prestamo=st.c_prestamo * masas["faltante_m3"],
                 botadero=st.c_botadero * masas["excedente_m3"])
    costo["total"] = float(sum(costo.values()))
    tierras = costo["corte"] + costo["relleno"]
    if tierras > RELACION_TIERRAS_ALERTA * costo["via"]:
        avisos.append(f"Movimiento de tierras (corte + relleno) = {tierras / costo['via']:.1f} veces el costo de la vía: "
                      "terreno muy empinado para estas especificaciones (¿curvas de retorno, muros u otra ruta?)")
    afect = np.array([a["d_chaflan"] + b["d_chaflan"] for a, b in zip(izq, der)])
    resumen = dict(longitud_m=L, estaciones=n, corte_m3=float(vc.sum()), relleno_compacto_m3=float(vf.sum()),
                   relleno_equiv_m3=float(vf.sum() * st.factor_terraplen),
                   transversal_m3=masas["transversal_m3"], acarreo_long_m3=masas["acarreo_long_m3"],
                   acarreo_m3m=masas["acarreo_m3m"], dma_m=masas["dma_m"], excedente_m3=masas["excedente_m3"],
                   faltante_m3=masas["faltante_m3"], ancho_afectacion_medio_m=float(afect.mean()),
                   ancho_afectacion_max_m=float(afect.max()), costo=costo,
                   estaciones_sin_chaflan=n_sin, cortes=tipo.count("corte"), terraplenes=tipo.count("terraplén"),
                   mixtas=tipo.count("mixta"), sobreancho_max_m=float(sob_arr.max()))
    if n_sin:
        avisos.append(f"{n_sin} chaflán(es) sin intersección dentro de {off[-1]:.0f} m: talud sin cerrar con el terreno "
                      "(¿talud más tendido o muro?)")
    if n_ladera:
        avisos.append(f"{n_ladera} sección(es) con terraplén sobre ladera > 20 %: se requieren escalones/obras especiales")
    # cobertura: el chaflán (y el eje) deben quedar dentro del área levantada; fuera la cota es una extrapolación
    pts_ch = np.vstack([np.column_stack([x + nrm[:, 0] * np.array([l["d_chaflan"] for l in izq]),
                                         y + nrm[:, 1] * np.array([l["d_chaflan"] for l in izq])]),
                        np.column_stack([x - nrm[:, 0] * np.array([l["d_chaflan"] for l in der]),
                                         y - nrm[:, 1] * np.array([l["d_chaflan"] for l in der])])])
    fuera = int((terr.tri.find_simplex(pts_ch) < 0).sum())
    if fuera:
        avisos.append(f"{fuera} chaflán(es) quedan fuera del área levantada: su cota es extrapolada (levantar más terreno)")
    # zonas sin datos: el TIN interpola con triángulos enormes y la cota allí no es medida
    vac, umbral_vac = terr.en_vacio(x, y)
    largo = float(np.sum(np.where(vac[:-1] | vac[1:], np.diff(s), 0.0))) if vac.any() else 0.0
    resumen["eje_en_vacio_m"] = largo
    if vac.any():
        avisos.append(f"El eje atraviesa ≈{largo:.0f} m ({100 * largo / s[-1]:.0f} % de su longitud) de zonas sin datos del "
                      f"levantamiento (a más de {umbral_vac:.0f} m de cualquier dato): la cota ahí es interpolada, no medida")
    visibilidad = _visibilidad(al, curvas_info, st, s, izq, der, off, avisos)
    return Secciones(s=s, x=x, y=y, zr=zr, zt_eje=np.interp(s, pf.s, pf.zt), cota_trabajo=cota_trab, Ac=Ac, Af=Af,
                     izq=izq, der=der, nrm=nrm, peralte_ext=per_ext, sobreancho=sob_arr, lado_sobreancho=lado_arr,
                     tipo=tipo, vol_corte=vc, vol_relleno=vf, masas=masas, resumen=resumen, avisos=avisos,
                     visibilidad=visibilidad, st=st, perfiles={"off": off, "ZL": ZL, "ZR": ZR})


def _visibilidad(al, curvas_info, st, s, izq, der, off, avisos, h_vis: float = 0.84):
    """Criterio simplificado de visibilidad de parada en curvas con corte en el lado interno: la línea de visión
    pasa a h_vis (≈ semisuma de ojo 1.08 y obstáculo 0.60) sobre el pavimento. M' = distancia entre la
    trayectoria del carril interior y el talud de corte a esa altura. [VERIFICAR contra el ejemplo del manual]."""
    out = []
    for c in curvas_info:
        centro = c["a0"] + c["Le"] + c["Lc"] / 2.0
        k = int(np.argmin(np.abs(s - centro)))
        Dp = normas.dp_tabla(normas.v_snap(c["V"]))
        Rc = radio_carril_interior(c["R"], st.carril, c["sob"]["S"])
        M = flecha_M(Rc, Dp)
        lado = izq[k] if c["sgn"] > 0 else der[k]
        info = dict(PI=c["idx"], R=c["R"], V=c["V"], Dp=Dp, M_req=M, S=c["sob"]["S"], Rc=Rc, M_disp=None, ok=True, corrimiento=0.0)
        if lado["modo"] == "corte" and lado["dB"] is not None:
            d_lane = (st.carril + c["sob"]["S"] / 2.0) / 2.0
            d_obst = lado["dB"] + st.tc_h * (h_vis + (lado["zE"] - lado["zB"]))
            M_disp = d_obst - d_lane
            info.update(M_disp=M_disp, ok=M_disp >= M, corrimiento=max(M - M_disp, 0.0))
            if not info["ok"]:
                avisos.append(f"Curva PI-{c['idx']}: visibilidad de parada limitada por el talud de corte; "
                              f"corrida del talud ≈ {info['corrimiento']:.1f} m (M req. {M:.1f} m, disp. {M_disp:.1f} m)")
        out.append(info)
    return out


# --------------------------------------------------------------------------- salidas
def tabla(sec: Secciones):
    filas = []
    for k in range(len(sec.s)):
        filas.append(dict(abscisa=sec.s[k], cota_terreno=sec.zt_eje[k], cota_rasante=sec.zr[k],
                          cota_trabajo=sec.cota_trabajo[k], tipo=sec.tipo[k], area_corte_m2=sec.Ac[k],
                          area_relleno_m2=sec.Af[k], chaflan_izq_m=sec.izq[k]["d_chaflan"],
                          chaflan_der_m=sec.der[k]["d_chaflan"], peralte_ext_pct=100 * sec.peralte_ext[k],
                          sobreancho_m=sec.sobreancho[k], masa_acumulada_m3=sec.masas["ordenada"][k]))
    return filas


def graficar(sec: Secciones, ruta: str, n_muestras: int = 4):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(14, 11))
    gs = fig.add_gridspec(4, 2, height_ratios=[0.9, 0.9, 1.0, 1.1])
    ks = []
    ks.append(int(np.argmax(sec.Ac)))                                  # mayor corte
    ks.append(int(np.argmax(sec.Af)))                                  # mayor relleno
    mix = [i for i, t in enumerate(sec.tipo) if t == "mixta"]
    ks.append(mix[len(mix) // 2] if mix else len(sec.s) // 2)
    cur = np.where(sec.sobreancho > 0)[0]
    ks.append(int(cur[len(cur) // 2]) if len(cur) else len(sec.s) // 3)
    titulos = ["Mayor área de corte", "Mayor área de relleno", "Sección mixta", "En curva (peralte y sobreancho)"]
    off, ZL, ZR = sec.perfiles["off"], sec.perfiles["ZL"], sec.perfiles["ZR"]
    for j, (k, tit) in enumerate(zip(ks, titulos)):
        ax = fig.add_subplot(gs[j // 2, j % 2])
        dL = sec.izq[k]["d_chaflan"]; dR = sec.der[k]["d_chaflan"]
        mL, mR = off <= dL + 6, off <= dR + 6
        ax.plot(-off[mL], ZL[k][mL], color="green", lw=1.4, label="Terreno")
        ax.plot(off[mR], ZR[k][mR], color="green", lw=1.4)
        P = sec.puntos_seccion(k)
        ax.plot([p[0] for p in P], [p[1] for p in P], "r-", lw=2, label="Plantilla (rasante, cuneta, taludes)")
        ax.plot([-dL, dR], [sec.izq[k]["z_chaflan"], sec.der[k]["z_chaflan"]], "ko", ms=5, label="Chaflanes")
        ax.set_title(f"K{sec.s[k] / 1000:.3f} · {tit}\nAc={sec.Ac[k]:.1f} m²  Af={sec.Af[k]:.1f} m²  "
                     f"peralte={100 * sec.peralte_ext[k]:+.1f} %  sobreancho={sec.sobreancho[k]:.2f} m", fontsize=9)
        ax.set_aspect("equal"); ax.grid(alpha=0.3)
        if j == 0:
            ax.legend(fontsize=7, loc="upper left")
    ax = fig.add_subplot(gs[2, :])
    ax.fill_between(sec.s, 0, sec.Ac, color="red", alpha=0.5, label="Área de corte")
    ax.fill_between(sec.s, 0, -sec.Af, color="blue", alpha=0.5, label="Área de relleno")
    ax.axhline(0, color="k", lw=0.6); ax.set_ylabel("m²"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title("Áreas de las secciones a lo largo del eje")
    ax = fig.add_subplot(gs[3, :], sharex=ax)
    ax.plot(sec.masas["s"], sec.masas["ordenada"], "k-", lw=1.8)
    ax.axhline(0, color="gray", lw=0.8)
    for lb in sec.masas["lobulos"]:
        ax.axvspan(lb["s0"], lb["s1"], color="orange" if lb["signo"] > 0 else "skyblue", alpha=0.2)
    r = sec.resumen
    ax.set_title(f"Diagrama de masas (f = {sec.st.factor_terraplen:.2f}) · acarreo longitudinal {r['acarreo_long_m3']:.0f} m³ · "
                 f"DMA {r['dma_m']:.0f} m · excedente {r['excedente_m3']:.0f} m³ · faltante {r['faltante_m3']:.0f} m³", fontsize=9)
    ax.set_xlabel("Abscisa (m)"); ax.set_ylabel("Ordenada (m³)"); ax.grid(alpha=0.3)
    plt.tight_layout(); fig.savefig(ruta, dpi=95); plt.close(fig)


def bordes(sec: Secciones):
    """Polilíneas en planta de la banca real (con sobreancho): borde de calzada, borde de berma, fondo de cuneta
    y chaflanes, por lado. Los chaflanes se parten donde cambia entre corte y relleno.
    Devuelve dict de listas de arrays Nx2."""
    c = sec.st.carril
    cL = c + np.where(sec.lado_sobreancho > 0, sec.sobreancho, 0.0)
    cR = c + np.where(sec.lado_sobreancho < 0, sec.sobreancho, 0.0)
    P = np.column_stack([sec.x, sec.y]); N = sec.nrm
    out = {"calzada_izq": P + N * cL[:, None], "calzada_der": P - N * cR[:, None],
           "berma_izq": P + N * np.array([l["dE"] for l in sec.izq])[:, None],
           "berma_der": P - N * np.array([l["dE"] for l in sec.der])[:, None]}
    chaf = {"corte": [], "relleno": []}
    for lado_, signo in ((sec.izq, 1.0), (sec.der, -1.0)):
        pts = [P[k] + signo * N[k] * lado_[k]["d_chaflan"] for k in range(len(P))]
        modos = [l["modo"] for l in lado_]
        ini = 0
        for k in range(1, len(P) + 1):
            if k == len(P) or modos[k] != modos[ini]:
                seg = np.array(pts[max(ini - 1, 0) if ini > 0 else 0:k])
                if len(seg) >= 2:
                    chaf[modos[ini]].append(seg)
                ini = k
    out["chaflan"] = chaf
    return out
