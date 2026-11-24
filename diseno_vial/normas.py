"""
Criterios normativos INVÍAS (Manual de Diseño Geométrico de Carreteras) y AASHTO Green Book.

FUENTE de cada dato:
  [DIAP]  leído de las diapositivas del curso (tablas 4.1, 4.4, 5.2, 5.4, casos de VECH, entretangencias)
  [AASHTO] fórmula o valor del Green Book (f_max métrico, Rmín = V²/127(e+f))
  [VERIFICAR] valor recordado/estimado: confirmar contra el manual antes de usar en un proyecto real
"""
from __future__ import annotations
import math

CATEGORIAS = ("primaria_2_calzadas", "primaria_1_calzada", "secundaria", "terciaria")
TERRENOS = ("plano", "ondulado", "montañoso", "escarpado")

# ---------------------------------------------------------------- Tabla 5.2 [DIAP]
# ANCHO_CALZADA[categoria][terreno] = {VTR: ancho_m}. Las claves existentes son
# las combinaciones categoría-terreno-VTR permitidas (los "-" de la tabla no existen).
def _fila(**kw):
    return {int(k[1:]): v for k, v in kw.items()}

ANCHO_CALZADA = {
    "primaria_2_calzadas": {
        "plano":      _fila(v80=7.30, v90=7.30, v100=7.30, v110=7.30),
        "ondulado":   _fila(v80=7.30, v90=7.30, v100=7.30, v110=7.30),
        "montañoso":  _fila(v70=7.30, v80=7.30, v90=7.30, v100=7.30),
        "escarpado":  _fila(v70=7.30, v80=7.30, v90=7.30)},
    "primaria_1_calzada": {
        "plano":      _fila(v80=7.30, v90=7.30, v100=7.30),
        "ondulado":   _fila(v70=7.30, v80=7.30, v90=7.30, v100=7.30),
        "montañoso":  _fila(v60=7.30, v70=7.30, v80=7.30, v90=7.30),
        "escarpado":  _fila(v60=7.00, v70=7.00, v80=7.00)},
    "secundaria": {
        "plano":      _fila(v60=7.30, v70=7.30, v80=7.30),
        "ondulado":   _fila(v50=7.00, v60=7.30, v70=7.30, v80=7.30),
        "montañoso":  _fila(v40=6.60, v50=7.00, v60=7.00, v70=7.00),
        "escarpado":  _fila(v40=6.00, v50=6.60, v60=7.00)},
    "terciaria": {
        "plano":      _fila(v40=6.00),
        "ondulado":   _fila(v30=6.00, v40=6.00),
        "montañoso":  _fila(v20=6.00, v30=6.00, v40=6.00),
        "escarpado":  _fila(v20=6.00, v30=6.00)},
}

# ---------------------------------------------------------------- Tabla 5.4 [DIAP]
# Berma (m). Dos calzadas: (derecha, izquierda). Terciaria: berma-cuneta.
BERMA = {
    "primaria_2_calzadas": {
        "plano":     _fila(v80=(2.5, 1.0), v90=(2.5, 1.0), v100=(2.5, 1.0), v110=(2.5, 1.0)),
        "ondulado":  _fila(v80=(2.0, 1.0), v90=(2.0, 1.0), v100=(2.5, 1.0), v110=(2.5, 1.0)),
        "montañoso": _fila(v70=(1.8, 0.5), v80=(1.8, 0.5), v90=(1.8, 0.5), v100=(2.0, 1.0)),
        "escarpado": _fila(v70=(1.8, 0.5), v80=(1.8, 0.5), v90=(1.8, 0.5))},
    "primaria_1_calzada": {
        "plano":     _fila(v80=2.0, v90=2.0, v100=2.5),
        "ondulado":  _fila(v70=1.8, v80=2.0, v90=2.0, v100=2.5),
        "montañoso": _fila(v60=1.5, v70=1.5, v80=1.8, v90=1.8),
        "escarpado": _fila(v60=1.5, v70=1.5, v80=1.8)},
    "secundaria": {
        "plano":     _fila(v60=1.0, v70=1.5, v80=1.8),
        "ondulado":  _fila(v50=1.0, v60=1.0, v70=1.5, v80=1.8),
        "montañoso": _fila(v40=0.5, v50=0.5, v60=1.0, v70=1.0),
        "escarpado": _fila(v40=0.5, v50=0.5, v60=0.5)},
    "terciaria": {
        "plano":     _fila(v40=1.0),
        "ondulado":  _fila(v30=0.5, v40=1.0),
        "montañoso": _fila(v20=0.5, v30=0.5, v40=0.5),
        "escarpado": _fila(v20=0.5, v30=0.5, v40=0.5)},
}

# ---------------------------------------------------------------- Tabla 4.1 [DIAP]
# Pendiente media máxima del corredor de ruta (%) según VTR.
PEND_MEDIA_MAX = {
    "primaria_2_calzadas": {70: 6, 80: 6, 90: 6, 100: 5, 110: 5},
    "primaria_1_calzada":  {60: 7, 70: 7, 80: 6, 90: 6, 100: 5},
    "secundaria":          {40: 7, 50: 7, 60: 7, 70: 7, 80: 6},
    "terciaria":           {20: 7, 30: 7, 40: 7},
}

# ---------------------------------------------------------------- Fricción lateral [AASHTO] (métrico, Green Book)
F_MAX = {20: 0.18, 30: 0.28, 40: 0.23, 50: 0.19, 60: 0.17, 70: 0.15,
         80: 0.14, 90: 0.13, 100: 0.12, 110: 0.11, 120: 0.09, 130: 0.08}

# Rampa relativa máxima de peralte Δ (%), por velocidad  [VERIFICAR contra tabla 3.x del manual]
RAMPA_PERALTE = {20: 0.80, 30: 0.75, 40: 0.70, 50: 0.65, 60: 0.60, 70: 0.55,
                 80: 0.50, 90: 0.47, 100: 0.44, 110: 0.41, 120: 0.38}

# Distancia de visibilidad de parada: calculada con la fórmula de las diapositivas
def dist_parada(V: float, i: float = 0.0, a: float = 3.4) -> float:
    """Dp = 0.695·V + V²/(254·(a/9.8 ± i))   [DIAP, diseño vertical]. i en m/m (+ ascenso)."""
    return 0.695 * V + V * V / (254.0 * (a / 9.8 + i))


def dp_tabla(V: int) -> float:
    """Distancia de parada de la tabla 4.4 de las diapositivas [DIAP]."""
    return {20: 20, 30: 35, 40: 50, 50: 65, 60: 85, 70: 105, 80: 130, 90: 160,
            100: 185, 110: 220, 120: 250, 130: 285}[int(V)]


# ---------------------------------------------------------------- Selección de parámetros
def clasificar_terreno(pend_transversal_pct: float) -> str:
    """Por pendiente transversal predominante (%). Límites INVÍAS [VERIFICAR]; el ejemplo
    de las diapositivas (85 % -> escarpado) es consistente."""
    p = pend_transversal_pct
    return "plano" if p <= 5 else "ondulado" if p <= 25 else "montañoso" if p <= 75 else "escarpado"


def rango_vtr(categoria: str, terreno: str) -> list[int]:
    """Velocidades de tramo homogéneo permitidas (tabla 5.2: celdas con valor)."""
    return sorted(ANCHO_CALZADA[categoria][terreno])


def ancho_calzada(categoria: str, terreno: str, vtr: int) -> float:
    try:
        return ANCHO_CALZADA[categoria][terreno][vtr]
    except KeyError:
        raise ValueError(f"VTR {vtr} km/h no permitida para {categoria}/{terreno}. "
                         f"Rango: {rango_vtr(categoria, terreno)}")


def berma(categoria: str, terreno: str, vtr: int):
    return BERMA[categoria][terreno][vtr]


def pend_media_max(categoria: str, vtr: int) -> float:
    """% -> m/m. Tabla 4.1 [DIAP]."""
    t = PEND_MEDIA_MAX[categoria]
    if vtr not in t:
        raise ValueError(f"Sin pendiente máx. en tabla 4.1 para {categoria} a {vtr} km/h")
    return t[vtr] / 100.0


# ---------------------------------------------------------------- Curvas horizontales
def r_min(V: float, e_max: float = 0.08) -> float:
    """Rmín = V²/(127·(e_max+f_max))  [AASHTO]. Con e=8 %: 70 km/h -> 168 m, 80 km/h -> 229 m
    (coincide con los ejemplos resueltos de las diapositivas)."""
    return V * V / (127.0 * (e_max + F_MAX[int(V)]))


def peralte(R: float, V: float, e_max: float = 0.08, bombeo: float = 0.02) -> float:
    """Peralte en función del radio: e = e_max·Rmín/R, acotado a [bombeo, e_max].
    Distribución lineal con la curvatura (AASHTO). Las tablas 3.x del manual usan el Método 5:
    [VERIFICAR] antes de uso definitivo."""
    e = e_max * r_min(V, e_max) / R
    return float(min(max(e, bombeo), e_max))


def requiere_peralte(R: float) -> bool:
    return R <= 3500          # [DIAP] R > 3500 m: el bombeo basta


def requiere_espiral(R: float) -> bool:
    return R <= 1000          # [DIAP] R > 1000 m: no requiere espirales


def longitud_espiral_min(V: float, R: float, e: float, ancho_carril: float,
                         bombeo: float = 0.02, J: float = 0.6) -> dict:
    """Longitud mínima de espiral = máx de:
       - rampa de peralte: L = a·(e + bombeo)/Δ    (a = distancia eje de giro -> borde; aquí el ancho del carril)
       - comodidad (Shortt): L = V³/(46.656·J·R)
       - recorrido de 2 s a V
       [VERIFICAR coeficientes J y criterios contra el manual INVÍAS]."""
    vv = min(RAMPA_PERALTE, key=lambda k: abs(k - V))
    delta = RAMPA_PERALTE[vv] / 100.0
    l_rampa = ancho_carril * (e + bombeo) / delta
    l_com = V ** 3 / (46.656 * J * R)
    l_2s = V / 3.6 * 2.0
    return {"rampa": l_rampa, "comodidad": l_com, "dos_seg": l_2s,
            "Le_min": max(l_rampa, l_com, l_2s)}


# ---------------------------------------------------------------- Velocidades específicas [DIAP]
def vech(vtrh: float, vech_ant: float | None, long_tangente: float, deflexion_deg: float) -> float:
    """Velocidad específica de curva horizontal según la longitud del segmento recto anterior."""
    if vech_ant is None:                       # primera curva del tramo
        vech_ant = vtrh
    if long_tangente < 150:
        v = vech_ant                           # caso 1
    elif long_tangente <= 400:                 # casos 2 y 3
        v = vech_ant if deflexion_deg < 45 else vech_ant - 10
    elif long_tangente <= 600:
        v = vtrh + 10                          # caso 4
    else:
        v = vtrh + 20                          # caso 5
    return float(min(max(v, vtrh - 10), vtrh + 20))


def entretangencia_min(V: float, espiralizada: bool = True, plano: bool = False) -> float:
    """Curvas espiralizadas: recorrido de 5 s a VETH. Circulares en terreno plano: 15 s [DIAP]."""
    seg = 15.0 if (plano and not espiralizada) else 5.0
    return V / 3.6 * seg


def entretangencia_max(V: float) -> float:
    return 15.0 * V           # 15 veces VETH (km/h) [DIAP]


# ---------------------------------------------------------------- Curvas verticales [DIAP]
def k_convexa(V: float, h1: float = 1.08, h2: float = 0.60) -> float:
    """K = Dp²/(100·(√2h1+√2h2)²). Reproduce la tabla 4.4 (p.ej. 80 km/h: 25.7 con h2=0.60)."""
    dp = dp_tabla(int(V))
    return dp * dp / (100.0 * (math.sqrt(2 * h1) + math.sqrt(2 * h2)) ** 2)


def k_concava(V: float) -> float:
    """K = Dp²/(120+3.5·Dp)  (criterio de faros) [DIAP]."""
    dp = dp_tabla(int(V))
    return dp * dp / (120.0 + 3.5 * dp)


def l_operacion(V: float) -> float:
    return 0.6 * V             # criterio de operación Lmín = 0.6·V [DIAP]


# ---------------------------------------------------------------- Tangente vertical [DIAP]
# Tabla 4.3: longitud mínima de la tangente vertical (m) según velocidad específica VTV
LTV_MIN = {20: 40, 30: 60, 40: 80, 50: 140, 60: 170, 70: 195, 80: 225, 90: 250,
           100: 280, 110: 305, 120: 335, 130: 360}
PEND_MIN = 0.005          # 0,3 % - 0,5 % para drenaje [DIAP]; se usa el valor alto

# Longitud crítica de pendiente (pérdida de 25 km/h del vehículo pesado, 180 kg/HP, entrada 90 km/h).
# Leída de forma APROXIMADA de la Fig. 4.2 de las diapositivas (±15 %). [VERIFICAR] -> solo avisos.
L_CRITICA_APROX = {3: 650, 4: 440, 5: 340, 6: 260, 7: 220, 8: 190, 9: 165}


def v_snap(V: float) -> int:
    """Velocidad al múltiplo de 10 válido más cercano (20..130) para entrar a las tablas."""
    return int(min(max(round(V / 10.0) * 10, 20), 130))


def ltv_min(V: float) -> float:
    return float(LTV_MIN[v_snap(V)])


def l_critica(pend: float) -> float:
    """Longitud crítica aproximada (m) para una pendiente m/m (interpola la lectura de la figura)."""
    g = abs(pend) * 100.0
    if g <= 3:
        return float("inf")
    xs = sorted(L_CRITICA_APROX)
    if g >= xs[-1]:
        return float(L_CRITICA_APROX[xs[-1]])
    for a, b in zip(xs[:-1], xs[1:]):
        if a <= g <= b:
            t = (g - a) / (b - a)
            return L_CRITICA_APROX[a] * (1 - t) + L_CRITICA_APROX[b] * t
    return float("inf")


def _s_convexa(h1=1.08, h2=0.60):
    return 200.0 * (math.sqrt(h1) + math.sqrt(h2)) ** 2


def l_min_convexa(A_pct: float, Dp: float, h1: float = 1.08, h2: float = 0.60) -> float:
    """Criterio de seguridad, curva convexa. Dp < L: L = A·Dp²/S ; Dp > L: L = 2·Dp − S/A  [DIAP]
    con S = 200·(√h1+√h2)². Ambas ramas coinciden en A = S/Dp."""
    S = _s_convexa(h1, h2)
    if A_pct <= 0:
        return 0.0
    return A_pct * Dp * Dp / S if A_pct >= S / Dp else max(2 * Dp - S / A_pct, 0.0)


def l_min_concava(A_pct: float, Dp: float) -> float:
    """Criterio de seguridad (faros), curva cóncava. Dp < L: L = A·Dp²/(120+3.5·Dp) ;
    Dp > L: L = 2·Dp − (120+3.5·Dp)/A  [DIAP]."""
    S = 120.0 + 3.5 * Dp
    if A_pct <= 0:
        return 0.0
    return A_pct * Dp * Dp / S if A_pct >= S / Dp else max(2 * Dp - S / A_pct, 0.0)


def l_curva_vertical(A_pct: float, tipo: str, V: float, paso: float = 10.0) -> dict:
    """Longitud de curva vertical simétrica: máx(seguridad, operación 0,6·V), redondeada hacia arriba a `paso` m.
    tipo: 'convexa' | 'concava'. A en % (diferencia algebraica de pendientes)."""
    Vs = v_snap(V)
    Dp = dp_tabla(Vs)
    seg = l_min_convexa(A_pct, Dp) if tipo == "convexa" else l_min_concava(A_pct, Dp)
    op = l_operacion(Vs)
    Lmin = max(seg, op)
    L = math.ceil(Lmin / paso) * paso
    return dict(seguridad=seg, operacion=op, Lmin=Lmin, L=float(L), Dp=Dp, V=Vs,
                K=L / A_pct if A_pct > 0 else float("inf"),
                Kmin=(k_convexa(Vs) if tipo == "convexa" else k_concava(Vs)))


# ---------------------------------------------------------------- Selección de ruta (Bruce) [DIAP]
def longitud_virtual(segmentos, k: float = 44.0, PR: float = 0.09) -> float:
    """X0 = X + k·Σy + Σ(p − PR)·l   (método de Bruce, diapositivas de línea de ceros).
    segmentos: lista de (longitud_m, desnivel_m) con signo (+ ascenso).
    Σy: suma de ascensos; el término de pendiente penaliza descensos con p > PR."""
    X = sum(l for l, _ in segmentos)
    sy = sum(dz for _, dz in segmentos if dz > 0)
    pen = sum((abs(dz) / l - PR) * l for l, dz in segmentos if dz < 0 and l > 0 and abs(dz) / l > PR)
    return X + k * sy + pen


def verificar():
    """Comprobaciones contra los ejemplos resueltos de las diapositivas."""
    assert round(r_min(70)) == 168, r_min(70)
    assert round(r_min(80)) == 229, r_min(80)
    assert abs(k_convexa(80) - 25.7) < 0.1, k_convexa(80)
    assert abs(k_convexa(70) - 16.8) < 0.1, k_convexa(70)
    assert abs(dist_parada(80) - 130) < 3 and abs(dist_parada(70) - 105) < 2
    assert clasificar_terreno(85) == "escarpado"
    assert rango_vtr("primaria_1_calzada", "escarpado") == [60, 70, 80]
    assert abs(72 - 5.0 / 0.07) < 1       # d = equidistancia/pendiente (diap. 17: 7 % y 5 m -> 72 m)
    # ejemplo de diseño vertical de las diapositivas: A = 11,53 %, Dp = 135,089 m -> Lmín = 354,938 m
    assert abs(l_min_concava(11.53, 135.089) - 354.938) < 0.05, l_min_concava(11.53, 135.089)
    # continuidad de las dos ramas del criterio de seguridad
    for Dp in (65, 130):
        a = _s_convexa() / Dp
        assert abs(l_min_convexa(a * 0.999, Dp) - l_min_convexa(a * 1.001, Dp)) < 0.5
        b = (120 + 3.5 * Dp) / Dp
        assert abs(l_min_concava(b * 0.999, Dp) - l_min_concava(b * 1.001, Dp)) < 0.5
    assert ltv_min(80) == 225 and ltv_min(70) == 195          # tabla 4.3
    return True


if __name__ == "__main__":
    verificar()
    print("normas.py: verificaciones contra las diapositivas OK")
    print(f"  Rmin(70)={r_min(70):.1f}  Rmin(80)={r_min(80):.1f}  K convexa(80)={k_convexa(80):.1f}  K cóncava(80)={k_concava(80):.1f}")
