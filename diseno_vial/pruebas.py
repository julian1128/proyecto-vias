"""Batería de pruebas: python -m diseno_vial.pruebas [-v]

Cada prueba verifica una pieza con una respuesta conocida o una propiedad que debe cumplirse siempre.
No requiere archivos externos (usa terrenos sintéticos)."""
from __future__ import annotations
import math, os, sys, tempfile, traceback
import numpy as np

from . import normas, topografia as T, linea_ceros as LC, alineamiento as A, perfil as PF, secciones as SC, optimizar as OPT, replanteo as RP, trazado as TR, exportar as X, validar as V

_RES = []
# Ruta raster real (337 puntos) capturada de un terreno empinado: con ella `ajustar_pis` devolvía PI duplicados.
_RUTA_HORQUILLAS = "100.00,150.00 97.09,161.46 97.09,169.45 97.09,177.45 97.09,185.44 97.09,193.43 97.09,201.43 97.09,209.42 97.09,217.41 97.09,225.41 97.09,233.40 97.09,241.39 97.09,249.39 97.09,257.38 97.09,265.37 97.09,273.37 105.09,289.35 113.08,305.34 121.07,321.33 129.07,337.31 137.06,345.31 145.05,353.30 145.05,345.31 145.05,337.31 153.05,345.31 161.04,353.30 169.03,361.29 177.03,369.29 185.02,377.28 201.01,385.27 216.99,393.27 232.98,401.26 224.99,393.27 240.98,401.26 256.96,409.26 272.95,417.25 280.94,417.25 272.95,409.26 280.94,409.26 288.94,409.26 296.93,409.26 312.92,417.25 328.90,425.24 344.89,433.24 360.88,441.23 368.87,441.23 376.86,441.23 384.86,441.23 400.84,449.22 392.85,441.23 408.84,449.22 424.82,457.22 432.82,457.22 424.82,449.22 432.82,449.22 424.82,441.23 432.82,441.23 440.81,441.23 432.82,433.24 440.81,433.24 448.80,433.24 456.80,433.24 464.79,433.24 472.78,433.24 480.78,433.24 496.76,441.23 512.75,449.22 528.74,457.22 536.73,457.22 544.72,457.22 552.72,457.22 560.71,457.22 568.70,457.22 576.70,457.22 584.69,457.22 600.68,465.21 608.67,465.21 616.66,465.21 624.66,465.21 632.65,465.21 640.64,465.21 656.63,473.20 664.62,481.20 672.62,473.20 688.60,465.21 696.60,473.20 712.58,481.20 728.57,473.20 736.56,473.20 744.56,457.22 752.55,441.23 760.54,441.23 768.54,441.23 784.52,449.22 800.51,457.22 808.50,473.20 808.50,481.20 808.50,489.19 816.50,489.19 824.49,489.19 824.49,497.18 840.48,489.19 840.48,497.18 832.49,513.17 840.48,513.17 832.49,529.16 832.49,537.15 848.47,529.16 840.48,545.14 832.49,561.13 848.47,553.14 840.48,569.12 856.47,561.13 872.45,553.14 864.46,569.12 856.47,585.11 856.47,593.10 856.47,601.10 856.47,609.09 864.46,593.10 864.46,601.10 864.46,609.09 872.45,593.10 872.45,601.10 872.45,609.09 880.45,593.10 888.44,585.11 896.43,577.12 904.43,569.12 920.41,561.13 912.42,577.12 904.43,593.10 920.41,585.11 936.40,577.12 952.39,569.12 968.37,561.13 984.36,553.14 976.37,561.13 992.35,553.14 984.36,561.13 976.37,569.12 968.37,577.12 960.38,585.11 976.37,577.12 968.37,585.11 960.38,593.10 976.37,585.11 968.37,593.10 984.36,585.11 976.37,593.10 968.37,609.09 976.37,601.10 992.35,593.10 984.36,609.09 992.35,601.10 984.36,617.08 1000.35,609.09 1016.33,601.10 1032.32,593.10 1048.31,585.11 1064.29,577.12 1080.28,569.12 1096.27,561.13 1104.26,561.13 1112.25,561.13 1128.24,553.14 1144.23,545.14 1152.22,545.14 1160.21,545.14 1168.21,545.14 1176.20,545.14 1184.19,545.14 1192.19,545.14 1200.18,545.14 1208.17,545.14 1224.16,537.15 1232.15,537.15 1240.15,537.15 1248.14,537.15 1256.13,537.15 1264.13,537.15 1272.12,537.15 1280.11,537.15 1288.11,537.15 1296.10,537.15 1304.09,537.15 1312.09,537.15 1320.08,537.15 1328.07,537.15 1336.07,537.15 1344.06,537.15 1352.05,537.15 1360.05,537.15 1368.04,537.15 1376.03,537.15 1384.03,537.15 1392.02,537.15 1400.02,537.15 1408.01,537.15 1416.00,537.15 1424.00,537.15 1431.99,537.15 1447.98,545.14 1455.97,545.14 1463.96,545.14 1479.95,553.14 1487.94,561.13 1495.94,569.12 1503.93,577.12 1511.92,585.11 1519.92,593.10 1527.91,609.09 1535.90,625.08 1535.90,633.07 1535.90,641.06 1535.90,649.06 1543.90,665.04 1551.89,681.03 1551.89,689.02 1551.89,697.02 1551.89,705.01 1551.89,713.00 1551.89,721.00 1551.89,728.99 1551.89,736.98 1551.89,744.98 1551.89,752.97 1551.89,760.96 1551.89,768.96 1551.89,776.95 1551.89,784.94 1559.88,800.93 1559.88,808.92 1559.88,816.92 1559.88,824.91 1559.88,832.90 1559.88,840.90 1559.88,848.89 1559.88,856.88 1559.88,864.88 1559.88,872.87 1559.88,880.86 1559.88,888.86 1559.88,896.85 1559.88,904.84 1559.88,912.84 1559.88,920.83 1559.88,928.82 1559.88,936.82 1559.88,944.81 1567.88,960.80 1567.88,968.79 1567.88,976.79 1567.88,984.78 1567.88,992.77 1567.88,1000.77 1567.88,1008.76 1567.88,1016.75 1567.88,1024.75 1567.88,1032.74 1567.88,1040.73 1567.88,1048.73 1567.88,1056.72 1567.88,1064.71 1567.88,1072.71 1567.88,1080.70 1567.88,1088.69 1567.88,1096.69 1567.88,1104.68 1567.88,1112.67 1567.88,1120.67 1567.88,1128.66 1567.88,1136.65 1567.88,1144.65 1567.88,1152.64 1567.88,1160.63 1567.88,1168.63 1567.88,1176.62 1575.87,1192.61 1583.86,1208.59 1583.86,1216.59 1583.86,1224.58 1583.86,1232.57 1583.86,1240.57 1583.86,1248.56 1583.86,1256.55 1583.86,1264.55 1583.86,1272.54 1583.86,1280.53 1583.86,1288.53 1583.86,1296.52 1583.86,1304.51 1583.86,1312.51 1583.86,1320.50 1583.86,1328.49 1583.86,1336.49 1583.86,1344.48 1591.86,1352.47 1607.84,1360.47 1615.84,1368.46 1623.83,1384.45 1615.84,1400.43 1615.84,1408.43 1615.84,1416.42 1615.84,1424.41 1623.83,1424.41 1631.82,1424.41 1623.83,1440.40 1615.84,1448.39 1607.84,1464.38 1607.84,1472.37 1607.84,1480.37 1623.83,1472.37 1639.82,1464.38 1655.80,1456.39 1647.81,1472.37 1663.80,1464.38 1655.80,1480.37 1671.79,1472.37 1663.80,1488.36 1679.78,1480.37 1671.79,1496.35 1687.78,1488.36 1695.77,1480.37 1687.78,1496.35 1703.76,1488.36 1695.77,1504.35 1687.78,1520.33 1695.77,1512.34 1711.76,1504.35 1727.74,1496.35 1743.73,1488.36 1735.74,1496.35 1727.74,1504.35 1743.73,1496.35 1735.74,1504.35 1750.00,1500.00"




def prueba(f):
    _RES.append(f)
    return f


# ------------------------------------------------------------------ terrenos sintéticos
def terreno_plano_inclinado(pend=0.03, L=2000.0, n=6000, seed=0):
    rng = np.random.default_rng(seed)
    x, y = rng.uniform(0, L, n), rng.uniform(0, L, n)
    return np.column_stack([x, y, 1000.0 + pend * x])


def terreno_ondulado(L=2000.0, n=20000, seed=1, ruido=0.03):
    rng = np.random.default_rng(seed)
    x, y = rng.uniform(0, L, n), rng.uniform(0, L, n)
    z = 1000 + 0.02 * x + 0.015 * y + 25 * np.sin(x / 450) * np.cos(y / 600) + 12 * np.sin(y / 300 + 1)
    return np.column_stack([x, y, z + rng.normal(0, ruido, n)])


# ------------------------------------------------------------------ normas
@prueba
def normas_contra_diapositivas():
    assert normas.verificar()


@prueba
def normas_tablas_consistentes():
    """Toda celda permitida de la tabla 5.2 (calzada) tiene berma en la 5.4. La 5.4 de las diapositivas trae
    una celda extra conocida (terciaria-escarpado a 40 km/h) que la 5.2 no permite: se documenta aquí."""
    discrepancias_conocidas = {("terciaria", "escarpado"): {40}}
    for cat in normas.CATEGORIAS:
        for ter in normas.TERRENOS:
            a, b = set(normas.ANCHO_CALZADA[cat][ter]), set(normas.BERMA[cat][ter])
            assert a <= b, ("celda de calzada sin berma", cat, ter, a - b)
            assert (b - a) == discrepancias_conocidas.get((cat, ter), set()), ("discrepancia nueva", cat, ter, b - a)
    for cat, tabla in normas.PEND_MEDIA_MAX.items():
        for v in tabla:
            assert any(v in normas.ANCHO_CALZADA[cat][t] for t in normas.TERRENOS), (cat, v)


@prueba
def normas_monotonia():
    """Rmín crece con la velocidad; la longitud mínima de tangente vertical también."""
    vs = sorted(normas.F_MAX)
    rm = [normas.r_min(v) for v in vs]
    assert all(b > a for a, b in zip(rm[:-1], rm[1:]))
    lt = [normas.LTV_MIN[v] for v in sorted(normas.LTV_MIN)]
    assert all(b > a for a, b in zip(lt[:-1], lt[1:]))
    # el peralte nunca excede e_max ni baja del bombeo
    for R in (50, 168, 500, 1000, 3000, 8000):
        e = normas.peralte(R, 70)
        assert 0.02 <= e <= 0.08


# ------------------------------------------------------------------ geometría vertical
@prueba
def parabola_geometria():
    s_pvi, z_pvi = np.array([0.0, 500.0, 1000.0]), np.array([80.0, 100.0, 85.0])     # +4 % y -3 %
    c = dict(i=1, s_pvi=500.0, L=200.0, A=-0.03 - 0.04)
    z = lambda s: float(PF.cota_diseno([s], s_pvi, z_pvi, [c])[0])
    assert abs(z(400) - 96.0) < 1e-9 and abs(z(600) - 97.0) < 1e-9          # PCV y PTV sobre las tangentes
    E = abs(c["A"]) * c["L"] / 8.0
    assert abs(z(500) - (100.0 - E)) < 1e-9 and abs(E - 1.75) < 1e-9        # externa E = A·L/8
    # continuidad de pendiente en PCV y PTV (derivada numérica)
    d = 1e-4
    for s0 in (400.0, 600.0):
        izq, der = (z(s0) - z(s0 - d)) / d, (z(s0 + d) - z(s0)) / d
        assert abs(izq - der) < 1e-3, (s0, izq, der)
    # coincide con la parábola analítica z = z_pcv + g1·x + A·x²/(2L)
    for x in (10, 60, 100, 150, 190):
        analitica = 96.0 + 0.04 * x + c["A"] * x * x / (2 * c["L"])
        assert abs(z(400 + x) - analitica) < 1e-9


@prueba
def curva_vertical_longitudes():
    for A in (1, 2, 4, 6, 8, 10):
        for tipo in ("convexa", "concava"):
            r = normas.l_curva_vertical(A, tipo, 80)
            assert r["L"] >= r["Lmin"] >= r["operacion"] and r["L"] % 10 == 0
            assert r["L"] - r["Lmin"] < 10.0 + 1e-9
    # más desnivel -> curva más larga
    L = [normas.l_curva_vertical(a, "convexa", 80)["L"] for a in (2, 4, 6, 8, 10)]
    assert all(b >= a for a, b in zip(L[:-1], L[1:]))


# ------------------------------------------------------------------ pipeline en terrenos de respuesta conocida
def _via_recta(pts, pend=0.07):
    terr = T.Terreno(pts, celda_tin=20.0)
    P, _ = LC.linea_de_ceros(terr, (100, 1000), (1900, 1000), pend, dh=5.0)
    al = A.disenar(P, 70, tol_pi=25)
    return terr, P, al


@prueba
def perfil_terreno_inclinado_recto():
    """Terreno plano inclinado al 3 %: rasante = terreno, sin curvas, corte/relleno ~ 0."""
    terr, P, al = _via_recta(terreno_plano_inclinado(0.03))
    assert len(al.curvas) == 0, "una vía recta no debe tener curvas horizontales"
    pf = PF.disenar_perfil(terr, al, "primaria_1_calzada", 70)
    assert pf.stats["n_pvi"] == 2 and pf.stats["n_curvas"] == 0, pf.stats
    assert np.max(np.abs(pf.zt - pf.zd)) < 0.15, np.max(np.abs(pf.zt - pf.zd))
    g = (pf.z_pvi[1] - pf.z_pvi[0]) / (pf.s_pvi[1] - pf.s_pvi[0])
    assert abs(g - 0.03) < 0.002, g
    assert not pf.violaciones


@prueba
def perfil_terreno_plano_horizontal():
    """Terreno horizontal: debe cumplir la pendiente mínima de drenaje (no rasante a 0 %)."""
    terr, P, al = _via_recta(terreno_plano_inclinado(0.0))
    pf = PF.disenar_perfil(terr, al, "primaria_1_calzada", 70)
    g = np.diff(pf.z_pvi) / np.diff(pf.s_pvi)
    assert np.all(np.abs(g) >= normas.PEND_MIN - 1e-6), g
    assert not pf.violaciones, pf.violaciones


@prueba
def perfil_cumple_restricciones_en_terreno_ondulado():
    terr = T.Terreno(terreno_ondulado(), celda_tin=20.0)
    P, _ = LC.linea_de_ceros(terr, (100, 150), (1850, 1700), 0.07, dh=5.0)
    al = A.disenar(P, 70, tol_pi=25)
    pf = PF.disenar_perfil(terr, al, "primaria_1_calzada", 70)
    assert not pf.violaciones, pf.violaciones
    g = np.diff(pf.z_pvi) / np.diff(pf.s_pvi)
    assert np.all(np.abs(g) <= pf.gmax + 1e-6) and np.all(np.abs(g) >= pf.gmin - 1e-6)
    # pendiente numérica de la rasante completa (incluye curvas) nunca supera la máxima
    gn = np.abs(np.diff(pf.zd) / np.diff(pf.s))
    assert gn.max() <= pf.gmax + 0.002, gn.max()
    # cada curva vertical cumple su longitud mínima y las tangentes entre curvas la tabla 4.3
    cv = [c for c in pf.curvas if c["L"] > 0]
    for c in cv:
        assert c["L"] >= c["Lmin"] - 1e-9
    for a, b in zip(cv[:-1], cv[1:]):
        gap = (b["s_pvi"] - b["L"] / 2) - (a["s_pvi"] + a["L"] / 2)
        assert gap >= normas.ltv_min(max(a["V"], b["V"])) - 1e-6 or (a["tipo"] != b["tipo"] and gap >= -1e-6), gap
    # el movimiento de tierras no puede ser un orden de magnitud peor que la cota inferior teórica
    # (rasante libre punto a punto con la misma pendiente máxima)
    tot = pf.stats["area_corte_m2"] + pf.stats["area_relleno_m2"]
    ratio = tot / pf.stats["area_total_cota_inferior_m2"]
    assert ratio < 6.0, f"movimiento de tierras {ratio:.1f}x la cota inferior"
    assert pf.stats["corte_max_m"] < 10 and pf.stats["relleno_max_m"] < 10, pf.stats


@prueba
def perfil_determinista():
    terr = T.Terreno(terreno_ondulado(), celda_tin=20.0)
    P, _ = LC.linea_de_ceros(terr, (100, 150), (1850, 1700), 0.07, dh=5.0)
    al = A.disenar(P, 70, tol_pi=25)
    a = PF.disenar_perfil(terr, al, "primaria_1_calzada", 70)
    b = PF.disenar_perfil(terr, al, "primaria_1_calzada", 70)
    assert np.allclose(a.z_pvi, b.z_pvi) and np.allclose(a.zd, b.zd)


# ------------------------------------------------------------------ alineamiento
@prueba
def alineamiento_continuidad_geometrica():
    terr = T.Terreno(terreno_ondulado(), celda_tin=20.0)
    P, _ = LC.linea_de_ceros(terr, (100, 150), (1850, 1700), 0.07, dh=5.0)
    al = A.disenar(P, 70, tol_pi=25)
    assert al.curvas, "se esperaban curvas en este terreno"
    S = al.muestrear(1.0)
    seg = np.diff(S[:, :2], axis=0)
    ln = np.hypot(*seg.T)
    assert ln.min() > 1e-9, "muestras duplicadas en el eje"
    assert ln.max() < 3.0, f"salto de {ln.max():.2f} m en la muestra: elementos discontinuos"
    ang = np.degrees(np.unwrap(np.arctan2(seg[:, 1], seg[:, 0])))      # unwrap en radianes, luego a grados
    salto = np.abs(np.diff(ang))
    # sin espirales/curvas el rumbo cambia suavemente: ningún cambio brusco entre muestras de ~2 m
    assert salto.max() < 4.0, f"cambio brusco de rumbo {salto.max():.2f}°"
    suma = sum(e["longitud"] for e in al.elementos)
    assert abs(suma - al.longitud) < 1e-6
    assert abs(S[-1, 2] - al.longitud) / al.longitud < 5e-3
    # el eje empieza y termina en BP y EP
    assert np.hypot(*(S[0, :2] - al.pis[0])) < 1e-6 and np.hypot(*(S[-1, :2] - al.pis[-1])) < 1e-6


@prueba
def alineamiento_normas_horizontales():
    terr = T.Terreno(terreno_ondulado(), celda_tin=20.0)
    P, _ = LC.linea_de_ceros(terr, (100, 150), (1850, 1700), 0.07, dh=5.0)
    al = A.disenar(P, 70, tol_pi=25)
    for c in al.curvas:
        assert c["R"] >= normas.r_min(c["V"]) - 1e-6, "radio menor que el mínimo"
        assert 0.02 <= c["e"] <= 0.08
        if normas.requiere_espiral(c["R"]):
            assert c["Le"] > 0 and c["Lc"] >= A.LC_MIN - 1e-6 or any("PI-%d" % c["idx"] in v for v in al.violaciones)
        # R, deflexión y longitudes son coherentes: Lc = R·(Δ − 2θs)
        assert abs(c["Lc"] - c["R"] * (c["delta"] - 2 * c["theta"])) < 1e-6
        # |CS - C| = |SC - C| = R
        assert abs(np.hypot(*(c["SC"] - c["C"])) - c["R"]) < 1e-3
        assert abs(np.hypot(*(c["CS"] - c["C"])) - c["R"]) < 1e-3


@prueba
def alineamiento_terreno_plano_sin_curvas():
    terr, P, al = _via_recta(terreno_plano_inclinado(0.0))
    assert len(al.curvas) == 0 and len(al.elementos) == 1 and al.elementos[0]["tipo"] == "linea"


# ------------------------------------------------------------------ robustez: errores esperados
@prueba
def errores_esperados():
    ok = False
    try:
        normas.ancho_calzada("primaria_1_calzada", "ondulado", 60)
    except ValueError as e:
        ok = "no permitida" in str(e)
    assert ok, "VTR inválida debe lanzar ValueError claro"
    terr = T.Terreno(terreno_plano_inclinado(0.06), celda_tin=20.0)
    ok = False
    try:
        LC.linea_de_ceros(terr, (100, 1000), (1900, 1000), 0.001, dh=5.0)   # compás de 5 km: no cabe en el terreno
    except RuntimeError as e:
        ok = "pendiente" in str(e)
    assert ok, "destino inalcanzable debe lanzar RuntimeError claro"


# ------------------------------------------------------------------ integración y archivos de salida
@prueba
def integracion_completa_y_archivos():
    from .main import disenar_via
    pts = terreno_ondulado()
    with tempfile.TemporaryDirectory() as d:
        r = disenar_via(pts, (100, 150), (1850, 1700), "primaria_1_calzada", vtr=70, terreno="ondulado",
                        equidistancia=5.0, salida=d, verbose=False, optimizar_rasante=False, replantear=True)
        assert r["trz"] is None or os.path.getsize(os.path.join(d, "optimizacion_trazado.csv")) > 0
        for f in ("via.xml", "via.dxf", "perfil.dxf", "secciones.dxf", "planta.png", "perfil.png", "secciones.png",
                  "curvas.csv", "pvi.csv", "secciones.csv", "resumen_movimiento_tierras.csv",
                  "replanteo_rutas.csv"):
            assert os.path.getsize(os.path.join(d, f)) > 0, f
        errores, info = V.validar_landxml(os.path.join(d, "via.xml"))
        assert not errores, errores
        assert info["perfil"]["PVI"] >= 2 and info["puntos"] > 100
        import ezdxf
        for f in ("via.dxf", "perfil.dxf", "secciones.dxf"):
            doc = ezdxf.readfile(os.path.join(d, f))
            assert len(doc.audit().errors) == 0, f
        # el perfil del LandXML termina donde termina el alineamiento
        assert abs(info["longitud"] - r["al"].longitud) < 0.01


@prueba
def validador_detecta_errores():
    """El validador debe rechazar un LandXML con el eje roto (prueba de que no da falsos OK)."""
    from .main import disenar_via
    with tempfile.TemporaryDirectory() as d:
        disenar_via(terreno_ondulado(), (100, 150), (1850, 1700), "primaria_1_calzada", vtr=70, terreno="ondulado",
                    salida=d, verbose=False, optimizar_rasante=False, replantear=False)
        ruta = os.path.join(d, "via.xml")
        txt = open(ruta, encoding="utf-8").read()
        i = txt.index("<Curve ")
        j = txt.index("<Start>", i)
        k = txt.index("</Start>", j)
        roto = txt[:j] + "<Start>1.0 1.0" + txt[k:]            # inicio de curva movido: discontinuidad
        ruta2 = os.path.join(d, "roto.xml")
        open(ruta2, "w", encoding="utf-8").write(roto)
        errores, _ = V.validar_landxml(ruta2)
        assert errores, "el validador no detectó un eje discontinuo"


# ------------------------------------------------------------------ regresiones de defectos corregidos
def _zigzag():
    xy = np.array([(0, 0), (200, 0), (200, 200), (400, 200), (400, 400)], float)
    return np.column_stack([xy, np.zeros(len(xy))])


@prueba
def alineamiento_detecta_solapes():
    """Giros de 90° cada 200 m a 100 km/h (R mín 394 m): las curvas no caben y deben marcarse inválidas."""
    al = A._disenar_una(_zigzag(), 100, tol_pi=5.0)
    assert not al.valido and al.solapes, "no se detectó el solape de curvas"
    assert al.longitud > 0


@prueba
def perfil_rechaza_eje_invalido():
    al = A._disenar_una(_zigzag(), 100, tol_pi=5.0)
    terr = T.Terreno(terreno_plano_inclinado(0.0, L=600.0, n=3000), celda_tin=10.0)
    try:
        PF.disenar_perfil(terr, al, "primaria_1_calzada", 100)
    except ValueError as e:
        assert "inválido" in str(e)
        return
    raise AssertionError("el perfil no debe diseñarse sobre un eje con curvas solapadas")


@prueba
def alineamiento_no_elige_recta_trampa():
    """Una recta sin curvas 'cumple' las normas por vacío: no debe elegirse si se aleja del terreno."""
    terr = T.Terreno(terreno_ondulado(), celda_tin=20.0)
    P, _ = LC.linea_de_ceros(terr, (100, 150), (1850, 1700), 0.07, dh=5.0)
    al = A.disenar(P, 70, tol_pi=25)
    extension = max(np.ptp(P[:, 0]), np.ptp(P[:, 1]))
    assert al.valido
    assert len(al.curvas) >= 1, "se eligió una recta que ignora la línea de ceros"
    assert al.desv_max <= max(50.0, 0.03 * extension) + 1e-6, al.desv_max


# ------------------------------------------------------------------ sección transversal, peralte, sobreancho, volúmenes
_CACHE = {}


def _caso_ondulado():
    """Terreno, eje y rasante de referencia (se calcula una sola vez)."""
    if "c" not in _CACHE:
        terr = T.Terreno(terreno_ondulado(), celda_tin=20.0)
        P, _ = LC.linea_de_ceros(terr, (100, 150), (1850, 1700), 0.07, dh=5.0)
        al = A.disenar(P, 70, tol_pi=25)
        pf = PF.disenar_perfil(terr, al, "primaria_1_calzada", 70)
        _CACHE["c"] = (terr, al, pf)
    return _CACHE["c"]


@prueba
def secciones_sobreancho_ejemplo_diapositivas():
    """Ejemplo resuelto: R=174 m, V=60, calzada 7.20 -> U=3.048, FA=0.121, Z=0.059, Ac=8.076, S=0.876."""
    r = SC.sobreancho(174, 60, 7.20, 2, "ejemplo_diapositiva")
    assert abs(r["U"] - 3.048) < 0.001 and abs(r["FA"] - 0.121) < 0.001 and abs(r["Z"] - 0.059) < 0.001
    assert abs(r["Ac"] - 8.076) < 0.002 and abs(r["S"] - 0.876) < 0.002, r
    assert abs(SC.radio_carril_interior(174, 3.60, 0.876) - 171.543) < 0.001
    # a mayor radio, menor sobreancho; nunca negativo
    S = [SC.sobreancho(R, 70, 7.30)["S"] for R in (150, 250, 400, 800, 2000)]
    assert all(b <= a for a, b in zip(S[:-1], S[1:])) and min(S) >= 0.0


@prueba
def secciones_flecha_M_ejemplo_diapositivas():
    """Ejemplo: R'c = 171.543 m, Dp = 75.042 m -> M = 4.088 m."""
    assert abs(SC.flecha_M(171.543, 75.042) - 4.088) < 0.005
    Dp = 0.695 * 60 + 60 ** 2 / (254 * (3.4 / 9.81 + 0.0785))
    assert abs(Dp - 75.042) < 0.01


@prueba
def secciones_peralte_transicion():
    m, e = 0.02, 0.08
    assert SC.pendientes_peralte(0.0, e, m) == (-m, -m)
    so, si = SC.pendientes_peralte(1.0, e, m)
    assert abs(so - e) < 1e-12 and abs(si + e) < 1e-12
    # el borde externo pasa por cero en f = m/(e+m) y el interno solo empieza a bajar cuando el externo supera +m
    f0 = m / (e + m)
    assert abs(SC.pendientes_peralte(f0, e, m)[0]) < 1e-12 and SC.pendientes_peralte(f0, e, m)[1] == -m
    prev = SC.pendientes_peralte(0.0, e, m)
    for f in np.linspace(0.01, 1.0, 100):
        so, si = SC.pendientes_peralte(f, e, m)
        assert so >= prev[0] - 1e-12 and si <= prev[1] + 1e-12          # externo sube, interno nunca sube
        prev = (so, si)


@prueba
def secciones_areas_formula_cerrada():
    """Áreas numéricas = fórmulas cerradas (terreno plano o inclinado, con/sin cuneta, con/sin espesor de pavimento)."""
    dd = 0.05
    d = np.arange(0, 80 + dd / 2, dd)
    base = dict(bombeo=0.0, berma_pend=0.0, e_pav=0.0, cuneta_ancho=0.0, berma=1.8)
    b = 3.65 + 1.8
    lado = lambda zt, zr, st: SC._lado(d, zt, zr, 0.0, 3.65, st)
    # A) terreno inclinado s=0.20 con el eje a nivel: parte bajo la plataforma s·b²/2 + triángulo del talud
    s = 0.20
    st = SC.SeccionTipo(tc_h=1.0, tf_h=1.5, **base)
    cut, fil = lado(+s * d, 0.0, st), lado(-s * d, 0.0, st)
    Ac = s * b * b / 2 + 0.5 * s * b * (b * s * st.tc_h / (1 - s * st.tc_h))
    Af = s * b * b / 2 + 0.5 * s * b * (b * s * st.tf_h / (1 - s * st.tf_h))
    assert cut["modo"] == "corte" and fil["modo"] == "relleno"
    assert abs(cut["Ac"] - Ac) < 0.01 and abs(fil["Af"] - Af) < 0.01, (cut["Ac"], Ac, fil["Af"], Af)
    assert abs(cut["d_chaflan"] - b / (1 - s * st.tc_h)) < 0.01 and abs(fil["d_chaflan"] - b / (1 - s * st.tf_h)) < 0.01
    h = 2.0
    # B) relleno sobre terreno horizontal
    r = lado(0 * d, h, SC.SeccionTipo(tf_h=1.5, **base))
    assert abs(r["Af"] - (h * b + 0.5 * h * 1.5 * h)) < 0.01 and r["Ac"] == 0.0
    # C) corte sin cuneta
    r = lado(0 * d, -h, SC.SeccionTipo(tc_h=1.0, **base))
    assert abs(r["Ac"] - (h * b + 0.5 * h * 1.0 * h)) < 0.01
    # D) corte con cuneta 0.9 m de ancho y talud interior 3:1 (profundidad 0.3 m)
    r = lado(0 * d, -h, SC.SeccionTipo(tc_h=1.0, **{**base, "cuneta_ancho": 0.9}))
    assert abs(r["Ac"] - (h * b + 0.9 * (h + 0.15) + 0.5 * (h + 0.3) * (h + 0.3))) < 0.01
    # E) relleno con espesor de pavimento 0.5 m: el movimiento de tierras se mide a subrasante
    r = lado(0 * d, h, SC.SeccionTipo(tf_h=1.5, **{**base, "e_pav": 0.5}))
    assert abs(r["Af"] - (1.5 * b + 1.5 * (0.5 * 1.5) + 0.5 * 1.5 * (1.5 * 1.5))) < 0.01


@prueba
def secciones_diagrama_masas_caso_manual():
    s_ = np.arange(0, 161, 20.0)
    Vc = np.array([100, 100, 100, 100, 0, 0, 0, 0.0]); Vf = np.array([0, 0, 0, 0, 100, 100, 100, 100.0])
    m = SC.diagrama_masas(s_, Vc, Vf, 1.0)
    assert np.allclose(m["ordenada"], [0, 100, 200, 300, 400, 300, 200, 100, 0])
    assert m["acarreo_long_m3"] == 400 and abs(m["acarreo_m3m"] - 32000) < 1e-6 and abs(m["dma_m"] - 80) < 1e-9
    assert m["excedente_m3"] == 0 and m["faltante_m3"] == 0 and len(m["lobulos"]) == 1
    m = SC.diagrama_masas(s_, Vc, Vf, 1.25)           # el terraplén compacto equivale a 500 m³ de corte
    assert abs(m["ordenada_final"] + 100) < 1e-9 and abs(m["faltante_m3"] - 100) < 1e-9
    # compensación transversal: corte y relleno en el mismo tramo
    m = SC.diagrama_masas(np.array([0, 10.0]), np.array([50.0]), np.array([20.0]), 1.0)
    assert m["transversal_m3"] == 20 and m["excedente_m3"] == 30
    # volúmenes por áreas extremas de un prisma
    assert abs(SC.volumenes_areas_extremas([0, 50, 100], [10, 10, 10]).sum() - 1000) < 1e-9


@prueba
def secciones_chaflan_sobre_terreno_y_talud():
    """Verificación independiente: cada chaflán cae sobre el TIN y sobre la recta del talud."""
    terr, al, pf = _caso_ondulado()
    st = SC.SeccionTipo(tc_h=1.0, tf_h=1.5)
    sec = SC.calcular(terr, al, pf, st)
    comprobadas = 0
    for k in range(0, len(sec.s), 3):
        for lado_, signo in ((sec.izq, +1.0), (sec.der, -1.0)):
            l = lado_[k]
            if l["sin_chaflan"]:
                continue
            px = sec.x[k] + signo * sec.nrm[k, 0] * l["d_chaflan"]
            py = sec.y[k] + signo * sec.nrm[k, 1] * l["d_chaflan"]
            assert abs(terr.z([px], [py])[0] - l["z_chaflan"]) < 0.06, (k, terr.z([px], [py])[0], l["z_chaflan"])
            if l["modo"] == "corte" and l["dB"] is not None and l["d_chaflan"] > l["dB"] + 1e-6:
                pend = (l["z_chaflan"] - l["zB"]) / (l["d_chaflan"] - l["dB"])
                assert abs(pend - 1 / st.tc_h) < 1e-6, pend
            if l["modo"] == "relleno" and l["d_chaflan"] > l["dE"] + 1e-6:
                pend = (l["zE"] - l["z_chaflan"]) / (l["d_chaflan"] - l["dE"])
                assert abs(pend - 1 / st.tf_h) < 1e-6, pend
            comprobadas += 1
    assert comprobadas > 100


@prueba
def secciones_peralte_y_sobreancho_coherentes_con_alineamiento():
    terr, al, pf = _caso_ondulado()
    st = SC.SeccionTipo()
    sec = SC.calcular(terr, al, pf, st)
    assert al.curvas
    en_curva = np.zeros(len(sec.s), bool)
    for c in al.curvas:
        a0, a3 = c["s_ini"], c["s_fin"]
        en_curva |= (sec.s >= a0) & (sec.s <= a3)
        # en el tramo circular el peralte es el pleno y el sobreancho el calculado
        k = int(np.argmin(np.abs(sec.s - (a0 + c["Le"] + c["Lc"] / 2))))
        if c["Lc"] > 10:
            assert abs(sec.peralte_ext[k] - c["e"]) < 1e-9, (sec.peralte_ext[k], c["e"])
            S = SC.sobreancho(c["R"], c["V"], st.calzada, st.n_carriles, st.vehiculo)["S"]
            assert abs(sec.sobreancho[k] - S) < 1e-9
            assert sec.lado_sobreancho[k] == (+1 if c["sgn"] > 0 else -1), "el sobreancho va por el lado interno"
            # dirección física del peralte: el borde EXTERNO sube y el INTERNO baja respecto al eje; el carril interno es más ancho
            zr_k = sec.zr[k]
            (cl, zl), (cd, zd_) = sec.izq[k]["puntos"][1], sec.der[k]["puntos"][1]
            if c["sgn"] > 0:                      # curva a la izquierda: interno = izquierdo, externo = derecho
                assert zd_ > zr_k > zl and cl > cd, "peralte volcado hacia el lado equivocado"
                assert abs((zd_ - zr_k) - c["e"] * cd) < 1e-9 and abs((zr_k - zl) - c["e"] * cl) < 1e-9
            else:                                 # curva a la derecha: interno = derecho, externo = izquierdo
                assert zl > zr_k > zd_ and cd > cl, "peralte volcado hacia el lado equivocado"
                assert abs((zl - zr_k) - c["e"] * cl) < 1e-9 and abs((zr_k - zd_) - c["e"] * cd) < 1e-9
        # la rampa de peralte no supera la pendiente relativa máxima (la espiral se dimensionó con ese criterio)
        delta = normas.RAMPA_PERALTE[min(normas.RAMPA_PERALTE, key=lambda v: abs(v - c["V"]))] / 100.0
        m_ = (sec.s >= a0) & (sec.s <= a0 + c["Le"])
        borde = sec.peralte_ext[m_] * st.carril
        if m_.sum() > 2:
            grad = np.abs(np.diff(borde) / np.diff(sec.s[m_]))
            assert grad.max() <= delta + 1e-6, (grad.max(), delta)
    assert np.allclose(sec.peralte_ext[~en_curva], -st.bombeo) and np.allclose(sec.sobreancho[~en_curva], 0.0)


@prueba
def secciones_volumenes_y_masas_consistentes():
    terr, al, pf = _caso_ondulado()
    sec = SC.calcular(terr, al, pf, SC.SeccionTipo(factor_terraplen=1.25))
    vc = np.sum((sec.Ac[:-1] + sec.Ac[1:]) / 2 * np.diff(sec.s))
    vf = np.sum((sec.Af[:-1] + sec.Af[1:]) / 2 * np.diff(sec.s))
    assert abs(vc - sec.resumen["corte_m3"]) < 1e-6 and abs(vf - sec.resumen["relleno_compacto_m3"]) < 1e-6
    assert abs(sec.masas["ordenada"][-1] - (vc - 1.25 * vf)) < 1e-6
    assert np.all(sec.Ac >= 0) and np.all(sec.Af >= 0) and sec.masas["acarreo_long_m3"] >= 0
    assert abs(sec.resumen["excedente_m3"] - max(vc - 1.25 * vf, 0)) < 1e-6
    t = SC.tabla(sec)
    assert len(t) == len(sec.s) and t[0]["abscisa"] == 0.0 and abs(t[-1]["abscisa"] - al.longitud) < 1e-9
    # convergencia de la integración numérica al afinar la malla transversal
    fino = SC.calcular(terr, al, pf, SC.SeccionTipo(factor_terraplen=1.25), dd=0.1)
    assert abs(fino.resumen["corte_m3"] - vc) / vc < 0.01 and abs(fino.resumen["relleno_compacto_m3"] - vf) / vf < 0.02


@prueba
def secciones_sensibilidad_a_condiciones_de_la_via():
    """El movimiento de tierras debe responder en el sentido físico a ancho, talud de corte y de terraplén."""
    terr, al, pf = _caso_ondulado()
    base = SC.calcular(terr, al, pf, SC.SeccionTipo(tc_h=1.0, tf_h=1.5)).resumen
    roca = SC.calcular(terr, al, pf, SC.SeccionTipo(tc_h=0.25, tf_h=1.5)).resumen
    blando = SC.calcular(terr, al, pf, SC.SeccionTipo(tc_h=1.5, tf_h=1.5)).resumen
    tendido = SC.calcular(terr, al, pf, SC.SeccionTipo(tc_h=1.0, tf_h=2.0)).resumen
    ancho = SC.calcular(terr, al, pf, SC.SeccionTipo(calzada=10.9, tc_h=1.0, tf_h=1.5)).resumen
    assert roca["corte_m3"] < base["corte_m3"] < blando["corte_m3"], "talud de corte más vertical -> menos corte"
    assert roca["ancho_afectacion_medio_m"] < base["ancho_afectacion_medio_m"] < blando["ancho_afectacion_medio_m"]
    assert tendido["relleno_compacto_m3"] > base["relleno_compacto_m3"], "talud de terraplén más tendido -> más relleno"
    assert ancho["corte_m3"] > base["corte_m3"] and ancho["relleno_compacto_m3"] > base["relleno_compacto_m3"]
    assert ancho["costo"]["total"] > base["costo"]["total"]


@prueba
def secciones_rechazos_esperados():
    zz = _zigzag()
    al = A._disenar_una(zz, 100, tol_pi=5.0)
    terr = T.Terreno(terreno_plano_inclinado(0.0, L=600.0, n=3000), celda_tin=10.0)
    try:
        SC.calcular(terr, al, None, SC.SeccionTipo())
    except ValueError:
        pass
    else:
        raise AssertionError("no debe calcular secciones sobre un eje inválido")
    from .main import disenar_via
    try:
        disenar_via(terreno_plano_inclinado(0.03), (100, 1000), (1900, 1000), "primaria_2_calzadas", vtr=100,
                    terreno="plano", salida=tempfile.mkdtemp(), verbose=False)
    except NotImplementedError as e:
        assert "Calzadas separadas" in str(e)
    else:
        raise AssertionError("dos calzadas aún no está soportado y debe decirlo")


# ------------------------------------------------------------------ rasante optimizada por volúmenes
def _fuerza_bruta(costo, W):
    import itertools
    n, M = costo.shape
    mejor, ruta_m = np.inf, None
    for ruta in itertools.product(range(M), repeat=n):
        if any(abs(ruta[k + 1] - ruta[k]) > W[k] for k in range(n - 1)):
            continue
        c = sum(costo[k, ruta[k]] for k in range(n))
        if c < mejor:
            mejor, ruta_m = c, ruta
    return mejor, ruta_m


@prueba
def optimizar_dp_malla_igual_fuerza_bruta():
    rng = np.random.default_rng(7)
    for _ in range(25):
        n, M = 5, 6
        costo = rng.uniform(0, 10, (n, M))
        costo[rng.random((n, M)) < 0.2] = np.inf                     # niveles no permitidos
        W = rng.integers(1, 3, n - 1)
        bruto, _ = _fuerza_bruta(costo, W)
        if not np.isfinite(bruto):
            continue
        ruta, total = OPT.dp_malla(costo, W)
        assert abs(total - bruto) < 1e-9, (total, bruto)
        assert all(abs(ruta[k + 1] - ruta[k]) <= W[k] for k in range(n - 1)), "ruta con salto de pendiente no permitido"
        assert abs(sum(costo[k, ruta[k]] for k in range(n)) - total) < 1e-9
        # la versión genérica (cotas explícitas) debe dar el mismo costo sobre la misma malla
        z = np.tile(np.arange(M, dtype=float), (n, 1))
        ruta2, total2 = OPT.dp_rasante(costo, z, W.astype(float))
        assert abs(total2 - total) < 1e-9


@prueba
def optimizar_dp_no_peor_que_solucion_factible():
    """Regresión del defecto original: el óptimo (λ=0) no puede costar más que la rasante base evaluada con el mismo
    criterio (corte y relleno de la sección completa). Con una malla relativa al terreno daba el doble."""
    terr, al, pf = _caso_ondulado()
    st = SC.SeccionTipo()
    base = SC.calcular(terr, al, pf, st).resumen
    c_base = st.c_corte * base["corte_m3"] + st.c_relleno * base["relleno_compacto_m3"]
    res = OPT.optimizar_rasante(terr, al, "primaria_1_calzada", 70, st, lambdas=[0.0], pulir=False)
    dp0 = [d for d in res["dp"] if abs(d["lam"]) < 1e-12][0]
    assert dp0["costo_nodo"] <= c_base * 1.005, (dp0["costo_nodo"], c_base)


@prueba
def optimizar_monotonia_del_multiplicador():
    """Propiedad de la optimización paramétrica exacta: el balance neto N(λ) de la solución óptima no crece con λ."""
    terr, al, pf = _caso_ondulado()
    res = OPT.optimizar_rasante(terr, al, "primaria_1_calzada", 70, SC.SeccionTipo(),
                                lambdas=np.linspace(-0.8, 1.0, 7), pulir=False)
    d = sorted(res["dp"], key=lambda x: x["lam"])
    N = [x["N"] for x in d]
    assert all(b <= a + 1e-6 for a, b in zip(N[:-1], N[1:])), N
    assert d[0]["Vc"] > d[-1]["Vc"] and d[0]["Vf"] < d[-1]["Vf"]       # λ alto: menos corte, más relleno


@prueba
def optimizar_no_empeora_y_es_consistente():
    terr, al, pf = _caso_ondulado()
    st = SC.SeccionTipo()
    if "opt" not in _CACHE:
        _CACHE["opt"] = OPT.optimizar_rasante(terr, al, "primaria_1_calzada", 70, st)
    res = _CACHE["opt"]
    base = res["candidatos"][0]
    mejor = res["mejor"]
    assert mejor["costo"] <= base["costo"] + 1e-9
    assert res["mejora_pct"] > 2.0, f"mejora {res['mejora_pct']:.2f} %: el pulido por costo real debería mejorar este caso"
    # el costo reportado coincide con una evaluación independiente del perfil final con todas las secciones
    indep = SC.calcular(terr, al, res["pf"], st).resumen["costo"]["total"]
    assert abs(indep - mejor["costo"]) < 1e-6, (indep, mejor["costo"])
    # y el perfil final cumple todo
    pfm = res["pf"]
    assert not pfm.violaciones, pfm.violaciones
    g = np.diff(pfm.z_pvi) / np.diff(pfm.s_pvi)
    assert np.all(np.abs(g) <= pfm.gmax + 1e-9) and np.all(np.abs(g) >= pfm.gmin - 1e-9)
    assert abs(np.max(np.abs(np.diff(pfm.zd) / np.diff(pfm.s)))) <= pfm.gmax + 0.002
    # coherencia del movimiento de tierras: se reduce el excedente respecto de la base
    assert mejor["excedente_m3"] < base["excedente_m3"]


@prueba
def optimizar_factibilidad_rechaza_configuraciones_invalidas():
    """El pulido solo puede aceptar perfiles que cumplan pendiente máx./mín., tangentes mínimas y curvas dentro del perfil."""
    terr, al, pf = _caso_ondulado()
    L_prev = np.array([0.0] + [c["L"] for c in pf.curvas] + [0.0])
    ok = OPT._factible(pf, pf.s_pvi.copy(), pf.z_pvi.copy(), al, 70, L_prev)
    assert ok is not None, "la configuración del propio perfil debe ser factible"
    # 1) un PVI casi encima del siguiente: las curvas verticales se pisan / la tangente queda corta
    sp = pf.s_pvi.copy(); sp[2] = sp[1] + 40.0
    assert OPT._factible(pf, sp, pf.z_pvi.copy(), al, 70, L_prev) is None
    # 2) pendiente por encima de la máxima
    zp = pf.z_pvi.copy(); zp[-1] += 80.0
    assert OPT._factible(pf, pf.s_pvi.copy(), zp, al, 70, L_prev) is None
    # 3) tangente más plana que la mínima de drenaje
    zp = pf.z_pvi.copy(); zp[1] = zp[0] + 0.001 * (pf.s_pvi[1] - pf.s_pvi[0])
    assert OPT._factible(pf, pf.s_pvi.copy(), zp, al, 70, L_prev) is None
    # 4) una curva que se sale del perfil (PVI pegado al extremo)
    sp = pf.s_pvi.copy(); sp[1] = 20.0
    assert OPT._factible(pf, sp, pf.z_pvi.copy(), al, 70, L_prev) is None


@prueba
def optimizar_factibilidad_detecta_conflictos_de_curvas_con_pendientes_legales():
    """Aísla la comprobación de curvas: todas las pendientes son legales (0,6 %-3,4 %) y aun así la configuración
    es inválida por geometría vertical (dos crestas a 200 m: tangente corta; PVI pegado al inicio: curva fuera)."""
    terr, al, pf = _caso_ondulado()
    L_prev = np.full(4, 100.0)
    z0 = 1024.3
    # dos curvas convexas seguidas separadas 200 m (mínimo de tangente a 90 km/h: 250 m)
    sp = np.array([0.0, 450.0, 650.0, 2381.0])
    zp = np.array([z0, z0 + 0.034 * 450, z0 + 0.034 * 450 + 0.006 * 200, 0.0])
    zp[3] = zp[2] - 0.025 * (2381.0 - 650.0)
    g = np.diff(zp) / np.diff(sp)
    assert np.all(np.abs(g) <= pf.gmax) and np.all(np.abs(g) >= pf.gmin), g          # el caso es legal en pendientes
    assert OPT._factible(pf, sp, zp, al, 70, L_prev) is None, "tangente de 100 m entre curvas debe rechazarse"
    # separadas lo suficiente (700 m) la misma geometría sí es factible
    sp2 = np.array([0.0, 450.0, 1150.0, 2381.0])
    zp2 = np.array([z0, z0 + 0.034 * 450, z0 + 0.034 * 450 + 0.006 * 700, 0.0])
    zp2[3] = zp2[2] - 0.025 * (2381.0 - 1150.0)
    assert OPT._factible(pf, sp2, zp2, al, 70, L_prev) is not None
    # PVI a 40 m del inicio con pendientes legales: su curva vertical se saldría del perfil
    sp3 = np.array([0.0, 40.0, 1500.0, 2381.0])
    zp3 = np.array([z0, z0 + 0.04 * 40, z0 + 0.04 * 40 - 0.02 * 1460, 0.0])
    zp3[3] = zp3[2] + 0.02 * (2381.0 - 1500.0)
    g3 = np.diff(zp3) / np.diff(sp3)
    assert np.all(np.abs(g3) <= pf.gmax) and np.all(np.abs(g3) >= pf.gmin), g3
    assert OPT._factible(pf, sp3, zp3, al, 70, L_prev) is None, "curva fuera del perfil debe rechazarse"


@prueba
def optimizar_flag_desactivable_y_resultado_integrado():
    from .main import disenar_via
    with tempfile.TemporaryDirectory() as d:
        r = disenar_via(terreno_ondulado(), (100, 150), (1850, 1700), "primaria_1_calzada", vtr=70, terreno="ondulado",
                        salida=d, verbose=False, optimizar_rasante=False, replantear=False)
        assert r["opt"] is None and r["rp"] is None and not os.path.exists(os.path.join(d, "optimizacion_rasante.csv"))
        assert r["pf"].modo == "eje"
    with tempfile.TemporaryDirectory() as d:
        r = disenar_via(terreno_ondulado(), (100, 150), (1850, 1700), "primaria_1_calzada", vtr=70, terreno="ondulado",
                        salida=d, verbose=False, optimizar_rasante=True, replantear=False)
        assert r["opt"] is not None
        for f in ("optimizacion_rasante.csv", "optimizacion_rasante.png"):
            assert os.path.getsize(os.path.join(d, f)) > 0, f
        # el perfil y las secciones exportados son los de la mejor alternativa
        assert abs(r["sec"].resumen["costo"]["total"] - r["opt"]["mejor"]["costo"]) < 1e-6
        assert r["sec"].resumen["costo"]["total"] <= r["opt"]["candidatos"][0]["costo"] + 1e-9
        errores, info = V.validar_landxml(os.path.join(d, "via.xml"))
        assert not errores, errores


# ------------------------------------------------------------------ replanteo del trazado según la vía
def terreno_colina(H, sig=250.0, L=2000.0, n=30000, seed=3):
    rng = np.random.default_rng(seed)
    x, y = rng.uniform(0, L, n), rng.uniform(0, L, n)
    z = 1000 + H * np.exp(-((x - 1000) ** 2 + (y - 1000) ** 2) / (2 * sig ** 2))
    return np.column_stack([x, y, z + rng.normal(0, 0.02, n)])


def _costo_camino(G, nodos):
    """Suma de pesos del grafo a lo largo de un camino de nodos (inf si falta una arista)."""
    tot = 0.0
    for a, b in zip(nodos[:-1], nodos[1:]):
        w = G[a, b]
        if w == 0:
            return np.inf
        tot += w
    return float(tot)


@prueba
def replanteo_tabla_costo_formula_cerrada():
    st = SC.SeccionTipo(bombeo=0.0, berma_pend=0.0, e_pav=0.0, cuneta_ancho=0.0, berma=1.8, tc_h=1.0, tf_h=1.5,
                        c_corte=1.0, c_relleno=1.4)
    sg, ct = LC.tabla_costo_transversal(st)
    b = 3.65 + 1.8
    for s in (0.05, 0.10, 0.20):
        cut = s * b * b / 2 + 0.5 * s * b * (b * s * st.tc_h / (1 - s * st.tc_h))
        fil = s * b * b / 2 + 0.5 * s * b * (b * s * st.tf_h / (1 - s * st.tf_h))
        esperado = st.c_corte * cut + st.c_relleno * fil
        i = int(round(s / 0.01))
        assert abs(ct[i] - esperado) / esperado < 0.01, (s, ct[i], esperado)
    fin = ct[np.isfinite(ct)]
    assert np.all(np.diff(fin) >= -1e-9), "el costo transversal debe crecer con la pendiente lateral"
    assert not np.isfinite(ct[int(round(0.70 / 0.01))]), "con talud de relleno 1.5H:1V no cierra a 70 % de pendiente"
    # secciones más anchas cuestan más a cualquier pendiente
    ang = LC.tabla_costo_transversal(SC.SeccionTipo(calzada=6.0, berma=0.5))[1]
    anc = LC.tabla_costo_transversal(SC.SeccionTipo(calzada=10.9))[1]
    assert np.all(ang[1:30] < anc[1:30])


@prueba
def replanteo_ruta_respeta_pendiente_y_alfa0_es_la_mas_corta():
    rng = np.random.default_rng(5)
    x, y = rng.uniform(0, 2000, 20000), rng.uniform(0, 2000, 20000)
    terr = T.Terreno(np.column_stack([x, y, 1000 + 0.15 * y]), celda_tin=20.0)
    st = SC.SeccionTipo()
    r0 = LC.ruta_costo_minimo(terr, (100, 200), (1900, 1900), 0.07, st, alfa=0.0, h=10.0, devolver_grafo=True)
    r2 = LC.ruta_costo_minimo(terr, (100, 200), (1900, 1900), 0.07, st, alfa=2.0, h=10.0, devolver_grafo=True)
    for P, info in (r0, r2):
        seg = np.hypot(*np.diff(P[:, :2], axis=0).T)
        pend = np.abs(np.diff(P[:, 2])) / np.maximum(seg, 1e-9)
        medios = pend[1:-1]                                           # los tramos extremos incluyen el ajuste al punto real
        assert medios.max() <= 0.07 * LC.HOLGURA_PENDIENTE * 1.03, f"pendiente longitudinal {medios.max():.3f} supera la máxima con holgura"
        assert info["longitud"] > np.hypot(1800, 1700) * 1.04, "la ruta recta (10,6 %) no es factible: debe alargarse"
    assert r0[1]["longitud"] <= r2[1]["longitud"] + 2 * 10.0, "α=0 minimiza la longitud"


@prueba
def replanteo_ruta_reacciona_al_ancho_cuando_hay_alternativas():
    """Colina empinada que la vía debe rodear: la sección ancha pasa más lejos (laderas más suaves) y cada ruta es la
    óptima para su propia sección (óptimo cruzado en el mismo grafo)."""
    terr = T.Terreno(terreno_colina(150), celda_tin=20.0)
    stA = SC.SeccionTipo(calzada=6.0, berma=0.5)
    stB = SC.SeccionTipo(calzada=10.9)
    PA, iA = LC.ruta_costo_minimo(terr, (100, 1000), (1900, 1000), 0.07, stA, alfa=1.0, h=10.0, devolver_grafo=True)
    PB, iB = LC.ruta_costo_minimo(terr, (100, 1000), (1900, 1000), 0.07, stB, alfa=1.0, h=10.0, devolver_grafo=True)
    dA = np.min(np.hypot(PA[:, 0] - 1000, PA[:, 1] - 1000))
    dB = np.min(np.hypot(PB[:, 0] - 1000, PB[:, 1] - 1000))
    assert dB > dA + 50, (dA, dB)
    assert iB["longitud"] > iA["longitud"] + 50, (iA["longitud"], iB["longitud"])
    # óptimo cruzado: el camino de A no es mejor que el de B bajo el modelo de B, y viceversa
    assert abs(_costo_camino(iA["grafo"], iA["nodos"]) - iA["costo_modelo"]) < 1e-6
    assert _costo_camino(iB["grafo"], iA["nodos"]) >= iB["costo_modelo"] - 1e-6
    assert _costo_camino(iA["grafo"], iB["nodos"]) >= iA["costo_modelo"] - 1e-6
    # y la ruta es una sola para cada sección: sin ancho de por medio el costo del modelo sube con el ancho
    assert iB["costo_modelo"] > iA["costo_modelo"]


@prueba
def alineamiento_variantes_validas_unicas_y_consistentes_con_disenar():
    terr = T.Terreno(terreno_ondulado(), celda_tin=20.0)
    P, _ = LC.linea_de_ceros(terr, (100, 150), (1850, 1700), 0.07, dh=5.0)
    vs = A.variantes(P, 70, tol_pi=25)
    assert vs and all(v.valido for v in vs)
    claves = [tuple(np.round(v.pis, 0).ravel()) + tuple(round(c["R"]) for c in v.curvas) for v in vs]
    assert len(set(claves)) == len(claves), "variantes duplicadas"
    al = A.disenar(P, 70, tol_pi=25)
    k = tuple(np.round(al.pis, 0).ravel()) + tuple(round(c["R"]) for c in al.curvas)
    assert k in set(claves), "el diseño por fidelidad debe ser una de las variantes"


@prueba
def replanteo_elige_menor_costo_con_prioridad_normativa():
    if "rp" not in _CACHE:
        terr = T.Terreno(terreno_ondulado(), celda_tin=20.0)
        _CACHE["rp"] = (terr, SC.SeccionTipo(), RP.replantear(terr, (100, 150), (1850, 1700), "primaria_1_calzada", 70, SC.SeccionTipo()))
    terr, st, rp = _CACHE["rp"]
    validas = [c for c in rp["candidatas"] if c["valida"]]
    assert len(rp["candidatas"]) >= 4 and validas
    mejor = rp["mejor"]
    assert all((mejor["obs"], mejor["costo"]) <= (c["obs"], c["costo"]) for c in validas)
    ref = rp["referencia"]
    if ref is not None:
        assert (mejor["obs"], mejor["costo"]) <= (ref["obs"], ref["costo"])
    # el costo reportado coincide con una evaluación independiente del diseño elegido
    indep = SC.calcular(terr, mejor["al"], mejor["pf"], st).resumen["costo"]["total"]
    assert abs(indep - mejor["costo"]) < 1e-6
    assert mejor["al"].valido and not mejor["pf"].violaciones
    assert all(c["variantes"] >= 1 for c in validas)


@prueba
def replanteo_integrado_en_pipeline():
    from .main import disenar_via
    with tempfile.TemporaryDirectory() as d:
        r = disenar_via(terreno_ondulado(), (100, 150), (1850, 1700), "primaria_1_calzada", vtr=70, terreno="ondulado",
                        salida=d, verbose=False, optimizar_rasante=False, replantear=True, optimizar_trazado=False)
        assert r["rp"] is not None and os.path.getsize(os.path.join(d, "replanteo_rutas.csv")) > 0
        assert abs(r["sec"].resumen["costo"]["total"] - r["rp"]["mejor"]["costo"]) < 1e-6
        errores, _ = V.validar_landxml(os.path.join(d, "via.xml"))
        assert not errores, errores


@prueba
def alineamiento_pi_colineales_y_muchos_quiebres_pequenos():
    """Regresión: ruta con decenas de quiebres casi colineales (raster) y PI exactamente colineales no deben
    romper el diseño (antes: ZeroDivisionError por deflexión 0 tras agotar las iteraciones de limpieza)."""
    # (a) muchos quiebres de ~0.5-1° que la limpieza debe eliminar todos aunque superen 40
    x = np.arange(0.0, 3000.0, 15.0)
    y = 0.2 * np.sin(np.arange(len(x)) * 2.1)                  # ondulación que con tol bajo deja cientos de vértices
    P = np.column_stack([x, y, np.zeros(len(x))])
    pis = A.ajustar_pis(P[:, :2], tol=0.01)
    if len(pis) > 2:
        d = A._deflexiones(pis)
        assert np.all(d >= 2.0 - 1e-9), f"quedaron quiebres de {d.min():.3f}° sin eliminar"
    # (b) PI exactamente colineales: el diseño no debe lanzar y no genera curvas
    Q = np.array([[0, 0, 0], [500, 0, 0], [1000, 0, 0], [1500, 0, 0], [2000, 0, 0.0]])
    al = A._disenar_una(Q, 70, tol_pi=1.0)
    assert al.valido and len(al.curvas) == 0
    # (c) una curva real entre tramos colineales sí se diseña
    R_ = np.array([[0, 0, 0], [800, 0, 0], [1600, 0, 0], [2200, 500, 0], [2700, 500, 0.0]])
    al2 = A._disenar_una(R_, 70, tol_pi=5.0)
    assert len(al2.curvas) >= 1


@prueba
def replanteo_ruta_no_usa_nodos_fuera_del_levantamiento():
    """Regresión: la ruta se pegaba al borde del terreno aprovechando cotas extrapoladas (sin datos)."""
    rng = np.random.default_rng(5)
    x, y = rng.uniform(0, 2000, 20000), rng.uniform(0, 2000, 20000)
    terr = T.Terreno(np.column_stack([x, y, 1000 + 0.15 * y]), celda_tin=20.0)
    st = SC.SeccionTipo()
    P, info = LC.ruta_costo_minimo(terr, (100, 200), (1900, 1900), 0.07, st, alfa=0.0, h=10.0)
    assert np.all(terr.tri.find_simplex(P[:, :2]) >= 0), "la ruta sale del casco del levantamiento"
    assert P[:, 0].max() < terr.xmax - 20 and P[:, 1].max() < terr.ymax - 20 and P[:, 0].min() > terr.xmin + 20
    try:
        LC.ruta_costo_minimo(terr, (100, 200), (1995, 1900), 0.07, st, h=10.0)
    except RuntimeError as e:
        assert "fuera del área levantada" in str(e)
    else:
        raise AssertionError("un destino pegado al borde debe rechazarse")


@prueba
def alineamiento_filtro_de_pi_colineales_se_activa():
    """Prueba directa del filtro defensivo: si el ajuste de PI entrega un PI colineal (deflexión 0), el diseño lo descarta
    en lugar de dividir por cero."""
    from unittest import mock
    pis = np.array([[0, 0], [800, 0], [1600, 0], [2200, 500], [2700, 500.0]])      # el PI (800,0) es colineal
    P = np.column_stack([pis, np.zeros(len(pis))])
    with mock.patch.object(A, "ajustar_pis", return_value=pis):
        al = A._disenar_una(P, 70, tol_pi=5.0)
    assert len(al.pis) == 4, "el PI colineal debe eliminarse"
    assert np.all(A._deflexiones(al.pis) > 0.05)
    assert len(al.curvas) >= 1


@prueba
def replanteo_no_oculta_errores_inesperados():
    """Una inviabilidad esperada (RuntimeError/ValueError) marca la candidata como no válida; cualquier otro error es un
    defecto y debe propagarse (antes un 'except Exception' escondía una división por cero como 'ruta no válida')."""
    from unittest import mock
    terr = T.Terreno(terreno_ondulado(), celda_tin=20.0)
    P = np.array([[100, 150, 0], [900, 600, 0], [1850, 1700.0, 0]])
    args = (terr, "x", P, "primaria_1_calzada", 70, SC.SeccionTipo(), 0.08, 3.65, False, 25.0)
    with mock.patch.object(A, "variantes", side_effect=RuntimeError("inviable")):
        fila = RP._evaluar(*args)
    assert not fila["valida"] and "inviable" in fila["error"]
    with mock.patch.object(A, "variantes", side_effect=ZeroDivisionError("defecto")):
        try:
            RP._evaluar(*args)
        except ZeroDivisionError:
            return
    raise AssertionError("un error inesperado no debe ocultarse como candidata no válida")


# ------------------------------------------------------------------ optimización de PI y radios por costo real
@prueba
def trazado_criterio_de_aceptacion_lexicografico():
    k = lambda o, c: dict(key=(o, c))
    assert TR._mejor(k(0, 200.0), k(1, 100.0)), "menos observaciones gana aunque cueste más"
    assert not TR._mejor(k(1, 50.0), k(0, 100.0)), "más observaciones nunca se acepta aunque cueste menos"
    assert TR._mejor(k(0, 99.0), k(0, 100.0)), "a igualdad de observaciones gana el menor costo"
    assert not TR._mejor(k(0, 100.0 - 1e-9), k(0, 100.0)), "una mejora despreciable no se acepta (evita ciclos)"
    assert not TR._mejor(None, k(0, 100.0))


@prueba
def trazado_reevaluacion_reproduce_el_costo_del_diseno_original():
    terr, al, pf = _caso_ondulado()
    st = SC.SeccionTipo()
    original = SC.calcular(terr, al, pf, st).resumen["costo"]["total"]
    P = np.column_stack([al.pis, np.zeros(len(al.pis))])
    r = TR.evaluar_configuracion(terr, al.pis, [c["R"] for c in al.curvas], P, "primaria_1_calzada", 70, st,
                                 0.08, 3.65, False, al.tol_pi)
    assert r is not None and abs(r["key"][1] - original) < 1e-6, (r["key"], original)
    assert r["key"][0] == len(al.violaciones) + len(pf.violaciones)


@prueba
def trazado_rechaza_configuraciones_inviables():
    """Cada caso inviable se diseña de forma que no lo frene otra cosa que la defensa que se quiere probar:
    (1) PI fuera del levantamiento con curvas perfectamente válidas; (2) curvas solapadas moderadamente dentro del
    terreno. El control positivo (la configuración original) debe aceptarse."""
    terr, al, pf = _caso_ondulado()
    st = SC.SeccionTipo()
    args = ("primaria_1_calzada", 70, st, 0.08, 3.65, False)
    R = [c["R"] for c in al.curvas]
    P = np.column_stack([al.pis, np.zeros(len(al.pis))])
    # control positivo
    ok = TR.evaluar_configuracion(terr, al.pis, R, P, *args)
    assert ok is not None and ok["al"].valido
    # (1) PI fuera del casco del TIN, con alineamiento válido
    pis1 = np.array([[100.0, 100.0], [terr.xmax + 50.0, 1000.0], [1950.0, 1900.0]])
    a1 = A.disenar_desde_pis(pis1, np.column_stack([pis1, np.zeros(3)]), 70, radios=[304.0])
    assert a1.valido and terr.tri.find_simplex(a1.pis[1:2])[0] < 0, "el caso debe ser válido y tener el PI fuera"
    assert TR.evaluar_configuracion(terr, pis1, [304.0], np.column_stack([pis1, np.zeros(3)]), *args) is None
    # (2) el segundo PI a 100 m del primero: las curvas se pisan unos 150 m, todo dentro del terreno
    pis2 = al.pis.copy()
    u = (pis2[2] - pis2[1]) / np.hypot(*(pis2[2] - pis2[1]))
    pis2[2] = pis2[1] + 100.0 * u
    a2 = A.disenar_desde_pis(pis2, P, 70, radios=R)
    assert (not a2.valido) and np.all(terr.tri.find_simplex(a2.pis) >= 0), "el caso debe ser un solape dentro del terreno"
    assert TR.evaluar_configuracion(terr, pis2, R, P, *args) is None


@prueba
def trazado_encuentra_la_mejora_en_colina_con_optimo_conocido():
    """PI sobre el flanco de una colina de 150 m: el costo baja mucho al alejarlo; la búsqueda debe encontrarlo."""
    terr = T.Terreno(terreno_colina(150), celda_tin=20.0)
    st = SC.SeccionTipo()
    pis = np.array([[100.0, 1000.0], [1000.0, 850.0], [1900.0, 1000.0]])
    P = np.column_stack([pis, np.zeros(3)])
    al0 = A.disenar_desde_pis(pis, P, 70)
    res = TR.optimizar_trazado(terr, al0, P, "primaria_1_calzada", 70, st, max_evals=80)
    assert res["key1"][0] <= res["key0"][0] and res["key1"][1] <= res["key0"][1]
    assert res["mejora_pct"] > 50.0, res["mejora_pct"]
    d0 = np.hypot(*(pis[1] - np.array([1000.0, 1000.0])))
    d1 = np.hypot(*(res["al"].pis[1] - np.array([1000.0, 1000.0])))
    assert d1 > d0 + 200.0, (d0, d1)
    assert np.allclose(res["al"].pis[[0, -1]], pis[[0, -1]]), "BP y EP no se mueven"
    # el costo reportado coincide con una evaluación independiente del diseño final
    indep = SC.calcular(terr, res["al"], res["pf"], st).resumen["costo"]["total"]
    assert abs(indep - res["key1"][1]) < 1e-6
    assert not res["al"].violaciones and not res["pf"].violaciones and res["al"].valido


@prueba
def trazado_no_empeora_en_terreno_ondulado_y_respeta_el_presupuesto():
    terr, al, pf = _caso_ondulado()
    st = SC.SeccionTipo()
    P = np.column_stack([al.pis, np.zeros(len(al.pis))])
    res = TR.optimizar_trazado(terr, al, P, "primaria_1_calzada", 70, st, max_evals=25)
    assert res["n_eval"] <= 25, res["n_eval"]
    assert res["key1"] <= res["key0"] and res["mejora_pct"] >= 0.0
    costos = [c for _, (_, c) in res["historial"]]
    assert all(b <= a + 1e-6 for a, b in zip(costos[:-1], costos[1:])), "el historial debe ser no creciente"
    assert res["al"].valido


@prueba
def trazado_integrado_en_pipeline():
    from .main import disenar_via
    kw = dict(vtr=70, terreno="ondulado", verbose=False, optimizar_rasante=False, replantear=False)
    with tempfile.TemporaryDirectory() as d0, tempfile.TemporaryDirectory() as d:
        sin = disenar_via(terreno_ondulado(), (100, 150), (1850, 1700), "primaria_1_calzada", salida=d0,
                          optimizar_trazado=False, **kw)
        con = disenar_via(terreno_ondulado(), (100, 150), (1850, 1700), "primaria_1_calzada", salida=d,
                          optimizar_trazado=True, **kw)
        assert len(sin["al"].pis) > 2 and con["trz"] is not None
        assert con["sec"].resumen["costo"]["total"] <= sin["sec"].resumen["costo"]["total"] + 1e-6, "mover PI no puede empeorar"
        assert os.path.getsize(os.path.join(d, "optimizacion_trazado.csv")) > 0
        import ezdxf
        capas = {l.dxf.name for l in ezdxf.readfile(os.path.join(d, "via.dxf")).layers}
        if con["trz"]["mejora_pct"] > 0:
            assert "EJE_ANTES_DE_MOVER_PI" in capas
        errores, _ = V.validar_landxml(os.path.join(d, "via.xml"))
        assert not errores, errores


@prueba
def optimizar_rasante_sin_solucion_dp_no_aborta():
    """Regresión (topografía real): si la programación dinámica no encuentra rasante factible, se conserva la base."""
    from unittest import mock
    terr, al, pf = _caso_ondulado()
    st = SC.SeccionTipo()
    base = (pf, SC.calcular(terr, al, pf, st))
    with mock.patch.object(OPT, "dp_malla", side_effect=RuntimeError("Sin rasante factible")):
        res = OPT.optimizar_rasante(terr, al, "primaria_1_calzada", 70, st, base=base, pulir=False)
    assert res["mejor"]["lam"] is None and res["mejora_pct"] == 0.0 and res["pf"] is pf
    assert len(res["candidatos"]) == 1, "solo queda la base cuando ningún candidato DP es factible"


@prueba
def trazado_cambio_de_radios_acorta_la_via_en_terreno_plano():
    """En terreno llano un radio mayor corta más la esquina (R·(2·tan(Δ/2) − Δ) menos de longitud): la búsqueda solo
    de radios debe aumentarlos mientras quepan, sin mover los PI."""
    rng = np.random.default_rng(11)
    x, y = rng.uniform(0, 2000, 15000), rng.uniform(0, 2000, 15000)
    terr = T.Terreno(np.column_stack([x, y, np.full(15000, 1000.0)]), celda_tin=20.0)
    st = SC.SeccionTipo()
    pis = np.array([[100.0, 1000.0], [1000.0, 1000.0], [1800.0, 1600.0]])          # deflexión ≈ 37°, tangentes largas
    P = np.column_stack([pis, np.zeros(3)])
    al0 = A.disenar_desde_pis(pis, P, 70, radios=[304.0])
    assert al0.valido and not al0.violaciones
    res = TR.optimizar_trazado(terr, al0, P, "primaria_1_calzada", 70, st, mover_pi=False, max_evals=40)
    assert np.allclose(res["al"].pis, pis), "con mover_pi=False los PI no se mueven"
    assert res["al"].curvas[0]["R"] > al0.curvas[0]["R"] * 1.2, (al0.curvas[0]["R"], res["al"].curvas[0]["R"])
    assert res["al"].longitud < al0.longitud - 1.0, (al0.longitud, res["al"].longitud)
    assert res["key1"][1] < res["key0"][1] and res["key1"][0] <= res["key0"][0]
    # y con ambas familias desactivadas no cambia nada
    quieto = TR.optimizar_trazado(terr, al0, P, "primaria_1_calzada", 70, st, mover_pi=False, mover_radios=False)
    assert quieto["mejora_pct"] == 0.0 and quieto["n_eval"] == 1


@prueba
def alineamiento_curva_de_retorno_es_inviable_y_no_revienta():
    """Regresión: una deflexión cercana a 180° (horquilla) hacía Ts = R·tan(Δ/2) ~ 1e17 m; el diseño intentaba muestrear esa
    geometría y agotaba la memoria (MemoryError). Debe declararse inviable con un mensaje claro."""
    # (fin, ¿se corta antes de construir la geometría?): ~177° da una tangente de ~10 km (se construye y queda inválida por
    # solape); ~179.7° y 180° dan tangentes de >10^5 m y se declaran inviables sin construir nada
    for fin, se_corta in (((0.0, 60.0), False), ((0.0, 5.0), True), ((-300.0, 0.0), True)):
        pis = np.array([[0.0, 0.0], [1000.0, 0.0], list(fin)])
        al = A.disenar_desde_pis(pis, np.column_stack([pis, np.zeros(3)]), 40)
        assert not al.valido and any("curva de retorno" in m for m in al.solapes), al.solapes
        if se_corta:
            assert al.elementos == [] and al.muestrear().shape == (0, 3)
        else:
            assert al.elementos and np.isfinite(al.longitud) and al.muestrear().shape[1] == 3
    # el resto de la cadena lo trata como inviable y no como defecto
    terr = T.Terreno(terreno_plano_inclinado(0.0, L=3000.0, n=8000), celda_tin=20.0)
    P = np.array([[100.0, 100.0, 0.0], [1100.0, 100.0, 0.0], [100.0, 160.0, 0.0]])
    fila = RP._evaluar(terr, "horquilla", P, "primaria_1_calzada", 70, SC.SeccionTipo(), 0.08, 3.65, False, 5.0)
    assert not fila["valida"]
    # la defensa de muestreo también protege si llega una geometría absurda por otra vía
    absurda = A._disenar_una(np.array([[0.0, 0, 0], [1000.0, 0, 0], [2000.0, 500.0, 0]]), 70, tol_pi=5.0)
    absurda.elementos = [dict(tipo="linea", ini=(0.0, 0.0), fin=(2e17, 0.0), s0=0.0, longitud=2e17)]
    try:
        absurda.muestrear()
    except ValueError as e:
        assert "absurda" in str(e)
    else:
        raise AssertionError("debe rechazar elementos rectos de longitud absurda")


@prueba
def alineamiento_ajuste_de_pi_no_produce_duplicados_ni_nan():
    """Regresión: en terreno empinado la división de deflexiones grandes dejaba PI consecutivos idénticos; la división por
    una longitud nula producía NaN y solo un filtro posterior los descartaba 'por casualidad'."""
    P = np.array([[float(v) for v in par.split(",")] for par in _RUTA_HORQUILLAS.split()])
    assert P.shape == (337, 2)
    with np.errstate(invalid="raise", divide="raise"):               # cualquier NaN/división inválida es un error
        for tol in (12.5, 25.0, 37.5, 50.0, 75.0, 100.0, 150.0, 200.0):
            pis = A.ajustar_pis(P, tol)
            assert np.all(np.isfinite(pis))
            assert np.all(np.hypot(*np.diff(pis, axis=0).T) >= A.PI_DUPLICADO), f"PI duplicados con tol={tol}"
            if len(pis) > 2:
                assert np.all(np.isfinite(A._deflexiones(pis)))
            A._quitar_pi_colineales(pis)


@prueba
def alineamiento_pi_coincidentes_son_invalidos_y_sin_nan():
    with np.errstate(invalid="raise", divide="raise"):
        pis = np.array([[0.0, 0.0], [500.0, 0.0], [500.0, 0.0], [1000.0, 400.0]])
        al = A.disenar_desde_pis(pis, np.column_stack([pis, np.zeros(4)]), 70)
        assert not al.valido and any("coincidentes" in m for m in al.solapes)
        assert A._unit(np.array([0.0, 0.0])).tolist() == [0.0, 0.0]
        assert np.all(np.isfinite(A._deflexiones(np.array([[0.0, 0.0], [0.0, 0.0], [10.0, 0.0], [10.0, 10.0]]))))
        # el filtro de colineales descarta los repetidos (también el que coincide con el punto final)
        q = A._quitar_pi_colineales(np.array([[0.0, 0.0], [500.0, 0.0], [500.0, 0.0], [900.0, 300.0], [900.0, 300.2]]))
        assert q.shape[0] == 3 and np.allclose(q[-1], [900.0, 300.2])


@prueba
def replanteo_alfa0_con_aristas_inconstruibles_sin_nan():
    """Regresión: α=0 multiplicaba 0·inf (aristas cuyo talud no cierra) y producía NaN. La arista inconstruible se excluye de
    forma explícita para cualquier α; si no hay ruta construible, el error es claro (no un FloatingPointError)."""
    rng = np.random.default_rng(4)
    x, y = rng.uniform(0, 2000, 20000), rng.uniform(0, 2000, 20000)
    terr = T.Terreno(np.column_stack([x, y, 1000 + 0.80 * y]), celda_tin=20.0)   # 80 %: ningún talud cierra
    st = SC.SeccionTipo()
    with np.errstate(invalid="raise", divide="raise"):
        for alfa in (0.0, 1.0):
            try:
                LC.ruta_costo_minimo(terr, (100, 1000), (1900, 1000), 0.07, st, alfa=alfa, h=10.0)
            except RuntimeError as e:
                assert "Sin ruta factible" in str(e)
            else:
                raise AssertionError("en una ladera de 80 % no debe haber ruta construible")
    # en un terreno moderado, α=0 sigue dando la ruta más corta construible y sin NaN
    terr2 = T.Terreno(terreno_colina(60), celda_tin=20.0)
    with np.errstate(invalid="raise", divide="raise"):
        P0, i0 = LC.ruta_costo_minimo(terr2, (100, 1000), (1900, 1000), 0.07, st, alfa=0.0, h=10.0)
    assert np.isfinite(i0["costo_modelo"]) and i0["longitud"] < 2200


@prueba
def secciones_avisa_cuando_el_movimiento_de_tierras_domina_el_costo():
    """Un diseño 'válido' en terreno muy empinado puede costar decenas de veces el de la vía: debe advertirse."""
    terr = T.Terreno(terreno_colina(150), celda_tin=20.0)
    pis = np.array([[100.0, 1000.0], [1000.0, 850.0], [1900.0, 1000.0]])        # PI sobre el flanco de la colina
    P = np.column_stack([pis, np.zeros(3)])
    al = A.disenar_desde_pis(pis, P, 70)
    pf = PF.disenar_perfil(terr, al, "primaria_1_calzada", 70)
    sec = SC.calcular(terr, al, pf, SC.SeccionTipo())
    c = sec.resumen["costo"]
    assert c["corte"] + c["relleno"] > SC.RELACION_TIERRAS_ALERTA * c["via"]
    assert any("Movimiento de tierras" in a for a in sec.avisos)
    # y el caso de referencia, razonable, no dispara el aviso
    terr0, al0, pf0 = _caso_ondulado()
    sec0 = SC.calcular(terr0, al0, pf0, SC.SeccionTipo())
    assert not any("Movimiento de tierras" in a for a in sec0.avisos)


@prueba
def topografia_detecta_zonas_sin_datos_sin_falsos_positivos():
    rng = np.random.default_rng(8)
    x, y = rng.uniform(0, 2000, 25000), rng.uniform(0, 2000, 25000)
    con_hueco = ~((x > 800) & (x < 1200) & (y > 800) & (y < 1200))
    pts = np.column_stack([x, y, 1000 + 0.02 * x])
    terr_h = T.Terreno(pts[con_hueco], celda_tin=20.0)
    vac, um = terr_h.en_vacio([1000.0, 300.0, 1700.0, 830.0], [1000.0, 300.0, 1700.0, 1000.0])
    assert vac.tolist() == [True, False, False, False] and um >= 50.0, (vac, um)   # (830,1000) está a 30 m del borde: cubierto
    terr_ok = T.Terreno(pts, celda_tin=20.0)                      # sin hueco: ningún punto debe marcarse
    gx, gy = np.meshgrid(np.linspace(100, 1900, 60), np.linspace(100, 1900, 60))
    assert not terr_ok.en_vacio(gx.ravel(), gy.ravel())[0].any(), "falsos positivos en un levantamiento completo"
    terr, al, pf = _caso_ondulado()
    assert not terr.en_vacio(*al.muestrear(5.0)[:, :2].T)[0].any()


@prueba
def topografia_vacios_en_datos_de_curvas_de_nivel_sin_falsos_positivos():
    """Curvas de nivel: densas a lo largo (1 m) y separadas 25 m entre sí. La longitud de las aristas del TIN marcaba como
    'sin datos' casi la mitad del levantamiento; el criterio de distancia al dato no."""
    filas = []
    for k, yc in enumerate(np.arange(0.0, 2001.0, 25.0)):
        xs = np.arange(0.0, 2001.0, 1.0)
        xs = xs[~((xs > 800) & (xs < 1200) & (yc > 800) & (yc < 1200))]            # hueco de 400 x 400 m
        filas.append(np.column_stack([xs, np.full(len(xs), yc), np.full(len(xs), 1000.0 + 1.25 * k)]))
    terr = T.Terreno(np.vstack(filas), celda_tin=2.0)
    gx, gy = np.meshgrid(np.linspace(50, 1950, 120), np.linspace(50, 1950, 120))
    Q = np.column_stack([gx.ravel(), gy.ravel()])
    en_hueco = (Q[:, 0] > 800) & (Q[:, 0] < 1200) & (Q[:, 1] > 800) & (Q[:, 1] < 1200)
    cercano_al_hueco = (Q[:, 0] > 750) & (Q[:, 0] < 1250) & (Q[:, 1] > 750) & (Q[:, 1] < 1250)
    vac = terr.en_vacio(Q[:, 0], Q[:, 1])[0]
    assert not vac[~cercano_al_hueco].any(), "falsos positivos fuera del hueco (zona cubierta por curvas)"
    assert vac[(Q[:, 0] > 900) & (Q[:, 0] < 1100) & (Q[:, 1] > 900) & (Q[:, 1] < 1100)].all(), "el centro del hueco debe marcarse"
    assert vac[en_hueco].mean() > 0.4


@prueba
def secciones_avisa_cuando_el_eje_cruza_una_zona_sin_datos():
    rng = np.random.default_rng(8)
    x, y = rng.uniform(0, 2000, 25000), rng.uniform(0, 2000, 25000)
    pts = np.column_stack([x, y, 1000 + 0.02 * x])
    hueco = pts[~((pts[:, 0] > 800) & (pts[:, 0] < 1200) & (pts[:, 1] > 800) & (pts[:, 1] < 1200))]
    P = np.array([[100.0, 1000.0, 0.0], [1900.0, 1000.0, 0.0]])
    avisos = {}
    for nombre, datos in (("con_hueco", hueco), ("completo", pts)):
        terr = T.Terreno(datos, celda_tin=20.0)
        al = A._disenar_una(P, 70, tol_pi=25.0)
        pf = PF.disenar_perfil(terr, al, "primaria_1_calzada", 70)
        avisos[nombre] = SC.calcular(terr, al, pf, SC.SeccionTipo()).avisos
    a = [m for m in avisos["con_hueco"] if "zonas sin datos" in m]
    assert a, avisos["con_hueco"]
    m = float(a[0].split("≈")[1].split(" m")[0])
    assert 250 < m < 420, f"el hueco es de 400 m y se marcan los que distan más de 50 m de un dato (~300 m): se reportó {m}"
    assert not any("zonas sin datos" in m for m in avisos["completo"])


# ------------------------------------------------------------------ datos reales: CSV, curvas de nivel y zonas sin datos
@prueba
def topografia_cargar_csv_con_y_sin_encabezado_y_otros_separadores():
    datos = np.array([[100.5, 200.25, 10.0], [101.5, 201.0, 11.5], [102.0, 202.5, 9.75], [103.0, 203.5, 8.0]])
    with tempfile.TemporaryDirectory() as d:
        casos = {"a.csv": "Este,Norte,Cota\n", "b.csv": "", "c.txt": "X;Y;Z\n", "d.csv": "cota,norte,este\n"}
        for nombre, enc in casos.items():
            sep = ";" if nombre == "c.txt" else ","
            cuerpo = "\n".join(sep.join(f"{v}" for v in fila) for fila in datos)
            if nombre == "d.csv":                                    # columnas en otro orden, reconocidas por nombre
                cuerpo = "\n".join(sep.join(f"{v}" for v in fila[::-1]) for fila in datos)
            ruta = os.path.join(d, nombre)
            open(ruta, "w").write(enc + cuerpo + "\n")
            out = T.cargar_puntos(ruta)
            assert out.shape == (4, 3) and np.allclose(out, datos), (nombre, out)


@prueba
def replanteo_ruta_esquiva_las_zonas_sin_datos():
    rng = np.random.default_rng(12)
    x, y = rng.uniform(0, 2000, 40000), rng.uniform(0, 2000, 40000)
    pts = np.column_stack([x, y, 1000 + 0.01 * x])
    hueco = pts[~((pts[:, 0] > 700) & (pts[:, 0] < 1300) & (pts[:, 1] > 700) & (pts[:, 1] < 1300))]   # 600 x 600 m
    terr = T.Terreno(hueco, celda_tin=20.0)
    st = SC.SeccionTipo()
    P, info = LC.ruta_costo_minimo(terr, (100, 1000), (1900, 1000), 0.07, st, alfa=1.0, h=10.0)
    assert not terr.en_vacio(P[:, 0], P[:, 1])[0].any(), "la ruta no debe pasar por una zona sin datos"
    # rodear el cuadrado sin datos efectivo (~500 m) exige al menos 2·√(650²+250²) + 500 ≈ 1.893 m; la malla añade un poco
    assert 1890 < info["longitud"] < 1800 * 1.12, info["longitud"]
    terr_ok = T.Terreno(pts, celda_tin=20.0)
    P2, i2 = LC.ruta_costo_minimo(terr_ok, (100, 1000), (1900, 1000), 0.07, st, alfa=1.0, h=10.0)
    assert i2["longitud"] < 1800 * 1.02 and info["longitud"] > i2["longitud"] + 80
    try:
        LC.ruta_costo_minimo(terr, (1000, 1000), (1900, 1000), 0.07, st, h=10.0)
    except RuntimeError as e:
        assert "zona sin datos" in str(e) or "fuera del área" in str(e)
    else:
        raise AssertionError("un origen en una zona sin datos debe rechazarse")


@prueba
def trazado_cruzar_una_zona_sin_datos_cuenta_como_observacion():
    rng = np.random.default_rng(12)
    x, y = rng.uniform(0, 2000, 40000), rng.uniform(0, 2000, 40000)
    pts = np.column_stack([x, y, 1000 + 0.01 * x])
    hueco = pts[~((pts[:, 0] > 700) & (pts[:, 0] < 1300) & (pts[:, 1] > 700) & (pts[:, 1] < 1300))]
    terr = T.Terreno(hueco, celda_tin=20.0)
    st = SC.SeccionTipo()
    args = ("primaria_1_calzada", 70, st, 0.08, 3.65, False)
    cruza = np.array([[100.0, 1000.0], [1000.0, 1000.0], [1900.0, 1100.0]])            # atraviesa el hueco
    rodea = np.array([[100.0, 1000.0], [1000.0, 500.0], [1900.0, 1000.0]])             # lo evita por el sur
    r1 = TR.evaluar_configuracion(terr, cruza, [None], np.column_stack([cruza, np.zeros(3)]), *args)
    r2 = TR.evaluar_configuracion(terr, rodea, [None], np.column_stack([rodea, np.zeros(3)]), *args)
    assert r1 is not None and r2 is not None
    assert r1["sec"].resumen["eje_en_vacio_m"] > 100 and r2["sec"].resumen["eje_en_vacio_m"] == 0.0
    assert r1["key"][0] >= 1 and r2["key"][0] == 0, (r1["key"], r2["key"])
    assert TR._mejor(r2, r1), "el trazado que no inventa terreno se prefiere aunque cueste más"


# ------------------------------------------------------------------ ejecución
def main():
    """python -m diseno_vial.pruebas [-k texto]   (-k: ejecuta solo las pruebas cuyo nombre contiene el texto)"""
    filtro = sys.argv[sys.argv.index("-k") + 1] if "-k" in sys.argv else None
    pruebas = [f for f in _RES if filtro is None or filtro in f.__name__]
    fallos = 0
    for f in pruebas:
        try:
            f()
            print(f"  OK     {f.__name__}")
        except Exception:
            fallos += 1
            print(f"  FALLO  {f.__name__}")
            print("         " + traceback.format_exc(limit=4).replace("\n", "\n         "))
    print(f"\n{len(pruebas) - fallos}/{len(pruebas)} pruebas OK")
    return 1 if fallos else 0


if __name__ == "__main__":
    sys.exit(main())
