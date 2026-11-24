"""Exportación para Civil 3D: LandXML 1.2 (superficie + alineamiento con espirales + perfil del terreno)
y DXF (curvas de nivel, línea de ceros, eje, PI, bordes, abscisado)."""
from __future__ import annotations
import datetime, math
from xml.sax.saxutils import escape
import numpy as np

from . import alineamiento as A


def _az(p, q):
    return A.azimut((q[0] - p[0], q[1] - p[1]))


def _ne(p):                      # LandXML: "norte este"
    return f"{p[1]:.4f} {p[0]:.4f}"


def _caras_con_datos(terr, max_arista: float = 400.0):
    """Triángulos del TIN con datos: se descartan los de aristas muy largas y los que cubren zonas sin datos
    (centroide a más de `terr.d_vacio` de cualquier punto medido). Sin esto, el casco convexo del TIN se extiende
    como superficie sobre lo que nadie levantó."""
    xy = terr.pts[:, :2]
    f = terr.tri.simplices
    a = xy[f]
    arista = np.max(np.linalg.norm(a - np.roll(a, 1, axis=1), axis=2), axis=1)
    cen = a.mean(axis=1)
    vacio = terr.en_vacio(cen[:, 0], cen[:, 1])[0]
    return f[(arista <= max_arista) & ~vacio]


def landxml_superficie(ruta, terr, nombre_superficie="Terreno", proyecto="Topografía", max_arista=400.0):
    """LandXML 1.2 SOLO con la superficie TIN (sin alineamiento): para importar la topografía extraída."""
    xy, z = terr.pts[:, :2], terr.pts[:, 2]
    s = _caras_con_datos(terr, max_arista)
    usados = np.unique(s)                                 # no se escriben los puntos que ya no pertenecen a ninguna cara
    nuevo = np.full(len(xy), -1, dtype=np.int64); nuevo[usados] = np.arange(len(usados))
    now = datetime.datetime.now()
    L = ['<?xml version="1.0" encoding="UTF-8"?>',
         '<LandXML xmlns="http://www.landxml.org/schema/LandXML-1.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
         'xsi:schemaLocation="http://www.landxml.org/schema/LandXML-1.2 http://www.landxml.org/schema/LandXML-1.2/LandXML-1.2.xsd" '
         f'date="{now:%Y-%m-%d}" time="{now:%H:%M:%S}" version="1.2" language="Spanish" readOnly="false">',
         '  <Units><Metric areaUnit="squareMeter" linearUnit="meter" volumeUnit="cubicMeter" temperatureUnit="celsius" '
         'pressureUnit="HPA" diameterUnit="millimeter" angularUnit="decimal degrees" directionUnit="decimal degrees"/></Units>',
         f'  <Project name="{escape(proyecto)}"/>',
         '  <Application name="diseno_vial" desc="Superficie desde curvas de nivel" manufacturer="-" version="0.7"/>',
         "  <Surfaces>", f'    <Surface name="{escape(nombre_superficie)}"><Definition surfType="TIN">', "      <Pnts>"]
    for i, k in enumerate(usados, start=1):
        L.append(f'        <P id="{i}">{xy[k, 1]:.4f} {xy[k, 0]:.4f} {z[k]:.4f}</P>')
    L.append("      </Pnts>\n      <Faces>")
    for t in s:
        L.append(f"        <F>{nuevo[t[0]] + 1} {nuevo[t[1]] + 1} {nuevo[t[2]] + 1}</F>")
    L.append("      </Faces>\n    </Definition></Surface>\n  </Surfaces>\n</LandXML>")
    with open(ruta, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    return dict(puntos=len(usados), triangulos=len(s), triangulos_descartados=len(terr.tri.simplices) - len(s))


def _elementos_xml(al) -> str:
    out = []
    for e in al.elementos:
        s0 = e["s0"]
        if e["tipo"] == "linea":
            if e["longitud"] < 1e-6:
                continue
            out.append(f'      <Line dir="{_az(e["ini"], e["fin"]):.6f}" length="{e["longitud"]:.4f}" staStart="{s0:.4f}">'
                       f'<Start>{_ne(e["ini"])}</Start><End>{_ne(e["fin"])}</End></Line>')
        elif e["tipo"] == "espiral":
            sp = e["sp"]
            rot = "ccw" if e["sgn"] > 0 else "cw"
            rs, re_ = ("INF", f'{e["R"]:.4f}') if e["entrada"] else (f'{e["R"]:.4f}', "INF")
            out.append(f'      <Spiral length="{e["longitud"]:.4f}" radiusEnd="{re_}" radiusStart="{rs}" rot="{rot}" '
                       f'spiType="clothoid" theta="{math.degrees(sp["theta"]):.6f}" totalY="{sp["Ys"]:.4f}" '
                       f'totalX="{sp["Xs"]:.4f}" tanLong="{sp["tan_larga"]:.4f}" tanShort="{sp["tan_corta"]:.4f}" '
                       f'staStart="{s0:.4f}"><Start>{_ne(e["ini"])}</Start><PI>{_ne(e["pi"])}</PI><End>{_ne(e["fin"])}</End></Spiral>')
        else:
            c = e["curva"]
            sgn, R, dc = e["sgn"], e["R"], e["delta"]
            u1 = c["u1"]
            az_ini = A.azimut(u1) - sgn * math.degrees(c["theta"])
            az_fin = az_ini - sgn * math.degrees(dc)
            rot = "ccw" if sgn > 0 else "cw"
            tan = R * math.tan(dc / 2) if dc > 1e-9 else 0.0
            dirv = np.array([math.sin(math.radians(az_ini)), math.cos(math.radians(az_ini))])
            pi = np.array(e["ini"]) + tan * dirv
            out.append(f'      <Curve rot="{rot}" chord="{2 * R * math.sin(dc / 2):.4f}" crvType="arc" delta="{math.degrees(dc):.6f}" '
                       f'dirStart="{az_ini % 360:.6f}" dirEnd="{az_fin % 360:.6f}" external="{R * (1 / math.cos(dc / 2) - 1):.4f}" '
                       f'length="{e["longitud"]:.4f}" midOrd="{R * (1 - math.cos(dc / 2)):.4f}" radius="{R:.4f}" '
                       f'tangent="{tan:.4f}" staStart="{s0:.4f}"><Start>{_ne(e["ini"])}</Start>'
                       f'<Center>{_ne(e["centro"])}</Center><End>{_ne(e["fin"])}</End><PI>{_ne(pi)}</PI></Curve>')
    return "\n".join(out)


def landxml(ruta, terr, al, nombre_via="Eje 1", nombre_superficie="Terreno", proyecto="Diseño vial",
            max_arista=400.0, ds_perfil=10.0, rasante=None):
    """Escribe LandXML 1.2. `rasante` (opcional): lista [(abscisa, cota)] de PVI para <ProfAlign>."""
    xy, z = terr.pts[:, :2], terr.pts[:, 2]
    s = _caras_con_datos(terr, max_arista)
    now = datetime.datetime.now()
    S = al.muestrear(ds_perfil)
    S_p = np.arange(0, al.longitud + 1e-6, ds_perfil)
    if S_p[-1] < al.longitud - 1e-6:
        S_p = np.append(S_p, al.longitud)
    xp, yp = np.interp(S_p, S[:, 2], S[:, 0]), np.interp(S_p, S[:, 2], S[:, 1])
    zp = terr.z(xp, yp)
    L = []
    L.append('<?xml version="1.0" encoding="UTF-8"?>')
    L.append('<LandXML xmlns="http://www.landxml.org/schema/LandXML-1.2" '
             'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
             'xsi:schemaLocation="http://www.landxml.org/schema/LandXML-1.2 http://www.landxml.org/schema/LandXML-1.2/LandXML-1.2.xsd" '
             f'date="{now:%Y-%m-%d}" time="{now:%H:%M:%S}" version="1.2" language="Spanish" readOnly="false">')
    L.append('  <Units><Metric areaUnit="squareMeter" linearUnit="meter" volumeUnit="cubicMeter" '
             'temperatureUnit="celsius" pressureUnit="HPA" diameterUnit="millimeter" '
             'angularUnit="decimal degrees" directionUnit="decimal degrees"/></Units>')
    L.append(f'  <Project name="{escape(proyecto)}"/>')
    L.append('  <Application name="diseno_vial" desc="Prototipo INVIAS/AASHTO" manufacturer="-" version="0.1"/>')
    L.append("  <Surfaces>")
    L.append(f'    <Surface name="{escape(nombre_superficie)}"><Definition surfType="TIN">')
    L.append("      <Pnts>")
    for i, (p, zz) in enumerate(zip(xy, z), start=1):
        L.append(f'        <P id="{i}">{p[1]:.4f} {p[0]:.4f} {zz:.4f}</P>')
    L.append("      </Pnts>\n      <Faces>")
    for t in s:
        L.append(f"        <F>{t[0] + 1} {t[1] + 1} {t[2] + 1}</F>")
    L.append("      </Faces>\n    </Definition></Surface>\n  </Surfaces>")
    L.append('  <Alignments name="Alineamientos">')
    L.append(f'    <Alignment name="{escape(nombre_via)}" length="{al.longitud:.4f}" staStart="0" desc="VTR {al.vtr:.0f} km/h">')
    L.append("      <CoordGeom>")
    L.append(_elementos_xml(al))
    L.append("      </CoordGeom>")
    L.append(f'      <Profile name="{escape(nombre_via)} - Terreno"><ProfSurf name="{escape(nombre_superficie)}">'
             "<PntList2D>" + " ".join(f"{a:.3f} {b:.3f}" for a, b in zip(S_p, zp)) + "</PntList2D></ProfSurf></Profile>")
    if rasante:
        L.append(f'      <Profile name="{escape(nombre_via)} - Rasante"><ProfAlign name="Rasante">')
        for item in rasante:
            sta, el = item[0], item[1]
            Lc = item[2] if len(item) > 2 else 0.0
            if Lc > 0:
                L.append(f'        <ParaCurve length="{Lc:.3f}">{sta:.3f} {el:.3f}</ParaCurve>')
            else:
                L.append(f"        <PVI>{sta:.3f} {el:.3f}</PVI>")
        L.append("      </ProfAlign></Profile>")
    L.append("    </Alignment>\n  </Alignments>\n</LandXML>")
    with open(ruta, "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    return dict(puntos=len(xy), triangulos=len(s))


def dxf(ruta, curvas, P_ceros, al, origen, destino, intervalo, ancho_calzada, berma, maestra=5, sec=None, ruta_elegida=None,
        al_antes=None):
    import ezdxf
    doc = ezdxf.new("R2018", setup=True); doc.units = 6
    msp = doc.modelspace()
    capas = {"CURVA_MENOR": 8, "CURVA_MAESTRA": 7, "TEXTO_CURVA": 7, "LINEA_CEROS": 1, "EJE": 2,
             "PI": 6, "PUNTOS_CURVA": 4, "BORDE_CALZADA": 3, "BORDE_BERMA": 5, "ABSCISADO": 7, "TEXTO": 7}
    for n, c in capas.items():
        doc.layers.add(n, color=c)
    doc.layers.get("LINEA_CEROS").dxf.linetype = "DASHED"
    for zc, P in curvas:
        mae = abs(zc / (intervalo * maestra) - round(zc / (intervalo * maestra))) < 1e-6
        msp.add_lwpolyline(P.tolist(), dxfattribs={"layer": "CURVA_MAESTRA" if mae else "CURVA_MENOR", "elevation": zc})
        if mae and len(P) > 12:
            m = len(P) // 2
            t = P[min(m + 1, len(P) - 1)] - P[m - 1]
            ang = math.degrees(math.atan2(t[1], t[0]))
            ang = ang - 180 if ang > 90 else ang + 180 if ang < -90 else ang
            msp.add_text(f"{zc:.0f}", height=intervalo * 0.8, dxfattribs={"layer": "TEXTO_CURVA", "insert": (*P[m], zc), "rotation": ang})
    msp.add_lwpolyline(P_ceros[:, :2].tolist(), dxfattribs={"layer": "LINEA_CEROS"})
    if ruta_elegida is not None and ruta_elegida is not P_ceros:
        doc.layers.add("RUTA_ELEGIDA", color=6)
        msp.add_lwpolyline(ruta_elegida[:, :2].tolist(), dxfattribs={"layer": "RUTA_ELEGIDA"})
    S = al.muestrear(2.0)
    if al_antes is not None and al_antes is not al:
        doc.layers.add("EJE_ANTES_DE_MOVER_PI", color=8)
        msp.add_lwpolyline(al_antes.muestrear(2.0)[:, :2].tolist(), dxfattribs={"layer": "EJE_ANTES_DE_MOVER_PI"})
    msp.add_lwpolyline(S[:, :2].tolist(), dxfattribs={"layer": "EJE"})
    tang = np.gradient(S[:, :2], axis=0); tang /= np.linalg.norm(tang, axis=1)[:, None]
    nrm = np.column_stack([-tang[:, 1], tang[:, 0]])
    if sec is not None:                                  # banca real: sobreancho, chaflanes
        from . import secciones as SCm
        for n_, c_ in {"CHAFLAN_CORTE": 5, "CHAFLAN_RELLENO": 4}.items():
            doc.layers.add(n_, color=c_)
        bd = SCm.bordes(sec)
        for k_, capa in (("calzada_izq", "BORDE_CALZADA"), ("calzada_der", "BORDE_CALZADA"),
                         ("berma_izq", "BORDE_BERMA"), ("berma_der", "BORDE_BERMA")):
            msp.add_lwpolyline(bd[k_].tolist(), dxfattribs={"layer": capa})
        for modo, capa in (("corte", "CHAFLAN_CORTE"), ("relleno", "CHAFLAN_RELLENO")):
            for seg in bd["chaflan"][modo]:
                msp.add_lwpolyline(seg.tolist(), dxfattribs={"layer": capa})
    else:
        for off, capa in ((ancho_calzada / 2, "BORDE_CALZADA"), (ancho_calzada / 2 + _b(berma), "BORDE_BERMA")):
            for sg in (1, -1):
                msp.add_lwpolyline((S[:, :2] + sg * off * nrm).tolist(), dxfattribs={"layer": capa})
    paso = max(intervalo * 4, 20.0)
    for st in np.arange(0, al.longitud, paso):
        i = int(np.argmin(np.abs(S[:, 2] - st)))
        a, b = S[i, :2] + nrm[i] * (ancho_calzada / 2 + _b(berma) + paso * 0.1), S[i, :2] - nrm[i] * (ancho_calzada / 2 + _b(berma) + paso * 0.1)
        msp.add_line(a.tolist(), b.tolist(), dxfattribs={"layer": "ABSCISADO"})
        if round(st) % int(paso * 5) == 0:
            msp.add_text(f"K{int(st // 1000)}+{int(st % 1000):03d}", height=paso * 0.2,
                         dxfattribs={"layer": "TEXTO", "insert": (S[i, :2] + nrm[i] * (ancho_calzada / 2 + _b(berma) + paso * 0.3)).tolist(),
                                     "rotation": math.degrees(math.atan2(nrm[i][1], nrm[i][0])) - 90})
    h = max(intervalo, 3.0)
    for k, pi in enumerate(al.pis):
        nombre = "BP" if k == 0 else "EP" if k == len(al.pis) - 1 else f"PI-{k}"
        msp.add_circle(pi, h, dxfattribs={"layer": "PI"})
        msp.add_text(nombre, height=h * 1.5, dxfattribs={"layer": "PI", "insert": (pi[0] + h * 1.5, pi[1] + h * 1.5)})
    for c in al.curvas:
        for nom in ("TE", "SC", "CS", "ST"):
            if c["Le"] == 0 and nom in ("SC", "CS"):
                continue
            msp.add_point(c[nom], dxfattribs={"layer": "PUNTOS_CURVA"})
            msp.add_text(f"{nom}{c['idx']}", height=h, dxfattribs={"layer": "PUNTOS_CURVA", "insert": (c[nom][0] + h, c[nom][1] - h)})
    doc.saveas(ruta)
    return doc


def _b(berma):
    return max(berma) if isinstance(berma, (tuple, list)) else berma


def dxf_perfil(ruta, pf, exag: float = 10.0, z_base: float | None = None):
    """DXF del perfil longitudinal en coordenadas (abscisa, cota) con exageración vertical `exag`
    (la cota dibujada es (z - z_base)*exag; los rótulos muestran la cota real)."""
    import ezdxf
    doc = ezdxf.new("R2018", setup=True); doc.units = 6
    msp = doc.modelspace()
    for n, c in {"TERRENO": 3, "RASANTE": 1, "PVI": 2, "CUADRICULA": 8, "TEXTO": 7, "CURVA_VERTICAL": 30}.items():
        doc.layers.add(n, color=c)
    zb = z_base if z_base is not None else math.floor(min(pf.zt.min(), pf.zd.min()) / 10.0) * 10.0 - 10.0
    zmax = math.ceil(max(pf.zt.max(), pf.zd.max()) / 10.0) * 10.0 + 10.0
    Y = lambda z: (np.asarray(z) - zb) * exag
    S0, S1 = float(pf.s[0]), float(pf.s[-1])
    msp.add_lwpolyline(list(zip(pf.s, Y(pf.zt))), dxfattribs={"layer": "TERRENO"})
    msp.add_lwpolyline(list(zip(pf.s, Y(pf.zd))), dxfattribs={"layer": "RASANTE"})
    z = zb
    while z <= zmax + 1e-6:
        msp.add_line((S0, float(Y(z))), (S1, float(Y(z))), dxfattribs={"layer": "CUADRICULA"})
        msp.add_text(f"{z:.0f}", height=4, dxfattribs={"layer": "TEXTO", "insert": (S0 - 40, float(Y(z)) - 2)})
        z += 10.0
    st = 0.0
    while st <= S1 + 1e-6:
        msp.add_line((st, 0), (st, float(Y(zmax))), dxfattribs={"layer": "CUADRICULA"})
        msp.add_text(f"K{int(st // 1000)}+{int(st % 1000):03d}", height=4,
                     dxfattribs={"layer": "TEXTO", "insert": (st - 12, -14), "rotation": 0})
        st += 100.0
    g = np.diff(pf.z_pvi) / np.diff(pf.s_pvi)
    for i, (a, b) in enumerate(zip(pf.s_pvi, pf.z_pvi)):
        msp.add_circle((float(a), float(Y(b))), 3, dxfattribs={"layer": "PVI"})
        msp.add_text(f"PVI{i} K{int(a // 1000)}+{int(a % 1000):03d} z={b:.2f}", height=3.5,
                     dxfattribs={"layer": "PVI", "insert": (float(a) + 4, float(Y(b)) + 6), "rotation": 90})
    for c in pf.curvas:
        if c["L"] > 0:
            x0, x1 = c["s_pvi"] - c["L"] / 2, c["s_pvi"] + c["L"] / 2
            msp.add_text(f"CV {c['tipo']} L={c['L']:.0f} A={abs(c['A']) * 100:.2f}% K={c['K']:.1f}", height=3.5,
                         dxfattribs={"layer": "CURVA_VERTICAL", "insert": (x0, float(Y(zmax)) + 8)})
            msp.add_line((x0, 0), (x0, float(Y(zmax))), dxfattribs={"layer": "CURVA_VERTICAL"})
            msp.add_line((x1, 0), (x1, float(Y(zmax))), dxfattribs={"layer": "CURVA_VERTICAL"})
    for i, gi in enumerate(g):
        xm = (pf.s_pvi[i] + pf.s_pvi[i + 1]) / 2
        msp.add_text(f"{gi * 100:+.2f}%", height=4, dxfattribs={
            "layer": "TEXTO", "insert": (float(xm), float(Y(np.interp(xm, pf.s_pvi, pf.z_pvi))) - 14)})
    doc.saveas(ruta)
    return doc


def dxf_secciones(ruta, sec, cada: float = 20.0, columnas: int = 4, ancho_celda: float = 70.0, alto_celda: float = 34.0):
    """Hoja de secciones transversales en DXF (1:1), cada `cada` m, en cuadrícula. La cota se dibuja relativa a la
    rasante del eje (eje de la rasante a y = 10 dentro de cada celda); los rótulos dan las cotas reales."""
    import ezdxf
    doc = ezdxf.new("R2018", setup=True); doc.units = 6
    msp = doc.modelspace()
    for n_, c_ in {"TERRENO_SEC": 3, "PLANTILLA_SEC": 1, "TEXTO_SEC": 7, "MARCO_SEC": 8, "EJE_SEC": 2}.items():
        doc.layers.add(n_, color=c_)
    off, ZL, ZR = sec.perfiles["off"], sec.perfiles["ZL"], sec.perfiles["ZR"]
    ks = [int(np.argmin(np.abs(sec.s - st_))) for st_ in np.arange(0, sec.s[-1] + 1e-6, cada)]
    for n, k in enumerate(ks):
        col, fila = n % columnas, n // columnas
        ox, oy = col * ancho_celda, -fila * alto_celda
        base = 10.0
        zr = sec.zr[k]
        dL, dR = sec.izq[k]["d_chaflan"], sec.der[k]["d_chaflan"]
        mL, mR = off <= dL + 3, off <= dR + 3
        tl = [(ox - d, oy + base + (z - zr)) for d, z in zip(off[mL], ZL[k][mL])][::-1]
        tr = [(ox + d, oy + base + (z - zr)) for d, z in zip(off[mR], ZR[k][mR])]
        msp.add_lwpolyline(tl + tr[1:], dxfattribs={"layer": "TERRENO_SEC"})
        P = sec.puntos_seccion(k)
        msp.add_lwpolyline([(ox + d, oy + base + (z - zr)) for d, z in P], dxfattribs={"layer": "PLANTILLA_SEC"})
        msp.add_line((ox, oy), (ox, oy + alto_celda - 4), dxfattribs={"layer": "EJE_SEC"})
        msp.add_lwpolyline([(ox - ancho_celda / 2 + 2, oy - 3), (ox + ancho_celda / 2 - 2, oy - 3),
                            (ox + ancho_celda / 2 - 2, oy + alto_celda - 4), (ox - ancho_celda / 2 + 2, oy + alto_celda - 4)],
                           close=True, dxfattribs={"layer": "MARCO_SEC"})
        msp.add_text(f"K{int(sec.s[k] // 1000)}+{int(sec.s[k] % 1000):03d}  rasante {zr:.2f}  terreno {sec.zt_eje[k]:.2f}",
                     height=1.2, dxfattribs={"layer": "TEXTO_SEC", "insert": (ox - ancho_celda / 2 + 3, oy + alto_celda - 7)})
        msp.add_text(f"Ac={sec.Ac[k]:.1f} m2  Af={sec.Af[k]:.1f} m2  {sec.tipo[k]}", height=1.2,
                     dxfattribs={"layer": "TEXTO_SEC", "insert": (ox - ancho_celda / 2 + 3, oy + alto_celda - 9)})
    doc.saveas(ruta)
    return doc
