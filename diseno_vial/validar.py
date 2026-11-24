"""Comprobaciones de consistencia del LandXML generado (no sustituye abrirlo en Civil 3D)."""
import math, sys
import xml.etree.ElementTree as ET

NS = {"l": "http://www.landxml.org/schema/LandXML-1.2"}


def _pt(el):
    n, e = map(float, el.text.split())
    return e, n                                   # (este, norte)


def validar_landxml(ruta, tol=0.05):
    raiz = ET.parse(ruta).getroot()
    errores, info = [], {}
    pts = {p.get("id"): tuple(map(float, p.text.split())) for p in raiz.iterfind(".//l:Surface//l:P", NS)}
    caras = [tuple(f.text.split()) for f in raiz.iterfind(".//l:Surface//l:F", NS)]
    info["puntos"], info["caras"] = len(pts), len(caras)
    if any(i not in pts for c in caras for i in c):
        errores.append("Hay caras que referencian puntos inexistentes")
    cg = raiz.find(".//l:Alignment/l:CoordGeom", NS)
    elems = list(cg)
    prev_end, prev_sta, longitud = None, 0.0, 0.0
    for k, el in enumerate(elems):
        tag = el.tag.split("}")[1]
        ini, fin = _pt(el.find("l:Start", NS)), _pt(el.find("l:End", NS))
        L, sta = float(el.get("length")), float(el.get("staStart"))
        if prev_end is not None:
            if math.dist(prev_end, ini) > tol:
                errores.append(f"Discontinuidad antes del elemento {k} ({tag}): {math.dist(prev_end, ini):.3f} m")
            if abs(sta - prev_sta) > tol:
                errores.append(f"Abscisa no encadena en elemento {k} ({tag}): {sta:.3f} vs {prev_sta:.3f}")
        if tag == "Line" and abs(math.dist(ini, fin) - L) > tol:
            errores.append(f"Línea {k}: longitud {L:.3f} vs distancia {math.dist(ini, fin):.3f}")
        if tag == "Curve":
            c = _pt(el.find("l:Center", NS)); R = float(el.get("radius"))
            for nom, p in (("inicio", ini), ("fin", fin)):
                if abs(math.dist(c, p) - R) > tol:
                    errores.append(f"Curva {k}: {nom} a {math.dist(c, p):.3f} m del centro, R={R:.3f}")
            if abs(R * math.radians(float(el.get("delta"))) - L) > tol:
                errores.append(f"Curva {k}: R·Δ ≠ longitud")
        if tag == "Spiral":
            Ls = float(el.get("length")); th = math.radians(float(el.get("theta")))
            R = float(el.get("radiusEnd") if el.get("radiusEnd") != "INF" else el.get("radiusStart"))
            if abs(Ls / (2 * R) - th) > 1e-6:
                errores.append(f"Espiral {k}: θ ≠ Ls/2R")
            if abs(math.hypot(float(el.get("totalX")), float(el.get("totalY"))) - math.dist(ini, fin)) > tol:
                errores.append(f"Espiral {k}: cuerda ≠ distancia inicio-fin")
        prev_end, prev_sta = fin, sta + L
        longitud = sta + L
    # ---- perfil de diseño (ProfAlign)
    pa = raiz.find(".//l:ProfAlign", NS)
    if pa is not None:
        items = []
        for el in pa:
            tag = el.tag.split("}")[1]
            sta, cota = map(float, el.text.split())
            items.append((tag, sta, cota, float(el.get("length", 0.0))))
        info["perfil"] = {"PVI": sum(1 for i in items if i[0] == "PVI"), "ParaCurve": sum(1 for i in items if i[0] == "ParaCurve")}
        st = [i[1] for i in items]
        if any(b <= a for a, b in zip(st[:-1], st[1:])):
            errores.append("Perfil: las abscisas de los PVI no son estrictamente crecientes")
        if items and (items[0][0] != "PVI" or items[-1][0] != "PVI"):
            errores.append("Perfil: el primero y el último elemento deben ser PVI")
        for k, (tag, sta, cota, Lc) in enumerate(items):
            if tag == "ParaCurve":
                if Lc <= 0:
                    errores.append(f"Perfil: curva en K{sta:.1f} con longitud no positiva")
                if k > 0 and items[k - 1][0] == "ParaCurve" and sta - Lc / 2 < items[k - 1][1] + items[k - 1][3] / 2 - tol:
                    errores.append(f"Perfil: las curvas verticales en K{items[k - 1][1]:.1f} y K{sta:.1f} se solapan")
                if sta - Lc / 2 < st[0] - tol or sta + Lc / 2 > st[-1] + tol:
                    errores.append(f"Perfil: la curva en K{sta:.1f} se sale del perfil")
        if items and abs(items[-1][1] - float(raiz.find(".//l:Alignment", NS).get("length"))) > 0.5:
            errores.append("Perfil: la abscisa final no coincide con la longitud del alineamiento")
    ps = raiz.find(".//l:ProfSurf/l:PntList2D", NS)
    if ps is not None and len(ps.text.split()) % 2:
        errores.append("Perfil del terreno: lista de puntos con número impar de valores")
    al = raiz.find(".//l:Alignment", NS)
    if abs(float(al.get("length")) - longitud) > tol:
        errores.append(f"Longitud del alineamiento {al.get('length')} ≠ suma de elementos {longitud:.3f}")
    info["elementos"] = {t: sum(1 for e in elems if e.tag.endswith(t)) for t in ("Line", "Spiral", "Curve")}
    info["longitud"] = longitud
    # el eje debe quedar dentro de la superficie
    xs = [p[1] for p in pts.values()]; ys = [p[0] for p in pts.values()]
    for nom, p in (("inicio", _pt(elems[0].find("l:Start", NS))), ("fin", _pt(elems[-1].find("l:End", NS)))):
        if not (min(xs) <= p[0] <= max(xs) and min(ys) <= p[1] <= max(ys)):
            errores.append(f"El {nom} del eje queda fuera de la superficie")
    return errores, info


if __name__ == "__main__":
    e, i = validar_landxml(sys.argv[1])
    print("Info:", i)
    print("OK: sin errores" if not e else "ERRORES:\n  " + "\n  ".join(e))
    sys.exit(1 if e else 0)
