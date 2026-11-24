#!/usr/bin/env python3
"""
Extrae topografía de un DXF (convertido desde DWG con LibreDWG: ver compilar_libredwg.sh) a un CSV Este,Norte,Cota.

  python extraer_topografia_dwg.py fuente.dxf --listar
  python extraer_topografia_dwg.py fuente.dxf --capas CVL_CURV_D CVL_CURV_G --salida puntos.csv

Entidades leídas con cota: LWPOLYLINE (grupo 38 = elevación), POLYLINE 3D (z de cada vértice), POINT, 3DFACE, LINE.
NO se pueden leer los objetos propietarios de Civil 3D (superficies TIN, alineamientos, etc.): LibreDWG los salta.
"""
import argparse, collections, sys
import numpy as np, pandas as pd, ezdxf


def vertices(e):
    t = e.dxftype()
    if t == "LWPOLYLINE":
        z = float(e.dxf.elevation)
        return [(x, y, z) for x, y, *_ in e.get_points("xy")]
    if t == "POLYLINE":
        return [(v.dxf.location.x, v.dxf.location.y, v.dxf.location.z) for v in e.vertices]
    if t == "POINT":
        p = e.dxf.location; return [(p.x, p.y, p.z)]
    if t == "3DFACE":
        return [tuple(e.dxf.get(f"vtx{i}", (0, 0, 0))) for i in range(4)]
    if t == "LINE":
        return [tuple(e.dxf.start), tuple(e.dxf.end)]
    return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dxf"); ap.add_argument("--capas", nargs="*"); ap.add_argument("--listar", action="store_true")
    ap.add_argument("--salida", default="puntos_topografia.csv")
    a = ap.parse_args()
    doc = ezdxf.readfile(a.dxf); msp = doc.modelspace()
    if a.listar:
        info = collections.defaultdict(lambda: [0, 0, set(), collections.Counter()])
        for e in msp:
            v = vertices(e)
            if not v: continue
            i = info[e.dxf.layer]; i[0] += 1; i[1] += len(v); i[3][e.dxftype()] += 1
            i[2].add(round(float(np.mean([p[2] for p in v])), 3))
        print(f"{'capa':<34}{'entidades':>10}{'vértices':>10}{'cotas distintas':>17}  rango de cota / tipos")
        for k, (n, nv, zs, tp) in sorted(info.items(), key=lambda kv: -kv[1][1]):
            zs = sorted(zs); print(f"{k:<34}{n:>10}{nv:>10}{len(zs):>17}  [{zs[0]:.2f}, {zs[-1]:.2f}]  {dict(tp)}")
        return
    filas = []
    for e in msp:
        if a.capas and e.dxf.layer not in a.capas: continue
        for x, y, z in vertices(e): filas.append((x, y, z))
    df = pd.DataFrame(filas, columns=["Este", "Norte", "Cota"])
    df.to_csv(a.salida, index=False)
    print(f"{len(df):,} vértices -> {a.salida} | cota [{df.Cota.min():.2f}, {df.Cota.max():.2f}] | "
          f"Este [{df.Este.min():.1f}, {df.Este.max():.1f}] Norte [{df.Norte.min():.1f}, {df.Norte.max():.1f}]")


if __name__ == "__main__":
    sys.exit(main())
