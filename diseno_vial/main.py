"""Prototipo v0.1: nube de puntos -> TIN -> línea de ceros -> alineamiento horizontal espiralizado
-> LandXML + DXF (Civil 3D).

Ejemplo:
  python -m diseno_vial.main puntos.xlsx --hoja Topografia --escala 10 \\
      --categoria primaria_1_calzada --vtr 60 --origen 100,150 --destino 1750,1500 --equidistancia 5 --salida salida/
"""
from __future__ import annotations
import argparse, csv, json, os, sys
import numpy as np

from . import normas, topografia as T, linea_ceros as LC, alineamiento as A, perfil as PF, secciones as SC, optimizar as OPT, replanteo as RP, trazado as TR, exportar as X


def disenar_via(puntos, origen, destino, categoria, vtr=None, terreno=None, equidistancia=5.0,
                e_max=0.08, celda_tin=None, tol_pi=None, salida="salida", verbose=True,
                material_corte="suelo_comun", seccion=None, optimizar_rasante=True, replantear=True, optimizar_trazado=True):
    """`seccion`: dict opcional con campos de SC.SeccionTipo (tc_h, tf_h, e_pav, cuneta_ancho, factor_terraplen,
    precios...). `material_corte` fija el talud de corte (roca_sana, roca_fracturada, suelo_comun, suelo_blando)
    salvo que `seccion` traiga tc_h. `optimizar_rasante`: busca la rasante de menor costo real (movimiento de
    tierras con balance de masas) partiendo de la rasante base; nunca queda peor que la base. `replantear`: genera
    rutas candidatas (compás y rutas de costo mínimo con la sección real), las diseña completas y elige la de menor
    costo con menos observaciones normativas; si es False se usa solo la línea de ceros por compás.
    `optimizar_trazado`: búsqueda local de la posición de los PI y de los radios por costo real (antes de optimizar
    la rasante); nunca empeora el diseño de partida."""
    if categoria == "primaria_2_calzadas":
        raise NotImplementedError("Calzadas separadas (separador central, dos calzadas) aún no están implementadas: "
                                  "el prototipo diseña una sola calzada. Usa primaria_1_calzada, secundaria o terciaria.")
    os.makedirs(salida, exist_ok=True)
    extension = max(np.ptp(puntos[:, 0]), np.ptp(puntos[:, 1]))
    celda_tin = celda_tin or max(extension / 100.0, 1.0)
    terr = T.Terreno(puntos, celda_tin=celda_tin)
    cl, pend_pred = T.clasificar(terr)
    terreno = terreno or cl
    rango = normas.rango_vtr(categoria, terreno)
    if vtr is None:
        vtr = rango[len(rango) // 2]
    ancho = normas.ancho_calzada(categoria, terreno, vtr)           # valida la combinación
    berma = normas.berma(categoria, terreno, vtr)
    pend = normas.pend_media_max(categoria, vtr)
    if verbose:
        print(f"Terreno: {terreno} (pendiente predominante {pend_pred:.1f} %)   VTR {vtr} km/h (rango {rango})")
        print(f"Calzada {ancho:.2f} m | berma {berma} m | pendiente media máx. {pend:.0%} (tabla 4.1) | d compás = {equidistancia / pend:.0f} m")
    kw = dict(calzada=ancho, n_carriles=2, berma=(max(berma) if isinstance(berma, (tuple, list)) else berma),
              tc_h=SC.MATERIALES_CORTE[material_corte])
    kw.update(seccion or {})
    st = SC.SeccionTipo(**kw)
    tol_pi = tol_pi or max(extension / 80.0, 5.0)
    info, rp, P_compas = None, None, None
    if replantear:
        rp = RP.replantear(terr, origen, destino, categoria, vtr, st, equidistancia=equidistancia, e_max=e_max,
                           ancho_carril=ancho / 2, terreno_plano=(terreno == "plano"), tol_pi=tol_pi)
        mejor = rp["mejor"]
        P, al, pf = mejor["P"], mejor["al"], mejor["pf"]
        P_compas = rp["referencia"]["P"] if rp["referencia"] is not None else None
        if verbose:
            print(f"Replanteo: {len(rp['candidatas'])} rutas candidatas diseñadas completas (rasante base, todas las variantes de PI):")
            for c in rp["candidatas"]:
                if not c["valida"]:
                    print(f"   {c['nombre']:<30} no válida: {c['error']}")
                else:
                    print(f"   {c['nombre']:<30} L={c['longitud']:.0f} m, {c['curvas']} curvas, obs {c['obs']}, costo {c['costo']:,.0f}"
                          f" ({c['variantes']} variantes)" + ("   <-- elegida" if c is mejor else ""))
            ref = rp["referencia"]
            if ref is not None and ref is not mejor:
                print(f"   la ruta elegida cuesta {100 * (ref['costo'] - mejor['costo']) / ref['costo']:.1f} % menos que la línea de ceros por compás")
    else:
        P, info = LC.linea_de_ceros(terr, origen, destino, pend, dh=equidistancia)
        P_compas = P
        if verbose:
            print(f"Línea de ceros: {info['longitud']:.0f} m, pendiente media {info['pend_media']:.1%}, ascensos {info['ascensos']:.0f} m, Bruce X0 = {info['bruce']:.0f} m")
        al = A.disenar(P, vtr, e_max=e_max, ancho_carril=ancho / 2, tol_pi=tol_pi, terreno_plano=(terreno == "plano"))
        if not al.valido:
            raise RuntimeError("No se encontró un alineamiento horizontal válido (las curvas se pisan): "
                               + "; ".join(al.solapes[:3]) + ". Prueba una VTR menor, otra pendiente/equidistancia "
                               "o un trazado con más distancia entre cambios de dirección.")
        pf = PF.disenar_perfil(terr, al, categoria, vtr)
    al_antes, trz = al, None
    if optimizar_trazado and len(al.pis) > 2:
        trz = TR.optimizar_trazado(terr, al, P, categoria, vtr, st, e_max=e_max, ancho_carril=ancho / 2,
                                   terreno_plano=(terreno == "plano"))
        if trz["pf"] is not None and trz["mejora_pct"] > 0:
            al, pf = trz["al"], trz["pf"]
        if verbose and trz["key0"] is not None:
            print(f"Optimización de PI y radios: {trz['n_eval']} evaluaciones; costo (rasante base) {trz['key0'][1]:,.0f} -> "
                  f"{trz['key1'][1]:,.0f} ({trz['mejora_pct']:.1f} % menos)")
    filas = A.tabla_curvas(al)
    if verbose:
        print(f"Alineamiento: {al.longitud:.0f} m, {len(al.curvas)} curvas, {len(al.violaciones)} observaciones normativas")
        for f in filas:
            print("   PI-{PI}: Δ={deflexion_deg:.1f}° V={V_kmh:.0f} R={R_m:.0f} e={peralte_pct:.1f}% Le={Le_m:.0f} Lc={Lc_m:.0f} Ts={Ts_m:.0f} (ajuste {desviacion_media_m:.1f} m)".format(**f))
        for v in al.violaciones:
            print("   OBSERVACIÓN:", v)
    sec = SC.calcular(terr, al, pf, st)
    opt = None
    if optimizar_rasante:
        base_costo = sec.resumen["costo"]["total"]
        opt = OPT.optimizar_rasante(terr, al, categoria, vtr, st, base=(pf, sec))
        pf_base, sec_base = pf, sec
        pf, sec = opt["pf"], opt["sec"]
        if verbose:
            m = opt["mejor"]
            print(f"Optimización de rasante: {len(opt['candidatos'])} alternativas evaluadas; mejor = {m['modo']} | "
                  f"costo {base_costo:,.0f} -> {m['costo']:,.0f} ({opt['mejora_pct']:.1f} % menos)")
            print(f"   excedente {sec_base.resumen['excedente_m3']:.0f} -> {sec.resumen['excedente_m3']:.0f} m³ | "
                  f"corte {sec_base.resumen['corte_m3']:.0f} -> {sec.resumen['corte_m3']:.0f} m³ | "
                  f"relleno {sec_base.resumen['relleno_compacto_m3']:.0f} -> {sec.resumen['relleno_compacto_m3']:.0f} m³")
            if opt["aviso_malla"]:
                print("   AVISO OPTIMIZACIÓN:", opt["aviso_malla"])
    if verbose:
        stp = pf.stats
        print(f"Perfil: {stp['n_pvi']} PVI, {stp['n_curvas']} curvas verticales, pendiente máx. usada {stp['pend_max_usada_pct']:.1f} % (límite {pf.gmax:.0%}), "
              f"corte máx. {stp['corte_max_m']:.1f} m, relleno máx. {stp['relleno_max_m']:.1f} m")
        for r in pf.tabla():
            if r["L_m"] > 0:
                print("   PVI-{PVI} K{k:.3f}: {tipo} A={A_pct:.2f}% V={V_kmh:.0f} L={L_m:.0f} m (Lmin {Lmin_m:.0f}) K={K:.1f}".format(k=r["abscisa"] / 1000, **r))
        for v in pf.violaciones:
            print("   VIOLACIÓN PERFIL:", v)
        for v in pf.avisos:
            print("   AVISO PERFIL:", v)
        r = sec.resumen
        print(f"Secciones: {r['estaciones']} (cada 10 m) | corte {r['corte_m3']:.0f} m³, relleno compacto {r['relleno_compacto_m3']:.0f} m³ "
              f"(equiv. {r['relleno_equiv_m3']:.0f} con f={st.factor_terraplen}), excedente {r['excedente_m3']:.0f} m³, faltante {r['faltante_m3']:.0f} m³")
        print(f"   acarreo longitudinal {r['acarreo_long_m3']:.0f} m³, DMA {r['dma_m']:.0f} m | afectación media {r['ancho_afectacion_medio_m']:.1f} m (máx {r['ancho_afectacion_max_m']:.1f} m) "
              f"| sobreancho máx {r['sobreancho_max_m']:.2f} m | cortes {r['cortes']} / terraplenes {r['terraplenes']} / mixtas {r['mixtas']}")
        print("   costo (unid. relativas): " + ", ".join(f"{k} {v:,.0f}" for k, v in r["costo"].items()))
        for v in sec.avisos:
            print("   AVISO SECCIONES:", v)
    curvas, tr = T.curvas_de_nivel(terr, equidistancia, max_arista=extension / 5)
    ruta_xml = os.path.join(salida, "via.xml")
    ruta_dxf = os.path.join(salida, "via.dxf")
    r = X.landxml(ruta_xml, terr, al, proyecto="Diseño vial prototipo", max_arista=extension / 5, rasante=pf.pvis())
    X.dxf(ruta_dxf, curvas, (P_compas if P_compas is not None else P), al, origen, destino, equidistancia, ancho, berma, sec=sec, ruta_elegida=P, al_antes=al_antes)
    X.dxf_secciones(os.path.join(salida, "secciones.dxf"), sec)
    SC.graficar(sec, os.path.join(salida, "secciones.png"))
    X.dxf_perfil(os.path.join(salida, "perfil.dxf"), pf)
    PF.graficar(pf, os.path.join(salida, "perfil.png"))
    with open(os.path.join(salida, "secciones.csv"), "w", newline="", encoding="utf-8") as f:
        ft = SC.tabla(sec)
        w = csv.DictWriter(f, fieldnames=list(ft[0].keys())); w.writeheader()
        w.writerows([{k: (round(float(v), 3) if isinstance(v, (float, np.floating)) else v) for k, v in r_.items()} for r_ in ft])
    with open(os.path.join(salida, "resumen_movimiento_tierras.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["concepto", "valor"])
        for k, v in sec.resumen.items():
            if k == "costo":
                for kk, vv in v.items():
                    w.writerow([f"costo_{kk}", round(float(vv), 2)])
            else:
                w.writerow([k, round(float(v), 3)])
    if trz is not None and trz["key0"] is not None:
        with open(os.path.join(salida, "optimizacion_trazado.csv"), "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f); w.writerow(["evaluacion", "observaciones", "costo_total_rasante_base"])
            for n_, (o_, c_) in trz["historial"]:
                w.writerow([n_, o_, round(c_, 1)])
    if rp is not None:
        with open(os.path.join(salida, "replanteo_rutas.csv"), "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["ruta", "valida", "variantes_PI", "longitud_m", "curvas", "observaciones", "corte_m3", "relleno_m3",
                        "excedente_m3", "costo_total_rasante_base", "costo_peor_variante", "separacion_media_a_elegida_m", "elegida"])
            for c in rp["candidatas"]:
                if not c["valida"]:
                    w.writerow([c["nombre"], "no", 0, "", "", "", "", "", "", "", "", "", c["error"]]); continue
                w.writerow([c["nombre"], "sí", c["variantes"], round(c["longitud"], 1), c["curvas"], c["obs"], round(c["corte"], 1),
                            round(c["relleno"], 1), round(c["excedente"], 1), round(c["costo"], 1),
                            round(c["costo_peor_variante"], 1), round(c.get("sep_mejor_media", 0.0), 1),
                            "sí" if c is rp["mejor"] else ""])
    if opt is not None:
        OPT.graficar_comparacion(opt, os.path.join(salida, "optimizacion_rasante.png"))
        with open(os.path.join(salida, "optimizacion_rasante.csv"), "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["alternativa", "lambda", "costo_total", "corte_m3", "relleno_m3", "excedente_m3", "faltante_m3",
                        "acarreo_long_m3", "dma_m", "violaciones_perfil", "elegida"])
            for c in opt["candidatos"]:
                w.writerow([c["modo"], "" if c["lam"] is None else round(c["lam"], 3), round(c["costo"], 1),
                            round(c["corte_m3"], 1), round(c["relleno_m3"], 1), round(c["excedente_m3"], 1),
                            round(c["faltante_m3"], 1), round(c["acarreo_long_m3"], 1), round(c["dma_m"], 1),
                            c["violaciones"], "sí" if c is opt["mejor"] else ""])
    with open(os.path.join(salida, "pvi.csv"), "w", newline="", encoding="utf-8") as f:
        filas_pf = pf.tabla()
        w = csv.DictWriter(f, fieldnames=list(filas_pf[0].keys())); w.writeheader()
        w.writerows([{k: (round(float(v), 3) if isinstance(v, (float, np.floating)) else v) for k, v in r_.items()} for r_ in filas_pf])
    with open(os.path.join(salida, "curvas.csv"), "w", newline="", encoding="utf-8") as f:
        campos = ["PI", "deflexion_deg", "V_kmh", "R_m", "peralte_pct", "Le_m", "Lc_m", "Ts_m", "desviacion_media_m"]
        w = csv.DictWriter(f, fieldnames=campos); w.writeheader()          # el encabezado se escribe aunque no haya curvas
        w.writerows([{k: (round(v, 3) if isinstance(v, float) else v) for k, v in r_.items()} for r_ in filas])
    _grafico(os.path.join(salida, "planta.png"), terr, tr, (P_compas if P_compas is not None else P), al, origen, destino, equidistancia, sec=sec, ruta_elegida=P, al_antes=al_antes)
    return dict(terr=terr, P=P, info=info, al=al, pf=pf, sec=sec, opt=opt, rp=rp, trz=trz, filas=filas, xml=ruta_xml, dxf=ruta_dxf, landxml=r,
                parametros=dict(categoria=categoria, terreno=terreno, vtr=vtr, ancho=ancho, berma=berma, pend=pend))


def _grafico(ruta, terr, tr, P, al, origen, destino, dh, sec=None, ruta_elegida=None, al_antes=None):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 9))
    z = terr.pts[:, 2]
    niv = np.arange(np.floor(z.min() / dh) * dh, z.max() + dh, dh)
    ax.tricontourf(tr, z, niv, cmap="terrain", alpha=0.45)
    ax.tricontour(tr, z, niv, colors="0.35", linewidths=0.3)
    ax.tricontour(tr, z, niv[::5], colors="0.1", linewidths=0.9)
    ax.plot(P[:, 0], P[:, 1], "r--", lw=1.2, label="Línea de ceros (compás)")
    if ruta_elegida is not None and ruta_elegida is not P:
        ax.plot(ruta_elegida[:, 0], ruta_elegida[:, 1], color="limegreen", ls=":", lw=1.8, label="Ruta elegida por el replanteo")
    S = al.muestrear(2.0)
    if al_antes is not None and al_antes is not al:
        Sa = al_antes.muestrear(2.0)
        ax.plot(Sa[:, 0], Sa[:, 1], color="0.35", ls="--", lw=1.2, label="Eje antes de mover PI y radios")
    ax.plot(S[:, 0], S[:, 1], "b-", lw=2.2, label="Eje con espirales")
    ax.plot(al.pis[1:-1, 0], al.pis[1:-1, 1], "m^", label="PI")
    if sec is not None:
        bd = SC.bordes(sec)
        for k_ in ("calzada_izq", "calzada_der"):
            ax.plot(bd[k_][:, 0], bd[k_][:, 1], color="k", lw=0.8)
        for k_, c_, et in (("corte", "darkorange", "Chaflán de corte"), ("relleno", "royalblue", "Chaflán de terraplén")):
            for j_, seg in enumerate(bd["chaflan"][k_]):
                ax.plot(seg[:, 0], seg[:, 1], color=c_, lw=1.3, label=et if j_ == 0 else None)
    ax.plot(*origen, "ko", ms=8); ax.plot(*destino, "ks", ms=8)
    ax.set_aspect("equal"); ax.legend(loc="upper left"); ax.set_title("Planta: curvas de nivel, línea de ceros y eje")
    plt.tight_layout(); fig.savefig(ruta, dpi=100); plt.close(fig)


def _xy(s):
    return tuple(float(v) for v in s.split(","))


def cli():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("archivo"); ap.add_argument("--hoja", default="0")
    ap.add_argument("--escala", type=float, default=1.0, help="multiplica X,Y,Z (solo para probar el terreno de ejemplo)")
    ap.add_argument("--categoria", default="primaria_1_calzada", choices=normas.CATEGORIAS)
    ap.add_argument("--terreno", default=None, choices=normas.TERRENOS)
    ap.add_argument("--vtr", type=int, default=None)
    ap.add_argument("--origen", required=True, type=_xy); ap.add_argument("--destino", required=True, type=_xy)
    ap.add_argument("--equidistancia", type=float, default=5.0)
    ap.add_argument("--emax", type=float, default=0.08)
    ap.add_argument("--material", default="suelo_comun", choices=list(SC.MATERIALES_CORTE), help="talud de corte según material")
    ap.add_argument("--talud-terraplen", type=float, default=None, help="H por 1V del terraplén (por defecto 1.5)")
    ap.add_argument("--e-pav", type=float, default=None, help="espesor de la estructura de pavimento (m)")
    ap.add_argument("--factor-terraplen", type=float, default=None, help="factor de terraplén compacto (1.25-1.30)")
    ap.add_argument("--sin-optimizar-rasante", action="store_true", help="usa la rasante base (más rápido)")
    ap.add_argument("--sin-replanteo", action="store_true", help="usa solo la línea de ceros por compás (más rápido)")
    ap.add_argument("--sin-optimizar-trazado", action="store_true", help="no mueve PI ni radios por costo (más rápido)")
    ap.add_argument("--celda-tin", type=float, default=None,
                    help="celda de adelgazado del TIN en m (por defecto extensión/100; con CURVAS DE NIVEL usar 2: "
                         "el adelgazado grueso destruye su precisión)")
    ap.add_argument("--salida", default="salida")
    a = ap.parse_args()
    pts = T.cargar_puntos(a.archivo, a.hoja) * a.escala
    sec_kw = {}
    if a.talud_terraplen is not None: sec_kw["tf_h"] = a.talud_terraplen
    if a.e_pav is not None: sec_kw["e_pav"] = a.e_pav
    if a.factor_terraplen is not None: sec_kw["factor_terraplen"] = a.factor_terraplen
    disenar_via(pts, a.origen, a.destino, a.categoria, a.vtr, a.terreno, a.equidistancia, a.emax, salida=a.salida,
                material_corte=a.material, seccion=sec_kw or None,
                celda_tin=a.celda_tin, optimizar_rasante=not a.sin_optimizar_rasante, replantear=not a.sin_replanteo, optimizar_trazado=not a.sin_optimizar_trazado)


if __name__ == "__main__":
    sys.exit(cli())
