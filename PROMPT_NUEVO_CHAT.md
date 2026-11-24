# Texto para pegar en el chat nuevo (adjuntar `diseno_vial_v0.8.zip`)

Continúo un proyecto: un programa en Python (paquete `diseno_vial`) que diseña una vía desde una nube de puntos o curvas de nivel
con INVÍAS/AASHTO y exporta LandXML/DXF para Civil 3D. Descomprime el zip y **lee primero `ESTADO.md`** (sección «RETOMAR EN UN CHAT NUEVO»):
ahí están el estado, el entorno, los datos reales de Felipe (`topografia_felipe/`) y los pendientes en orden.

Instrucciones de trabajo: trabaja con cuidado y verificando (cada cambio con su prueba; comprueba que `python -m diseno_vial.pruebas` da 64/64 antes y después);
lanza las corridas largas desacopladas (`setsid nohup ...`) y sondéalas, porque cada orden tiene un límite de 300 s; no borres pruebas;
guarda siempre lo importante en `/mnt/user-data/outputs` (el entorno se reinicia y borra el resto).

Primera tarea: ejecutar el flujo completo sobre `topografia_felipe/puntos_curvas_de_nivel.csv` (secundaria 40 km/h terreno montañoso, origen y destino =
extremos de su línea de ceros, `--celda-tin 2 --equidistancia 1`), revisar el resultado y seguir con el pendiente 2 (primera curva que no cabe).
