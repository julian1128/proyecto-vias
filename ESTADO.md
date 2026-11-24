# RETOMAR EN UN CHAT NUEVO (leer primero)

**Versión 0.8** — 64 pruebas pasan (`python -m diseno_vial.pruebas`, ~2 min; `-k texto` filtra). Proyecto: diseñar una vía desde una nube de
puntos o curvas de nivel (INVÍAS/AASHTO, diapositivas del curso EIA) y exportar LandXML/DXF para Civil 3D.

## Cómo arrancar
1. `pip install numpy scipy pandas openpyxl matplotlib ezdxf` (numpy ≥ 2). Descomprimir el zip y trabajar dentro de la carpeta.
2. Diseño completo: `python -m diseno_vial.main PUNTOS.csv --categoria secundaria --vtr 40 --terreno montañoso --origen=E,N --destino=E,N --equidistancia 1 --celda-tin 2 --salida salida/`
   (`--origen=` con signo igual si la coordenada es negativa; **con curvas de nivel usar siempre `--celda-tin 2`**; opciones `--sin-replanteo`, `--sin-optimizar-trazado`, `--sin-optimizar-rasante` aceleran).
3. **Las corridas largas (más de ~4 min) deben lanzarse desacopladas** (`setsid nohup sh -c '...' &`) y sondearse: cada orden tiene un límite de 300 s.

## ⚠ Entorno
El entorno de ejecución se reinició 3 veces y borró `/home/claude` y `/tmp` (LibreDWG compilado, DXF convertidos, scripts). Solo persiste `/mnt/user-data`.
Todo lo importante está en la carpeta del proyecto. Para leer un DWG: `sh herramientas/compilar_libredwg.sh` (~5 min, deja `dwg2dxf`) →
`dwg2dxf -y -o fuente.dxf archivo.dwg` → `python herramientas/extraer_topografia_dwg.py fuente.dxf --listar` / `--capas ... --salida puntos.csv`.
LibreDWG no lee los objetos propietarios de Civil 3D (superficies, alineamientos): solo geometría estándar (curvas como polilíneas con cota, etc.).

## Datos reales de Felipe (`topografia_felipe/`, origen: `01_LineaCeros_FELIPE.dwg`)
- `puntos_curvas_de_nivel.csv`: 243.215 vértices de 621 curvas de nivel (1 m, cotas 83-148 m, coordenadas locales: Este −2900…−1440, Norte 2740…4614, metros).
  `superficie_terreno.xml`: TIN LandXML listo para Civil 3D (119.908 puntos, 239.566 caras, recortado en zonas sin datos).
  `lineas_de_referencia_dwg.csv`: su `Linea_Ceros_Felipe` (184 vértices, 2.566 m, de (−2338.4, 2947.5) a (−2033.7, 4524.2), ~7 % de pendiente) y otra polilínea (capa 0).
  `drenaje_quebrada_y_cano.csv`: quebrada (19 polilíneas) y caño.
- **Precisión de la topografía**: con solo las curvas pares (2 m) el TIN predice las impares con RMS 0,32 m y error máx. 1,5 m (mejor que RBF y que el TIN adelgazado). Usar el TIN con todos los vértices (celda 2 m).
- **Terreno**: gradiente mediano 23 % (p90 39 %) en la zona levantada; las curvas cubren una franja irregular (~40 % del casco convexo; el resto es zona sin datos).
- **Hallazgos del primer diseño sobre datos reales**: a 70 km/h es inviable (la ladera no admite radios de 168-304 m; casi todas las variantes «curvas solapadas»).
  Con 20-40 km/h en terreno montañoso sí salen alineamientos válidos pero con 1-5 observaciones normativas (ninguno con 0). La línea de Felipe es muy tortuosa
  (deflexión mediana 17,7° por vértice, hasta 156°): no es un alineamiento sino una guía; el diseñador la simplifica con tangentes largas.
- **La última corrida completa se cayó en el optimizador de rasante** (`Sin rasante factible…`: el terreno sube más rápido que el 7 %); ya está corregido (degrada a la rasante base),
  pero **todavía no hay un resultado completo sobre los datos reales**. Es lo primero que hay que hacer en el chat nuevo (secundaria 40 km/h montañoso y/o terciaria 30).

## Pendiente, en orden de prioridad
1. Correr el flujo completo sobre la topografía real (desacoplado, ~6-10 min), revisar observaciones, avisos y costo, y compararlo con la línea de Felipe evaluada con el mismo modelo.
2. **Primera/última curva no cabe entre el punto inicial y el primer PI** (error típico «tangente 0: las curvas se pisan»): permitir empezar dentro de una curva o desplazar el PI.
3. **Rutas tortuosas**: añadir variantes con suavizado de la ruta (promedio móvil por abscisa a 50-300 m) antes de simplificar a PI, evaluadas por costo real como el resto de variantes.
4. **Curvas de volteo (horquillas)**: hoy una deflexión ≥ ~150° se declara inviable con mensaje claro; en terreno montañoso colombiano hacen falta (diseñar por arco, no por PI).
5. Calzadas separadas; secciones y corredor en Civil 3D (no viajan en el LandXML actual); más vehículos de diseño; tablas de talud por material; valores `[VERIFICAR]` contra el manual.
6. Probar `via.xml` / `superficie_terreno.xml` en Civil 3D (*Insertar → Importar LandXML*) y reportar mensajes.

## Lecciones de esta sesión (para no repetirlas)
- **Nunca reemplazar un bloque de pruebas por índices de texto** hasta el marcador final: borró 13 pruebas sin avisar (se vio porque el recuento bajó de 64 a 52). Comprobar siempre el número de pruebas.
- Cada arreglo lleva su prueba y una **prueba de mutaciones** (romper el código a propósito y comprobar que alguna prueba falla). Las mutaciones «equivalentes» (defensa en profundidad) se documentan.
- Antes de afirmar generalidad, barrido de casos distintos (terrenos empinados, huecos, ruido, extensión) y datos reales: encontró 6 defectos que las pruebas sintéticas no mostraban.


---

# Diseño vial automático — estado del proyecto

**Versión:** 0.8 (+ topografía real desde DWG, zonas sin datos por distancia, holgura de pendiente) · **Normas:** INVÍAS (Manual de Diseño Geométrico) + AASHTO Green Book 2011
**Guía:** diapositivas del curso de Diseño de Vías (EIA) — línea de ceros, velocidad de diseño, radios, espirales, entretangencias, vertical, sección transversal.

## Cómo ejecutarlo
```
pip install numpy scipy pandas openpyxl matplotlib ezdxf        # numpy >= 2 (usa np.trapezoid)
python -m diseno_vial.normas        # verifica las normas contra los ejemplos de las diapositivas
python -m diseno_vial.pruebas       # batería completa (59 pruebas, ~85 s, sin archivos externos; `-k texto` ejecuta solo las que contienen ese texto)
python -m diseno_vial.main PUNTOS.xlsx --hoja Topografia --categoria primaria_1_calzada --vtr 70 \
       --origen E,N --destino E,N --equidistancia 5 --salida salida/ \
       [--material roca_sana|roca_fracturada|suelo_comun|suelo_blando] [--talud-terraplen 1.5] [--e-pav 0.55] [--factor-terraplen 1.25] [--sin-optimizar-rasante] [--sin-replanteo] [--sin-optimizar-trazado]
python -m diseno_vial.validar salida/via.xml
```
Entrada: nube de puntos (X,Y,Z / Este,Norte,Cota; xlsx, csv, txt, las). Salida: `via.xml` (LandXML 1.2: superficie + alineamiento + perfil del terreno + rasante), `via.dxf` (planta con banca real y chaflanes), `perfil.dxf`, `secciones.dxf` (hoja de secciones cada 20 m), `planta.png`, `perfil.png`, `secciones.png`, `curvas.csv`, `pvi.csv`, `secciones.csv`, `resumen_movimiento_tierras.csv`, `optimizacion_rasante.csv/.png` (alternativas evaluadas y comparación base vs optimizada), `replanteo_rutas.csv` (rutas candidatas y la elegida), `optimizacion_trazado.csv` (historial de la búsqueda de PI y radios). En el DXF, la capa `RUTA_ELEGIDA` es la ruta del replanteo, `LINEA_CEROS` la del compás y `EJE_ANTES_DE_MOVER_PI` el eje previo a la optimización de PI y radios.
Civil 3D: *Insertar → Importar LandXML* (superficie, alineamiento con espirales, perfil del terreno y rasante con curvas verticales). El DXF se abre directo y se guarda como DWG.

## Hecho (v0.1)
| Módulo | Qué hace |
|---|---|
| `normas.py` | Tablas 4.1 (pend. media máx.), 5.2 (calzada), 5.4 (bermas) leídas de las diapositivas; Rmín AASHTO (reproduce 168 m a 70 km/h y 229 m a 80 km/h); peralte, espiral mínima, VECH (5 casos), entretangencias, K vertical (reproduce tabla 4.4), longitud virtual de Bruce |
| `topografia.py` | Nube → TIN (Delaunay) → curvas de nivel exactas; clasificación automática del terreno por pendiente predominante |
| `linea_ceros.py` | Método del compás (d = equidistancia/pendiente) con búsqueda A*: menor longitud, menores deflexiones, penaliza contrapendiente |
| `alineamiento.py` | Línea de ceros → tangentes (Douglas-Peucker + ajuste) → PI → radio que mejor ajusta a la línea respetando Rmín, espacio de tangentes, entretangencia y Lc ≥ 20 m → espirales (clotoide) |
| `exportar.py` | LandXML (TIN + Line/Spiral/Curve + perfil del terreno) y DXF (curvas, línea de ceros, eje, bordes, abscisado, PI, TE/SC/CS/ST) |
| `perfil.py` | Rasante: LP de mínimo corte/relleno con pendiente máx. (tabla 4.1) y mín. (0,5 %), PVI por simplificación, cotas de PVI por LP, curvas verticales (seguridad por Dp + operación 0,6·V, redondeo a 10 m), VECV según curvas horizontales, tangente mínima (tabla 4.3), avisos de longitud crítica |
| `secciones.py` | Sección tipo (calzada, bermas, bombeo, cuneta 3:1, taludes de corte/terraplén, espesor de pavimento), **peralte** con transición en 3 etapas a lo largo de la espiral, **sobreancho** (fórmula de las diapositivas, por el lado interno), chaflanes sobre el TIN, áreas por integración, volúmenes por áreas extremas, **diagrama de masas** con factor de terraplén compacto (acarreo longitudinal, DMA, compensación transversal, excedente/faltante), costos y avisos (terraplén sobre ladera > 20 %, visibilidad en curva) |
| `replanteo.py` + `linea_ceros.ruta_costo_minimo` | **Replanteo en planta según la sección**: tabla de costo transversal por metro (misma plantilla que las secciones: calzada, bermas, cuneta, taludes) → ruta de costo mínimo (Dijkstra, 16 direcciones, pendiente máx., solo nodos dentro del levantamiento) con distinto peso α del movimiento de tierras; cada ruta se diseña completa con **todas las variantes de ajuste de PI** (eligiendo por costo real, no por fidelidad geométrica) y se elige la de menos observaciones normativas y menor costo |
| `trazado.py` | **Optimización local de PI y radios por costo real**: parte del alineamiento válido, mueve cada PI ±paso en x e y (pasos decrecientes de 2 % a 0,25 % de la extensión) y cambia el radio de cada curva (×1,25 y ÷1,25); cada configuración se evalúa con el diseño completo (curvas con espirales, rasante, secciones, volúmenes, masas). Criterio lexicográfico (menos observaciones, luego menor costo). Rechaza curvas solapadas, PI fuera del levantamiento y chaflanes con cota extrapolada. BP y EP no se mueven. Presupuesto de evaluaciones configurable (`max_evals`) |
| `optimizar.py` | **Rasante optimizada por costo real**: (1) áreas reales de corte/relleno para cada cota candidata en cada estación; (2) programación dinámica exacta sobre una malla de cotas absoluta (restricción de pendiente) con multiplicador λ para el balance de masas → límite teórico y rasantes objetivo; (3) ajuste a PVI y curvas verticales normativas; (4) **pulido por búsqueda local** de cotas y abscisas de PVI sobre el costo total real (corte, relleno, acarreo, botadero, préstamo) manteniendo pendientes, tangentes mínimas y curvas. Elige la mejor alternativa, incluida la base: nunca queda peor |
| `validar.py` | Comprueba continuidad, abscisas y geometría del LandXML (eje y perfil) |
| `pruebas.py` | 59 pruebas (+ óptimo de la programación dinámica contra fuerza bruta, monotonía en λ, regresión del defecto original, factibilidad del pulido): normas, parábola vertical, ejemplos resueltos de las diapositivas (sobreancho, flecha M), áreas contra fórmulas cerradas, diagrama de masas calculado a mano, chaflán sobre TIN y talud, dirección física del peralte, sensibilidad a ancho/taludes, errores esperados, integración y regresiones |

Prueba base (terreno sintético ×10, primaria 1 calzada, 70 km/h): línea de ceros 2206 m, eje 2193 m, 2 curvas, 0 observaciones; perfil con 6 PVI y 4 curvas verticales, pendiente máx. 7 %, corte máx. 3,2 m, relleno máx. 1,8 m. LandXML válido; DXF sin errores de auditoría.

Secciones del caso base (calzada 7,30 m, bermas 1,80 m, corte 1H:1V, terraplén 1,5H:1V, f = 1,25): corte ≈ 29.100 m³, relleno compacto ≈ 6.400 m³, afectación media 15 m, sobreancho máx. 0,94 m, acarreo longitudinal ≈ 6.400 m³ con DMA ≈ 350 m, excedente ≈ 21.100 m³.

Replanteo del caso base: 6 rutas candidatas (compás + 5 de costo mínimo) con 2-3 variantes de PI cada una; la elegida cuesta **8,9 % menos** que la línea de ceros por compás con la misma rasante base (140.779 → 128.317). Combinado con la optimización de rasante: 140.779 → **113.279 (−19,5 %)**, excedente 21.085 → 368 m³. Tarda ~40 s en total.

Cadena completa del caso base (replanteo → mover PI y radios → rasante optimizada), ~36 s: **140.779 → 128.317 → 112.443 → 107.171 (−23,9 %)**. El eje pasa de 2.193 a 2.270 m (más largo) pero el corte baja de 29.088 a 7.638 m³ y el excedente a 387 m³; 0 observaciones normativas.

Optimización de rasante del caso base (misma sección): costo total 140.779 → 127.802 (**−9,2 %**), corte 29.088 → 18.182 m³, relleno compacto 6.403 → 14.638 m³, excedente 21.085 → 0 m³, acarreo longitudinal 6.438 → 2.686 m³. En el terreno ondulado de las pruebas: 146.046 → 135.955 (−6,9 %). Tarda unos 15-20 s; `--sin-optimizar-rasante` la omite.

### Defectos encontrados y corregidos (v0.2 a v0.7)
1. **Curvas solapadas pasaban en silencio**: en un terreno ondulado el eje medía 4777 m (línea de ceros: 2456 m) y el perfil se calculaba sobre él. Ahora se detectan los solapes (`al.valido`), el perfil se niega a diseñarse sobre un eje inválido y `disenar_via` lanza un error claro.
2. **Una recta «cumplía» las normas por vacío**: al escalar la tolerancia de PI, la solución con 0 curvas ganaba aunque se alejara 135 m de la línea de ceros. Ahora solo cuentan diseños a menos de max(50 m, 3 % de la extensión) de la línea de ceros; si ninguno cumple, se avisa.
3. Tablas 5.2 y 5.4 de las diapositivas se contradicen en una celda (terciaria-escarpado, 40 km/h): se conserva la transcripción y queda documentado en las pruebas.
4. **Prueba de mutaciones (v0.3):** se rompió el código a propósito en 8 puntos; 7 los detectaban las pruebas y 1 no (peralte volcado hacia el lado equivocado: solo se comprobaba su magnitud). Se reforzó la prueba con la geometría física (borde externo sube, interno baja, carril interno más ancho) y ahora se detectan las variantes a izquierda y derecha.
6. **Optimizador, primer intento inválido (v0.4):** la programación dinámica con malla de cotas *relativa al terreno* (paso 0,25 m) costaba el doble que la rasante base: el incremento máximo alcanzable por estación quedaba hasta 0,25 m por debajo del límite de pendiente, así que en un tramo sostenido al 7 % la rasante se hundía hasta el borde de la malla. Se detectó porque un óptimo exacto no puede costar más que una solución factible. Corrección: malla **absoluta** de 0,05 m con un nivel de holgura (el ajuste a PVI restituye el límite exacto). Además, el rango de λ permitía excavar «gratis» (349.000 m³): ahora λ ∈ (−c_corte, c_relleno/f).
7. Prueba de mutaciones del optimizador: 6 roturas intencionales; 5 se detectaban y 1 no (el pulido ignorando conflictos de curvas). Un primer intento de prueba tampoco la cubría (sus casos «inválidos» violaban en realidad la pendiente máxima, que se rechaza antes). Se añadió una prueba con pendientes legales y curvas en conflicto; ahora se detectan las 6.
8. **Replanteo, rutas pegadas al borde (v0.5):** la ruta de costo mínimo se pegaba al borde del terreno aprovechando cotas extrapoladas (sin datos), con pendientes falsas de hasta 11,5 % frente al 7 % permitido. Ahora solo se usan nodos dentro del casco del TIN con margen de 3 celdas, un origen/destino pegado al borde se rechaza con mensaje claro, y las secciones avisan si un chaflán cae fuera del área levantada.
9. **Ajuste de PI con deflexión 0 (v0.5):** con rutas raster la limpieza de quiebres casi rectos se detenía a las 40 iteraciones y el diseño dividía por una deflexión cero (`ZeroDivisionError`). Ahora se eliminan todos los quiebres colineales por pasada y hay un filtro defensivo. El error estaba oculto: un `except Exception` lo convertía en «ruta no válida»; ahora solo se toleran los errores de inviabilidad esperados.
10. Mutaciones del replanteo y los arreglos 8 y 9: 6 roturas intencionales, todas detectadas (una primera tanda excedió los 300 s al correr la batería completa en cada mutación; se añadió el filtro `-k`).
11. **Prueba de configuraciones inviables mal diseñada (v0.6):** los dos casos «inviables» producían solapes enormes (3.339 m y 71.626 m), así que solo los frenaba la validez del alineamiento y nunca se ejercitaban las demás defensas. Se rediseñaron: PI fuera del terreno con curvas válidas y solape moderado (148 m) dentro del terreno.
12. **Mutaciones equivalentes por defensa en profundidad:** quitar una sola comprobación (casco de los PI, validez del alineamiento) no cambia el comportamiento porque otra capa rechaza lo mismo (cobertura de chaflanes; perfil y secciones). No son detectables por ninguna prueba de conducta; las mutaciones combinadas (quitar todas las capas redundantes a la vez) sí se detectan. Las demás mutaciones del optimizador de trazado (aceptar peores, ignorar observaciones, no aplicar cambios de radio) también.
13. **Entorno reiniciado** durante la sesión: se perdieron `ezdxf` y los datos de ejemplo de `/home/claude`; se reinstalaron dependencias y se recuperó el terreno desde el zip original. El proyecto en `/mnt/user-data/outputs` no se vio afectado. Toda la batería se volvió a pasar tras el reinicio.
14. **Curvas de retorno (horquillas) agotaban la memoria (v0.7):** con una deflexión cercana a 180° la tangente `T = R·tan(Δ/2)` llega a ~10¹⁷ m; el diseño intentaba muestrear esa geometría (`MemoryError`). Ahora se declara inviable con un mensaje claro («curva de retorno, requiere una curva de volteo»), `muestrear` rechaza elementos absurdos y el replanteo descarta esa candidata.
15. **PI duplicados y NaN silenciosos en `ajustar_pis`:** en terreno empinado 7 de 48 llamadas devolvían PI consecutivos idénticos; el NaN resultante solo se «arreglaba» porque un filtro posterior los descartaba por casualidad. Ahora se eliminan de raíz, `_unit` y `_deflexiones` no pueden dar NaN, el filtro de colineales descarta repetidos y `disenar_desde_pis` rechaza PI coincidentes. La prueba de regresión usa la ruta real capturada (337 puntos), porque ninguna entrada sintética (1.600 pruebas) lo reproducía.
16. **`0·∞ = NaN` en la ruta de costo mínimo con α = 0:** el resultado era correcto por casualidad; ahora «arista no construible» es una exclusión explícita para cualquier α.
17. **Diseños «válidos» con costo absurdo pasaban sin aviso:** ahora se avisa si corte + relleno superan 5 veces el costo de la vía (casos con pendientes naturales de 20 % o más).
18. **Zonas sin datos ignoradas:** el eje podía cruzar un hueco de 300 m del levantamiento (21 % de su longitud) diseñando sobre terreno interpolado sin avisar. Ahora se detectan por triángulos del TIN de más de 8 veces la arista mediana y se avisa con los metros afectados (subestima la extensión de huecos pequeños; no bloquea rutas).
19. Nota sobre el script antiguo `via_linea_de_ceros.py`: su costo para elegir la ruta omite el área bajo la plataforma en terreno inclinado (s·b²/2 por lado). No afecta a este paquete (la evaluación final usa secciones reales), pero sus costos de búsqueda están subestimados.

## ¿Es general o específica de un caso? (barrido de 13 casos, v0.7)
Arquitectura: general (nube de puntos arbitraria, origen/destino, categoría, VTR, sección y costos como parámetros; sin coordenadas fijas en el paquete). **Validación: limitada a terrenos sintéticos de 2-6 km**, nunca con topografía real.

| Caso | Resultado |
|---|---|
| Desarrollo (2 km, primaria 70), sentido inverso, otra diagonal | funciona; invertir el sentido cambia el costo ≈ 3 % (heurísticas dependientes del sentido) |
| Mismo terreno como secundaria 50 y terciaria 30 | funciona |
| Relieve ×2,5 (≈ 20 %), primaria 70 / secundaria 40 montañosa | funciona; costos ≈ 10× (aviso de movimiento de tierras dominante) |
| Relieve ×5 (≈ 40 %), terciaria 30 escarpada | **inviable**: «ninguna ruta candidata produjo un alineamiento válido» (necesita curvas de volteo) |
| Extensión de 6 km | funciona; 31 s sin optimizaciones (la resolución de la búsqueda de rutas es extensión/250) |
| Ruido de 0,4 m + hueco de 300×300 m | funciona; avisa de zonas sin datos |
| Colinas aisladas de 60 y 150 m | funciona; la ruta las rodea |
| Ladera uniforme de 20 % a subir | «funciona» con costo de 10,7 millones (aviso) |
| Terreno casi llano (1 %) | funciona (VTR 80, rango de la tabla 5.2) |

Antes de este barrido fallaban dos casos reales (puntos 14 y 15 de arriba). **No** se ha probado con datos reales (LiDAR de millones de puntos, ríos, vías existentes, discontinuidades), con extensiones de decenas de km ni con el sistema de coordenadas MAGNA-SIRGAS.

## Pendiente (orden sugerido)
0. **Curvas de volteo (horquillas)**: imprescindibles en terreno montañoso y escarpado colombiano; hoy el programa solo dice «inviable». Requiere diseñar por arco (centro y radio) en lugar de por PI.
0b. **Prueba con topografía real** (un levantamiento o LiDAR de tu proyecto): es lo que realmente demostraría la generalidad. Medir rendimiento con millones de puntos (TIN, malla de búsqueda).
1. **Comparador de escenarios** (ancho, material de corte, VTR) en lote, con tabla y gráfica de costo, volúmenes, desplazamiento de la ruta y del eje. Ya se comprobó que la ruta reacciona al ancho cuando el terreno ofrece alternativas.
2. **Búsqueda conjunta** de trazado y rasante: hoy los PI/radios se optimizan con la rasante base y la rasante se optimiza después; el eje óptimo para la rasante final puede ser otro. Valorar también añadir/quitar PI y tener en cuenta la visibilidad en curva (hoy solo avisos).
3. **Calzadas separadas** (primaria de dos calzadas): hoy se rechaza con un error claro.
4. **Civil 3D**: ensamblaje/corredor y secciones no viajan en el LandXML actual. Opciones: `CrossSects` en LandXML (formato por validar), script (.NET/LISP/Dynamo) o superficie de diseño.
5. Más vehículos de diseño para el sobreancho; tablas de talud por material y altura del manual; transición de peralte en R entre 1000 y 3500 m; criterio de drenaje de curvas verticales; parámetros por JSON/interfaz; curvas espiral-espiral; inserción/eliminación de PVI en el pulido; huecos en el levantamiento.
6. Prueba en Civil 3D por el usuario y ajustes según lo que reporte.

## Limitaciones conocidas
- **No probado dentro de Civil 3D** (no tengo acceso): validar importación y reportar cualquier mensaje.
- El terreno de ejemplo es sintético; coordenadas locales (sin MAGNA-SIRGAS).
- Si el destino cae en una zona más empinada que la pendiente máxima, la línea hace un gancho (se vio con destino 1900,1850): falta restricción dura de deflexión.
- Con tangentes > 600 m, la regla de VECH de las diapositivas sube la velocidad específica (VTR+20), lo que agranda los radios.
- No hay curvas espiral-espiral; si no cabe una curva se reporta como observación normativa.
- En trazados con giros cerrados y poca distancia entre ellos (p. ej. 26° y 64° separados 315 m a 80–90 km/h) el programa entrega el mejor eje válido **con observaciones** (entretangencia insuficiente) en vez de ocultarlas.
- Longitud crítica de pendiente: leída de forma aproximada de la gráfica 4.2 (±15 %); solo genera avisos.
- La pendiente máx. de tangente se toma de la tabla 4.1 (la misma que la media del corredor, como en las diapositivas).
- Generalidad: ver la sección «¿Es general o específica de un caso?». Los umbrales dependen de la extensión (tolerancia de PI = extensión/80, malla de rutas = extensión/250, celda del TIN = extensión/100): en vías de decenas de km la resolución puede ser insuficiente. La detección de zonas sin datos usa 8× la arista mediana (no detecta huecos menores de ~230 m con celdas de 20 m).
- Mover PI y radios: búsqueda local (puede quedar en un óptimo local), con presupuesto de ~120 evaluaciones (≈ 40 s en el caso base). El eje final puede alejarse de la línea de ceros 30-45 m (el costo real manda sobre la fidelidad geométrica). El criterio no incluye los avisos de visibilidad en curva: el caso base conserva 2 avisos (talud de corte, corrida ≈ 2,4-2,7 m).
- Replanteo: el costo no es monótono en α (rutas parecidas dan costos distintos por cómo cierran las curvas); se elige la mejor de las candidatas evaluadas, que es el costo real más bajo encontrado, no necesariamente el óptimo global. El modelo de costo por metro asume terreno plano inclinado en cada arista y eje a nivel del terreno (sin balance de masas ni acarreo); la evaluación final sí usa secciones reales con acarreo.
- Optimización de rasante: el límite teórico de la programación dinámica (sin curvas verticales normativas) es mucho más bajo que lo alcanzable con la geometría real (en el caso base ≈ 10.000 frente a ≈ 40.000 de costo de tierras): casi todo lo perdido se debe a curvas verticales largas (V = 90 → L ≈ 300 m) y tangentes mínimas de 250 m. La mejora real (6-9 %) depende de los costos unitarios (botadero 0,5, préstamo 2,0, acarreo 0,002 por m³·m son placeholders): con otros precios el óptimo cambia.
- El pulido es búsqueda local (puede quedar en un óptimo local); evalúa unas cientos de combinaciones de cotas/abscisas de PVI.
- Sección transversal: los taludes por tipo de material (`MATERIALES_CORTE`), la pendiente de berma (4 %), la cuneta (0,80 m, 3:1), el espesor de pavimento (0 por defecto) y los costos unitarios son **valores típicos/placeholder**, no del manual: ajustarlos a tu proyecto. El factor de terraplén compacto (1,25) es el del taller N°5.
- El criterio de visibilidad en curva (M requerida vs. disponible) es simplificado (línea de visión a 0,84 m sobre el pavimento); la fórmula de M reproduce el ejemplo de las diapositivas, pero la comprobación con el talud es aproximada.
- Transiciones de peralte y sobreancho solo en curvas espiralizadas (R ≤ 1000 m). El volumen se calcula por áreas extremas cada 10 m (sin corrección por curvatura).

## Valores a VERIFICAR contra el manual (marcados `[VERIFICAR]` en el código)
Límites de clasificación de terreno por pendiente transversal · tipos de talud por material · pendiente de berma · costos unitarios · rampa relativa de peralte Δ por velocidad · criterio de comodidad J=0,6 en la espiral · peralte vs. radio (aquí e = e_max·Rmín/R; el manual usa el Método 5 de AASHTO) · longitud mínima de curvas con deflexión < 6° · e_max = 8 %.
