# -*- coding: utf-8 -*-
import osmnx as ox
from osmnx.projection import project_gdf
import networkx as nx
from itertools import islice
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import LineString
from shapely.affinity import translate
import matplotlib.lines as mlines
from geopy.geocoders import Nominatim


def obtener_coordenadas(nombre_lugar):
    geolocalizador = Nominatim(user_agent="ruta_saludable")
    ubicacion = geolocalizador.geocode(nombre_lugar)
    if ubicacion:
        return (ubicacion.latitude, ubicacion.longitude)
    else:
        raise ValueError(f"No se encontró el lugar: {nombre_lugar}")



# ============ 1) Grafo peatonal (área algo más grande) ============
ox.settings.use_cache = True
ox.settings.log_console = False

# Centro aproximado Miraflores
lat, lon = -12.1209, -77.0290
G = ox.graph_from_point((lat, lon), dist=1500, network_type='walk')  # 1.5 km
G_simple = nx.Graph(G)

# ============ 2) Puntos seleccionados ============
# Opción F:
"""
origen_nombre = input("¿Desde dónde quieres salir? (ej. Parque Kennedy): ")
destino_nombre = input ("¿A dónde quieres ir? (ej. Larcomar): ")


origen = obtener_coordenadas(origen_nombre)
destino = obtener_coordenadas(destino_nombre)
"""


def pedir_modo_entrada():
    while True:
        modo = input("¿Cómo quieres ingresar los lugares? Escribe 'nombre' o 'coordenadas': ").strip().lower()
        if modo in ['nombre', 'coordenadas']:
            return modo
        print("⚠️ Opción no válida. Escribe 'nombre' o 'coordenadas'.")

modo_entrada = pedir_modo_entrada()

def pedir_coordenadas(mensaje):
    while True:
        nombre = input(mensaje)
        try:
            return obtener_coordenadas(nombre)
        except ValueError as e:
            print(f"⚠️ {e}")
            print("Por favor, intenta con otro nombre más específico o revisa la ortografía.")

# Entrada robusta
"""
origen = pedir_coordenadas("¿Desde dónde quieres salir? (ej. Parque Kennedy): ")
destino = pedir_coordenadas("¿A dónde quieres ir? (ej. Larcomar): ")
"""
if modo_entrada == 'nombre':
    def pedir_coordenadas_por_nombre(mensaje):
        while True:
            nombre = input(mensaje)
            try:
                return obtener_coordenadas(nombre)
            except ValueError as e:
                print(f"⚠️ {e}")
                print("Intenta con otro nombre más específico.")
    origen = pedir_coordenadas_por_nombre("¿Desde dónde quieres salir? (ej. Parque Kennedy): ")
    destino = pedir_coordenadas_por_nombre("¿A dónde quieres ir? (ej. Larcomar): ")

else:
    def pedir_coordenadas_directas(mensaje):
        while True:
            entrada = input(mensaje)
            try:
                lat, lon = map(float, entrada.split(","))
                return (lat, lon)
            except:
                print("⚠️ Formato incorrecto. Usa 'latitud,longitud' (ej. -12.1209,-77.0290)")
    origen = pedir_coordenadas_directas("Ingresa las coordenadas de origen (latitud,longitud): ")
    destino = pedir_coordenadas_directas("Ingresa las coordenadas de destino (latitud,longitud): ")



# Nodos más cercanos
nO = ox.distance.nearest_nodes(G, origen[1],  origen[0])
nD = ox.distance.nearest_nodes(G, destino[1], destino[0])

# ============ 3) Funciones auxiliares ============
def ruta_dist_m(Gmulti, ruta):
    """Suma la longitud (m) de una ruta en un MultiDiGraph."""
    d = 0.0
    for u, v in zip(ruta[:-1], ruta[1:]):
        data = Gmulti.get_edge_data(u, v)
        if not data:
            continue
        attrs = next(iter(data.values()))
        d += float(attrs.get('length', 0.0))
    return d

def jaccard_edges(Gmulti, r1, r2):
    """Similitud Jaccard entre conjuntos de aristas de dos rutas."""
    e1 = {(u, v) for u, v in zip(r1[:-1], r1[1:])}
    e2 = {(u, v) for u, v in zip(r2[:-1], r2[1:])}
    inter = len(e1 & e2)
    union = len(e1 | e2)
    return inter / union if union else 1.0

# ============ 4) Generar k rutas distintas ============
candidatas = islice(nx.shortest_simple_paths(G_simple, nO, nD, weight='length'), 12)

rutas = []
for r in candidatas:
    if all(jaccard_edges(G, r, r_exist) < 0.85 for r_exist in rutas):
        rutas.append(r)
    if len(rutas) == 3:
        break

if not rutas:
    rutas = [nx.shortest_path(G_simple, nO, nD, weight='length')]

print("\nVERIFICACIÓN DE RUTAS GENERADAS:")
for i, r in enumerate(rutas):
    print(f"Ruta {i+1}: {len(r)} nodos | distancia ≈ {ruta_dist_m(G, r):.1f} m")



# Distancias y tiempos
dist_m = [ruta_dist_m(G, r) for r in rutas]
v_ms = 1.4  # velocidad promedio al caminar (m/s)
tiempos_min = [d / v_ms / 60 for d in dist_m]

# ============ 5) Parques cercanos ============
tags = {"leisure": ["park", "garden"]}
parques_wgs = ox.features_from_point((lat, lon), dist=1500, tags=tags)

def ruta_linestring_wgs(Gmulti, r):
    coords = [(Gmulti.nodes[n]['x'], Gmulti.nodes[n]['y']) for n in r]
    return LineString(coords)

rutas_wgs = [ruta_linestring_wgs(G, r) for r in rutas]

# Proyectar para medir distancias (m)
parques_proj = project_gdf(parques_wgs)
rutas_proj = gpd.GeoSeries(rutas_wgs, crs="EPSG:4326").to_frame("geometry")
rutas_proj = project_gdf(rutas_proj)


# Desplazar visualmente cada ruta para diferenciarlas
offset_px = 1  # ajusta este valor si quieres más separación
rutas_wgs_offset = [translate(r, xoff=offset_px*i, yoff=offset_px*i) for i, r in enumerate(rutas_wgs)]

# Buffer de sombra (40 m alrededor de parques)
buffer_m = 40
dist_min_park_m = []
pct_en_sombra = []

parques_union = parques_proj.union_all()
parques_buffer = parques_union.buffer(buffer_m)

for geom in rutas_proj.geometry:
    dmin = geom.distance(parques_union)
    dist_min_park_m.append(dmin)
    inter = geom.intersection(parques_buffer)
    pct = (inter.length / geom.length) if geom.length > 0 else 0
    pct_en_sombra.append(float(pct))

# Índices
idx_salud = max(range(len(rutas)), key=lambda i: (pct_en_sombra[i], -dist_min_park_m[i], -dist_m[i]))
idx_corta = min(range(len(rutas)), key=lambda i: dist_m[i])

# ============ 6) Visualización ============

colores = ['#E53935','#FDD835', '#1E88E5']  # rojo, amarillo, azul
labels  = ['Más corta', 'Ruta saludable', 'Otra']

fig, ax = ox.plot_graph(G, show=False, close=False, bgcolor='black', node_size=0, edge_color='#777777', edge_linewidth=0.6)

# Dibujar parques
if not parques_wgs.empty:
    parques_wgs.plot(ax=ax, facecolor='#4CAF50', edgecolor='#2E7D32', alpha=0.35)

# Dibujar rutas
#orden = [idx_corta] + [idx_salud] + [i for i in range(len(rutas)) if i not in {idx_corta, idx_salud}]

# Asegura que los índices sean únicos y en orden deseado
orden = []
if idx_corta != idx_salud:
    orden = [idx_corta, idx_salud]
else:
    orden = [idx_corta]
# Agrega los demás sin repetir
orden += [i for i in range(len(rutas)) if i not in orden]

for j, i in enumerate(orden):
    ox.plot_graph_route(G, rutas[i], ax=ax, route_linewidth=4, route_color=colores[j], orig_dest_size=30, show=False, close=False)

# Leyenda con colores reales
line_corta = mlines.Line2D([], [], color=colores[0], linewidth=4, label='Más corta')
line_salud = mlines.Line2D([], [], color=colores[1], linewidth=4, label='Ruta saludable')
line_otra  = mlines.Line2D([], [], color=colores[2], linewidth=4, label='Otra')
plt.legend(handles=[line_corta, line_salud, line_otra], loc='lower right', facecolor='white', edgecolor='black')



# ============ 7) Resumen de rutas ============
def fmt(m): return f"{m/1000:.2f} km"

print("—"*60)
print("RESUMEN DE RUTAS")
for i, (d, t, pct, dpark) in enumerate(zip(dist_m, tiempos_min, pct_en_sombra, dist_min_park_m), 1):
    etiqueta = ""
    if i-1 == idx_corta:
        etiqueta = "🟥 Más corta"
    elif i-1 == idx_salud:
        etiqueta = "🟨 Ruta saludable"
    else:
        etiqueta = "🟦 Otra"
    print(f"{etiqueta}: {fmt(d)} | {t:.1f} min | cerca/parques: {dpark:.0f} m | sombra≈ {pct*100:.0f}%")

delta_min = tiempos_min[idx_salud] - tiempos_min[idx_corta]
if idx_salud != idx_corta:
    if delta_min <= 3:
        msg_tiempo = f"solo ~{delta_min:.1f} min más que la más corta"
    else:
        msg_tiempo = f"~{delta_min:.1f} min más que la más corta"
    print("\nRECOMENDACIÓN:")
    print(f"✅ Toma la **Ruta saludable** (🟨): {fmt(dist_m[idx_salud])}, {tiempos_min[idx_salud]:.1f} min, {msg_tiempo}.")
    print("   Pasa junto a parques/jardines (más sombra), lo que **reduce tu exposición directa al sol** ☀️🌳.")
else:
    print("\nRECOMENDACIÓN:")
    print("✅ La ruta más corta también es la más sombreada cercana a parques. Buenas noticias para tu bienestar ☀️🌳")

plt.show()