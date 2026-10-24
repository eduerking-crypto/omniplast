"""Modelo de dispersión de microplásticos (advección + difusión + sedimentación Stokes).

Modelo lagrangiano simple de partículas para humedales costeros: viento, corriente
uniforme o gradiente, caminata aleatoria turbulenta y velocidad de hundimiento por
ley de Stokes corregida (Kowalski et al. 2016; Waldschläger & Schüttrumpf 2019).
"""
import numpy as np
import pandas as pd


def stokes_settling_velocity(diameter_mm, polymer_density, water_density=1000,
                             dynamic_viscosity=0.001, g=9.81):
    """Velocidad terminal (m/s) con correccion empirica para particulas mayores."""
    d = diameter_mm / 1000
    dp = polymer_density - water_density
    if dp <= 0:
        return 0.0
    v_stokes = g * d ** 2 * dp / (18 * dynamic_viscosity)
    reynolds = water_density * abs(v_stokes) * d / dynamic_viscosity
    if reynolds < 0.1:
        return float(v_stokes)
    drag = 1 + 0.15 * reynolds ** 0.687
    return float(v_stokes / drag)


def simulate_drift(origin_lat, origin_lon, hours=24, dt=600, n_particles=1000,
                   wind_speed_ms=3.0, wind_direction_deg=45, wind_factor=0.02,
                   current_speed_ms=0.15, current_direction_deg=90,
                   diffusion_m2s=0.05, diameter_mm=0.5, polymer_density=960,
                   depth_m=2.0, bbox=None, seed=42):
    """Simulacion Monte Carlo de trayectorias; devuelve posiciones, densidad y destinos."""
    rng = np.random.default_rng(seed)
    n_steps = int(hours * 3600 / dt)
    wind_dir = np.radians(wind_direction_deg)
    curr_dir = np.radians(current_direction_deg)
    vx = current_speed_ms * np.cos(curr_dir) + wind_factor * wind_speed_ms * np.cos(wind_dir)
    vy = current_speed_ms * np.sin(curr_dir) + wind_factor * wind_speed_ms * np.sin(wind_dir)
    sink = stokes_settling_velocity(diameter_mm, polymer_density)
    m_per_deg_lat = 110574.0
    m_per_deg_lon = 111320.0 * np.cos(np.radians(origin_lat))
    lat = np.full(n_particles, origin_lat, float)
    lon = np.full(n_particles, origin_lon, float)
    z = rng.uniform(0, depth_m, n_particles)
    states = np.array(["flotando"] * n_particles, dtype=object)
    history = []
    for step in range(n_steps):
        active = states == "flotando"
        if not active.any():
            break
        if step % 3 == 0:
            history.append({"paso": step, "min": step * dt / 60,
                            "lat_media": float(lat[active].mean()),
                            "lon_media": float(lon[active].mean()),
                            "activos": int(active.sum())})
        sigma = np.sqrt(2 * diffusion_m2s * dt)
        d_lat_m = rng.normal(0, sigma, n_particles) + vy * dt
        d_lon_m = rng.normal(0, sigma, n_particles) + vx * dt
        lat[active] += d_lat_m[active] / m_per_deg_lat
        lon[active] += d_lon_m[active] / m_per_deg_lon
        z[active] -= sink * dt
        beaching = active & np.zeros(n_particles, bool)
        if bbox:
            out = (lon < bbox["lon_min"]) | (lon > bbox["lon_max"]) | \
                  (lat < bbox["lat_min"]) | (lat > bbox["lat_max"])
            beaching = active & out
        settled = active & (z <= 0) & (sink > 0)
        states[beaching] = "varado_en_costa"
        states[settled] = "sedimentado"
        z[settled] = 0
    counts = pd.Series(states).value_counts()
    displacement = np.sqrt((lat - origin_lat) * m_per_deg_lat) ** 0
    disp_m = np.hypot((lat - origin_lat) * m_per_deg_lat, (lon - origin_lon) * m_per_deg_lon)
    return {"test": "Modelo de dispersión lagrangiano (Monte Carlo)",
            "posiciones_finales": pd.DataFrame({"lat": lat, "lon": lon, "estado": states,
                                                "z_m": z}),
            "destinos": counts.to_dict(),
            "historia": pd.DataFrame(history),
            "desplazamiento_m": {"media": float(disp_m.mean()),
                                 "p95": float(np.percentile(disp_m, 95)),
                                 "max": float(disp_m.max())},
            "parametros": {"horas": hours, "n": n_particles, "dt_s": dt,
                           "v_ms": (wind_speed_ms, wind_direction_deg),
                           "c_ms": (current_speed_ms, current_direction_deg),
                           "hundimiento_m_s": sink, "densidad_polimero": polymer_density},
            "interpretacion": (f"Tras {hours} h: {counts.get('flotando', 0)} flotando, "
                               f"{counts.get('sedimentado', 0)} sedimentadas, "
                               f"{counts.get('varado_en_costa', 0)} varadas; "
                               f"desplazamiento medio {disp_m.mean():.0f} m. "
                               f"Vstokes={sink*1e6:.1f} µm/s (Ø {diameter_mm} mm, "
                               f"ρ {polymer_density} kg/m³)")}
