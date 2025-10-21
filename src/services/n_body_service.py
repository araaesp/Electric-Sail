from n_body.n_body_solver import NBodySolver, CollisionException
from utils.converter import kepler_to_cartesian, probe_kepler_to_cartesian
from model.electric_probe import ElectricSailProbe
import pandas as pd
import json


def solve(data):
    bodies, time_span = process_data(data)

    solver = NBodySolver(bodies)
    try:
        y = solver.solve(time_span=time_span, max_step=time_span / 1000)
    except CollisionException as exception:
        y = exception.simulation
        print(exception.message)

    y = format_result(y)
    df = pd.DataFrame(y, columns=create_cols_name(y))

    return json.loads(df.to_json(orient='split', index=False))


def process_data(data):
    corpos_finais = data.get('bodies', [])

    if 'keplerian_bodies' in data:
        corpos_finais.extend(data['keplerian_bodies'])
    
    if 'solar_probes' in data:
        corpos_finais.extend(data['solar_probes'])

    if 'keplerian_solar_probes' in data:
        corpos_finais.extend(data['keplerian_solar_probes'])

    if 'electric_sail_probes' in data:
        corpos_finais.extend(data['electric_sail_probes'])

    start_time_str = data.get('start_time')
    end_time_str = data.get('end_time')
    
    from datetime import datetime
    start_time = datetime.strptime(start_time_str, "%Y-%m-%d %H:%M:%S")
    end_time = datetime.strptime(end_time_str, "%Y-%m-%d %H:%M:%S")
    time_span = (end_time - start_time).total_seconds()
    
    return corpos_finais, time_span

def format_result(y_t):
    cols = int((y_t.shape[1] - 1) / 2)
    file_cols = [i for i in range(cols + 1)]

    # if remove_z_axis:
    #     file_cols = [i for i in range(cols) if i % 3 != 0 or i == 0]

    return y_t[:, file_cols]


def create_cols_name(data):
    number_of_bodies = int(data.shape[1] / 3)
    cols = ['T']
    for body_number in range(number_of_bodies):
        for letter in ['X{}', 'Y{}', 'Z{}']:
            cols.append(letter.format(body_number))

    return cols
