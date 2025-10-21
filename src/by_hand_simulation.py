import json
import logging
import concurrent.futures
from datetime import datetime

import numpy as np

from services import n_body_service
from n_body.n_body_solver import NBodySolver, CollisionException
from utils.decorators import log_time
from utils.logger import create_logger, LOGGER_NAME
from utils.schemas import RequestDto

OUTPUT_FILE_PATH = '../simulation_results/{}-{}.txt'

JSON_FILE_PATH = '../resources/electric_sail_probe.json'
#JSON_FILE_PATH = '../resources/demonstration.json'

logger = logging.getLogger(LOGGER_NAME)

with open(JSON_FILE_PATH) as json_file:
    input_json = json.load(json_file)


@log_time
def write_on_file(y_t, theta, remove_z_axis=False):
    cols = int(y_t.shape[1])

    if remove_z_axis:
        file_cols = [i for i in range(cols) if i % 3 != 0 or i == 0]
    else:
        file_cols = [i for i in range(cols)]

    np.savetxt(OUTPUT_FILE_PATH.format(theta, datetime.now().strftime('%Y-%m-%d_%H-%M-%S.%f')), y_t[:, file_cols],
               delimiter=',')


@log_time
def run_synchronous_simulation():
    logger.info("Starting synchronous process")

    # thetas = range(0, 360, 30)

    # for true_anomaly in thetas:
    run_simulation(0)


@log_time
def run_parallel_simulation():
    logger.info("Starting process")

    alpha = range(-90, 91, 30)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(run_simulation, alpha)


def run_simulation(alpha):
    data = RequestDto().load(input_json)
    #keplerian_solar_probes = data.get('keplerian_solar_probes')
    #if 90 < alpha < 270:
    #    alpha = alpha - 180
    #keplerian_solar_probes[0].alpha = np.deg2rad(alpha)

    bodies, time_span = n_body_service.process_data(data)

    logger.info(
        f'Simulation with {len(bodies)} bodies. Time span = {time_span / 86400:.2f} days. Alpha = {alpha:.2f} degrees')

    solver = NBodySolver(bodies)

    try:
        y = solver.solve(time_span=time_span, max_step=86400 / 2)
    except CollisionException as exception:
        y = exception.simulation
        print(exception.message)

    write_on_file(y, alpha, False)

if __name__ == '__main__':
    create_logger()
    run_synchronous_simulation()
    #run_parallel_simulation()
    exit(0)