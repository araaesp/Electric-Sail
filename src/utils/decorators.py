import time
import logging

from functools import wraps
from utils.logger import LOGGER_NAME

logger = logging.getLogger(LOGGER_NAME)


def log_time(func):
    @wraps(func)
    def wrapper_decorator(*args, **kwargs):
        start_time = time.perf_counter()
        value = func(*args, **kwargs)
        end_time = time.perf_counter()
        run_time = end_time - start_time
        logger.info(f'{func.__name__!r} was executed in {run_time:.6f} seconds')
        return value

    return wrapper_decorator
