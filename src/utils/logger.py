import logging

LOGGER_NAME = 'app_logger'


def create_logger():
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    file_handler = logging.FileHandler('app.log')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.DEBUG)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
