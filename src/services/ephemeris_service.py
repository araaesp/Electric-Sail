from skyfield.api import Topos, load
import json
from pytz import timezone
import datetime

ts = load.timescale()
planets = load('de421.bsp')


def get_body_data(body_name, json_date):
    date = datetime.datetime.strptime(json_date, '%Y-%m-%dT%H:%M:%S.%f')
    time = ts.utc(timezone('UTC').localize(date))

    with open('../resources/solar-system.json') as json_file:
        data = json.load(json_file)
        astra = data[body_name]

    ephemeris = planets[body_name].at(time)

    astra['position'] = ephemeris.position.km.tolist()
    astra['velocity'] = ephemeris.velocity.km_per_s.tolist()

    return astra


def get_available_bodies():
    with open('../resources/solar-system.json') as json_file:
        data = json.load(json_file)
        return list(data.keys())
