from flask import Flask, request
from flask_restful import Api, Resource, reqparse
from services import n_body_service, ephemeris_service
from utils.logger import create_logger
from utils.schemas import RequestDto

app = Flask(__name__)
api = Api(app)

parser = reqparse.RequestParser()
parser.add_argument('time')
parser.add_argument('bodies')
parser.add_argument('satellites')


class SolverController(Resource):
    def get(self):
        return 'welcome to my api'

    def post(self):
        body_json = request.get_json()
        data = RequestDto().load(body_json)

        y = n_body_service.solve(data)

        return y


class BodyInfo(Resource):
    def get(self, body_name, time):
        return ephemeris_service.get_body_data(body_name, time)


class Bodies(Resource):
    def get(self):
        return ephemeris_service.get_available_bodies()


api.add_resource(SolverController, '/')
api.add_resource(BodyInfo, '/body-info/<string:body_name>/<string:time>')
api.add_resource(Bodies, '/bodies')

if __name__ == '__main__':
    create_logger()
    app.run(debug=True)
