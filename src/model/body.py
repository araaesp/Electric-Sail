from numpy import array


class Body:
    def __init__(self, name, mass, position, velocity, radius):
        """
        Class constructor

        Parameters:
            mass: float. Body mass in kg
            position: list or numpy array. [X, Y, Z] position of the body in km
            velocity: list or numpy array. [VX, VY, VZ] velocity of the body in km/s
            radius: float. Body radius in km
        """
        self.name = name
        self.mass = mass
        self.position = array(position)
        self.velocity = array(velocity)
        self.radius = radius

    def __str__(self):
        return "{}\nMass: {}\nPosition: {}\nVelocity: {}\nRadius: {}\n"\
            .format(self.name,self.mass, self.position, self.velocity, self.radius)
