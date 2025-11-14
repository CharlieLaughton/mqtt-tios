import os

class Configuration():
    def __init__(self):
        self.broker = os.getenv('TIOSBROKER', 'localhost')
        self.port = int(os.environ.get('TIOSPORT', 1883))
        self.username = os.environ.get('TIOSUSERNAME')
        self.password = os.environ.get('TIOSPASSWORD')

config = Configuration()
