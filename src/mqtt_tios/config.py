""" config.py: load tios-related environment variables """
import os

class Configuration():
    def __init__(self):
        self.broker = os.getenv('TIOSBROKER', 'localhost')
        self.port = int(os.getenv('TIOSPORT', 1883))
        self.username = os.getenv('TIOSUSERNAME')
        self.password = os.getenv('TIOSPASSWORD')

config = Configuration()
