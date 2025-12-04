# src/mqtt_tios/control.py
""" control.py: graceful shutdown handling for tios CLI tools """
import signal


class GracefulKiller:
    kill_now = False

    def __init__(self):
        signal.signal(signal.SIGINT, self.exit_gracefully)
        signal.signal(signal.SIGTERM, self.exit_gracefully)

    def exit_gracefully(self, signum, frame):
        self.kill_now = True


killer = GracefulKiller()
