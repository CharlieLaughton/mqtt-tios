import zlib
import numpy as np
from mdtraj.utils import box_vectors_to_lengths_and_angles
from mdtraj.formats import NetCDFTrajectoryFile, XTCTrajectoryFile
from mqttutils import MqttReader


class TiosXTCWriter():
    def __init__(self, broker_address, sim_id, xtcfilename, port=1883,
                 timeout=60):
        # Set the broker address and port
        self.broker_address = broker_address
        self.port = port
        self.subscription = f"tios/{sim_id}/state"
        self.xtcfilename = xtcfilename

        self._reader = MqttReader(broker_address, self.subscription,
                                  port=port, timeout=timeout)
        self.xtcfile = None
        self.framebuffer = None
        self.saved_frames = 0

    def write_frame(self):
        msg = self._reader.readmessage()
        if msg is None:
            print('End of transmission')
            return
        value = msg.payload
        self.framebuffer = np.frombuffer(
            zlib.decompress(value),
            dtype=np.float32).reshape((-1, 3))
        if self.xtcfile is None:
            self.xtcfile = XTCTrajectoryFile(self.xtcfilename, 'w')
        xyz = self.framebuffer[4:]
        box = self.framebuffer[1:4]
        t = self.framebuffer[0, 0]
        self.saved_frames += 1
        self.xtcfile.write(xyz, time=t, step=self.saved_frames, box=box)

    def timedout(self):
        return self._reader.timedout

    def close(self):
        self.xtcfile.close()
        self._reader.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class TiosNCWriter():
    def __init__(self, broker_address, sim_id, ncfilename, port=1883,
                 timeout=60):
        # Set the broker address and port
        self.broker_address = broker_address
        self.port = port
        self.subscription = f"tios/{sim_id}/state"
        self.ncfilename = ncfilename

        self._reader = MqttReader(broker_address, self.subscription,
                                  port=port, timeout=timeout)
        self.ncfile = None
        self.framebuffer = None
        self.saved_frames = 0

    def write_frame(self):
        msg = self._reader.readmessage()
        if msg is None:
            print('End of transmission')
            return
        value = msg.payload
        self.framebuffer = np.frombuffer(
            zlib.decompress(value),
            dtype=np.float32).reshape((-1, 3))
        if self.ncfile is None:
            self.ncfile = NetCDFTrajectoryFile(self.ncfilename, 'w')
        xyz = self.framebuffer[4:]
        box = self.framebuffer[1:4]
        a, b, c, alpha, beta, gamma = box_vectors_to_lengths_and_angles(*box)
        t = self.framebuffer[0, 0]
        self.saved_frames += 1
        self.ncfile.write(xyz, time=t,
                          cell_lengths=(a, b, c),
                          cell_angles=(alpha, beta, gamma))

    def timedout(self):
        return self._reader.timedout

    def close(self):
        self.ncfile.close()
        self._reader.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
