""" clients.py: Mqtt-based clients for Tios """
import paho.mqtt.client as mqtt
from queue import SimpleQueue
import string
import random
from .config import config
from time import time


def random_id():
    alphabet = string.ascii_lowercase + string.digits
    return ''.join(random.choices(alphabet, k=6))


ONLINE = 'online'.encode('utf-8')
OFFLINE = 'offline'.encode('utf-8')
EOT = 'EOT'.encode('utf-8')


class MyDict(dict):
    ''' Provides ORM-like behavior '''
    def __init__(self, client, topic, *args):
        self._client = client
        self._topic = topic
        super().__init__(*args)
        
    def __setitem__(self, key, value):
        if not (super().__contains__(key)) or \
                super().__getitem__(key) != value:
            print(f'sending {value} to {key}/{self._topic}')
            self._client.publish(f'tios/{key}/{self._topic}',
                                 value, retain=True)
        super().__setitem__(key, value)

    def __delitem__(self, key):
        if not super().__contains__(key):
            return
        print(f'deleting {key}/{self._topic}')
        self._client.publish(f'tios/{key}/{self._topic}', b'', retain=True)
        super().__delitem__(key)

    def set(self, key, value):
        ''' avoids 'echo' when the data is from the broker '''
        super().__setitem__(key, value)

        
class TiosSimDict():
    ''' An ORM-like object representing data on a Tios broker.

    A monitor has four attributes (state, status, summary, and checkpoint)
    each of which is a dictionary keyed by Tios simulation Id.

    Setting dictionary values leads to a publishing event.

    '''
    def __init__(self,
                 simIds=None,
                 broker_address=None,
                 port=None,
                 username=None,
                 password=None,
                 verbose=False):
        
        self.broker_address = broker_address or config.broker
        self.port = port or config.port
        self.verbose = verbose

        self.topics = []
        if simIds is None:
            simIds = ['+']
            
        if not isinstance(simIds, list):
            simIds = [simIds]
        for simId in simIds:
            self.topics.append((f'tios/{simId}/status', 2))
            self.topics.append((f'tios/{simId}/summary', 2))
            self.topics.append((f'tios/{simId}/checkpoint', 2))
            self.topics.append((f'tios/{simId}/state', 0))
                
        self._client = mqtt.Client(mqtt.CallbackAPIVersion.VERSION2,
                                   client_id='tios_monitor' + random_id())
        self._client.on_connect = self._on_connect
        self._client.on_subscribe = self._on_subscribe
        self._client.on_message = self._on_message
        username = username or config.username
        password = password or config.password
        self._client.username_pw_set(username=username, password=password)
        
        self._status = MyDict(self._client, 'status', {})
        self._summary = MyDict(self._client, 'summary', {})
        self._checkpoint = MyDict(self._client, 'checkpoint', {})
        self._state = MyDict(self._client, 'state', {})
        
        # Connect to the broker
        if self.verbose:
            print(f'Connecting to {broker_address}:{port}')
        self._client.connect(broker_address, port, 60)
        self._client.loop_start()
        
    def _on_connect(self, client, userdata, mid, reason_code, properties):
        if reason_code != 0:
            raise ConnectionError(
                f'Error - connection failed, reason code={reason_code}')
        
        self._client.subscribe(self.topics)
        if self.verbose:
            print(f"Connected to MQTT broker at {self.broker_address}:{self.port}")

    def _on_subscribe(self, client, userdata, mid,
                      reason_code_list, properties):
        for reason_code in reason_code_list:
            if reason_code == 128:
                raise ConnectionError(
                    f'Error - subscription failed, reason={reason_code}')
            if self.verbose:
                print(f'Subscription succesful {reason_code}')

    def _on_message(self, client, userdata, msg):
        if self.verbose:
            print(f"Received message on topic '{msg.topic}':"
                  f" {len(msg.payload)} bytes") 
        _, simId, tag = msg.topic.split('/')
        if tag == 'status':
            self.status.set(simId, msg.payload)
        
        elif tag == 'summary':
            self.summary.set(simId, msg.payload)

        elif tag == 'checkpoint':
            self.checkpoint.set(simId, msg.payload)

        elif tag == 'state':
            self.state.set(simId, msg.payload)

    @property
    def status(self):
        return self._status
        
    @property
    def state(self):
        return self._state

    @property
    def checkpoint(self):
        return self._checkpoint

    @property
    def summary(self):
        return self._summary
        
    def delete(self, simId):
        ''' delete a simulation from the database '''
        del self.status[simId]
        del self.checkpoint[simId]
        del self.summary[simId]
        del self.state[simId]
        
    def close(self):
        """Close the MQTT connection."""
        self._client.disconnect()
        self._client.loop_stop()


class TiosSim():
    ''' An ORM-like object representing a simulation on a Tios broker.

    A simulation has four attributes (state, status, summary, and checkpoint)
    state is a SimpleQueue, the others are of Type bytes.


    Setting attribute values leads to a publishing event.

    '''
    def __init__(self,
                 simId,
                 broker_address=None,
                 port=None,
                 username=None,
                 password=None,
                 verbose=False):
        
        self.simId = simId
        self.broker_address = broker_address or config.broker
        self.port = port or config.port
        self.verbose = verbose

        self.topics = [(f'tios/{simId}/status', 2),
                       (f'tios/{simId}/summary', 2),
                       (f'tios/{simId}/checkpoint', 2),
                       (f'tios/{simId}/state', 0)]
                
        self._client = mqtt.Client(mqtt.CallbackAPIVersion.VERSION2,
                                   client_id='tios_sim' + random_id())
        self._client.on_connect = self._on_connect
        self._client.on_subscribe = self._on_subscribe
        self._client.on_message = self._on_message
        username = username or config.username
        password = password or config.password
        self._client.username_pw_set(username=username, password=password)
        
        self._status = None
        self._summary = None
        self._checkpoint = None
        self._state = SimpleQueue()
        self._poke = None

        self.exists = False
        
        # Connect to the broker
        if self.verbose:
            print(f'Connecting to {broker_address}:{port}')
        self._client.connect(broker_address, port, 60)
        self._client.loop_start()
        
    def _on_connect(self, client, userdata, mid, reason_code, properties):
        if reason_code != 0:
            raise ConnectionError(
                f'Error - connection failed, reason code={reason_code}')
        
        self._client.subscribe(self.topics)
        if self.verbose:
            print(f"Connected to MQTT broker at {self.broker_address}:{self.port}")

    def _on_subscribe(self, client, userdata, mid, reason_code_list, properties):
        for reason_code in reason_code_list:
            if reason_code == 128:
                raise ConnectionError(
                    f'Error - subscription failed, reason={reason_code}')
            if self.verbose:
                print(f'Subscription succesful {reason_code}')

    def _on_message(self, client, userdata, msg):
        self.exists = True
        if self.verbose:
            print(f"Received message on topic '{msg.topic}':"
                  f" {len(msg.payload)} bytes") 
        _, simId, tag = msg.topic.split('/')
        if tag == 'status':
            self._status = msg.payload
        
        elif tag == 'summary':
            self._summary = msg.payload

        elif tag == 'checkpoint':
            self._checkpoint = msg.payload

        elif tag == 'state':
            self._state.put(msg.payload)

        elif tag == 'poke':
            self._poke = msg.payload

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self, value):
        # self._status = value
        self._client.publish(f'tios/{self.simId}/status', value, 2,
                             retain=True)
        
    @property
    def state(self):
        return self._state

    @state.setter
    def state(self, value):
        # self._state.put(value)
        self._client.publish(f'tios/{self.simId}/state', value, 0,
                             retain=True)

    @property
    def checkpoint(self):
        return self._checkpoint

    @checkpoint.setter
    def checkpoint(self, value):
        # self._checkpoint = value
        self._client.publish(f'tios/{self.simId}/checkpoint', value, 2,
                             retain=True)

    @property
    def summary(self):
        return self._summary

    @summary.setter
    def summary(self, value):
        # self._summary = value
        self._client.publish(f'tios/{self.simId}/summary', value, 2,
                             retain=True)
        
    @property
    def poke(self):
        return self._poke
    
    @poke.setter
    def poke(self):
        self._client.publish(f'tios/{self.simId}/poke', time.now())
        
    def delete(self):
        ''' delete a simulation from the database '''
        self.status = b''
        self.checkpoint = b''
        self.summary = b''
        self.state = b''
        
    def close(self):
        """Close the MQTT connection."""
        self._client.disconnect()
        self._client.loop_stop()