""" clients.py: Mqtt-based clients for Tios """
import paho.mqtt.client as mqtt
import time
from queue import SimpleQueue
import string
import random
from .config import config


def random_id():
    alphabet = string.ascii_lowercase + string.digits
    return ''.join(random.choices(alphabet, k=6))


ONLINE = 'online'.encode('utf-8')
OFFLINE = 'offline'.encode('utf-8')
EOT = 'EOT'.encode('utf-8')


class TiosPublisher():
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

        self._status_topic = f'tios/{simId}/status'
        self._state_topic = f'tios/{simId}/state'
        self._summary_topic = f'tios/{simId}/summary'
        self._checkpoint_topic = f'tios/{simId}/checkpoint'
        self._subscribed_topic = f'tios/{simId}/subscribed'

        self._client = mqtt.Client(mqtt.CallbackAPIVersion.VERSION2,
                                   client_id='tios_publisher_' + random_id())
        self._client.on_connect = self._on_connect
        self._client.on_subscribe = self._on_subscribe
        self._client.on_message = self._on_message
        username = username or config.username
        password = password or config.password
        self._client.username_pw_set(username, password)
        self._client.will_set(self._status_topic, payload=OFFLINE, qos=1,
                              retain=True)

        self._status = None
        self._summary = None
        self._checkpoint = None
        self._state = None
        self.has_subscribers = None

        # Connect to the broker
        try:
            self._client.connect(self.broker_address, self.port, 60)
        except Exception as e:
            print(f'Failed to connect to {self.broker_address}:{self.port}')
            raise e
        self._client.loop_start()

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self, new_status):
        if not isinstance(new_status, bytes):
            raise ValueError('Error - value must be bytes')
        self._client.publish(self._status_topic,
                             payload=new_status,
                             qos=1, retain=True)
        self._status = new_status
        if self.verbose:
            print(f'Status updated to {self._status}')

    @property
    def summary(self):
        return self._summary

    @summary.setter
    def summary(self, new_summary):
        if not isinstance(new_summary, bytes):
            raise ValueError('Error - value must be bytes')
        self._client.publish(self._summary_topic,
                             payload=new_summary,
                             qos=1, retain=True)
        self._summary = new_summary
        if self.verbose:
            print(f'Summary updated to {self._summary}')

    @property
    def checkpoint(self):
        return self._checkpoint

    @checkpoint.setter
    def checkpoint(self, new_checkpoint):
        if not isinstance(new_checkpoint, bytes):
            raise ValueError('Error - value must be bytes')
        self._client.publish(self._checkpoint_topic,
                             payload=new_checkpoint,
                             qos=1, retain=True)
        self._checkpoint = new_checkpoint
        if self.verbose:
            print('Checkpoint saved')

    @property
    def state(self):
        return self._state

    @state.setter
    def state(self, new_state):
        if not isinstance(new_state, bytes):
            raise ValueError('Error - value must be bytes')
        self._client.publish(self._state_topic,
                             payload=new_state,
                             qos=1, retain=False)
        self._state = new_state
        if self.verbose:
            print('New state published')

    def _on_connect(self, client, userdata, mid, reason_code, properties):
        if reason_code != 0:
            raise ConnectionError(
                f'Error - connection failed, reason code={reason_code}')

        self.status = ONLINE
        self._client.subscribe(self._subscribed_topic)
        if self.verbose:
            print(f"Connected to broker at {self.broker_address}:{self.port}")

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
            print(f"Received message on topic '{msg.topic}': "
                  f"{len(msg.payload)} bytes")
        tag = msg.topic.split('/')[2]
        if tag == 'subscribed':
            self.has_subscribers = msg.payload.decode('utf-8') == 'True'
            if self.verbose:
                print(f'Received subscription status: {self.has_subscribers} ')

    def close(self):
        """Close the MQTT connection."""
        self.status = OFFLINE
        self.state = EOT
        time.sleep(1)
        self._client.disconnect()
        self._client.loop_stop()
        if self.verbose:
            print('Client closed')


class TiosSubscriber():
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

        self._status_topic = f'tios/{simId}/status'
        self._state_topic = f'tios/{simId}/state'
        self._summary_topic = f'tios/{simId}/summary'
        self._checkpoint_topic = f'tios/{simId}/checkpoint'
        self._subscribed_topic = f'tios/{simId}/subscribed'

        self._client = mqtt.Client(mqtt.CallbackAPIVersion.VERSION2,
                                   client_id='tios_subscriber_' + random_id())
        self._client.on_connect = self._on_connect
        self._client.on_subscribe = self._on_subscribe
        self._client.on_message = self._on_message
        self._client.username_pw_set(username=username or config.username,
                                     password=password or config.password)
        self._client.will_set(self._subscribed_topic, payload=b'False',
                              qos=1, retain=True)

        self.status = None
        self.summary = None
        self.checkpoint = None
        self.state = None
        self.states = SimpleQueue()
        self._closing = False

        # Connect to the broker
        self._client.connect(self.broker_address, self.port, 60)
        self._client.loop_start()

    def _on_connect(self, client, userdata, mid, reason_code, properties):
        if reason_code != 0:
            raise ConnectionError(
                f'Error - connection failed, reason code={reason_code}')

        self._client.publish(self._subscribed_topic,
                             payload=b'True',
                             qos=1, retain=True)

        self._client.subscribe([(self._status_topic, 2),
                                (self._summary_topic, 2),
                                (self._checkpoint_topic, 2),
                                (self._state_topic, 2),
                                (self._subscribed_topic, 2)])
        if self.verbose:
            print(f"Connected to broker at {self.broker_address}:{self.port}")

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
            print(f"Received message on topic '{msg.topic}': "
                  f"{len(msg.payload)} bytes")
        tag = msg.topic.split('/')[2]
        if tag == 'status':
            self.status = msg.payload
        elif tag == 'summary':
            self.summary = msg.payload
        elif tag == 'checkpoint':
            self.checkpoint = msg.payload
        elif tag == 'state':
            self.state = msg.payload
            self.states.put(self.state)
        elif tag == 'subscribed':
            if msg.payload.decode('utf-8') == 'False':
                if not self._closing:
                    self._client.publish(self._subscribed_topic,
                                         payload=b'True',
                                         qos=1, retain=True)

    def close(self):
        """Close the MQTT connection."""
        self._closing = True
        self._client.publish(self._subscribed_topic,
                             payload=b'False',
                             qos=1, retain=True)
        time.sleep(1)
        self._client.disconnect()
        self._client.loop_stop()
        if self.verbose:
            print('Client closed')


class TiosMonitor():
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
            self.topics.append((f'tios/{simId}/subscribed', 2))

        self._client = mqtt.Client(mqtt.CallbackAPIVersion.VERSION2,
                                   client_id='tios_monitor_' + random_id())
        self._client.on_connect = self._on_connect
        self._client.on_subscribe = self._on_subscribe
        self._client.on_message = self._on_message
        self._client.username_pw_set(username=username or config.username,
                                     password=password or config.password)

        self.status = {}
        self.summary = {}
        self.subscribed = {}

        # Connect to the broker
        self._client.connect(self.broker_address, self.port, 60)
        self._client.loop_start()

    def _on_connect(self, client, userdata, mid, reason_code, properties):
        if reason_code != 0:
            raise ConnectionError(
                f'Error - connection failed, reason code={reason_code}')

        self._client.subscribe(self.topics)
        if self.verbose:
            print(f"Connected to broker at {self.broker_address}:{self.port}")

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
            print(f"Received message on topic '{msg.topic}': "
                  f"{len(msg.payload)} bytes")
        _, simId, tag = msg.topic.split('/')
        if tag == 'status':
            self.status[simId] = msg.payload
        elif tag == 'summary':
            self.summary[simId] = msg.payload
        elif tag == 'subscribed':
            self.subscribed[simId] = msg.payload.decode('utf-8') == 'True'

    def delete(self, simId):
        self._client.publish(f'tios/{simId}/status', b'', retain=True)
        self._client.publish(f'tios/{simId}/checkpoint', b'', retain=True)
        self._client.publish(f'tios/{simId}/summary', b'', retain=True)
        self._client.publish(f'tios/{simId}/subscribed', b'', retain=True)
        self._client.publish(f'tios/{simId}/state', b'', retain=True)
        if self.verbose:
            print(f'Simulation {simId} deleted')

    def refresh(self):
        for topic in self.topics:
            self._client.unsubscribe(topic[0])
        self.status = {}
        self.summary = {}
        self.subscribed = {}
        for topic in self.topics:
            self._client.subscribe(topic)

    def close(self):
        """Close the MQTT connection."""
        self._client.disconnect()
        self._client.loop_stop()
        if self.verbose:
            print('Client closed')
