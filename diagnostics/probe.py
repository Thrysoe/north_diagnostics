import json
from diagnostics.diagnostic import Diagnostic
from utils.file_reader import read_probe_data
import os

class Probe(Diagnostic):
    TOTAL_PROBES = 50

    _position_cache = {}
    _mapping_cache = {}
    _data_cache = {}  # Format: {shot: (time, bias_voltage, {channel: I})}

    @classmethod
    def get_position_cache(cls):
        return cls._position_cache

    @classmethod
    def get_mapping_cache(cls):
        return cls._mapping_cache

    @classmethod
    def get_data_cache(cls):
        return cls._data_cache

    def __init__(self, path: str, shot: int, number: int, caching: bool = False):
        super().__init__(path, shot)
        self.number = number
        self.caching = caching

        utils_path = os.path.join(os.path.dirname(__file__), '../utils')
        self.position_path = os.path.join(utils_path, 'probe_positions.json')
        self.mapping_path = os.path.join(utils_path, 'probe_mappings_new.json') if self.shot > 9831 else os.path.join(utils_path, 'probe_mappings__deprecated.json')

        self._active = None
        self._position = None
        self._channel = None
        self._data_loaded = False

        if self.caching:
            self._load_position_cache()
            self._load_mapping_cache()
            self._load_data_cache()  # Now loads multiple channels

    def _load_position_cache(self):
        if self.shot not in Probe._position_cache:
            print(f"🔄 Loading position data for Shot #{self.shot}")
            with open(self.position_path, 'r') as file:
                Probe._position_cache[self.shot] = json.load(file)["probes"]

    def _load_mapping_cache(self):
        if self.shot not in Probe._mapping_cache:
            print(f"🔄 Loading mapping data for Shot #{self.shot}")
            with open(self.mapping_path, 'r') as file:
                Probe._mapping_cache[self.shot] = json.load(file)

    def _load_data_cache(self):
        """ Load all active probe channels into cache if caching is enabled. """
        if self.shot not in Probe._data_cache:
            print(f"🔄 Bulk loading probe data for Shot #{self.shot}")
            channels = []
            mapping = Probe._mapping_cache[self.shot]
            for probe_id, channel in mapping.items():
                try:
                    channels.append(int(channel))
                except (TypeError, ValueError):
                    continue  # Skip unmapped or non-numeric entries
            print(f"read channels: {channels}")
            try:
                time, bias_voltage, current_dict = read_probe_data(self.path, self.shot, channels)
                Probe._data_cache[self.shot] = (time, bias_voltage, current_dict)
            except FileNotFoundError:
                print(f"⚠️  Data file not found for shot {self.shot} at path {self.path}")
                Probe._data_cache[self.shot] = (None, None, {})

    @property
    def active(self):
        if self._active is None:
            self.position  # This will load the position and set _active
        return self._active

    @property
    def position(self):
        if self._position is None:
            if self.caching:
                probe_info = Probe._position_cache[self.shot].get(str(self.number))
            else:
                print(f"🔄 Loading configuration for Shot #{self.shot} : Probe {self.number}")
                with open(self.position_path, 'r') as file:
                    positions = json.load(file)["probes"]
                probe_info = positions.get(str(self.number))

            if probe_info:
                self._active = probe_info["active"]
                x, y, z = probe_info["position"]
                self._position = {"x": x, "y": y, "z": z, "r": (x**2 + y**2)**0.5}
            else:
                self._position = None
                self._active = False
        return self._position

    @property
    def channel(self):
        if self._channel is None:
            if self.caching:
                self._channel = Probe._mapping_cache[self.shot].get(str(self.number))
            else:
                print(f"🔄 Loading channel for Shot #{self.shot} : Probe {self.number}")
                with open(self.mapping_path, 'r') as file:
                    mapping = json.load(file)
                self._channel = mapping.get(str(self.number))

            if self._channel is None:
                print(f"⚠️  Probe {self.number} does not have a mapped channel.")
        return self._channel

    @property
    def time(self):
        if not self._data_loaded:
            self._load_data()
        return self._time

    @property
    def bias_voltage(self):
        if not self._data_loaded:
            self._load_data()
        return self._bias_voltage

    @property
    def current(self):
        if not self._data_loaded:
            self._load_data()
        return self._current

    def _load_data(self):
        if self.active and not self._data_loaded:
            if self.caching:
                print(f"🔄 Accessing cached data for Shot #{self.shot} : Probe {self.number}")
                time, bias_voltage, current_dict = Probe._data_cache[self.shot]
                ch = int(self.channel) if self.channel is not None else None
                self._time = time
                self._bias_voltage = bias_voltage
                self._current = current_dict.get(ch)
                if self._current is None:
                    print(f"⚠️  Channel {ch} not found in cached data for shot {self.shot}")
            else:
                print(f"🔄 Loading data for Shot #{self.shot} : Probe {self.number}")
                try:
                    self._time, self._bias_voltage, self._current = read_probe_data(self.path, self.shot, self.channel)
                except FileNotFoundError:
                    print(f"⚠️  Data file not found for shot {self.shot} at path {self.path}")
                    self._time = None
                    self._bias_voltage = None
                    self._current = None
            self._data_loaded = True
