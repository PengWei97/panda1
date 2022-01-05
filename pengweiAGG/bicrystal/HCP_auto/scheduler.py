import optuna
import logging

class SimulationKernel(object):
    def __init__(self) -> None:
        super().__init__()

    def prepare(self):
        pass

    def simulate(self):
        pass

    def post_processing(self):
        pass

class ObservationKernel(object):
    def __init__(self) -> None:
        super().__init__()

    def prepare(self):
        pass

    def observe(self):
        pass

    def post_processing(self):
        pass

class CostKernel(object):

    def __init__(self) -> None:
        super().__init__()

    def criterion(self, sim_obj, ob_obj):
        pass

    def stop_criterion(self):
        pass


class Scheduler(object):
    def __init__(self, simkernel=None, obkernel=None, costkerkel=None, postkernel=None) -> None:
        super().__init__()

        self.simkernel = simkernel
        self.obkernel = obkernel
        self.costkerkel = costkerkel
        self.postkernel = postkernel

    def prepare(self):
        pass

    def run(self):
        pass

    def post_processing(self):
        pass