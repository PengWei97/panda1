import sys
sys.path.append("/home/pengwei/projects/panda/pengweiAGG/bicrystal/HCP") # 修改为当前路径
import os
import hashlib
import shutil
import json
import scheduler as qpkernel
# import quams.utilities as tools
import numpy as np
import optuna
import logging
import scipy.stats as st
import pandas as pd
import jinja2

class SimulationKernel(qpkernel.SimulationKernel):

    def __init__(self, base=None, executable=None, cmd=None, script_filename_1=None, script_filename_2=None) -> None:
        super().__init__()
        self.base = base
        self.executable = executable
        self.cmd = cmd
        self.script_filename_1 = script_filename_1
        self.script_filename_2 = script_filename_2

    def prepare(self, trial):

        executable = os.path.join(self.base, "data", self.executable)
        config = os.path.join(self.base, "data", self.script_filename_1)

        info_str = json.dumps(trial)
        sha = hashlib.sha1(info_str.encode("utf-8"))
        task_folder_name = "tasks"

        if not os.path.exists(os.path.join(self.base, task_folder_name)):
            os.mkdir(os.path.join(self.base, task_folder_name))
            print("[Done] %s created" % task_folder_name)
    
        shutil.copyfile(executable, os.path.join(self.base, task_folder_name, self.executable))
        os.system("chmod +x %s" % os.path.join(self.base, task_folder_name, self.executable))

        with open(config, "r") as f1:
            config_content = "".join([line for line in f1])

        # 
        my_loading = trial["my_loading"]
        Euler_angle = trail["Euler_angle"]
        template = jinja2.Template(config_content)

        with open(os.path.join(self.base, task_folder_name, self.script_filename_1), "w") as f:
            config_content = template.render(my_GBmob0=my_loading, my_filename='FCC_Cu_phi45_A{0}_L{1}'.format(Euler_angle, my_loading))
            f.write(config_content)
        print("Input file_1 has been modified")
        
        with open(os.path.join(self.base, "tasks", task_folder_name, self.script_filename_2), "w", encoding="utf-8") as f1:
            f1.write("Texture File" + "\n" + "\n" + "File generated from MATLAB"+ "\n" + "B 2" +"\n" + "    0.00   %s   0.00   1.00"%Euler_angle + "\n" + "    0.00   0.00   0.00   1.00")
        print("Input file_2 has been modified")

        return task_folder_name, os.path.join(self.base, task_folder_name)

    def simulate(self, task_folder_path):
        pwd = os.getcwd()
        os.chdir(task_folder_path)
        os.system(self.cmd)
        os.chdir(pwd)

    def post_processing(self, task_folder_path):
        resultfile = "copper_graingrowth/simu_copper_graingrowth.csv"
        resultfile = os.path.join(task_folder_path, resultfile)
        if not os.path.exists(resultfile):
            return 5

        data = pd.read_csv(resultfile)
    
        return {
            "time": data["time"],
            "grainarea": 50.0 - data["gr0_area"] / 100.0
        }

class ObervationKernel(qpkernel.ObservationKernel):

    def __init__(self, base=None) -> None:
        super().__init__()
        self.base = base

    def post_processing(self):
        expfile = os.path.join(self.base, "data", "exp.txt")
        data = pd.read_csv(expfile)

        return {
            "time": data["time"] / 60.0,
            "grainarea": data["gr0_area"]
        }

class CostKerkel(qpkernel.CostKernel):

    def __init__(self) -> None:
        super().__init__()

    def criterion(self, sim_obj, ob_obj):

        predicted_time = np.interp(ob_obj["grainarea"], sim_obj["grainarea"], sim_obj["time"])
        
        return np.sum((predicted_time - ob_obj["time"]) ** 2)

    def stop_criterion(self):
        pass

class Scheduler(qpkernel.Scheduler):

    def __init__(self, simkernel=None, obkernel=None, costkerkel=None, postkernel=None, study_name="default") -> None:
        super().__init__(simkernel=simkernel, obkernel=obkernel, costkerkel=costkerkel, postkernel=postkernel)

        optuna.logging.get_logger("optuna").addHandler(logging.StreamHandler(sys.stdout))
        # study_name = "example-study"  # Unique identifier of the study.
        storage_name = "sqlite:///{}.db".format(study_name)
        # search_space = {"my_loading": list(np.linspace(0.0, 50.0, 6)), "Euler_angle": list(np.linspace(0.0, 180, 19))}
        search_space = {"my_loading": [10], "Euler_angle": [10,20]} # 正交
        self.study = optuna.create_study(study_name=study_name, storage=storage_name, sampler=optuna.samplers.GridSampler(search_space))

    def objective(self, trial_):

        trial = self.prepare(trial_)
        
        task_folder_name, task_folder_path = self.simkernel.prepare(trial)
        self.simkernel.simulate(task_folder_path)
        # simresult = self.simkernel.post_processing(task_folder_path)
        # expresult = self.obkernel.post_processing()

        return 0

    def prepare(self, trial):
        my_loading = trial.suggest_loguniform("my_loading", 0, 50)
        Euler_angle = trial.suggest_loguniform("Euler_angle", 0, 50)

        return {"my_loading": my_loading, "Euler_angle": Euler_angle}

    def run(self, ntrial=500):

        self.study.optimize(self.objective, n_trials=ntrial)

    def post_processing(self):
        pass

    def best_result(self):
        params = self.study.best_trial.params
        info_str = json.dumps(params)
        sha = hashlib.sha1(info_str.encode("utf-8"))
        task_folder_name = sha.hexdigest()
        return params, task_folder_name