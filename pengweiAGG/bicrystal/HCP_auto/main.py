from graingrowth import Scheduler, SimulationKernel, ObervationKernel, CostKerkel

"""
export LD_LIBRARY_PATH=/home/zhongjing/work/moose-projects/mystudy/lib:$LD_LIBRARY_PATH
mpirun -np 20 ./panda-opt -i paramaterCopper.i
"""

if __name__ == "__main__":

    base = "./"
    cmd = "mpirun -np 20 ./panda-opt -i bicrystal_HCP.i > bicrystal_HCP_Ti_phi45.log"

    simkernel = SimulationKernel(base, cmd=cmd, executable="panda-opt", script_filename_1="bicrystal_HCP.i", script_filename_2="grn_2_rand_2D.tex")
    print("[Done] Simulation kernel is ready")
    # trial = {"sigma": 1.0}
    # task_folder_name, task_folder_path = simkernel.prepare(trial=trial)
    # simkernel.simulate(task_folder_path=task_folder_path)
    # obkernel = ObervationKernel(base)
    # print("[Done] Observation kernel is ready")

    # sim_obj = simkernel.post_processing(task_folder_path=task_folder_path)
    # ob_obj = obkernel.post_processing()

    # costkernel = CostKerkel()
    # print("[Done] Cost kernel is ready")
    # costvalue = costkernel.criterion(sim_obj, ob_obj)

    # print("cost value = %.4f" % costvalue)

    scheduler = Scheduler(simkernel=simkernel)
    print("[Done] Scheduler is ready")

    scheduler.run()
