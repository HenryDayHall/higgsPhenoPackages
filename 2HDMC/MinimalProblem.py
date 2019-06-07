# need to turn off display so that pyplot works
import matplotlib
matplotlib.use('Agg')
import subprocess
import os
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

def MC_run(number_tries):
    """ Generate points and run MinimalProblem on them to get their validity and their Br(H->bb)
    

    Parameters
    ----------
    number_tries : int
        number of random points to check

    Returns
    -------
    input_points : numpy array of floats
        the random positions tried

    checks: numpy array of bools
       is each point permited

    brHbb: numpy array of floats
        the branching ratio of each point

    """
    # set of parmeters for chosing points
    mH_min = 150.; mH_max = 800.
    mA_min = 250.; mA_max = 900.
    sba_value = 1.
    tb_min = 0.; tb_max = 25.
    m12_2_min = 0.; m12_2_max = (mA_max**2)*tb_max/(1 + tb_max**2)
    # chose the points
    mH = np.random.uniform(mH_min, mH_max, number_tries)
    mA = np.random.uniform(mA_min, mA_max, number_tries)
    mC = mA
    sba = np.full(number_tries, sba_value)
    tb = np.random.uniform(tb_min, tb_max, number_tries)
    m12_2 = np.random.uniform(m12_2_min, m12_2_max, number_tries)
    input_vals = np.vstack((mH, mA, mC, sba, tb, m12_2)).T
    string_inputs = input_vals.astype(str)
    # set of arrays we are going to create
    checks = np.empty((number_tries, 4), dtype=bool)
    brHbb = np.empty((number_tries, 1), dtype=float)
    # other instructions for the program
    program_name = "./MinimalProblem"
    yukawa = '1'
    number_outputs = 5
    for i, inps in enumerate(string_inputs):
        command = [program_name] + list(inps) + [yukawa]
        out = subprocess.check_output(command).split()[-number_outputs:]
        checks[i] = np.array(out[:4], dtype=int).astype(bool)
        brHbb[i] = float(out[-1])
    return input_vals, checks, brHbb


def main():
    num_points = 1000
    input_vals, checks, brHbb = MC_run(num_points)
    np.savetxt("input_values.csv", input_vals,
               header="Chosen input values; mH, mA, sba, tb, m12_2")
    np.savetxt("checks.csv", checks,
               header="Validity checks; paramerters_set stability unitarity perturbativity")
    np.savetxt("BrHbb.csv", brHbb, header="Branching ratio H -> b ~b")
    # allowed is the combination of all checks
    allowed = np.all(checks, axis=1)
    # now make the probelms plots
    # mH vs BrHbb
    plt.scatter(input_vals[allowed, 0], brHbb[allowed])
    plt.xlabel("m_H"); plt.ylabel("Br(H -> b b)")
    plt.title("Type I")
    plt.savefig("mH_vs_BrHbb.png")
    plt.clf()
    # mH vs mA
    plt.scatter(input_vals[allowed, 0], input_vals[allowed, 1],
                c=input_vals[allowed, 3])
    plt.xlabel("m_H"); plt.ylabel("m_A")
    plt.colorbar(label="tan(beta)")
    plt.title("Type I")
    plt.savefig("mH_vs_mA.png")

if __name__ == '__main__':
    main()

