# In[]
from random import seed
from Python_module_Quang import *
from scipy.stats import qmc

import numpy as np
import os
import pandas as pd
import string
# import pandas.DataFrame as df

# In[]
# @staticmethod
# def load_file(directory, filename, globals=None):
#     if os.path.splitext(filename)[1] == "":
#         filename = filename + ".txt"
#     if globals is None:
#         globals = dict()
#     globals.update({"__builtins__": None})
#     with open(os.path.join(str(directory), filename), "r") as infile:
#         return eval(infile.read(), globals, {})

# def rev_min_max_func(scaled_val, min, max):
#     """ inverse scaled from min to max """
#     # print("np.shape(scaled_val): ", np.shape(scaled_val))
#     y = np.zeros((np.shape(scaled_val)[0], np.shape(scaled_val)[1]))
#     # print("shape of y:", np.shape(y))
#     # print("range(np.shape(scaled_val)[0]):",
#     #       list(range(np.shape(scaled_val)[0])))
#     # print("range(np.shape(scaled_val)[1]):",
#     #       list(range(np.shape(scaled_val)[1])))
#     # row: number of samples
#     for i in range(np.shape(scaled_val)[0]):
#         # column: number of parameters
#         for j in range(np.shape(scaled_val)[1]):
#             # print("i,j: ", i, " and ", j)
#             # print("y before: ", y)
#             y[i, j] = scaled_val[i, j] * min[j] + (max[j] - min[j])
#             # print("y afer: ", y)
#     return y


def sample_set_RBniCS_formater(inv_sample):
    sample_set_string_i = ""
    sample_set_string_j = ""
    for i in range(np.shape(inv_sample)[0]):
        # print("np.shape(inv_sample)[0]", np.shape(inv_sample)[0])
        for j in range(np.shape(inv_sample)[1]):
            # print("np.shape(inv_sample)[0]", np.shape(inv_sample)[1])
            sample_set_string_j = sample_set_string_j + f"{inv_sample[i,j]}, "
            # print(f"sample_set_string_j: {sample_set_string_j}")
        sample_set_string_j = sample_set_string_j[:-2:]
        # if i == int(np.shape(inv_sample)[0]):
        #     sample_set_string_i = f"{sample_set_string_i}, ({sample_set_string_j})"
        # else:

        if np.shape(inv_sample)[1] == 1:
            # sample_set_string_i = f"{sample_set_string_i} ({sample_set_string_j},),"
            sample_set_string_i = sample_set_string_i + " (" + sample_set_string_j + ",),"
        else:
            sample_set_string_i = f"{sample_set_string_i} ({sample_set_string_j}),"

        # print(f"sample_set_string_i: {sample_set_string_i}")
        sample_set_string_j = ""
    # sample_set_string_j = sample_set_string_j[:0:]
    sample_set_string_i = sample_set_string_i[:-1:]
    sample_set_string_i = f"[{sample_set_string_i}]"
    sample_set_string_i = ''.join(sample_set_string_i.split(' ', 1))
    return sample_set_string_i


# In[] Latin Hypercube sampling

# def Latin_Hypercube_sampling(
#     mu_range,
#     data_path,
#     random_state=0,
#     sample_size=40,
#     model_type="Elastic",  # "Elastic", "Heat"
#     type_of_sample_set="testing_set",
# ):
#     # print(
#     #     f"Latin Hypercube sampling {type_of_sample_set}: {configuration}_w{load_case}"
#     # )

#     dimension, _ = np.shape(mu_range)
#     # print(dimension)
#     sampler = qmc.LatinHypercube(d=dimension,
#                                  seed=random_state)  # number of parameters
#     sample = sampler.random(n=sample_size)  # number of samples
#     # print(mu_range[0][0])
#     print(f"mu_range: {mu_range}")
#     l_bounds = mu_range[0][0]
#     u_bounds = mu_range[0][1]
#     sample = qmc.scale(sample, l_bounds, u_bounds)
#     # print(f"min(sample): {min(sample)}")
#     # print(f"max(sample): {max(sample)}")
#     # print("sample: ", sample)
#     # inv_sample = rev_min_max_func(
#     #     scaled_val=sample,
#     #     min=[mu_range[i][0] for i in range(np.shape(mu_range)[0])],
#     #     max=[mu_range[i][1] for i in range(np.shape(mu_range)[0])])
#     sample_save_string = sample_set_RBniCS_formater(sample)

#     with open(f'{data_path}/{type_of_sample_set}.txt', 'w') as f:
#         f.write(sample_save_string)

# def check_create_dir(data_RB_path):
#     if os.path.exists(data_RB_path) == False:
#         os.mkdir(data_RB_path)


def _Latin_Hypercube_sampling(
    mu_range,
    data_path,
    random_state=None,
    sample_size=40,
    model_type="Elastic",  # "Elastic", "Heat"
    type_of_sample_set="testing_set",
):
    dimension, _ = np.shape(mu_range)
    sampler = qmc.LatinHypercube(d=dimension,
                                 seed=random_state)  # number of parameters
    sample = sampler.random(n=sample_size)  # number of samples

    print(f"mu_range: {mu_range}")
    l_bounds = [mu_range[i][0] for i in range(np.shape(mu_range)[0])]
    u_bounds = [mu_range[i][1] for i in range(np.shape(mu_range)[0])]
    sample = qmc.scale(sample, l_bounds, u_bounds)

    df = pd.DataFrame(
        sample, columns=[f'mu_{i}' for i in range(np.shape(mu_range)[0])])
    df.to_csv(f'{data_path}/{type_of_sample_set}.csv', index=False)

    sample_save_string = sample_set_RBniCS_formater(sample)

    with open(f'{data_path}/{type_of_sample_set}.txt', 'w') as f:
        f.write(sample_save_string)


# def check_create_dir(data_RB_path):
#     if os.path.exists(data_RB_path) == False:
#         os.mkdir(data_RB_path)


# In[]
# if __name__ == "__main__":
def Latin_Hypercube_sampling_func(
        modelName,
        mu_range,
        sample_size_RB,
        type_of_sample_set="testing_set",
        model_type="Elastic",  # "Elastic", "Heat"
        # load_cases = [11]  # 11, 22, 12
):
    print(f"Generate Latin Hypercube sampling for {type_of_sample_set}")
    for load_case in [
            11,
            # 22,
            # 12,
    ]:  # 11, 22, 12
        if type_of_sample_set == "training_set":
            # sample_size_RB = 250
            # sample_size_SCM = 250
            random_state_RB = None
            random_state_SCM = 2
        elif type_of_sample_set == "testing_set":
            # sample_size_RB = 40
            # sample_size_SCM = 40
            random_state_RB = None
            random_state_SCM = 20

        data_RB_path = f"{modelName}/{type_of_sample_set}"
        data_RB_path = os.path.join(os.getcwd(), data_RB_path)
        check_create_dir(data_RB_path)

        # RB sample set
        _Latin_Hypercube_sampling(
            mu_range=mu_range,
            random_state=random_state_RB,
            sample_size=sample_size_RB,
            type_of_sample_set=type_of_sample_set,
            model_type=model_type,  # "Elastic", "Heat"
            data_path=data_RB_path)

        # data_SCM_path = f"{modelName}/scm/{type_of_sample_set}"
        # data_SCM_path = os.path.join(os.getcwd(), data_SCM_path)
        # check_create_dir(data_SCM_path)

        # # SCM sample set
        # Latin_Hypercube_sampling(
        #     configuration=configuration,  # "Macedo", "Oliveira"
        #     load_case=load_case,  # 11, 22, 12
        #     mu_range=mu_range,
        #     random_state=random_state_SCM,
        #     sample_size=sample_size_SCM,
        #     type_of_sample_set=type_of_sample_set,
        #     model_type=model_type,  # "Elastic", "Heat"
        #     data_path=data_SCM_path)


# %%
if __name__ == "__main__":
    mu_range = [
        (0.31, 0.46),  # mu0 - radius
        # (0.8, 1.2),  # mu1 - Young modulus, E, fiber
    ]
    sample_size_RB_training = 100
    sample_size_RB_testing = 40
    for sample_size_RB in [
            sample_size_RB_training,
            # sample_size_RB_testing,
    ]:
        if sample_size_RB == sample_size_RB_training:
            type_of_sample_set = "training_set"
        else:
            type_of_sample_set = "testing_set"
        corrector = 11
        N = 30
        modelName = f"ElasticBlock_EIM_{corrector}_N{N}"
        Latin_Hypercube_sampling_func(
            modelName=modelName,
            sample_size_RB=sample_size_RB,
            mu_range=mu_range,
            type_of_sample_set=type_of_sample_set,
            model_type="Elastic",  # "Elastic", "Heat"
            # load_cases = [11]  # 11, 22, 12
        )
# %%
