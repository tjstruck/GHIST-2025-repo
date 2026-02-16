import yaml
import numpy as np

# submissionfile = 'final_submissions/growth/growth_3559911.yaml'
# goldstandardfile = 'groundtruth/GHIST_2025_growth_final_goldstandard.yaml'

def relative_root_mean_squared_error(true, pred):
    n = len(true) # update
    squared_error = np.square((true - pred) / true)
    rrmse = np.sqrt(np.sum(squared_error))
    return rrmse

# with open(submissionfile) as stream:
#     try:
#         submission = yaml.safe_load(stream)
#     except yaml.YAMLError as exc:
#         print(exc)

# with open(goldstandardfile) as stream:
#     try:
#         goldstandard = yaml.safe_load(stream)
#     except yaml.YAMLError as exc:
#         print(exc)

# keys = list(goldstandard['parameters'].keys())
# keys.sort()
# RRMSE = relative_root_mean_squared_error(np.array([goldstandard['parameters'][key] for key in keys]), np.array([submission['parameters'][key] for key in keys]))

# # print(RRMSE)