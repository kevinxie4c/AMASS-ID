import numpy as np

a = np.load("data/body_models/smplh/male/model.npz")
np.savetxt("data/parent_of_vertex.tx", np.argmax(a["weights"], 1))
