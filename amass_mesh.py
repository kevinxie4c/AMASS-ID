import torch
import numpy as np

from human_body_prior.tools.omni_tools import copy2cpu as c2c
from human_body_prior.body_model.body_model import BodyModel

class AMASSmesh:

    def __init__(self, amass_npz_fname):
        #print(self)
        #print(amass_npz_fname)
        #comp_device = torch.device('cpu')  # will throw an "invalid pointer" exception if embedded in C++
        bdata = np.load(amass_npz_fname)
        subject_gender = bdata['gender']
        print('Data keys available:%s'%list(bdata.keys()))
        bm_fname = 'data/body_models/smplh/{}/model.npz'.format(subject_gender)
        dmpl_fname = 'data/body_models/dmpls/{}/model.npz'.format(subject_gender)

        num_betas = 16
        num_dmpls = 8

        bm = BodyModel(bm_fname=bm_fname, num_betas=num_betas, num_dmpls=num_dmpls, dmpl_fname=dmpl_fname)
        faces = c2c(bm.f)

        frame_length = len(bdata['trans'])

        #body_parms = {
        #        'root_orient': torch.Tensor(bdata['poses'][:, :3]).to(comp_device), # controls the global root orientation
        #        'pose_body': torch.Tensor(bdata['poses'][:, 3:66]).to(comp_device), # controls the body
        #        'pose_hand': torch.Tensor(bdata['poses'][:, 66:]).to(comp_device), # controls the finger articulation
        #        'trans': torch.Tensor(bdata['trans']).to(comp_device), # controls the global body position
        #        'betas': torch.Tensor(np.repeat(bdata['betas'][:num_betas][np.newaxis], repeats=frame_length, axis=0)).to(comp_device), # controls the body shape. Body shape is static
        #        'dmpls': torch.Tensor(bdata['dmpls'][:, :num_dmpls]).to(comp_device) # controls soft tissue dynamics
        #        }
        body_parms = {
                'root_orient': torch.Tensor(bdata['poses'][:, :3]), # controls the global root orientation
                'pose_body': torch.Tensor(bdata['poses'][:, 3:66]), # controls the body
                'pose_hand': torch.Tensor(bdata['poses'][:, 66:]), # controls the finger articulation
                'trans': torch.Tensor(bdata['trans']), # controls the global body position
                'betas': torch.Tensor(np.repeat(bdata['betas'][:num_betas][np.newaxis], repeats=frame_length, axis=0)), # controls the body shape. Body shape is static
                'dmpls': torch.Tensor(bdata['dmpls'][:, :num_dmpls]) # controls soft tissue dynamics
                }

        print('Body parameter vector shapes: \n{}'.format(' \n'.join(['{}: {}'.format(k,v.shape) for k,v in body_parms.items()])))
        print('frame_length = {}'.format(frame_length))
        body_pose_beta = bm(**{k:v for k,v in body_parms.items() if k in ['pose_body', 'betas', 'pose_hand', 'dmpls', 'trans', 'root_orient']})
        #body_pose_beta = bm(**{k:v for k,v in body_parms.items() if k in ['pose_body', 'betas',]})
        self.faces = faces
        self.vertices = body_pose_beta.v.numpy()
