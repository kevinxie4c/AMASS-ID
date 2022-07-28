import sys
import os
import argparse
import torch
import numpy as np

from human_body_prior.tools.omni_tools import copy2cpu as c2c
from human_body_prior.body_model.body_model import BodyModel
from numpy import linalg

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--threshold', action='store', default=0.1, type=float, help='velocity threshold for detecting contact points (default: 0.1)')
parser.add_argument('-o', '--outdir', action='store', default="output",  help='output folder (default: output)')
parser.add_argument('pose_file')
args = parser.parse_args()

bindir = os.path.abspath(os.path.dirname(__file__))

amass_npz_fname = args.pose_file
outdir = args.outdir
threshold = args.threshold

comp_device = torch.device('cpu')
#amass_npz_fname = 'SFU/0005/0005_Walking001_poses.npz'
bdata = np.load(amass_npz_fname)
subject_gender = bdata['gender']

print('Data keys available:%s'%list(bdata.keys()))
bm_fname = '{}/data/body_models/smplh/{}/model.npz'.format(bindir, subject_gender)
dmpl_fname = '{}/data/body_models/dmpls/{}/model.npz'.format(bindir, subject_gender)

num_betas = 16
num_dmpls = 8

bm = BodyModel(bm_fname=bm_fname, num_betas=num_betas, num_dmpls=num_dmpls, dmpl_fname=dmpl_fname).to(comp_device)
faces = c2c(bm.f)

frame_length = len(bdata['trans'])

body_parms = {
        'root_orient': torch.Tensor(bdata['poses'][:, :3]).to(comp_device), # controls the global root orientation
        'pose_body': torch.Tensor(bdata['poses'][:, 3:66]).to(comp_device), # controls the body
        'pose_hand': torch.Tensor(bdata['poses'][:, 66:]).to(comp_device), # controls the finger articulation
        'trans': torch.Tensor(bdata['trans']).to(comp_device), # controls the global body position
        'betas': torch.Tensor(np.repeat(bdata['betas'][:num_betas][np.newaxis], repeats=frame_length, axis=0)).to(comp_device), # controls the body shape. Body shape is static
        'dmpls': torch.Tensor(bdata['dmpls'][:, :num_dmpls]).to(comp_device) # controls soft tissue dynamics
        }

print('Body parameter vector shapes: \n{}'.format(' \n'.join(['{}: {}'.format(k,v.shape) for k,v in body_parms.items()])))
print('frame_length = {}'.format(frame_length))

body_pose_beta = bm(**{k:v for k,v in body_parms.items() if k in ['pose_body', 'betas', 'trans', 'root_orient']})

# the frame rate of the SFU mocap data is 120 fps
dt = 1/120

par_vtx = np.loadtxt(bindir + '/data/parent_of_vertex.txt')

low_pts = []

for i in range(0, frame_length - 1):
    dx = body_pose_beta.v[i+1] - body_pose_beta.v[i]
    vel = dx / dt
    vel_norm = linalg.norm(vel, axis=1)
    indices = np.where(vel_norm < threshold)[0]
    if len(indices) > 0:
        pts = body_pose_beta.v[i][indices]
        min_idx = np.argmin(pts[:,2])   # gravity force is in -z direction
        low_pts.append(pts[min_idx,:])

A = np.stack(low_pts, axis=0)
p = np.linalg.pinv(A).dot(np.ones(A.shape[0]))
print("plane equation: {}x+{}y+{}z-1=0".format(p[0], p[1], p[2]))
#n = p / np.linalg.norm(p)
#d = -1 / np.linalg.norm(p)
n = np.array([0.0, 0.0, 1.0]) # why the commented code above didn't work?
d = 0

f_node = open(outdir + '/contact_nodes.txt', 'w')
f_pt_idx = open(outdir + '/contact_point_indices.txt', 'w')

for i in range(0, frame_length - 1):
    dx = body_pose_beta.v[i+1] - body_pose_beta.v[i]
    pts = body_pose_beta.v[i].numpy()
    vel = dx / dt
    vel_norm = linalg.norm(vel, axis=1)
    lFoot = []
    rFoot = []
    indices = np.where(np.logical_and(vel_norm < threshold, np.abs(pts.dot(n) + d) < 0.02))[0]
    f_pt_idx.write(' '.join(str(x) for x in indices) + "\n")
    for j in indices:
        x = par_vtx[j]
        if int(x) == 4 or int(x) == 7 or int(x) == 10:
            lFoot.append(pts[j])
        if int(x) == 5 or int(x) == 8 or int(x) == 11:
            rFoot.append(pts[j])
    #contact_nodes = []
    if len(lFoot) > 0:
        #v = np.mean(np.stack(lFoot), axis=0)
        #contact_nodes.append('lFoot {} {} {}'.format(v[0], v[1], v[2]))
        f_node.write('lFoot {} '.format(len(lFoot)))
        for v in lFoot:
            f_node.write('{} {} {} '.format(v[0], v[1], v[2]))
    if len(rFoot) > 0:
        #v = np.mean(np.stack(rFoot), axis=0)
        #contact_nodes.append('rFoot {} {} {}'.format(v[0], v[1], v[2]))
        f_node.write('rFoot {} '.format(len(rFoot)))
        for v in rFoot:
            f_node.write('{} {} {} '.format(v[0], v[1], v[2]))
    #f_node.write(' '.join(contact_nodes) + "\n")
    f_node.write("\n")

f_node.close()
f_pt_idx.close() 
