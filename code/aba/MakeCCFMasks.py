# download from https://github.com/benfulcher/AllenSDK/blob/master/MakeCCFMasks.py
# modified by EJC on 4/12/20

import numpy as np
import csv
import os
import nrrd
import scipy.io as sio
import pandas as pd
import sys

# import matplotlib.pyplot as plt
# %matplotlib inline
#-------------------------------------------------------------------------------
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.structure_tree import StructureTree
from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi
from allensdk.config.manifest import Manifest
from allensdk.core.reference_space import ReferenceSpace
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Global parameters:
resolution = 25
adultMouseStructureGraphID = 1 # 1 is the id of the adult mouse structure graph:

# Set this to what structures to make a mask of:
# Set input/output filenames:
basedir = '/Users/Eli/Dropbox/Neurodegeneration/TauSpread/tau-spread/'
structIDSource = 'data/aba/atlas/structIDs.csv'
outputFilename = 'data/aba/atlas/mask_aba.mat'
basedir = str(sys.argv[1])
structIDSource = str(sys.argv[2])
outputFilename = str(sys.argv[3])
os.chdir(basedir)
print("Making a mask for Oh structures as %s" % outputFilename)

# Set max number of voxels:
maxVoxels = 0; # (0: no max)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
oapi = OntologiesApi()
structure_graph = oapi.get_structures_with_sets([adultMouseStructureGraphID])
# Removes some unused fields returned by the query:
structure_graph = StructureTree.clean_structures(structure_graph)
tree = StructureTree(structure_graph)

# Example:
# tree.get_structures_by_name(['Dorsal auditory area'])
# The annotation download writes a file, so we will need somwhere to put it
annotation_dir = os.path.dirname(structIDSource)
Manifest.safe_mkdir(annotation_dir)
annotation_path = os.path.join(annotation_dir,'annotation.nrrd')

#-------------------------------------------------------------------------------
# Use the connectivity API:
mcapi = MouseConnectivityApi()
# The name of the latest ccf version (a string):
annotation_version = mcapi.CCF_VERSION_DEFAULT
if not os.path.exists(annotation_path):
    mcapi.download_annotation_volume(annotation_version,resolution,annotation_path)
annotation,meta = nrrd.read(annotation_path)


# Build a reference space from a StructureTree and annotation volume, the third argument is
# the resolution of the space in microns
rsp = ReferenceSpace(tree,annotation,[resolution,resolution,resolution])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# So now we're ready to go through structures, and extract their coordinates
structureID_df = pd.read_csv(structIDSource)
structureIDs = structureID_df['ids'].to_numpy()
print("Retrieved %u structures from %s..." % (len(structureIDs),structIDSource))

# A complete mask for one structure
# midPoint = 227 # Middle of brain (z-coordinate)
# coOrds = np.zeros((len(structureIDs),3))
# Assign labels to each voxel according to the list of structureIDs
for ind,sID in enumerate(structureIDs):
    structure_mask = rsp.make_structure_mask([sID])
    if ind==0:
        whole_cortex_mask = np.zeros(structure_mask.shape,dtype=np.uint16)
        # whole_cortex_mask = structure_mask

    # Filter a subset, maxVoxels voxels from the mask
    i = np.nonzero(structure_mask)
    if maxVoxels>0 and maxVoxels<len(i):
        rp = np.random.permutation(len(i))
        rp = rp[:maxVoxels]
        i_filter = i(rp)
    else:
        i_filter = i
    whole_cortex_mask[i_filter] = ind+1
    print("%u / %u: Set %u pixels to %u, %s" % (ind+1,structureIDs.shape[0],i_filter[0].shape[0],ind+1,structureID_df['names'].iloc[ind]))

# np.unique(whole_cortex_mask)

#-------------------------------------------------------------------------------
# Write to h5 file:
sio.savemat(outputFilename,{'mask':whole_cortex_mask,'region_inds':np.arange(len(structureIDs))+1,'region_names':structureID_df['ids'].to_numpy()})
print("Saved mask to %s" % outputFilename)

# View in coronal section
# fig, ax = plt.subplots(figsize=(10, 10))
# plt.imshow(whole_cortex_mask[150, :], interpolation='none', cmap=plt.cm.afmhot)