################################################################################
##                                                                            ##
##    CONTAINS: Methods related to associating segments and/or trajectories   ##
##              corresponding to the same particle. ALL METHODS CURRENTLY     ##
##              RELY ON TRUTH INFORMATION (07/11/2023).                       ##
##    NOTES: "ghdr" refers to a genie_hdr dataset extracted from an hdf5      ##
##           file as seen in common/file_parsing.py                           ##
##           "traj" refers to trajectories dataset extracted from an hdf5     ##
##           file as seen in common/file_parsing.py                           ##
##           Modified from: https://github.com/edhinkle/mesonless_numubarCC/  ##
##           blob/main/common/singleParticleAssociation_methods.py            ##
##                                                                            ##
################################################################################

import sys
sys.path.append('../../')
import truth_studies.util.truth_ixn_methods as truth_ixn_methods
import numpy as np


####------------------- CONNECT TRAJECTORIES USING PDG ---------------------####

# Method below used to connect trajectories of a single particle which have been  
# assigned different traj_ids. If a particle scatters or reinteracts in a way such 
# that multiple child particles are produced, the lineage stops, even if one of 
# those child particles has the same PDG ID as the original particle. We do this by 
# ensuring that the original particle/traj_id (and any other traj_id accepted as 
# being the same particle) does not have the same parent as any other 
# particles/traj_ids (is an "only child") and also does not have multiple children.
def same_pdg_connected_trajectories(track_pdg, track_id, vertex_assoc_traj,\
                                    traj, ghdr):
    traj_id_set = {track_id} # initialize a set of traj_ids associated with a single particle
    
    ## WALK UP THE FAMILY TREE
    this_pdg = track_pdg
    this_track_id = track_id

    while this_pdg==track_pdg: # stop if a member of an older generation has a different PDG ID than the original track/particle
        particle_mask = vertex_assoc_traj['traj_id'] == this_track_id # mask to find trajectories assoc. w/ current track
        parent_track_id = vertex_assoc_traj[particle_mask]['parent_id'] # find parent ID of current track
        parent_mask = vertex_assoc_traj['parent_id'] == parent_track_id # mask to find trajectories w/ current track AS PARENT TRACK
        this_generation =  vertex_assoc_traj[parent_mask] # find all tracks with the same parent as the current track
        if len(this_generation) == 1: # only move forward if current track is an "only child"
            this_pdg = truth_ixn_methods.find_parent_pdg(vertex_assoc_traj[particle_mask]['parent_id'],
                                   vertex_assoc_traj[particle_mask]['vertex_id'],
                                   traj, ghdr) # get parent PDG ID of current track
            if len(this_pdg) != 1: # if parent PDG ID is not a single value, break while loop
                break
            elif this_pdg==track_pdg: # if parent PDG ID of track matches current/original track's PDG ID, add parent to track id set
                this_track_id = vertex_assoc_traj[particle_mask]['parent_id'].tolist()[0] # also makes parent track the new "current" track
                traj_id_set.add(this_track_id)
        else: break # break while loop if current track/particle is not an "only child"

    ## WALK DOWN THE FAMILY TREE
    this_pdg = track_pdg
    this_track_id = track_id

    while this_pdg==track_pdg: # stop if/when a child has a different PDG ID than the original track
        particle_mask = vertex_assoc_traj['parent_id'] == this_track_id  # mask to find trajectories w/ current track AS PARENT TRACK
        child_particle=vertex_assoc_traj[particle_mask] # find all tracks which are children of the current track
        if len(child_particle)==1: # only move forward if current track has only one child
            this_pdg = child_particle['pdg_id']
            if child_particle['pdg_id']==track_pdg: # if the child's PDG ID matches the original track PDG ID, add child track id to set
                this_track_id = child_particle['traj_id'].tolist()[0] # also makes child track the new "current" track
                traj_id_set.add(this_track_id)
        else: break # break while loop if current track/particle has more than one child particle

    return traj_id_set

# Method below used to connect trajectories associated with a single pi0 decay product  
# assigned different traj_ids. If a particle scatters or reinteracts in a way such 
# that non-photons/electron/positrons are produced, the lineage stops. We do this by 
# ensuring that the original particle/traj_id comes from a pi0 parent
# Method above "connect_pi0_decay_prod_shower_traj" is used to recursively find children within
# the pi0 decay product shower.
def recursive_shower_association(valid_pdgs, start_traj_list, traj_id_set, final_traj_id_set, vertex_assoc_traj):

    #print("Start Trajectory List:", start_traj_list)
    #print("Final Trajectory ID set:", final_traj_id_set)
    final_traj_id_set.update(traj_id_set)
    child_track_ids = []

    if len(start_traj_list) == 0:
        return final_traj_id_set
    elif len(start_traj_list) == 1:
        start_traj_id = start_traj_list[0]
        particle_mask = vertex_assoc_traj['traj_id'] == start_traj_id
        start_pdg = vertex_assoc_traj[particle_mask]['pdg_id']
        if start_pdg in valid_pdgs:
            final_traj_id_set.add(start_traj_id)
            child_particle_mask = vertex_assoc_traj['parent_id'] == start_traj_id
            child_particles = vertex_assoc_traj[child_particle_mask]
            child_track_ids = child_particles['traj_id'].tolist()
            #print("Child Trajectory IDs:", child_track_ids)
            if len(child_track_ids) == 0:
                return final_traj_id_set
            else:
                final_traj_id_set.update(recursive_shower_association(valid_pdgs, child_track_ids, traj_id_set, final_traj_id_set, vertex_assoc_traj))
                return final_traj_id_set
        else:
            return final_traj_id_set
    elif len(start_traj_list) > 1:
        for start_traj_id in start_traj_list:
            particle_mask = vertex_assoc_traj['traj_id'] == start_traj_id
            start_pdg = vertex_assoc_traj[particle_mask]['pdg_id']
            if start_pdg in valid_pdgs:
                final_traj_id_set.add(start_traj_id)
                child_particle_mask = vertex_assoc_traj['parent_id'] == start_traj_id
                child_particles = vertex_assoc_traj[child_particle_mask]
                child_track_ids.extend(child_particles['traj_id'].tolist())
            else:
                #print("Child particle with invalid PDG ID found. Children still added to list.")
                #print("Child PDG ID:", start_pdg)
                child_particle_mask = vertex_assoc_traj['parent_id'] == start_traj_id
                child_particles = vertex_assoc_traj[child_particle_mask]
                child_track_ids.extend(child_particles['traj_id'].tolist())
        #print("Child Trajectory IDs:", child_track_ids)
        final_traj_id_set.update(recursive_shower_association(valid_pdgs, child_track_ids, traj_id_set, final_traj_id_set, vertex_assoc_traj))
        return final_traj_id_set
    else:
        print("Unforeseen case. Exiting.")
        print("Length of start_traj_list:", len(start_traj_list))
        return final_traj_id_set


def connect_pi0_decay_prod_shower_traj(track_pdg, track_id, vertex_assoc_traj,\
                                    traj, ghdr):

    traj_id_set = {track_id} # initialize a set of traj_ids associated with a single particle

    # Check that input track has a pi0 parent
    particle_mask = vertex_assoc_traj['traj_id'] == track_id
    #print("Particles:", vertex_assoc_traj[particle_mask]['pdg_id'])
    parent_pdg = truth_ixn_methods.find_parent_pdg(vertex_assoc_traj[particle_mask]['parent_id'],
                                            vertex_assoc_traj[particle_mask]['vertex_id'],
                                            traj, ghdr)
    #print("Parent PDG ID:", parent_pdg)
    if abs(parent_pdg)!=111:
        print("Input track does not have a pi0 parent. No shower trajectories to connect.")
        return traj_id_set # initialize a set of traj_ids associated with a single particle

    ## WALK DOWN THE FAMILY TREE
    this_pdg = track_pdg
    this_track_id = track_id
    child_particle_mask = vertex_assoc_traj['parent_id'] == this_track_id
    child_particles = vertex_assoc_traj[child_particle_mask]
    valid_pdgs = np.array([22, 11, -11]) # valid PDG IDs for pi0 decay products
    final_traj_id_set = set()

    traj_id_set.update(recursive_shower_association(valid_pdgs, child_particles['traj_id'], traj_id_set, final_traj_id_set, vertex_assoc_traj))
    #print("FINAL Trajectory ID set:", traj_id_set)
    return traj_id_set

def shower_pdg_fraction(traj_id_set, vertex_assoc_traj, traj, ghdr):

    total_traj = len(traj_id_set)
    photons = 0.
    electrons = 0.
    for traj in traj_id_set:
        particle_mask = vertex_assoc_traj['traj_id'] == traj
        pdg = vertex_assoc_traj[particle_mask]['pdg_id']
        if pdg == 22:
            photons += 1.
        elif pdg == 11 or pdg == -11:
            electrons += 1.
        else:
            print("Unexpected PDG ID found in shower. PDG ID:", pdg)
            break
    return electrons/total_traj
    

####-------- FIND TRAJECTORY AT END OR BEGINNING OF PARTICLE TRACK ---------####

def find_trajectory_at_vertex(traj_id_set, vertex_assoc_traj,traj, ghdr):
    
    is_prim = truth_ixn_methods.is_primary_particle(traj_id_set, vertex_assoc_traj,traj, ghdr)

    if is_prim == False:
        print("Track ID set does not represent a primary particle. Therefore, \
              there does not exist a trajectory coming from the vertex.")
        return 
    else:
        traj_id_at_vertex = 0 # initialize traj_id_at_vertex variable
        for tid in traj_id_set: # loop through trajectories in set to find trajectory with muon (anti)neutrino parent
            particle_mask = vertex_assoc_traj['traj_id'] == tid
            parent_pdg = truth_ixn_methods.find_parent_pdg(vertex_assoc_traj[particle_mask]['parent_id'],
                                         vertex_assoc_traj[particle_mask]['vertex_id'],
                                         traj, ghdr) 
            
            if abs(parent_pdg)==14:
                traj_id_at_vertex = tid
                break
            else: continue

    return traj_id_at_vertex


def find_forward_primary_particle_end_trajectory(traj_id_set, vertex_assoc_traj,traj, ghdr):
    
    tid_at_vertex = find_trajectory_at_vertex(traj_id_set, vertex_assoc_traj,traj, ghdr)
    particle_mask = vertex_assoc_traj['traj_id'] == tid_at_vertex
    start_pt = vertex_assoc_traj[particle_mask]['xyz_start']

    traj_id_at_end = tid_at_vertex # initialize traj_id_at_vertex variable
    end_z = start_pt[2] # initialize end z value
    for tid in traj_id_set: # loop through trajectories in set to find trajectory with largest z value
        traj_mask = vertex_assoc_traj['traj_id'] == tid
        tid_end = vertex_assoc_traj[traj_mask]['xyz_end']
        if tid_end[2]>end_z:
            end_z = tid_end[2]
            traj_id_at_end = tid
        else: continue

    return traj_id_at_end

# Not all particles are forward-going
def distance_between_points(p1, p2):
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)

def find_primary_particle_end_trajectory(traj_id_set, vertex_assoc_traj,traj, ghdr):
    
    tid_at_vertex = find_trajectory_at_vertex(traj_id_set, vertex_assoc_traj,traj, ghdr)
    particle_mask = vertex_assoc_traj['traj_id'] == tid_at_vertex
    start_pt = vertex_assoc_traj[particle_mask]['xyz_start']
    if len(start_pt) > 1:
        print("Multiple start points found for trajectory ID:", tid_at_vertex)
    else:
        start_pt = start_pt[0]

    traj_id_at_end = tid_at_vertex # initialize traj_id_at_vertex variable
    max_diff = 0.
    for tid in traj_id_set: # loop through trajectories in set to find trajectory with largest z value
        traj_mask = vertex_assoc_traj['traj_id'] == tid
        tid_end = vertex_assoc_traj[traj_mask]['xyz_end']
        if len(tid_end) > 1:
            print("Multiple end points found for trajectory ID:", tid)
            break
        else:
            tid_end = tid_end[0]
        if distance_between_points(tid_end, start_pt)>max_diff:
            traj_id_at_end = tid
        else: continue

    return traj_id_at_end

def find_secondary_particle_end_trajectory(traj_id_set, vertex_assoc_traj, start_traj_id):
    
    particle_mask = vertex_assoc_traj['traj_id'] == start_traj_id
    start_pt = vertex_assoc_traj[particle_mask]['xyz_start']
    if len(start_pt) > 1:
        print("Multiple start points found for trajectory ID:", start_traj_id)
    else:
        start_pt = start_pt[0]


    traj_id_at_end = start_traj_id # initialize traj_id_at_vertex variable
    max_diff = 0.
    for tid in traj_id_set: # loop through trajectories in set to find trajectory with largest z value
        traj_mask = vertex_assoc_traj['traj_id'] == tid
        tid_end = vertex_assoc_traj[traj_mask]['xyz_end']
        if len(tid_end) > 1:
            print("Multiple end points found for trajectory ID:", tid)
            break
        else:
            tid_end = tid_end[0]
        #print("End Point:", tid_end)
        #print("Start Point:", start_pt)
        if distance_between_points(tid_end, start_pt)>max_diff:
            traj_id_at_end = tid
        else: continue

    return traj_id_at_end

def find_forward_secondary_particle_end_trajectory(traj_id_set, vertex_assoc_traj, start_traj_id):
    
    particle_mask = vertex_assoc_traj['traj_id'] == start_traj_id
    start_pt = vertex_assoc_traj[particle_mask]['xyz_start']
    if len(start_pt) > 1:
        print("Multiple start points found for trajectory ID:", start_traj_id)
    else:
        start_pt = start_pt[0]

    traj_id_at_end = start_traj_id # initialize traj_id_at_vertex variable
    end_z = start_pt[2] # initialize end z value
    for tid in traj_id_set: # loop through trajectories in set to find trajectory with largest z value
        traj_mask = vertex_assoc_traj['traj_id'] == tid
        tid_end = vertex_assoc_traj[traj_mask]['xyz_end']
        if len(tid_end) > 1:
            print("Multiple end points found for trajectory ID:", tid)
        else:
            tid_end = tid_end[0]
        if tid_end[2]>end_z:
            end_z = tid_end[2]
            traj_id_at_end = tid
        else: continue

    return traj_id_at_end