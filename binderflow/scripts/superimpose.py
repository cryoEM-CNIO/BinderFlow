import pyrosetta
from pyrosetta import rosetta
from types import ModuleType

prc: ModuleType = pyrosetta.rosetta.core
prp: ModuleType = pyrosetta.rosetta.protocols
pru: ModuleType = pyrosetta.rosetta.utility
prn: ModuleType = pyrosetta.rosetta.numeric
prs: ModuleType = pyrosetta.rosetta.std 
pr_conf: ModuleType = pyrosetta.rosetta.core.conformation
pr_scoring: ModuleType = pyrosetta.rosetta.core.scoring
pr_res: ModuleType = pyrosetta.rosetta.core.select.residue_selector
pr_options: ModuleType = pyrosetta.rosetta.basic.options


def superpose_pose_by_chain(pose, ref, chain: str, strict: bool=False) -> float:
    """
    superpose by PDB chain letter

    :param pose:
    :param ref:
    :param chain:
    :return:
    """
    atom_map = prs.map_core_id_AtomID_core_id_AtomID()
    chain_sele: pr_res.ResidueSelector = pr_res.ChainSelector(chain)
    r: int  # reference pose residue number (Fortran)
    m: int  # mobile pose residue number (Fortran)
    for r, m in zip(pr_res.selection_positions(chain_sele.apply(ref)),
                    pr_res.selection_positions(chain_sele.apply(pose))
                    ):
        if strict:
            assert pose.residue(m).name3() == ref.residue(r).name3(), 'Mismatching residue positions!'
        ref_atom = pyrosetta.AtomID(ref.residue(r).atom_index("CA"), r)

        mobile_atom = pyrosetta.AtomID(pose.residue(m).atom_index("CA"), m)
        atom_map[mobile_atom] = ref_atom
    return pr_scoring.superimpose_pose(mod_pose=pose, ref_pose=ref, atom_map=atom_map)