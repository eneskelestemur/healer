'''
    Retrosynthesis decomposition tree structure.
'''
import random
from typing import List, Tuple

from rdkit import Chem

import healer.utils.utils as utils
from healer.domain.composition import RetroStep, CompositionPath
from healer.domain.reaction_template import ReactionTemplate21
from healer.utils.plotting import plot_retrosynthesis_tree


class SplitNode:
    '''
        A node in the retrosynthesis decomposition tree, holding one molecule
        and zero-or-more outgoing SplitBranch entries.
    '''
    def __init__(self, molecule: Chem.Mol) -> None:
        self.molecule: Chem.Mol = molecule
        self.children: List[SplitBranch] = []

class SplitBranch:
    '''
        Associates a RetroStep with its two child SplitNode instances.
    '''
    def __init__(self, step: RetroStep, children: Tuple[SplitNode, SplitNode]) -> None:
        self.step: RetroStep = step
        self.children: Tuple[SplitNode, SplitNode] = children


class RetrosynthesisTree:
    '''
        Builds a retrosynthesis tree for a given molecule using specified
        reactions, up to a maximum depth and with optional filters:  
        - fragment-size

        Time & space complexity
        ------------------------
        Let  
        R = number of reactions tried per node  
        S = average number of splits returned per reaction  
        b = R x S    (average branching factor)  
        d = maximum tree depth  

        Then the total number of nodes is  
        N = 1 + b + b^2 + … + b^d = O(b^d)  

        - Time complexity: O(b^d)  
        - Space complexity: O(b^d)  
    '''
    def __init__(
        self,
        root_molecule: Chem.Mol,
        reactions: List[ReactionTemplate21],
        max_depth: int=1,
        min_heavy_atoms: int = 3,
    ) -> None:
        '''
            Initialize the retrosynthesis tree with a root molecule and 
            a list of reactions.

            Args:
                root_molecule: The starting molecule to decompose.
                reactions: List of ReactionTemplate21 objects to use for retrosynthesis.
                max_depth: Maximum depth of the retrosynthesis tree.
                min_heavy_atoms: Minimum number of heavy atoms in reactants to consider them valid.
        '''
        self.root: SplitNode = SplitNode(root_molecule)
        self.reactions: List[ReactionTemplate21] = reactions
        self.max_depth: int = max_depth
        self.min_heavy_atoms: int = min_heavy_atoms

    def build(self) -> None:
        '''
            Recursively expand the tree from the root.
        '''
        self._expand(self.root, depth=0)

    def _expand(self, node: SplitNode, depth: int) -> None:
        '''
            Recursively expand the retrosynthesis tree from the given node.
            This function applies all reactions to the current molecule and
            creates child nodes for each valid split.

            Args:
                node: The current SplitNode to expand.
                depth: Current depth in the retrosynthesis tree.
        '''
        if depth >= self.max_depth:
            return

        for rxn in self.reactions:
            for (r1, r2) in rxn.run_retro(node.molecule):
                if not self._passes_filters(r1, r2):
                    continue

                step = RetroStep(
                    product=node.molecule,
                    reaction=rxn,
                    reactants=(r1, r2),
                )
                child1 = SplitNode(r1)
                child2 = SplitNode(r2)
                branch = SplitBranch(step, (child1, child2))
                node.children.append(branch)

                self._expand(child1, depth + 1)
                self._expand(child2, depth + 1)

    def _passes_filters(self, r1: Chem.Mol, r2: Chem.Mol) -> bool:
        '''
            Check whether the given pair of reactants meets all criteria.
            Currently filters on minimum heavy atom count but can be extended.
        '''
        # if (r1.GetNumHeavyAtoms() < self.min_heavy_atoms or r2.GetNumHeavyAtoms() < self.min_heavy_atoms):
        #     return False
        if any(not utils.sanitize_mol(r)[0] for r in (r1, r2)):
            return False
        if any(r.GetNumHeavyAtoms() < self.min_heavy_atoms for r in (r1, r2)):
            return False
        return True
    
    def get_composition_paths(self, randomize: bool=False, random_seed: int=-1) -> List[CompositionPath]:
        '''
            Traverse the built tree and return all retrosynthetic composition paths. 
            This will deduplicate CompositionPaths on the fly. The final list will 
            be sorted by the number of fragments in each path (ascending).

            NOTE: Fragment order matters in CompositionPaths, for examples:
                - (A, B) is different from (B, A)

            Args:
                randomize: If True, shuffle the order of branches at each node.
                random_seed: Seed for randomization (if enabled). If < 0, no seed is set.

            Returns:
                List of CompositionPath instances representing retrosynthetic paths.
                Each CompositionPath includes:
                    - steps: the sequence of RetroStep instances applied
                    - fragments: the final fragment tuple
        '''
        results: List[CompositionPath] = []
        seen: set[Tuple[str, ...]] = set()

        def _recurse(node: SplitNode) -> List[Tuple[List[Chem.Mol], List[RetroStep]]]:
            all_res: List[Tuple[List[Chem.Mol], List[RetroStep]]] = [([node.molecule], [])]
            for branch in node.children:
                left = _recurse(branch.children[0])
                right = _recurse(branch.children[1])
                for (fragL, stepsL) in left:
                    for (fragR, stepsR) in right:
                        frags = fragL + fragR
                        steps = [branch.step] + stepsL + stepsR
                        all_res.append((frags, steps))
            return all_res

        for (frags, steps) in _recurse(self.root):
            if not steps:
                continue
            sig = tuple(Chem.MolToSmiles(m) for m in frags)
            if sig in seen:
                continue
            seen.add(sig)
            results.append(CompositionPath(steps=tuple(steps), fragments=tuple(frags)))

        results.sort(key=lambda x: len(x.fragments))
        if randomize:
            if random_seed >= 0:
                random.seed(random_seed)
            else:
                random.seed()
            random.shuffle(results)
        return results

    def display_tree(self, **kwargs) -> None:
        '''
            Visualize the retrosynthesis tree.
        '''
        plot_retrosynthesis_tree(self.root, **kwargs)

