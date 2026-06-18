'''
    Helper functions for plotting.
'''
from typing import TYPE_CHECKING, Tuple, Optional

import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from rdkit.Chem import Draw
from rdkit import Chem

if TYPE_CHECKING:
    from healer.application.tree_builder import SplitNode


def plot_retrosynthesis_tree(
    root: 'SplitNode',
    figsize: Tuple[float, float] = (6, 4),
    name_prop: str = '_Name',  # 'image', 'smiles', or Mol property
    font_size: int = 8,
    save_path: Optional[str] = None,
    save_dpi: int = 300
) -> None:
    '''
        Simple Matplotlib-only tree, using either images or text labels per name_prop.

        Args:
            root: the root SplitNode of the retrosynthesis tree
            figsize: size of the figure in inches
            name_prop: property to use for node labels; 'image' for images, 'smiles' for SMILES,
                       or any valid Mol property (e.g. '_Name')
            font_size: font size for text labels
            save_path: optional path to save the figure
            save_dpi: DPI for saving the figure
    '''
    positions: dict[object, tuple[float, float]] = {}   # maps objects to (x, y) positions
    def compute_positions(obj, depth=0) -> int:
        from healer.application.tree_builder import SplitNode, SplitBranch
        if isinstance(obj, SplitNode) and not obj.children:
            x = compute_positions.counter
            compute_positions.counter += 1
            positions[obj] = (x, -depth)
            return 1
        if isinstance(obj, SplitNode):
            total = 0
            xs = []
            for branch in obj.children:
                width = compute_positions(branch, depth+1)
                total += width
                xs.append(positions[branch][0])
            positions[obj] = (sum(xs)/len(xs), -depth)
            return total
        if isinstance(obj, SplitBranch):
            total = 0
            xs = []
            for child in obj.children:
                width = compute_positions(child, depth+1)
                total += width
                xs.append(positions[child][0])
            positions[obj] = (sum(xs)/len(xs), -depth)
            return total
        raise TypeError

    compute_positions.counter = 0
    compute_positions(root)

    fig, ax = plt.subplots(figsize=figsize)
    for obj, (x, y) in positions.items():
        from healer.application.tree_builder import SplitNode
        if isinstance(obj, SplitNode):
            # edges to branches
            for branch in obj.children:
                bx, by = positions[branch]
                ax.plot([x, bx], [y, by], 'k-')
            # node label/image
            if name_prop == 'image':
                img = Draw.MolToImage(obj.molecule, size=(150,100))
                im = OffsetImage(img, zoom=0.5)
                ab = AnnotationBbox(im, (x,y), frameon=False)
                ax.add_artist(ab)
            else:
                if name_prop == 'smiles':
                    label = Chem.MolToSmiles(obj.molecule)
                else:
                    label = obj.molecule.GetProp(name_prop) if obj.molecule.HasProp(name_prop) else f'node_{id(obj)%1000}'
                ax.text(x, y, label, ha='center', va='center', fontsize=font_size,
                        bbox=dict(boxstyle='round', fc='white', ec='black', lw=0.5))
        else:
            # SplitBranch edges to split-nodes
            for child in obj.children:
                cx, cy = positions[child]
                ax.plot([x, cx], [y, cy], 'k--')
            rxn_name = getattr(obj.step.reaction, 'name', '')
            ax.text(x, y, rxn_name, ha='center', va='center', fontsize=font_size,
                    bbox=dict(boxstyle='ellipse', fc='lightgray', ec='black', lw=0.5))
    ax.axis('off')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=save_dpi)
    plt.show()

