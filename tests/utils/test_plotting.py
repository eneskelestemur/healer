"""Tests for retrosynthesis tree plotting."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from healer.application.tree_builder import RetrosynthesisTree  # noqa: E402
from healer.utils.plotting import plot_retrosynthesis_tree  # noqa: E402


class TestPlotting:
    def test_a_built_tree_renders(self, penicillin, all_reactions):
        tree = RetrosynthesisTree(penicillin, all_reactions, max_depth=1)
        tree.build()

        plot_retrosynthesis_tree(tree.root)
        assert plt.gcf().get_axes()
        plt.close("all")

    def test_an_unexpanded_root_renders(self, penicillin, all_reactions):
        tree = RetrosynthesisTree(penicillin, all_reactions, max_depth=0)
        tree.build()

        plot_retrosynthesis_tree(tree.root)
        plt.close("all")

    def test_the_tree_object_exposes_a_display_helper(self, penicillin, all_reactions):
        tree = RetrosynthesisTree(penicillin, all_reactions, max_depth=1)
        tree.build()

        tree.display_tree()
        plt.close("all")
