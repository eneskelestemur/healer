'''
    Reaction class to define chemical reacitons using rdkit.
'''

import os
import json

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


class ReactionTemplate21:
    '''
        Wraps rdkit reaction functions to define chemical reactions.
        This class is specifically designed to handle reactions with
        2 reactants and 1 product. If the reaction has more than 2 reactants
        or more than 1 product, the functions may not work as expected.
    '''
    def __init__(self, name=None, **kwargs):
        '''
            Constructor for the Reaction class.

            Args:
                name: str, name of the reaction.
                kwargs: additional properties of the reaction.
                    reaction_smarts: str, SMARTS string for the reaction. Same as syn_smarts.
                    display_smarts: str, SMARTS string for the reaction.
                    descriptions: str, description of the reaction.
                    long_name: str, long name of the reaction.
                    retro_smarts: str, SMARTS string for the retro reaction.
                    rhs_classes: list of int, reaction classes.
                    tags: list of str, tags for the reaction.
                    tier: int, tier of the reaction.
        '''
        self._name = name
        self.sanitized = False

        # update properties if kwargs are provided
        self.display_smarts = None
        self.descriptions = None
        self.long_name = None
        self.retro_smarts = None
        self.rhs_classes = []
        self.tags = []
        self.tier = None
        for key, value in kwargs.items():
            if key == 'syn_smarts' or key == 'reaction_smarts':
                self.reaction_smarts = value
            setattr(self, key, value)

    def get_reaction_smarts(self):
        '''
            Returns the reaction SMARTS string.
        '''
        return self._reaction_smarts
    
    def set_reaction_smarts(self, reaction_smarts):
        '''
            Sets the reaction SMARTS string.
        '''
        self._reaction_smarts = reaction_smarts
        self._reaction = AllChem.ReactionFromSmarts(reaction_smarts)
        try:
            flags = AllChem.SanitizeRxn(self._reaction)
            self._reaction.RemoveUnmappedReactantTemplates(0.1)
            self._reaction.RemoveUnmappedProductTemplates(0.1)
            if len(self.get_reactants()) == 2 and len(self.get_products()) == 1:
                self.sanitized = True
        except:
            self.sanitized = False

    def get_reaction(self):
        '''
            Returns the rdkit reaction object.
        '''
        return self._reaction
    
    def set_reaction(self, reaction):
        '''
            Can't set the reaction directly. Use set_reaction_smarts instead.
        '''
        raise ValueError("Can't set the reaction directly. Use set_reaction_smarts instead.")

    # TODO: Update the setters to check if the reaction is valid.
    def get_name(self): return self._name
    def set_name(self, name): self._name = name
    
    def get_display_smarts(self): return self._display_smarts
    def set_display_smarts(self, display_smarts): self._display_smarts = display_smarts

    def get_descriptions(self): return self._descriptions
    def set_descriptions(self, descriptions): self._descriptions = descriptions

    def get_long_name(self): return self._long_name
    def set_long_name(self, long_name): self._long_name = long_name

    def get_retro_smarts(self): return self._retro_smarts
    def set_retro_smarts(self, retro_smarts): self._retro_smarts = retro_smarts

    def get_rhs_classes(self): return self._rhs_classes
    def set_rhs_classes(self, rhs_classes): self._rhs_classes = rhs_classes

    def get_tags(self): return self._tags
    def set_tags(self, tags): self._tags = tags

    def get_tier(self): return self._tier
    def set_tier(self, tier): self._tier = tier

    # Properties of the Reaction class:
    reaction_smarts = property(get_reaction_smarts, set_reaction_smarts)
    reaction = property(get_reaction, set_reaction)
    name = property(get_name, set_name)
    display_smarts = property(get_display_smarts, set_display_smarts)
    descriptions = property(get_descriptions, set_descriptions)
    long_name = property(get_long_name, set_long_name)
    retro_smarts = property(get_retro_smarts, set_retro_smarts)
    rhs_classes = property(get_rhs_classes, set_rhs_classes)
    tags = property(get_tags, set_tags)
    tier = property(get_tier, set_tier)

    def __str__(self):
        return self.name
    
    def __repr__(self):
        return self.display_smarts
    
    def __hash__(self):
        return hash(AllChem.ReactionToSmiles(self.reaction, canonical=True))
    
    def get_reactants(self):
        '''
            Returns the reactants of the reaction sorted by molecular weight 
            in descending order.
        '''
        return sorted(list(self.reaction.GetReactants()), key=lambda x: Descriptors.MolWt(x), reverse=True)
    
    def get_products(self):
        '''
            Returns the products of the reaction sorted by molecular weight 
            in descending order.
        '''
        return sorted(list(self.reaction.GetProducts()), key=lambda x: Descriptors.MolWt(x), reverse=True)

    def get_reactants_smarts(self):
        '''
            Returns the SMARTS of the reactants.
        '''
        return [Chem.MolToSmarts(reactant) for reactant in self.get_reactants()]
    
    def get_products_smarts(self):
        '''
            Returns the SMARTS of the products.
        '''
        return [Chem.MolToSmarts(product) for product in self.get_products()]
    
    def get_reactants_smiles(self):
        '''
            Returns the SMILES of the reactants.
        '''
        return [Chem.MolToSmiles(reactant) for reactant in self.get_reactants()]
    
    def get_products_smiles(self):
        '''
            Returns the SMILES of the products.
        '''
        return [Chem.MolToSmiles(product) for product in self.get_products()]
    
    def is_valid(self):
        '''
            Returns True if the reaction template is valid.
        '''
        return self.sanitized
    
    def is_reactant(self, mol):
        '''
            Returns True if the molecule is a reactant in the reaction.
        '''
        return self.reaction.IsMoleculeReactant(mol)
    
    def is_product(self, mol):
        '''
            Returns True if the molecule is a product in the reaction.
        '''
        return self.reaction.IsMoleculeProduct(mol)
    
    def run(self, *reactants):
        '''
            Runs the reaction on the reactants and returns the products.
        '''
        return self.reaction.RunReactants(list(reactants), maxProducts=10)

    def run_retro(self, *products):
        '''
            Runs the retro reaction on the products and returns the reactants.
        '''
        retro_reaction = AllChem.ReactionFromSmarts(self.retro_smarts)
        return retro_reaction.RunReactants(list(products), maxProducts=10)

