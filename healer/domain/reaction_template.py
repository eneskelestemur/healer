'''
    Reaction class to define chemical reacitons using rdkit.
'''
import inspect
from itertools import chain

from rdkit import Chem
from rdkit.Chem import Descriptors, rdChemReactions

from healer.domain.building_block import BuildingBlock


class ReactionTemplate21:
    '''
        Wraps rdkit reaction functions to define chemical reactions.
        This class is specifically designed to handle reactions with
        2 reactants and 1 product. If the reaction has more than 2 reactants
        or more than 1 product, the functions may not work as expected.
    '''
    def __init__(
            self,
            name: str,
            reaction_smarts: str,
            retro_smarts: str,
            display_smarts: str = None,
            rhs_classes: list[int] = tuple(),
            tags: list[str] = tuple(),
            description: str = "",
            long_name: str = None,
            tier: int = None
    ):
        '''
            Constructor for the Reaction class.

            Args:
                name: str, name of the reaction.
                kwargs: additional properties of the reaction.
                    reaction_smarts: str, SMARTS string for the reaction. Same as syn_smarts.
                    retro_smarts: str, SMARTS string for the retro reaction.
                    display_smarts: str, SMARTS string for the reaction.
                    description: str, description of the reaction.
                    long_name: str, long name of the reaction.
                    rhs_classes: list of int, reaction classes.
                    tags: list of str, tags for the reaction.
                    tier: int, tier of the reaction.
        '''
        self.name = name
        self.reaction_smarts = reaction_smarts
        self.retro_smarts = retro_smarts
        self.display_smarts = display_smarts
        self.rhs_classes = rhs_classes
        self.tags = tags
        self.description = description
        self.long_name = long_name
        self.tier = tier

        self.sanitized_ = False
        self._reaction = rdChemReactions.ReactionFromSmarts(reaction_smarts)
        self._reaction.Initialize()

        # check if a reaction is valid
        try:
            _ = rdChemReactions.SanitizeRxn(self._reaction)
            self._reaction.RemoveUnmappedReactantTemplates(0.1)
            self._reaction.RemoveUnmappedProductTemplates(0.1)
            if len(self.get_reactants()) == 2 and len(self.get_products()) == 1:
                self.sanitized_ = True
        except:
            self.sanitized_ = False

    def __str__(self):
        return self.name
    
    def __repr__(self):
        return self.name
    
    def __hash__(self):
        return hash(rdChemReactions.ReactionToSmiles(self._reaction, canonical=True))
    
    @classmethod
    def from_reaction_json(cls, name: str, reaction_json: dict):
        cls_parameters = [
            p.name for p in inspect.signature(cls.__init__).parameters.values()
            if p.name != "self"
        ]

        valid_params = {key: val for key, val in reaction_json.items() if key in cls_parameters}

        # check for the fact that reaction_smarts could be called syn_smarts
        if "syn_smarts" in reaction_json.keys():
            valid_params["reaction_smarts"] = reaction_json["syn_smarts"]

        # add in the name
        valid_params["name"] = name

        return ReactionTemplate21(**valid_params)

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
        self._reaction = rdChemReactions.ReactionFromSmarts(reaction_smarts)
        try:
            flags = rdChemReactions.SanitizeRxn(self._reaction)
            self._reaction.RemoveUnmappedReactantTemplates(0.1)
            self._reaction.RemoveUnmappedProductTemplates(0.1)
            if len(self.get_reactants()) == 2 and len(self.get_products()) == 1:
                self.sanitized_ = True
        except:
            self.sanitized_ = False

    def get_rdkit_reaction_object(self):
        return self._reaction
    
    def get_reactants(self, sort_by_mw=True):
        '''
            Returns the reactants of the reaction.
            If sort_by_mw is True, the reactants are sorted by molecular weight 
            in descending order. Otherwise, the reactants are returned in the order
            they are defined in the reaction.
        '''
        if not sort_by_mw:
            return list(self._reaction.GetReactants())
        return sorted(list(self._reaction.GetReactants()), key=lambda x: Descriptors.MolWt(x), reverse=True)
    
    def get_products(self):
        '''
            Returns the products of the reaction sorted by molecular weight 
            in descending order.
        '''
        return sorted(list(self._reaction.GetProducts()), key=lambda x: Descriptors.MolWt(x), reverse=True)

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
    
    def get_reactant_index(self, mol: Chem.Mol):
        '''
            Returns the index of the reactant in the reaction.
        '''
        reactants = self.get_reactants(False)
        return [i for i, reactant in enumerate(reactants) if mol.HasSubstructMatch(reactant)]
    
    def is_valid(self):
        '''
            Returns True if the reaction template is valid.
        '''
        return self.sanitized_
    
    def is_reactant(self, mol):
        '''
            Returns True if the molecule is a reactant in the reaction.
        '''
        return self._reaction.IsMoleculeReactant(mol)
    
    def is_product(self, mol):
        '''
            Returns True if the molecule is a product in the reaction.
        '''
        return self._reaction.IsMoleculeProduct(mol)
    
    def run_syn(self, *reactants, ):
        '''
            Runs the reaction on the reactants and returns all possible products.
        '''
        assert len(reactants) == 2, "Reaction must have exactly 2 reactants."

        reactants_ordered = self._order_reactants_by_annotations(reactants)
        products = []
        for reactant_pair in reactants_ordered:
            products += self._reaction.RunReactants(reactant_pair, maxProducts=10)
        
        return list(chain(*products))

    def run_retro(self, product):
        '''
            Runs the retro reaction on the product and returns the possible tuples of
            reactants. Each tuple contains the reactants that can form the product
            using the reaction template.
        '''
        retro_reaction = rdChemReactions.ReactionFromSmarts(self.retro_smarts)
        return list(retro_reaction.RunReactants([product], maxProducts=100))

    def _order_reactants_by_annotations(self, reactants):
        '''
            Orders the reactants based on the 'rxn_annotations' property.
            This is used to ensure that the reactants are in the same order
            as they are defined in the reaction template.
        '''
        ann0 = (reactants[0].get_parsed_prop('rxn_annotations').get(self.name, [])
                if isinstance(reactants[0], BuildingBlock) else self.get_reactant_index(reactants[0]))
        ann1 = (reactants[1].get_parsed_prop('rxn_annotations').get(self.name, [])
                if isinstance(reactants[1], BuildingBlock) else self.get_reactant_index(reactants[1]))

        orderings = []
        for p0 in ann0:
            for p1 in ann1:
                if p0 == p1:
                    continue
                out = [None, None]
                out[p0] = reactants[0]
                out[p1] = reactants[1]
                orderings.append(out)

        return orderings

