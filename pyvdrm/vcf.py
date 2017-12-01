'''
Classes for dealing with amino acid mutation sets
'''
import re

class Mutation(object):
    '''Mutation has optional wildtype, position, and call'''

    def __init__(self, pos, variant, wildtype=None):

        self.wildtype = wildtype
        self.pos = int(pos)
        self.variant = variant

    @classmethod
    def from_string(cls, residue):
        '''Parse mutation from string'''
        match = re.match(r"([A-Z]?)(\d+)([idA-Z])", residue)
        if match is None:
            # raise a proper pyvcf exception
            raise ValueError

        if match.group(0) != residue:
            # user probably supplied ambiguous variant def
            raise ValueError

        wildtype = None if not match.group(1) else match.group(1)

        return cls(match.group(2), match.group(3), wildtype)

    def setWt(self, seq):
        '''I really don't like this; please don't actually use this'''
        self.wildtype = seq[self.pos - 1]

    def __repr__(self):
        if self.wildtype:
            return "{}{}{}".format(self.wildtype, self.pos, self.variant)
        return "{}{}".format(self.pos, self.variant)

    def __eq__(self, other):
        if self.pos != other.pos:
            return False

        if self.wildtype is not None and other.wildtype is not None:
            # if the wt is specified for wt and variant, they must match
            # otherwise the user is doing something weird
            if self.wildtype != other.wildtype:
                raise ValueError

        # now that we agree on the wt and position
        return (self.pos, self.variant) == (other.pos, other.variant)

    def __hash__(self):
        return hash((self.pos, self.variant))


class MutationSet(object):
    '''Handle sets of mutations at a position'''

    def __init__(self, mutations):

        # test that mutations is non-empty before this
        pos = None
        wildtype = None
        for mutation in mutations:
            if pos is None:
                pos = mutation.pos
            else:
                if pos != mutation.pos:
                    # mismatched positions in input mutations
                    raise ValueError
            if wildtype is None:
                wildtype = mutation.wildtype
            else:
                if mutation.wildtype is not None and wildtype is not None:
                    if mutation.wildtype != wildtype:
                        raise ValueError

        self.wildtype = wildtype
        self.pos = pos
        self.mutations = set(mutations)

    @classmethod
    def from_string(cls, residue):
        '''Syntactic sugar for specifying multiple variants
        '''

        match = re.match(r"([A-Z]?)(\d+)([idA-Z]*)", residue)
        if match is None:
            # todo: raise a pyvcf parse error
            raise ValueError

        wildtype = None if not match.group(1) else match.group(1)

        return cls([Mutation(match.group(2), amino, wildtype) for amino in
                    match.group(3)])

    def __len__(self):
        return len(self.mutations)

    def __contains__(self, call):
        for mutation in self.mutations:
            if call == mutation:
                return True
        return False

    def __iter__(self):
        self._mu_iter = list(self.mutations).__iter__()
        return self._mu_iter

    def __next__(self):
        return self._mu_iter.__next__()

    def __reversed__(self):
        return list(self.mutations).__reversed__()

    def __repr__(self):
        rep = str(self.pos)
        if self.wildtype:
            rep = self.wildtype + rep

        for mutation in self.mutations:
            rep += mutation.variant

        return rep
