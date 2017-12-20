"""
Classes for dealing with amino acid mutation sets
"""
import re
from operator import attrgetter


class VariantCalls:
    def __init__(self, text=None, reference=None, sample=None):
        """ Construct a set of Mutations given two aligned amino acid sequences

        :param str reference: the wild-type reference
        :param sample: amino acids present at each position, either a string or
        a list of strings
        """
        self.reference = reference
        if text is not None:
            terms = text.split()
            self.mutation_sets = frozenset(MutationSet(term) for term in terms)
        else:
            if len(reference) != len(sample):
                raise ValueError(
                    'Reference length was {} and sample length was {}.'.format(
                        len(reference),
                        len(sample)))

            self.mutation_sets = {MutationSet(pos=i, variants=alt, wildtype=ref)
                                  for i, (alt, ref) in enumerate(zip(sample, reference),
                                                                 1)
                                  if ref != alt}

    def __str__(self):
        return ' '.join(map(str, sorted(self.mutation_sets,
                                        key=attrgetter('pos'))))

    def __repr__(self):
        text = str(self)
        return 'VariantCalls({!r})'.format(text)

    def __eq__(self, other):
        return self.mutation_sets == other.mutation_sets

    def __hash__(self):
        return hash(self.mutation_sets)


class Mutation(object):
    """Mutation has optional wildtype, position, and call"""

    def __init__(self, text=None, wildtype=None, pos=None, variant=None):
        """ Initialize.

        :param str text: will be parsed for wildtype (optional), position,
            and variant
        :param str wildtype: amino acid abbreviation for wild type
        :param str|int pos: position
        :param str variant: single amino acid abbreviation, or 'i' for
            insertion, or 'd' for deletion
        """
        if text is not None:
            match = re.match(r"([A-Z]?)(\d+)([idA-Z])", text)
            if match is None:
                raise ValueError('Mutation text expects wild type (optional), '
                                 'position, and one variant.')

            if match.group(0) != text:
                # user probably supplied ambiguous variant def
                raise ValueError('Mutation text only allows one variant.')

            wildtype, pos, variant = match.groups()
        self.wildtype = wildtype or None
        self.pos = int(pos)
        self.variant = variant

    def extract_wildtype(self, seq):
        """I really don't like this; please don't actually use this"""
        self.wildtype = seq[self.pos - 1]

    def __repr__(self):
        text = str(self)
        return "Mutation({!r})".format(text)

    def __str__(self):
        text = self.wildtype or ''
        text += '{}{}'.format(self.pos, self.variant)
        return text

    def __eq__(self, other):
        if self.pos != other.pos:
            return False

        if self.wildtype is not None and other.wildtype is not None:
            # if the wt is specified for wt and variant, they must match
            # otherwise the user is doing something weird
            if self.wildtype != other.wildtype:
                message = 'Wild type mismatch between {} and {}.'.format(self,
                                                                         other)
                raise ValueError(message)

        # now that we agree on the wt and position
        return (self.pos, self.variant) == (other.pos, other.variant)

    def __hash__(self):
        return hash((self.pos, self.variant))


class MutationSet(object):
    """Handle sets of mutations at a position"""

    def __init__(self,
                 text=None,
                 wildtype=None,
                 pos=None,
                 variants=None,
                 mutations=None):
        """ Initialize

        :param str text: will be parsed for wildtype (optional), position,
            and variants
        :param str wildtype: amino acid abbreviation for wild type
        :param int|str pos: position
        :param str variants: zero or more amino acid abbreviations, or 'i' for
            insertion, or 'd' for deletion
        :param mutations: a sequence of Mutation objects, with matching
            positions and wild types
        """
        if text:
            match = re.match(r"([A-Z]?)(\d+)([idA-Z]*)", text)
            if match is None:
                raise ValueError

            wildtype, pos, variants = match.groups()

        if variants:
            mutations = frozenset(Mutation(wildtype=wildtype,
                                           pos=pos,
                                           variant=variant)
                                  for variant in variants)
        else:
            mutations = frozenset(mutations)
            positions = {mutation.pos for mutation in mutations}
            wildtypes = {mutation.wildtype for mutation in mutations}
            if pos is not None:
                positions.add(pos)
            wildtypes.add(wildtype)
            wildtypes.discard(None)
            if len(positions) > 1:
                message = 'Multiple positions found: {}.'.format(
                    ', '.join(map(str, sorted(positions))))
                raise ValueError(message)
            if not wildtypes and not mutations:
                raise ValueError('No wildtype and no variants.')
            if not positions:
                raise ValueError('No position and no variants.')
            if len(wildtypes) > 1:
                message = 'Multiple wildtypes found: {}.'.format(
                    ', '.join(sorted(wildtypes)))
                raise ValueError(message)
            pos = positions.pop()
            if wildtypes:
                wildtype = wildtypes.pop()
        self.wildtype = wildtype or None
        self.pos = int(pos)
        self.mutations = mutations

    def __len__(self):
        return len(self.mutations)

    def __contains__(self, call):
        for mutation in self.mutations:
            if call == mutation:
                return True
        return False

    def __eq__(self, other):
        if self.pos != other.pos:
            return False
        if self.wildtype is not None and other.wildtype is not None:
            # if the wt is specified for wt and variant, they must match
            # otherwise the user is doing something weird
            if self.wildtype != other.wildtype:
                message = 'Wild type mismatch between {} and {}.'.format(self,
                                                                         other)
                raise ValueError(message)
        return self.mutations == other.mutations

    def __hash__(self):
        return hash((self.pos, self.mutations))

    def __iter__(self):
        self._mu_iter = list(self.mutations).__iter__()
        return self._mu_iter

    def __next__(self):
        return self._mu_iter.__next__()

    def __reversed__(self):
        return list(self.mutations).__reversed__()

    def __str__(self):
        text = self.wildtype or ''
        text += str(self.pos)
        text += ''.join(sorted(mutation.variant for mutation in self.mutations))
        return text

    def __repr__(self):
        text = str(self)
        return 'MutationSet({!r})'.format(text)
