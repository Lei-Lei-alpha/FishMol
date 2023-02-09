import numpy as np
from recordclass import make_dataclass, dataobject
# from recordclass import make_dataclass, dataobject, astuple, asdict
import itertools
from .data import elements

class Atom(np.ndarray):
    """
    Array representing atom symbol and coordinates in real space.
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Attributes
    cell : the Cell object, the unit cell associated to the coordinates.
    basis : {'Cartesian', 'Crystal'}, Describes whether the array contains crystal or cartesian coordinates.
    """
    def __new__(cls, symb, pos, cell = None, basis='Cartesian', pbc = (1,1,1)):
        dc_cell = make_dataclass("Cell", 'lattice')
        crys_names = ['Crystal', 'Crys', 'Cr', 'S']
        cart_names = ['Cartesian', 'Cart', 'Ca', 'R']
        coord_names = make_dataclass("Coord_names", 'crys_names cart_names')
        dt = np.dtype([('symb', np.unicode_, 2), ('pos', np.float64, (3,))])
        obj = np.asarray([(symb, pos)], dtype=dt).view(cls)
        obj._symb = symb
        obj._pos = pos
        obj._basis = basis
        obj._coord_names = coord_names(crys_names, cart_names)
        if isinstance(cell, dataobject):
            obj._cell = cell
        else:
            obj._cell = dc_cell(np.asarray(cell))
        obj._pbc = pbc
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._cell = getattr(obj, '_cell', None)
        self._basis = getattr(obj, '_basis', None)
        self._pbc = getattr(obj, '_pbc', (0,0,0))

    @property
    def symb(self):
        return self._symb
    @property
    def pos(self):
        return self._pos
    @property
    def cell(self):
        return self._cell
    @property
    def basis(self):
        return self._basis
    @property
    def coord_names(self):
        return self._coord_names
    @property
    def pbc(self):
        return self._pbc

    def to_cart(self):
        """
        Converts the coordinates to Cartesian and return a new Atom object.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Returns
        out : Atom, new Atom object insured to have basis='Cartesian'.
        """
        if self.basis in self.coord_names.cart_names:
            return self
        else:
            pos = xys2cart(self.pos, self.cell)
            return Atom(symb = self.symb, pos=pos, cell=self.cell, basis=self.coord_names.cart_names[0], pbc = self.pbc)

    def to_crys(self):
        """
        Converts the coordinates to Crystal and return a new Atom object.
        Returns
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        out : Atom, new Atom object insured to have basis='Crystal'.
        """
        if self.basis in self.coord_names.crys_names:
            return self
        else:
            pos = cart2xys(self.pos, self.cell)
            return Atom(symb = self.symb, pos=pos, cell=self.cell, basis=self.coord_names.crys_names[0], pbc = self.pbc)

    def to_basis(self, basis):
        """
        Converts the coordinates to the desired basis and return a new object.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Parameters
        basis : {'Cartesian', 'Crystal'}, basis to which the coordinates are converted.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Returns
        out : Atom, new Atom object insured to have basis=basis.
        """
        
        if basis in self.coord_names.cart_names:
            return self.to_crys()
        elif basis in self.coord_names.cart_names:
            return self.to_cart()
        else:
            raise NameError("Trying to convert to an unknown basis")

    def vec(self, other, mic = False):
        """
        Calculate the vector connecting two Atom.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Parameters
        other : Atom
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Returns
        out : Atom, vector connecting self and other with the same basis as self.
        """
        if mic:
            a2b = other.to_crys().pos - self.to_crys().pos
            for i in range(3):
                a2b[i] = a2b[i] - round(a2b[i])
            if self.basis in self.coord_names.cart_names:
                return xys2cart(a2b, self.cell)
            else:
                return a2b
        else:
            a2b = other.to_crys().pos - self.to_crys().pos
            if self.basis in self.coord_names.cart_names:
                return xys2cart(a2b, self.cell)
            else:
                return a2b
            
    def dist(self, other, mic = False):
        """
        Calculate the distance between two Atom objects.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Parameters
        other : Atom
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Returns
        out : float, the distance between self and other Atoms in angstrom.
        """
        if mic:
            a2b = other.to_crys().pos - self.to_crys().pos
            for i in range(3):
                a2b[i] = a2b[i] - round(a2b[i])
            a2b = xys2cart(a2b, self.cell)
            return np.sqrt(np.dot(a2b, a2b))
        else:
            a2b = other.to_cart().pos - self.to_cart().pos
            return np.sqrt(np.dot(a2b, a2b))

class Atoms(np.ndarray):
    """
    Array representing coordinates in real space under periodic boundary conditions.
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Attributes
    cell : the Cell object, the unit cell associated to the coordinates.
    basis : {'Cartesian', 'Crystal'}, Describes whether the array contains crystal or cartesian coordinates.
    """
    def __new__(cls, symbs, pos, cell = None, basis='Cartesian', pbc = (1,1,1)):
        dt = np.dtype([('symbol', np.unicode_, 2), ('position', np.float64, (3,))])
        dc_cell = make_dataclass("Cell", 'lattice')
        crys_names = ['Crystal', 'Crys', 'Cr', 'S']
        cart_names = ['Cartesian', 'Cart', 'Ca', 'R']
        coord_names = make_dataclass("Coord_names", 'crys_names cart_names')
        if isinstance(symbs, list):
            symbs = np.array(symbs, dtype = "<U2")
        if isinstance(symbs, str):
            a = []
            for i, letter in enumerate(symbs):
                if letter.isdigit():
                    if i == 0:
                        raise Exception("Wrong formmat, a valid chemical formula should start with a chemical symbol!")
                    else:
                        a += [a[i-1] for n in range(int(letter)-1)]
                else:
                    a += [letter]
            symbs = np.array(a, dtype = "<U2")
        pos =  np.array(pos, dtype = np.float64)
        obj = np.asarray([np.array((symb, position), dtype=dt)] for symb, position in zip(symbs, pos)).view(cls)
        obj._symbs = symbs
        obj._pos = pos
        obj._basis = basis
        obj._coord_names = coord_names(crys_names, cart_names)
        if isinstance(cell, dataobject):
            obj._cell = cell
        else:
            obj._cell = dc_cell(np.asarray(cell))
        obj._pbc = pbc
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._cell = getattr(obj, '_cell', None)
        self._basis = getattr(obj, '_basis', None)
        self._pbc = getattr(obj, '_pbc', (0,0,0))

    @property
    def symbs(self):
        return self._symbs
    @property
    def pos(self):
        return self._pos
    @property
    def cell(self):
        return self._cell
    @property
    def basis(self):
        return self._basis
    @property
    def coord_names(self):
        return self._coord_names
    @property
    def pbc(self):
        return self._pbc
    @pos.setter
    def pos(self, pos_lst):
        self._pos = pos_lst
        
    def __len__(self):
        return len(self.symbs)
    
    def __getitem__(self, n):
        """
        Enable slicing of Atoms object.
        """
        if isinstance(n, int):
            if abs(n) >= len(self):
                raise IndexError("The index (%d) is out of range."%n)
            if n < 0:
                n += len(self)
            select = n
            
        elif isinstance(n, str):
            select = [x for x in range(len(self)) if self.symbs[x] == n]
        
        elif isinstance(n, tuple) or isinstance(n, list):
            if all(isinstance(x, str) for x in n):
                select = [x for x in range(len(self)) if self.symbs[x] in n]

            elif all(isinstance(x, int) for x in n):
                select = n
            else:
                raise Exception("Use indices only or chemical symbols only enumerates, mixed indices-chemical symbol selection is not enabled yet!")
        
        elif isinstance(n, slice):
            start = n.start
            stop = n.stop
            step = n.step
            if start is None:
                start = 0
            if stop is None:
                stop = len(self)
            if step is None:
                step = 1
            select = list(range(start,stop,step))
        
        else:
            raise TypeError("Invalid argument type.")
        
        return self.__class__(symbs = self.symbs[select], pos = self.pos[select], cell = self.cell, basis = self.basis, pbc = self.pbc)
    
    def to_cart(self):
        """
        Converts the coordinates to Cartesian and return a new Atoms object.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Returns
        out : Atoms, new Atoms object insured to have basis='Cartesian'.
        """
        if self.basis in self.coord_names.cart_names:
            return self
        else:
            if self.cell is None:
                print("Cell \U0001F35E not defined, cannot convert to absolute atom positions!")
            else:
                self.pos = xys2cart(self.pos, self.cell)
            return self

    def to_crys(self):
        """
        Converts the coordinates to Crystal and return a new Atom object.
        Returns
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        out : Atom, new Atom object insured to have basis='Crystal'.
        """
        if self.basis in self.coord_names.crys_names:
            return self
        else:
            if self.cell is None:
                print("\U0001F35E Cell not defined, cannot convert to fractional atom positions!")
            else:
                self.pos = cart2xys(self.pos, self.cell)
            return self

    def to_basis(self, basis):
        """
        Converts the coordinates to the desired basis and return a new object.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Parameters
        basis : {'Cartesian', 'Crystal'}, basis to which the coordinates are converted.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Returns
        out : Atom, new Atom object insured to have basis=basis.
        """
        if basis in Atoms.crys_names:
            return self.to_crys()
        elif basis in Atoms.cart_names:
            return self.to_cart()
        else:
            raise NameError("Trying to convert to an unknown basis")

    def vec(self, a, b, normalise = False, absolute = True, mic = False):
        """
        Calculate the vector connecting two atoms in the Atoms object with indices m and n.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Parameters
        a, b : int, indices of two atoms in the Atoms object. Required.
        absolute : boolean, if true, calculates the absolute vector in angstrom, otherwise calculates the vector with fractional positions. Default: True. Optional.
        normalise: boolean, if true, calculated the normalised vector, otherwise calculates the vector with a length equals to the distance between the two atoms m and n. Default: False. Optional.
        mic: boolean, if true, wraps the atom within the cell and using minimum image convension. Default: False. Optional.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Returns
        out : Array, vector connecting self and other with the same basis as self.
        """
        pos = self.pos[[a,b]]    
        if mic:
            if self.basis in self.coord_names.cart_names:
                pos = cart2xys(pos, self.cell)
            a2b = pos[1] - pos[0]
            for i in range(3):
                a2b[i] = a2b[i] - round(a2b[i])
            if absolute:
                a2b = xys2cart(a2b, self.cell)
        else:
            a2b = pos[1] - pos[0]
            if not absolute:
                a2b = cart2xys(a2b, self.cell)
        if normalise:
                    a2b = a2b/np.linalg.norm(a2b)
        return a2b
    
    def vecs(self, a, b, normalise = False, absolute = True, mic = False):
        if any([isinstance(a, int),  isinstance(b, int)]):
            combs = itertools.product(a, b)
        elif all([isinstance(a, list), isinstance(b, list)]):
            if len(a) != len(b):
                combs = itertools.product(a, b)
            elif len(a) == len(b):
                combs = zip(a, b)
        a2bs = np.asarray([self.vec(*comb, normalise = normalise, absolute = absolute, mic = mic) for comb in combs])
        return a2bs
            
    def dist(self, a, b, mic = False):
        """
        Calculate the distance between two Atom objects.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Parameters
        a, b : int, the indcices of atoms to be calculated.
        mic : boolean, whether or not to use the minimum image convension.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Returns
        out : float, the distance between a and b Atoms.
        """
        if mic:
            a2b = self.pos[b] - self.pos[a]
            a2b = cart2xys(a2b, self.cell)
            for i in range(3):
                a2b[i] = a2b[i] - round(a2b[i])
            a2b = xys2cart(a2b, self.cell)
        else:
            a2b = self.pos[b] - self.pos[a]
        return np.linalg.norm(a2b)
    
    def dists(self, at_g1, at_g2, cutoff = None, mic = False):
        """
        Calculates the distaces between two atom groups.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Parameters
        at_g1, at_g2 : lists of int, the indcices of atoms groups to be calculated.
        mic : boolean, whether or not to use the minimum image convension.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Returns
        out : float, the distance between a and b Atoms.
        """
        pairs = list(itertools.product(at_g1, at_g2))
        # pairs = [[pair[0], pair[1]] for pair in pairs if not pair[0] == pair[1]]
        distances = np.zeros(len(pairs))
        
        for i, pair in enumerate(pairs):
            distances[i] = self.dist(pair[0], pair[1], mic = mic)

        if cutoff is not None:
            mask = np.where(distances <= cutoff)
            distances = distances[mask]
            pairs = np.array(pairs)[mask]
        return pairs, distances
    
    def angle(self, a, b, c, mic = False):
        """
        Calculates the angle formed by atoms a, b and c, or the angle between vectors ba and bc.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Parameters
        a, b : int, the indcices of atoms to be calculated.
        mic : boolean, whether or not to use the minimum image convension.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Returns
        out : float, the distance between a and b Atoms.
        """
        if mic:
            b2a = self[a].to_crys().pos - self[b].to_crys().pos
            b2c = self[c].to_crys().pos - self[b].to_crys().pos
            for i in range(3):
                b2a[i] = b2a[i] - round(b2a[i])
                b2c[i] = b2c[i] - round(b2c[i])
            b2a = xys2cart(b2a, self.cell)
            b2c = xys2cart(b2c, self.cell)
        else:
            b2a = self[a].to_cart().pos - self[b].to_cart().pos
            b2c = self[c].to_cart().pos - self[b].to_cart().pos
        return np.rad2deg(np.arccos(b2a.dot(b2c)/(np.linalg.norm(b2a)*np.linalg.norm(b2c))))
    
    def angles(self, at_g1, at_g2, at_g3, mic = False):
        """
        Calculates the angle formed by atoms a, b and c, or the angle between vectors ba and bc.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Parameters
        a, b : int, the indcices of atoms to be calculated.
        mic : boolean, whether or not to use the minimum image convension.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Returns
        out : np.array, the angles between atom groups at_g1, at_g2, at_g3.
        """
        pairs = list(itertools.product(at_g1, at_g2, at_g3))
        # pairs = [pair for pair in pairs if len(pair) == len(set(pair))]
        angles = np.zeros(len(pairs))

        for i, pair in enumerate(pairs):
                angles[i] = self.angle(pair[0],pair[1], pair[2], mic = mic)
        return pairs, angles
    
    def dihedral(self, idx, mic = False):
        """
        Calculates the angle formed by atoms a, b and c, or the angle between vectors ba and bc.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Parameters
        idx : list or tuple, the indcices of atoms to be calculated.
        mic : boolean, whether or not to use the minimum image convension.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Returns
        out : np.array, the dihedral angle.
        The dihedral angle is defined by three vectors, if the input idx has four elements, the dihedral angle is defined by 0->1, 1->2, 2->3; if the input idx has six elements, the dihedral angle is defined by 0->1, 2->3, 4->5
        """
        if len(idx) == 3:
            vec1, vec2, vec3 = idx
        elif len(idx) == 4:
            vec1, vec2, vec3 = self.vec(idx[0],idx[1], mic = mic), self.vec(idx[1],idx[2], mic = mic), self.vec(idx[2],idx[3], mic = mic)
        elif len(idx) == 6:
            vec1, vec2, vec3 = self.vec(idx[0],idx[1], mic = mic), self.vec(idx[2],idx[3], mic = mic), self.vec(idx[4],idx[5], mic = mic)
        else:
            raise Exception("Wrong length of indices, please specify exactly four or six atoms.")
        numerator = np.dot(np.cross(vec1, vec3), np.cross(vec2, vec3).T)
        denominator = np.linalg.norm(np.cross(vec1, vec3)) * np.linalg.norm(np.cross(vec2, vec3))
        delta = np.arccos(numerator/denominator)
        return np.rad2deg(delta)
    
    def dihedrals(self, at_g, mic = False):
        """
        Calculates the angle formed by atoms a, b and c, or the angle between vectors ba and bc.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Parameters
        idx : list or tuple, the indcices of atoms to be calculated.
        mic : boolean, whether or not to use the minimum image convension.
        -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Returns
        out : np.array, the dihedral angle.
        The dihedral angle is defined by three vectors, if the input idx has four elements, the dihedral angle is defined by 0->1, 1->2, 2->3; if the input idx has six elements, the dihedral angle is defined by 0->1, 2->3, 4->5
        """
        if len(at_g) == 4:
            pairs = list(itertools.product(at_g[0], at_g[1], at_g[2], at_g[3]))
            # pairs = [pair for pair in pairs if len(pair) == len(set(pair))]
            vec1 = np.asarray([self.vec(idx[0], idx[1], mic = mic) for idx in pairs])
            vec2 = np.asarray([self.vec(idx[1], idx[2], mic = mic) for idx in pairs])
            vec3 = np.asarray([self.vec(idx[2], idx[3], mic = mic) for idx in pairs])
        elif len(at_g) == 6:
            pairs = list(itertools.product(at_g[0], at_g[1], at_g[2], at_g[3], at_g[4], at_g[5]))
            # pairs = [pair for pair in pairs if len(pair) == len(set(pair))]
            vec1 = np.asarray([self.vec(idx[0], idx[1], mic = mic) for idx in pairs])
            vec2 = np.asarray([self.vec(idx[2], idx[3], mic = mic) for idx in pairs])
            vec3 = np.asarray([self.vec(idx[4], idx[5], mic = mic) for idx in pairs])
        else:
            raise Exception("Wrong length of atom groups, please specify exactly four or six atom groups.")
        numerator = np.diagonal(np.dot(np.cross(vec1, vec3), np.cross(vec2, vec3).T)) # drop  np.cross(a, c)[0].dot(np.cross(b, c)[1])
        denominator = np.linalg.norm(np.cross(vec1, vec3), axis = 1) * np.linalg.norm(np.cross(vec2, vec3), axis = 1)
        delta = np.arccos(numerator/denominator)
        return pairs, np.rad2deg(delta)
    
    def calc_com(self):
        """
        Calculates the centre of mass of atoms.
        """
        masses = np.array([data.elements[symb] for symb in self.symbs])
        com = np.dot(masses, self.pos)/masses.sum()
        return com
    
    def calc_cog(self):
        """
        Calculate the centre of geometry.
        """
        return self.pos.mean(axis=0)
        
    def wrap_pos(self, center=(0.5, 0.5, 0.5), pretty_translation=False, eps=1e-7):

        if isinstance(self.pbc[0], int):
            pbc = [True if x == 1 else False for x in self.pbc]
        
        shift = np.asarray(center) - 0.5 - eps

        # Don't change coordinates when pbc is False
        shift[np.logical_not(pbc)] = 0.0

        assert np.asarray(self.cell.lattice)[np.asarray(pbc)].any(axis=1).all(), (self.cell.lattice, pbc)

        fractional = np.linalg.solve(self.cell.lattice.T, np.asarray(self.pos).T).T - shift

        if pretty_translation:
            fractional = translate_pretty(fractional, pbc)
            shift = np.asarray(center) - 0.5
            shift[np.logical_not(pbc)] = 0.0
            fractional += shift
        else:
            for i, periodic in enumerate(pbc):
                if periodic:
                    fractional[:, i] %= 1.0
                    fractional[:, i] += shift[i]
        
        self.pos = xys2cart(fractional, self.cell)

        return self
    
    def at_sel(self, n, inverse_select = False):
        """
        Select by chemical symbols
        """
        if isinstance(n, int):
            if abs(n) >= len(self):
                raise IndexError("The index (%d) is out of range."%n)
            if n < 0:
                n += len(self)
            if not inverse_select:
                select = n
            else:
                select = [x for x in list(range(len(self))) if x != n]

        elif isinstance(n, str):
            if not inverse_select:
                select = [x for x in range(len(self)) if self.symbs[x] == n]
            else:
                select = [x for x in range(len(self)) if self.symbs[x] != n]

        elif isinstance(n, tuple) or isinstance(n, list):
            if all(isinstance(x, str) for x in n):
                if inverse_select:
                    select = [x for x in range(len(self)) if self.symbs[x] not in n]
                else:
                    select = [x for x in range(len(self)) if self.symbs[x] in n]

            elif all(isinstance(x, int) for x in n):
                if inverse_select:
                    select = [x for x in range(len(self)) if x not in n]
                else:
                    select = n
            else:
                raise Exception("Use indices only or chemical symbols only enumerates, mixed indices-chemical symbol selection is not enabled yet!")

        elif isinstance(n, slice):
            start = n.start
            stop = n.stop
            step = n.step
            if start is None:
                start = 0
            if stop is None:
                stop = len(self)
            if step is None:
                step = 1
            select = list(range(start,stop,step))

        else:
            raise TypeError("Invalid argument type.")
        try:
            atoms = self.__class__(symbs = self.symbs[select], pos = self.pos[select], cell = self.cell, basis = self.basis, pbc = self.pbc)
        except IndexError:
                atoms = self.__class__(symbs = self.symbs[list(itertools.chain.from_iterable(select))],
                                       pos = self.pos[list(itertools.chain.from_iterable(idx))], cell = self.cell, basis = self.basis, pbc = self.pbc)
        return select, atoms
        
    # def autobond(self):
    
    # def at_group(self):
        
def cart2xys(pos, cell):
    """
    Cartesian (absolute) position in angstrom to fractional position (scaled position in lattice).
    """
    pos = np.asarray(pos)
    bg = np.linalg.inv(cell.lattice)
    xyzs = np.tensordot(bg, pos.T, axes=([-1], 0)).T
    return xyzs

def xys2cart(pos, cell):
    """
    Fractional position (scaled position in lattice) to cartesian (absolute) position in angstrom.
    """
    pos = np.asarray(pos)
    xyzr = np.tensordot(cell.lattice, pos.T, axes=([-1], 0)).T
    return xyzr

def translate_pretty(fractional, pbc):
    """Translates atoms such that fractional positions are minimized."""

    for i in range(3):
        if not pbc[i]:
            continue

        indices = np.argsort(fractional[:, i])
        sp = fractional[indices, i]

        widths = (np.roll(sp, 1) - sp) % 1.0
        fractional[:, i] -= sp[np.argmin(widths)]
        fractional[:, i] %= 1.0
    return fractional
