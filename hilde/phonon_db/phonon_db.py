import os
import numpy as np
import operator
import os
import re
import warnings
import collections
import functools
import numbers
# Import ase
from ase.utils import Lock, basestring, PurePath
from ase.db.core import Database, now, lock, parse_selection, check, str_represents
from ase.calculators.calculator import all_properties, all_changes
from ase.parallel import world, DummyMPI, parallel_function, parallel_generator
from ase.data import atomic_numbers
from ase.atoms import Atoms, symbols2numbers, string2symbols

# Import Hilde
from hilde.phonon_db.row import PhononRow
from hilde.materials_fp.material_fingerprint import get_phonon_bs_fingerprint_phononpy, get_phonon_dos_fingerprint_phononpy
from hilde.structure.structure import pAtoms

# File largely copied from ase.db.core modified to use PhononRows over AtomsRow
T2000 = 946681200.0  # January 1. 2000
YEAR = 31557600.0  # 365.25 days

default_key_descriptions = {
    'id': ('ID', 'Uniqe row ID', ''),
    'age': ('Age', 'Time since creation', ''),
    'formula': ('Formula', 'Chemical formula', ''),
    'user': ('Username', '', ''),
    'calculator': ('Calculator', 'ASE-calculator name', ''),
    'energy': ('Energy', 'Total energy', 'eV'),
    'fmax': ('Maximum force', '', 'eV/Ang'),
    'smax': ('Maximum stress', '', '`\\text{eV/Ang}^3`'),
    'pbc': ('PBC', 'Periodic boundary conditions', ''),
    'charge': ('Charge', '', '|e|'),
    'mass': ('Mass', '', 'au'),
    'magmom': ('Magnetic moment', '', 'au'),
    'unique_id': ('Unique ID', '', ''),
    'volume': ('Volume', 'Volume of unit-cell', '`\\text{Ang}^3`')}

reserved_keys = set(all_properties +
                    all_changes +
                    list(atomic_numbers) +
                    ['id', 'unique_id', 'ctime', 'mtime', 'user',
                     'momenta', 'constraints', 'natoms', 'formula', 'age',
                     'calculator', 'calculator_parameters',
                     'key_value_pairs', 'data'])


seconds = {'s': 1,
           'm': 60,
           'h': 3600,
           'd': 86400,
           'w': 604800,
           'M': 2629800,
           'y': YEAR}

longwords = {'s': 'second',
             'm': 'minute',
             'h': 'hour',
             'd': 'day',
             'w': 'week',
             'M': 'month',
             'y': 'year'}

ops = {'<': operator.lt,
       '<=': operator.le,
       '=': operator.eq,
       '>=': operator.ge,
       '>': operator.gt,
       '!=': operator.ne}

invop = {'<': '>=', '<=': '>', '>=': '<', '>': '<=', '=': '!=', '!=': '='}

word = re.compile('[_a-zA-Z][_0-9a-zA-Z]*$')

reserved_keys = set(all_properties +
                    all_changes +
                    list(atomic_numbers) +
                    ['id', 'unique_id', 'ctime', 'mtime', 'user',
                     'momenta', 'constraints', 'natoms', 'formula', 'age',
                     'calculator', 'calculator_parameters',
                     'key_value_pairs', 'data'])

numeric_keys = set(['id', 'energy', 'magmom', 'charge', 'natoms', 'natoms_in_Sc'])
numeric_keys = set(['id', 'energy', 'magmom', 'charge', 'natoms', 'natoms_in_Sc'])

def str_represents(value, t=int):
    try:
        t(value)
    except ValueError:
        return False
    return True

def check(key_value_pairs):
    for key, value in key_value_pairs.items():
        if not word.match(key) or key in reserved_keys:
            raise ValueError('Bad key: {}'.format(key))
        try:
            string2symbols(key)
        except ValueError:
            pass
        else:
            warnings.warn(
                'It is best not to use keys ({0}) that are also a '
                'chemical formula.  If you do a "db.select({0!r})",'
                'you will not find rows with your key.  Instead, you wil get '
                'rows containing the atoms in the formula!'.format(key))
        if not isinstance(value, (numbers.Real, basestring, np.bool_)):
            print(value)
            raise ValueError('Bad value for {!r}: {}'.format(key, value))
        if isinstance(value, basestring):
            for t in [int, float]:
                if str_represents(value, t):
                    raise ValueError(
                        'Value ' + value + ' is put in as string ' +
                        'but can be interpreted as ' +
                        '{}! Please convert '.format(t.__name__) +
                        'to {} using '.format(t.__name__) +
                        '{}(value) before '.format(t.__name__) +
                        'writing to the database OR change ' +
                        'to a different string.')

def connect(name, type='extract_from_name', create_indices=True, use_lock_file=True, append=True, serial=False):
    """Create connection to database.
    Modified to link to PhononDatabase types
    name: str
        Filename or address of database.
    type: str
        One of 'json', 'db', 'postgresql',
        (JSON, SQLite, PostgreSQL).
        Default is 'extract_from_name', which will guess the type
        from the name.
    use_lock_file: bool
        You can turn this off if you know what you are doing ...
    append: bool
        Use append=False to start a new database.
    """

    if type == 'extract_from_name':
        if name is None:
            type = None
        elif not isinstance(name, basestring):
            type = 'json'
        elif (name.startswith('postgresql:/') or
              name.startswith('postgres:/')):
            type = 'postgresql'
        else:
            type = os.path.splitext(name)[1][1:]
            if type == '':
                raise ValueError('No file extension or database type given')

    if type is None:
        return Database()

    if not append and world.rank == 0 and os.path.isfile(name):
        os.remove(name)

    if isinstance(name, PurePath):
        name = str(name)

    if type != 'postgresql' and isinstance(name, basestring):
        name = os.path.abspath(name)

    if type == 'json':
        from hilde.phonon_db.phonon_jsondb import PhononJSONDatabase
        return PhononJSONDatabase(name, use_lock_file=use_lock_file, serial=serial)
    if type == 'db':
        from hilde.phonon_db.phonon_sqlitedb import PhononSQLite3Database
        return PhononSQLite3Database(name, create_indices, use_lock_file,
                               serial=serial)
    if type == 'postgresql':
        from hilde.phonon_db.phonon_postgresqldb import PhononPostgreSQLDatabase
        return PhononPostgreSQLDatabase(name)
    raise ValueError('Unknown database type: ' + type)


class PhononDatabase(Database):
    """Base class for all databases."""
    def __init__(self, filename=None, create_indices=True, use_lock_file=False, serial=False):
        """Database object.
        Args:
            filename: str
                Filename of the database
            create_indices: bool

            serial: bool
                Let someone else handle parallelization.  Default behavior is
                to interact with the database on the master only and then
                distribute results to all slaves.
        """
        if isinstance(filename, basestring):
            filename = os.path.expanduser(filename)
        self.filename = filename
        self.create_indices = create_indices
        if use_lock_file and isinstance(filename, basestring):
            self.lock = Lock(filename + '.lock', world=DummyMPI())
        else:
            self.lock = None
        self.serial = serial
        self._metadata = None  # decription of columns and other stuff

    # def _write(self, atoms, key_value_pairs, data, id=None):
    #     check(key_value_pairs)
    #     return 1

    def get_phonon(self, selection=None, **kwarg):
        '''
        Gets a phonopy object from a database row
        Args:
            selection: selection criteria for the database query
            kwargs   : additional selection criteria not stored in selection
        Returns:
            the phonopy object of the row
        '''
        row = self.get(selection, **kwarg)
        return row.to_phonon()

    @parallel_function
    @lock
    def update(self, id, phonon=None, delete_keys=[], data=None, **add_key_value_pairs):
        """Update and/or delete key-value pairs of row(s).

        id: int
            ID of row to update.
        phonon: Phonopy object
            Optionally update the Phononpy data (positions, cell, ...).
        data: dict
            Data dict to be added to the existing data.
        delete_keys: list of str
            Keys to remove.

        Use keyword arguments to add new key-value pairs.

        Returns number of key-value pairs added and removed.
        """
        if not isinstance(id, numbers.Integral):
            if isinstance(id, list):
                err = ('First argument must be an int and not a list.\n'
                       'Do something like this instead:\n\n'
                       'with db:\n'
                       '    for id in ids:\n'
                       '        db.update(id, ...)')
                raise ValueError(err)
            raise TypeError('id must be an int')

        check(add_key_value_pairs)

        row = self._get_row(id)

        if phonon:
            oldrow = row
            row = PhononRow(phonon)

            # Copy over data, kvp, ctime, user and id
            row._data = oldrow._data
            kvp = oldrow.key_value_pairs
            row.__dict__.update(kvp)
            row._keys = list(kvp)
            row.ctime = oldrow.ctime
            row.user = oldrow.user
            row.id = id

        kvp = row.key_value_pairs

        n = len(kvp)
        for key in delete_keys:
            kvp.pop(key, None)
        n -= len(kvp)
        m = -len(kvp)
        kvp.update(add_key_value_pairs)
        m += len(kvp)

        moredata = data
        data = row.get('data', {})
        if moredata:
            data.update(moredata)
        if not data:
            data = None

        self._write(row, kvp, data, row.id)

        return m, n

    @parallel_generator
    def select(self, selection=None, filter=None, explain=False, verbosity=1, limit=None, offset=0, sort=None, include_data=True, columns='all', **kwargs):
        """Select rows.

        Return PhononRow iterator with results.  Selection is done
        using key-value pairs and the special keys:

            formula, age, user, calculator, natoms, energy, magmom
            and/or charge.

        selection: int, str or list
            Can be:

            * an integer id
            * a string like 'key=value', where '=' can also be one of
              '<=', '<', '>', '>=' or '!='.
            * a string like 'key'
            * comma separated strings like 'key1<value1,key2=value2,key'
            * list of strings or tuples: [('charge', '=', 1)].
        filter: function
            A function that takes as input a row and returns True or False.
        explain: bool
            Explain query plan.
        verbosity: int
            Possible values: 0, 1 or 2.
        limit: int or None
            Limit selection.
        offset: int
            Offset into selected rows.
        sort: str
            Sort rows after key.  Prepend with minus sign for a decending sort.
        include_data: bool
            Use include_data=False to skip reading data from rows.
        columns: 'all' or list of str
            Specify which columns from the SQL table to include.
            For example, if only the row id and the energy is needed,
            queries can be speeded up by setting columns=['id', 'energy'].
        """
        if sort:
            if sort == 'age':
                sort = '-ctime'
            elif sort == '-age':
                sort = 'ctime'
            elif sort.lstrip('-') == 'user':
                sort += 'name'
        keys, cmps = parse_selection(selection, **kwargs)
        for row in self._select(keys, cmps, explain=explain,
                                verbosity=verbosity,
                                limit=limit, offset=offset, sort=sort,
                                include_data=include_data,
                                columns=columns):
            if filter is None or filter(row):
                yield row