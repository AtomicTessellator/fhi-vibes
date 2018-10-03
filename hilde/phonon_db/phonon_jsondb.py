from __future__ import absolute_import, print_function
import os
import sys

import numpy as np

from ase.db.core import ops, lock, now
from ase.db.jsondb import JSONDatabase
from ase.io.jsonio import encode, decode
from ase.parallel import world, parallel_function
from ase.utils import basestring

from hilde.phonon_db.phonon_db import PhononDatabase
from hilde.phonon_db.row import PhononRow


class PhononJSONDatabase(PhononDatabase, JSONDatabase, object):
    def _write(self, phonon, key_value_pairs, data, id):
        PhononDatabase._write(self, phonon, key_value_pairs, data)
        bigdct = {}
        ids = []
        nextid = 1
        if (isinstance(self.filename, basestring) and
            os.path.isfile(self.filename)):
            try:
                bigdct, ids, nextid = self._read_json()
            except (SyntaxError, ValueError):
                pass
        mtime = now()
        if isinstance(phonon, PhononRow):
            row = phonon
        else:
            row = PhononRow(phonon)
            row.ctime = mtime
            row.user = os.getenv('USER')
        dct = {}
        for key in row.__dict__:
            if key[0] == '_' or key in row._keys or key == 'id':
                continue
            dct[key] = row[key]
        dct['mtime'] = mtime
        if key_value_pairs:
            dct['key_value_pairs'] = key_value_pairs
        if data:
            dct['data'] = data
        constraints = row.get('constraints')
        if constraints:
            dct['constraints'] = constraints
        if id is None:
            id = nextid
            ids.append(id)
            nextid += 1
        else:
            assert id in bigdct
        bigdct[id] = dct
        self._write_json(bigdct, ids, nextid)
        return id

    def _get_row(self, id):
        bigdct, ids, nextid = self._read_json()
        if id is None:
            assert len(ids) == 1
            id = ids[0]
        dct = bigdct[id]
        dct['id'] = id
        return PhononRow(dct)

    def _select(self, keys, cmps, explain=False, verbosity=0, limit=None, offset=0, sort=None, include_data=True, columns='all'):
        if explain:
            yield {'explain': (0, 0, 0, 'scan table')}
            return

        if sort:
            if sort[0] == '-':
                reverse = True
                sort = sort[1:]
            else:
                reverse = False

            def f(row):
                return row.get(sort, missing)

            rows = []
            missing = []
            for row in self._select(keys, cmps):
                key = row.get(sort)
                if key is None:
                    missing.append((0, row))
                else:
                    rows.append((key, row))

            rows.sort(reverse=reverse, key=lambda x: x[0])
            rows += missing

            if limit:
                rows = rows[offset:offset + limit]
            for key, row in rows:
                yield row
            return

        try:
            bigdct, ids, nextid = self._read_json()
        except IOError:
            return

        if not limit:
            limit = -offset - 1

        cmps = [(key, ops[op], val) for key, op, val in cmps]
        n = 0
        for id in ids:
            if n - offset == limit:
                return
            dct = bigdct[id]
            if not include_data:
                dct.pop('data', None)
            row = PhononRow(dct)
            row.id = id
            for key in keys:
                if key not in row:
                    break
            else:
                for key, op, val in cmps:
                    if(key == 'supercell_matrix'):
                        val = list(val.flatten())
                    if isinstance(key, int):
                        value = np.equal(row.numbers, key).sum()
                    else:
                        value = row.get(key)
                        if key == 'pbc':
                            assert op in [ops['='], ops['!=']]
                            value = ''.join('FT'[x] for x in value)
                    if value is None or not op(value, val):
                        break
                else:
                    if n >= offset:
                        yield row
                    n += 1