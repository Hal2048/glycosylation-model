"""CFG glycan name → tree parser.

Direct port of GlycoGlyph `src/modules/glycantoJSON.js` (glycantojson).
Takes an IUPAC-style CFG name and returns the nested dict tree used by the
drawing code.
"""

from __future__ import annotations

import re
from typing import Any, Dict

from .monoslist import MONOS_SET


def glycan_to_json(glycan_name: str) -> Dict[str, Any]:
    """Parse a CFG glycan name into a nested node dict.

    Returns the root node with keys:
        name, linkage, monosaccharide, substituent, linknum, anomer,
        uniqueNodeID, children
    """
    name = glycan_name

    # remove whitespace
    name = re.sub(r' ', '', name)
    # replace sulfated position brackets
    name = re.sub(r'\((?=\d[a-zA-Z]\))', '[', name)
    name = re.sub(r'(\[\d[a-zA-Z])\)', r'\1]', name)
    # α/β ahead of a digit/? + hyphen → a/b
    name = re.sub(r'α(?=(\d|\?)\-)', 'a', name)
    name = re.sub(r'β(?=(\d|\?)\-)', 'b', name)

    # GlcN(Gc) → GlcNGc etc.
    if 'GlcN(Gc)' in name:
        name = name.replace('GlcN(Gc)', 'GlcNGc')
    if 'GalN(Gc)' in name:
        name = name.replace('GalN(Gc)', 'GalNGc')
    if 'ManN(Gc)' in name:
        name = name.replace('ManN(Gc)', 'ManNGc')

    arr = _glystr_to_array(name)
    tree = _parse_tree(arr)
    return tree


def _glystr_to_array(s: str):
    # split keeping separators: ')' '(' '-X' '-<letter>'
    parts = re.split(r'(\)|\(|-[\d,?]|-(?=[a-zA-z]))', s)
    arr = [p for p in parts if p != '']

    newarr = []
    for i in range(len(arr) - 1):
        if re.search(r'\-[\d,?]', arr[i + 1]):
            newarr.append(arr[i] + arr[i + 1])
        elif arr[i].find('-') < 0:
            newarr.append(arr[i])
    newarr.append(arr[-1])
    return newarr


def _parse_tree(arr):
    result = []
    for idx, ch in enumerate(arr):
        result = _parse_tree_reducer(result, ch, idx, arr)
    return result[0]


def _parse_tree_reducer(array, ch, idx, src):
    # find last unmatched '(' index
    index_bracket = -1
    for i in range(len(array) - 1, -1, -1):
        if array[i] == '(':
            index_bracket = i
            break

    if ch == '(':
        array.append(ch)
    elif ch == ')':
        if index_bracket >= 0:
            array.pop(index_bracket)
    elif ch != '-':
        link = None
        mono = None
        sub = None
        linknum = None

        last_close = ch.rfind(']')
        chnosub = ch[last_close + 1:] if last_close >= 0 else ch
        if chnosub != ch:
            sub = ch[:len(ch) - len(chnosub)]

        for i in range(len(chnosub)):
            candidate = chnosub[:len(chnosub) - i]
            if candidate in MONOS_SET:
                link = chnosub[len(chnosub) - i:]
                linknum = int(link[-1]) if link and link[-1].isdigit() else None
                mono = candidate
                break

        # special cases
        if re.search(r'Neu5\,9Ac2?', ch):
            mono = 'Neu5Ac'
            sub = '[9Ac]'
            link = re.sub(r'Neu5\,9Ac2?', '', ch)
            linknum = int(link[-1]) if link and link[-1].isdigit() else None

        if re.search(r'(GlcNGc)|(GalNGc)|(ManNGc)', ch):
            mono = mono.replace('Gc', '') if mono else mono
            sub = '[2Gc]'
            link = re.sub(r'(GlcNGc)|(GalNGc)|(ManNGc)', '', ch)
            linknum = int(link[-1]) if link and link[-1].isdigit() else None

        if link is None:
            link = ''
        if mono is None:
            mono = ''
        if sub is None:
            sub = ''
        if linknum is None:
            linknum = 0

        anomer = link[0] if link else ''
        if anomer.isdigit():
            anomer = '?'
        if link.endswith('-'):
            link = link[:-1]

        # collect children = items after the last unmatched '('
        children_source = array[index_bracket + 1:] if index_bracket >= 0 else array[:]
        # remove those items from array
        if index_bracket >= 0:
            del array[index_bracket + 1:]
        else:
            array.clear()
        children = list(reversed(children_source))

        array.append({
            'name': ch,
            'linkage': link,
            'monosaccharide': mono,
            'substituent': sub,
            'linknum': linknum,
            'anomer': anomer,
            'uniqueNodeID': idx,
            'children': children,
        })
    return array
