#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

import requests
from Bio.KEGG.KGML import KGML_parser
from Bio.KEGG.KGML.KGML_pathway import Relation

from .. import config
from ..lib.sutils import *

HOST = 'http://rest.kegg.jp'

DATA_ROOT = config.DATA_DIR / 'kegg'

safe_create_dir(DATA_ROOT)


@timeit
def get_all_pathways(force=False):
    list_url = f'{HOST}/list/pathway/hsa'
    path = DATA_ROOT / 'kegg-all-pathways.txt'
    if force or not path.exists():
        r = requests.get(list_url)
        if not r.status_code == 200:
            print('Failed to get pathway:', path, r.text)
            sys.exit(-1)
        with open(path, 'w') as f:
            f.write(r.text)
            print('Saved KGML file at:', path)

    with open(path) as f:
        lines = [line for line in f.readlines() if line.strip() != '']
        all_pathways = [line.split('\t')[0] for line in lines]
        print('Number of pathways:', len(all_pathways))
        return all_pathways


def get_pathway_kgml(pathway_id='hsa04151'):
    print('Reading pathway:', pathway_id)
    pw_url = f'{HOST}/get/{pathway_id}/kgml'

    # get pathway data if not exists
    path = DATA_ROOT / f'{pathway_id}.kgml'
    if not path.exists():
        print('Could not find hsa file getting from kegg:', pw_url)
        r = requests.get(pw_url)
        if not r.status_code == 200:
            print('Failed to get pathway:', path, r.text)
            sys.exit(-1)
        with open(path, 'w') as f:
            f.write(r.text)
            print('Saved KGML file at:', path)

    # parse data
    with open(path, 'r') as f:
        pathway = KGML_parser.read(f)
    print('\n', pathway, sep='')
    print('  entry:', len(pathway.entries))
    print('  reaction:', len(pathway.reactions))
    print('  relation:', len(pathway.relations))
    entries, relations = prune_kgml(pathway)
    print('Finished reading:', pathway_id)
    print(' entry:', len(entries.keys()), 'relation:', len(relations), 'new_relation:',
          len([1 for r in relations if hasattr(r, '_pamogk')]))
    return entries, relations


def prune_kgml(pathway):
    """
    prune_KGML prune given kgml pathway
      This function prunes the given entries and relations by deleting
      empty entries and corresponding relations, and empty relations and
      corresponding entries. Further pruning is done for our case. We delete
      inner pathways, complex components, PCrel (by inserting edge among
      nodes it connects), maplink, orthogonals
      Inputs:
        pathway:
            entries: index map of entries
            relations: list of relations
      Outputs:
          entries: pruned entries
          relations: pruned relations
          new_relations: new relations added in place of removed compounds
    """

    # copy relations as list and entries as dict
    entries = dict([(e.id, e) for e in pathway.entries.values()])
    relations = [r for r in pathway.relations]
    new_relations = []

    # replace groups with their components and relations between components
    for e in [e for e in entries.values() if e.type == 'group']:
        comp_eids = [comp.id for comp in list(e.components)]
        for e1_id in comp_eids:
            for e2_id in comp_eids:
                if can_add_relation(e1_id, e2_id, relations, new_relations):
                    new_relations.append((e1_id, e2_id))

    # prune orthologs, groups, other, map in entries
    # if any "other" type of entry is included, relations that has
    # these entries are removed from the relations
    type_filter = ['ortholog', 'group', 'other', 'map']
    relations = [r for r in relations if r.entry1.type not in type_filter and r.entry2.type not in type_filter]
    # remove orthologs, groups, and others from entries
    for e in [e for e in entries.values() if e.type in type_filter]:
        del entries[e.id]

    # % prune compounds in entries
    # retrieve gene ids
    gene_ids = set([e.id for e in entries.values() if e.type == 'gene'])

    # for every compound
    for cid in [e.id for e in entries.values() if e.type == 'compound']:
        # find indexes of compounds in relations as source and
        # derive the related ids from it
        source_ids = set([r.entry1.id for r in relations if r.entry2.id == cid])
        source_set = source_ids & gene_ids  # take just genes, ignore others

        # find indexes of compounds in relations as destination
        # and derive the related ids from it
        dest_ids = set([r.entry2.id for r in relations if r.entry1.id == cid])
        dest_set = dest_ids & gene_ids  # take just genes, ignore others

        # create the relation among source and destination if
        # it is not self-relation and if the relation does not already
        # exist in the relations
        # counter = 0
        if not source_ids or not dest_ids:
            continue
        print('src =', source_ids, 'dst =', dest_ids)
        for tmp_s in source_set:
            for tmp_d in dest_set:
                # check if the relation already exists
                if can_add_relation(tmp_s, tmp_d, relations, new_relations):
                    new_relations.append((tmp_s, tmp_d))
                else:
                    print([tmp_s, tmp_d], 'already in relations')
        # removed relations containing this compound entry
        relations = [r for r in relations if r.entry1.id != cid and r.entry2.id != cid]

        # remove compounds from entries
        del entries[cid]

    # turn new relations to relation objects
    new_rel = []
    for s, d in new_relations:
        r = Relation()
        r._entry1 = entries[s]
        r._entry2 = entries[d]
        r._pamogk = True  # to mark relation as added by PAMOGK
        new_rel.append(r)
    relations += new_rel

    return entries, relations


def can_add_relation(s, d, relations, new_relations):
    """
    returns true if (s, d) is a valid new relation
    """
    if s == d:
        return False
    for r in relations:
        if r.entry1.id == s and r.entry2.id == d:
            return False
        if r.entry2.id == s and r.entry1.id == d:
            return False
    for ns, nd in new_relations:
        if ns == s and nd == d:
            return False
        if nd == s and ns == d:
            return False
    return True


if __name__ == '__main__':
    for pw_id in get_all_pathways():
        print(pw_id)
        entry_map, rels = get_pathway_kgml(pw_id)
        print(pw_id, len(entry_map.keys()), len(rels), '\n')
