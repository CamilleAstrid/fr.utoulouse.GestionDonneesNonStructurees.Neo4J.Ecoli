#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pprint import pprint
#import Graph as gr # Graph library from part 1 of the project
import re
import numpy as np
import sys


def create_graph(directed = True, weighted = False, weight_attribute = None): # TP1
	"""
	create a dictionnary representing a graph and returns it.
	"""
	g = { 'nodes': {}, 'edges': {}, 'nb_edges': 0, 'directed': directed, 'weighted': weighted, 'weight_attribute': weight_attribute }
	return g

def add_node(g, n, attributes = None): # TP1
	"""
	add a node n (node id provided as a string or int) to the graph g.
	attributes on the node can be provided by a dict.
	returns the node n attributes.
	"""
	if n not in g['nodes']: # ensure node does not already exist
		if attributes is None: # create empty attributes if not provided
			attributes = {}
		g['nodes'][n] = attributes
		g['edges'][n] = {} # init outgoing edges
	return g['nodes'][n] # return node attributes

def add_edge(g, n1, n2, attributes = None, n1_attributes = None, n2_attributes = None): # TP1
	# create nodes if they do not exist
	if n1 not in g['nodes']: add_node(g, n1, n1_attributes) # ensure n1 exists
	if n2 not in g['nodes']: add_node(g, n2, n2_attributes) # ensure n2 exists
	# add edge(s) only if they do not exist
	if n2 not in g['edges'][n1]:
		if attributes is None: # create empty attributes if not provided
			attributes = {}
		g['edges'][n1][n2] = attributes
		if not g['directed']:
			g['edges'][n2][n1] = g['edges'][n1][n2] # share the same attributes as n1->n2
		g['nb_edges'] += 1
	return g['edges'][n1][n2] # return edge attributes

def load_OBO(filename):
	"""
	parse the OBO file and returns the graph
	obsolete terms are discarded
	only is_a and part_of relationships are loaded

	Extract of a file to be parsed:
	[Term]
	id: GO:0000028
	name: ribosomal small subunit assembly
	namespace: biological_process
	def: "The aggregation, arrangement and bonding together of constituent RNAs and proteins to form the small ribosomal subunit." [GOC:jl]
	subset: gosubset_prok
	synonym: "30S ribosomal subunit assembly" NARROW [GOC:mah]
	synonym: "40S ribosomal subunit assembly" NARROW [GOC:mah]
	is_a: GO:0022618 ! ribonucleoprotein complex assembly
	relationship: part_of GO:0042255 ! ribosome assembly
	relationship: part_of GO:0042274 ! ribosomal small subunit biogenesis
	"""
	def parseTerm(lines):
		# search for obsolete
		for l in lines:
			if l.startswith('is_obsolete: true'): return
		# otherwise create node
		go_id = re_go_id.match(lines.pop(0)).group(1)
		go_node = add_node(g,go_id)
		go_node['id'] = go_id
		go_node['type'] = 'GOTerm'
		for line in lines:
			if re_go_name.match(line): go_node['name'] = re_go_name.match(line).group(1)
			elif re_go_namespace.match(line): go_node['namespace'] = re_go_namespace.match(line).group(1)
			elif re_go_def.match(line): go_node['def'] = re_go_def.match(line).group(1)
			# relationships
			elif re_go_is_a.match(line): 
				parent_id = re_go_is_a.match(line).group(1)
				e = add_edge(g,go_id, parent_id)
				e['type'] = 'is_a'
			elif re_go_part_of.match(line): 
				parent_id = re_go_part_of.match(line).group(1)
				e = add_edge(g, go_id, parent_id)
				e['type'] = 'part_of'
	# instantiate directed graph and additionnal graph attributes
	g=create_graph(directed=True, weighted=False)
	g['alt_id'] = {} # alternate GO ids
	# regexp to parse term lines
	re_go_id          = re.compile(r'^id:\s+(GO:\d+)\s*$')
	re_go_name        = re.compile(r'^name:\s+(.+)\s*$')
	re_go_namespace   = re.compile(r'^namespace:\s+(.+)\s*$')
	re_go_def         = re.compile(r'^def:\s+"(.+)"\s.*$')
	re_go_alt_id      = re.compile(r'^alt_id:\s+(GO:\d+)\s*$')
	re_go_is_a        = re.compile(r'^is_a:\s+(GO:\d+)\s')
	re_go_xref        = re.compile(r'^xref:\s+(\S+)\s*$')
	re_go_part_of      = re.compile(r'^relationship:\s+part_of\s+(GO:\d+)\s')
	with open(filename) as f: 
		line = f.readline().rstrip()
		# skip header
		while not line.startswith('[Term]'): 
			line = f.readline().rstrip()
		buff = []  
		line = f.readline()
		stop = False
		while line and not stop:
			line = line.rstrip()
			# new Term
			if line.startswith('[Term]'):
				parseTerm(buff)
				buff=[]
			# last Term
			elif line.startswith('[Typedef]'):
				parseTerm(buff)
				stop=True
			# or append to buffer
			else:
				buff.append(line)
			line = f.readline()
	return g

if __name__ == "__main__": # argv[1] → go-basic.obo, argv[2] → nodes|edges→is_a/part_of
	go = load_OBO(sys.argv[1])
	#pprint(go)
	if sys.argv[2]=='nodes':
		print("id\tdesc\tdef\tnamespace")
		for n in go['nodes'].values():
			print(f"{n['id']}\t{n['name']}\t{n['def']}\t{n['namespace']}")
	elif sys.argv[2]=='edges':
		print("term1\tterm2\trelationship")
		for t1 in go['edges']:
			for t2 in go['edges'][t1]:
				rel = go['edges'][t1][t2]['type']
				if sys.argv[3]==rel:
					print(f"{t1}\t{t2}\t{rel}")
