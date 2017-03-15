#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from openpyxl import load_workbook, worksheet
import sys, os
import cgitb
import json
import cherrypy
import getopt

class _list(list):
		
	def last(self):
		return self[len(self)-1]

	def empty(self):
		return len(self) == 0


class _dict(dict):
		
	def minimize(self):
		_min = (self.keys()[0], self[self.keys()[0]])
		for k in self.keys():
			if self[k] < _min[1]:
				_min = (k, self[k])

		return _min[0]
				

class Insight:
		
	def __init__(self, name, _id, date, insight, source):
		self.name = name
		self.insight = insight
		self._id = _id
		self.date = date
		self.source = source

	def print_insight(self):
		print ("NAME:" + str(self.name) + "\n" + "INSIGHT:" + str(self.insight) + "\n" + "ID:" + str(self._id) + "\n" + "DATE:" + str(self.date) + "\n" + "SOURCE:" + str(self.source) + "\n")

class Process:
	
	def __init__(self, name, _id):
		self.name = name
		self._id = _id
		self.seq = {}
		self.descrip = {}

	def add_name_to_step(self, step, name):
		if step in self.seq:
			self.seq[step].append(name)
		else:
			self.seq[step] = [name]

	def add_des_to_step(self, step, des):
		self.descrip[step] = des

	def add_seq_arr(self, ind, arr):
		self.seq[ind] = arr

	def add_des_arr(self, ind, arr):
		self.descrip[ind] = arr



class Node:
	
	def __init__(self, name, _id, mu, hierarchy, pos, hyperstring):
		self.name = name
		self.hierarchy = hierarchy
		self.mu = mu
		self._id = _id
		self.hyperstring = hyperstring
		self.pos = pos

class Edge:
	
	def __init__(self, source, target, layer, sim):
		self.source = source
		self.target = target
		self.layer = layer
		self.sim = sim

class Group:

	def __init__(self, name, mask):
		self.name = name
		self.mask = mask
		self.fte = 0
		self.budget = 0

	def define_fte(self, val):
		self.fte = val

	def define_budget(self, val):
		self.budget = val

class Network:
	
	def __init__(self):
		self.edges = []
		self.nodes = []
		self.groups = []
		self.processes = []
		self.insights = []

	
	def add_edge(self, e):
		for ed in self.edges:
			if e.target == ed.target and e.source == ed.target:
				return 
			if e.source == ed.target and e.target == ed.source:
				return
		self.edges.append(e)

	def add_node(self, n):
		for no in self.nodes:
			if n._id == no._id:
				return
		self.nodes.append(n)

	def add_insight(self, i):
		self.insights.append(i)

	def check_for_nodes(self, _id):
		if len(self.nodes) == 0:
			return True
		for n in self.nodes:
			if _id == n._id:
				return False
		return True

	def add_to_hyperstring(self, name, val):
		for n  in self.nodes:
			if name == n.name:
				n.hyperstring += val
				return

	def compute_id(self, name, pos):
		_id = 0
		for n in self.nodes:
			if n.pos == pos and n.name == name:
				return n._id
			if n._id >= _id:
				_id = n._id + 1
		return _id

	def find_id(self, name):
		for n in self.nodes:
			if n.name == name:
				return n._id

	def find_id_by_pos(self, pos):
		for n in self.nodes:
			if n.pos == pos:
				return n._id

	def count_layers(self):
		curr_layer = 0
		for n in self.nodes:
			if n.hyperstring > curr_layer:
				curr_layer = n.hyperstring

		num = 0
		while curr_layer != 0:
			curr_layer = curr_layer >> 1
			num +=1 
		return num
	
	def add_group(self, group):
		for g in self.groups:
			if group.name == g.name:
				return
		else:
			self.groups.append(group)
		
	def find_next_mask(self, name):
		b = 1
		for k in self.groups:
			if name == k.name:
				return k.mask
			if k.mask >= b:
				b = k.mask * 2
		return b

	def add_process(self, p):
		self.processes.append(p)
	
	
#@app.route("/main")
def main(argv):
	n = Network()
	g_active = False
	p_active = False
	f_active = False
	b_active = False
	i_active = False
	input_g = ""
	input_p = ""
	input_b = ""
	input_i = ""

	_nextG = False
	_nextP = False
	_nextF = False	
	_nextB = False
	_nextI = False

	for a in argv:
		if a == '-g':
			_nextG = True
			continue
		elif a == "-f":
			_nextF = True
			continue
		elif a == '-p':
			_nextP = True
			continue
		elif a == '-b':
			_nextB = True
		elif a == '-i':
			_nextI = True
			continue
		elif _nextG:
			input_g = a
			g_active = True
		elif _nextF:
			input_g = a
			f_active = True
		elif _nextP:
			input_p = a
			p_active = True
		elif _nextB:
			input_r = a
			b_active = True
		elif _nextI:
			input_i = a
			i_active = True	


	# If there's already a data file, add nodes/links to network
	if os.path.isfile("./data/network.json"):
		data_file = open("./data/network.json")
		d = json.load(data_file)
		nodes = d['nodes']
		links = d['links']
		groups = d['groups']
		ps = d['processes']
		insights = d['insights']
		for no in nodes:
			nn = Node(no['name'], no['id'], no['mu'], no['hierarchy'], no['position'], no['hyperstring'])
			n.add_node(nn)
		for l in links:
			nl = Edge(l['source'], l['target'], l['layer'], l['sim'])
			n.add_edge(nl)
		for g in groups:
			ng = Group(g['name'], g['mask'])
			ng.define_fte(g['fte'])
			ng.define_budget(g['budget'])
			n.add_group(ng)
		for p in ps:
			np = Process(p['name'], p['id'])
			for k in p:
				if k == "name" or k == "id":
					continue
				elif 'seq' in k:
					np.add_seq_arr(k, p[k])
				elif 'des' in k:
					np.add_des_arr(k, p[k])
			n.add_sim(np)
	
		for i in insights:
			ni = Insight(i['name'], i['id'], i['date'], i['insight'], i['source'])
			n.add_insight(ni)

	net = open("./data/network.json", 'w')
	if g_active:
		wb = load_workbook(input_g, data_only=True, use_iterators=True)
		n_ws = wb.get_sheet_by_name("nodes")
		e_ws = wb.get_sheet_by_name("edges")

		bits = {}
		bit_ind = n.count_layers()
		print bit_ind
		curr_sub = _list([])
	
		for row in n_ws.iter_rows():
			if row[0].value == 'Name':
				continue
			
			# Find node information	
			name = row[0].value
			pos = row[1].value
			_id = n.compute_id(name, pos)
			if n.check_for_nodes(_id):
				pos = row[1].value
				mu = row[2].value
				hierarchy = row[3].value

				# Add node to network
				nn = Node(name, _id, mu, hierarchy, pos, pow(2, bit_ind))
				n.add_node(nn)
			else:
				n.add_to_hyperstring(name, pow(2, bit_ind))
			curr_sub.append(_id)
		#Don't add subnetwork links, makes the graph too convoluted
		
		ng = input_g.split(".xlsx")
		name = ng[0].split("data/")[1]
		g = Group(name, n.find_next_mask(name))
		n.add_group(g)

	if f_active:
		wb = load_workbook(input_g, data_only=True, use_iterators=True)
		n_ws = wb.get_sheet_by_name("nodes")
		e_ws = wb.get_sheet_by_name("edges")

		curr_sub = _list([])
		for row in n_ws.iter_rows():
			if row[0].value == 'Name':
				continue
			
			# Find node information	
			name = row[0].value
			pos = row[1].value
			_id = row[4].value
			if n.check_for_nodes(_id):
				mu = row[2].value
				hierarchy = row[3].value

				# Add node to network
				nn = Node(name, _id, mu, hierarchy, pos, 0)
				n.add_node(nn)
			else:
				n.add_to_hyperstring(name, 0)
			curr_sub.append(_id)
	
		# If foundation file, construct edges by id's	
		for row in e_ws.iter_rows():
			if row[0].value == "Source":
				continue
	
			source = row[2].value
			target = row[3].value
			
			
			ne = Edge(source, target, 0, 0)
			n.add_edge(ne)
		
		ng = input_g.split(".xlsx")
		name = ng[0].split("data/")[1]
		g = Group(name, 0)
		n.add_group(g)
															
	
	if p_active:
		wb = load_workbook(input_p, data_only=True, use_iterators=True)
		p_ws = wb.get_sheet_by_name("process")
		d_ws = wb.get_sheet_by_name("description")
	
		f = input_p.split(".xlsx")
		name = f[0].split("data/")[1]

		np = Process(name.lower(), len(n.processes) + 1)
		for row in p_ws.iter_rows():
			if type(row[0].value) == unicode:
				continue
			
			# find person's id & add to proper sim dict
			_id = int(row[2].value)
			_step = int(row[0].value)
			np.add_name_to_step(_step, _id)

		for row in d_ws.iter_rows():
			if type(row[0].value) == unicode:
				continue
			
			_step = int(row[0].value)
			des = row[1].value
			np.add_des_to_step(_step, des)

		n.add_process(np)

		# Now add the simulation's edges 
		#i = 0
		#while i < len(ns.seq):
		#	if i < len(ns.seq) - 1:
		#		_source = ns.seq[i]
		#		_target = ns.seq[i+1]
		#	else:
		#		_source = ns.seq[i]
		#		_target = ns.seq[0]
		#	ne = Edge(_source, _target, 0, ns._id)
		#	n.add_edge(ne)
		#	i += 1 

	if b_active:
		wb = load_workbook(input_r, data_only=True, use_iterators=True)
		r_ws = wb.get_sheet_by_name("resources")
		
		for row in r_ws.iter_rows():
			if row[0].value.lower() == 'group':
				continue
			
			name = row[0].value
			found = False
			for g in n.groups:
				if g.name.lower() == name.lower():
					g.define_fte(row[2].value)
					g.define_budget(row[1].value)
					found = True
			if not found:
				ng = Group(name, n.find_next_mask(name))
				ng.define_fte(row[2].value)
				ng.define_budget(row[1].value)
				n.add_group(ng)
				
	if i_active:
		wb = load_workbook(input_i, data_only=True, use_iterators=True)
		i_ws = wb.get_sheet_by_name("insights")

		for row in i_ws.iter_rows():
			
			if row[0].value.lower() == "name":
				continue

			name = row[0].value
			pos = row[1].value
			date = row[2].value
			insight = row[3].value
			_id = n.compute_id(name, pos)
			source = row[4].value

			ni = Insight(name, _id, date, insight, source)
			n.add_insight(ni)

	
	# Write nodes to file	
	net.write('{\n\t"nodes": [')
	for node in n.nodes:
		net.write("\n\t\t{\n")
		net.write('\t\t\t"id": ' + str(node._id) + ",\n")
		net.write('\t\t\t"position": "' + str(node.pos) + '",\n')
		net.write('\t\t\t"mu": "' + str(node.mu) + '",\n')
		net.write('\t\t\t"hierarchy": ' + str(node.hierarchy) + ',\n')
		net.write('\t\t\t"name" : "' + str(node.name) + '",\n')
		net.write('\t\t\t"hyperstring": ' + str(node.hyperstring) + '\n')
		net.write('\t\t},')

	if len(n.nodes) > 0:
		net.seek(-1, os.SEEK_CUR)
	
	
	net.write('\n\t],\n\t"links": [')	
	for link in n.edges:
		_source = link.source
		_tar = link.target
		_sim = link.sim
		net.write("\n\t\t{\n")
		net.write('\t\t\t"source": ' + str(_source) + ",\n")
		net.write('\t\t\t"target": ' + str(_tar) + ",\n")
		net.write('\t\t\t"value": ' + str(1) + ",\n")
		net.write('\t\t\t"layer": ' + str(link.layer) + ",\n")
		net.write('\t\t\t"sim": ' + str(_sim) + '\n')
		net.write('\t\t},')

	if len(n.edges) > 0:
		net.seek(-1, os.SEEK_CUR)


	net.write('\n\t],\n\t"groups": [')
	for g in n.groups:
		_name = g.name
		_mask = g.mask
		_fte = g.fte
		_budget = g.budget
		net.write("\n\t\t{\n")
		net.write('\t\t\t"name": "' + str(_name) + '",\n')
		net.write('\t\t\t"mask": ' + str(_mask) + ',\n')
		net.write('\t\t\t"fte": ' + str(_fte) + ',\n')
		net.write('\t\t\t"budget": ' + str(_budget) + '\n')
		net.write('\t\t},')

	if len(n.groups) > 0:
		net.seek(-1, os.SEEK_CUR)
	
		
		
	net.write('\n\t],\n\t"processes": [')
	for p in n.processes:
		_name = p.name
		_id = p._id
		net.write("\n\t\t{\n")
		net.write('\t\t\t"name": "' + str(_name) + '",\n')
		net.write('\t\t\t"id": ' + str(_id) + ',\n')
		for k in p.seq:
			name = "seq_" + str(k)
			net.write('\t\t\t"' + name + '": ' + str(p.seq[k]) + ',\n')
		
		i = 1
		for k in p.descrip:
			name = "des_" + str(k)
			if i == len(p.seq):
				net.write('\t\t\t"' +  name + '": "' + str(p.descrip[k]) + '"\n')
			else:
				net.write('\t\t\t"' + name + '": "' + str(p.descrip[k]) + '",\n')
			i += 1
			
		net.write('\t\t},')

	if len(n.processes) > 0:
		net.seek(-1, os.SEEK_CUR)


	net.write('\n\t],\n\t"insights": [')
	for i in n.insights:
		_name = i.name
		_id = i._id
		_date = i.date
		_insight = i.insight
		_source = i.source

		net.write("\n\t\t{\n")
		net.write('\t\t\t"name": "' + str(_name) + '",\n')
		net.write('\t\t\t"id": ' + str(_id) + ',\n')
		net.write('\t\t\t"date": "' + str(_date) + '",\n')
		net.write('\t\t\t"insight": "' + str(_insight) + '",\n')
		net.write('\t\t\t"source": "' + str(_source) + '"\n')
		net.write('\t\t},')
			
	
	if len(n.insights) > 0:	
		net.seek(-1, os.SEEK_CUR)
	
	net.write('\n\t]\n}')

#cherrypy.quickstart(main(sys.argv))
main(sys.argv)
