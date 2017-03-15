if (window.File && window.FileReader && window.FileList && window.Blob) {
	console.log("Great Success!")
} else {
	alert('The File APIs are not fully supported in this browser')
}

var netstr = "";

Array.prototype.last = function() {
	return this[this.length-1];
};

Array.prototype.empty = function() {
	return this.length == 0;
};

Array.prototype.contains = function(i) {
	for (n in this) {
		if (i == n) {
			return true;
		}
	}
	return false;
};

function minimize_obj (obj) {
	_min  = [Object.keys(obj)[0], obj[Object.keys(obj)[0]]];
	for (k in this) {
		if (obj[k] < _min[1]) {
			_min = [k, obj[k]];
		}
	}
	return _min[0];
};


class Node {
	
	constructor(name, _id, mu, hierarchy, pos, hyperstring) {
		this.name = name;
		this._id = _id;
		this.mu = mu;
		this.hierarchy = hierarchy;
		this.pos = pos;
		this.hyperstring = hyperstring;
	}
}

class Edge { 

	constructor(source, target, layer) {
		this.source = source;
		this.target = target;
		this.layer = layer;
	}
}

class Network {
	
	constructor() {
		this.edges = [];
		this.nodes = [];
		this.groups = [];
	}

}

Network.prototype.add_edge = function(e) {
	for (ed in this.edges) {
		if (e.target == ed.target && e.source == ed.target) { return; }
		if (e.source == ed.target && e.target == ed.source) { return; }	
	}
	this.edges.push(e);
};

Network.prototype.add_node = function(n) {
	for (no in this.nodes) {
		if (n.name == no.name) { return; }
	}
	this.nodes.push(n);
}

Network.prototype.check_for_nodes = function(name) {
	if (this.nodes.length == 0) {
		return true;
	}
	for (n in this.nodes) {
		if (name == n.name) { return false; }
	}
	return true;
};

Network.prototype.add_to_hyperstring = function (name, val) {
	for (n in this.nodes) {
		if (name == n.name) {
			n.hyperstring += val;
			return;
		}
	}
};

Network.prototype.dijkstras = function(group) {
	traversal = [group[0]];
	while (traversal.length != group.length) {
		v = traversal.last();
		distances = {};
		for (i = 0; i < group.length; i++) {
			var n = group[i];
			if (n !=v && !traversal.contains(n)) { distances[n] = this.compute_distance(v, n); } 
		}
		_next = minimize_obj(distances);
		if (!traversal.contains(_next)) { traversal.push(_next); }
	}
	return traversal;
};

Network.prototype.compute_distance = function(v, t) {
	visited = [];
	neigh = this.neighbors(v, visited);
	dist = 0;
	while (!neigh.empty()) {
		if (neigh.contains(t)) { return dist; }
		new_neigh = [];
		dist += 1;
		for (n in neigh) {
			visited.append(n);
			for (k in this.neighbors(n, visited)) {	new_neigh.push(k); }
		}
		neigh = new_neigh;
	}
};

Network.prototype.neighbors = function(v, visited) {
	adj = [];
	for (e in this.edges) {
		if (adj.contains(e.target) || adj.contains(e.source)) { continue; }
		if (e.target == v && !visited.contains(e.source)) {	adj.push(e.source); }
		if (e.source == v && !visited.contains(e.target)) { adj.push(e.target); }
	}
	return adj;
};
	
Network.prototype.find_id = function(name) {
	_id = 0;
	this.nodes.forEach (function(n) {
		if (n.name == name) { return n._id; }
		if (n._id >= _id) { _id = n._id + 1; }
	});
	return _id;
};

Network.prototype.count_layers = function() {
	curr_layer = 0;
	for (n in this.nodes) { 
		if (n.hyperstring > curr_layer) { curr_layer = n.hyperstring; }
	}
			
	num = 0;
	while (curr_layer != 0) {
		curr_layer = curr_layer >> 1;
		num++;
	}
	return num;
};

function add_to_network(f) {
	var n = new Network();
	var exists = false;

	function check_exists(url) {
		var http = new XMLHttpRequest();
		http.open('HEAD', url, false);
		http.send();
		return http.status != 404;
	}

	exists = check_exists("/data/network.json");
	if (exists) {
		d3.json("/data/network.json", function(d) {
			var nodes = d['nodes'];
			var links = d['links'];
			for (no in nodes) {
				var nn = new Node(no['name'], no['id'], no['mu'], no['hierarchy'], no['position'], no['hyperstring']);
				n.add_node(nn);
			}
			for (l in links) {
				var nl = new Edge(l['source'], l['target'], l['layer']);
				n.add_edge(nl);
			}
		});
	}
	
	var net_url = "/data/network.json";
	
	var bits = {};
	var bit_ind = n.count_layers();
	var curr_net = [];
	var data = JSON.parse(f);
	data.forEach(function(no) {
		var name = no['Name'];
		var _id = n.find_id(name);
		var mu = no["MU"];
		var hierarchy = no["Hierarchy"];
		var pos = no["Position"];
		if (n.check_for_nodes(name)) {
			var nn = new Node(name, _id, mu, hierarchy, pos, Math.pow(2, bit_ind));
			n.add_node(nn);
		} else {
			n.add_to_hyperstring(name, Math.pow(2, bit_ind));
		}
		curr_net.push(_id);
	});
	
	netstr += '{\n\t"nodes": [';
	n.nodes.forEach(function(node) {
		netstr += '\n\t\t{\n';
		netstr += '\t\t\t"id": ' + node._id.toString() + ",\n";
		netstr += '\t\t\t"position": "' + node.pos.toString() + '",\n';
		netstr += '\t\t\t"mu": "' + node.mu.toString() + '",\n';
		netstr += '\t\t\t"hierarchy": ' + node.hierarchy.toString() + ',\n';
		netstr += '\t\t\t"name": "' + node.name.toString() + '",\n';
		netstr += '\t\t\t"hyperstring": ' + node.hyperstring.toString() + '\n';
		netstr += '\t\t},';
	});

	netstr = netstr.slice(0, -1);
	var traversal = n.dijkstras(curr_net);
	var ind = 0;
	while (ind < traversal.length) {
		if (ind < traversal.length - 1) {
			var ne = new Edge(traversal[ind], traversal[ind+1], Math.pow(2, bit_ind));
		} else {
			var ne = new Edge(traversal[ind], traversal[0], Math.pow(2, bit_ind));
		}
		n.add_edge(ne);
		ind++;
	}
	
	netstr += '\n\t],\n\t"links": [';
	n.edges.forEach(function(e) {
		var _s = e.source;
		var _t = e.target;
		netstr += "\n\t\t{\n";
		netstr += '\t\t\t"source": ' + _s.toString() + ",\n";
		netstr += '\t\t\t"target": ' + _t.toString() + ",\n";
		netstr += '\t\t\t"value": 1,\n';
		netstr += '\t\t\t"layer": ' + e.layer.toString() + ",\n";
		netstr += '\t\t\t"sim": 0\n';
		netstr += '\t\t},';
	});
		
	netstr = netstr.slice(0, -1);
	netstr += '\n\t]\n}';
}


function ExceltoJson(wb) {
	var result = "";
	wb.SheetNames.forEach(function(sheetname) {
		var roa = XLSX.utils.sheet_to_row_object_array(wb.Sheets[sheetname]);
		if (roa.length > 0) { result = roa; }
	});
	return result;
}


function handleFileSelect(evt) {
	var f = evt.target.files[0];
	var reader = new FileReader();
	var deferred = $.Deferred();
	reader.onload = (function(e) {
		var d = e.target.result;
		var data = new Uint8Array(d);
		var arr = new Array();
		for (i = 0; i != data.length; ++i) { arr[i] = String.fromCharCode(data[i]); }
		var bstr = arr.join("");
		var wb = XLSX.read(bstr, {type: "binary"});
		var	json = JSON.stringify(ExceltoJson(wb), 2, 2);
		add_to_network(json);
	});
	reader.readAsArrayBuffer(f);
	reader.onloadend = function(e) {
		document.getElementById("vis").data  = JSON.parse(netstr);
	};
	return f.name;
}
