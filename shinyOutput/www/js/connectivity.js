function compute_connectivity(data) {
	data.nodes.forEach(function(d, i) {
		d.connectivity = 0;
		data.links.forEach(function(l, v) {
			if (l.source.name == d.name || l.target.name == d.name) {
				d.connectivity += 1;
			}
		});
	});
}
