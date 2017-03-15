function person() {
	this.name = "";
	this.id = 0;
	this.influence = 0;
	this.group_membership = 0;
	this.connectivity = 0;
	this.hierarchy = 0;
}

function compute_radii(data) {
	var people = [];

	data.nodes.forEach(function(d, i) {
		var np = new person();
		np.name = d.name;
		np.hierarchy = d.hierarchy;
		np.group_membership = 0;
		np.connectivity = 0;
		np.id = d.id;
		var groups = d.hyperstring
		while (groups > 0) {
			if (groups & 1) {
				np.group_membership += 1;
			}
			groups = groups >> 1;
		}
		people.push(np);
	});
	
	for (j=0; j < people.length; j++) {
		var filtered = data.links.filter(function (l) {
			return l.source == people[j].id || l.target == people[j].id;
			});
		people[j].connectivity = filtered.length;
	}

	for (i = 0; i < people.length; i++) {
		
		function compute_inf(p) {
			inf = p.group_membership + p.connectivity - p.hierarchy;
			if (inf > 0) { 
				return parseInt(5 + Math.log10(p.group_membership + p.connectivity - p.hierarchy));
			} else {
				return 3;
			}
		};
		people[i].influence = compute_inf(people[i]);
	}

	return people;
}
