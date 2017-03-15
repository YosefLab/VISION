function person() {
	this.name = "";
	this.num_groups = 0;
	this.hierarchy = 0;
	this.mu = "";
	this.connectivity = 0;
}

function group() {
	this.name = "";
	this.budget = 0;
	this.fte = 0;
	this.mu = "";
	this.mask = 0;
}

function arr_contains(arr, val) {
	for (i = 0; i < arr.length; i++) {
		if (arr[i] == val) {
			return true;
		}
	}
	return false;
}

function object_contains(obj, val) {
	var keys = [];
	for (var k in obj) { keys.push(k); }
	for (var i = 0; i < keys.length; i++) {
		if (keys[i] == val) {
			return true;
		}
	}
	return false;
}

function object_size(obj) {
	keys = [];
	for (k in obj) {
		keys.push(k);
	}
	return keys.length;
}

function SummaryBox() {

	function addSummaryBoxContent(layers, layer_mask,  data, insights, active_p, p_des, p_seq, curr, event) {
		sbid = "#vis-summary";
		content = "<h3>Network Overview</h3>";
		if (object_size(insights) == 0 && !active_p) {
			master_mask = 0;
			if (layers.length > 0) {
				for (i = 0; i < layers.length; i++) {
					master_mask += layer_mask[layers[i]]
				}
			
				var g = [];
		
				data.groups.forEach(function(d, i) {
					var ng = new group();
					ng.budget = d.budget;
					ng.fte = d.fte;
					ng.mask = d.mask;
					var words = d.name.split(" ");
					var _name = "";
					var i = 0;
					while (i < words.length) {
						var w = words[i].toLowerCase();
						if (i == words.length - 1) {
							_name += w;
						} else{
							_name += w + "-";
						}
						i += 1;
					}
					ng.name = _name;
					g.push(ng);
				});



				var people = [];
				var mus = {};	
				data.nodes.forEach(function(d, i) {
					if ((master_mask & d.hyperstring) == master_mask) {
						var np = new person();
						np.name = d.name;
						np.hierarchy = d.hierarchy;
						np.mu = d.mu

						groups = d.hyperstring;
						while (groups > 0) {
							if (groups & 1 == 1) {
								np.num_groups += 1;
							}
							groups = groups >> 1;
						}
						people.push(np);
					}
				});
				
			
				var avg_num_groups = 0;
				var avg_connect = 0;
				if (people.length == 0) {
					content += '<p>' + "None" + '</span></p>';
					content += '<p>' + "Average number of groups:\tN/A" + '</span></p>';
				} 
				else {
					content += "<table>";
					content += "<tr><th>Name</th><th>MU</th></tr>";
					for (j = 0; j < people.length; j++) {
						if ((layers.length == 1 && layer_mask[layers[0]] != 0) || layers.length > 1) {
							if (mus[people[j].mu]) {
								mus[people[j].mu] += 1;
							} else {
								mus[people[j].mu] = 1;
							}
							content += "<tr>";
							content += '<td>' + people[j].name + '</td>';
							content += '<td>' + people[j].mu + '</td>';
							content += '</tr>';
						}
						avg_num_groups += people[j].num_groups;
						data.links.forEach(function(l, i) {
							if (l.source.name == people[j].name || l.target.name == people[j].name) { avg_connect += 1; }
						});
						avg_connect += people[j].num_groups;
					}
					content += "</table>";

					avg_connect /= people.length;
					avg_num_groups /= people.length;
					content += '<p class="sub">' + "Average group membership: \t" + Number(avg_num_groups.toFixed(1)) + '</span></p>';
					content += '<p class="sub">' + "Average connectivity: \t" + Number(avg_connect.toFixed(1)) + '</span></p>';
					content += "<table>";
					content += "<tr><th>MU</th><th>Percent</th></tr>";
					for (var k in mus) {
						content += '<td>'+ k + '</td>';
						content += '<td>' + Math.round(mus[k] / people.length * 100, 3) + "%</td>";
						content += "</tr>"
					}
					content += "</table>";
				}
				
				
				var total_budget = 0;
				var agg_budget = 0;
				var total_fte = 0;
				for (j = 0; j < g.length; j++) {
					if (g[j].name.toLowerCase() == "total") {
						total_budget = g[j].budget;
					}
					if (layers.length == 1 && layer_mask[layers[0]] == 0) {
						agg_budget = total_budget;
					}
					if (arr_contains(layers, g[j].name.toLowerCase()) && !(layers.length == 1 && layer_mask[layers[0]] == 0)) {
						if (g[j].name.toLowerCase() != "total") {
							agg_budget += g[j].budget;
							total_fte += g[j].fte;
						}
					}
				}
				total_fte = total_fte.toFixed(2);
				perc_budget = (agg_budget / total_budget * 100).toFixed(2);
				if (agg_budget == 0) {
					perc_budget = 0.00;
				}
				content += '<p class="sub">Total FTE:\t' + total_fte + '</span></p>';
				content += '<p class="sub">Cumulative Budget:\t' + "$" + agg_budget + '</span></p>';
				content += '<p class="sub">Percent of Total Budget:\t' + perc_budget + '%</span></p>';
			}		
		} else if (object_size(insights) > 0 && !active_p) {
			content += "<table>";
			content += "<tr><th>Date</th><th>Insight</th><th>Source</th></tr>";
			for (var k in insights) {
				for (var j = 0; j < insights[k].length; j++) {
					content += "<tr>";
					content += "<td>" + insights[k][j][0] + "</td>";
					content += "<td>" + insights[k][j][1] + "</td>";
					content += "<td>" + insights[k][j][2] + "</td>";
					content += "</tr>"
				}
			}
			content += "</table>"; 
			
			content += '<p class="sub">Total Number of Insights:\t' + object_size(insights) + "</p>";
		} else if (active_p != "") {
			content += '<table>';
			content += "<tr><th>Process Step Description - " + (curr+1) + "</th></tr>";
			content += "<tr><td>";
			content += p_des[active_p][curr];
			content += "</td></tr>";
			content += "</table>";

			content += "<table>";
			content += "<tr><th>Name</th><th>MU</th></tr>";
							
			data.nodes.forEach(function(n, i) {
				if (arr_contains(p_seq[active_p][curr], n.id)) {
					content += "<tr>";
					content += "<td>" + n.name + "</td>";
					content += "<td>" + n.mu + "</td>";
					content += "</tr>";
				}
			});
			content += "</table>";

		}
		$(sbid).html(content);
		
		var xOffset = 20;
		var yOffset = 10;

	}

	function showNodeOverview(data, n, insights, active_p, p_seq, curr, event) {
		_id = "#text-container-vis-node";
		content = "";
		var hyperstring = n.hyperstring;
		var p = new person();
		p.name = n.name;	
		while (hyperstring > 0) {
			if ((hyperstring & 1) == 1) {
				p.num_groups += 1;
			}
			hyperstring = hyperstring >> 1;
		}
	
		data.links.forEach(function(l, i) {
			if (l.source.id == n.id || l.target.id == n.id) {
				p.connectivity += 1 
			} 
		});
		p.connectivity += p.num_groups;
		content += '<p>' + p.name + "</p>";	
		content += '<p>Group Membership:\t' + p.num_groups + "</p>";
		content += '<p>Connectivity:\t' + p.connectivity + "</p>";
		if (object_size(insights) > 0) {
			content += '<table class="insight-table">';
			content += '<tr><th>Date</th><th>Insight</th><th>Source</th></tr>'
			if (object_contains(insights, n.id)) {
				for (var i = 0; i < insights[n.id].length; i++) {
					content += '<tr>';
					content += '<td id="insight">' + insights[n.id][i][0] + "</td>";
					content += '<td id="insight">' + insights[n.id][i][1] + "</td>";
					content += '<td id="insight">' + insights[n.id][i][2] + "</td>";
					content += '</tr>';
				}
			}
			content += "</table>";
		} else if (active_p && arr_contains(p_seq[active_p][curr], n.id)) {
			content += "<br><br";	
			content += "<p class='sub'>INSERT PERSON'S RESPONSIBILITIES HERE</p>";
		}

		$(_id).html(content);	
	}
	
	function hideNodeOverview(data, event) {
		_id = "#text-container-vis-node";
		content = "";
		$(_id).html(content);
	}

	return {
		addSummaryBoxContent: addSummaryBoxContent,
		showNodeOverview: showNodeOverview,
		hideNodeOverview: hideNodeOverview,
	}
}
