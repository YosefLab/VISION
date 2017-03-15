root = exports ? this

active_layers = []
layer_mask = {}
curr_mask = 0
layers = []

curr_pulse = ""
p_active = false
p_des = {}
p_seq = {}
active_p = ""
ps = []

b_date = ""
e_date = ""
insight_active = false
insight_nodes = {}

Array.prototype.remove = (args...) ->
	output = []
	for arg in args
		index = @indexOf arg
		output.push @splice(index, 1) if index isnt -1
	output = output[0] if args.length is 1
	output

Array.prototype.contains = (val) ->
	this.forEach (a) ->
		if a == val	
			return true

	return false


RadialPlacement = () ->
	values = d3.map()
	increment = 20
	radius = 200
	center = {"x": 0, "y": 0}
	start = -120
	current = start

	radialLocation = (center, angle, radius) ->
		x = (center.x + radius * Math.cos(angle * Math.PI / 180))
		y = (center.y + radius * Math.sin(angle * Math.PI / 180))
		{"x": x, "y": y}

	placement = (key) ->
		value = values.get(key)
		if !values.has(key)
			value = place(key)
		value

	place = (key) ->
		value = radialLocation(center, current, radius)
		values.set(key, value)
		current += increment
		value

	setKeys = (keys) ->
		values = d3.map()
		increment = 360 / keys.length
		keys.forEach (k) -> place(k)


	placement.keys = (_) ->
		if !arguments.length
			return d3.keys(values)
		setKeys(_)
		placement
	
	placement.center = (_) ->
		if !arguments.length
			return center
		center = _
		placement

	placement.radius = (_) ->
		if !arguments.length
			return radius
		radius = _
		placement

	placement.start = (_) ->
		if !arguments.length
			return start
		current = start
		placement

	placement.increment = (_) ->
		if !arguments.length
			return increment
		increment = _
		placement

	return placement

Network = () ->

	width = 900
	height = $(window).height()

	allData = []
	curLinksData = []
	curNodesData = []
	linkedByIndex = {}

	nodesG = null
	linksG = null

	node = null
	link = null

	
	layout = "force"
	filter = "all"
	sort = "mu"
	groupCenters = null	
	
	force = d3.layout.force()
	nodeColors = d3.scale.category20()
	layerColors = d3.scale.category10()
	simColors = ["green", "red", "purple"]
	tooltip = Tooltip("vis-tooltip", 230)
	summary = SummaryBox()

	charge = (node) -> Math.pow(node.size, 1) / 4


	network = (selection, data) ->
		allData = setupData(data)


		vis = d3.select(selection).append("svg")
			.attr("width", width)
			.attr("height", height)
		linksG = vis.append("g").attr("id", "links")
		nodesG = vis.append("g").attr("id", "nodes")
		

		force.size([width, height])
		setLayout("force")
		setFilter("all")

		update()

	update = () ->
		curNodesData = filterNodes(allData.nodes)
		curLinksData = filterLinks(allData.links, curNodesData)
		updateNodes()
		force.nodes(curNodesData)
	
		if layout == "radial"
			people = sortByMu(curNodesData, curLinksData)
			updateCenters(people)

			
		if layout == "force"
			force.links(curLinksData)
			updateLinks()
		else
			if link
				link.data([]).exit().remove()
				link = null	

		legend = d3.select("svg").selectAll("g.legend")
			.data(nodeColors)
			.enter().append("g")
			.attr("class", "legend")
		
		l_y = $(window).height() - 50
		l_x = 50		
		incr_x = 55
		incr_y = 20

		allMU = []
		i = 0

		allData.nodes.forEach (n) ->
			if n.mu not in allMU
				allMU.push n.mu
				curr_x = l_x + (i * incr_x)
				legend.append("rect")
					.attr("x", () -> curr_x)
					.attr("y", l_y)
					.attr("width", incr_x)
					.attr('height', incr_y)
					.style('fill', () -> nodeColors(n.mu))
					.style("opacity", 0.8)
					.on("mouseover", () ->
						node.each (d) ->
							elem = d3.select(this)
							if d.mu != n.mu
								elem.style("opacity", 0.2)
							if d.mu == n.mu
								elem.style("opacity", 1.0).style("fill", (d) -> nodeColors(d.mu)))
					.on("mouseout", () ->
						#node.each (d) ->
						#	elem = d3.select(this)
						#	elem.style("opacity", 1.0))
						network.layer())
				legend.append("text")
					.attr("x", curr_x + 3)
					.attr("y", l_y + 15)
					.attr("class", 'legend-subs')
					.text(n.mu) 
					.on("mouseover", () ->
						node.each (d) ->
							elem = d3.select(this)
							if d.mu != n.mu
								elem.style("opacity", 0.2)
							if d.mu == n.mu
								elem.style("opacity", 1.0).style("fill", (d) -> nodeColors(d.mu)))
					.on("mouseout", () ->
						network.layer())
				
				i += 1	
	
		# Start force graph
		force.start()
		

	network.toggleLayout = (newLayout) ->
		force.stop()
		setLayout(newLayout)
		network.layer()
		update()

	network.add_layer = (link) ->
		force.stop()
		if link in active_layers
			active_layers.remove(link)
		else
			active_layers.push link

		update()

	network.layer = () ->
		
		node.each (d) -> d.layered = false
		if link
			link.each (l) -> l.layered = false
		if not insight_active	
			if "org" in active_layers and active_layers.length == 1
				node.each (d) ->
					elem = d3.select(this)
					elem.style("fill", (d) -> nodeColors(d.mu))
						.style("stroke-width", 1.0)
						.style("stroke", "#ddd")
						.style("opacity", 1)
						.style("r", (d) -> d.radius)
					d.layered = false
				summary.addSummaryBoxContent(active_layers, layer_mask, allData, insight_nodes, active_p, p_des, p_seq, curr_pulse, d3.event)
				update()
				return 
		
			master_mask = 0
			for l in active_layers
				master_mask += layer_mask[l]	
		
			node.each (d) ->
				for l in active_layers
					if (d.hyperstring & master_mask) == master_mask
						elem = d3.select(this)
						elem.style("fill", (e) -> nodeColors(e.mu))
							.style("stroke", "black")
							.style("stroke-width", 2.0)
							.style("opacity", 1)
							.style("r", (d) -> d.radius + 2)
						d.layered = true
					

			unlayered = node.filter (n) -> !n.layered

			unlayered.each (u) ->
				elem = d3.select(this)
				elem.style("opacity", 0.2)
					.style("stroke-width", 1.0)
					.style("stroke", "#ddd")
					.style('fill', "#FFF")

			summary.addSummaryBoxContent(active_layers, layer_mask, allData, insight_nodes, active_p, p_des, p_seq, curr_pulse, d3.event)
			update()


	network.pulse = () ->
		
		if p_active
			curr_nodes = p_seq[active_p][curr_pulse]
			node.each (n) -> 
				elem = d3.select(this)
				if n.id in curr_nodes
					elem.classed("active-node", true)
						.style("r", 5)
						.style("fill", (d) -> nodeColors(d.mu))
				else
					elem.classed("active-node", false)
						.style("r", (d) -> d.radius)
						.style("opacity", 0.2)
						.style("fill", "#FFF")
			summary.addSummaryBoxContent(active_layers, layer_mask, allData, insight_nodes, active_p, p_des, p_seq, curr_pulse, d3.event)


		else
			node.each (n) ->
				elem = d3.select(this)
				elem.classed("active-node", false)
					.style("r", (d) -> d.radius)
					.style("opacity", 1.0)
					.style("fill", (d) -> nodeColors(d.mu))
	
		

	network.updateData = (newData) ->
		allData = setupData(newData)
		link.remove()
		node.remove()
		update()

	network.nodeOverview = (d, i) ->
		node.each (n) ->
			elem = d3.select(this)
			elem.style('stroke', (d) -> if d.layered then 'black' else "#ddd")
				.style("stroke-width", (d) ->  if d.layered then 2.0 else 1.0)
				.style("opacity", 0.2)
			if n.overview and n.id != d.id
				n.active = false
				summary.hideNodeOverview(allData)
		if d.overview
			node.each (n) ->
				d3.select(this).style("opacity", (d) ->
					if (active_layers.length > 0 and layer_mask[active_layers[0]] != 0 ) or active_layers.length > 1 
						if d.layered then 1.0 else 0.2
					else
						1.0)
					.style("fill", (d) ->
						if (active_layers.length > 0 and layer_mask[active_layers[0]] != 0) or active_layers.length > 1
							if d.layered then nodeColors(d.mu) else "#FFF"
						else
							nodeColors(d.mu))
			summary.hideNodeOverview(allData)
			d.overview= false
		else
			d3.select(this).style('stroke', 'black')
				.style('stroke-width', 2.0)
				.style('opacity', 1.0)
				.style('fill', (d) -> nodeColors(d.mu))
			summary.showNodeOverview(allData, d, insight_nodes, active_p, p_seq, curr_pulse)
			d.overview = true

	network.viewInsights = () ->
		insight_nodes = {}
		if insight_active
			allData.insights.forEach (i) ->
				date = (i.date.split " ")[0]
			
				day = parseInt((date.split '-')[2])
				month = parseInt((date.split '-')[1])
				year = parseInt((date.split '-')[0])

				b_day = parseInt((b_date.split '-')[2])
				b_month = parseInt((b_date.split '-')[1])
				b_year = parseInt((b_date.split '-')[0])

				e_day = parseInt((e_date.split '-')[2])
				e_month = parseInt((e_date.split '-')[1])
				e_year = parseInt((e_date.split '-')[0])
				if year >= b_year and year < e_year
					ins = [(i.date.split " ")[0], i.insight, i.source]
					if insight_nodes[i.id]
						insight_nodes[i.id].push ins
					else
						insight_nodes[i.id] = [ins]
				if year >= b_year and year == e_year
					if month >= b_month and month <= e_month
						if day >= b_day and day <= e_day
							ins = [(i.date.split " ")[0], i.insight, i.source]
							if insight_nodes[i.id]
								insight_nodes[i.id].push ins
							else	
								insight_nodes[i.id] = [ins]
			node.each (n) ->
				elem = d3.select(this)
				insight_ids = Object.keys(insight_nodes)
				if n.id.toString() in insight_ids
					elem.classed('active-node', true)
						.style('r', 10)
						.style('fill', (d) -> nodeColors(d.mu))
				else
					elem.classed('active-node', false)
						.style('r', (d) -> d.radius)
						.style("opacity", 0.2)
						.style("fill", "#FFF")
		else
			node.classed("active-node", false)
				.style("r", (d) -> d.radius)
			network.layer()

		summary.addSummaryBoxContent(active_layers, layer_mask, allData, insight_nodes, active_p, p_des, p_seq, curr_pulse, d3.event)
	
	updateHTML = (name, layer) ->
		if "_" in name
			r_name = ""
			words = name.split "_" 
			i = 0
			while i < words.length
				w = words[i]
				if i == words.lenght - 1
					r_name += w
				else
					r_name += w + " "
				i += 1
			capitalized = r_name.charAt(0).toUpperCase() + r_name.slice(1)
		else
			capitalized = name.charAt(0).toUpperCase() + name.slice(1)
	
		if layer
			words = name.split " "
		else
			words = name.split "_"
		
		i = 0
		_id = ""
		while i < words.length
			w = words[i].toLowerCase()
			if i == words.length - 1
				_id += w
			else
				_id += w + "-"
			i += 1
		
		a = document.createElement('i')
		a.setAttribute("id", _id)
		a.setAttribute('class', 'icon-ok-outline')
		a.innerHTML = capitalized
		a.addEventListener "click", (d) -> 
			newTab = d3.select(this).attr("id")
			if layer
				network.add_layer(newTab)
				activate_layers(newTab)
				network.layer()
			else
				num = ""
				ps.forEach (p) ->
					if p == newTab
						num = p
				if active_p == num
					active_p = ""
					p_active = false
				else
					active_p = num
					p_active = true
				activate_processes(num, newTab)
				curr_pulse = 0
				network.pulse()
		if layer
			if layers.length != 0
				b = document.createElement("br")
				document.getElementById("layer-sub").appendChild(b)	
			document.getElementById("layer-sub").appendChild(a)
		else
			document.getElementById("process-sub").appendChild(a)
				
	setupData = (data) ->
		data.groups.forEach (g) ->
			if g.name.toLowerCase() != "total"
				updateHTML(g.name, true)
			words = g.name.split " "
			i = 0
			_id = ""
			while i < words.length
				w = words[i].toLowerCase()
				if i == words.length - 1
					_id += w
				else 
					_id += w + "-"
				i += 1
			if _id.toLowerCase() != "total"
				layers.push _id
				layer_mask[_id] = g.mask

		data.processes.forEach (p) ->
			updateHTML(p.name, false)
			words = p.name.split "_"
			i = 0
			_id = ""
			while i < words.length
				w = words[i].toLowerCase()
				if i == words.length - 1
					_id += w
				else
					_id += w + "-"
				i += 1
			ps.push _id

			seq_arr = []
			des_arr = []
			for k of p
				parts = k.split "_"
				if parts.length > 1
					if parts[0] == 'seq'
						seq_arr.push p[k]
					if parts[0] == 'des'
						des_arr.push p[k]
			p_seq[_id] = seq_arr
			p_des[_id] = des_arr

		influence_arr = compute_radii(data)

		data.nodes.forEach (n) ->
			# set initial x/y to values within the width/height
			# of the visualization
			n.x = Math.floor(width * n.hierarchy)
			n.y = Math.floor(height * n.hierarchy)
			
			find_radius = (n, influence_arr) ->	
				for p in influence_arr
					if p.name == n.name
						return p.influence
			n.radius = find_radius(n, influence_arr)

		nodesMap = mapNodes(data.nodes)
		data.links.forEach (l) ->
			l.source = nodesMap.get(l.source)
			l.target = nodesMap.get(l.target)
			linkedByIndex["#{l.source.id}, #{l.target.id}"] = 1
		data
	
	# Helper function to map node id's to node objects.
	# Returns d3.map of ids -> nodes
	mapNodes = (nodes) ->
		nodesMap = d3.map()
		nodes.forEach (n) ->
			nodesMap.set(n.id, n)
		nodesMap

	nodeCounts = (nodes, attr) ->
		counts = {}
		nodes.forEach (d) ->
			counts[d[attr]] ?= 0
			counts[d[attr]] += 1
		counts

	neighboring = (a,b) ->
		linkedByIndex[a.id + "," + b.id] or
			linkedByIndex[b.id + "," + a.id]

	filterNodes = (allNodes) ->
		filteredNodes = []
		allNodes.forEach (a) -> a.active = false
		
		if active_layers.length > 0

			filteredNodes = allNodes.filter (n) ->
				if !n.active
					for l in active_layers
						if ((n.hyperstring & layer_mask[l]) == layer_mask[l]) or "org" in active_layers
							n.active == true
							return true
				return false
			
		filteredNodes


	sortByMu = (nodes, links) ->
		people  = []
		final = []
		counts = nodeCounts(nodes, "mu")
		people = d3.entries(counts).sort (a, b) ->
			b.value - a.value
		people.forEach (v) -> 
			node.each (n) ->
				if n.mu == v.key
					final.push n.name
		final

	updateCenters = (people) ->
		if layout == "radial"
			groupCenters = RadialPlacement().center({"x":width/2, "y":height/2 - 50})
				.radius(600).increment(18).keys(people)
	
	filterLinks = (allLinks, curNodes) ->
		curNodes = mapNodes(curNodes)
		filtered_links = []

		filtered_links = allLinks.filter (l) ->
			if l.layer == 0
				return curNodes.get(l.source.id) and curNodes.get(l.target.id)
		return filtered_links
		
			
	updateNodes = () ->
		node = nodesG.selectAll("circle.node")
			.data(curNodesData, (d) -> d.id)

		
		node.enter().append("circle")
			.attr("class", "node")
			.attr("cx", (d) -> d.x)
			.attr("cy", (d) -> d.y)
			.attr("r", (d)-> d.radius)
			.style("fill", (d) -> nodeColors(d.mu))
			.style("stroke-width", 1.0)
			.call(force.drag)

		node.exit().remove()


		node.on("mouseover", showDetails)
			.on("mouseout", hideDetails)
			.on("click", network.nodeOverview)
		

	updateLinks = () ->
		link = linksG.selectAll("line.link")
			.data(curLinksData, (d) -> "#{d.source.id}_#{d.target.id}")
		link.enter().append("line")
			.attr("class", "link")
			.attr("stroke", "#555")
			.attr("stroke-opacity", 0)
			.attr("x1", (d) -> d.source.x)
			.attr("y1", (d) -> d.source.y)
			.attr("x2", (d) -> d.target.x)
			.attr("y2", (d) -> d.target.y)

		link.exit().remove()

	setLayout = (newLayout) ->
		layout = newLayout
		if layout == "force"
			force.on("tick", forceTick)
				.charge(-20)
				.linkDistance(10)
		else if layout == "radial"
			force.on("tick", radialTick)
				.charge(charge)

	setFilter = (newFilter) ->
		filter = newFilter	

	forceTick = (e) ->
		node
			.attr("cx", (d) -> d.x)
			.attr("cy", (d) -> d.y)
				
		link
			.attr("x1", (d) -> d.source.x)
			.attr("y1", (d) -> d.source.y)
			.attr("x2", (d) -> d.target.x)
			.attr("y2", (d) -> d.target.y)

	radialTick = (e) ->
		node.each(moveToRadialLayout(e.alpha))

		node
			.attr("cx", (d) -> d.x)
			.attr("cy", (d) -> d.y)

		if e.alpha < 0.009
			force.stop()
			updateLinks()

	moveToRadialLayout = (alpha) ->
		k = alpha * 0.1
		(d) ->
			centerNode = groupCenters(d.name)
			mask = 0
			for l in active_layers
				mask += layer_mask[l]
			d.x += (centerNode.x - d.x) * k
			d.y += (centerNode.y - d.y) * k
	
	strokeFor = (d) ->
		d3.rgb(nodeColors(d.mu)).darker().toString()


	showDetails = (d,i) ->
		content = '<p class="main title">' + d.name + '</span></p>'
		content += '<p class="main">' + d.position + '</span></p>'
		content += '<p class="main">' + d.mu + '</span></p>'
		tooltip.showTooltip(content, d3.event)

		if link	
			link.style("stroke-opacity", (l) ->
					if layout == "force"
						if (l.source.id == d.id or l.target.id == d.id) then 1.0 else 0.0
					else if layout == "radial"
						if l.layer == 0
							if (l.source.id == d.id or l.target.id == d.id) then 1.0 else 0.0
					else
						0.0
				)
	
		node.style("stoke", (n) -> 
			if (n.layered or n.overview) then "black" else "#ddd")
			.style("stroke-width", (n) ->
				if (n.overview or n.layered) then 2.0 else 1.0)
		

		d3.select(this).style("stroke", "black")
			.style("stroke-width", 2.0)
	
	hideDetails = (d,i) ->
		tooltip.hideTooltip()
		node_summary = false
		node.forEach (n) ->
			if n.overview
				node_summary = true
		node.style("stroke", (n) -> if (n.overview or n.layered) then 'black' else '#ddd')
			.style("stroke-width", (n) -> if (n.overview or n.layered) then 2.0 else 1.0)
		if link
			link.style("stroke-opacity", (l) ->
				if l.layer == "" and l.seq < curr_pulse
					1
				else if (l.layered and layout=="force")
					1 
				else 
					0
			)

	return network


activate_layers = (group) ->
	d3.selectAll("#layers i").classed("active", false).classed('icon-ok', false).classed('icon-ok-outline', true)
	for l in active_layers
		d3.select("#layers ##{l}").classed("active", true).classed('icon-ok', true).classed('icon-ok-outline', false)

activate_processes = (num, p) ->
	d3.selectAll("#process i").classed("active", false).classed("icon-ok", false).classed("icon-ok-outline", true)
	if num == active_p
		d3.select("#process ##{p}").classed('active', true).classed('icon-ok', true).classed('icon-ok-outline', false)

activate_layout = (link) ->
	d3.selectAll("#layouts i").classed("active", false).classed("icon-ok", false).classed("icon-ok-outline", true)
	d3.select("#layouts ##{link}").classed("active", true).classed('icon-ok', true).classed('icon-ok-outline', false)

activate_insights = () ->
	if insight_active
		insight_active = false
		d3.select("#start").classed("active", false).classed("icon-ok", false).classed("icon-ok-outline", true)
	else
		insight_active = true
		d3.select("#start").classed("active", true).classed('icon-ok', true).classed("icon-ok-outline", true)

activate = (group, link) ->
	d3.selectAll("##{group} a").classed("active", false)
	d3.select("##{group} ##{link}").classed("active", true)

$ ->
	myNetwork = Network()
	
	d3.selectAll("#layouts i").on "click", (d) ->
		newLayout = d3.select(this).attr("id")
		activate_layout(newLayout)
		myNetwork.toggleLayout(newLayout)

	$(document).on "keyup", (d) ->
		if active_p == ""
			return
		if d.which == 39
			if curr_pulse < p_seq[active_p].length - 1
				curr_pulse += 1
				console.log(curr_pulse)
				myNetwork.pulse()
		if d.which == 37
			if curr_pulse == 0
				return
			curr_pulse -= 1
			myNetwork.pulse()

	$("#begin-date").on 'change', (d) -> 
		b_date = this.value

	$("#end-date").on 'change', (d) ->
		e_date = this.value

	d3.selectAll("#start").on "click", (d) ->
		activate_insights()
		myNetwork.viewInsights()

	d3.json "data/network.json", (json) ->
		myNetwork("#vis", json)	
