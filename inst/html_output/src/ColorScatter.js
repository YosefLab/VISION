/*
   Initializes a zoomable scatter plot in the element "parent"
   parent = , for example, "#chart_div"
*/
function ColorScatter(parent, colorbar, legend)
{
    if(colorbar === undefined) { this.colorbar_enabled = false;}
    else {this.colorbar_enabled = colorbar;}

    if (legend == undefined) { this.legend_enabled = false; }
    else {this.legend_enabled = legend;}

    var colorbar_height = 20;
    var colorbar_width = 200;

    var self = this;
    var xdomain = [-2, 2];
    var ydomain = [-2, 2];

    var margin = {right: 20, left: 40, top: 20, bottom: 20};
    this.width = $(parent).width() - margin.left - margin.right;
    this.height = $(parent).height() - margin.top - margin.bottom;

    this.circle_radius = 4.0

    this.x = d3.scale.linear()
        .domain(xdomain)
        .range([0, self.width]);

    this.y = d3.scale.linear()
        .domain(ydomain)
        .range([self.height, 0]);

    this.xAxis = d3.svg.axis()
        .scale(self.x)
        .orient("bottom")
        .tickSize(-self.height);

    this.yAxis = d3.svg.axis()
        .scale(self.y)
        .orient("left")
        .ticks(5)
        .tickSize(-self.width);

    this.zoom = d3.behavior.zoom()
        .x(self.x)
        .y(self.y)
        .scaleExtent([0.2, 32])
        .on("zoom", self.redraw());

    //This gets overwritten when you set data
    this.colorScale = null

    this.tip = d3.tip()
        .attr('class', 'd3-tip')
        .offset([-10, 0])
        .html(function(d){ return d[3]+": "+d[2]; });


    this.svg = d3.select(parent).append("svg")
        .attr("height", self.height + margin.top + margin.bottom)
        .attr("width", self.width + margin.right + margin.left)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
        .call(self.zoom)
        .call(self.tip);

    this.svg.append("rect")
        .attr("width", self.width)
        .attr("height", self.height);


    this.svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + self.height + ")")
        .call(self.xAxis);

    this.svg.append("g")
        .attr("class", "y axis")
        .call(self.yAxis);
    

    if(this.colorbar_enabled){
        this.colorbar = this.svg.append("g")

        this.colorbar.append("rect")
            .attr("x", this.width - colorbar_width - 10)
            .attr("y", 10)
            .attr("width", colorbar_width)
            .attr("height", colorbar_height);
    }

    if (this.legend_enabled) {
        this.legend = this.svg.append("g")
            .style("stroke", "black")
            .style("stroke-width", "1px")
    }


    this.points = [];

    this.selected = -1;
    this.selectedPoints = {};
    this.selected_links = [];

    this.hovered = -1;

    this.last_event = -1;

    this.tree_points = [];
    this.tree_adj = [];
    this.curr_tree_node = {};
}

ColorScatter.prototype.clearData = function() {
    this.points = []
    this.tree_points = []
    this.tree_adj = []
}

ColorScatter.prototype.setData = function(points, isFactor,
    full_color_range, selected_cells, diverging_colormap)
{

    if(full_color_range === undefined){
        full_color_range = false
    }

    if(selected_cells === undefined){
        selected_cells = points.map(x => x[3]) //select all cells
    }

    if(diverging_colormap === undefined){
        diverging_colormap = true
    }

    selected_cells = _.keyBy(selected_cells, x => x)

    this.points = points;

    var cvals = points
        .filter(e => e[3] in selected_cells) // filter for selected cells
        .map(e => e[2])                      // extract 3rd column

    //Adjust circle size based on number of points
    this.circle_radius = this.pointsToRadius(points.length)

    if(isFactor)
    {
        //Find unique values
        var unique = d3.set(cvals).values();
        if(unique.length <= 10) { this.colorScale = d3.scale.category10();}
        else { this.colorScale = d3.scale.category20();}

        this.colorScale.domain(unique);
        this.setColorBar();


        this.setLegend(this.colorScale, unique);

    } else if(cvals[0] !== null) {

        cvals.sort(d3.ascending); // Needed for quantile
        var low, high, mid;
        if(full_color_range){
            low = d3.min(cvals)
            high = d3.max(cvals)
            mid = (low + high)/2
        } else {
            low = d3.quantile(cvals, 0.02)
            high = d3.quantile(cvals, 0.98)
            mid = (low + high)/2
        }

        if(diverging_colormap){
            this.colorScale = d3.scale.linear()
                .domain([low, mid, high])
                .range(["blue", "green", "red"]);
        } else {
            this.colorScale = d3.scale.linear()
                .domain([low, mid, high])
                .range(["#d8d8d8", "#395252", "#000000"]);
        }

        // Format the bound labels
        var label_low = parseFloat(low.toFixed(2)).toString();
        var label_high = parseFloat(high.toFixed(2)).toString();

        this.setLegend();
        this.setColorBar(this.colorScale.range(),
            label_low,
            label_high);

    } else {

        this.colorScale = function(){return "steelblue"}

        this.setLegend();
        this.setColorBar();
    }

    // Compute colors for all points, add to points[n][4]
    var self = this
    this.points.forEach(function(x) {
        x[4] = self.colorScale(x[2])
    })
}

ColorScatter.prototype.pointsToRadius = function(n_points)
{
    // Pick a point size based on the number of scatter
    // plot points

    // Using a piecewise function for this
    var a, b


    var scale_factor = Math.min(this.width, this.height)/600

    if(n_points <= 10000){
        a = 116500
        b = .284
    } else {
        a = 19800
        b = .748
    }

    return Math.pow(a/n_points, b) * scale_factor

}

ColorScatter.prototype.setTreeData = function(tree_points, tree_adj)
{

    this.tree_points = tree_points
    this.tree_adj = tree_adj
}

/*
 * Automatically adjusts the scale if the points are too wide
 * or too narrow
 */
ColorScatter.prototype.autoZoom = function() {
    var self = this;
    var xvals = self.points.map(function(e){return e[0];});
    var yvals = self.points.map(function(e){return e[1];});

    xvals.sort(d3.ascending); // Needed for quantile
    yvals.sort(d3.ascending); // Needed for quantile

    var xmax = d3.quantile(xvals, 0.99);
    var ymax = d3.quantile(yvals, 0.99);

    var xmin = d3.quantile(xvals, 0.01);
    var ymin = d3.quantile(yvals, 0.01);

    var current_x_range = self.x.domain()
    var current_y_range = self.y.domain()

    var do_x_rescale = false;
    var do_y_rescale = false;

    // These cover if more than X % of the data is outside the viewport on
    // either side
    if(xmin < current_x_range[0]) {do_x_rescale = true;}
    if(xmax > current_x_range[1]) {do_x_rescale = true;}

    if(ymin < current_y_range[0]) {do_y_rescale = true;}
    if(ymax > current_y_range[1]) {do_y_rescale = true;}

    // These cover if the data is too scrunched up
    var x_coverage = (xmax-xmin) / (current_x_range[1] - current_x_range[0])
    if(x_coverage < .20) {do_x_rescale = true;}

    var y_coverage = (ymax-ymin) / (current_y_range[1] - current_y_range[0])
    if(y_coverage < .20) {do_y_rescale = true;}

    if(do_x_rescale || do_y_rescale){
        var x_pad = (xmax-xmin)*0.30;
        self.x.domain([xmin-x_pad, xmax+x_pad]);

        self.zoom
            .x(self.x)
    }

    if(do_y_rescale || do_x_rescale){
        var y_pad = (ymax-ymin)*0.30;
        self.y.domain([ymin-y_pad, ymax+y_pad]);

        self.zoom
            .y(self.y)
    }
    this.redraw(true)();
}

ColorScatter.prototype.setLegend = function(colors, values) {
	
    if (!this.legend_enabled) { return; }

    var self = this;

    this.legend.selectAll("text").remove();
    this.legend.selectAll("rect").remove();

    if (colors !== undefined) {

        var square_size = 20
        var l_x = self.width - square_size -10;
        var l_y = 10;
        var i = 0;
        values.forEach(function(n) {
            var curr_y = l_y + (i * (square_size+5));
            self.legend.append("rect")
                .attr("x", l_x)
                .attr("y", curr_y)
                .attr("width", square_size)
                .attr("height", square_size)
                .style("fill", colors(n));
            self.legend.append("text")
                .attr("x", l_x - 5)
                .attr("y", curr_y + square_size/2)
                .attr("class", "legend-subs")
                .attr("text-anchor", "end")
                .attr("alignment-baseline", "middle")
                .text(n);
            i += 1;
        });
    }
	
}	

ColorScatter.prototype.setColorBar = function(colors, label_low, label_high)
{
    if(!this.colorbar_enabled){ return; }

    this.colorbar.selectAll("text").remove();
    this.colorbar.select("defs").remove();

    if(colors === undefined){
        // Clear color bar - useful for factors
        var colorbar_rect = this.colorbar.select('rect')

        colorbar_rect
            .style("fill", "rgba(0,0,0,0)")
            .style("stroke", "none")
    }
    else
    {
        var gradient = this.colorbar.append("svg:defs")
            .append("svg:linearGradient")
            .attr("id", "gradient")
            .attr("x1", "0%")
            .attr("y1", "0%")
            .attr("x2", "100%")
            .attr("y2", "0%")
            .attr("spreadMethod", "pad");

        for(var i = 0; i < colors.length; i++){
            var stop_percentage = Math.floor((100 / (colors.length-1))*i);

            gradient.append("svg:stop")
                .attr("offset", stop_percentage.toString() + "%")
                .attr("stop-color", colors[i])
                .attr("stop-opacity", 0.8);
        }
        var colorbar_rect = this.colorbar.select('rect')

        colorbar_rect
            .style("fill", "url(#gradient)")
            .style("stroke", "black")
            .style("stroke-width", "1px");

        var label_low_x = parseInt(colorbar_rect.attr("x"))
        var label_high_x = label_low_x +
            parseInt(colorbar_rect.attr("width")) - 10

        this.colorbar.append("text")
            .attr("text-anchor", "middle")
            .attr("x", label_low_x)
            .attr("y", 40)
            .attr("font-size", "10px")
            .text(label_low)
            .attr("fill", "black");

        this.colorbar.append("text")
            .attr("text-anchor", "middle")
            .attr("x", label_high_x)
            .attr("y", 40)
            .attr("font-size", "10px")
            .text(label_high)
            .attr("fill", "black");
    }
};

ColorScatter.prototype.setSelected = function(selected_key, event_id)
{

    var self = this;

    if (typeof(selected_key) === "string"){
        selected_key = [selected_key]
    }

    self.selectedPoints = _.keyBy(selected_key, x => x)

    var circles = self.svg.selectAll("circle")
        .data(self.points);

    if(_.isEmpty(self.selectedPoints)){
        circles
            .classed("selected", false)
            .classed("not-selected", false)

    } else {

        circles
            .classed("selected", function(d){return self.selectedPoints[d[3]] !== undefined})
            .classed("not-selected", function(d){return self.selectedPoints[d[3]] === undefined})
    }

    var event = new CustomEvent('select-cells', {detail: self.selectedPoints })

    window.dispatchEvent(event);

};

ColorScatter.prototype.setHovered_TreeNode = function(node, clusters, event_id) {
	var self = this;
	self.curr_tree_node = node;

	if (event_id === undefined) {
		event_id = Math.random();
	}

	if (this.last_event !== event_id) {
		this.last_event = event_id;
		if (self.curr_tree_node == null) {
			this.svg.selectAll("circle")
				.classed("point-faded", false)
				.classed("point-hover", false);
		}
		else {
			this.svg.selectAll("circle")
				.classed("point-faded", true)
				.classed("point-hover", function(d, i) {
					var sample_name = d[3];
					return (clusters[sample_name] == self.curr_tree_node.name); 
				});
		}
	}
	return;

}

ColorScatter.prototype.hover_cells = function(cell_ids)
{

    //Needed to prevent infinite loops with linked hover and select events
    if(_.isEmpty(cell_ids)){
        //Clear the hover
        this.svg.selectAll("circle")
            .classed("point-faded", false)
            .classed("point-hover", false);
    }
    else{
        this.svg.selectAll("circle")
            .classed("point-faded", true)
            .classed("point-hover", function (d) {
                return cell_ids.indexOf(d[3]) > -1;
            });
    }

};

ColorScatter.prototype.toggleLasso = function(enable) {

    var self = this;
    if (enable) { 

        var lasso_start = function() {
            self.lasso.items()
                .classed({"not-possible": true})
                .style("fill", null)

            self.setSelected([])
        };

        var lasso_draw = function() {
            // Style possible dots

            self.lasso.items().filter(function(d) { return d.possible===true})
                .classed({"not-possible":false, "possible":true});

            // Style not possible dots
            self.lasso.items().filter(function(d) { return d.possible===false})
                .classed({"not-possible": true, "possible": false});
        };

        var lasso_end = function() {

            self.lasso.items()
                .classed({"not-possible":false, "possible":false})

            var cellNames = self.lasso.items()
                .filter(function(d) { return d.selected === true })
                .data()
                .map(d => d[3])

            self.setSelected(cellNames)
            self.redraw(false)();

        };

        self.lasso_area = self.svg.append("rect")
            .attr("width", self.width)
            .attr("height", self.height)
            .style("opacity", 0);

        self.lasso = d3.lasso()
            .closePathDistance(100000)
            .closePathSelect(true)
            .hoverSelect(true)
            .area(self.lasso_area)
            .on("start", lasso_start)
            .on("draw", lasso_draw)
            .on("end", lasso_end);

        self.svg.call(self.lasso);
        self.svg
            .on("mousedown.zoom", null)
            .on("touchstart.zoom", null)
            .on("touchmove.zoom", null)
            .on("touchend.zoom", null)

        self.lasso.items(self.svg.selectAll("circle.scatter"));
    } else {
        self.lasso_area.remove();
        self.svg.call(self.zoom);
    }

}

ColorScatter.prototype.redraw = function(performTransition) {
    var self = this;
    return function(){
        self.svg.select(".x.axis").call(self.xAxis);
        self.svg.select(".y.axis").call(self.yAxis);

        self.svg.call(self.tip);

        var circles = self.svg.selectAll("circle.scatter")
            .data(self.points);


        circles.enter().append("circle")
            .classed("scatter", true);

        circles.style("fill", function(d){return d[4];})
            .attr("r", self.circle_radius * Math.pow(self.zoom.scale(), .5))
            .on("click", function(d){self.setSelected(d[3]);})
            .on("mouseover", function(d,i){self.tip.show(d,i);})
            .on("mouseout", function(d,i){self.tip.hide(d,i);})

        var tree_circles = self.svg.selectAll("circle.tree")
            .data(self.tree_points);

        tree_circles.enter().append("circle")
            .style("fill", '#05ff65')
            .classed("tree", true);

        tree_circles
            .attr("r", 5 * Math.pow(self.zoom.scale(), .5))

        var tree_edges = self.svg.selectAll("line.tree")
            .data(self.tree_adj);

        tree_edges.enter().append("line")
            .attr("stroke", "#05ff65")
            .attr("stroke-width", 2)
            .attr("fill", 'none')
            .classed("tree", true);

        tree_edges
            .attr("stroke-width", 2 * Math.pow(self.zoom.scale(), .5))


        if(performTransition !== undefined && performTransition === true)
        {
            circles
                .transition()
                .duration(1000)
                .attr("cx", function(d){return self.x(d[0]);})
                .attr("cy", function(d){return self.y(d[1]);});

            tree_circles
                .transition()
                .duration(1000)
                .attr("cx", function(d){return self.x(d[0]);})
                .attr("cy", function(d){return self.y(d[1]);});

            tree_edges
                .transition()
                .duration(1000)
                .attr("x1", function(d){return self.x(self.tree_points[d[0]][0])})
                .attr("y1", function(d){return self.y(self.tree_points[d[0]][1])})
                .attr("x2", function(d){return self.x(self.tree_points[d[1]][0])})
                .attr("y2", function(d){return self.y(self.tree_points[d[1]][1])})
        }
        else
        {
            circles
                .attr("cx", function(d){return self.x(d[0]);})
                .attr("cy", function(d){return self.y(d[1]);});

            tree_circles
                .attr("cx", function(d){return self.x(d[0]);})
                .attr("cy", function(d){return self.y(d[1]);});

            tree_edges
                .attr("x1", function(d){return self.x(self.tree_points[d[0]][0])})
                .attr("y1", function(d){return self.y(self.tree_points[d[0]][1])})
                .attr("x2", function(d){return self.x(self.tree_points[d[1]][0])})
                .attr("y2", function(d){return self.y(self.tree_points[d[1]][1])})

        }

        circles.exit().remove();
        tree_circles.exit().remove();
        tree_edges.exit().remove();
    };
};


ColorScatter.prototype.selectCellRange = function(points, lower, upper) {
	
	var self = this;
	var allCellNames = points.map(function(e) { return e[3]; });
    var mappedData = points.map(function(e){return [e[2], e[3]];}); //extract 2nd and 3rd columns

	var filteredData = mappedData.filter(function(e) {
		return e[0] >= lower && e[0] <= upper
	});

	var cellNames = filteredData.map(function(e){return e[1];});
	
	cellNames.forEach(function(e) {
		self.selectedPoints.push(allCellNames.indexOf(e));
	});

	
    var circles = this.svg.selectAll("circle")
        .data(points);
    circles
		.classed("point-selected", function(d, i){return self.selectedPoints.indexOf(i) != -1;})
		.attr("opacity", function(d,i) {
			if (self.selectedPoints.indexOf(i) == -1) {
				return .4
			}
		});
	return cellNames
}

ColorScatter.prototype.getSelected = function() {
	var self = this;

	var allCellNames = this.points.map(function(e) { return e[3]; });
	
	var selectedCellNames = [];
	this.selectedPoints.forEach(function(d) {
		selectedCellNames.push(allCellNames[d]);
	});

	return selectedCellNames;

}

ColorScatter.prototype.unselect = function() {
	var self = this;
	this.selectedPoints.length = 0;	
    var circles = this.svg.selectAll("circle")
        .data(this.points);
    circles
		.classed("point-selected", function(d, i){return self.selectedPoints.indexOf(i) != -1;})
		.attr("opacity", 1.0);
}

