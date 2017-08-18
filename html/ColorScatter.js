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

	this.legend_height = 20;
	this.legend_width = 200;

    var self = this;
    var xdomain = [-2, 2];
    var ydomain = [-2, 2];

    this.margin = {top: 20, right: 20, bottom: 15, left: 40};
    this.width = $(parent).width();
    this.height = $(parent).height();

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
    this.colorScale = d3.scale.linear()
        .domain([0,0.5,1])
        .range(["blue", "green", "red"]);

    this.tip = d3.tip()
        .attr('class', 'd3-tip')
        .offset([-10, 0])
        .html(function(d){ return d[3]+": "+d[2]; });


    this.svg = d3.select(parent).append("svg")
        .attr("height", self.height + self.margin.top + self.margin.bottom)
        .attr("width", self.width)
        .append("g")
        //.call(self.zoom)
        .call(self.tip);

    this.svg.append("rect")
        .attr("width", self.width)
        .attr("height", self.height);


    this.svg.append("g")
        .attr("class", "x axis")
        .call(self.xAxis);

    this.svg.append("g")
        .attr("class", "y axis")
        .call(self.yAxis);
    

    if(this.colorbar_enabled){
        this.colorbar_svg = this.svg
            .attr("width", this.width + this.margin.left + this.margin.right)
            .attr("height", colorbar_height+20);

        this.colorbar = this.colorbar_svg.append("rect")
            //.attr("x", this.width + this.margin.left + this.margin.right - colorbar_width - 120)
            .attr("x", this.width - colorbar_width - 10)
            .attr("y", 10)
            .attr("width", colorbar_width)
            .attr("height", colorbar_height);
    }

    if (this.legend_enabled) {
		this.legend_svg = this.svg
			.attr("width", this.width + this.margin.left + this.margin.right)
			.attr("height", this.legend_height + 20);

		this.legend = this.legend_svg.append("g")
			.attr("x", this.width - this.legend_width - 10)
			.attr("y", 10)
			.attr("width", this.legend_width)
			.attr("height", this.legend_height);
	}


    this.points = [];

    this.selected = -1;
    this.selectedPoints = [];
    this.selected_links = [];

    this.hovered = -1;
    this.hovered_links = [];

    this.last_event = -1;

    this.curr_tree_node = {};
}

ColorScatter.prototype.setData = function(points, isFactor)
{
    this.points = points;
    var cvals = points.map(function(e){return e[2];}); //extract 3rd column

    if(isFactor)
    {
        //Find unique values
        var unique = d3.set(cvals).values();
        if(unique.length <= 10) { this.colorScale = d3.scale.category10();}
        else { this.colorScale = d3.scale.category20();}

        this.colorScale.domain(unique);
        this.setColorBar();


		this.setLegend(this.colorScale, unique);
    }
    else
    {
        cvals.sort(d3.ascending); // Needed for quantile
        var low = d3.quantile(cvals, 0.1);
        var high = d3.quantile(cvals, 0.9);
        var mid = d3.mean(cvals);

        this.colorScale = d3.scale.linear()
            .domain([low, mid, high])
            .range(["#48D1CC", "#7CFC00", "red"]);

        // Format the bound labels
        var label_low = parseFloat(low.toFixed(2)).toString();
        var label_high = parseFloat(high.toFixed(2)).toString();

        this.setLegend();
        this.setColorBar(this.colorScale.range(),
               label_low,
               label_high);

    }

    this.redraw(true)();
};

ColorScatter.prototype.setLegend = function(colors, values) {
	
	if (!this.legend_enabled) { return; } 

	var self = this;

	this.legend.style("display", "block");
	this.legend_svg.selectAll("text").remove();
	//this.legend_svg.selectAll("defs").remove();


	if (colors == undefined) {
		this.legend.style("display", "none");
		this.legend
			.style("fill", "rgba(0,0,0,0)")
			.style("stroke", "none");
	} else {
        this.legend
            .style("stroke", "black")
            .style("stroke-width", "1px")

		incr_x = self.legend_width / values.length;
		l_x = self.width - self.legend_width - 10;
		l_y = 10; 
		i = 0;
		values.forEach(function(n) {
			curr_x = l_x + (i * incr_x);
			self.legend.append("rect")
				.attr("x", curr_x)
				.attr("y", l_y)
				.attr("width", incr_x)
				.attr("height", self.legend_height)
				.style("fill", colors(n));
			self.legend.append("text")
				.attr("x", curr_x + 3)
				.attr("y", l_y + 15)
				.attr("class", "legend-subs")
				.text(n);
			i += 1;
		});
	}


	
}	

ColorScatter.prototype.setColorBar = function(colors, label_low, label_high)
{
    if(!this.colorbar_enabled){ return; }

    this.colorbar_svg.selectAll("text").remove();
    this.colorbar_svg.select("defs").remove();

    if(colors === undefined){
        // Clear color bar - useful for factors
        this.colorbar
            .style("fill","rgba(0,0,0,0)")
            .style("stroke","none");
    }
    else
    {
        var gradient = this.colorbar_svg.append("svg:defs")
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

        this.colorbar
            .style("fill", "url(#gradient)")
            .style("stroke", "black")
            .style("stroke-width", "1px");

        var label_low_x = this.colorbar.attr("x");
        var label_high_x = (parseInt(this.colorbar.attr("x")) +
            parseInt(this.colorbar.attr("width")-10)).toString();

        this.colorbar_svg.append("text")
            .attr("text-anchor", "middle")
            .attr("x", label_low_x)
            .attr("y", 40)
            .attr("font-size", "10px")
            .text(label_low)
            .attr("fill", "white");

        this.colorbar_svg.append("text")
            .attr("text-anchor", "middle")
            .attr("x", label_high_x)
            .attr("y", 40)
            .attr("font-size", "10px")
            .text(label_high)
            .attr("fill", "white");
    }
};

ColorScatter.prototype.setSelected = function(selected_index, event_id)
{

	var self = this;
    if(event_id === undefined){
        event_id = Math.random();
    }

    //Needed to prevent infinite loops with linked hover and select events
    if(this.last_event !== event_id) {
        this.last_event = event_id;
        this.selected = selected_index;
        this.redraw()();
        this.selected_links.forEach(function (e) {
            e.setSelected(selected_index, event_id);
        });
	}

	if (this.selectedPoints.indexOf(selected_index) == -1) {
		this.selectedPoints.push(selected_index);
	} else {
		this.selectedPoints.splice(this.selectedPoints.indexOf(selected_index), 1);
	}

    var circles = this.svg.selectAll("circle")
        .data(this.points);
    circles
		.classed("point-selected", function(d, i){return self.selectedPoints.indexOf(i) != -1;})
		.attr("opacity", function(d,i) {
			if (self.selectedPoints.indexOf(i) == -1) {
				return .4
			}
		});
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

ColorScatter.prototype.setHovered = function(hovered_indices, event_id)
{
    if(event_id === undefined){
        event_id = Math.random();
    }

    //test for single index, and wrap in list
    if(typeof(hovered_indices) === "number"){hovered_indices = [hovered_indices];}
    if(hovered_indices.length === 0){hovered_indices = [-1];}

    //Needed to prevent infinite loops with linked hover and select events
    if(this.last_event !== event_id) {
        this.last_event = event_id;
        this.hover_col = hovered_indices;
        if(hovered_indices.length === 1 && hovered_indices[0] === -1){
             //Clear the hover
            this.svg.selectAll("circle")
                .classed("point-faded", false)
                .classed("point-hover", false);
        }
        else{
            this.svg.selectAll("circle")
                .classed("point-faded", true)
                .classed("point-hover", function (d, i) {
                    return hovered_indices.indexOf(i) > -1;
                });
        }

        this.hovered_links.forEach(function (e) {
            e.setHovered(hovered_indices, event_id);
        });
    }
};

ColorScatter.prototype.toggleLasso = function(enable) {

	var self = this;
	if (enable) { 
		self.lasso_start = function() {
			self.lasso.items()
				.classed({"not_possible": true, "selected": false});
		};

		self.lasso_draw = function() {
			// Style possible dots
		
	
			self.lasso.items().filter(function(d) { return d.possible===true})
				.classed({"not_possible":false, "possible":true});

			// Style not possible dots
			self.lasso.items().filter(function(d) { return d.possible===false})
				.classed({"not_possible": true, "possible": false});
		};
	
		self.lasso_end = function() {

			self.lasso.items().filter(function(d) { return d.selected===true})
				.classed({"not_possible":false, "possible":false})
				.classed("point-selected", true);

			self.lasso.items().filter(function(d) { return d.selected===false })
				.classed({"not_possible": false, "possible": false});

			var cellNames = self.points.map(function(e){return e[3];}); //extract 2nd and 3rd columns

			var selected = self.lasso.items().filter(function(d) { return d.selected === true });
			selected[0].forEach(function(d) {
				self.selectedPoints.push(cellNames.indexOf(d.__data__[3]));
			});	

			self.lasso.items()
				.attr("opacity", function(d,i) {
					if (self.selectedPoints.indexOf(i) == -1) {
						return .4
					}
				});

		};

		self.lasso_area = self.svg.append("rect")
				.attr("width", self.width)
				.attr("height", self.height)
				.style("opacity", 0);

		self.lasso = d3.lasso()
			.closePathDistance(75)
			.closePathSelect(true)
			.hoverSelect(true)
			.area(self.lasso_area)
			.on("start", self.lasso_start)
			.on("draw", self.lasso_draw)
			.on("end", self.lasso_end);

		self.svg.call(self.lasso);

		self.lasso.items(self.svg.selectAll("circle"));
	} else {
		self.lasso_area.remove();
	}

}

ColorScatter.prototype.redraw = function(performTransition) {
    var self = this;
    return function(){
    self.svg.select(".x.axis").call(self.xAxis);
    self.svg.select(".y.axis").call(self.yAxis);

	self.svg.call(self.tip);

    var circles = self.svg.selectAll("circle")
        .data(self.points);

    circles.enter().append("circle").attr("r",4.0);
    circles.style("fill", function(d){return self.colorScale(d[2]);})
        .on("click", function(d,i){self.setSelected(i);})
        .on("mouseover", function(d,i){self.tip.show(d,i); self.setHovered(i);})
        .on("mouseout", function(d,i){self.tip.hide(d,i); self.setHovered(-1);})



    if(performTransition !== undefined && performTransition === true)
    {
        circles
            .transition()
            .duration(1000)
            .attr("cx", function(d){return self.x(d[0]);})
            .attr("cy", function(d){return self.y(d[1]);});
    }
    else
    {
        circles
            .attr("cx", function(d){return self.x(d[0]);})
            .attr("cy", function(d){return self.y(d[1]);});
    }

    circles.exit().remove();
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
