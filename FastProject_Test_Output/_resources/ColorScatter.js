/*
   Initializes a zoomable scatter plot in the element "parent"
   parent = , for example, "#chart_div"
*/
function ColorScatter(parent, colorbar)
{
    if(colorbar === undefined) { this.colorbar_enabled = false;}
    else {this.colorbar_enabled = colorbar;}

    var colorbar_height = 20;
    var colorbar_width = 200;

    var self = this;
    var xdomain = [-2, 2];
    var ydomain = [-2, 2];

    this.margin = {top: 20, right: 20, bottom: 15, left: 40};
    this.width = $(parent).width() - this.margin.right - this.margin.left;
    this.height = $(parent).height() - this.margin.top - this.margin.bottom;

    if(this.colorbar_enabled){
        this.height -= colorbar_height;
        this.height -= 15;
    }

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
        .attr("width", self.width + self.margin.left + self.margin.right)
        .attr("height", self.height + self.margin.top + self.margin.bottom)
        .append("g")
        .attr("transform", "translate(" + self.margin.left + "," + self.margin.top + ")")
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
        this.colorbar_svg = d3.select(parent).append("svg")
            .attr("width", this.width + this.margin.left + this.margin.right)
            .attr("height", colorbar_height+20);

        this.colorbar = this.colorbar_svg.append("rect")
            .attr("x", this.width + this.margin.left + this.margin.right - colorbar_width - 50)
            .attr("y", 0)
            .attr("width", colorbar_width)
            .attr("height", colorbar_height);
    }

    this.points = [];

    this.selected = -1;
    this.selected_links = [];

    this.hovered = -1;
    this.hovered_links = [];

    this.last_event = -1;
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
    }
    else
    {
        cvals.sort(d3.ascending); // Needed for quantile
        var low = d3.quantile(cvals, 0.1);
        var high = d3.quantile(cvals, 0.9);
        var mid = d3.mean(cvals);

        this.colorScale = d3.scale.linear()
            .domain([low, mid, high])
            .range(["blue", "green", "red"]);

        // Format the bound labels
        var label_low = parseFloat(low.toFixed(2)).toString();
        var label_high = parseFloat(high.toFixed(2)).toString();

        this.setColorBar(this.colorScale.range(),
               label_low,
               label_high);
    }

    this.redraw(true)();
};

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
            parseInt(this.colorbar.attr("width"))).toString();

        this.colorbar_svg.append("text")
            .attr("text-anchor", "middle")
            .attr("x", label_low_x)
            .attr("y", 30)
            .attr("font-size", "10px")
            .text(label_low);

        this.colorbar_svg.append("text")
            .attr("text-anchor", "middle")
            .attr("x", label_high_x)
            .attr("y", 30)
            .attr("font-size", "10px")
            .text(label_high);
    }
};

ColorScatter.prototype.setSelected = function(selected_index, event_id)
{
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
};

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

ColorScatter.prototype.redraw = function(performTransition) {
    var self = this;
    return function(){
    self.svg.select(".x.axis").call(self.xAxis);
    self.svg.select(".y.axis").call(self.yAxis);
    var circles = self.svg.selectAll("circle")
        .data(self.points);

    circles.enter().append("circle").attr("r",4.0);
    circles.style("fill", function(d){return self.colorScale(d[2]);})
        .on("click", function(d,i){self.setSelected(i);})
        .on("mouseover", function(d,i){self.tip.show(d,i); self.setHovered(i);})
        .on("mouseout", function(d,i){self.tip.hide(d,i); self.setHovered(-1);})
        .classed("point-selected", function(d, i){return i === self.selected;});

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
