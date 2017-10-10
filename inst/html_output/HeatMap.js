/**
 * Created by David DeTomaso on 6/24/2015.
 */

function HeatMap(parent)
{
    var self = this;
    this.h = 1;  //height of row
    this.w = 4;  //width of column

    this.width = 600;
    this.height = 450;

    this.svg = d3.select(parent).append("svg")
        .attr("width", self.width)
        .attr("height", self.height);

    var offset = 0;
        
    this.labels = this.svg.append("g");
    
    this.labels.append("text")
        .classed("col_label", true)
        .attr("x", 0)
        .attr("y", offset + 20)
        .attr("font-size", "20px");
        
    this.labels.append("text")
        .classed("row_label", true)
        .attr("x", 0)
        .attr("y", offset + 40)
        .attr("font-size", "20px");

    this.labels.append("rect")
        .attr("x", this.width-40)
        .attr("y", offset)
        .attr("width", 40)
        .attr("height", 40)
        .style("fill", "white");

    this.labels.append("text")
        .classed("rect_label", true)
        .attr("text-anchor", "middle")
        .attr("x", this.width-20)
        .attr("y", offset+23)
        .attr("font-size", "10px");
        
    offset += 40;
    offset += 10; //Some margin
        
    
    this.cluster_bar = this.svg.append("g");
    this.cluster_bar.on("mouseleave", function(){
        self.setHovered(-1);
    });
    this.cluster_bar
        .attr("transform", "translate(0," +
        (offset)+")");

    this.cluster_bar_props = {
        height: 10,
        colors: d3.scale.category10().domain(d3.range(10))
    };
    
    offset += this.cluster_bar_props.height;
    offset += 10; //Some margin

    this.grid = this.svg.append("g");
    this.grid
        .attr("transform", "translate(0," +
        (offset)+")");
    
    offset += this.height;
    this.svg.attr("height", offset);

    //define a color scale using the min and max expression values
    this.colorScale = d3.scale.linear()
        .domain([-0.6, 0, 0.6])
        .range(["steelblue", "white", "lightcoral"]);

    this.data = [];
    this.selected = -1;
    this.selected_links = [];

    this.hover_cols = -1;
    this.hovered_links = [];
    
    this.hover_rows = -1;

    this.last_event = 0;

    this.col_clusters = null; //Cluster assignment for each column in data matrix
    this.col_order = null;  //Position for each column in data matrix
    
    this.row_labels = [];
    this.col_labels = [];

}

HeatMap.prototype.cluster_columns = function(assignments)
{
    var self = this;
    this.col_clusters = assignments;

    //Argsort the col_clusters
    var rr = d3.range(0,this.col_clusters.length);
    rr.sort(function(a,b){return self.col_clusters[a] - self.col_clusters[b];});
    //Argsort again to get rank
    var rr2 = d3.range(0,rr.length);
    rr2.sort(function(a,b){return rr[a] - rr[b];});
    
    this.col_order = rr2;
    this.redraw()();
};

HeatMap.prototype.setData = function(data, render)
{
    var formatted_data = data.map(function(row,row_i){
        return row.map(function(value, col_i){
            return {"row":row_i, "col":col_i, "value":value};
        });
    });
            
    this.data = formatted_data;
    var N_ROWS = data.length;
    var N_COLS = data[0].length;

    this.w = Math.floor(this.width/N_COLS);
    if(this.w === 0) {this.w = 1;}
    this.h = Math.floor(this.height/N_ROWS);
    if(this.h === 0) {this.h = 1;}

    this.col_clusters = Array.apply(null, Array(N_COLS)).map(Number.prototype.valueOf,0);
    this.col_order = d3.range(0, N_COLS);  //Initial sorting
    
    if(render){
        this.redraw(true)();
    }
};

HeatMap.prototype.setSelected = function(selected_index, event_id)
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

HeatMap.prototype.setHovered = function(hovered_indices, event_id)
{
    if(event_id === undefined){
        event_id = Math.random();
    }

    //test for single index, and wrap in list
    if(typeof(hovered_indices) === "number"){hovered_indices = [hovered_indices];}

    //Needed to prevent infinite loops with linked hover and select events
    if(this.last_event !== event_id) {
        this.last_event = event_id;
        this.hover_cols = hovered_indices;
        this.grid.selectAll("g").selectAll("rect")
            .classed("heatmap-hover", function (d, i) {
                return hovered_indices.indexOf(i) > -1;
            });
        this.hovered_links.forEach(function (e) {
            e.setHovered(hovered_indices, event_id);
        });
        if(this.hover_cols.length === 1)
        {
            var ii = this.hover_cols[0];
            this.labels.select(".col_label").text(this.col_labels[ii]);
        }
    }
};

HeatMap.prototype.setHoveredRow = function(hovered_row_indices)
{
    if(typeof(hovered_row_indices) === "number"){
        hovered_row_indices = [hovered_row_indices];
    }
    
    this.hover_rows = hovered_row_indices;
    
    if(this.hover_rows.length === 1)
    {
        var ii = this.hover_rows[0];
        this.labels.select(".row_label").text(this.row_labels[ii]);
    }
    
};

//Sets the text and color for the square upper-right indicator
HeatMap.prototype.setHoveredIndicator = function(data_val)
{
    if(data_val !== undefined)
    {
        this.labels.select("rect")
            .style("fill", this.colorScale(data_val));
        this.labels.select(".rect_label")
            .text(data_val);
    }
    else
    {
        this.labels.select("rect")
            .style("fill", "white");
        this.labels.select(".rect_label")
            .text("");
    }

};

HeatMap.prototype.redraw = function(performTransition) {
    var self = this;
    return function(){
        //self.svg.select(".x.axis").call(self.xAxis);
        //self.svg.select(".y.axis").call(self.yAxis);

        //Draw cluster-bar
        var clusterRects = self.cluster_bar.selectAll("rect")
            .data(self.col_clusters);

        clusterRects.enter()
            .append("rect")
            .attr('height', self.cluster_bar_props.height)
            .attr('y', self.cluster_bar_props.offset);

        clusterRects
            .attr('width', self.w)
            .style('fill',function(d) {
                return self.cluster_bar_props.colors(d);})
            .attr('x', function(d,i){
                return (self.col_order[i] * self.w); })
            .on("mouseover", function(d) {
                var ii = d3.range(self.col_clusters.length);
                var selected_i = ii.filter(function(e,j){
                    return self.col_clusters[j] === d;});
                self.setHovered(selected_i);
                });

        clusterRects.exit().remove();


        //generate heatmap rows
        var heatmapRow = self.grid.selectAll("g")
            .data(self.data);

        heatmapRow.enter()
            .append("g");

        //heatmapRow.attr("transform", function(d, j){
        //        return "translate(0," + (j*self.h)+ ")"});
                
        heatmapRow.exit().remove();

        //generate heatmap columns
        var heatmapRects = heatmapRow
            .selectAll("rect")
            .data(function(d) {
                return d;
            });

        heatmapRects.enter().append("rect")
            .on("mouseover", function(d){self.setHovered(d.col); self.setHoveredRow(d.row); self.setHoveredIndicator(d.value);});

        self.svg
            .on("mouseleave", function(d){self.setHovered(-1); self.setHoveredRow(-1); self.setHoveredIndicator();});

        heatmapRects.style('fill',function(d) {
                return self.colorScale(d.value);})
            .attr('width',self.w)
            .attr('height',self.h)
            .attr('y',function(d){ return d.row*self.h;})
            .attr('x', function(d) {
                return (self.col_order[d.col] * self.w);});

        heatmapRects.exit().remove();

    };
};
