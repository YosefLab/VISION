/**
 * Created by David DeTomaso on 6/24/2015.
 */

function HeatMap(parent)
{
    var self = this;
    this.h = 1;  //height of row

    this.width = $(parent).width();
    this.height = $(parent).height();

    var otherHeight = 0;
    //Subtract height of anything else
    $(parent).children().each(function(i,e){ otherHeight += $(e).outerHeight(true);});
    this.height -= otherHeight;

    this.svg = d3.select(parent).append("svg")
        .attr("width", self.width)
        .attr("height", self.height);

    var offset = 0;
        
    this.labels = this.svg.append("g");
        
    this.labels.append("text")
        .classed("row_label", true)
        .attr("x", 0)
        .attr("y", offset + 28)
        .attr("font-size", "20px");

    this.labels.append("rect")
        .attr("x", this.width-40)
        .attr("y", offset)
        .attr("width", 40)
        .attr("height", 30)
        .style("fill", "white");

    this.labels.append("text")
        .classed("rect_label", true)
        .attr("text-anchor", "middle")
        .attr("x", this.width-20)
        .attr("y", offset+17)
        .attr("font-size", "10px");
        
    offset += 30;
    offset += 10; //Some margin
        
    this.xlabel_margin = 30;
    this.heat_height = this.height - offset - this.xlabel_margin;
    this.grid_gap = 10; //# of pixels between the plus and minus heat maps
    this.grid_start = offset;
    this.grid_xoffset = 50; //# of pixels to the left of the grids for the + and - labels
    this.grid_xoffsetr = 50; //# of pixels to the right of the grids for the ylabel

    this.grid_plus = this.svg.append("g");
    this.grid_minus = this.svg.append("g");

    this.grid_plus_label = this.svg.append("text")
        .attr("text-anchor", "middle")
        .attr("x", this.grid_xoffset/2)
        .attr("y", 0)
        .attr("font-size", "25px")
        .text("");

    this.grid_minus_label = this.svg.append("text")
        .attr("text-anchor", "middle")
        .attr("x", this.grid_xoffset/2)
        .attr("y", 0)
        .attr("font-size", "25px")
        .text("");

    this.grid_x_label = this.svg.append("text")
        .attr("text-anchor", "middle")
        .attr("x", this.grid_xoffset+(this.width - this.grid_xoffsetr)/2)
        .attr("y", this.grid_start + this.heat_height + this.xlabel_margin*.8)
        .attr("font-size", "25px")
        .text("Sample Clusters");

    this.grid_y_label = this.svg.append("text")
        .attr("text-anchor", "middle")
        .attr("x", this.width - this.grid_xoffsetr/2)
        .attr("y", this.height/2)
        .attr("font-size", "25px")
        .style("writing-mode", "tb")
        .text("Genes");
    

    //define a color scale using the min and max expression values
    this.colorScale = d3.scale.linear()
        .domain([0, 4])
        .range(["white", "lightcoral"]);

    this.data_plus = [];
    this.data_minus = [];
    this.selected = -1;
    this.selected_links = [];

    this.hover_cols = -1;
    this.hovered_links = [];
    
    this.hovered_row = "";

    this.last_event = 0;

    this.col_labels = [];

    this.cluster_assignments = []; //Used in hoverCol.  Denotes which cluster each sample is assigned to
}

HeatMap.prototype.setLabels = function(xlabel, ylabel)
{
    this.grid_x_label.text(xlabel);
    this.grid_y_label.text(ylabel);
};


function dist_mat(data)
{
    //creates a matrix of size data.length x data.length
    //distance between entries(rows) of data
    //each item in data (which is an array) is an array of numbers
    
    //Create the array
    var output = [];
    for(var i = 0; i < data.length; i++)
    {
        output.push([]);
    }

    //Calculate distances
    for(var i = 0; i < data.length; i++)
    {
        for(var j = i; j < data.length; j++)
        {
            var dd = dist(data[i], data[j]);
            output[i][j] = dd;
            output[j][i] = dd;
        }
    }

    return output;
}

function dist(x,y)
{
    //x and y are arrays of numbers of the same length
    
    //Euclidean
    var total = 0;
    for(var i=0; i < x.length; i++)
    {
        var temp = x[i] - y[i];
        total = total + temp * temp;
    }
    return Math.sqrt(total);

    //1-Pearsons
    var stdx = d3.deviation(x);
    var stdy = d3.deviation(y);
    var meanx = d3.mean(x);
    var meany = d3.mean(y);
    var total = 0;
    for(var i = 0; i < x.length; i++)
    {
        total = total + (x[i] - meanx) * (y[i] - meany);
    }
    total = total / (stdx * stdy);
    return 1-total;
}

function cluster_order(data)
{

    var distances = dist_mat(data);
    var clusters = [];

    //Convert array of rows into array of objects
    //row_data = row vector
    //index = index in the original matrix
    
    for(var i=0; i < data.length; i++)
    {
        clusters.push([{'index':i, 'row_data':data[i]}]);
    }

    //clusters is a list of clusters
    //each cluster is a list of row vectors
    
    while(clusters.length > 1)
    {
        clusters = merge_closest_cluster_pair(clusters, distances);
    }

    var leaves_indices = [];
    for(var i = 0; i < clusters[0].length; i++)
    {
        leaves_indices.push(clusters[0][i].index);
    }

    return leaves_indices;
}

function merge_closest_cluster_pair(clusters, distances)
{
    var best_pair = [0,1];
    var closest_distance = 1e99;

    for(var i = 0; i < clusters.length; i++)
    {
        for(var j = i+1; j < clusters.length; j++)
        {
            var distance = distance_between_clusters(clusters[i], clusters[j], distances);
            if(distance < closest_distance)
            {
                closest_distance = distance;
                best_pair = [i,j];
            }
        }
    }

    var clusterA = clusters[best_pair[0]];
    var clusterB = clusters[best_pair[1]];
    var merged_cluster = merge_clusters(clusterA, clusterB, distances);

    var new_clusters = [];
    for(var i = 0; i < clusters.length; i++)
    {
        if(i != best_pair[0] && i != best_pair[1])
        {
            new_clusters.push(clusters[i]);
        }
    }

    new_clusters.push(merged_cluster);
    return new_clusters;
}

function merge_clusters(clusterA, clusterB, distances)
{
    //Merge can happen in two orders
    //Compare first/last clusters in either order to determine the best order
    
    //AB
    var a_last = clusterA[clusterA.length-1];
    var b_first = clusterB[0];
    var distAB = distances[a_last.index][b_first.index];

    //BA
    var a_first = clusterA[0];
    var b_last = clusterB[clusterB.length-1];
    var distBA = distances[a_first.index][b_last.index];

    var new_cluster;
    if(distBA < distAB)
    {
        new_cluster = clusterB.concat(clusterA);
    }
    else
    {
        new_cluster = clusterA.concat(clusterB);
    }

    return new_cluster;
}

function distance_between_clusters(clusterA, clusterB, distances)
{
    //Compute the average distance
    var total_dist = 0;
    for(var i = 0; i < clusterA.length; i++)
    {
        for(var j = 0; j < clusterB.length; j++)
        {
            var rowA = clusterA[i];
            var rowB = clusterB[j];
            total_dist = total_dist + distances[rowA.index][rowB.index];
        }
    }
    total_dist = total_dist / clusterA.length / clusterB.length;
    return total_dist;
    
//    var min_dist = 1e99;
//    for(var i = 0; i < clusterA.length; i++)
//    {
//        for(var j = 0; j < clusterB.length; j++)
//        {
//            var rowA = clusterA[i];
//            var rowB = clusterB[j];
//            var curr_dist = distances[rowA.index][rowB.index];
//            if(curr_dist < min_dist)
//            {
//                min_dist = curr_dist;
//            }
//        }
//    }
//
//    return min_dist;
}

function order_cluster_group_rows(cluster_group)
{
    if(cluster_group[0].data.length === 0)
    {
        return;
    }
    //create data matrix from cluster_group
    var data = [];
    for(var i = 0; i < cluster_group.length; i++)
    {
        var col_obj = cluster_group[i];
        var row = col_obj.data.map(function(x){return x.value;});
        data.push(row);
    }

    data = d3.transpose(data);

    // Get optimal leaf orderings from heirarchical clustering
    var leaf_order = cluster_order(data);

    // Reorder the cluster_group in place
    for(var i = 0; i < cluster_group.length; i++)
    {
        var col_obj = cluster_group[i];
        var old_data = col_obj.data;
        var new_data = [];
        for(var j = 0; j < leaf_order.length; j++)
        {
            new_data.push(old_data[leaf_order[j]]);
        }
        col_obj.data = new_data;
    }
}


function order_cluster_group_cols(self, cluster_group_plus, cluster_group_minus, TOTAL_SAMPLES)
{
    // Cluster the columns in the heatmap
    // First, merge the + and - sets into a single column
    // Then, use that column to calculate a new column order
    // Lastly, use the new column order to modify the cluster_group_plus(minus) objects
    
    if(cluster_group_plus[0].data.length === 0 && cluster_group_minus[0].data.length === 0)
    {
        return [cluster_group_plus, cluster_group_minus];
    }
    //create data matrix from cluster_group
    var data = [];
    for(var i = 0; i < cluster_group_plus.length; i++)
    {
        var col_obj_p = cluster_group_plus[i];
        var col_obj_m = cluster_group_minus[i];
        var row = col_obj_p.data.map(function(x){return x.value;});
        row = row.concat(col_obj_m.data.map(function(x){return x.value;}));
        data.push(row);
    }

    // Get optimal leaf orderings from heirarchical clustering
    var leaf_order = cluster_order(data);

    cluster_group_plus = leaf_order.map(function(e){return cluster_group_plus[e];});
    cluster_group_minus = leaf_order.map(function(e){return cluster_group_minus[e];});


    // Change the x_offset for each cluster
    var x_offset = self.grid_xoffset;
    for (var j = 0; j < cluster_group_plus.length; j++)
    {
        var clust_p = cluster_group_plus[j];
        var width = clust_p.weight / TOTAL_SAMPLES * (self.width - self.grid_xoffset - self.grid_xoffsetr);

        for (var k = 0; k < clust_p.data.length; k++)
        {
            clust_p.data[k].x = x_offset;
        }
        
        x_offset = x_offset + width;
    }

    // Same thing for the minus clusters
    var x_offset = self.grid_xoffset;
    for (var j = 0; j < cluster_group_minus.length; j++)
    {
        var clust_m = cluster_group_minus[j];
        var width = clust_m.weight / TOTAL_SAMPLES * (self.width - self.grid_xoffset - self.grid_xoffsetr);

        for (var k = 0; k < clust_m.data.length; k++)
        {
            clust_m.data[k].x = x_offset;
        }
        
        x_offset = x_offset + width;
    }

    return [cluster_group_plus, cluster_group_minus];

}


HeatMap.prototype.setData = function(data, cluster_assignments, gene_labels, gene_signs, sample_labels)
{
    //Data is an array of rows, each containing an array of values for each col
    //cluster_assignments is an array of numbers indicating assignment
    
    this.col_labels = sample_labels;

    this.cluster_assignments = cluster_assignments;
    var dataT = d3.transpose(data);
    var TOTAL_SAMPLES = cluster_assignments.length;
    var clusters = {};
    for(var j = 0; j < cluster_assignments.length; j++)
    {
        if(clusters[cluster_assignments[j]] === undefined){
            clusters[cluster_assignments[j]] = {'weight':0, 'data':[]};
        }

        var clust = clusters[cluster_assignments[j]];
        clust.weight = clust.weight + 1;
        clust.data.push(dataT[j]);
    }

    //Now clusters is a dict of "cluster index" -> ('weight' -> # of cols in cluster, 'data' -> N_cols x N_row for cluster)
    for (var clust_i in clusters)
    {
        var clust = clusters[clust_i];
        var cdT = d3.transpose(clust.data);
        clust.data = cdT.map(function(e){return d3.mean(e);});
    }
    //Now, each clusters data is just a length N_genes array of values for that cluster.

    //Convert to a list for easier d3 binding
    var cluster_list = [];
    for (var clust_i in clusters)
    {
        var clust = clusters[clust_i];
        clust.index = parseInt(clust_i);
        cluster_list.push(clust);
    }

    //Add in information on the width of each cluster
    var x_offset = this.grid_xoffset;
    for (var j = 0; j < cluster_list.length; j++)
    {
        var clust = cluster_list[j];
        var width = clust.weight / TOTAL_SAMPLES * (this.width - this.grid_xoffset - this.grid_xoffsetr);
        
        clust.data = clust.data.map(function(e,i){
            return {"value":e, "x":x_offset, "width":width, "index": clust.index, "gene": gene_labels[i]};
        });

        x_offset = x_offset + width;
    }

    var N_ROWS = cluster_list[0].data.length;

    //this.h = Math.floor((this.heat_height - this.grid_gap)/N_ROWS);
    this.h = (this.heat_height - this.grid_gap)/N_ROWS;
    if(this.h === 0) {this.h = (this.heat_height - this.grid_gap)/N_ROWS;}

    //Split data into data_plus and data_minus
    //Unsigned sigs (sign = 0) go in data_plus

    var cluster_list_plus = [];
    var cluster_list_minus = [];
    for(var i = 0; i < cluster_list.length; i++)
    {
        var clust = cluster_list[i];

        //get positive rows
        var clust_plus = {};
        clust_plus.data = clust.data.filter(function(e,i){
            return gene_signs[i] === 0 || gene_signs[i] === 1;
        });

        //get negative rows
        var clust_minus = {};
        clust_minus.data = clust.data.filter(function(e,i){
            return gene_signs[i] === -1;
        });

        //copy over all other properties
        for(var key in clust){
            if(key !== "data"){
                clust_plus[key] = clust[key];
                clust_minus[key] = clust[key];
            }
        }

        cluster_list_plus.push(clust_plus);
        cluster_list_minus.push(clust_minus);

    }

    // Cluster within data_minus and data_plus
    
    order_cluster_group_rows(cluster_list_plus);
    order_cluster_group_rows(cluster_list_minus);

    var out = order_cluster_group_cols(this, cluster_list_plus, cluster_list_minus, TOTAL_SAMPLES);
    cluster_list_plus = out[0];
    cluster_list_minus = out[1];

    this.data_plus = cluster_list_plus;
    this.data_minus = cluster_list_minus;

    this.setColormap();

    this.redraw()();
};

HeatMap.prototype.setColormap = function()
{
    var self = this;
    var vals_p = self.data_plus.map(function(x){
        return x.data.map(function(y){
            return y.value;}
            );
    });
    var vals_m = self.data_minus.map(function(x){
        return x.data.map(function(y){
            return y.value;}
            );
    });
    var vals = vals_p.concat(vals_m);

    vals = d3.merge(vals);
    vals = vals.sort(d3.ascending);
    var upper = d3.quantile(vals, 0.99);
    var low = d3.quantile(vals, 0.01);
    var mid = (low + upper)/2;

    self.colorScale = d3.scale.linear()
        .domain([low, mid, upper])
        .range(["steelblue", "white", "lightcoral"]);

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

};

HeatMap.prototype.setHoveredCol = function(hovered_col_index)
{
    var hovered_indices = [];
    if(hovered_col_index !== undefined) {
        for (var i = 0; i < this.cluster_assignments.length; i++) {
            if (this.cluster_assignments[i] === hovered_col_index) {
                hovered_indices.push(i);
            }
        }
    }

    this.hovered_links.forEach(function (e) {
        e.setHovered(hovered_indices);
    });
};

HeatMap.prototype.setHoveredRowLabel = function(row_label)
{

    this.hovered_row = row_label;
    this.labels.select(".row_label").text(this.hovered_row);
    
};

//Sets the text and color for the square upper-right indicator
HeatMap.prototype.setHoveredIndicator = function(data_val)
{
    if(data_val !== undefined)
    {
        this.labels.select("rect")
            .style("fill", this.colorScale(data_val));
        this.labels.select(".rect_label")
            .text(data_val.toFixed(3));
    }
    else
    {
        this.labels.select("rect")
            .style("fill", "white");
        this.labels.select(".rect_label")
            .text("");
    }

};

HeatMap.prototype.redraw = function() {
    var self = this;
    return function(){

        //generate heatmap columns for plus grid
        
        var pos_grid_start = self.grid_start;
        self.grid_plus.attr("transform",
                "translate(0," + (pos_grid_start)+")"
                );

        var heatmapCols_plus = self.grid_plus.selectAll("g")
            .data(self.data_plus);

        heatmapCols_plus.enter()
            .append("g");

        heatmapCols_plus.exit().remove();

        //generate heatmap rows
        var heatmaRects_plus = heatmapCols_plus
            .selectAll("rect")
            .data(function(d) {
                return d.data;
            });

        heatmaRects_plus.enter().append("rect")
            .on("mouseover", function(d,i){ self.setHoveredRowLabel(d.gene); self.setHoveredCol(d.index); self.setHoveredIndicator(d.value);});

        heatmaRects_plus.style('fill',function(d) {
                return self.colorScale(d.value);})
            .attr('width',function(d){return d.width;})
            .attr('height',self.h)
            .attr('y',function(d,i){ return i*self.h;})
            .attr('x', function(d) {
                return d.x;});

        heatmaRects_plus.exit().remove();

        var num_pos_rects;
        if(self.data_plus.length === 0)
        {
            num_pos_rects = 0;
        } else {
            num_pos_rects = self.data_plus[0].data.length;
        }

        var neg_grid_start = 0;
        neg_grid_start += pos_grid_start; // Start location for all grids
        neg_grid_start += self.h*num_pos_rects; // Offset for space taken up by positive grid
        if(num_pos_rects > 0)
        {
            neg_grid_start += self.grid_gap;  // If there's a positive grid, introduce a gap between them
        }

        self.grid_minus.attr("transform",
                "translate(0," + (neg_grid_start) + ")"
                );

        //generate heatmap columns for minus grid
        var heatmapCols_minus = self.grid_minus.selectAll("g")
            .data(self.data_minus);

        heatmapCols_minus.enter()
            .append("g");

        heatmapCols_minus.exit().remove();

        //generate heatmap rows
        var heatmapRects_minus = heatmapCols_minus
            .selectAll("rect")
            .data(function(d) {
                return d.data;
            });

        heatmapRects_minus.enter().append("rect")
            .on("mouseover", function(d,i){ self.setHoveredRowLabel(d.gene); self.setHoveredCol(d.index); self.setHoveredIndicator(d.value);});

        self.svg
            .on("mouseleave", function(d){ self.setHoveredRowLabel(""); self.setHoveredCol(); self.setHoveredIndicator();});

        heatmapRects_minus.style('fill',function(d) {
                return self.colorScale(d.value);})
            .attr('width',function(d){return d.width;})
            .attr('height',self.h)
            .attr('y',function(d,i){ return i*self.h;})
            .attr('x', function(d) {
                return d.x;});

        heatmapRects_minus.exit().remove();

        var num_neg_rects;
        if(self.data_minus.length === 0)
        {
            num_neg_rects = 0;
        } else {
            num_neg_rects = self.data_minus[0].data.length;
        }


        //Plus and minus labels for the grids
        if(num_pos_rects > 0){
            var pos_grid_center = Math.floor(pos_grid_start + (self.h*num_pos_rects) / 2);
            self.grid_plus_label
                .attr("y", pos_grid_center+9) // Y lines up with baseline, need to offset to vertically center
                .text("+");                   // Some browsers have an attribute to do this, doesn't work on FireFox
        }
        else
        {
            self.grid_plus_label.text("");
        }

        if(num_neg_rects > 0){
            var neg_grid_center = Math.floor(neg_grid_start + (self.h*num_neg_rects) / 2);
            self.grid_minus_label
                .attr("y", neg_grid_center+9) // See note for plus label above
                .text("-");
        }
        else
        {
            self.grid_minus_label.text("");
        }
    };
};
