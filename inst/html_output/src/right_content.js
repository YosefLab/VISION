function Right_Content()
{
    this.scatter = {}
}

Right_Content.prototype.init = function()
{
    var self = this;
    self.dom_node = $("#right-content");
    self.scatter = new ColorScatter("#scatter-div", true, true);

    self.scatterColorOptions = $(self.dom_node).find("input[name='scatterColorButtons']")

    $(self.scatterColorOptions).on('change', function(){
        self.update({ 'colorScatterOption': '' }) // No need to send value
    })

    //Enable Toggling of Lasso Select
    $(self.dom_node).find("#lasso-select").on("click", function() {
        var tog = $(this).html();
        if (tog == "Enable Lasso Select") {
            self.scatter.toggleLasso(true);
            $(this).html("Disable Lasso Select");
        } else {
            self.scatter.toggleLasso(false);
            $(this).html("Enable Lasso Select");
        }
    });

    $(self.dom_node)
        .find("#export-button")
        .on("click", function () {
            self.exportSigProj()
        });


}

Right_Content.prototype.update = function(updates)
{
    var self = this;
    
    var needsUpdate = ('main_vis' in updates) ||
        ('plotted_item' in updates) ||
        ('plotted_item_type' in updates) ||
        ('plotted_projection' in updates) ||
        ('filter_group' in updates) ||
        ('plotted_pc' in updates) ||
        ('colorScatterOption' in updates);

    if (!needsUpdate) return;

    var main_vis = get_global_status('main_vis');

    if(main_vis === 'sigvp'){
        self.draw_sigvp();

    } else if (main_vis === "tree") {
        self.draw_tree();

    } else if (main_vis === "pcannotator") {
        self.draw_pca();

    } else {
        throw "Bad main_vis value!";
    }

}

Right_Content.prototype.draw_sigvp = function() {

    var self = this;

    $(self.dom_node).find("#plotted-value-option").show()

    var item_key = get_global_status('plotted_item');
    var item_type = get_global_status('plotted_item_type');
    var proj_key = get_global_status('plotted_projection');

    var projection = get_global_data('sig_projection_coordinates')
    var values = get_global_data('plotted_values')

    var isFactor;
    if(item_type !== "gene"){
        isFactor = get_global_data('sig_info').isFactor
    } else {
        isFactor = false;
    }

    if(self.getScatterColorOption() == "rank" && !isFactor){
        values = self.rank_values(values)
    }


    $('#plot-title').text(proj_key);
    $('#plot-subtitle').text(item_key);

    var points = [];
    var sample_labels = Object.keys(values).sort()

    _.each(sample_labels, (sample_label) => {
        var x = projection[sample_label][0]
        var y = projection[sample_label][1]
        var sig_score = values[sample_label]
        points.push([x, y, sig_score, sample_label]);
    })

    self.scatter.setData(points, isFactor);

}


Right_Content.prototype.draw_tree = function() {

    var self = this;

    $(self.dom_node).find("#plotted-value-option").show()

    var item_key = get_global_status('plotted_item');
    var item_type = get_global_status('plotted_item_type');
    var proj_key = get_global_status('plotted_projection');
    var filter_group = get_global_status('filter_group');

    var tree_points = api.tree.tree_points(filter_group, proj_key);
    var tree_adjlist = api.tree.tree(filter_group)

    var projection = get_global_data('tree_projection_coordinates')
    var values = get_global_data('plotted_values')

    var isFactor;
    if(item_type !== "gene"){
        isFactor = get_global_data('sig_info').isFactor
    } else {
        isFactor = false;
    }

    if(self.getScatterColorOption() == "rank" && !isFactor){
        values = self.rank_values(values)
    }

    return $.when(tree_points, tree_adjlist) // Runs when both are completed
        .then(function(treep, treel){

            // Massage treep for easier D3 binding

            tree_points = []

            $('#plot-title').text(proj_key);
            $('#plot-subtitle').text(item_key);

            var points = [];
            var sample_labels = Object.keys(values).sort()

            _.each(sample_labels, (sample_label) => {
                var x = projection[sample_label][0]
                var y = projection[sample_label][1]
                var sig_score = values[sample_label]
                points.push([x, y, sig_score, sample_label]);
            })

            var tree_points = [];
            for (var i = 0; i < treep[0].length; i++) {
                var x = treep[0][i];
                var y = treep[1][i];
                tree_points.push([x, y, "Node " + i]);
            }

            // Change tree adjacency list into a list of pairs
            var tree_adj = []

            for (var i = 0; i < treel.length; i++) {
                for (var j = i+1; j < treel[i].length; j++) {
                    if (treel[i][j] == 1) {
                        tree_adj.push([i, j])
                    }
                }
            }

            self.scatter.setData(points, isFactor, tree_points, tree_adj);

        });
}

Right_Content.prototype.draw_pca = function() {

    var self = this;

    $(self.dom_node).find("#plotted-value-option").hide()

    var item_key = get_global_status('plotted_item');
    var item_type = get_global_status('plotted_item_type');

    var pc_key = get_global_status('plotted_pc');

    var pca = get_global_data('pca_projection_coordinates')
    var values = get_global_data('plotted_values')

    var isFactor;
    if(item_type === "gene"){
        isFactor = get_global_data('sigIsPrecomputed')[item_key];
    } else {
        isFactor = false;
    }

    $("#plot-title").text("PC: ".concat(pc_key));
    $("#plot-subtitle").text(item_key);

    var points = []
    var sample_labels = Object.keys(values).sort()

    _.each(sample_labels, (sample_label) => {
        var x = pca[sample_label][pc_key-1]
        var y = values[sample_label]
        var sig_score = null
        points.push([x, y, sig_score, sample_label]);
    })

    self.scatter.setData(points, isFactor);

}

// Called when the window is resized
Right_Content.prototype.resize = function() {

    $('#scatter-div').children().remove();
    this.scatter = new ColorScatter("#scatter-div", true, true);
    this.update({'plotted_item': ''}) // Tricks into replotting

}

Right_Content.prototype.getScatterColorOption = function() {
    return this.scatterColorOptions.filter(":checked").val()
}


//Returns the rank of each value in values (object)
//Averages ranks that are ties
Right_Content.prototype.rank_values = function(values)
{
    var pairs = _.toPairs(values)
    pairs.sort(function(a, b) { return a[1] - b[1];})
    var ranks = {}

    var current_group_start = 0
    var current_group_end;
    var last_val = pairs[0][1]
    var current_val;
    var current_group_rank;
    var i, j;

    for(i = 1; i < pairs.length; i++)
    {
        current_val = pairs[i][1]
        if(current_val !== last_val){
            current_group_end = i-1;
            current_group_rank = (current_group_end + current_group_start)/2
            for(j = current_group_start; j <= current_group_end; j++)
            {
                ranks[pairs[j][0]] =  current_group_rank
            }
            current_group_start = i
        }
        last_val = current_val
    }

    // Need to wrap up the final group
    current_group_end = pairs.length-1;
    current_group_rank = (current_group_end + current_group_start)/2
    for(j = current_group_start; j <= current_group_end; j++)
    {
        ranks[pairs[j][0]] =  current_group_rank
    }

    return ranks
}

/*
Exports a zip with data in it
 */
Right_Content.prototype.exportSigProj = function()
{
    var self = this;
    var zip = new JSZip();

    var main_vis = get_global_status('main_vis')

    var type = get_global_status('plotted_item_type')
    var plotted_item = get_global_status('plotted_item')
    var values = get_global_data('plotted_values')

    //Convert the data that's in the scatter plot to a tab-delimited table

    var proj;
    if (main_vis === 'sigvp') {
        proj = get_global_data('sig_projection_coordinates')
    } else if (main_vis ==='pcannotator') {
        proj = get_global_data('pca_projection_coordinates')
    } else if (main_vis ==='tree') {
        proj = get_global_data('tree_projection_coordinates')
    }

    var table;
    if (main_vis === 'sigvp' || main_vis === 'tree') {

        table = _.map(proj, (value, key) => {
            return [key, proj[key][0], proj[key][1], values[key]]
        });

        table = [["Cell", "X", "Y", plotted_item]].concat(table);

    } else if (main_vis ==='pcannotator') {
        var plotted_pc = get_global_status('plotted_pc')
        table = _.map(proj, (value, key) => {
            return [key, proj[key][plotted_pc-1], values[key]]
        });

        table = [["Cell", "PC: "+plotted_pc, plotted_item]].concat(table);
    }

    table = table.map(function(x){ return x.join("\t");});
    var scatter_csv_str = table.join("\n");
    zip.file("Scatter.txt", scatter_csv_str);

    //Get the scatter plot and convert to a PNG
    var svg = $(self.dom_node).find('#scatter-div').children('svg');
    svg.attr("version", 1.1)
        .attr("xmlns", "http://www.w3.org/2000/svg");
    var svg2 = svgCopy(svg.get(0));

    var html_data = svg2.parentNode.innerHTML;
    zip.file("Scatter.svg", html_data);

    var imgsrc = "data:image/svg+xml;base64," + btoa(html_data);

    var image = new Image();
    image.onload = function()
    {
        var canvas = document.createElement("canvas");
        canvas.width = image.width;
        canvas.height = image.height;
        var context = canvas.getContext("2d");
        context.drawImage(image, 0,0);

        var canvasdata = canvas.toDataURL("image/png");
        //Strip off the data URI portion
        var scatter_png = canvasdata.substring(canvasdata.indexOf(",")+1);

        //Take the result and stick it into a zip

        zip.file("Scatter.png", scatter_png   , {base64: true});

        var zip_uri = "data:application/zip;base64," + zip.generate({type:"base64"});

        var a = document.createElement("a");
        var proj_name;
        if(main_vis === 'pcannotator'){
            proj_name = 'PC' + get_global_status('plotted_pc')
        } else {
            proj_name = get_global_status('plotted_projection')
        }

        a.download = plotted_item+"_"+proj_name+".zip";
        a.href = zip_uri;
        a.click();

    };

    image.setAttribute("src", imgsrc)

}
