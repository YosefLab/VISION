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
    this.init()
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
