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

    self.setLoadingStatus = createLoadingFunction(self.dom_node);

    var proj_promise = api.projections.list()
        .then(function(proj_names) {

            var projSelect = self.dom_node.find('#SelectProjScatter')
            projSelect.children().remove()

            _.each(proj_names, function (proj) {
                projSelect.append(
                    $('<option>', {
                        value: proj,
                        text: proj
                    }));
            });

            projSelect.chosen({
                'width': '110px',
                'disable_search_threshold': 99,
            })
                .off('change')
                .on('change', function () {
                    set_global_status({
                        'plotted_projection':$(this).val(),
                    });
                })
                .trigger('chosen:updated')

        });

    var treeproj_promise = api.tree.list()
        .then(function(proj_names) {

            var projSelect = self.dom_node.find('#SelectTrajectoryProjScatter')
            projSelect.children().remove()

            _.each(proj_names, function (proj) {
                projSelect.append(
                    $('<option>', {
                        value: proj,
                        text: proj
                    }));
            });

            projSelect.chosen({
                'width': '110px',
                'disable_search_threshold': 99,
            })
                .off('change')
                .on('change', function () {
                    set_global_status({
                        'plotted_trajectory':$(this).val(),
                    });
                })
                .trigger('chosen:updated')

        });

    return $.when(proj_promise, treeproj_promise)


}

Right_Content.prototype.update = function(updates)
{
    var self = this;

    var needsUpdate = ('main_vis' in updates) ||
        ('plotted_item' in updates) ||
        ('plotted_item_type' in updates) ||
        ('plotted_projection' in updates) ||
        ('plotted_trajectory' in updates) ||
        ('plotted_pc' in updates) ||
        ('colorScatterOption' in updates) ||
        ('selected_cluster' in updates);

    if (!needsUpdate) return;

    var main_vis = get_global_status('main_vis');

    var autoZoom = false
    if('plotted_projection' in updates ||
       'plotted_trajectory' in updates ||
       'main_vis' in updates) {

        autoZoom = true
    }

    if(main_vis === 'clusters' || main_vis === "pcannotator"){
        self.draw_sigvp(autoZoom);

    } else if (main_vis === "tree") {
        self.draw_tree(autoZoom);

    } else {
        throw "Bad main_vis value!";
    }

    // Update the dropdown if plotted projection changes elsewhere
    if('plotted_projection' in updates) {
        var proj_key = get_global_status('plotted_projection');
        var projSelect = self.dom_node.find('#SelectProjScatter')
        projSelect.val(proj_key)
        projSelect.trigger('chosen:updated') // Changes shown item, but doesn't fire update event
    }

    // Update the dropdown if plotted trajectory changes elsewhere
    if('plotted_trajectory' in updates) {
        var proj_key = get_global_status('plotted_trajectory');
        var projSelect = self.dom_node.find('#SelectTrajectoryProjScatter')
        projSelect.val(proj_key)
        projSelect.trigger('chosen:updated') // Changes shown item, but doesn't fire update event
    }

    if('main_vis' in updates){
        $('#plot-subtitle-latent').hide()
        $('#plot-subtitle-trajectory').hide()
        if(main_vis === 'tree'){
            $('#plot-subtitle-trajectory').show()
        } else {
            $('#plot-subtitle-latent').show()
        }


    }

}

Right_Content.prototype.select_default_proj = function()
{
    var proj = $(this.dom_node.find('#SelectProjScatter')).val()
    var traj = $(this.dom_node.find('#SelectTrajectoryProjScatter')).val()

    var update = {}
    update['plotted_projection'] = proj
    update['plotted_trajectory'] = traj
    return update;
}

Right_Content.prototype.draw_sigvp = function(autoZoom) {

    var self = this;

    var item_key = get_global_status('plotted_item');
    var item_type = get_global_status('plotted_item_type');
    var proj_key = get_global_status('plotted_projection');

    var projection = get_global_data('sig_projection_coordinates')
    var values = get_global_data('plotted_values')

    var isFactor = (typeof(_.values(values)[0]) === 'string') &&
                   (_.values(values)[0] !== "NA")

    var full_color_range, diverging_colormap
    if(item_type === "gene"){
        $(self.dom_node).find("#plotted-value-option").hide()
        full_color_range = true
        diverging_colormap = false
    } else if(item_type === "meta"){
        $(self.dom_node).find("#plotted-value-option").hide()
        full_color_range = false
        diverging_colormap = true
    } else {
        $(self.dom_node).find("#plotted-value-option").show()
        full_color_range = false
        diverging_colormap = true
    }

    if(self.getScatterColorOption() == "rank"
        && !isFactor && item_type !== "gene"
        && item_type !== "meta") {

        values = self.rank_values(values)
    }


    $('#plot-title').text(item_key);

    var points = [];
    var sample_labels = Object.keys(values).sort()

    _.each(sample_labels, (sample_label) => {
        var x = projection[sample_label][0]
        var y = projection[sample_label][1]
        var sig_score = values[sample_label]
        points.push([x, y, sig_score, sample_label]);
    })

    // Get selected cells
    var selected_cluster = get_global_status('selected_cluster')
    var clusters = get_global_data('clusters')

    var selected_cells
    if (selected_cluster !== ''){
        selected_cells = _(clusters)
            .toPairs(clusters)
            .filter(x => x[1] === selected_cluster)
            .map(x => x[0])
            .value()
    } else {
        selected_cells === undefined
    }

    self.scatter.clearData()
    self.scatter.setData(points, isFactor, full_color_range, selected_cells, diverging_colormap);

    if (autoZoom){
        self.scatter.autoZoom();
    } else {
        self.scatter.redraw(true)();
    }
}


Right_Content.prototype.draw_tree = function(autoZoom) {

    var self = this;

    var item_key = get_global_status('plotted_item');
    var item_type = get_global_status('plotted_item_type');
    var proj_key = get_global_status('plotted_trajectory');

    var milestonePromise = api.tree.milestones(proj_key);

    var projection = get_global_data('tree_projection_coordinates')
    var values = get_global_data('plotted_values')

    var isFactor = (typeof(_.values(values)[0]) === 'string') &&
                   (_.values(values)[0] !== "NA")

    var full_color_range, diverging_colormap
    if(item_type === "gene"){
        $(self.dom_node).find("#plotted-value-option").hide()
        full_color_range = true
        diverging_colormap = false
    } else if(item_type === "meta"){
        $(self.dom_node).find("#plotted-value-option").hide()
        full_color_range = false
        diverging_colormap = true
    } else {
        $(self.dom_node).find("#plotted-value-option").show()
        full_color_range = false
        diverging_colormap = true
    }

    if(self.getScatterColorOption() == "rank"
        && !isFactor && item_type !== "gene"
        && item_type !== "meta") {

        values = self.rank_values(values)
    }

    return $.when(milestonePromise) // Runs when both are completed
        .then(function(milestoneCoordinates){

            var treep = milestoneCoordinates[0]
            var treel = milestoneCoordinates[1]

            // Massage treep for easier D3 binding

            tree_points = []

            $('#plot-title').text(item_key);

            var points = [];
            var sample_labels = Object.keys(values).sort()

            _.each(sample_labels, (sample_label) => {
                var x = projection[sample_label][0]
                var y = projection[sample_label][1]
                var sig_score = values[sample_label]
                points.push([x, y, sig_score, sample_label]);
            })

            var tree_points = [];
            for (var i = 0; i < treep.length; i++) {
                var x = treep[i][0];
                var y = treep[i][1];
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

            // Get selected cells
            var selected_cluster = get_global_status('selected_cluster')
            var clusters = get_global_data('clusters')

            var selected_cells
            if (selected_cluster !== ''){
                selected_cells = _(clusters)
                    .toPairs(clusters)
                    .filter(x => x[1] === selected_cluster)
                    .map(x => x[0])
                    .value()
            } else {
                selected_cells === undefined
            }

            self.scatter.clearData()
            self.scatter.setData(points, isFactor, full_color_range, selected_cells, diverging_colormap)
            self.scatter.setTreeData(tree_points, tree_adj)

            if (autoZoom){
                self.scatter.autoZoom();
            } else {
                self.scatter.redraw(true)();
            }

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
        isFactor = (typeof(Object.values(values)[0]) === "string") &&
                   (_.values(values)[0] !== "NA")
    } else {
        isFactor = false;
    }

    $("#plot-title").text(item_key);

    var points = []
    var sample_labels = Object.keys(values).sort()

    _.each(sample_labels, (sample_label) => {
        var x = pca[sample_label][pc_key-1]
        var y = values[sample_label]
        var sig_score = null
        points.push([x, y, sig_score, sample_label]);
    })

    self.scatter.clearData()
    self.scatter.setData(points, isFactor);
    self.scatter.redraw(true)();

}

// Called when the window is resized
Right_Content.prototype.resize = function() {

    $('#scatter-div').children().remove();
    this.scatter = new ColorScatter("#scatter-div", true, true);
    this.update({
        'plotted_item': '',
        'main_vis': '',
    }) // Tricks into replotting

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
    if (main_vis === 'clusters' || main_vis === 'pcannotator') {
        proj = get_global_data('sig_projection_coordinates')
    } else if (main_vis ==='tree') {
        proj = get_global_data('tree_projection_coordinates')
    } else {
        throw "Bad main_vis value!";
    }

    var table;
    if (main_vis === 'tree' || main_vis === 'clusters') {

        table = _.map(proj, (value, key) => {
            return [key, proj[key][0], proj[key][1], values[key]]
        });

        table = [["Cell", "X", "Y", plotted_item]].concat(table);

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
        if (main_vis === 'tree'){
            proj_name = get_global_status('plotted_trajectory')
        } else {
            proj_name = get_global_status('plotted_projection')
        }

        a.download = plotted_item+"_"+proj_name+".zip";
        a.href = zip_uri;
        a.click();

    };

    image.setAttribute("src", imgsrc)

}

Right_Content.prototype.hover_cells = function(cell_ids)
{
    this.scatter.hover_cells(cell_ids);
}

Right_Content.prototype.getSelectedCells = function() { 
    
    return this.scatter.getSelected();

}	
