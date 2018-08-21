var global_status = {};
global_status.main_vis = "clusters"; // Selected from 4 options on top
/* Options are:
    'clusters'
    'pcannotator'
    'tree'
*/

// Determine the projected coordinates
global_status.plotted_projection = "";
global_status.plotted_trajectory = "";
global_status.plotted_pc = 1;

// Indicate whether or not we have pooled data
global_status.pooled = false;

// Determine projected values
global_status.plotted_item = "";  // name of signature, meta or gene that is plotted
global_status.plotted_item_type = ""; // either 'signature', 'meta', 'gene', or 'signature-gene'

global_status.cluster_var = ""; // which cluster variable are we using

global_status.selected_cell = ""; // which cell(s) is/are currently selected
global_status.selection_type = "none"; // either 'cell', or 'cells', or 'pool', or 'pools', or 'none'
global_status.selection_name = "Selection"; // For plot legends

var global_data = {};

global_data.sig_projection_coordinates = {};
global_data.pca_projection_coordinates = {};
global_data.tree_projection_coordinates = {};
global_data.plotted_values = {}; // Holds gene expression, signature scores/ranks, etc...
global_data.sig_info = {};  // Holds the information for the last plotted signature
global_data.clusters = {};  // Maps cell ID to cluster ID
global_data.cluster_variables = [];
global_data.meta_sigs = []; // Holds list of sig names that are meta-data

var lower_left_content;
var upper_left_content;
var right_content;

function get_global_data(key){
    return global_data[key];
}

function get_global_status(key){
    return global_status[key];
}

function set_global_status(update){

    // Delete all updates that are not changes
    var toDelete = []
    _.forIn(update, function(value, key) {
        if (global_status[key] === value)
        {
            toDelete.push(key);
        }
    });
    _.forEach(toDelete, function(key) {
        delete update[key]
    });

    // Now all updates represent changes

    _.forIn(update, function(value, key) {
        global_status[key] = value;
    });

    // Now all global_status are current, and the 'update' dict
    // holds keys for entries that just changed
    //
    var all_promises = [];
    var right_content_promises = [];
    var lower_left_content_promises = [];
    var upper_left_content_promises = [];

    right_content.setLoadingStatus(true);
    lower_left_content.setLoadingStatus(true);
    upper_left_content.setLoadingStatus(true);

    // Just get all the PC's using one call
    if('main_vis' in update && get_global_status('main_vis') === 'pcannotator' &&
        _.isEmpty(global_data.pca_projection_coordinates))
    {
        var pc_promise = api.pc.coordinates()
            .then(function(projection){
                global_data.pca_projection_coordinates = projection;
            });

        all_promises.push(pc_promise);
        right_content_promises.push(pc_promise);
        lower_left_content_promises.push(pc_promise);
    }

    if(('plotted_trajectory' in update && get_global_status('main_vis') === 'tree') ||
       ('main_vis' in update && get_global_status('main_vis') === 'tree')
    ){
        var proj_key = get_global_status('plotted_trajectory');
        var proj_promise = api.tree.coordinates(proj_key)
            .then(function(projection){
                global_data.tree_projection_coordinates = projection;
            });

        all_promises.push(proj_promise);
        right_content_promises.push(proj_promise);
        upper_left_content_promises.push(proj_promise);
    }

    if('plotted_projection' in update || 'main_vis' in update){
        if( get_global_status('main_vis') === 'clusters' ||
            get_global_status('main_vis') === 'pcannotator'){

            var proj_key = get_global_status('plotted_projection');
            var proj_promise = api.projections.coordinates(proj_key)
                .then(function(projection){
                    global_data.sig_projection_coordinates = projection;
                });

            all_promises.push(proj_promise);
            right_content_promises.push(proj_promise);
            upper_left_content_promises.push(proj_promise);

        }
    }


    // Updates plotted item values
    if('plotted_item' in update || 'plotted_item_type' in update) {

        var type = get_global_status('plotted_item_type');
        var sig_key = get_global_status('plotted_item');

        var val_promise;

        if(type === 'signature'){
            val_promise = api.signature.scores(sig_key)
        } else if (type === 'meta') {
            val_promise = api.signature.meta(sig_key)
        } else if (type === 'gene' || type === 'signature-gene') {
            val_promise = api.expression.gene(sig_key)
        } else {
            throw 'Bad "plotted_item_type"!';
        }

        var val_promise_after = val_promise.then(function(values)
        {
            global_data.plotted_values = values;
        });

        all_promises.push(val_promise_after);
        right_content_promises.push(val_promise_after);
        lower_left_content_promises.push(val_promise_after);

    }

    // Updates signature info if we're plotting a new signature
    if('plotted_item' in update || 'plotted_item_type' in update) {

        var type = get_global_status('plotted_item_type');
        var new_sig = get_global_status('plotted_item');
        var sig_info = get_global_data('sig_info');

        if((type === 'signature') &&
            (_.isEmpty(sig_info) || sig_info.name !== new_sig)
        )
        {
            var sig_info_promise = api.signature.info(new_sig)
                .then(new_info => {
                    global_data.sig_info = new_info;
                })

            all_promises.push(sig_info_promise);
            lower_left_content_promises.push(sig_info_promise);

        }

    }

    if('cluster_var' in update) {
        var cluster_var = get_global_status('cluster_var')
        var cell_clusters_promise = api.clusters.cells(cluster_var)
            .then( data => {
                global_data.clusters = data
            })
        all_promises.push(cell_clusters_promise);
        right_content_promises.push(cell_clusters_promise);
        lower_left_content_promises.push(cell_clusters_promise);
        upper_left_content_promises.push(cell_clusters_promise);

    }

    $.when.apply($, right_content_promises).then(
        function() {
            right_content.update(update).then(() => {
                right_content.setLoadingStatus(false);
            })
        }
    );

    $.when.apply($, lower_left_content_promises).then(
        function() {
            lower_left_content.update(update).then(() => {
                lower_left_content.setLoadingStatus(false);
            })
        }
    );

    $.when.apply($, upper_left_content_promises).then(
        function() {
            upper_left_content.update(update).then(() => {
                upper_left_content.setLoadingStatus(false);
            })
        }
    );


}

window.onresize = function()
{
    right_content.resize()
    lower_left_content.resize()

    /*
    if($('#heatmap-div').is(":visible"))
    {
        $('#heatmap-div').find('svg').remove();
        global_heatmap = new HeatMap('#heatmap-div');
    }

    if ($('#tree_div').is(":visible"))
    {
        $('#tree_div').find('svg').remove();
    }
    */

    //Render
    //drawHeat();

};

window.onload = function()
{

    lower_left_content = new Lower_Left_Content()
    upper_left_content = new Upper_Left_Content()
    right_content = new Right_Content()


    var lower_left_promise = lower_left_content.init();
    var right_promise = right_content.init();

    // Get the cluster assignments for cells
    var cellClustersPromise = api.clusters.list()
        .then(function(data) {
            global_data.cluster_variables = data;
            global_status.cluster_var = global_data.cluster_variables[0];

            return upper_left_content.init();
        }).then(function() {

            var cluster_var = get_global_status('cluster_var')
            var cell_clusters_promise = api.clusters.cells(cluster_var)
                .then( data => {
                    global_data.clusters = data
                })
            return cell_clusters_promise
        })

    var sessionInfoPromise = api.sessionInfo().then(info => {
        if(info.name.length > 0){
            $('#SampleNameSpan').text(' - ' + info.name)
        }

        if(info.has_tree){
            $('#nav-bar')
                .find(".nav-link[data-main-vis='tree']")
                .removeClass('disabled')
        }

        global_data.meta_sigs = info.meta_sigs
        global_status.pooled = info.pooled
        global_status.ncells = info.ncells
    });

    // When it's all done, run this
    $.when(right_promise, lower_left_promise,
        cellClustersPromise, sessionInfoPromise)
        .then(function(){
            var update0 = {'main_vis': 'clusters'}
            var update1 = upper_left_content.select_default_sig();
            var update2 = right_content.select_default_proj();
            var update = Object.assign({}, update0, update1, update2); // Merge
            set_global_status(update)
        });

    // Enable the nav-bar functionality
    $('#nav-bar').find('.nav-link').click(function(){
        if( $(this).hasClass('disabled') ) { return; }
        $('#nav-bar').find('.nav-link').removeClass('active')
        $(this).addClass('active')
        set_global_status({
            'main_vis': $(this).data('main-vis')
        })
    });

    window.addEventListener('hover-cells', function(e) {

        // Disable hover events if there are too many cells
        // This is necessary for performance
        var ncells = get_global_status('ncells')
        if( ncells < 10000){
            var list_of_cell_ids = e.detail
            right_content.hover_cells(list_of_cell_ids)
            upper_left_content.hover_cells(list_of_cell_ids)
            lower_left_content.hover_cells(list_of_cell_ids)
        }
    });

    window.addEventListener("select-cells", function(e) {

        var cells = e.detail.cells;
        var update = {}

        if (cells.length == 0){
            update['selection_type'] = 'none'
            update['selected_cell'] = cells
            set_global_status(update)
            return;
        }

        if(get_global_data('pooled')){
            if(cells.length == 1){
                update['selected_cell'] = cells
                update['selection_type'] = 'pool'
            } else {
                update['selected_cell'] = cells
                update['selection_type'] = 'pools'
            }
        } else {
            if(cells.length == 1){
                update['selected_cell'] = cells
                update['selection_type'] = 'cell'
            } else {
                update['selected_cell'] = cells
                update['selection_type'] = 'cells'
            }
        }

        var name = 'Selection'
        if('name' in e.detail){
            name = e.detail.name
        }
        update['selection_name'] = name;

        set_global_status(update)
    });

};
