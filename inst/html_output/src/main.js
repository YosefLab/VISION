var global_status = {};
global_status.main_vis = "cells"; // Selected from 4 options on top
/* Options are:
    'cells'
    'genes'
    'tree' (this is being phased out)
*/

// Determine the projected coordinates
global_status.plotted_projectionX = [null, null];  // [Name of projection dataframe, Name of column]
global_status.plotted_projectionY = [null, null];
global_status.plotted_trajectory = "";

// Indicate whether or not we have pooled data
global_status.pooled = false;

// Is trajectory data available?
global_status.has_tree = false;
global_status.has_sigs = true;
global_status.has_mods = true;
global_status.has_lca = false;

// Determine projected values
global_status.plotted_item = "";  // name of signature, meta or gene that is plotted
global_status.plotted_item_type = ""; // either 'signature', 'protein', 'meta', 'gene', or 'signature-gene'

global_status.cluster_var = ""; // which cluster variable are we using

global_status.selected_cell = ""; // which cell(s) is/are currently selected
global_status.selection_type = "none"; // either 'cell', or 'cells', or 'pool', or 'pools', or 'none'
global_status.selection_name = "Selection"; // For plot legends

// handling for scatterplot of modulescors vs signature scores
global_status.enrichment = false;
global_status.enrichment_module= "";


var global_data = {};

global_data.sig_projection_coordinatesX = {};
global_data.sig_projection_coordinatesY = {};
global_data.tree_projection_coordinates = {};
global_data.plotted_values = {}; // Holds gene expression, signature scores/ranks, etc...
global_data.sig_info = {};  // Holds the information for the last plotted signature
global_data.meta_sigs = []; // Holds list of sig names that are meta-data
global_data.extra_plotted_values = {}; // Holds extra info like module hotspot scores
global_data.mod_gene_list = {}; // Holds extra info like module hotspot scores

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
        if (global_status[key] === value && !key.endsWith('_changed'))
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

    // Updates plotted item values
    if('plotted_item' in update || 'plotted_item_type' in update) {

        var type = get_global_status('plotted_item_type');
        var sig_key = get_global_status('plotted_item');

        var val_promise;

        if(type === 'signature'){
            val_promise = api.signature.scores(sig_key)
        } else if (type === 'protein') {
            val_promise = api.protein.values(sig_key)
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
        upper_left_content_promises.push(val_promise_after);
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
    
    if('enrichment' in update || 'enrichment_module' in update) {
        var enriched = get_global_status("enrichment")
        if(enriched) {
            var mod_sig_promise = api.modules.hotspot_score(get_global_status("enrichment_module"))
            .then(new_info => {
                global_data.extra_plotted_values = new_info;
            })
        
            all_promises.push(mod_sig_promise);
            lower_left_content_promises.push(mod_sig_promise);
            
            
            var mod_gene_list_promise = api.modules.mod_gene_list(get_global_status("enrichment_module"))
            .then(new_info => {
                global_data.mod_gene_list = new_info;
                console.log(new_info);
            })
        
            all_promises.push(mod_gene_list_promise);
            lower_left_content_promises.push(mod_gene_list_promise);
            
            
        }
        
    } else if ('plotted_item' in update || 'plotted_item_type' in update) {
        global_status["enrichment"] = false;
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


    var right_promise = right_content.init();

    var sessionInfoPromise = api.sessionInfo().then(info => {
        if(info.name.length > 0){
            $('#SampleNameSpan').text(' - ' + info.name)
            document.title = 'VISION - ' + info.name
        }

        if(info.has_tree){
            $('#nav-bar')
                .find(".nav-link[data-main-vis='tree']")
                .removeClass('disabled')
        }

        if(info.has_lca){
            $('#nav-bar')
                .find(".nav-link[data-main-vis='genes']")
                .removeClass('disabled')
        }

        global_data.meta_sigs = info.meta_sigs
        global_status.pooled = info.pooled
        global_status.ncells = info.ncells
        global_status.has_tree = info.has_tree
        global_status.has_lca = info.has_lca
        global_status.has_sigs = info.has_sigs
        global_status.has_mods = info.has_mods
        global_status.has_proteins = info.has_proteins
        global_status.has_dendrogram = info.has_dendrogram
        global_status.dendrogram = info.dendrogram

        global_data.clusters_unique = info.clusters
    })

    var allCellClustersPromise = api.clusters.meta_levels()
        .then(function(data) {
            global_data.meta_levels = data;
        })


    var combinedPromise = $.when(sessionInfoPromise, allCellClustersPromise)
        .then(function(){
            return upper_left_content.init();
        })

    var combinedPromise2 = $.when(sessionInfoPromise, allCellClustersPromise)
        .then(function(){
            return lower_left_content.init();
        })

    // When it's all done, run this
    $.when(right_promise, combinedPromise, combinedPromise2, allCellClustersPromise)
        .then(function(){
            var has_tree = get_global_status('has_tree')
            var update0 = {'main_vis': (has_tree ? 'tree': 'cells')}
            var update1 = upper_left_content.select_default_sig();
            var update2 = right_content.select_default_proj();
            var update = Object.assign({}, update0, update1, update2); // Merge
            set_global_status(update)

            if (has_tree){ // need to change styling for navbar button to select Traj
                $('#nav-bar').find('.nav-link').removeClass('active')
                $('a[data-main-vis="tree"]').addClass('active')
            }
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

    // Pre-cache loading spinner
    var img = $('<img />', {
        src: 'css/loading.svg',
        alt: 'loading-spinner'
    });

};
