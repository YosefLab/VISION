var global_stack = [];

var global_status = {};
global_status.main_vis = "sigvp"; // Selected from 3 options on top
/* Options are:
    'sigvp'
    'pcannotator'
    'tree'
*/

// Determine the projected coordinates
global_status.plotted_projection = "";
global_status.plotted_pc = 1;

// Determine projected values
global_status.plotted_item = "";  // name of signature, meta or gene that is plotted
global_status.plotted_item_type = ""; // either 'signature', 'meta', or 'gene'


global_status.upper_range = "";
global_status.lower_range = "";
global_status.pc1 = "";
global_status.pc2 = "";
global_status.subset = [];
global_status.subset_criteria = "Rank";


var global_data = {};
global_data.sigIsMeta = {};

global_data.sig_projection_coordinates = {};
global_data.pca_projection_coordinates = {};
global_data.tree_projection_coordinates = {};
global_data.plotted_values = {}; // Holds gene expression, signature scores/ranks, etc...
global_data.sig_info = {};  // Holds the information for the last plotted signature

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

    if(('plotted_projection' in update && get_global_status('main_vis') === 'tree') ||
       ('main_vis' in update && get_global_status('main_vis') === 'tree')
    ){
        var proj_key = get_global_status('plotted_projection');
        var proj_promise = api.tree.coordinates(proj_key)
            .then(function(projection){
                global_data.tree_projection_coordinates = projection;
            });

        all_promises.push(proj_promise);
        right_content_promises.push(proj_promise);
        upper_left_content_promises.push(proj_promise);
    }

    if(('plotted_projection' in update && get_global_status('main_vis') === 'sigvp') ||
       ('main_vis' in update && get_global_status('main_vis') === 'sigvp')
    ){
        var proj_key = get_global_status('plotted_projection');
        var proj_promise = api.projection.coordinates(proj_key)
            .then(function(projection){
                global_data.sig_projection_coordinates = projection;
            });

        all_promises.push(proj_promise);
        right_content_promises.push(proj_promise);
        upper_left_content_promises.push(proj_promise);
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
        } else if (type === 'gene') {
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

    $.when.apply($, right_content_promises).then(
        function() {
            right_content.update(update);
        }
    );

    $.when.apply($, lower_left_content_promises).then(
        function() {
            lower_left_content.update(update);
        }
    );

    $.when.apply($, upper_left_content_promises).then(
        function() {
            upper_left_content.update(update);
        }
    );


}

$(window).resize(function()
{
    right_content.resize()

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

    //Link the scatter/heatmap
    //global_scatter.hovered_links.push(global_heatmap);
    //global_heatmap.hovered_links.push(global_scatter);
    //global_tree.hovered_links.push(global_tree);

    //Render
    //drawHeat();


});

window.onload = function()
{


    //global_heatmap = new HeatMap("#heatmap-div");

    //Link the scatter/heatmap
    //global_scatter.hovered_links.push(global_heatmap);
    //global_heatmap.hovered_links.push(global_scatter);
    //global_tree.hovered_links.push(global_scatter);
    

    lower_left_content = new Lower_Left_Content()
    upper_left_content = new Upper_Left_Content()
    right_content = new Right_Content()


    var lower_left_promise = lower_left_content.init();
    var upper_left_promise = upper_left_content.init();
    var right_promise = right_content.init();

    // Get the 'isMeta' vector for signatures
    var sigIsMetaPromise = api.signature.listMeta()
        .then(function(sigIsMeta) {
            global_data.sigIsMeta = sigIsMeta;
        });

    // When it's all done, run this
    $.when(upper_left_promise, right_promise,
        sigIsMetaPromise, lower_left_promise
    )
        .then(function(){
            upper_left_content.select_default();
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


    api.sessionInfo().then(info => {
        if(info.name.length > 0){
            $('#SampleNameSpan').text(' - ' + info.name)
        }

        if(info.has_tree){
            $('#nav-bar')
                .find(".nav-link[data-main-vis='tree']")
                .removeClass('disabled')
        }
    });

    window.addEventListener('hover-cells', function(e) {
        var list_of_cell_ids = e.detail
        right_content.hover_cells(list_of_cell_ids)
        upper_left_content.hover_cells(list_of_cell_ids)
        lower_left_content.hover_cells(list_of_cell_ids)
    });

    /*
    $("#subset-criteria").change(function() {
        global_status.subset_criteria = $(this).val();
    });


    // Set Listeners for Cell Subset Analysis
    var upperRange = $("#upper-input");
    var lowerRange = $("#lower-input");

    upperRange.on("input", function() {
        global_status.upper_range = this.value; 
    });

    lowerRange.on("input", function() {
        global_status.lower_range = this.value;
    });


    // Set Listeners for PC Analysis
    var pc1 = $("#pc1-input");
    var pc2 = $("#pc2-input");

    pc1.on("input", function() {
        global_status.pc1 = this.value;
    });

    pc2.on("input", function() {
        global_status.pc2 = this.value;
    });


    //Define color option (for scatter) change function
    $('input[name=scatterColorButtons]').change(function(){
        var val = $('input[name=scatterColorButtons]:checked').val();
        global_status.scatterColorOption = val;
        right_content.update()();
    });

    $("#reload_heatmap").on("click", function() {
        console.log('here');
        drawHeat(); 
    });

    // Create listeners for main visualization modes
    $("#proj_tab").on("click", function() {
        global_status.main_vis = "sigvp";
        createTableFromData();
        right_content.update()();
    });

    $("#pc_tab").on("click", function() {

        global_status.main_vis = "pcannotator"
        createTableFromData();
        right_content.update()();


    });

    $("#tree_tab").on("click", function() {
        global_status.main_vis = "tree";
        createTableFromData();
        right_content.update()();
    });

    var criteriaList = ["Rank", "Value"];
    for (var i = 0; i < criteriaList.length; i++) {
        var criteria = criteriaList[i];
        var option = $(document.createElement("option"));
        option.text(criteria).val(criteria);
        $("#subset-criteria").append(option);
    }

    */

};

function goBack(d) {

    if (global_stack.length == 1) {
        return;
    }

    curr_status = global_stack.pop();
    global_status = global_stack.pop();
    createTableFromData();
    right_content.update()();
    drawHeat();
}

