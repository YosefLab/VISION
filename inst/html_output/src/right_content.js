function Right_Content()
{
}

Right_Content.prototype.init = function()
{
    var self = this;
    self.dom_node = $("#right-content");

    self.scatterColorOptions = $(self.dom_node).find("input[name='scatterColorButtons']")
    self.scatterLayoutOptions = $(self.dom_node).find("input[name='scatterLayoutButtons']")

    self.selectedLayout = 'full'

    // Holds cached plot data -- needed when redrawing all plots
    self.layoutPlotData = {
        'full': [{}],
        'splitH': [{}, {}],
        'splitV': [{}, {}],
        'split4': [{}, {}, {}, {}],
    }

    self.layoutMeta = {
        'full': {selectedPlot: 0, initialized: true, plotDivs: ['#scatter-div']},
        'splitH': {selectedPlot: 0, initialized: false,
            plotDivs: [
                "#scatter-splitH-div0", "#scatter-splitH-div1",
            ]},
        'splitV': {selectedPlot: 0, initialized: false,
            plotDivs: [
                "#scatter-splitV-div0", "#scatter-splitV-div1",
            ]},
        'split4': {selectedPlot: 0, initialized: false,
            plotDivs: [
                "#scatter-split-div0", "#scatter-split-div1",
                "#scatter-split-div2", "#scatter-split-div3",
            ]},
    }

    self.layoutPlotData['full'][0]['scatter'] = new ColorScatter("#scatter-div");

    $(self.scatterColorOptions).on('change', function(){
        self.update({ 'colorScatterOption': '' }) // No need to send value
    })

    $(self.scatterLayoutOptions).on('change', function(eventObj){
        var selectedValue = eventObj.currentTarget.value

        $(self.dom_node).find('.scatter-layout-container').hide()

        if(selectedValue === 'full'){
            $(self.dom_node).find('#scatter-div').show()
        } else if(selectedValue === 'splitH') {
            $(self.dom_node).find('#scatter-splitH-div').show()
        } else if(selectedValue === 'splitV') {
            $(self.dom_node).find('#scatter-splitV-div').show()
        } else if(selectedValue === 'split4') {
            $(self.dom_node).find('#scatter-split-div').show()
        }

        self.changeLayout(selectedValue)
    })

    var _scatter_relayout = function(e) {
        // Used to synchronize zoom events between visible scatter plots
        // First need to remove the listener to prevent infinite loop

        self.dom_node.get(0).removeEventListener('scatter_relayout', _scatter_relayout)

        var promises = []

        var visiblePlotData = self.getVisiblePlotData()
        for (var i = 0; i < visiblePlotData.length; i++){
            if (e.detail.origin !== visiblePlotData[i]['scatter']){
                promises.push(
                    visiblePlotData[i]['scatter'].relayout(e.detail.newLayout)
                )
            }
        }

        Promise.all(promises).then(function(){
            self.dom_node.get(0).addEventListener('scatter_relayout', _scatter_relayout)
        })
    }

    self.dom_node.get(0).addEventListener('scatter_relayout', _scatter_relayout)

    // Allow for clicking to specify selected plot div
    $(self.dom_node).find(".scatter-split-plot-div")
        .on('click', function(e){
            var targetId = e.currentTarget.id
            var layoutMeta = self.layoutMeta[self.selectedLayout]
            layoutMeta['selectedPlot'] = parseInt(targetId[targetId.length-1])

            // Need to clear current selection, but only inside the current scatter-layout-container
            $(e.currentTarget).parents('.scatter-layout-container')
                .find('.scatter-split-plot-div')
                .removeClass('active')

            $(e.currentTarget).addClass('active')

            // Fire the update so the lower-left area changes too
            //
            var plottedData = self.getSelectedPlotData()

            var update = {}
            update['plotted_item_type'] = plottedData['item_type'];
            update['plotted_item'] = plottedData['item_key'];

            set_global_status(update);
        })

    // Set up events for saving/loading selections
    var saveSelectionModal = $('#saveSelectionModal')
    saveSelectionModal.on('shown.bs.modal', function() {
        // Update selection name if it has a name already
        var selectionName = get_global_status('selection_name')
        var selectionInput = $(this).find('#saveSelectionName')
        if (selectionName !== "Selection"){
            selectionInput.val(selectionName)
        } else {
            selectionInput.val("")
        }
        selectionInput.focus()
    }).on('keyup', function(event){
        if(event.key != "Enter") return;
        $(this).find('.confirm-button').click()
        event.preventDefault();
    })

    saveSelectionModal.find(".confirm-button").on("click", function() {
        // Save selection and close modal
        var selectionName = saveSelectionModal.find('#saveSelectionName').val()

        if (selectionName === "") return;

        var selectedCells = get_global_status('selected_cell')

        api.cells.saveSelection(selectionName, selectedCells)

        set_global_status({
            'selection_name': selectionName
        })
        saveSelectionModal.modal('hide')

    })

    var loadSelectionModal = $('#loadSelectionModal')
    var lSM_loading = createLoadingFunction(loadSelectionModal.find('.modal-body'))
    loadSelectionModal.on('show.bs.modal', function() {
        // Load selections and display

        var modalList = loadSelectionModal.find('.selection-list')
        var emptyList = loadSelectionModal.find('.empty-list')
        lSM_loading(true)

        api.cells.listSelections().then(data => {
            modalList.empty()

            var confirmButton = loadSelectionModal.find('.confirm-button')
            if(data.length === 0){

                modalList.hide()
                emptyList.show()
                confirmButton.attr('disabled', true)
                return
            }
            emptyList.hide()
            modalList.show()

            confirmButton.attr('disabled', false)

            var itemClickFun = function(){
                modalList.find('.loadable-selection-option').removeClass('selected')
                $(this).addClass('selected')
            }

            _.each(data, (name, i) => {
                var e = document.createElement('div')
                e.innerText = name
                e.classList.add('loadable-selection-option')
                if (i == 0){ // Select first element
                    e.classList.add('selected')
                }
                e.addEventListener('click', itemClickFun)
                modalList.append(e)
            })

        }).always(() => {
            lSM_loading(false)
        })

    }).on('keyup', function(event){
        if(event.key != "Enter") return;
        $(this).find('.confirm-button').click()
        event.preventDefault();
    })

    loadSelectionModal.find(".confirm-button").on("click", function() {
        // Load selection and close modal
        self.setLoadingStatus(true)

        var selectionToLoad = loadSelectionModal.find('.loadable-selection-option.selected').text()

        api.cells.getSelection(selectionToLoad).then(cells => {
            var update = {
                selected_cell: cells,
                selection_type: get_global_status('pooled')? 'pools' : 'cells',
                selection_name: selectionToLoad,
            }
            set_global_status(update)
        })


        loadSelectionModal.modal('hide')
    })

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

    var plotData = self.getSelectedPlotData()

    // Update the selection status in ALL plots
    if('selected_cell' in updates){
        var visiblePlotData = this.getVisiblePlotData()
        for (var i = 0; i < visiblePlotData.length; i++){

            // Need to check in case it's not initialized
            if ('scatter' in visiblePlotData[i]) {
                visiblePlotData[i]['scatter'].updateSelection()
            }
        }

        var selection_type = get_global_status('selection_type')
        var selection_button = $(this.dom_node).find('#save-selection-button')
        if (selection_type === 'cells' || selection_type === 'pools'){
            selection_button.attr('disabled', false)
        } else {
            selection_button.attr('disabled', true)
        }

    }

    var needsUpdate = ('main_vis' in updates) ||
        ('plotted_item' in updates) ||
        ('plotted_item_type' in updates) ||
        ('plotted_projection' in updates) ||
        ('plotted_trajectory' in updates) ||
        ('colorScatterOption' in updates);

    if (!needsUpdate) return $.Deferred().resolve().promise()

    var main_vis = get_global_status('main_vis');

    var autoZoom = false

    if('plotted_projection' in updates ||
       'plotted_trajectory' in updates ||
       'main_vis' in updates) {

        autoZoom = true
        var projection;
        if (get_global_status('main_vis') === 'tree'){
            projection = get_global_data('tree_projection_coordinates')
        } else {
            projection = get_global_data('sig_projection_coordinates')
        }

        // If the projection is changing, then we need to change all the plots
        var allPlotData = this.getAllPlotData()
        for (var i = 0; i < allPlotData.length; i++){
            allPlotData[i]['projection'] = projection
            allPlotData[i]['needsUpdate'] = true
        }
    }

    if ('colorScatterOption' in updates) {

        // Mark all plots for update
        var allPlotData = this.getAllPlotData()
        for (var i = 0; i < allPlotData.length; i++){
            allPlotData[i]['needsUpdate'] = true
        }
    }

    var item_key = get_global_status('plotted_item');
    var item_type = get_global_status('plotted_item_type');
    var values = get_global_data('plotted_values')

    // If this, then we just need to change the selected visualization
    if('plotted_item' in updates ||
       'plotted_item_type' in updates ||
       'plotted_values' in updates) {

        plotData['item_key'] = item_key
        plotData['item_type'] = item_type
        plotData['values'] = values
        plotData['needsUpdate'] = true

    }

    // Update all visible plots if they need it
    var visiblePlotData = this.getVisiblePlotData()
    for (i = 0; i < visiblePlotData.length; i++){
        if (visiblePlotData[i]['needsUpdate']){
            this.draw_scatter(visiblePlotData[i], autoZoom)
        }
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

    // Show either the Projection or Tree lists
    if('main_vis' in updates){
        $('#plot-subtitle-latent').hide()
        $('#plot-subtitle-trajectory').hide()
        if(main_vis === 'tree'){
            $('#plot-subtitle-trajectory').show()
        } else {
            $('#plot-subtitle-latent').show()
        }
    }

    // Need to return a promise
    return $.Deferred().resolve().promise()

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

/*
 * scatter: Instance of ColorScatter - which plot to update
 * item_key: string - name of plotted item
 * item_type: string - one of 'signature', 'meta', 'gene', or 'signature-gene'
 * values: sample name (str) -> int/str - values to plot
 * projection: sample name (str) -> [x, y] - coordinates for each sample
 * autoZoom: bool - whether or not to re-zoom the plot
 */
Right_Content.prototype.draw_sigvp = function(scatter, item_key, item_type, values, projection, autoZoom) {

    var self = this;

    var sample_value = _.values(values)[0]
    var isFactor = (typeof(sample_value) === 'string') &&
                   (sample_value !== "NA")

    var full_color_range, diverging_colormap
    if(item_type === "gene" || item_type === 'signature-gene'){
        full_color_range = true
        diverging_colormap = false
    } else if(item_type === "meta"){
        full_color_range = false
        diverging_colormap = true
    } else {
        full_color_range = false
        diverging_colormap = true
    }

    if(self.getScatterColorOption() == "rank" && item_type === 'signature') {

        values = self.rank_values(values)
    }


    var points = [];
    var sample_labels = Object.keys(values).sort()
    var selected_cells = get_global_status('selected_cell')
    var selected_cells_map = _.keyBy(selected_cells, x => x)
    if(selected_cells.length == 1){ // Just a single cell
        selected_cells_map = {} // Don't style anything
    }

    _.each(sample_labels, (sample_label) => {
        var x = projection[sample_label][0]
        var y = projection[sample_label][1]
        var sig_score = values[sample_label]
        var selected = sample_label in selected_cells_map

        points.push({
            x: x, y: y,
            value: sig_score, label: sample_label,
            selected: selected,
        });
    })

    scatter.setData({
        points: points,
        isFactor: isFactor,
        full_color_range: full_color_range,
        diverging_colormap: diverging_colormap,
        title: item_key,
        autozoom: autoZoom,
    });

}


/*
 * scatter: Instance of ColorScatter - which plot to update
 * item_key: string - name of plotted item
 * item_type: string - one of 'signature', 'meta', 'gene', or 'signature-gene'
 * values: sample name (str) -> int/str - values to plot
 * projection: list
 *             [0]: sample name (str) -> [x, y] - coordinates for each sample
 *             [1]: Milestone coordiantes (list), [vData, adjMat]
 * autoZoom: bool - whether or not to re-zoom the plot
 */
Right_Content.prototype.draw_tree = function(scatter, item_key, item_type, values, projection, autoZoom) {

    var self = this;

    var isFactor = (typeof(_.values(values)[0]) === 'string') &&
                   (_.values(values)[0] !== "NA")

    var full_color_range, diverging_colormap
    if(item_type === "gene" || item_type === "signature-gene"){
        full_color_range = true
        diverging_colormap = false
    } else if(item_type === "meta"){
        full_color_range = false
        diverging_colormap = true
    } else {
        full_color_range = false
        diverging_colormap = true
    }

    if(self.getScatterColorOption() == "rank"
        && !isFactor && item_type !== "gene"
        && item_type !== "meta") {

        values = self.rank_values(values)
    }

    var milestoneCoordinates = projection[1]
    var projection = projection[0]

    var treep = milestoneCoordinates[0]
    var treel = milestoneCoordinates[1]

    // Massage treep for easier D3 binding

    tree_points = []

    var points = [];
    var sample_labels = Object.keys(values).sort()
    var selected_cells = get_global_status('selected_cell')
    var selected_cells_map = _.keyBy(selected_cells, x => x)
    if(selected_cells.length == 1){ // Just a single cell
        selected_cells_map = {} // Don't style anything
    }

    _.each(sample_labels, (sample_label) => {
        var x = projection[sample_label][0]
        var y = projection[sample_label][1]
        var sig_score = values[sample_label]
        var selected = sample_label in selected_cells_map

        points.push({
            x: x, y: y,
            value: sig_score, label: sample_label,
            selected: selected,
        });
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

    scatter.setData({
        points: points,
        isFactor: isFactor,
        full_color_range: full_color_range,
        diverging_colormap: diverging_colormap,
        tree_points: tree_points,
        tree_adj: tree_adj,
        title: item_key,
        autozoom: autoZoom,
    });
}


// Called when the window is resized
Right_Content.prototype.resize = function() {

    var visiblePlots = this.getVisiblePlotData()

    for (var i = 0; i < visiblePlots.length; i++) {
        visiblePlots[i]['scatter'].resize()
    }
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

Right_Content.prototype.hover_cells = function(cell_ids)
{
    var visiblePlots = this.getVisiblePlotData()

    for (var i = 0; i < visiblePlots.length; i++) {
        visiblePlots[i]['scatter'].hover_cells(cell_ids)
    }
}

Right_Content.prototype.changeLayout = function(newLayout)
{
    // Old plot data
    var refPlotData = this.getSelectedPlotData()
    var refVisiblePlotData = this.getVisiblePlotData()

    // Change layout here
    this.selectedLayout = newLayout
    var newLayoutMeta = this.layoutMeta[this.selectedLayout]

    var splitInitialized = newLayoutMeta['initialized']

    var plotData, visiblePlotData
    if(splitInitialized){

        plotData = this.getSelectedPlotData()

        for (var prop in refPlotData){
            if(prop === 'scatter') {continue; }
            if(prop === 'needsUpdate') {continue; }
            plotData[prop] = refPlotData[prop]
        }

    } else {
        // New plot data
        visiblePlotData = this.getVisiblePlotData()
        var plotDivs = newLayoutMeta['plotDivs']

        // Loop through visible splits assigning their content from exists splits
        for (var i = 0; i < visiblePlotData.length; i++){

            plotData = visiblePlotData[i]
            refPlotData = refVisiblePlotData[i % refVisiblePlotData.length]

            plotData['scatter'] = new ColorScatter(plotDivs[i])
            plotData['needsUpdate'] = true

            for (prop in refPlotData){
                if (prop === 'scatter') {continue; }
                if(prop === 'needsUpdate') {continue; }
                plotData[prop] = refPlotData[prop]
            }

        }

        newLayoutMeta['initialized'] = true
    }

    visiblePlotData = this.getVisiblePlotData()

    // Need to purge all non-visible plots to reduce WebGL context count
    // Then need to re-create all visible plots
    var allPlotData = this.getAllPlotData()
    for (var i = 0; i < allPlotData.length; i++) {
        if ('scatter' in allPlotData[i]) {
            allPlotData[i]['scatter'].clear()
        }
    }

    var autoZoom = false
    for (var i = 0; i < visiblePlotData.length; i++) {
        $.extend(visiblePlotData[i]['scatter'].currentZoom, refPlotData['scatter'].currentZoom)
        this.draw_scatter(visiblePlotData[i], autoZoom)
    }

}

Right_Content.prototype.draw_scatter = function(plotData, autoZoom)
{
    var main_vis = get_global_status('main_vis');

    if(main_vis === 'clusters' || main_vis === "pcannotator"){
        this.draw_sigvp(plotData['scatter'], plotData['item_key'], plotData['item_type'], plotData['values'], plotData['projection'], autoZoom);

    } else if (main_vis === "tree") {
        this.draw_tree(plotData['scatter'], plotData['item_key'], plotData['item_type'], plotData['values'], plotData['projection'], autoZoom);

    } else {
        throw "Bad main_vis value!";
    }

    plotData['needsUpdate'] = false

}

Right_Content.prototype.getSelectedPlotData = function()
{
    var layoutMeta = this.layoutMeta[this.selectedLayout]
    var visiblePlotData = this.layoutPlotData[this.selectedLayout]

    return visiblePlotData[layoutMeta.selectedPlot]
}

Right_Content.prototype.getVisiblePlotData = function()
{
    var visiblePlotData = this.layoutPlotData[this.selectedLayout]
    return visiblePlotData
}

Right_Content.prototype.getAllPlotData = function()
{
    return _.flatMap(this.layoutPlotData, x => _.values(x))
}

Right_Content.prototype.getAllLayoutMeta = function()
{
    return _.flatMap(this.layoutMeta, x => _.values(x))
}
