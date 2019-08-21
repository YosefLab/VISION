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
    self.syncCheck = $(this.dom_node).find('#sync-check').get(0);
    self.treeView = false

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
            if (e.detail.origin !== visiblePlotData[i]['scatter'] &&
                e.detail.projKeyX[0] === visiblePlotData[i]['projection_keyX'][0] && // Only sync zoom if same plot
                e.detail.projKeyX[1] === visiblePlotData[i]['projection_keyX'][1] &&
                e.detail.projKeyY[0] === visiblePlotData[i]['projection_keyY'][0] &&
                e.detail.projKeyY[1] === visiblePlotData[i]['projection_keyY'][1]
            ){
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
            var plottedData = self.getSelectedPlotData()

            self.syncSelectedProj()

            // Fire the update so the lower-left area changes too
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

            self.proj_names_scatter = proj_names
            var projSelect = self.dom_node.find('#SelectProjScatter')
            var projSelectX = self.dom_node.find('#SelectProjScatterX')
            var projSelectY = self.dom_node.find('#SelectProjScatterY')
            projSelect.children().remove()

            _.each(proj_names, function (val, key) {
                projSelect.append(
                    $('<option>', {
                        value: key,
                        text: key
                    }));
            });

            var update_axes_dropdowns = function(){
                var selected_proj = $(projSelect).val()

                projSelectX.children().remove()
                projSelectY.children().remove()

                _.each(proj_names[selected_proj], function (val) {
                    projSelectX.append(
                        $('<option>', {
                            value: val,
                            text: val
                        }));
                    projSelectY.append(
                        $('<option>', {
                            value: val,
                            text: val
                        }));
                });

                $(projSelectX).val(proj_names[selected_proj][0])
                $(projSelectY).val(proj_names[selected_proj][1])
                projSelectX.trigger('chosen:updated')
                projSelectY.trigger('chosen:updated')
            }

            projSelectX.chosen({
                'width': '150px',
                'disable_search_threshold': 99,
            })
                .off('change')
                .on('change', function () {
                    set_global_status({
                        'plotted_projectionX': [$(projSelect).val(), $(projSelectX).val()],
                        'plotted_projectionY': [$(projSelect).val(), $(projSelectY).val()],
                    });
                })

            projSelectY.chosen({
                'width': '150px',
                'disable_search_threshold': 99,
            })
                .off('change')
                .on('change', function () {
                    set_global_status({
                        'plotted_projectionX': [$(projSelect).val(), $(projSelectX).val()],
                        'plotted_projectionY': [$(projSelect).val(), $(projSelectY).val()],
                    });
                })

            projSelect.chosen({
                'width': '150px',
                'disable_search_threshold': 99,
            })
                .off('change')
                .on('change', function () {
                    update_axes_dropdowns()
                    set_global_status({
                        'plotted_projectionX': [$(this).val(), $(projSelectX).val()],
                        'plotted_projectionY': [$(this).val(), $(projSelectY).val()],
                    });
                })
                .trigger('chosen:updated')

            update_axes_dropdowns()

        });

    var treeproj_promise = api.tree.list()
        .then(function(proj_names) {

            self.proj_names_trajectory = proj_names
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
        var visiblePlotData = self.getVisiblePlotData()
        for (var i = 0; i < visiblePlotData.length; i++){

            // Need to check in case it's not initialized
            if ('scatter' in visiblePlotData[i]) {
                visiblePlotData[i]['scatter'].updateSelection()
            }
        }

        var selection_type = get_global_status('selection_type')
        var selection_button = $(self.dom_node).find('#save-selection-button')
        if (selection_type === 'cells' || selection_type === 'pools'){
            selection_button.attr('disabled', false)
        } else {
            selection_button.attr('disabled', true)
        }

    }

    // Show either the Projection or Tree lists
    var treeViewChanged = false

    if('main_vis' in updates){

        var main_vis = get_global_status('main_vis');

        $('#plot-subtitle-latent').hide()
        $('#plot-subtitle-trajectory').hide()
        if(main_vis === 'tree'){
            $('#plot-subtitle-trajectory').show()
            treeViewChanged = self.treeView === false
            self.treeView = true
        } else {
            $('#plot-subtitle-latent').show()
            treeViewChanged = self.treeView === true
            self.treeView = false
        }
    }

    if ('colorScatterOption' in updates) {

        // Mark all plots for update
        var allPlotData = self.getAllPlotData()
        for (var i = 0; i < allPlotData.length; i++){
            allPlotData[i]['needsUpdate'] = true
        }
    }

    // If this, then we just need to change the selected visualization
    if('plotted_item' in updates ||
       'plotted_item_type' in updates ||
       'plotted_values' in updates) {

        var item_key = get_global_status('plotted_item');
        var item_type = get_global_status('plotted_item_type');
        var values = get_global_data('plotted_values')

        plotData['item_key'] = item_key
        plotData['item_type'] = item_type
        plotData['values'] = values
        plotData['needsUpdate'] = true

    }

    var autoZoom = false

    var update_promise;

    if('plotted_projectionX' in updates ||
       'plotted_projectionY' in updates ||
       'plotted_trajectory' in updates) {

        autoZoom = true
        var proj_promise;
        var projection;
        var projection_type;

        if (
            ('plotted_projectionX' in updates || 'plotted_projectionY' in updates) ||
            (treeViewChanged && !self.treeView)
        ){

            var proj_keyX = get_global_status('plotted_projectionX');
            var proj_keyY = get_global_status('plotted_projectionY');
            var proj_promiseX = api.projections.coordinates(proj_keyX[0], proj_keyX[1])
            var proj_promiseY = api.projections.coordinates(proj_keyY[0], proj_keyY[1])
            proj_promise = $.when(proj_promiseX, proj_promiseY)
                .then(function(projectionX, projectionY){
                    var keys = _.keys(projectionX)
                    var values = _.map(keys, key => [projectionX[key], projectionY[key]])
                    var projection = _.zipObject(keys, values)
                    return {
                        'projection': projection,
                        'projection_type': 'default',
                        'projection_keyX': proj_keyX,
                        'projection_keyY': proj_keyY,
                    }
                });
        } else if (
            'plotted_trajectory' in updates ||
            (treeViewChanged && self.treeView)
        ) {
            var proj_key = get_global_status('plotted_trajectory');
            proj_promise = api.tree.coordinates(proj_key)
                .then(function(projection){
                    return {
                        'projection': projection,
                        'projection_type': 'tree',
                        'projection_keyX': [proj_key, ''],
                        'projection_keyY': [proj_key, ''],
                    }
                });
        }

        update_promise = proj_promise.then(function(data) {
            var projection = data['projection']
            var projection_type = data['projection_type']
            var projection_keyX = data['projection_keyX']
            var projection_keyY = data['projection_keyY']

            // If the projection is changing, then we need to change all the plots
            var synchronized = $(self.syncCheck).is(':checked');
            if (synchronized || treeViewChanged) {
                var allPlotData = self.getAllPlotData()
                for (var i = 0; i < allPlotData.length; i++){
                    allPlotData[i]['projection'] = projection
                    allPlotData[i]['projection_type'] = projection_type
                    allPlotData[i]['projection_keyX'] = projection_keyX
                    allPlotData[i]['projection_keyY'] = projection_keyY
                    allPlotData[i]['needsUpdate'] = true
                }
            } else {
                var selectedPlotData = self.getSelectedPlotData()
                selectedPlotData['projection'] = projection
                selectedPlotData['projection_type'] = projection_type
                selectedPlotData['projection_keyX'] = projection_keyX
                selectedPlotData['projection_keyY'] = projection_keyY
                selectedPlotData['needsUpdate'] = true
            }
        })
    } else {
        update_promise = $.Deferred().resolve().promise()
    }

    return update_promise.then(function() {
        // Update all visible plots if they need it
        var visiblePlotData = self.getVisiblePlotData()
        for (i = 0; i < visiblePlotData.length; i++){
            if (visiblePlotData[i]['needsUpdate']){
                self.draw_scatter(visiblePlotData[i], autoZoom)
            }
        }
    })
}

Right_Content.prototype.select_default_proj = function()
{
    var proj_name = $(this.dom_node.find('#SelectProjScatter')).val()
    var proj_X = $(this.dom_node.find('#SelectProjScatterX')).val()
    var proj_Y = $(this.dom_node.find('#SelectProjScatterY')).val()
    var traj = $(this.dom_node.find('#SelectTrajectoryProjScatter')).val()

    var update = {}
    update['plotted_projectionX'] = [proj_name, proj_X]
    update['plotted_projectionY'] = [proj_name, proj_Y]
    update['plotted_trajectory'] = traj
    return update;
}

/*
 * This is used to make sure that changing the selected projection
 * updates the dropdown menus up top accordingly
 */
Right_Content.prototype.syncSelectedProj = function()
{
    var self = this
    var plottedData = self.getSelectedPlotData()
    var proj_keyX = plottedData['projection_keyX']
    var proj_keyY = plottedData['projection_keyY']

    if (proj_keyX[0] in self.proj_names_scatter) {
        var projSelect = self.dom_node.find('#SelectProjScatter')
        var projSelectX = self.dom_node.find('#SelectProjScatterX')
        var projSelectY = self.dom_node.find('#SelectProjScatterY')
        $(projSelect).val(proj_keyX[0])

        projSelectX.children().remove()
        projSelectY.children().remove()

        _.each(self.proj_names_scatter[proj_keyX[0]], function (val) {
            projSelectX.append(
                $('<option>', {
                    value: val,
                    text: val
                }));
            projSelectY.append(
                $('<option>', {
                    value: val,
                    text: val
                }));
        });

        $(projSelectX).val(proj_keyX[1])
        $(projSelectY).val(proj_keyY[1])
        projSelect.trigger('chosen:updated')
        projSelectX.trigger('chosen:updated')
        projSelectY.trigger('chosen:updated')
    } else {
        var projSelect = self.dom_node.find('#SelectTrajectoryProjScatter')
        $(projSelect).val(proj_keyX[0])
        projSelect.trigger('chosen:updated')
    }

}

/*
 * scatter: Instance of ColorScatter - which plot to update
 * item_key: string - name of plotted item
 * item_type: string - one of 'signature', 'meta', 'gene', or 'signature-gene'
 * values: sample name (str) -> int/str - values to plot
 * projection: sample name (str) -> [x, y] - coordinates for each sample
 * autoZoom: bool - whether or not to re-zoom the plot
 */
Right_Content.prototype.draw_sigvp = function(plotData, autoZoom) {

    var self = this;
    var scatter = plotData['scatter']
    var item_key = plotData['item_key']
    var item_type = plotData['item_type']
    var values = plotData['values']
    var projection = plotData['projection']
    var proj_keyX = plotData['projection_keyX']
    var proj_keyY = plotData['projection_keyY']

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
        proj_keyX: proj_keyX,
        proj_keyY: proj_keyY,
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
Right_Content.prototype.draw_tree = function(plotData, autoZoom) {

    var self = this;

    var scatter = plotData['scatter']
    var item_key = plotData['item_key']
    var item_type = plotData['item_type']
    var values = plotData['values']
    var projection = plotData['projection']
    var proj_keyX = plotData['projection_keyX']
    var proj_keyY = plotData['projection_keyY']

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
        proj_keyX: proj_keyX,
        proj_keyY: proj_keyY,
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

    if (plotData['projection_type'] === 'default') {
        this.draw_sigvp(plotData, autoZoom);

    } else if (plotData['projection_type'] === 'tree') {
        this.draw_tree(plotData, autoZoom);

    } else {
        throw "Bad projection_type value!";
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
