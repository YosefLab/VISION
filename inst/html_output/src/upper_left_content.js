/*
 * This script contains several components
 *
 * Upper_Left_Content
 *    - Signature_Table
 *    - Meta_Table
 *    - Gene_Select
 *
 */

function Upper_Left_Content()
{
    this.children = []
    this.sig_table = {}
    this.pc_table = {}
    this.gene_select = {}
    this.dom_node = document.getElementById("upper-left-content");
}

Upper_Left_Content.prototype.init = function()
{
    var sig_table = new Signature_Table()
    this.children.push(sig_table)
    this.sig_table = sig_table

    var sig_table_promise = sig_table.init();

    var pc_table = new Meta_Table()
    this.children.push(pc_table)
    this.pc_table = pc_table

    var pc_table_promise = pc_table.init();

    var gene_select = new Gene_Select()
    this.children.push(gene_select)
    this.gene_select = gene_select

    var gene_select_promise = gene_select.init()

    this.setLoadingStatus = createLoadingFunction(
        document.getElementById("upper-left-content")
    );

    // Initialize cluster dropdown in the top row
    var clust_dropdown = $(this.dom_node).find('#cluster-group-select');
    var cluster_variables = get_global_data('cluster_variables')

    clust_dropdown.empty();
    for(var i=0; i<cluster_variables.length; i++){
        clust_dropdown.append($("<option />")
            .val(cluster_variables[i])
            .text(cluster_variables[i]));
    }
    clust_dropdown
        .on('change', function () {
            set_global_status({
                'cluster_var':$(this).val(),
            });
        })


    return $.when(sig_table_promise, pc_table_promise, gene_select_promise);

}

Upper_Left_Content.prototype.update = function(updates)
{
    // Updates passed to children components
    _.each(this.children, function(child){
        child.update(updates)
    });

}

Upper_Left_Content.prototype.select_default_sig = function()
{
    var update = {}
    update['plotted_item'] = this.sig_table.get_top_sig()
    update['plotted_item_type'] = 'signature'
    return update;
}

Upper_Left_Content.prototype.hover_cells = function()
{
}


function Signature_Table()
{
    this.dom_node = document.getElementById("table-div-container");
    this.matrix = {}
    this.clusters = {}
    this.sorted_column = 'Consistency'
    this.filterSig = $(this.dom_node).find('#sig_filt_input')
    this.is_filtering = false
    this.is_collapsed = {} // cluster -> boolean, holds whether a cluster is collapsed
    this.tooltip = {} // will contain a Popper.js object
}


Signature_Table.prototype.init = function()
{

    var self = this;

    var debounced = _.debounce(function(){self.doneTyping()}, 500)
    self.filterSig.on('input', debounced)


    var popper = document.querySelector('body > #sig-tooltip')
    var arrow = popper.querySelector('.arrow')
    var initialNode = document.querySelector('body > #hidden-node')

    var outPopper = new Popper(initialNode, popper, {placement: 'top',
        modifiers: {
            arrow: {
                enabled: true,
                element: arrow,
            },
            preventOverflow: {
                enabled: false,
            },
            hide: {
                enabled: false,
            },
        }
    })

    this.tooltip = outPopper;

    // Set up the scrolling behavior
    var header_div = this.dom_node.querySelector('.sig-tables-header')
    var tables_div = this.dom_node.querySelector('.sig-tables-wrapper')
    $(tables_div).on('scroll', function() {
        $(header_div).scrollLeft($(this).scrollLeft())
    })

    var update_promise = self.update({})
    return update_promise;

}

Signature_Table.prototype.render = function()
{
    var self = this;
    var matrix = self.matrix;
    var main_vis = get_global_status('main_vis');

    var cluster_ids = _(self.clusters).values().uniq().sort(x => x).value()

    // If the filtering is enabled, we're just going to put everything into one cluster
    if (self.is_filtering){
        cluster_ids = [1]
    }

    // Create the Header row
    var header_row = $(self.dom_node).children(".sig-tables-header").find(".proj-row");
    header_row.find("th:not(:first-child)").remove()
    _.each(matrix.proj_labels, function(proj_label){
        var new_cell = $("<th>")
        var new_item = $("<div>")
        new_item.html("<div>" + proj_label + "</div>")
        new_cell.on("click", function() {
            var col_name = $(this).text()
            self.sorted_column = col_name
            self.render()
        });
        new_cell.append(new_item)
        header_row.append(new_cell)
    });

    // Add padding for scrollbar width of table below it
    $(self.dom_node).children(".sig-tables-header")
        .css("padding-right", detect_browser_scrollbar_width() + "px")

    var clusterTableDiv = $(self.dom_node).children('.sig-tables-wrapper').first()
    var new_table_divs = []

    // Remove old tables
    clusterTableDiv.children('.sig-table-div').remove()

    var curr_cl;
    for (var i = 0; i < cluster_ids.length; i++) {
        curr_cl = cluster_ids[i];

        if (!(curr_cl in self.is_collapsed)){
            self.is_collapsed[curr_cl] = true;
        }

        // Create new table and add to table_div
        var new_table_div = document.createElement("div");
        new_table_div.setAttribute("class", "sig-table-div");

        var new_table = document.createElement("table");
        new_table.setAttribute("id", "table"+ curr_cl);
        new_table.setAttribute("class", "sig-table");
        $(new_table).data("cluster", curr_cl)

        var thead = document.createElement("thead");

        var tbody = document.createElement("tbody");

        new_table.appendChild(thead);
        new_table.appendChild(tbody);


        new_table_div.appendChild(new_table);
        new_table_divs.push(new_table_div);

        if (typeof(matrix.sig_labels) == "string") {
            matrix.sig_labels = [matrix.sig_labels];
        }

        // Format cell data for better d3 binding
        var sig_labels = matrix.sig_labels.filter(
            function(x) { return self.is_filtering || self.clusters[x] == curr_cl; }
        );

        var zscores = []; // rows of z-scores
        var pvals = []; // rows of p-values
        var sig_indices = [] // signature index corresponding to each row

        for (var ind = 0; ind < sig_labels.length; ind++) {
            var sig = sig_labels[ind];
            var sig_index = matrix.sig_labels.indexOf(sig)

            pvals.push(matrix.pvals[sig_index]);
            zscores.push(matrix.zscores[sig_index]);
            sig_indices.push(sig_index);
        }

        var formatted_data_matrix = _.zip(sig_indices, zscores, pvals)
            .map(x => {
                var sig_index = x[0]
                var zscores_row = x[1]
                var pvals_row = x[2]
                return _.zip(zscores_row, pvals_row).map((x, j)  => {
                    return {
                        "row": sig_index, "col": j,
                        "zscore": x[0], "pval": x[1]
                    }
                })
            })

        var formatted_data_w_row_labels = d3.zip(sig_labels, formatted_data_matrix);

        // Sort data if necessary

        var sort_col = matrix.proj_labels.indexOf(self.sorted_column);
        if(sort_col > -1){
            var sortFun = function(a,b){
                return a[1][sort_col].zscore - b[1][sort_col].zscore; // ascending order
            };
            formatted_data_w_row_labels.sort(sortFun);
        }

        $(new_table_div).data('table-sort-val', 0)
        if (cluster_ids.length > 1) {
            // Pull out the cluster leader and put it first

            var row_vals;
            if (main_vis === "pcannotator") {
                row_vals = _.map(formatted_data_w_row_labels, x => {
                    return _(x[1]).map(y => y.pval).sum()
                })
            } else {
                row_vals = _.map(formatted_data_w_row_labels, x => {
                    return (x[1][0].zscore) // sort by zscore (geary-c actually)
                })
            }

            var leader_row_i = row_vals.indexOf(_.min(row_vals))
            var leader_row = formatted_data_w_row_labels[leader_row_i];
            formatted_data_w_row_labels.splice(leader_row_i, 1);
            formatted_data_w_row_labels.splice(0, 0, leader_row);

            if (sort_col > -1) {
                $(new_table_div).data('table-sort-val', leader_row[1][sort_col].zscore)
            }
        }

        if (main_vis == "pcannotator") {
            var colorScale = d3.scale.linear()
                .domain([-.5,0,.5])
                .range(["steelblue", "white", "lightcoral"])
                .clamp(true);
        } else {
            var colorScale = d3.scale.linear()
                .domain([1,.5,0])
                .range(["steelblue","white", "lightcoral"])
                .clamp(true);
        }
        var colorScaleCluster = d3.scale.linear()
            .domain([-10,0,10])
            .range(["steelblue","white", "lightcoral"])
            .clamp(true);

        var content_rows = d3.select(new_table).select('tbody').selectAll('tr')
            .data(formatted_data_w_row_labels);
        content_rows.enter().append('tr');
        content_rows.exit().remove();

        var content_row = content_rows.selectAll("td")
            .data(function(d){return [d[0]].concat(d[1]);})

        content_row.enter().append('td');
        content_row.exit().remove();

        if (main_vis === 'pcannotator') {
            content_row
                .filter(function(d,i) { return i > 0;})
                .style('background-color', function(d){return colorScale(d.zscore);})
                .on("click", function(d){tableClickFunction_PC(matrix.sig_labels[d.row], 'signature')})
                .append('div');

        } else {
            content_row
                .filter(function(d,i) { return i > 0;})
                .on("click", function(d){tableClickFunction_clusters(matrix.sig_labels[d.row], matrix.proj_labels[d.col], 'signature')})

            content_row
                .filter(function(d,i) { return i == 1;})
                .style('background-color', function(d){return colorScale(d.zscore);})
                .append('div')
                .text(function(d){ return d.zscore.toPrecision ? d.zscore.toPrecision(2) : d.zscore;})

            content_row
                .filter(function(d,i) { return i > 1;})
                .style('background-color', function(d){return colorScaleCluster(d.zscore);})
                .append('div')
                .text(function(d){ return d.pval < -1.301 ? '*' : '';})
        }

        // Hover actions
        var rowHoverFunc = function(header_row){
            return function(d, i){
                var tooltip_str;
                if(main_vis === 'pcannotator'){
                    tooltip_str = "corr = " + d.zscore.toFixed(2)
                } else {
                    if(main_vis === 'clusters' && i > 0)
                        tooltip_str = "p<" + _pval_format(d.pval)
                    else
                        tooltip_str = "z=" + d.zscore.toFixed(2) + ", p<" + _pval_format(d.pval)
                }
                createTooltip(self.tooltip, this, tooltip_str)
                hoverRowCol(header_row, this, matrix.proj_labels[d.col])
            }
        }(header_row, new_table);
        var rowUnHoverFunc = function(header_row){
            return function(d){
                destroyTooltip(self.tooltip)
                unhoverRowCol(header_row, this, matrix.proj_labels[d.col])
            }
        }(header_row, new_table);

        content_row
            .filter(function(d,i) { return i > 0;})
            .on("mouseenter", rowHoverFunc)
            .on("mouseleave", rowUnHoverFunc)

        // Add text for signature names
        content_row.filter(function(d,i) { return i == 0;})
            .html(function(d){return '<div class="sig-row-label"><div>'+d+'</div></div>';})

        // Add cluster expand/collapse behavior ONLY if more than one cluster
        if (cluster_ids.length > 1) {

            $(new_table).addClass('collapsible')

            if (self.is_collapsed[curr_cl]){
                $(new_table).addClass('collapsed')
            }

            $(new_table).children('tbody').children('tr:first-child').children("td:first-child")
                .on("click", function() { self.clickSummaryRow(this); });
        }

    }

    new_table_divs = _.sortBy(new_table_divs, d => {
        return $(d).data('table-sort-val')
    });

    _.forEach(new_table_divs, new_table_div => {
        clusterTableDiv.append(new_table_div);
    });

    // Apply compact styling for pcannotator
    if(main_vis === "pcannotator" || main_vis === "clusters" || main_vis === "tree"){
        $(self.dom_node).find('table').addClass('compact')
    } else {
        $(self.dom_node).find('table').removeClass('compact')
    }
    if(main_vis === "clusters" || main_vis === "tree"){
        $(self.dom_node).find('table').addClass('first-col')
    } else {
        $(self.dom_node).find('table').removeClass('first-col')
    }

    // Trigger existinging filter
    self.filterSig.trigger('input');
}

Signature_Table.prototype.clickSummaryRow = function(d){

    var table_id = $(d).closest('table')
    var cluster_id = table_id.data('cluster')

    if (this.is_collapsed[cluster_id]){
        this.is_collapsed[cluster_id] = false;
        table_id.removeClass('collapsed')
    } else {
        this.is_collapsed[cluster_id] = true;
        table_id.addClass('collapsed')
    }
}


Signature_Table.prototype.update = function(updates)
{
    var self = this;

    var clusters_promise;
    if (_.isEmpty(self.clusters)){

        clusters_promise = api.signature.clusters(false)
            .then(function(cls){
                self.clusters = cls;
                return true;
            })
    } else {
        clusters_promise = false
    }

    var matrix_promise;
    if( 'main_vis' in updates || _.isEmpty(self.matrix) || 'cluster_var' in updates){
        var main_vis = get_global_status('main_vis');
        var cluster_var = get_global_status('cluster_var');

        if (main_vis == "clusters") {
            matrix_promise = api.clusters.sigProjMatrix(cluster_var, false);
        } else if (main_vis == "tree") {
            matrix_promise = api.tree.sigProjMatrix(false);
        } else {
            matrix_promise = api.filterGroup.pCorr(false);
        }

        matrix_promise = matrix_promise
            .then(function(matrix) {
                self.matrix = matrix;
                return true;
            });
    } else {
        matrix_promise = false
    }

    return $.when(matrix_promise, clusters_promise)
        .then(function(matrix_update, clusters_update){
            if(matrix_update || clusters_update){
                self.render();
            }
        });

}

Signature_Table.prototype.get_top_sig = function()
{
    var matrix = this.matrix;
    var s_i = matrix.pvals.map(function(e){return e[0];}).argSort();

    return matrix.sig_labels[s_i[0]]
}


function Meta_Table()
{
    this.dom_node = document.getElementById("meta-table-div");
    this.matrix = {}
    this.sorted_column = 'All'
    this.is_collapsed = {} // cluster -> boolean, holds whether a cluster is collapsed
    this.tooltip = {} // will contain a Popper.js object
}

Meta_Table.prototype.init = function()
{
    var self = this;

    var popper = document.querySelector('body > #sig-tooltip')
    var arrow = popper.querySelector('.arrow')
    var initialNode = document.querySelector('body > #hidden-node')

    var outPopper = new Popper(initialNode, popper, {placement: 'top',
        modifiers: {
            arrow: {
                enabled: true,
                element: arrow,
            },
            preventOverflow: {
                enabled: false,
            },
            hide: {
                enabled: false,
            },
        }
    })

    this.tooltip = outPopper;

    // Set up the scrolling behavior
    var header_div = this.dom_node.querySelector('.sig-tables-header')
    var tables_div = this.dom_node.querySelector('.sig-tables-wrapper')
    $(tables_div).on('scroll', function() {
        $(header_div).scrollLeft($(this).scrollLeft())
    })

    var update_promise = self.update({})

    return update_promise;
}

Meta_Table.prototype.update = function(updates)
{
    var self = this;

    var clusters_promise;
    if (_.isEmpty(self.clusters)){

        clusters_promise = api.signature.clusters(true)
            .then(function(cls){
                self.clusters = cls;
                return true;
            })
    } else {
        clusters_promise = false
    }

    var matrix_promise;
    if('main_vis' in updates || _.isEmpty(self.matrix) || 'cluster_var' in updates){
        var main_vis = get_global_status('main_vis');
        var cluster_var = get_global_status('cluster_var');

        if (main_vis == "clusters") {
            matrix_promise = api.clusters.sigProjMatrix(cluster_var, true);
        } else if (main_vis === "tree") {
            matrix_promise = api.tree.sigProjMatrix(true);
        } else {
            matrix_promise = api.filterGroup.pCorr(true);
        }

        matrix_promise = matrix_promise.then(
            function(matrix){
                self.matrix = matrix
                return true
            });
    } else {
        matrix_promise = false
    }

    return $.when(matrix_promise, clusters_promise)
        .then(function(matrix_update, clusters_update){
            if(matrix_update || clusters_update){
                self.render();
            }
        });
}

Meta_Table.prototype.render = function()
{
    var self = this;
    var matrix = self.matrix;
    var main_vis = get_global_status('main_vis');

    var cluster_ids = _(self.clusters).values().uniq().sort(x => x).value()

    // Create the Header row
    var header_row = $(self.dom_node).children(".sig-tables-header").find(".proj-row");
    header_row.find("th:not(:first-child)").remove()
    _.each(matrix.proj_labels, function(proj_label){
        var new_cell = $("<th>")
        var new_item = $("<div>")
        new_item.html("<div>" + proj_label + "</div>")
        new_cell.on("click", function() {
            var col_name = $(this).text()
            self.sorted_column = col_name
            self.render()
        });
        new_cell.append(new_item)
        header_row.append(new_cell)
    });

    header_row
        .find("th:first-child")
        .on("click", function() {
            self.sorted_column = "Name"
            self.render()
        });

    // Add padding for scrollbar width of table below it
    $(self.dom_node).children(".sig-tables-header")
        .css("padding-right", detect_browser_scrollbar_width() + "px")

    if (typeof(matrix.sig_labels) == "string") {
        matrix.sig_labels = [matrix.sig_labels];
    }

    var clusterTableDiv = $(self.dom_node).children('.sig-tables-wrapper').first()
    var new_table_divs = []

    // Remove old tables
    clusterTableDiv.children('.sig-table-div').remove()

    var curr_cl;
    for (var i = 0; i < cluster_ids.length; i++) {
        curr_cl = cluster_ids[i];

        if (!(curr_cl in self.is_collapsed)){
            self.is_collapsed[curr_cl] = true;
        }

        // Create new table and add to table_div
        var new_table_div = document.createElement("div");
        new_table_div.setAttribute("class", "sig-table-div");

        var new_table = document.createElement("table");
        new_table.setAttribute("id", "table"+ curr_cl);
        new_table.setAttribute("class", "sig-table");
        $(new_table).data("cluster", curr_cl)

        var thead = document.createElement("thead");

        var tbody = document.createElement("tbody");

        new_table.appendChild(thead);
        new_table.appendChild(tbody);


        new_table_div.appendChild(new_table);
        new_table_divs.push(new_table_div);

        // Format cell data for better d3 binding
        var sig_labels = matrix.sig_labels.filter(
            function(x) { return self.is_filtering || self.clusters[x] == curr_cl; }
        );


        var zscores = []; // rows of z-scores
        var pvals = []; // rows of p-values
        var sig_indices = [] // signature index corresponding to each row

        for (var ind = 0; ind < sig_labels.length; ind++) {
            var sig = sig_labels[ind];
            var sig_index = matrix.sig_labels.indexOf(sig)

            pvals.push(matrix.pvals[sig_index]);
            zscores.push(matrix.zscores[sig_index]);
            sig_indices.push(sig_index);
        }

        var formatted_data_matrix = _.zip(sig_indices, zscores, pvals)
            .map(x => {
                var sig_index = x[0]
                var zscores_row = x[1]
                var pvals_row = x[2]
                return _.zip(zscores_row, pvals_row).map((x, j)  => {
                    return {
                        "row": sig_index, "col": j,
                        "zscore": x[0], "pval": x[1]
                    }
                })
            })

        var formatted_data_w_row_labels = d3.zip(sig_labels, formatted_data_matrix);

        // Sort data if necessary
        var sortFun;
        if(self.sorted_column === "Name"){
            sortFun = function(a, b){
                if(a < b) { return -1; }
                if(a > b) { return 1; }
                return 0;
            }
            formatted_data_w_row_labels.sort(sortFun);
        } else {
            var sort_col = matrix.proj_labels.indexOf(self.sorted_column);
            if(sort_col > -1){
                var sortFun = function(a,b){
                    return a[1][sort_col].zscore - b[1][sort_col].zscore; // Descending order
                };
                formatted_data_w_row_labels.sort(sortFun);
            }
        }


        $(new_table_div).data('table-sort-val', 0)
        if (cluster_ids.length > 1) {
            // Pull out the cluster leader and put it first
            // Just use the shortest name (works for factor clusters)
            var row_vals = _.map(formatted_data_w_row_labels, x => {
                return x[0].length
            })

            var leader_row_i = row_vals.indexOf(_.min(row_vals))
            var leader_row = formatted_data_w_row_labels[leader_row_i];
            formatted_data_w_row_labels.splice(leader_row_i, 1);
            formatted_data_w_row_labels.splice(0, 0, leader_row);

            if (sort_col > -1) {
                $(new_table_div).data('table-sort-val', leader_row[1][sort_col].zscore)
            }
        }

        var colorScale;
        if (main_vis == "pcannotator") {
            colorScale = d3.scale.linear()
                .domain([-.5,0,.5])
                .range(["steelblue", "white", "lightcoral"])
                .clamp(true);
        } else {
            colorScale = d3.scale.linear()
                .domain([1,.5,0])
                .range(["steelblue","white", "lightcoral"])
                .clamp(true);
        }
        var colorScaleCluster = d3.scale.linear()
            .domain([-10,0,10])
            .range(["steelblue","white", "lightcoral"])
            .clamp(true);


        var content_rows = d3.select(new_table).select('tbody').selectAll('tr')
            .data(formatted_data_w_row_labels);
        content_rows.enter().append('tr');
        content_rows.exit().remove();

        var content_row = content_rows.selectAll("td")
            .data(function(d){return [d[0]].concat(d[1]);})

        content_row.enter().append('td');
        content_row.exit().remove();


        if (main_vis === "pcannotator") {
            content_row
                .filter(function(d,i) { return i > 0;})
                .style('background-color', function(d){return colorScale(d.zscore);})
                .on("click", function(d){tableClickFunction_PC(matrix.sig_labels[d.row], 'meta')})
                .append('div');

        } else {
            content_row
                .filter(function(d,i) { return i > 0;})
                .on("click", function(d){tableClickFunction_clusters(matrix.sig_labels[d.row], matrix.proj_labels[d.col], 'meta')})

            content_row
                .filter(function(d,i) { return i == 1;})
                .style('background-color', function(d){return colorScale(d.zscore);})
                .append('div')
                .text(function(d){ return d.zscore.toPrecision ? d.zscore.toPrecision(2) : d.zscore;})

            content_row
                .filter(function(d,i) { return i > 1;})
                .style('background-color', function(d){return colorScaleCluster(d.zscore);})
                .append('div')
                .text(function(d){ return d.pval < -1.301 ? '*' : '';})
        }

        // Hover actions
        var rowHoverFunc = function(header_row){
            return function(d){
                var tooltip_str;
                if(main_vis === 'pcannotator'){
                    tooltip_str = "corr = " + d.zscore.toFixed(2)
                } else {
                    if(main_vis === 'clusters' && i > 0)
                        tooltip_str = "p<" + _pval_format(d.pval)
                    else
                        tooltip_str = "z=" + d.zscore.toFixed(2) + ", p<" + _pval_format(d.pval)
                }
                createTooltip(self.tooltip, this, tooltip_str)
                hoverRowCol(header_row, this, matrix.proj_labels[d.col])
            }
        }(header_row, new_table);
        var rowUnHoverFunc = function(header_row){
            return function(d){
                destroyTooltip(self.tooltip)
                unhoverRowCol(header_row, this, matrix.proj_labels[d.col])
            }
        }(header_row, new_table);

        content_row
            .filter(function(d,i) { return i > 0;})
            .on("mouseenter", rowHoverFunc)
            .on("mouseleave", rowUnHoverFunc)

        // Add text for signature names
        content_row.filter(function(d,i) { return i == 0;})
            .html(function(d){return '<div class="sig-row-label"><div>'+d+'</div></div>';})

        // Add cluster expand/collapse behavior ONLY if more than one cluster
        if (cluster_ids.length > 1 && curr_cl !== 1) {

            $(new_table).addClass('collapsible')

            if (self.is_collapsed[curr_cl]){
                $(new_table).addClass('collapsed')
            }

            $(new_table).children('tbody').children('tr:first-child').children("td:first-child")
                .on("click", function() { self.clickSummaryRow(this); });
        }
    }

    new_table_divs = _.sortBy(new_table_divs, d => {
        return $(d).data('table-sort-val')
    });

    _.forEach(new_table_divs, new_table_div => {
        clusterTableDiv.append(new_table_div);
    });


    // Apply compact styling for pcannotator
    if(main_vis === "pcannotator" || main_vis === "clusters" || main_vis === "tree"){
        $(self.dom_node).find('table').addClass('compact')
    } else {
        $(self.dom_node).find('table').removeClass('compact')
    }
    if(main_vis === "clusters" || main_vis === "tree"){
        $(self.dom_node).find('table').addClass('first-col')
    } else {
        $(self.dom_node).find('table').removeClass('first-col')
    }

}

Meta_Table.prototype.clickSummaryRow = function(d){

    var table_id = $(d).closest('table')
    var cluster_id = table_id.data('cluster')

    if (this.is_collapsed[cluster_id]){
        this.is_collapsed[cluster_id] = false;
        table_id.removeClass('collapsed')
    } else {
        this.is_collapsed[cluster_id] = true;
        table_id.addClass('collapsed')
    }
}

function Gene_Select()
{
    this.dom_node = document.getElementById("genes-table");
    this.recent_genes = [];
}

Gene_Select.prototype.init = function()
{
    var gene_promise = api.expression.genes.list()
        .then(function(genes){

            var geneSelect = $('#SelectGene');

            _.each(genes, function (gene) {
                geneSelect.append(
                    $('<option>', {
                        value: gene,
                        text: gene
                    }));
            });

            geneSelect.chosen({
                'width': '150px'
            })
                .on('change', function () {
                    set_global_status({
                        'plotted_item':$(this).val(),
                        'plotted_item_type': 'gene'
                    });
                });
        });


    this.render_recent_genes()

    return $.when(gene_promise);

}

Gene_Select.prototype.update = function(updates)
{
    // Update the 'recent-genes' list
    if('plotted_item' in updates){
        var gene = get_global_status('plotted_item')
        var plotted_item_type = get_global_status('plotted_item_type')
        if(plotted_item_type === 'gene'){
            if(this.recent_genes.indexOf(gene) === -1){
                this.recent_genes.unshift(gene)

                if(this.recent_genes.length > 10){
                    this.recent_genes.pop();
                }

                this.render_recent_genes()
            }
        }
    }
}

Gene_Select.prototype.render_recent_genes = function()
{
    var recent_genes_div = $(this.dom_node).find('#recent-genes-div')
    var recent_genes_list = $(this.dom_node).find('#recent-genes-list')

    if(this.recent_genes.length == 0){
        recent_genes_div.hide()
        return
    }

    recent_genes_list.children().remove()

    this.recent_genes.forEach(function (gene) {
        recent_genes_list.append(
            $('<div>',{'class': 'recent-gene'}).text(gene)
        )
    })

    recent_genes_list.children().on('click', function() {
        var gene = $(this).text()
        set_global_status({
            'plotted_item_type': 'gene',
            'plotted_item': gene
        })
    });

    recent_genes_div.show()
}

Signature_Table.prototype.doneTyping = function()
{
    var self = this;
    var val = self.filterSig.val();

    if (val.length > 0 && !self.is_filtering) {
        // Just started filtering, need to re-render table
        self.is_filtering = true
        self.render()
        return // No need to finish - this function is called at the end of render()
    }

    if (val.length == 0 && self.is_filtering) {
        // Just ended filtering, need to re-render table
        self.is_filtering = false
        self.render()
        return // No need to finish - this function is called at the end of render()
    }

    var vals = val.split(",");
    vals = vals.map(function(str){return str.trim();})
        .filter(function(str){ return str.length > 0;})
        .map(function(str){ return str.toLowerCase();})

    var tablerows = $(self.dom_node).find('.sig-table-div').find('tr')
    tablerows.removeClass('filtered');

    var posvals = vals.filter(function(str){ return str[0] != '!';});
    var negvals = vals.filter(function(str){ return str[0] == '!';})
        .map(function(str){ return str.slice(1);})
        .filter( function(str){return str.length > 0;});

    if(posvals.length > 0){
        tablerows.filter(function(i, element){
            var sig_text = $(element).children('td').first().html().toLowerCase();
            for(var j = 0; j < posvals.length; j++)
            {
                if(sig_text.indexOf(posvals[j]) > -1)
                {
                    return false;
                }
            }
            return true;
        }).addClass('filtered');
    }

    if(negvals.length > 0){
        tablerows.filter(function(i, element){
            var sig_text = $(element).children('td').first().html().toLowerCase();
            for(var j = 0; j < negvals.length; j++)
            {
                if(sig_text.indexOf(negvals[j]) > -1)
                {
                    return true;
                }
            }
            return false;
        }).addClass('filtered');
    }
}


// Function that's triggered when clicking on table cell
function tableClickFunction_PC(row_key, item_type)
{
    var update = {}
    update['plotted_item_type'] = item_type;
    update['plotted_item'] = row_key;

    set_global_status(update);
}

function tableClickFunction_clusters(row_key, col_key)
{
    var update = {}

    var meta_sigs = get_global_data('meta_sigs')

    var item_type
    if (meta_sigs.indexOf(row_key) > -1) {
        item_type = 'meta'
    } else {
        item_type = 'signature'
    }

    update['plotted_item_type'] = item_type;
    update['plotted_item'] = row_key;

    if (col_key === 'Consistency')
    {
        col_key = '' // This is used to indicate 'no cluster'
    }
    update['selected_cluster'] = col_key;

    set_global_status(update);
}


function hoverRowCol(header_row, node, col){
    $(header_row).find('th').filter((i, e) => $(e).text() == col).addClass('highlight')
    $(node).siblings('td:first-child').addClass('highlight')


    var hovered_cells;
    if (col === 'Consistency'){
        hovered_cells = [];
    } else {
        var clusters = get_global_data('clusters'); // clusters maps cell_id to cluster
        hovered_cells = _(clusters)
            .pickBy(val => val === col)
            .keys()
            .value();
    }

    var event = new CustomEvent('hover-cells',
        { detail: hovered_cells, bubbles: true } )

    node.dispatchEvent(event)
}

function unhoverRowCol(header_row, node, col){
    $(header_row).find('th').filter((i, e) => $(e).text() == col).removeClass('highlight')
    $(node).siblings('td:first-child').removeClass('highlight')

    var hovered_cells = [];
    var event = new CustomEvent('hover-cells',
        { detail: hovered_cells, bubbles: true } )

    node.dispatchEvent(event)
}

function createTooltip(popper, node, text) {
    var popper_node = popper.popper
    var inner_node = popper_node.querySelector('.inner')
    inner_node.textContent = text

    popper.reference = node;
    popper.update();
}

function destroyTooltip(popper) {
    var node = document.querySelector('body #hidden-node')
    popper.reference = node;
    popper.update();
}

function _pval_format(x) {
    var p = Math.pow(10, x)
    if(p < 0.01){
        return p.toExponential(1)
    } else {
        return p.toFixed(2)
    }
}
