/*
 * This script contains several components
 *
 * Upper_Left_Content
 *    - Signature_Table
 *    - Meta_Table
 *    - Gene_Select
 *    - DE_Select
 *
 */

function Upper_Left_Content()
{
    this.children = []
    this.sig_table = {}
    this.protein_table = {}
    this.pc_table = {}
    this.gene_select = {}
    this.dendrogram = {}
    this.mod_table = {}
    this.dom_node = document.getElementById("upper-left-content");
}

Upper_Left_Content.prototype.init = function()
{
    var self = this

    // Initialize cluster dropdown in the top row
    // Must happen before initializing child components
    var clust_dropdown = $(this.dom_node).find('#cluster-group-select');
    var cluster_variables = Object.keys(get_global_data('meta_levels'))

    clust_dropdown.empty();
    for(var i=0; i<cluster_variables.length; i++){
        clust_dropdown.append($("<option />")
            .val(cluster_variables[i])
            .text(cluster_variables[i]));
    }
    clust_dropdown
        .on('change', function () {
            self.update({
                'cluster_var':$(this).val(),
            });
            $(this).blur()
        })

    // 2-level variables are boring
    var meta_levels = get_global_data('meta_levels')
    for (var i=0; i < cluster_variables.length; i++) {
        var cv = cluster_variables[i]
        if (meta_levels[cv].length > 2) {
            clust_dropdown.val(cv);
            break;
        }
    }


    var sig_table = new Signature_Table()
    this.children.push(sig_table)
    this.sig_table = sig_table

    var sig_table_promise = sig_table.init();

    var protein_table = new Protein_Table()
    this.children.push(protein_table)
    this.protein_table = protein_table

    var protein_table_promise = protein_table.init();

    var pc_table = new Meta_Table()
    this.children.push(pc_table)
    this.pc_table = pc_table

    var pc_table_promise = pc_table.init();

    var gene_select = new Gene_Select()
    this.children.push(gene_select)
    this.gene_select = gene_select

    var gene_select_promise = gene_select.init()


    // De select Promise
    var de_select = new DE_Select()
    this.children.push(de_select)
    this.de_select = de_select

    var de_select_promise = de_select.init()
    
    var has_dend = get_global_status('has_dendrogram')
    // Dend Promise
    if (has_dend) {
      var dend = new Dend()
      this.children.push(dend)
      this.dend = dend
  
      var dend_promise = dend.init()
    }
    
    
    var mod_table = new Modules_Table()
    this.children.push(mod_table)
    this.mod_table = mod_table

    var mod_table_promise = mod_table.init();
    
    this.setLoadingStatus = createLoadingFunction(
        document.getElementById("upper-left-content")
    );

    // If no signatures, need to hide the 'Signatures' nav item
    // If no proteins, need to hide the 'Proteins' nav item
    // If no dendrogram, hide that nav item
    // same for modules

    var has_sigs = get_global_status('has_sigs')
    var has_mods = get_global_status('has_mods')
    var has_proteins = get_global_status('has_proteins')
    var nav = $(this.dom_node).find('.section-nav')
    var dendrogramTab = $('#dend-section')
    var sigTableTab = nav.find('#sig-table-tab')
    var modTableTab = nav.find('#modules-table-tab')
    var proteinTableTab = nav.find('#protein-table-tab')
    var metaTableTab = nav.find('#meta-table-tab')
    
    var has_dendrogram = get_global_status('has_dendrogram')
    
    if(!has_sigs){
        sigTableTab.parent().hide()
        if(!has_proteins){
            metaTableTab.click()
        } else {
            proteinTableTab.click()
        }
    }
    if(!has_proteins){
        proteinTableTab.parent().hide()
    }
    if(!has_mods){
        modTableTab.parent().hide()
    }
    if(!has_dendrogram){
        dendrogramTab.hide()
    }

    $(this.dom_node).find('.nav-link')
        .on('click', function(e) {
            if (e.target.id === "genes-table-tab" || e.target.id === "de-table-tab") {
                clust_dropdown.prop('disabled', true)
            } else {
                clust_dropdown.prop('disabled', false)
            }
        });


    return $.when(sig_table_promise, protein_table_promise, pc_table_promise, gene_select_promise, de_select_promise, dend_promise, mod_table_promise);

}

var _get_cluster_var = function()
{
    // Kind of a hack
    var clust_dropdown = $('#cluster-group-select');
    return clust_dropdown.val();
}

Upper_Left_Content.prototype.update = function(updates)
{
    // Updates passed to children components

    var child_promises = []

    _.each(this.children, function(child){
        child_promises.push(child.update(updates))
    });

    return $.when.apply($, child_promises)

}

Upper_Left_Content.prototype.select_default_sig = function()
{
    var update = {}

    var has_sigs = get_global_status('has_sigs')
    var has_proteins = get_global_status('has_proteins')

    var plotted_item
    var item_type
    if (has_sigs){
        plotted_item = this.sig_table.get_top_sig()
        item_type = 'signature'
    } else if (has_proteins) {
        plotted_item = this.protein_table.get_top_sig()
        item_type = 'protein'
    } else {
        plotted_item = this.pc_table.get_top_sig()
        item_type = 'meta'
    }

    if (item_type == "signature") {
        // Could also be meta-data because we mix those in
        var meta_sigs = get_global_data('meta_sigs')
        if (meta_sigs.indexOf(plotted_item) > -1) {
            item_type = 'meta'
        }
    }

    update['plotted_item']  = plotted_item
    update['plotted_item_type'] = item_type

    return update;
}


function Signature_Table()
{
    this.dom_node = document.getElementById("table-div-container");
    this.matrix = {}
    this.clusters = {}
    this.sorted_column = 'Score'
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
                return b[1][sort_col].zscore - a[1][sort_col].zscore; // Descending order
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
                    return (x[1][0].zscore*-1) // sort by zscore (geary-c actually)
                })
            }

            var leader_row_i = row_vals.indexOf(_.min(row_vals))
            var leader_row = formatted_data_w_row_labels[leader_row_i];
            formatted_data_w_row_labels.splice(leader_row_i, 1);
            formatted_data_w_row_labels.splice(0, 0, leader_row);

            if (sort_col > -1) {
                $(new_table_div).data('table-sort-val', leader_row[1][sort_col].zscore*-1)
            }
        }

        if (main_vis == "pcannotator") {
            var colorScale = d3.scale.linear()
                .domain([-.5,0,.5])
                .range(["steelblue", "white", "lightcoral"])
                .clamp(true);
        } else {
            var colorScale = d3.scale.linear()
                .domain([0,.5,1])
                .range(["steelblue","white", "lightcoral"])
                .clamp(true);
        }
        var colorScaleCluster = d3.scale.linear()
            .domain([0,0.5,1])
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
                .on("click", function(d){tableClickFunction_clusters(matrix.sig_labels[d.row], 'signature')})

            content_row
                .filter(function(d,i) { return i == 1;})
                .style('background-color', function(d){return colorScale(d.zscore);})
                .append('div')
                .text(function(d){ return d.zscore.toPrecision ? d.zscore.toPrecision(2) : d.zscore;})

            content_row
                .filter(function(d,i) { return i > 1;})
                .style('background-color', function(d){return colorScaleCluster(d.zscore);})
                .append('div')
                .text(function(d){ return d.pval < 0.05 ? '*' : '';})
        }

        // Hover actions
        var rowHoverFunc = function(header_row){
            return function(d, i){
                var tooltip_str;
                if(main_vis === 'pcannotator'){
                    tooltip_str = "corr = " + d.zscore.toFixed(2)
                } else {
                    if(main_vis === 'clusters' && i > 0)
                        tooltip_str = "AUC=" + _auc_format(d.zscore) + " p<" + _pval_format(d.pval)
                    else
                        tooltip_str = "C'=" + d.zscore.toFixed(2) + ", p<" + _pval_format(d.pval)
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

    // If there are no signatures, don't update
    var has_sigs = get_global_status('has_sigs')
    if(!has_sigs){
        return $.when(null); // Empty promise
    }

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
        var cluster_var = _get_cluster_var();

        var fix_col_label = true
        if (main_vis == "clusters") {
            matrix_promise = api.clusters.sigProjMatrix(cluster_var, false);
        } else if (main_vis == "tree") {
            matrix_promise = api.tree.sigProjMatrix(false);
        } else {
            matrix_promise = api.pCorr.signatures();
            fix_col_label = false
        }

        matrix_promise = matrix_promise
            .then(function(matrix) {
                if (fix_col_label){
                    matrix.proj_labels[0] = "Score";
                }
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
    var s_i = matrix.zscores.map(function(e){return e[0];}).argMax();

    return matrix.sig_labels[s_i]
}


function Protein_Table()
{
    this.dom_node = document.getElementById("protein-div-container");
    this.matrix = {}
    this.clusters = {}
    this.sorted_column = 'Score'
    this.filterProtein = $(this.dom_node).find('#protein_filt_input')
    this.is_filtering = false
    this.is_collapsed = {} // cluster -> boolean, holds whether a cluster is collapsed
    this.tooltip = {} // will contain a Popper.js object
}


Protein_Table.prototype.init = function()
{

    var self = this;

    var debounced = _.debounce(function(){self.doneTyping()}, 500)
    self.filterProtein.on('input', debounced)


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

Protein_Table.prototype.render = function()
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
                return b[1][sort_col].zscore - a[1][sort_col].zscore; // Descending order
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
                    return (x[1][0].zscore*-1) // sort by zscore (geary-c actually)
                })
            }

            var leader_row_i = row_vals.indexOf(_.min(row_vals))
            var leader_row = formatted_data_w_row_labels[leader_row_i];
            formatted_data_w_row_labels.splice(leader_row_i, 1);
            formatted_data_w_row_labels.splice(0, 0, leader_row);

            if (sort_col > -1) {
                $(new_table_div).data('table-sort-val', leader_row[1][sort_col].zscore*-1)
            }
        }

        if (main_vis == "pcannotator") {
            var colorScale = d3.scale.linear()
                .domain([-.5,0,.5])
                .range(["steelblue", "white", "lightcoral"])
                .clamp(true);
        } else {
            var colorScale = d3.scale.linear()
                .domain([0,.5,1])
                .range(["steelblue","white", "lightcoral"])
                .clamp(true);
        }
        var colorScaleCluster = d3.scale.linear()
            .domain([0,0.5,1])
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
                .on("click", function(d){tableClickFunction_PC(matrix.sig_labels[d.row], 'protein')})
                .append('div');

        } else {
            content_row
                .filter(function(d,i) { return i > 0;})
                .on("click", function(d){tableClickFunction_clusters(matrix.sig_labels[d.row], 'protein')})

            content_row
                .filter(function(d,i) { return i == 1;})
                .style('background-color', function(d){return colorScale(d.zscore);})
                .append('div')
                .text(function(d){ return d.zscore.toPrecision ? d.zscore.toPrecision(2) : d.zscore;})

            content_row
                .filter(function(d,i) { return i > 1;})
                .style('background-color', function(d){return colorScaleCluster(d.zscore);})
                .append('div')
        }

        // Hover actions
        var rowHoverFunc = function(header_row){
            return function(d, i){
                var tooltip_str;
                if(main_vis === 'pcannotator'){
                    tooltip_str = "corr = " + d.zscore.toFixed(2)
                } else {
                    if(main_vis === 'clusters' && i > 0)
                        tooltip_str = "AUC=" + _auc_format(d.zscore) + " p<" + _pval_format(d.pval)
                    else
                        tooltip_str = "C'=" + d.zscore.toFixed(2)
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
    self.filterProtein.trigger('input');
}

Protein_Table.prototype.clickSummaryRow = function(d){

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


Protein_Table.prototype.update = function(updates)
{
    var self = this;

    // If there are no proteins, don't update
    var has_proteins = get_global_status('has_proteins')
    if(!has_proteins){
        return $.when(null); // Empty promise
    }

    var clusters_promise;
    if (_.isEmpty(self.clusters)){

        clusters_promise = api.protein.clusters()
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
        var cluster_var = _get_cluster_var()

        var fix_col_label = true
        if (main_vis == "clusters") {
            matrix_promise = api.clusters.proteinMatrix(cluster_var);
        } else if (main_vis == "tree") {
            matrix_promise = api.tree.proteinMatrix();
        } else {
            matrix_promise = api.pCorr.proteins();
            fix_col_label = false
        }

        matrix_promise = matrix_promise
            .then(function(matrix) {
                if (fix_col_label){
                    matrix.proj_labels[0] = "Score";
                }
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

Protein_Table.prototype.get_top_sig = function()
{
    var matrix = this.matrix;
    var s_i = matrix.zscores.map(function(e){return e[0];}).argMax();

    return matrix.sig_labels[s_i]
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
        var cluster_var = _get_cluster_var()

        var fix_col_label = true
        if (main_vis == "clusters") {
            matrix_promise = api.clusters.sigProjMatrix(cluster_var, true);
        } else if (main_vis === "tree") {
            matrix_promise = api.tree.sigProjMatrix(true);
        } else {
            matrix_promise = api.pCorr.meta();
            fix_col_label = false
        }

        matrix_promise = matrix_promise.then(
            function(matrix){
                if (fix_col_label){
                    matrix.proj_labels[0] = "Score";
                }
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

        // Sort data within each table group
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
                    return b[1][sort_col].zscore - a[1][sort_col].zscore; // Descending order
                };
                formatted_data_w_row_labels.sort(sortFun);
            }
        }


        // Sort the row groups themselves
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
                $(new_table_div).data('table-sort-val', leader_row[1][sort_col].zscore*-1)
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
                .domain([0,.5,1])
                .range(["steelblue","white", "lightcoral"])
                .clamp(true);
        }
        var colorScaleCluster = d3.scale.linear()
            .domain([0,0.5,1])
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
                .on("click", function(d){tableClickFunction_clusters(matrix.sig_labels[d.row], 'meta')})

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
                    tooltip_str = "p<" + _pval_format(d.pval)
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

Meta_Table.prototype.get_top_sig = function()
{
    var matrix = this.matrix;
    var s_i = matrix.zscores.map(function(e){return e[0];}).argMax();

    return matrix.sig_labels[s_i]
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
                'width': '150px',
                'max_shown_results': 1000,
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

                if(this.recent_genes.length > 20){
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



function DE_Select()
{
    //this.recent_genes = [];
    this.dom_node = document.getElementById("de-table");
    this.filterSig = $(this.dom_node).find('#sig_filt_input')
    this.numControl = $(this.dom_node).find('#num_control');
    this.numSelect = $(this.dom_node).find('#num');
    this.denomSelect = $(this.dom_node).find('#denom');
    this.submit_de = $(this.dom_node).find('#submit_de');
    this.new_de = $(this.dom_node).find('#new_de');
    this.de_error = $(this.dom_node).find('#de-error');
    this.de_results = $(this.dom_node).find('#de-results');
    this.de_table = $(this.dom_node).find('#de-results-table');
    this.de_controls = $(this.dom_node).find('.controls');

    this.min_cells = $(this.dom_node).find('#de_min_cells');
    this.sub_cells_n = $(this.dom_node).find('#de_sub_cells_n');
    this.sub_cells = $(this.dom_node).find('#de_sub_cells');
}

DE_Select.prototype.init = function()
{

    var self = this;

    var numControl = this.numControl;
    var numSelect = this.numSelect;
    var denomSelect = this.denomSelect;
    var submit_de = this.submit_de;
    var new_de = this.new_de;
    var de_error = this.de_error;
    var de_results = this.de_results;
    var de_controls = this.de_controls;

    var clusters = Object.keys(get_global_data('meta_levels'))

    this.setLoadingStatus = createLoadingFunction(this.dom_node);

    var table = this.de_table.DataTable(
        {
            'columns': [
                {'title': 'Feature'},
                {'title': 'Type'},
                {'title': 'logFC'},
                {'title': 'AUC'},
                {'title': 'FDR'},
            ],
            "pageLength": 500,
            "scrollY": '15vh',
            "order": [[4, "asc"]],
            "buttons": {
                dom: {
                    button: {
                        className: ''
                    }
                },
                buttons: [
                    {
                        extend: 'csvHtml5',
                        text: 'Export',
                        className: 'btn btn-secondary',
                    }
                ]
            },
            "pagingType": "simple",
            "dom": '<"fdr-filter">ftip'
        }
    );

    table.buttons().container().prependTo($(this.dom_node).find('.button-container'))

    $(this.dom_node).find('.fdr-holder').remove().appendTo(
        $(this.dom_node).find('.fdr-filter')
    )

    $.fn.dataTable.ext.search.push(
        function( settings, data) {
            var min = 0;
            var max = parseFloat( $('#max').val(), 10 );
            var fdr = parseFloat( data[4] ) || 0; // use data for the age column
            if ( ( isNaN( min ) && isNaN( max ) ) ||
                 ( isNaN( min ) && fdr <= max ) ||
                 ( min <= fdr   && isNaN( max ) ) ||
                 ( min <= fdr   && fdr <= max ) )
            {
                return true;
            }
            return false;
        }
    );

    $('#max').off();
    $('#max').keyup( function() {
        $('#de-results-table').DataTable().draw();
    } );


    this.de_table.on('click', "tr", function() {
        var cells = $(this).find("td")
        var feature = cells.eq(0).text()
        var type = cells.eq(1).text()

        if (type === 'Gene'){
            set_global_status({
                'plotted_item_type': 'gene',
                'plotted_item': feature,
            })
        } else {
            set_global_status({
                'plotted_item_type': 'protein',
                'plotted_item': feature,
            })
        }
        //Plotly.restyle("scatter-div", {selectedpoints: [null]});
    });



    numControl.append(
        $('<option>', {
            value: "selections",
            text: "Selections"
        }));

    _.each(clusters, name => {
        numControl.append(
            $('<option>', {
                value: name,
                text: name,
            }));
    });


    numSelect.chosen({
        'width': '150px',
        'max_shown_results': 1000,
        'placeholder_text_single': 'N/A',
        'disable_search': true,
    });

    denomSelect.chosen({
        'width': '150px',
        'max_shown_results': 1000,
        'placeholder_text_single': 'N/A',
        'disable_search': true,
    });

    function addClusters(select, cluster, num) {
        var cluster_vals = get_global_data("meta_levels")[cluster];
        var data = Array.from(new Set(Object.values(cluster_vals)));

        if (!num) {
            select.append(
                $('<option>', {
                    value: 'remainder',
                    text: 'Remainder',
                })
            );
        }

        _.each(data, name => {
            select.append(
                $('<option>', {
                    value: name,
                    text: name
                }));
        });

        select.trigger('chosen:updated')
    }

    submit_de.on('click', function () {

        de_error.hide();

        var type_n, group_num;
        var type_d, group_denom;

        var subtype_n = ""
        var subtype_d = ""

        var min_cells = self.min_cells.val();
        var sub_cells_n = self.sub_cells_n.val();
        var sub_cells = self.sub_cells.is(":checked")

        if (!min_cells) {
            min_cells = 0;
        }

        if (numControl.val() === 'selections') {

            if(numSelect.val() === 'current') {
                type_n = 'current';
                group_num = get_global_status("selected_cell");
                if (group_num.length === 0){
                    de_error.html("No cells currently selected<br>Make a selection using the lasso tool on the scatter plot");
                    de_error.show();
                    return;
                }
            } else {

                type_n = 'saved_selection';
                group_num = numSelect.val();
                if (group_num === null){
                    de_error.html("No saved selections<br>To save a selection, first select cells on the scatter plot using the lasso tool.  Then click 'Save' above the plot and give the selection a name.");
                    de_error.show();
                    return;
                }
            }

        } else {
            type_n = 'meta';
            subtype_n = numControl.val();
            group_num = numSelect.val();
        }

        if (denomSelect.val() == 'remainder') {
            type_d = 'remainder';
            group_denom = '';
        } else if (numControl.val() === 'selections') {
            type_d = 'saved_selection';
            group_denom = denomSelect.val();
            if (group_denom === null){
                de_error.html("No saved selections<br>To save a selection, first select cells on the scatter plot using the lasso tool.  Then click 'Save' above the plot and give the selection a name.");
                de_error.show();
                return;
            }
        } else {
            type_d = 'meta';
            subtype_d = subtype_n;
            group_denom = denomSelect.val();
        }


        self.setLoadingStatus(true, 0);

        api.de(
            type_n, subtype_n, group_num,
            type_d, subtype_d, group_denom,
            min_cells, sub_cells, sub_cells_n
        ).then(data => {
            set_global_status({"de":data})
        });

    });

    new_de.on('click', function () {
        de_results.hide()
        de_controls.show()
    })


    numControl.chosen({
        'width': '150px',
        'max_shown_results': 1000,
        'disable_search': true,
    }).on('change', function () {

        de_error.hide();

        var val = $(this).val();

        numSelect = self.numSelect
        numSelect.children().remove()

        denomSelect = self.denomSelect
        denomSelect.children().remove()

        if (val === "selections") {
            // add selection values
            self.addSelections(numSelect, denomSelect);
        }  else {
            // genotype
            addClusters(numSelect, val, true)
            addClusters(denomSelect, val, false)
        }
    });

    this.sub_cells.on('change', function() {
        var checked = $(this).is(":checked")
        if (checked) {
            $(this).parent().siblings('.max-control').removeClass('disabled')
        } else {
            $(this).parent().siblings('.max-control').addClass('disabled')
        }
    })


    this.addSelections(numSelect, denomSelect);

}

DE_Select.prototype.addSelections = function(selectNum, selectDenom) {

    return api.cells.listSelections().then(data => {

        selectNum.append(
            $('<option>', {
                value: 'current',
                text: 'Current',
            })
        );
        selectDenom.append(
            $('<option>', {
                value: 'remainder',
                text: 'Remainder',
            })
        );

        _.each(data, name => {
            selectNum.append(
                $('<option>', {
                    value: name,
                    text: name
                }));
        });

        _.each(data, name => {
            selectDenom.append(
                $('<option>', {
                    value: name,
                    text: name
                }));
        });

        selectNum.trigger('chosen:updated')
        selectDenom.trigger('chosen:updated')

    });

}


DE_Select.prototype.update = function(updates)
{
    // Update the 'recent-genes' list
    if('de' in updates){
        this.de_controls.hide();
        this.render_de();
    }

    if ('selected_cell' in updates || 'selection_type' in updates){
        this.de_error.hide();
    }

    if ('selections_changed' in updates){
        var numControl = this.numControl;
        var val = $(numControl).val();

        if (val === "selections") {
            // add selection values

            var numSelect = this.numSelect
            numSelect.children().remove()

            var denomSelect = this.denomSelect
            denomSelect.children().remove()
            this.addSelections(numSelect, denomSelect);
        }
    }
}

DE_Select.prototype.render_de = function()
{

    var features = get_global_status('de')["Feature"]
    var types = get_global_status('de')["Type"]
    var pvals = get_global_status('de')["pval"]
    var stats = get_global_status('de')["stat"]
    var logFCs = get_global_status('de')["logFC"]

    var dataSet = _.map(features, (f, i) => [features[i], types[i], logFCs[i], stats[i], pvals[i]]);

    this.de_results.show();

    this.de_table.DataTable().clear()
        .rows.add(dataSet)
        .columns.adjust() // Needed or else column headers don't align
        .draw()


    // Render the title
    var group = this.numControl.find('option:selected').text()
    var numText = this.numSelect.find('option:selected').text()
    var denomText = this.denomSelect.find('option:selected').text()

    var title = group + ": " + numText + " vs. " + denomText

    $(this.de_results).find('#de-results-title').text(title)

    this.setLoadingStatus(false);
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

Protein_Table.prototype.doneTyping = function()
{
    var self = this;
    var val = self.filterProtein.val();

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



function Dend()
{
    //this.recent_genes = [];
    this.dom_node = document.getElementById("dend");
    this.cluster_var_idx = 0;

}

Dend.prototype.init = function()
{

    var self = this;

    this.newick = get_global_status("dendrogram")
    this.tree = null


}


Dend.prototype.update = function(updates)
{
    if ('plotted_item' in updates){
        this.plotted_item = updates['plotted_item']

        this.render_dend();
    }
    if ('selected_cell' in updates) {
      this.update_tree_selection()
    }
}

Dend.prototype.update_tree_selection = function() {
    var cells = get_global_status("selected_cell")
    this.tree.selection_callback(function (nodes) {});
    var nodes = cells
    var all_nodes =  this.tree.get_nodes()
    
    if (cells.length > 0) {
        var node_selector = function(node) {
          if (!d3.layout.phylotree.is_leafnode(node)) {
            return true
          }
          var selected = nodes.includes(node.name)

          if (selected) {
            nodes.push(node.name)
          }
          return false
        }
       all_nodes = all_nodes.filter(node_selector)
    
        while (true) {
          var node_selector = function(node) {
              var in_nodes = (n) => nodes.includes(n.name)
              var selected = node.children.every(in_nodes)
              if (selected) {
                nodes.push(node.name)
              }
              return !selected
          }
          
          var all_nodes_new = all_nodes.filter(node_selector)
          if (all_nodes_new.length === all_nodes.length) {
            break
          }
          all_nodes = all_nodes_new
        }
    }
    
    var node_selector_final = function(node) {
          var selected = nodes.includes(node.target.name)
          return selected
      }
    
    this.tree.modify_selection(node_selector_final)

    this.tree.selection_callback(function (nodes) {
        var cells = nodes.map(function(node) {
          if (d3.layout.phylotree.is_leafnode(node)) {
            return node.name;
          }
        })
        
        if (cells.length > 0) {
          set_global_status({"selected_cell":cells, "selection_type":"cells"})
        }
      });
}

Dend.prototype.render_dend = function()
{
    this.setLoadingStatus = createLoadingFunction(this.dom_node);
    this.setLoadingStatus(true, 0);
    
    $("#dend_layout").prop('checked', true);
    
    var plotted = get_global_data("plotted_values");
    var plotted_item = this.plotted_item;
    
    function edgeStyler(dom_element, edge_object) {
      dom_element.style("stroke-width", "1pt");
    }
    
    var num_cells = Object.keys(plotted).length;
    var width = 720;
    var height = 720
    var container_id = '#tree_container';

    var s = num_cells*4;
    $(container_id).empty();
    var svg = d3.select(container_id).append("svg")
    .attr("width", width)
    .attr("height", height);

    var tree = d3.layout.phylotree()
    // create a tree layout object
    .svg(svg)
    .radial(true)
    .options({
      "brush": false,
      "zoom": true,
      "show-scale": false,
      "selectable": true,
      "minimum-per-node-spacing": 5,
      "minimum-per-level-spacing": 5,
      "max-radius": s/2,
      'left-right-spacing': 'fit-to-size'
      //'top-bottom-spacing': 'fit-to-size',
      //'max-radius': 800
    }).size([s, s])

    tree(this.newick)

    var tree_attributes = {};
    var scatter = right_content.layoutPlotData["full"][0].scatter;
    var plotly_defaults_colors = ["1F77B4", "FF7F0E", "2CA02C", "D62728", "9467BD", "8C564B","E377C2", "7F7F7F", "BCBD22", "17BECF"]
    
    if(get_global_data('meta_levels')[plotted_item] == null)  {
        var full_color_range = scatter.full_color_range;
        var diverging_colormap = scatter.diverging_colormap
        var c = Object.values(plotted);
        var cvals = c.filter(cv => cv !== "NA")
        cvals.sort(d3.ascending); // Needed for quantile
        var low, high;
        if(full_color_range){
            low = d3.min(cvals)
            high = d3.max(cvals)
        } else {
            low = d3.quantile(cvals, 0.02)
            high = d3.quantile(cvals, 0.98)
        }
        if(get_global_status("plotted_item_type") !== "gene") {
            var color_range = [
              '#440154', '#48186a', '#472d7b', '#424086', '#3b528b', '#33638d', '#2c728e', '#26828e', 
            '#21918c', '#1fa088', '#28ae80', '#3fbc73', '#5ec962', '#84d44b', '#addc30',
            '#d8e219', '#fde725'];
            var color_domain = [
              0, 0.06274509803921569, 0.12549019607843137, 0.18823529411764706, 0.25098039215686274, 
              0.3137254901960784, 0.3764705882352941, 0.4392156862745098, 0.5019607843137255, 0.5647058823529412, 
              0.6274509803921569, 0.6901960784313725, 0.7529411764705882, 0.8156862745098039, 0.8784313725490196, 0.9411764705882353, 1];
        } else {
            color_range = ['#d8d8d8','#395252','#000000'];
            color_domain = [0, 0.5, 1]
        }
       
        attribute_to_color = function(d) {
          if (d === undefined) {
            d = 0.5
          }
          if (high === low) {
            var x = 0.5
          } else {
             var x = (d - low) / (high-low);
          }
         
          return d3.scale.linear().domain(color_domain).range(color_range)(x);
        }
        
    } else {
        var labels = get_global_data('meta_levels')[plotted_item].sort();
        if (labels.length > 10) {
            attribute_to_color = function(d) {
                var i = labels.indexOf(d);
                var step = 8
                var L = Math.ceil(labels.length/step)*step
                var hue = Math.round(
                ( (i*L/step)%L + Math.floor(i/step) ) / L * 360 // Helps scatter values
                );
                return hsluv.hsluvToHex([hue, 90, 55])
             }
          } else {
            attribute_to_color = function(d) {
              var c = d3.scale.ordinal().domain(get_global_data('meta_levels')[plotted_item].sort()).range(plotly_defaults_colors)(d);
              return c
            }
      }
    } 
    
    var standard_label = tree.branch_name();
    
    function edgeStyler(element, node_data) {
      element.style ("stroke-width", 1)  
    }
    
    tree.style_edges(edgeStyler);
    
    function nodeStyler(element, node_data, is_radial) {
      if (d3.layout.phylotree.is_leafnode(node_data)) {
        var node_label = element.select("text");
        var level = plotted[node_data.name];
        
        var annotation = element.selectAll ("rect").data([level]);
        annotation.enter().append ("rect");
        if (is_radial) {
          annotation.attr ("height", 3).attr ("width", 30).attr ("x", -15).attr("y", -1.5)
            .style ("fill", function(d, i) {
              return attribute_to_color (d);
             }).style("stroke-width", 0);
           if (node_data.text_angle && is_radial) {
            annotation.attr("transform", function(d) {
              return "rotate(" + node_data.text_angle + ")";
            });
          } 
        } else {
            annotation.attr ("height", 4).attr ("width", 30).attr("y", -2)
              .style ("fill", function(d, i) {
                return attribute_to_color (d);
               }).style("stroke-width", 0);
            annotation.attr("transform", function(d) {
              return "rotate(0)";
              });
            }
      }
    }

    tree.style_edges(edgeStyler).style_nodes(function(element, node_data) { return nodeStyler(element, node_data, true)}).branch_name(function() {
        return ""
      });

    
    // sort nodes from phylotree example
    function sort_nodes (asc) {
      tree.traverse_and_compute (function (n) {
              var d = 1;
              if (n.children && n.children.length) {
                  d += d3.max (n.children, function (d) { return d["count_depth"];});
              }
              n["count_depth"] = d;
          });
          tree.resort_children (function (a,b) {
              return (a["count_depth"] - b["count_depth"]) * (asc ? 1 : -1);
          });
    }
    sort_nodes(false);
    
    var max_count_depth = Math.max(...tree.get_nodes().map(n => n["count_depth"]))
    // now change branch lengths
    tree.branch_length (function (n) {
      var parent_depth = n["parent"]["count_depth"];
      var length = 210 * (parent_depth - n["count_depth"]);
      if(n["count_depth"] == 1) {
        length = length + 120
      }
      return length
    });
    
    // try to collapse tree downwards
    
    
    
    tree.layout()
    tree.placenodes().update();
    
    var svgHeight = parseInt(svg.style("height").replace("px", ""));
    var svgWidth = parseInt(svg.style("width").replace("px", ""));
    
    
    var zoomFactor = 1000/(Math.min(svgHeight, svgWidth)*1.25);
    
    var horizontalOffset = (-1 * svgWidth / 2) * (1-zoomFactor);
    var verticalOffset = (-1 * svgHeight / 2) * (1-zoomFactor);
    
    
    svg.attr("transform", "translate(" +  horizontalOffset +"," + verticalOffset + ")" + "scale(" + zoomFactor + ")")
    
    
    
    $("#dend_layout").on("click", function(e) {
      var is_radial = $(this).prop("checked")
      tree.radial(is_radial)
      if (is_radial) {
        tree.size([s, s])
        tree.style_nodes(function(element, node_data) {return nodeStyler(element, node_data, true)}).placenodes().update();
         svg.attr("transform", "translate(" +  horizontalOffset +"," + verticalOffset + ")" + "scale(" + zoomFactor + ")")
      } else {
        tree.size([100, 770]);
        tree.style_nodes(function(element, node_data) {return nodeStyler(element, node_data, false)}).placenodes().update();
        svg.attr("transform", "translate(" +  0 +"," + 0 + ")") 
      }
      
    });
    
    //var collapseTo = 10;
    /* tree.traverse_and_compute(function(n){ 
      if(n["depth"] === collapseTo) {
        tree.collapse_node(n);
        console.log(n)
      }
    }); */
    
    
    this.tree = tree
    this.update_tree_selection()
    this.setLoadingStatus(false);
}


// START MODULES


function Modules_Table()
{
    this.dom_node = document.getElementById("modules-table-div-container");
    this.matrix = {}
    this.modules_table = $(this.dom_node).find('#modules-table');
    // this.clusters = {}
    // this.sorted_column = 'Score'
    // this.filterSig = $(this.dom_node).find('#sig_filt_input')
    // this.is_filtering = false
    // this.is_collapsed = {} // cluster -> boolean, holds whether a cluster is collapsed
    // this.tooltip = {} // will contain a Popper.js object
}


Modules_Table.prototype.init = function()
{

    var self = this;

    // var debounced = _.debounce(function(){self.doneTyping()}, 500)
    // self.filterSig.on('input', debounced)
    //
    //
    // var popper = document.querySelector('body > #sig-tooltip')
    // var arrow = popper.querySelector('.arrow')
    // var initialNode = document.querySelector('body > #hidden-node')
    //
    // var outPopper = new Popper(initialNode, popper, {placement: 'top',
    //     modifiers: {
    //         arrow: {
    //             enabled: true,
    //             element: arrow,
    //         },
    //         preventOverflow: {
    //             enabled: false,
    //         },
    //         hide: {
    //             enabled: false,
    //         },
    //     }
    // })
    //
    // this.tooltip = outPopper;
    //
    // // Set up the scrolling behavior
    // var header_div = this.dom_node.querySelector('.sig-tables-header')
    // var tables_div = this.dom_node.querySelector('.sig-tables-wrapper')
    // $(tables_div).on('scroll', function() {
    //     $(header_div).scrollLeft($(this).scrollLeft())
    // })
    //
    var update_promise = self.update({})
    return update_promise;

}

Modules_Table.prototype.render = function()
{
    var self = this;
    var matrix = self.matrix;
    var main_vis = get_global_status('main_vis');

    //var cluster_ids = _(self.clusters).values().uniq().sort(x => x).value()
    var table = this.modules_table.DataTable(
        {
            'columns': [
                {'title': 'Module'},
                {'title': 'C'},
                {'title': 'P'},
                {'title': 'FDR'},
            ],
            "pageLength": 10,
            "scrollY": '15vh',
            "order": [[1, "desc"]],
            "bFilter": false,
            "buttons": {
                dom: {
                    button: {
                        className: ''
                    }
                },
                buttons: [
                    {
                        extend: 'csvHtml5',
                        text: 'Export',
                        className: 'btn btn-secondary',
                    }
                ]
            },
            "pagingType": "simple",
        }
    );
    


    var dataSet = _.map(this.matrix, (f, i) => [f._row, f.C, f.pValue, f.FDR]);


    this.modules_table.DataTable().clear()
        .rows.add(dataSet)
         // Needed or else column headers don't align
        .draw()
        
    
    $('a[data-toggle="tab"]').on('shown.bs.tab', function(e){
      table.columns.adjust();
    });
    
    
    this.modules_table.on('click', "tr", function() {
        var cells = $(this).find("td")
        var module = cells.eq(0).text()
        
        set_global_status({
            'plotted_item_type': 'signature',
            'plotted_item': module,
        })
    });
    

    // // If the filtering is enabled, we're just going to put everything into one cluster
    // if (self.is_filtering){
    //     cluster_ids = [1]
    // }
    //
    // // Create the Header row
    // var header_row = $(self.dom_node).children(".sig-tables-header").find(".proj-row");
    // header_row.find("th:not(:first-child)").remove()
    // _.each(matrix.proj_labels, function(proj_label){
    //     var new_cell = $("<th>")
    //     var new_item = $("<div>")
    //     new_item.html("<div>" + proj_label + "</div>")
    //     new_cell.on("click", function() {
    //         var col_name = $(this).text()
    //         self.sorted_column = col_name
    //         self.render()
    //     });
    //     new_cell.append(new_item)
    //     header_row.append(new_cell)
    // });
    //
    // // Add padding for scrollbar width of table below it
    // $(self.dom_node).children(".sig-tables-header")
    //     .css("padding-right", detect_browser_scrollbar_width() + "px")
    //
    // var clusterTableDiv = $(self.dom_node).children('.sig-tables-wrapper').first()
    // var new_table_divs = []
    //
    // // Remove old tables
    // clusterTableDiv.children('.sig-table-div').remove()
    //
    // var curr_cl;
    // for (var i = 0; i < cluster_ids.length; i++) {
    //     curr_cl = cluster_ids[i];
    //
    //     if (!(curr_cl in self.is_collapsed)){
    //         self.is_collapsed[curr_cl] = true;
    //     }
    //
    //     // Create new table and add to table_div
    //     var new_table_div = document.createElement("div");
    //     new_table_div.setAttribute("class", "sig-table-div");
    //
    //     var new_table = document.createElement("table");
    //     new_table.setAttribute("id", "table"+ curr_cl);
    //     new_table.setAttribute("class", "sig-table");
    //     $(new_table).data("cluster", curr_cl)
    //
    //     var thead = document.createElement("thead");
    //
    //     var tbody = document.createElement("tbody");
    //
    //     new_table.appendChild(thead);
    //     new_table.appendChild(tbody);
    //
    //
    //     new_table_div.appendChild(new_table);
    //     new_table_divs.push(new_table_div);
    //
    //     if (typeof(matrix.sig_labels) == "string") {
    //         matrix.sig_labels = [matrix.sig_labels];
    //     }
    //
    //     // Format cell data for better d3 binding
    //     var sig_labels = matrix.sig_labels.filter(
    //         function(x) { return self.is_filtering || self.clusters[x] == curr_cl; }
    //     );
    //
    //     var zscores = []; // rows of z-scores
    //     var pvals = []; // rows of p-values
    //     var sig_indices = [] // signature index corresponding to each row
    //
    //     for (var ind = 0; ind < sig_labels.length; ind++) {
    //         var sig = sig_labels[ind];
    //         var sig_index = matrix.sig_labels.indexOf(sig)
    //
    //         pvals.push(matrix.pvals[sig_index]);
    //         zscores.push(matrix.zscores[sig_index]);
    //         sig_indices.push(sig_index);
    //     }
    //
    //     var formatted_data_matrix = _.zip(sig_indices, zscores, pvals)
    //         .map(x => {
    //             var sig_index = x[0]
    //             var zscores_row = x[1]
    //             var pvals_row = x[2]
    //             return _.zip(zscores_row, pvals_row).map((x, j)  => {
    //                 return {
    //                     "row": sig_index, "col": j,
    //                     "zscore": x[0], "pval": x[1]
    //                 }
    //             })
    //         })
    //
    //     var formatted_data_w_row_labels = d3.zip(sig_labels, formatted_data_matrix);
    //
    //     // Sort data if necessary
    //
    //     var sort_col = matrix.proj_labels.indexOf(self.sorted_column);
    //     if(sort_col > -1){
    //         var sortFun = function(a,b){
    //             return b[1][sort_col].zscore - a[1][sort_col].zscore; // Descending order
    //         };
    //         formatted_data_w_row_labels.sort(sortFun);
    //     }
    //
    //     $(new_table_div).data('table-sort-val', 0)
    //     if (cluster_ids.length > 1) {
    //         // Pull out the cluster leader and put it first
    //
    //         var row_vals;
    //         if (main_vis === "pcannotator") {
    //             row_vals = _.map(formatted_data_w_row_labels, x => {
    //                 return _(x[1]).map(y => y.pval).sum()
    //             })
    //         } else {
    //             row_vals = _.map(formatted_data_w_row_labels, x => {
    //                 return (x[1][0].zscore*-1) // sort by zscore (geary-c actually)
    //             })
    //         }
    //
    //         var leader_row_i = row_vals.indexOf(_.min(row_vals))
    //         var leader_row = formatted_data_w_row_labels[leader_row_i];
    //         formatted_data_w_row_labels.splice(leader_row_i, 1);
    //         formatted_data_w_row_labels.splice(0, 0, leader_row);
    //
    //         if (sort_col > -1) {
    //             $(new_table_div).data('table-sort-val', leader_row[1][sort_col].zscore*-1)
    //         }
    //     }
    //
    //     if (main_vis == "pcannotator") {
    //         var colorScale = d3.scale.linear()
    //             .domain([-.5,0,.5])
    //             .range(["steelblue", "white", "lightcoral"])
    //             .clamp(true);
    //     } else {
    //         var colorScale = d3.scale.linear()
    //             .domain([0,.5,1])
    //             .range(["steelblue","white", "lightcoral"])
    //             .clamp(true);
    //     }
    //     var colorScaleCluster = d3.scale.linear()
    //         .domain([0,0.5,1])
    //         .range(["steelblue","white", "lightcoral"])
    //         .clamp(true);
    //
    //     var content_rows = d3.select(new_table).select('tbody').selectAll('tr')
    //         .data(formatted_data_w_row_labels);
    //     content_rows.enter().append('tr');
    //     content_rows.exit().remove();
    //
    //     var content_row = content_rows.selectAll("td")
    //         .data(function(d){return [d[0]].concat(d[1]);})
    //
    //     content_row.enter().append('td');
    //     content_row.exit().remove();
    //
    //     if (main_vis === 'pcannotator') {
    //         content_row
    //             .filter(function(d,i) { return i > 0;})
    //             .style('background-color', function(d){return colorScale(d.zscore);})
    //             .on("click", function(d){tableClickFunction_PC(matrix.sig_labels[d.row], 'signature')})
    //             .append('div');
    //
    //     } else {
    //         content_row
    //             .filter(function(d,i) { return i > 0;})
    //             .on("click", function(d){tableClickFunction_clusters(matrix.sig_labels[d.row], 'signature')})
    //
    //         content_row
    //             .filter(function(d,i) { return i == 1;})
    //             .style('background-color', function(d){return colorScale(d.zscore);})
    //             .append('div')
    //             .text(function(d){ return d.zscore.toPrecision ? d.zscore.toPrecision(2) : d.zscore;})
    //
    //         content_row
    //             .filter(function(d,i) { return i > 1;})
    //             .style('background-color', function(d){return colorScaleCluster(d.zscore);})
    //             .append('div')
    //             .text(function(d){ return d.pval < 0.05 ? '*' : '';})
    //     }
    //
    //     // Hover actions
    //     var rowHoverFunc = function(header_row){
    //         return function(d, i){
    //             var tooltip_str;
    //             if(main_vis === 'pcannotator'){
    //                 tooltip_str = "corr = " + d.zscore.toFixed(2)
    //             } else {
    //                 if(main_vis === 'clusters' && i > 0)
    //                     tooltip_str = "AUC=" + _auc_format(d.zscore) + " p<" + _pval_format(d.pval)
    //                 else
    //                     tooltip_str = "C'=" + d.zscore.toFixed(2) + ", p<" + _pval_format(d.pval)
    //             }
    //             createTooltip(self.tooltip, this, tooltip_str)
    //             hoverRowCol(header_row, this, matrix.proj_labels[d.col])
    //         }
    //     }(header_row, new_table);
    //     var rowUnHoverFunc = function(header_row){
    //         return function(d){
    //             destroyTooltip(self.tooltip)
    //             unhoverRowCol(header_row, this, matrix.proj_labels[d.col])
    //         }
    //     }(header_row, new_table);
    //
    //     content_row
    //         .filter(function(d,i) { return i > 0;})
    //         .on("mouseenter", rowHoverFunc)
    //         .on("mouseleave", rowUnHoverFunc)
    //
    //     // Add text for signature names
    //     content_row.filter(function(d,i) { return i == 0;})
    //         .html(function(d){return '<div class="sig-row-label"><div>'+d+'</div></div>';})
    //
    //     // Add cluster expand/collapse behavior ONLY if more than one cluster
    //     if (cluster_ids.length > 1) {
    //
    //         $(new_table).addClass('collapsible')
    //
    //         if (self.is_collapsed[curr_cl]){
    //             $(new_table).addClass('collapsed')
    //         }
    //
    //         $(new_table).children('tbody').children('tr:first-child').children("td:first-child")
    //             .on("click", function() { self.clickSummaryRow(this); });
    //     }
    //
    // }
    //
    // new_table_divs = _.sortBy(new_table_divs, d => {
    //     return $(d).data('table-sort-val')
    // });
    //
    // _.forEach(new_table_divs, new_table_div => {
    //     clusterTableDiv.append(new_table_div);
    // });
    //
    // // Apply compact styling for pcannotator
    // if(main_vis === "pcannotator" || main_vis === "clusters" || main_vis === "tree"){
    //     $(self.dom_node).find('table').addClass('compact')
    // } else {
    //     $(self.dom_node).find('table').removeClass('compact')
    // }
    // if(main_vis === "clusters" || main_vis === "tree"){
    //     $(self.dom_node).find('table').addClass('first-col')
    // } else {
    //     $(self.dom_node).find('table').removeClass('first-col')
    // }
    //
    // // Trigger existinging filter
    // self.filterSig.trigger('input');
}

// Modules_Table.prototype.clickSummaryRow = function(d){
//
//     var table_id = $(d).closest('table')
//     var cluster_id = table_id.data('cluster')
//
//     if (this.is_collapsed[cluster_id]){
//         this.is_collapsed[cluster_id] = false;
//         table_id.removeClass('collapsed')
//     } else {
//         this.is_collapsed[cluster_id] = true;
//         table_id.addClass('collapsed')
//     }
// }


Modules_Table.prototype.update = function(updates)
{
    var self = this;

    // If there are no modules, don't update
    var has_mods = get_global_status('has_mods')
    if(!has_mods){
        return $.when(null); // Empty promise
    }
    

    var lcs_promise;
    if (_.isEmpty(self.matrix)){

        lcs_promise = api.modules.lcs().then(function(mat){
                self.matrix = mat;
                return true;
            })
    } else {
        lcs_promise = false
    }
    /**
    var matrix_promise;
    if( 'main_vis' in updates || _.isEmpty(self.matrix) || 'cluster_var' in updates){
        var main_vis = get_global_status('main_vis');
        var cluster_var = _get_cluster_var();

        var fix_col_label = true
        if (main_vis == "clusters") {
            matrix_promise = api.clusters.sigProjMatrix(cluster_var, false); // YANAY
        } else if (main_vis == "tree") {
            matrix_promise = api.tree.sigProjMatrix(false);
        } else {
            matrix_promise = api.pCorr.signatures();
            fix_col_label = false
        }

        matrix_promise = matrix_promise
            .then(function(matrix) {
                if (fix_col_label){
                    matrix.proj_labels[0] = "Score";
                }
                self.matrix = matrix;
                return true;
            });
    } else {
        matrix_promise = false
    }
    **/
    return $.when(lcs_promise)
        .then(function(lcs_update){
            if(lcs_update){
                self.render();
            }
        });

}

// END MODULES


// Function that's triggered when clicking on table cell
function tableClickFunction_PC(row_key, item_type)
{
    var update = {}

    if (item_type == "signature") {
        // Could also be meta-data because we mix those in
        var meta_sigs = get_global_data('meta_sigs')
        if (meta_sigs.indexOf(row_key) > -1) {
            item_type = 'meta'
        }
    }

    update['plotted_item_type'] = item_type;
    update['plotted_item'] = row_key;

    set_global_status(update);
}

function tableClickFunction_clusters(row_key, item_type)
{
    var update = {}

    if (item_type == "signature") {
        // Could also be meta-data because we mix those in
        var meta_sigs = get_global_data('meta_sigs')
        if (meta_sigs.indexOf(row_key) > -1) {
            item_type = 'meta'
        }
    }

    update['plotted_item_type'] = item_type;
    update['plotted_item'] = row_key;

    set_global_status(update);
}


function hoverRowCol(header_row, node, col){
    $(header_row).find('th').filter((i, e) => $(e).text() == col).addClass('highlight')
    $(node).siblings('td:first-child').addClass('highlight')
}

function unhoverRowCol(header_row, node, col){
    $(header_row).find('th').filter((i, e) => $(e).text() == col).removeClass('highlight')
    $(node).siblings('td:first-child').removeClass('highlight')
}

function createTooltip(popper, node, text) {
    var popper_node = popper.popper
    var inner_node = popper_node.querySelector('.inner')
    inner_node.textContent = text

    popper.reference = node;
    popper.scheduleUpdate();
}

function destroyTooltip(popper) {
    var node = document.querySelector('body #hidden-node')
    popper.reference = node;
    popper.scheduleUpdate();
}

function _pval_format(p) {
    if(p < 0.01){
        return p.toExponential(1)
    } else {
        return p.toFixed(2)
    }
}

function _auc_format(x) {
    // x is on range 0 to 1
    // Convert to positive 0.5-1 AUC
    var auc = Math.max(x, 1-x)
    return auc.toFixed(2)
}
