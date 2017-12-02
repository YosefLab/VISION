/*
 * This script contains several components
 *
 * Upper_Left_Content
 *    - Signature_Table
 *    - Precomputed_Table
 *    - Gene_Select
 *
 */

function Upper_Left_Content()
{
    this.children = []
    this.sig_table = {}
    this.pc_table = {}
    this.gene_select = {}
}

Upper_Left_Content.prototype.init = function()
{
    var sig_table = new Signature_Table()
    this.children.push(sig_table)
    this.sig_table = sig_table

    var sig_table_promise = sig_table.init();

    var pc_table = new Precomputed_Table()
    this.children.push(pc_table)
    this.pc_table = pc_table

    var pc_table_promise = pc_table.init();

    var gene_select = new Gene_Select()
    this.children.push(gene_select)
    this.gene_select = gene_select

    var gene_select_promise = gene_select.init()

    return $.when(sig_table_promise, pc_table_promise, gene_select_promise);

}

Upper_Left_Content.prototype.update = function(updates)
{
    // Updates passed to children components
    _.each(this.children, function(child){
        child.update(updates)
    });

}

Upper_Left_Content.prototype.select_default = function()
{
    var main_vis = get_global_status('main_vis');
    var new_projection = 'PCA: 1,2'

    var update = {}
    update['plotted_projection'] = new_projection
    if(main_vis === 'sigvp'){
        update['plotted_item'] = this.sig_table.get_top_sig(new_projection)
    }
    update['plotted_item_type'] = 'signature'

    set_global_status(update)
}


function Signature_Table()
{
    this.dom_node = document.getElementById("table-div-container");
    this.matrix = {}
    this.clusters = {}
    this.sorted_column = 'PCA: 1,2'
    this.filterSig = $(this.dom_node).find('#sig_filt_input')
}


Signature_Table.prototype.init = function()
{

    var self = this;

    var debounced = _.debounce(function(){self.doneTyping()}, 500)
    self.filterSig.on('input', debounced)

    var update_promise = self.update({})

    return update_promise;

}

Signature_Table.prototype.render = function()
{
    var self = this;
    var matrix = self.matrix;
    var main_vis = get_global_status('main_vis');
    var clusarr = Object.keys( self.clusters ).map(
        function ( key ) {return self.clusters[key];}
    );

    var clusmax = Math.max.apply(null, clusarr);

    // Create the Header row
    var header_row = $(self.dom_node).children("#header-div").find("#proj_row");
    header_row.find("th:not(:first-child)").remove()
    _.each(matrix.proj_labels, function(proj_label){
        var new_cell = $("<th>")
        var new_item = $("<div>")
        new_item.html(proj_label)
        new_item.on("click", function() {
            var col_name = $(this).html()
            self.sorted_column = col_name
            self.render()
        });
        new_cell.append(new_item)
        header_row.append(new_cell)
    });

    var clusterTableDiv = $(self.dom_node).children('#cluster-tables-div').first()

    // Remove old tables
    clusterTableDiv.children('.sig-table-div').remove()

    for (var curr_cl = 1; curr_cl <= clusmax; curr_cl++) {
        // Create new table and add to table_div
        var new_table_div = document.createElement("div");
        new_table_div.setAttribute("style", "height=calc((100vh - 88px) / 2)");
        new_table_div.setAttribute("style", "overflow: hidden");
        new_table_div.setAttribute("class", "sig-table-div");

        var new_table = document.createElement("table");
        new_table.setAttribute("id", "table"+ curr_cl);
        new_table.setAttribute("class", "sig-table");

        var thead = document.createElement("thead");

        var tbody = document.createElement("tbody");

        new_table.appendChild(thead);
        new_table.appendChild(tbody);


        new_table_div.appendChild(new_table);
        clusterTableDiv.append(new_table_div);
    }


    for (var curr_cl = 1; curr_cl <= clusmax; curr_cl++) {

        if (typeof(matrix.sig_labels) == "string") {
            matrix.sig_labels = [matrix.sig_labels];
        }

        // Format cell data for better d3 binding
        var sig_labels = matrix.sig_labels.filter(
            function(x) { return self.clusters[x] == curr_cl; }
        );

        var data = [];
        for (var ind = 0; ind < sig_labels.length; ind++) {
            var sig = sig_labels[ind];
            data.push(matrix.data[matrix.sig_labels.indexOf(sig)]);
        }

        if (typeof(sig_labels) == "string") {
            sig_labels = [sig_labels];
        }

        var formatted_data_matrix = sig_labels.map(function(row, i){
            return data[i].map(function(val, j){
                return {"val":val, "row":matrix.sig_labels.indexOf(sig_labels[i]), "col":j}
            });
        });


        var formatted_data_w_row_labels = d3.zip(sig_labels, formatted_data_matrix);

        // Sort data if necessary

        var sort_col = matrix.proj_labels.indexOf(self.sorted_column);
        if(sort_col > -1){
            var sortFun = function(a,b){
                var a_precomp = global_data.sigIsPrecomputed[a[0]];
                var b_precomp = global_data.sigIsPrecomputed[b[0]];
                if(a_precomp && b_precomp || !a_precomp && !b_precomp){
                    if (main_vis == "sigvp") {
                        return a[1][sort_col].val - b[1][sort_col].val;
                    } else {
                        return b[1][sort_col].val - a[1][sort_col].val;
                    }	
                }
                else if (a_precomp) { return -1;}
                else {return 1;}
            };
            formatted_data_w_row_labels.sort(sortFun);
        }

        if (main_vis == "pcannotator") {
            var colorScale = d3.scale.linear()
                .domain([0,0.4,.8])
                .range(["steelblue", "white", "lightcoral"])
                .clamp(true);
        } else {
            var colorScale = d3.scale.linear()
                .domain([0,-3,-50])
                .range(["steelblue","white", "lightcoral"])
                .clamp(true);
        }

        var content_rows = d3.select('#table'+curr_cl).select('tbody').selectAll('tr')
            .data(formatted_data_w_row_labels);
        content_rows.enter().append('tr');
        content_rows.exit().remove();

        var content_row = content_rows.selectAll("td")
            .data(function(d){return [d[0]].concat(d[1]);})

        content_row.enter().append('td');
        content_row.exit().remove();

        if (main_vis == "pcannotator") {
            content_row
                .filter(function(d,i) { return i > 0;})
                .text(function(d){
                    if(d.val > .8) { return "> .8";}
                    else if(d.val < .8) { return d.val.toFixed(2);}
                    else {return d.val.toPrecision(2);}
                })
                .style('background-color', function(d){return colorScale(d.val);})
                .on("click", function(d){tableClickFunction_PC(matrix.sig_labels[d.row], matrix.proj_labels[d.col], 'signature')});

        } else {
            content_row
                .filter(function(d,i) { return i > 0;})
                .text(function(d){
                    if(d.val < -50) { return "< -50";}
                    else if(d.val > -1) { return d.val.toFixed(2);}
                    else {return d.val.toPrecision(2);}
                })
                .style('background-color', function(d){return colorScale(d.val);})
                .on("click", function(d){tableClickFunction(matrix.sig_labels[d.row], matrix.proj_labels[d.col], 'signature')});
        }

        // Make signature names click-able
        content_row.filter(function(d,i) { return i == 0;})
            .text(function(d){return d;})

        $('#table'+curr_cl).children('tbody').children("tr:not(:first-child)").addClass('collapsed')
        $("#table"+curr_cl).children("tbody").children("tr:first-child").children("td:first-child").attr("id", "sigclust_" + curr_cl);
        $("#table"+curr_cl).children('tbody').children('tr:first-child').children("td:first-child")
            .on("click", function() { clickSummaryRow(this); });

        // Add '>' sign to top sig name to indiciate expandability
        $("#table"+curr_cl).children("tbody").children("tr:first-child").children("td:first-child")
            .text(function(i, origText) {
                return origText + " \u25B6 ";
            });

    }


    $(".sigclust").on("mouseover", function(d) {
        tooltip.showTooltip("Click To Toggle Cluster Display", d);
    })
        .on("mouseout", function() {
            tooltip.hideTooltip();
        });

    self.filterSig.trigger('input'); 
}

Signature_Table.prototype.update = function(updates)
{

    var clusters_promise;
    var self = this;
    if (_.isEmpty(self.clusters)){

        clusters_promise = api.signature.clusters(false, "1")
            .then(function(cls){
                self.clusters = cls;
                return true;
            })
    } else {
        clusters_promise = false
    }

    var matrix_promise;
    if( 'main_vis' in updates || 'filter_group' in updates || _.isEmpty(self.matrix)){
        var self = this;
        var main_vis = get_global_status('main_vis');
        var filter_group = get_global_status('filter_group');

        if (main_vis == "sigvp") {
            matrix_promise = api.filterGroup.sigProjMatrixP(filter_group, false);
        } else if (main_vis == "tree") {
            matrix_promise = api.filterGroup.treeSigProjMatrixP(filter_group, false);
        } else {
            matrix_promise = api.filterGroup.pCorr(filter_group, false);
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

Signature_Table.prototype.get_top_sig = function(projection)
{
    var matrix = this.matrix;
    var j = matrix.proj_labels.indexOf(projection);
    var s_i = matrix.data.map(function(e){return e[j];}).argSort();

    return matrix.sig_labels[s_i[0]]
}


function Precomputed_Table()
{
    this.dom_node = document.getElementById("precomp-table-div");
    this.matrix = {}
    this.sorted_column = 'PCA: 1,2'
}

Precomputed_Table.prototype.init = function()
{
    var self = this;

    var update_promise = self.update({})

    return update_promise;
}

Precomputed_Table.prototype.render = function()
{
    var self = this;
    var matrix = self.matrix;
    var main_vis = get_global_status('main_vis');

    var header_row = d3.select(self.dom_node).select('#precomp-table').select("thead").select("#proj_row").selectAll("th")
        .data([""].concat(matrix.proj_labels));

    header_row.enter().append('th');
    header_row.html(function(d){return "<div>" + d+"</div>";})
        .filter(function(d,i) {return i > 0;})
        .on("click", function(col_name) { 
            self.sorted_column = col_name;
            self.render();
        });
    header_row.exit().remove();


    if (typeof(matrix.sig_labels) == "string") {
        matrix.sig_labels = [matrix.sig_labels];
    }

    // Format cell data for better d3 binding
    var sig_labels = matrix.sig_labels;
    var data = [];
    for (var ind = 0; ind < sig_labels.length; ind++) {
        var sig = sig_labels[ind];
        data.push(matrix.data[matrix.sig_labels.indexOf(sig)]);
    }

    if (typeof(sig_labels) == "string") {
        sig_labels = [sig_labels];
    }

    var formatted_data_matrix = sig_labels.map(function(row, i){
        return data[i].map(function(val, j){
            return {"val":val, "row":matrix.sig_labels.indexOf(sig_labels[i]), "col":j}
        });
    });


    var formatted_data_w_row_labels = d3.zip(sig_labels, formatted_data_matrix);

    // Sort data if necessary

    var sort_col = matrix.proj_labels.indexOf(self.sorted_column);
    if(sort_col > -1){
        var sortFun = function(a,b){
            var a_precomp = global_data.sigIsPrecomputed[a[0]];
            var b_precomp = global_data.sigIsPrecomputed[b[0]];
            if(a_precomp && b_precomp || !a_precomp && !b_precomp){
                if (main_vis == "sigvp") {
                    return a[1][sort_col].val - b[1][sort_col].val;
                } else {
                    return b[1][sort_col].val - a[1][sort_col].val;
                }	
            }
            else if (a_precomp) { return -1;}
            else {return 1;}
        };
        formatted_data_w_row_labels.sort(sortFun);
    }

    if (main_vis == "pcannotator") {
        var colorScale = d3.scale.linear()
            .domain([0,0.2,0.5])
            .range(["steelblue", "white", "lightcoral"])
            .clamp(true);
    } else {
        var colorScale = d3.scale.linear()
            .domain([0,-3,-50])
            .range(["steelblue","white", "lightcoral"])
            .clamp(true);
    }


    var content_rows = d3.select(self.dom_node).select('#precomp-table').select('tbody').selectAll('tr')
        .data(formatted_data_w_row_labels);
    content_rows.enter().append('tr');
    content_rows.exit().remove();

    var content_row = content_rows.selectAll("td")
        .data(function(d){return [d[0]].concat(d[1]);})

    content_row.enter().append('td');
    content_row.exit().remove();

    if (main_vis == "pcannotator") {
        content_row
            .filter(function(d,i) { return i > 0;})
            .text(function(d){
                if(d.val > .5) { return "> .5";}
                else if(d.val < .5) { return d.val.toFixed(2);}
                else {return d.val.toPrecision(2);}
            })
            .style('background-color', function(d){return colorScale(d.val);})
            .on("click", function(d){tableClickFunction_PC(matrix.sig_labels[d.row], matrix.proj_labels[d.col], 'meta')});

    } else {
        content_row
            .filter(function(d,i) { return i > 0;})
            .text(function(d){
                if(d.val < -50) { return "< -50";}
                else if(d.val > -1) { return d.val.toFixed(2);}
                else {return d.val.toPrecision(2);}
            })
            .style('background-color', function(d){return colorScale(d.val);})
            .on("click", function(d){tableClickFunction(matrix.sig_labels[d.row], matrix.proj_labels[d.col], 'meta')});
    }

    // Make signature names 
    content_row.filter(function(d,i) { return i == 0;})
        .text(function(d){return d;});

}

Precomputed_Table.prototype.update = function(updates)
{
    var self = this;
    var matrix_promise;

    if('main_vis' in updates || 'filter_group' in updates || _.isEmpty(self.matrix)){
        var filter_group = get_global_status('filter_group');
        var main_vis = get_global_status('main_vis');

        if (main_vis === "sigvp") {
            matrix_promise = api.filterGroup.sigProjMatrixP(filter_group, true);
        } else if (main_vis === "tree") {
            matrix_promise = api.filterGroup.treeSigProjMatrixP(filter_group, true);
        } else {
            matrix_promise = api.filterGroup.pCorr(filter_group, true);
        }

        matrix_promise = matrix_promise.then(
            function(matrix){
                self.matrix = matrix
                return true
            });
    } else {
        matrix_promise = $.when(false);
    }

    return matrix_promise.then(function(matrix_updated) {
        if(matrix_updated){
            self.render()
        }
    });

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


    // Get a list of the projection names
    var proj_names_promise = api.filterGroup.listProjections("fano")
        .then(function(proj_names) {

            var projSelect = $('#SelectProj');

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
                .on('change', function () {
                    set_global_status({
                        'plotted_projection':$(this).val(),
                    });
                });

        });

    this.render_recent_genes()

    return $.when(gene_promise, proj_names_promise);

}

Gene_Select.prototype.update = function(updates)
{
    if('plotted_projection' in updates)
    {
        var plotted_projection = get_global_status('plotted_projection')
        var projSelect = $('#SelectProj');
        projSelect.val(plotted_projection).trigger('chosen:updated')
    }

    // Update the 'recent-genes' list
    if('plotted_item' in updates){
        var gene = get_global_status('plotted_item')
        var plotted_item_type = get_global_status('plotted_item_type')
        if(plotted_item_type === 'gene'){
            this.recent_genes = this.recent_genes.filter(function(e) {
                return e !== gene;
            })
            this.recent_genes.unshift(gene)

            if(this.recent_genes.length > 10){
                this.recent_genes.pop();
            }

            this.render_recent_genes()
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

    var vals = val.split(",");
    vals = vals.map(function(str){return str.trim();})
        .filter(function(str){ return str.length > 0;});

    var tablerows = $(self.dom_node).find('table').find('tr:not(:first-child)');
    tablerows.removeClass('filtered');

    var posvals = vals.filter(function(str){ return str[0] != '!';});
    var negvals = vals.filter(function(str){ return str[0] == '!';})
        .map(function(str){ return str.slice(1);})
        .filter( function(str){return str.length > 0;});

    if(posvals.length > 0){
        tablerows.filter(function(i, element){
            if(i == 0){return false;} // Don't filter out the header row
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
            if(i == 0){return false;} // Don't filter out the header row
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

    /*tablerows.removeClass('altRow')
            .not('.filtered').filter(':odd').addClass('altRow');
            */
}


function clickSummaryRow(d) {
    var clust = d.id.split("_")[1];

    var table_id = $("#table"+clust).children("tbody");
    if (!table_id.children("tr:not(:first-child)").hasClass('collapsed')) {
        table_id.children("tr:not(:first-child)").addClass('collapsed');

        table_id.children("tr:first-child").children("td:first-child")
            .text(function(i, origText) {
                origText = origText.split(" ")[0];
                return origText + " \u25B6 ";
            });
    } else {
        table_id.children("tr:not(:first-child)").removeClass('collapsed');

        table_id.children("tr:first-child").children("td:first-child")
            .text(function(i, origText) {
                origText = origText.split(" ")[0];
                return origText + " \u25BC";
            });
    }
}

// Function that's triggered when clicking on table cell
function tableClickFunction(row_key, col_key, item_type)
{

    var update = {}
    update['plotted_item_type'] = item_type;
    update['plotted_item'] = row_key;
    update['plotted_projection'] = col_key;

    set_global_status(update);
}

function tableClickFunction_PC(row_key, col_key, item_type)
{
    var update = {}
    update['plotted_item_type'] = item_type;
    update['plotted_item'] = row_key;
    update['plotted_pc'] = col_key;

    set_global_status(update);
}
