/*
 * This script contains several components
 *
 * Lower_Left_Content
 *    - Sig_Info
 *      - Sig_Heatmap
 *      - Cell_SigInfo
 *    - Gene_Info
 *    	- 
 *    - Meta_Info
 *    - Cell_Info
 *
 */

function Lower_Left_Content()
{
    this.children = []
    this.sig_info = {}
    this.gene_info = {}
    this.meta_info = {}
}

Lower_Left_Content.prototype.init = function()
{
    var sig_info = new Sig_Info()
    this.children.push(sig_info)
    this.sig_info = sig_info

    var sig_info_promise = sig_info.init();

    var gene_info = new Gene_Info()
    this.children.push(gene_info)
    this.gene_info = gene_info

    var gene_info_promise = gene_info.init();

    var meta_info = new Meta_Info()
    this.children.push(meta_info)
    this.meta_info = meta_info

    var meta_info_promise = meta_info.init();

    //var cell_info = new Cell_Info();
    //this.children.push(cell_info);
    //this.cell_info = cell_info

    //var cell_info_promise = cell_info.init();

    this.setLoadingStatus = createLoadingFunction(
        document.getElementById("lower-left-content")
    );

    return $.when(sig_info_promise, gene_info_promise, meta_info_promise);
}

Lower_Left_Content.prototype.update = function(updates)
{
    // Updates passed to children components
    _.each(this.children, function(child){
        child.update(updates)
    });

    if('plotted_item_type' in updates){
        var item_type = get_global_status('plotted_item_type')
        if(item_type === 'gene'){
            $(this.sig_info.dom_node).hide()
            $(this.meta_info.dom_node).hide()
            $(this.gene_info.dom_node).show()
        } else if (item_type === 'signature') {
            $(this.gene_info.dom_node).hide()
            $(this.meta_info.dom_node).hide()
            $(this.sig_info.dom_node).show()
        } else if (item_type === 'meta') {
            $(this.sig_info.dom_node).hide()
            $(this.gene_info.dom_node).hide()
            $(this.meta_info.dom_node).show()
        } else if (item_type == "cell") {
	    $(this.sig_info.dom_node).hide()
	    $(this.gene_info.dom_node).hide()
	    $(this.meta_info.dom_node).hide()
	}
    }

}

Lower_Left_Content.prototype.hover_cells = function()
{
}

function Sig_Info()
{
    this.dom_node = document.getElementById("sig-info");
    this.title = $(this.dom_node).find('#data-analysis-title');
    this.source = $(this.dom_node).find('#sig-source');
    this.content = $(this.dom_node).find('#sig-table-wrapper');
    this.bound_sig = ""
    this.sig_info_default = $(this.dom_node).find('#sig-info-default')
    this.sig_info_cluster = $(this.dom_node).find('#sig-info-cluster')
    this.sig_info_cell = $(this.dom_node).find("#sig-info-cell")

    this.heatmap = null
    this.plotted_heatmap = {
        'sig_key': '',
        'proj_key': '',
    }

    this.cell_info = null
    this.plotted_cellinfo = {
        'sig_key': '',
        'cells': [],
    }

    this.CellInfo = null;
}

Sig_Info.prototype.init = function()
{
    var self = this;
    var dt = $(this.content).find('#sig-info-table')
    dt.DataTable( {
        columns: [
            {
                'title': 'Gene',
                'render': function(data)
                {
                    return "<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=" + data + " target='_blank'>" + data + "</a>";
                }
            }, 
            {'title': 'Sign', 'className': 'dt-center'}],
        'paging': false,
        'info': true,
        'scrollY': '15vh',
        'scrollCollapse': true,
    })

    $(self.dom_node).find('#HeatmapButton').on('click', function()
    {
        self.drawHeat()
    })

    $(self.dom_node).find("#CellInfoButton").on('click', function()
    {
	self.addCellInfo()
    });

}

Sig_Info.prototype.update = function(updates)
{

    var sig_info = get_global_data('sig_info');
    if(sig_info.name === this.bound_sig || sig_info.isMeta)
    {
        return;
    }

    $(self.dom_node).find('#SigInfoButton').click();

    this.bound_sig = sig_info.name

    $(this.title).text(sig_info.name)

    // cut down the sig-info source
    var source = sig_info.source

    if(source.lastIndexOf('/') !== -1){
        source = source.slice(source.lastIndexOf('/')+1)
    }

    if(source.length > 30){
        source = source.slice(0, 27) + '...';
    }


    $(this.source).text(source);

    // Render the gene table

    var dt = $(this.content).find('#sig-info-table')

    var dataSet = _.map(sig_info.sigDict, function (value, key){
        if(value > 0){
            return [key, '+'];
        } else {
            return [key, '-'];
        }
    })

    dt.DataTable().clear()
        .rows.add(dataSet)
        .draw()

    if(this.sig_info_cluster.hasClass('active') && 'plotted_item' in updates){
        this.drawHeat();
    }

    if (this.sig_info_cell.hasClass("active")) {
	this.addCellInfo();
    }	
}

Sig_Info.prototype.build_cluster_dropdown_param = function()
{
    // Rebuild the 'param' if the first dropdown is changed
    var self = this;
    var vals = self.cluster_options[$(self.dom_node).find('#cluster_select_method').val()]
    var clust_dropdown_param = $(self.dom_node).find('#cluster_select_param');
    var old_val = clust_dropdown_param.val()

    clust_dropdown_param.empty();
    for(var i=0; i<vals.length; i++){
        clust_dropdown_param.append($("<option />").val(vals[i]).text(vals[i]));
    }

    if(vals.indexOf(old_val) > -1){
        clust_dropdown_param[0].selectedIndex = vals.indexOf(old_val)
    }
}

// Draw the heatmap
Sig_Info.prototype.drawHeat = function(){

    // Need to create it if it isn't there
    var self = this;

    var heatmap_div = $(self.dom_node).find('#heatmap-div')

    if( heatmap_div.children('svg').length === 0)
    {
        var heatmap_width = heatmap_div.parent().parent().width();
        var heatmap_height = heatmap_div.parent().parent().height()-40;

        self.heatmap = new HeatMap('#heatmap-div', heatmap_width, heatmap_height);
    }

    var sig_key = get_global_status('plotted_item'); // assume it's a signature

    var sig_info = get_global_data('sig_info');

    // plotted heatmap is based on sig_key
    // check if we are already showing the right heatmap and don't regenerate
    var need_plot = false;

    if (self.plotted_heatmap['sig_key'] !== sig_key){
        self.plotted_heatmap['sig_key'] = sig_key
        need_plot = true;
    }

    if(!need_plot){
        return $.when(true);
    }

    heatmap_div.addClass('loading')

    return $.when(
        api.signature.expression(sig_key))
        .then(function(sig_expression){


            //Construct data matrix
            // TODO: sort genes

            var dataMat = sig_expression.data;
            var gene_labels = sig_expression.gene_labels;
            var sample_labels = sig_expression.sample_labels;
            var clusters = get_global_data('clusters')

            var gene_signs = gene_labels.map(function(e){
                return sig_info.sigDict[e]
            });

            var assignments = sample_labels.map(sample => clusters[sample]);

            self.heatmap.setData(dataMat,
                assignments,
                gene_labels,
                gene_signs,
                sample_labels);

        }).always(function() {
            heatmap_div.removeClass('loading');
        });
}

// Add Cell info for signature data
Sig_Info.prototype.addCellInfo = function(){

    // Need to create it if it isn't there
    var self = this;

    var cellinfo_div = $(self.dom_node).find('#cell-dist-div')
    cellinfo_div.show()

    if (self.CellInfo == null) { 
        self.CellInfo = new Cell_Info("signature");
    } 

    var sig_key = get_global_status('plotted_item'); // assume it's a signature
    var cells = get_global_status("selected_cells");

    // plotted heatmap is based on sig_key
    // check if we are already showing the right heatmap and don't regenerate
    var need_plot = false;

    if (self.plotted_cellinfo['sig_key'] !== sig_key){
        self.plotted_cellinfo['sig_key'] = sig_key
        need_plot = true;
    }

    if (self.plotted_cellinfo['cells'] !== cells) { 
        self.plotted_cellinfo['cells'] = cells;
        need_plot = true;
    }

    if(!need_plot){
        return $.when(true);
    }

    cellinfo_div.addClass('loading')

    self.CellInfo.update()

    cellinfo_div.removeClass("loading");

}

function Gene_Info()
{
    this.dom_node = document.getElementById("gene-info");
    this.title = $(this.dom_node).find('#gene-analysis-title').get(0);
    this.chart = $(this.dom_node).find('#gene-dist-div').get(0);
    this.gene_info_default = $(this.dom_node).find('#gene-info-default');
    this.gene_info_cell = $(this.dom_node).find("#gene-info-cell");
    this.bound_gene = ""

    this.cell_info = null
    this.plotted_cellinfo = {
	'gene_key': '',
	'cells': [], 
    }

    this.CellInfo = null;
}

Gene_Info.prototype.init = function()
{

    var self = this;

    $(self.dom_node).find("#Gene-CellInfoButton").on('click', function()
    {
	self.addCellInfo()
    });



}

Gene_Info.prototype.update = function(updates)
{
    var item_type = get_global_status('plotted_item_type')
    var gene = get_global_status('plotted_item')

    if(item_type !== 'gene' || gene === this.bound_gene){
        return;
    }

    this.bound_gene = gene

    $(self.dom_node).find("GeneInfoButton").click();

    this.bound_gene = gene

    $(this.title).html(this.bound_gene)

    var urllink = "http://www.genecards.org/cgi-bin/carddisp.pl?gene=" + this.bound_gene;
    $(this.title).html(
        "<a href=" + urllink + " target='_blank'>" + this.bound_gene + "</a>"
    )

    var gene_exp = _.values(get_global_data('plotted_values'))

    drawDistChart(this.chart, gene_exp, 'Expression')

    if (this.gene_info_cell.hasClass('active')) { 
	this.addCellInfo();
    }

}

Gene_Info.prototype.addCellInfo = function() {

    // Need to create it if it isn't there
    var self = this;

    var cellinfo_div = $(self.dom_node).find('#cell-dist-div')
    cellinfo_div.show()

    if (self.CellInfo == null) { 
        self.CellInfo = new Cell_Info("gene");
    } 

    var gene_key = get_global_status('plotted_item'); // assume it's a gene
    var cells = get_global_status("selected_cells");

    // plotted heatmap is based on sig_key
    // check if we are already showing the right heatmap and don't regenerate
    var need_plot = false;

    if (self.plotted_cellinfo['gene_key'] !== gene_key){
        self.plotted_cellinfo['gene_key'] = gene_key
        need_plot = true;
    }

    if (self.plotted_cellinfo['cells'] !== cells) { 
        self.plotted_cellinfo['cells'] = cells;
        need_plot = true;
    }

    if(!need_plot){
        return $.when(true);
    }

    cellinfo_div.addClass('loading')

    self.CellInfo.update()

    cellinfo_div.removeClass("loading");

}

function Cell_Info(data_type)
{
    if (data_type == "signature") {
	this.dom_node = document.getElementById("sig-info");
    } else if (data_type == "meta") { 
	this.dom_node = document.getElementById("meta-info");
    } else if (data_type == "gene") { 
	this.dom_node = document.getElementById("gene-info");
    }

    this.data_type = data_type;
    this.title = $(this.dom_node).find("#cell-analysis-title");
    this.source = $(this.dom_node).find('#cell-source');
    this.chart = $(this.dom_node).find('#cell-dist-div').get(0);
    this.bound_cell = ""
}

Cell_Info.prototype.init = function()
{
}

Cell_Info.prototype.update = function()
{

    var self = this;
    if (this.data_type == "signature") { 

	var name = get_global_data('sig_info').name;

    } else if (this.data_type == "meta") { 

	var name = get_global_status('plotted_item')

    } else if (this.data_type == "gene") { 
	
	var name = get_global_status('plotted_item')

    } 
	
    var item_type = get_global_status('plotted_item_type')
    var cells = get_global_status('selected_cells')
    
     if (Object.keys(cells).length > 1) { 
	this.bound_cell = "Selected Subset";
    } else {
	this.bound_cell = cells[0]
    }

    $(this.title).text(name)
    $(this.source).text(this.bound_cell);

    var exp = get_global_data('plotted_values')

    var poolstatus = get_global_status("pooled")

    f_exp = []

    if (poolstatus && (this.data_type == "meta" || this.data_type == "gene")) { 
	var promises = [];
	exp = Object.keys(exp).forEach(function(k) { 
	    if (cells.indexOf(k) > -1) {
		var vals_promise = api.pool.values(k, self.data_type, name);
		promises.push(vals_promise);
	    }
	});
	
	return Promise.all(promises).then(function(vs) { 
	    console.log(vs);
	    for (var i = 0; i < vs.length; i++) {
		for (var j = 0; j < vs[i].length; j++) { 
		    f_exp.push(vs[i][j]);
		}
	    }
	    
	    console.log(f_exp);

	    drawDistChart(self.chart, f_exp, "Values")
	});
    } else {
	exp = Object.keys(exp).forEach(function(k) {
	    if (cells.indexOf(k) > -1) {
		f_exp.push(exp[k]);
	   } 
	})
	console.log(f_exp);



        drawDistChart(this.chart, f_exp, 'Values')
    }


}

function Meta_Info()
{
    this.dom_node = document.getElementById("meta-info");
    this.title = $(this.dom_node).find('#meta-analysis-title').get(0);
    this.chart = $(this.dom_node).find('#meta-dist-div').get(0);
    this.meta_info_default = $(this.dom_node).find('#meta-info-default');
    this.meta_info_cell = $(this.dom_node).find("#meta-info-cell");
    this.bound_meta = ""

    this.cell_info = null
    this.plotted_cellinfo = {
	'meta_key': '',
	'cells': [], 
    }

    this.CellInfo = null;
}

Meta_Info.prototype.init = function()
{
    var self = this;

    $(self.dom_node).find("#Meta-CellInfoButton").on('click', function()
    {
	self.addCellInfo()
    });


}

Meta_Info.prototype.update = function(updates)
{
    var item_type = get_global_status('plotted_item_type')
    var meta = get_global_status('plotted_item')

    if(item_type !== 'meta' || meta === this.bound_gene){
        return;
    }

    $(self.dom_node).find("MetaInfoButton").click();

    this.bound_meta = meta

    $(this.title).html(this.bound_meta)

    var meta_vals = _.values(get_global_data('plotted_values'))

    drawDistChart(this.chart, meta_vals)

    if (this.meta_info_cell.hasClass('active')) { 
	this.addCellInfo();
    }

}

Meta_Info.prototype.addCellInfo = function() {

    // Need to create it if it isn't there
    var self = this;

    var cellinfo_div = $(self.dom_node).find('#cell-dist-div')
    cellinfo_div.show()

    if (self.CellInfo == null) { 
        self.CellInfo = new Cell_Info("meta");
    } 

    var meta_key = get_global_status('plotted_item'); // assume it's meta infoe
    var cells = get_global_status("selected_cells");

    // plotted heatmap is based on sig_key
    // check if we are already showing the right heatmap and don't regenerate
    var need_plot = false;

    if (self.plotted_cellinfo['meta_key'] !== meta_key){
        self.plotted_cellinfo['meta_key'] = meta_key
        need_plot = true;
    }

    if (self.plotted_cellinfo['cells'] !== cells) { 
        self.plotted_cellinfo['cells'] = cells;
        need_plot = true;
    }

    if(!need_plot){
        return $.when(true);
    }

    cellinfo_div.addClass('loading')

    self.CellInfo.update()

    cellinfo_div.removeClass("loading");
}

//Draw Dist Scatter
function drawDistChart(parent_div, data, title) {

    var hist = create_dist(data)
    var isFactor = (typeof(data[0]) === "string") &&
                   (data[0] !== "NA")

    var x_vals = hist['centers']
    var counts = hist['counts']

    var x_axis_params;
    if(true) {
        x_axis_params = {
            type: 'category',
            categories: x_vals,
            tick: {
                rotate: 75,
                width: 100,
            },
            height: 100,
        }
    } else {
        x_axis_params = {
            type: 'indexed',
            tick: {
                rotate: 75,
                format: d3.format('.2n')
            },
            height: 100,
        }
    }

    var c3_params = {
        bindto: parent_div,
        data: {
            x: 'x',
            columns: [
                ['x'].concat(x_vals),
                ['y'].concat(counts)
            ],
            type: 'bar'
        },
        bar: {
            width: {
                ratio: 0.8
            }
        },
        axis: {
            x: x_axis_params,
            y: {
                type: 'indexed',
            }
        },
        legend: {
            show: false
        },
        size: {
            width: 400
        },
        padding: {
            right: 30,
        },
    }

    if(title !== undefined)
    {
        c3_params['title'] = {text: title}
    }

    c3.generate(c3_params)
}


function create_dist(data) {

    var counts;
    var centers;

    if((typeof(data[0]) === "string") && (data[0] !== "NA")) // Then it's a factor
    {
        var count_hist = {}
        data.forEach(function(x){
            count_hist[x] = 0;
        })
        data.forEach(function(x){
            count_hist[x] = count_hist[x] + 1;
        })
        counts = []
        centers = []
        _.forEach(count_hist, function(value, key){
            counts.push(value);
            centers.push(key);
        })

    } else {
        var num_values = 10

        // Need to filter out NA
        var data_filtered = _.filter(data, x => x !== "NA")
        var na_count = data.length - data_filtered.length
        data = data_filtered

        var data_true_min = Math.min.apply(null, data)
        var data_min = Math.max(data_true_min, 0)
        var data_max = Math.max.apply(null, data)

        var bin_width = (data_max - data_min)/num_values

        counts = Array.apply(Math, Array(num_values)).map(function() { return 0 });
        centers = Array.apply(Math, Array(num_values)).map(function() { return 0 });

        var low, high;
        var formatFn = d3.format('.2n')
        for (var i=0; i < num_values; i++) {
            low = bin_width*i + data_min
            high = bin_width*(i+1) + data_min

            data.forEach(function(d) {
                if (d >= low && d < high) {
                    counts[i] += 1;
                }
            });

            centers[i] = '[' + formatFn(low) + ', ' + formatFn(high) + ')';
        }


        if(data_true_min < 0){ // Add a "< 0" category
            low = -1e99
            high = -1e-10
            i = 0
            var lessThenZeroCounts = 0;

            data.forEach(function(d) {
                if (d >= low && d < high) {
                    lessThenZeroCounts += 1;
                }
            });

            counts = [lessThenZeroCounts].concat(counts)
            centers = ["< 0"].concat(centers)
        }

        if(na_count > 0) { // Add a "NA" category
            counts = [na_count].concat(counts)
            centers = ["NA"].concat(centers)
        }
    }

    return {'counts': counts, 'centers': centers}

}


function unselectRange() {
    global_scatter.unselect();

    global_status.subset = [];
}

function runSubsetAnalysis() {

    global_status.subset = global_scatter.getSelected();
    if (global_status.subset.length > 0) {
        api.analysis.run(global_status.subset)
    }

}

function runPCAnalysis() {

    var sig_key = global_status.plotted_signature;
    var pc1 = global_status.pc1;
    var pc2 = global_status.pc2;

    if (global_status.subset_criteria == "Rank") {
        //var sig_promise = api.signature.ranks(sig_key);
    } else {
        var sig_promise = api.signature.scores(sig_key);
    }

    var sig_info_promise = api.signature.info(sig_key);
    var proj_promise = api.pc.versus(pc1, pc2);

    return $.when(proj_promise, sig_promise, sig_info_promise)
        .then(function(projection, signature, sig_info) {

            $("#plot-title").text("PC " + pc1 + "vs. PC " + pc2);
            $("#plot-subtitle").text(sig_key);

            var points = [];
            for (var sample_label in signature) {
                var x = projection[sample_label][0];
                var y = projection[sample_label][1];
                var sig_score = signature[sample_label];
                points.push([x, y, sig_score, sample_label]);
            }

            global_scatter.setData(points, sig_info.isFactor);
        });
}

function selectRange() {
    var sig_key = global_status.plotted_signature;
    var proj_key = global_status.plotted_projection;

    var proj_promise = api.projections.coordinates(proj_key);

    if (global_status.subset_criteria == "Rank") {
        //var sig_promise = api.signature.ranks(sig_key);
    } else {
        var sig_promise = api.signature.scores(sig_key);
    }

    return $.when(sig_promise, proj_promise)
        .then(function(signature, projection) {

            var points = [];
            for (var sample_label in signature) {
                var x = projection[sample_label][0];
                var y = projection[sample_label][1];
                var sig_score = signature[sample_label];
                points.push([x, y, sig_score, sample_label]);
            }

            var cells = global_scatter.selectCellRange(points, global_status.lower_range, global_status.upper_range);
            global_status.subset = cells;

        });

}


function addToPCAnalysisBox(sig_info) {
    var dbid = "#pc-analysis-content";
    if (sig_info.isMeta) {
        $(dbid).html("");
        return;
    }

    var pcnum = global_status.plotted_pc.split(" ")[1];
    var loading_promise_pos = api.filterGroup.loadings_pos(pcnum);
    var loading_promise_neg = api.filterGroup.loadings_neg(pcnum);

    return $.when(loading_promise_pos, loading_promise_neg)
        .then(function(pos_loadings, neg_loadings) {


            var content = '<h4 style="font-size:16px; font-weight:bold">Top 5 Informative Postive Loading Vectors</h4>';
            content += "<table style='width:100%'>";
            content += "<tr><th>Gene Name</th><th>Normalized Loading</th></tr>";

            var toShow = Object.keys(pos_loadings).slice(0, 6);

            toShow.forEach(function(l) {
                //		urllink = "http://www.genecards.org/cgi-bin/carddisp.pl?gene=" + l;
                content += "<tr><td class='gene-cell'>" + l + "</td>";
                //		content += "<a href=" + urllink + " target='_blank'>" + l + "</a></td>";
                content += "<td>" + ((pos_loadings[l]*100).toPrecision(3)) + "%</td></tr>";
            });

            content += "</table><br>";

            content += '<h4 style="font-size:16px; font-weight:bold">Top 5 Informative Negative Loading Vectors</h4>';
            content += "<table style='width:100%'>";
            content += "<tr><th>Gene Name</th><th>Normalized Loading</th></tr>";


            toShow = Object.keys(neg_loadings).slice(0, 6);

            toShow.forEach(function(l) {
                var urllink = "http://www.genecards.org/cgi-bin/carddisp.pl?gene=" + l;
                content += "<tr><td class='gene-cell'>" + l + "</td>";
                //content += "<a href=" + urllink + " target='_blank'>" + l + "</a></td>";
                content += "<td>" + ((neg_loadings[l]*100).toPrecision(3)) + "%</td></tr>";
            });

            content += "</table>";
            $(dbid).html(content);
        });

}

