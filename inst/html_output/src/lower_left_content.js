/*
 * This script contains several components
 *
 * Lower_Left_Content
 *    - Sig_Info
 *      - Sig_Heatmap
 *    - Gene_Info
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
        }
    }

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

    this.cluster_options = { // Note: param values should be strings
        "KMeans": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"],
        //"PAM": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    }

    this.heatmap = null
    this.plotted_heatmap = {
        'sig_key': '',
        'proj_key': '',
        'cluster_method': '',
        'cluster_param': '',
    }

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
        'scrollY': '22vh',
        'scrollCollapse': true,
    })

    $(self.dom_node).find('#HeatmapButton').on('click', function()
    {
        self.drawHeat()
    })

    // Define cluster dropdown 
    var clust_dropdown = $(self.dom_node).find('#cluster_select_method')
    var clust_param = $(self.dom_node).find('#cluster_select_param')

    clust_dropdown.empty();
    $.each(self.cluster_options, function(name){
        clust_dropdown.append($("<option />").val(name).text(name));
    });
    clust_dropdown[0].selectedIndex = 0;

    // Call it now to initially populate second dropdown
    self.build_cluster_dropdown_param()

    //Define cluster dropdown's change function
    clust_dropdown.change(function(){
        self.build_cluster_dropdown_param()
        self.drawHeat();

    });

    //Define cluster dropdown's change function
    clust_param.change(function(){
        self.drawHeat();
    });


}

Sig_Info.prototype.update = function(updates)
{

    var sig_info = get_global_data('sig_info');
    if(sig_info.name === this.bound_sig || sig_info.isMeta)
    {
        // Needed to switch to Signature view if we plot a new projection
        // but stay on same signature
        if('plotted_projection' in updates){
            $(self.dom_node).find('#SigInfoButton').click();
        }

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
    var proj_key = get_global_status('plotted_projection');

    var cluster_method = $(self.dom_node).find('#cluster_select_method').val();
    var cluster_param = $(self.dom_node).find('#cluster_select_param').val();

    var sig_info = get_global_data('sig_info');

    // plotted heatmap is based on sig_key, proj_key, cluster_method, and cluster_param
    // check if we are already showing the right heatmap and don't regenerate
    var need_plot = false;

    if (self.plotted_heatmap['sig_key'] !== sig_key){
        self.plotted_heatmap['sig_key'] = sig_key
        need_plot = true;
    }
    if (self.plotted_heatmap['proj_key'] !== proj_key){
        self.plotted_heatmap['proj_key'] = proj_key
        need_plot = true;
    }
    if (self.plotted_heatmap['cluster_method'] !== cluster_method){
        self.plotted_heatmap['cluster_method'] = cluster_method
        need_plot = true;
    }
    if (self.plotted_heatmap['cluster_param'] !== cluster_param){
        self.plotted_heatmap['cluster_param'] = cluster_param
        need_plot = true;
    }

    if(!need_plot){
        return $.when(true);
    }

    heatmap_div.addClass('loading')

    return $.when(
        api.signature.expression(sig_key),
        api.projection.clusters(proj_key, cluster_method, cluster_param))
        .then(function(sig_expression, cluster){


            //Construct data matrix
            // TODO: sort genes

            var dataMat = sig_expression.data;
            var gene_labels = sig_expression.gene_labels;
            var sample_labels = sig_expression.sample_labels;

            var gene_signs = gene_labels.map(function(e){
                return sig_info.sigDict[e]
            });

            //var assignments = data.Clusters[proj_key][choice];
            var assignments = sample_labels.map(sample => cluster['data'][sample]);

            self.heatmap.setData(dataMat,
                assignments,
                gene_labels,
                gene_signs,
                sample_labels);

        }).always(function() {
            heatmap_div.removeClass('loading');
        });
}


function Gene_Info()
{
    this.dom_node = document.getElementById("gene-info");
    this.title = $(this.dom_node).find('#gene-analysis-title').get(0);
    this.chart = $(this.dom_node).find('#gene-dist-div').get(0);
    this.bound_gene = ""
}

Gene_Info.prototype.init = function()
{
}

Gene_Info.prototype.update = function(updates)
{
    var item_type = get_global_status('plotted_item_type')
    var gene = get_global_status('plotted_item')
    if(item_type !== 'gene' || gene === this.bound_gene){
        return;
    }

    this.bound_gene = gene

    var urllink = "http://www.genecards.org/cgi-bin/carddisp.pl?gene=" + this.bound_gene;
    $(this.title).html(
        "<a href=" + urllink + " target='_blank'>" + this.bound_gene + "</a>"
    )

    var gene_exp = _.values(get_global_data('plotted_values'))

    drawDistChart(this.chart, gene_exp, 'Expression')

}

function Meta_Info()
{
    this.dom_node = document.getElementById("meta-info");
    this.title = $(this.dom_node).find('#meta-analysis-title').get(0);
    this.chart = $(this.dom_node).find('#meta-dist-div').get(0);
    this.bound_meta = ""
}

Meta_Info.prototype.init = function()
{
}

Meta_Info.prototype.update = function(updates)
{
    var item_type = get_global_status('plotted_item_type')
    var meta = get_global_status('plotted_item')
    if(item_type !== 'meta' || meta === this.bound_gene){
        return;
    }

    this.bound_meta = meta

    $(this.title).html(this.bound_meta)

    var meta_vals = _.values(get_global_data('plotted_values'))

    drawDistChart(this.chart, meta_vals)

}

//Draw Dist Scatter
function drawDistChart(parent_div, data, title) {

    var hist = create_dist(data)
    var isFactor = typeof(data[0]) === "string"

    var x_vals = hist['centers']
    var counts = hist['counts']

    var x_axis_params;
    if(isFactor) {
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

    if(typeof(data[0]) === "string") // Then it's a factor
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
        var data_min = Math.min.apply(null, data)
        var data_max = Math.max.apply(null, data)

        var bin_width = (data_max - data_min)/num_values

        counts = Array.apply(Math, Array(num_values)).map(function() { return 0 });
        centers = Array.apply(Math, Array(num_values)).map(function() { return 0 });

        for (var i=0; i < num_values; i++) {
            var low = bin_width*i + data_min
            var high = bin_width*(i+1) + data_min

            data.forEach(function(d) {
                if (d >= low && d < high) {
                    counts[i] += 1;
                }
            });

            centers[i] = (high + low)/2
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

    var proj_promise = api.projection.coordinates(proj_key);

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

