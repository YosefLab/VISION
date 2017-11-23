var global_stack = [];

var global_status = {};
global_status.plotted_projection = "";
global_status.plotted_signature = "";
global_status.plotted_pc = "PC 1";
global_status.sorted_column = "";
global_status.signature_filter = "";
global_status.scatterColorOption = "rank";
global_status.filter_group = "";
global_status.filter_group_genes = []
global_status.selected_gene = "";
global_status.upper_range = "";
global_status.lower_range = "";
global_status.pc1 = "";
global_status.pc2 = "";
global_status.subset = [];
global_status.subset_criteria = "Rank";
global_status.main_vis = "sigvp";

var global_data = {};
global_data.sigIsPrecomputed = {};

var global_scatter = {};
var global_heatmap = {};

// Keys are cluster methods
// Values are list of allowed method parameter
var cluster_options = { // Note: param values should be strings
    "KMeans": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"],
    //"PAM": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
}

$(window).resize(function()
{
    $('#scatter_div').children().remove();
    global_scatter = new ColorScatter('#scatter_div', true, true);

    if($('#heatmap_div').is(":visible"))
    {
        $('#heatmap_div').find('svg').remove();
        global_heatmap = new HeatMap('#heatmap_div');
    } 

    if ($('#tree_div').is(":visible"))
    { 
        $('#tree_div').find('svg').remove();
    }

    //Link the scatter/heatmap
    global_scatter.hovered_links.push(global_heatmap);
    global_heatmap.hovered_links.push(global_scatter);
    //global_tree.hovered_links.push(global_tree);

    //Render
    drawChart();
    //drawHeat();


});

function doneTyping()
{
    var val = global_status.signature_filter.toLowerCase();
    var vals = val.split(",");
    vals = vals.map(function(str){return str.trim();})
        .filter(function(str){ return str.length > 0;});

    var tablerows = $('#table_div_container table').find('tr');
    tablerows.removeClass('hidden');

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
        }).addClass('hidden');
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
        }).addClass('hidden');
    }

    /*tablerows.removeClass('altRow')
            .not('.hidden').filter(':odd').addClass('altRow');
            */
}

function doneTyping_Gene()
{
    var val = global_status.gene_filter.toLowerCase();
    var vals = val.split(",");
    vals = vals.map(function(str){return str.trim();})
        .filter(function(str){ return str.length > 0;});

    var tablerows = $('#gene-table').find('tr');
    tablerows.removeClass('hidden');


    var posvals = vals.filter(function(str){ return str[0] != '!';});
    var negvals = vals.filter(function(str){ return str[0] == '!';})
        .map(function(str){ return str.slice(1);})
        .filter( function(str){return str.length > 0;});

    if(posvals.length > 0){
        tablerows.filter(function(i, element){
            if(i == 0){return false;} // Don't filter out the header row
            var gene_text = $(element).children('td').first().html().toLowerCase();
            for(var j = 0; j < posvals.length; j++)
            {
                if(gene_text.indexOf(posvals[j]) > -1)
                {
                    return false;
                }
            }
            return true;
        }).addClass('hidden');
    }

    if(negvals.length > 0){
        tablerows.filter(function(i, element){
            if(i == 0){return false;} // Don't filter out the header row
            var gene_text = $(element).children('td').first().html().toLowerCase();
            for(var j = 0; j < negvals.length; j++)
            {
                if(gene_text.indexOf(negvals[j]) > -1)
                {
                    return true;
                }
            }
            return false;
        }).addClass('hidden');
    }

    /*tablerows.removeClass('altRow')
            .not('.hidden').filter(':odd').addClass('altRow');
            */
}

window.onload = function()
{


    //Define some globals
    global_scatter = new ColorScatter("#scatter_div", true, true);
    global_heatmap = new HeatMap("#heatmap_div");

    //Link the scatter/heatmap
    global_scatter.hovered_links.push(global_heatmap);
    global_heatmap.hovered_links.push(global_scatter);
    //global_tree.hovered_links.push(global_scatter);

    global_status.precomputed = false;


    var clusters_promise = api.signature.clusters(global_status.precomputed, "1")
        .then(function(cls) {

            var clusarr = Object.keys( cls ).map(function ( key ) {return cls[key];});
            var clusmax = Math.max.apply(null, clusarr);

            for (var curr_cl = 1; curr_cl <= clusmax; curr_cl++) {
                // Create new table and add to table_div
                var new_table_div = document.createElement("div");
                new_table_div.setAttribute("style", "height=calc((100vh - 88px) / 2)");
                new_table_div.setAttribute("class", "table_div");
                new_table_div.setAttribute("style", "overflow: hidden");

                var table_div_container = document.getElementById("table_div_container");

                var new_table = document.createElement("table");
                new_table.setAttribute("id", "table"+ curr_cl);
                new_table.setAttribute("class", "sig-cluster-table");

                var thead = document.createElement("thead");
                if (curr_cl == 1) {
                    var tr = document.createElement("tr");
                    tr.setAttribute("id", "proj_row");
                    thead.appendChild(tr);
                }

                var tbody = document.createElement("tbody");

                new_table.appendChild(thead);
                new_table.appendChild(tbody);


                new_table_div.appendChild(new_table);
                table_div_container.appendChild(new_table_div);
            }

        });

    $("#filter_dropdown").change(function() {
        global_status.filter_group = $(this).val();

        api.filterGroup.genes(global_status.filter_group)
            .then(function(genes) {
                global_status.filter_group_genes = genes;
                updateMenuBar();
            });

        createTableFromData().then(function() {
            var pcnum = global_status.plotted_pc.split(" ")[1];
            $("#data-analysis-title").text(global_status.plotted_signature);
            $("#pc-analysis-title").text("Principal Component " + pcnum);

            drawChart();
            drawHeat();
            drawTree();
        });

    });	



    $("#subset-criteria").change(function() {
        global_status.subset_criteria = $(this).val();
    });

    var filterSig = $('#sig_filt_input');
    var filterSigTimer;
    var filterSigTimer_Timeout = 500;

    filterSig.on('input', function(){
        global_status.signature_filter = this.value;
        clearTimeout(filterSigTimer);
        filterSigTimer = setTimeout(doneTyping, filterSigTimer_Timeout);
    });

    var filterGene = $("#gene_filt_input");
    var filterGeneTimer;
    var filterGeneTimer_Timeout = 500;

    filterGene.on("input", function() {
        global_status.gene_filter = this.value;
        clearTimeout(filterGeneTimer);
        filterGeneTimer = setTimeout(doneTyping_Gene, filterGeneTimer_Timeout);
    })

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

    // Define cluster dropdown 
    var clust_dropdown = $('#cluster_select_method');
    clust_dropdown.empty();
    $.each(cluster_options, function(name){
        clust_dropdown.append($("<option />").val(name).text(name));
    });
    clust_dropdown[0].selectedIndex = 0; // Restore index on model change

    var build_cluster_dropdown_param = function()
    {
        // Rebuild the 'param' if the first dropdown is changed
        var vals = cluster_options[$('#cluster_select_method').val()]
        var clust_dropdown_param = $('#cluster_select_param');
        var old_val = clust_dropdown_param.val()

        clust_dropdown_param.empty();
        for(var i=0; i<vals.length; i++){
            clust_dropdown_param.append($("<option />").val(vals[i]).text(vals[i]));
        }

        if(vals.indexOf(old_val) > -1){
            clust_dropdown_param[0].selectedIndex = vals.indexOf(old_val)
        }
    }

    build_cluster_dropdown_param() // Call it now to initially populate it

    //Define cluster dropdown's change function
    $('#cluster_select_method').change(function(){
        build_cluster_dropdown_param()
        drawHeat();

    });

    //Define cluster dropdown's change function
    $('#cluster_select_param').change(function(){
        drawHeat();
    });

    //Define color option (for scatter) change function
    $('input[name=scatterColorButtons]').change(function(){
        var val = $('input[name=scatterColorButtons]:checked').val();
        global_status.scatterColorOption = val;
        drawChart();
    });

    $("#reload_heatmap").on("click", function() {
        console.log('here');
        drawHeat(); 
    });

    //Enable Toggling of Lasso Select
    $("#lasso-select").on("click", function() {
        var tog = document.getElementById("lasso-select").innerHTML;
        if (tog == "Enable Lasso Select") {
            global_scatter.toggleLasso(true);
            document.getElementById("lasso-select").innerHTML = "Disable Lasso Select";
        } else {
            global_scatter.toggleLasso(false);
            document.getElementById("lasso-select").innerHTML = "Enable Lasso Select";
        }
    });

    // Create listeners for main visualization modes
    $("#proj_tab").on("click", function() {
        global_status.main_vis = "sigvp";
        createTableFromData();
        drawChart();
    });

    $("#pc_tab").on("click", function() {

        global_status.main_vis = "pcannotator"
        createTableFromData();
        drawChart();


    });

    $("#tree_tab").on("click", function() {
        global_status.main_vis = "tree";
        createTableFromData();
        drawChart();
    });

    $("#gene_tab").on("click", function() {
        $('#gene-analysis-title').text("Gene Analysis");
        addToGeneBox();

        var filterGene = $('#gene_filt_input');
        filterGene.detach();

        // (Re)Create filter signature box
        var th = $('#filter_genes_container');
        $(th).append(filterGene);
        filterGene.show();
        filterGene.trigger('input'); 

    });

    // Make some service calls here
    // Get the list of filter groups
    var filterGroupPromise = api.filterGroup.list()
        .then(function(filters){

            for(var i = 0; i < filters.length; i++){
                var filter = filters[i]
                var option = $(document.createElement("option"));
                option.text(filter).val(filter);
                $('#filter_dropdown').append(option);
            }

        });

    var criteriaList = ["Rank", "Value"];
    for (var i = 0; i < criteriaList.length; i++) {
        var criteria = criteriaList[i];
        var option = $(document.createElement("option"));
        option.text(criteria).val(criteria);
        $("#subset-criteria").append(option);
    }

    // Get the 'isPrecomputed' vector for signatures
    var sigIsPrecomputedPromise = api.signature.listPrecomputed()
        .then(function(sigIsPrecomputed) {
            global_data.sigIsPrecomputed = sigIsPrecomputed;
        });

    // When it's all done, run this
    $.when(filterGroupPromise, sigIsPrecomputedPromise)
        .then(function(){
            $("#filter_dropdown").change() // Change event updates the table
        });

    drawHeat();

};

Element.prototype.remove = function() {
    this.parentElement.removeChild(this);
}

function addSigClusterDivs() {

    $(".table_div").remove()

    if (global_status.precomputed) {
        var new_table_div = document.createElement("div");
        new_table_div.setAttribute("style", "height=calc((100vh-88px) / 2)");

        new_table_div.setAttribute("id", "new_table_div");
        new_table_div.setAttribute("class", "table_div");


        var table_div_container = document.getElementById("table_div_container");

        var new_table = document.createElement("table");
        new_table.setAttribute("id", "precomp-table");
        new_table.setAttribute("class", "sig-cluster-table");

        var thead = document.createElement("thead");
        var tr = document.createElement("tr");
        tr.setAttribute("id", "proj_row");
        thead.appendChild(tr);

        var tbody = document.createElement("tbody");

        new_table.appendChild(thead);
        new_table.appendChild(tbody);


        new_table_div.appendChild(new_table);
        table_div_container.appendChild(new_table_div);

        return;
    } 

    var num_clusters_promise = api.signature.clusters(global_status.precomputed, global_status.filter_group)
        .then(function(cls) {

            var clusarr = Object.keys( cls ).map(function ( key ) {return cls[key];});
            var clusmax = Math.max.apply(null, clusarr);

            for (var curr_cl = 1; curr_cl <= clusmax; curr_cl++) {
                // Create new table and add to table_div
                var new_table_div = document.createElement("div");
                new_table_div.setAttribute("style", "height=calc((100vh - 88px) / 2)");
                new_table_div.setAttribute("id", "new_table_div");
                new_table_div.setAttribute("class", "table_div");


                var table_div_container = document.getElementById("table_div_container");

                var new_table = document.createElement("table");
                new_table.setAttribute("id", "table"+ curr_cl);
                new_table.setAttribute("class", "sig-cluster-table");

                var thead = document.createElement("thead");
                if (curr_cl == 1) {
                    var tr = document.createElement("tr");
                    tr.setAttribute("id", "proj_row");
                    thead.appendChild(tr);
                }

                var tbody = document.createElement("tbody");

                new_table.appendChild(thead);
                new_table.appendChild(tbody);


                new_table_div.appendChild(new_table);
                table_div_container.appendChild(new_table_div);
            }

            //table_div_container.appendChild(document.createElement("br"));

        });

}

function updateCurrentSelections(matrix)
{

    if (global_status.main_vis == "sigvp") {
        // If no sort specified, sort by PCA: 1,2
        if(matrix.proj_labels.indexOf(global_status.sorted_column) == -1 )
        {
            if(matrix.proj_labels.indexOf("PCA: 1,2" ) != -1){
                global_status.sorted_column = "PCA: 1,2";
            }
            else
            {
                global_status.sorted_column = "";
            }
        }	

        // If no projection specified, select the PCA: 1,2 projection
        if(matrix.proj_labels.indexOf(global_status.plotted_projection) == -1)
        {
            if(matrix.proj_labels.indexOf("PCA: 1,2" ) != -1){
                global_status.plotted_projection = "PCA: 1,2";
            }
            else
            {
                global_status.plotted_projection = matrix.proj_labels[0];
            }
        }
    } else if (global_status.main_vis == "pcannotator") {

        if (matrix.proj_labels.indexOf(global_status.sorted_column) == -1) {
            if (matrix.proj_labels.indexOf("PC 1") != -1) {
                global_status.sorted_column = "PC 1";
            } else {
                global_status.sorted_column = "";
            }
        }

        if (matrix.proj_labels.indexOf(global_status.plotted_pc) == -1) {
            if (matrix.proj_labels.indexOf("PC 1") != -1) {
                global_status.plotted_pc = "PC 1";
            } else {
                global_status.plotted_pc = matrix.proj_labels[0];
            }
        }
    }		

    // If no signature specified, selected the top signature
    if(matrix.sig_labels.indexOf(global_status.plotted_signature) == -1)
    {
        //Select top signature of sorted projection by default
        var j = matrix.proj_labels.indexOf(global_status.plotted_projection);
        var s_i = matrix.data.map(function(e){return e[j];}).argSort();

        if (typeof(matrix.sig_labels) == "string") {
            global_status.plotted_signature = matrix.sig_labels;
        } else {
            global_status.plotted_signature = matrix.sig_labels[s_i[0]];
        }
    }

}

function addToPCAnalysisBox(sig_info) {
    var dbid = "#pc-analysis-content";
    if (sig_info.isPrecomputed) {
        $(dbid).html("");
        return;
    }

    var pcnum = global_status.plotted_pc.split(" ")[1];
    var loading_promise_pos = api.filterGroup.loadings_pos(global_status.filter_group, pcnum);
    var loading_promise_neg = api.filterGroup.loadings_neg(global_status.filter_group, pcnum);

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

function addToGeneBox() {

    var dbid = "#gene-analysis-content";

    var gene_promise = api.expression.genes.list()

    return $.when(gene_promise
        .then(function(genes) {


            var gene_title = document.getElementById("selected-gene-title");
            gene_title.innerHTML = "Selected Gene: " + global_status.selected_gene;

            var gene_card = document.getElementById("view-gene-card");
            var urllink = "http://www.genecards.org/cgi-bin/carddisp.pl?gene=" + global_status.selected_gene;
            gene_card.innerHTML = "<a href=" + urllink + " target='_blank'>" + "\tView GeneCard" + "</a></td>";


            var content = "<table id='gene-table' style='width:75%; margin-left: 5%;'>";
            content += "<tr><th>Gene Name</th</tr>";

            for (var i = 0; i < genes.length; i++) {
                content += "<tr><td class='gene-cell'>" + genes[i] + "</td></tr>";
            }

            content += "</table>";

            $(dbid).html(content);

            $(".gene-cell").on("click", function(d) {
                global_status.selected_gene = d.target.innerHTML;
                global_status.scatterColorOption = "gene";
                addToGeneBox();
                drawChart();
                drawDistChart();
            });

        }));


}

function addToDataAnalysisBox(sig_info) {
    var dbid = "#data-analysis-content";
    if (sig_info.isPrecomputed) {
        $(dbid).html("");
        return;
    }


    var content = '<p> Source: ' + sig_info.source + '<br>';
    content += "<table style='width:100%'>";
    content += "<tr><th>Gene Name</th><th>Sign</th></tr>";


    Object.keys(sig_info.sigDict).forEach(function(key) {
        if (global_status.selected_gene == "") {
            global_status.selected_gene = key;
        }
        var urllink = "http://www.genecards.org/cgi-bin/carddisp.pl?gene=" + key;
        content += "<tr><td>";
        content += "<a href=" + urllink + " target='_blank'>" + key + "</a></td>";
        if (sig_info.sigDict[key] == 1) {
            content += "<td>+</td></tr>";
        } else {
            content += "<td>-</td></tr>";
        }
    });

    content += "</table>";

    $(dbid).html(content);
}

function updateMenuBar()  // Updates the ### Genes button in the menu bar
{
    var button = $('#show_genes_button');
    var num_genes = global_status.filter_group_genes.length
    button.text(num_genes + " Genes");
}

// Function that's triggered when clicking on table header column
function sortByColumn(col_name)  
{
    global_status.sorted_column = col_name;
    createTableFromData();
}

// Function that's triggered when clicking on table cell
function tableClickFunction(row_key, col_key)
{
    global_status.plotted_signature  = row_key;
    if (global_status.main_vis == "pcannotator") {
        global_status.plotted_pc = col_key;
    } else {
        global_status.plotted_projection = col_key;
    }
    $('#data-analysis-title').text(global_status.plotted_signature);
    drawChart();
    //drawHeat();
    drawTree();
}

//Draw Dist Scatter
function drawDistChart() {

    var selected_gene = global_status.selected_gene;
    var val_promise = api.expression.gene(selected_gene);

    return val_promise.then(function(values) {

        var expr_data = _.values(values)

        var hist = create_dist(expr_data);

        var x_vals = hist['centers']
        var counts = hist['counts']

        c3.generate({
            bindto: '#gene_dist_div',
            data: {
                x: 'x',
                columns: [
                    ['x'].concat(x_vals),
                    [selected_gene].concat(counts)
                ],
                type: 'bar'
            },
            bar: {
                width: {
                    ratio: 0.8
                }
            },
            axis: {
                x: {
                    type: 'indexed',
                    tick: {
                        rotate: 75,
                        format: d3.format('.2n')
                    },
                },
                y: {
                    type: 'indexed',
                }
            },
            legend: {
                show: false
            }
        })
    });
}


// Draw the scatter plot
function drawChart() {

    var status_copy = $.extend({}, global_status);
    global_stack.push(status_copy);

    var sig_key = global_status.plotted_signature;
    var proj_key = global_status.plotted_projection;
    var pc_key = global_status.plotted_pc;
    var filter_group = global_status.filter_group;


    var val_promise;
    if(global_status.scatterColorOption == "value" || 
        global_data.sigIsPrecomputed[sig_key])
    {val_promise = api.signature.scores(sig_key)}

    if(global_status.scatterColorOption == "rank") 
    {val_promise = api.signature.ranks(sig_key)}

    if (global_status.scatterColorOption == "gene") {
        val_promise = api.expression.gene(global_status.selected_gene)
    }

    if (global_status.main_vis == "sigvp") {

        if(sig_key.length == 0 && proj_key.length == 0){
            $('#plot-title').text("");
            $('#plot-subtitle').text("");
            global_scatter.setData([], false);
            return $().promise()
        }	

        var proj_promise = api.projection.coordinates(filter_group, proj_key);

        var sig_info_promise = api.signature.info(sig_key)

        return $.when(proj_promise, val_promise, sig_info_promise) // Runs when both are completed
            .then(function(projection, values, sig_info){

                $('#plot-title').text(proj_key);
                $('#plot-subtitle').text(sig_key);

                var points = [];

                for(var sample_label in values){
                    var x = projection[sample_label][0]
                    var y = projection[sample_label][1]
                    var sig_score = values[sample_label]
                    points.push([x, y, sig_score, sample_label]);
                }

                global_scatter.setData(points, sig_info.isFactor);

                addToDataAnalysisBox(sig_info);

                addToPCAnalysisBox(sig_info);

            });
    } else if (global_status.main_vis == "tree") {

        if(sig_key.length == 0 && proj_key.length == 0){
            $('#plot-title').text("");
            $('#plot-subtitle').text("");
            global_scatter.setData([], false);
            return $().promise()
        }	

        proj_promise = api.tree.coordinates(filter_group, proj_key);

        sig_info_promise = api.signature.info(sig_key);

        var tree_points = api.tree.tree_points(filter_group, proj_key);
        var tree_adjlist = api.tree.tree(filter_group)

        return $.when(proj_promise, val_promise, sig_info_promise, tree_points, tree_adjlist) // Runs when both are completed
            .then(function(projection, values, sig_info, treep, treel){


                // Massage treep for easier D3 binding

                tree_points = []

                $('#plot-title').text(proj_key);
                $('#plot-subtitle').text(sig_key);

                var points = [];

                for(var sample_label in values){
                    var x = projection[sample_label][0]
                    var y = projection[sample_label][1]
                    var sig_score = values[sample_label]
                    points.push([x, y, sig_score, sample_label]);
                }

                for (var i = 0; i < treep[0].length; i++) {
                    var x = treep[0][i];
                    var y = treep[1][i];
                    var score = 0
                    points.push([x, y, score, "Node " + i, "Tree"]);
                }

                //global_scatter.setData(points, sig_info.isFactor);

                global_scatter.addTree(points, treel);

                addToDataAnalysisBox(sig_info);

                addToPCAnalysisBox(sig_info);

            });
    } else if (global_status.main_vis == "pcannotator") {

        var sig_info_promise = api.signature.info(sig_key);
        var pc_key = global_status.plotted_pc.split(" ")[1];

        if (sig_key.length == 0 && pc_key.length == 0) {
            $("#plot-title").text("");
            $("#plot-subtitle").text("");
            global_scatter.setData([], false);
            return $().promise();
        }

        var pc_promise = api.pc.coordinates(filter_group, sig_key, pc_key);
        var sig_info_promise = api.signature.info(sig_key);

        return $.when(pc_promise, val_promise, sig_info_promise)
            .then(function(projection, values, sig_info) {


                $("#plot-title").text("Principal Component ".concat(pc_key));
                $("#plot-subtitle").text(sig_key);

                var points = []
                for(var sample_label in values){
                    var x = projection[sample_label][0]
                    var y = projection[sample_label][1]
                    var sig_score = values[sample_label]
                    points.push([x, y, sig_score, sample_label]);
                }

                global_scatter.setData(points, sig_info.isFactor);

                addToDataAnalysisBox(sig_info);

                addToPCAnalysisBox(sig_info);
            });

    } else {
        return;
    }

}

function create_dist(data) {

    var num_values = 10
    var data_min = Math.min.apply(null, data)
    var data_max = Math.max.apply(null, data)

    var bin_width = (data_max - data_min)/num_values

    var counts = Array.apply(Math, Array(num_values)).map(function() { return 0 });
    var centers = Array.apply(Math, Array(num_values)).map(function() { return 0 });

    for (var i=0; i < num_values; i++) {
        var low = bin_width*i
        var high = bin_width*(i+1)

        data.forEach(function(d) {
            if (d >= low && d < high) {
                counts[i] += 1;
            }
        });

        centers[i] = (high + low)/2
    }

    return {'counts': counts, 'centers': centers}

}

// Draw the Tree
function drawTree() {
    var sig_key = global_status.plotted_signature;
    var proj_key = global_status.plotted_projection;
    var filter_group = global_status.filter_group;

    if (!$('#tree-options').hasClass('active')) {
        $("#tree_div").hide();
        $('#instructions').show();
        return;
    }
    $("#instructions").hide();
    if (sig_key.length == 0) {
        $('#tree_div').hide();
        return $().promise();
    }
    $("#tree_div rect").show();
    return $.when(api.signature.info(sig_key),
        api.signature.expression(sig_key),
        api.projection.tree(filter_group, proj_key))
        .then(function(sig_info, sig_expression, tree_data) {

            if (sig_info.isPrecomputed) {
                return;
            }

            if (!$('#tree_div').is(":visible"))
            {
                $("#tree_div").find("svg").remove();
                $("#tree_div").show();
                //global_tree.hovered_links.push(global_scatter);
            }

            var sample_labels = sig_expression.sample_labels;

            //return global_tree.setData(tree_data[1], tree_data[3], sample_labels);
        });
}


// Draw the heatmap
function drawHeat(){
    var sig_key = global_status.plotted_signature;
    var proj_key = global_status.plotted_projection;
    var filter_group = global_status.filter_group
    var cluster_method = $('#cluster_select_method').val();
    var cluster_param = $('#cluster_select_param').val();

    if(sig_key.length == 0){
        $('#heatmap_div').hide();
        return $().promise();
    } else if (global_status.main_vis == "pcannotator") {
        $("#heatmap_div").hide();
        return $().promise();
    }

    $("#heatmap_div rect").show();
    return $.when(api.signature.info(sig_key),
        api.signature.expression(sig_key),
        api.projection.clusters(filter_group, proj_key, cluster_method, cluster_param))
        .then(function(sig_info, sig_expression, cluster){

            if(sig_info.isPrecomputed){
                $('#heatmap_div').hide();
                return
            }

            // Heatmap doesn't show for precomputed sigs
            // Need to recreate it if it isn't there
            if( !$('#heatmap_div').is(":visible"))
            {
                $('#heatmap_div').find('svg').remove();
                $('#heatmap_div').show();
                global_heatmap = new HeatMap('#heatmap_div');
                global_scatter.hovered_links.push(global_heatmap);
                global_heatmap.hovered_links.push(global_scatter);
            }


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

            global_heatmap.setData(dataMat,
                assignments,
                gene_labels,
                gene_signs,
                sample_labels);

        });
}


function createTableFromData()
{


    addSigClusterDivs();
    var matrix_promise;
    if (global_status.main_vis == "sigvp") {
        matrix_promise = api.filterGroup.sigProjMatrixP(global_status.filter_group, global_status.precomputed);
    } else if (global_status.main_vis == "tree") {
        matrix_promise = api.filterGroup.treeSigProjMatrixP(global_status.filter_group, global_status.precomputed);
    } else {
        matrix_promise = api.filterGroup.pCorr(global_status.filter_group, global_status.precomputed);
    }


    if (global_status.precomputed) {

        return $.when(matrix_promise
            .then(function(matrix) {

                var header_row = d3.select("#precomp-table").select("thead").select("#proj_row").selectAll("th")
                    .data([""].concat(matrix.proj_labels));
                header_row.enter().append('th');
                header_row.html(function(d){return "<div>" + d+"</div>";})
                    .filter(function(d,i) {return i > 0;})
                    .on("click", function(d) { sortByColumn(d);});
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

                var sort_col = matrix.proj_labels.indexOf(global_status.sorted_column);
                if(sort_col > -1){
                    var sortFun = function(a,b){
                        var a_precomp = global_data.sigIsPrecomputed[a[0]];
                        var b_precomp = global_data.sigIsPrecomputed[b[0]];
                        if(a_precomp && b_precomp || !a_precomp && !b_precomp){
                            if (global_status.main_vis == "sigvp") {
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

                if (global_status.main_vis == "pcannotator") {
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


                var content_rows = d3.select('#precomp-table').select('tbody').selectAll('tr')
                    .data(formatted_data_w_row_labels);
                content_rows.enter().append('tr');
                content_rows.exit().remove();

                var content_row = content_rows.selectAll("td")
                    .data(function(d){return [d[0]].concat(d[1]);})

                content_row.enter().append('td');
                content_row.exit().remove();

                if (global_status.main_vis == "pcannotator") {
                    content_row
                        .filter(function(d,i) { return i > 0;})
                        .text(function(d){
                            if(d.val > .5) { return "> .5";}
                            else if(d.val < .5) { return d.val.toFixed(2);}
                            else {return d.val.toPrecision(2);}
                        })
                        .style('background-color', function(d){return colorScale(d.val);})
                        .on("click", function(d){tableClickFunction(matrix.sig_labels[d.row], matrix.proj_labels[d.col])});

                } else {
                    content_row
                        .filter(function(d,i) { return i > 0;})
                        .text(function(d){
                            if(d.val < -50) { return "< -50";}
                            else if(d.val > -1) { return d.val.toFixed(2);}
                            else {return d.val.toPrecision(2);}
                        })
                        .style('background-color', function(d){return colorScale(d.val);})
                        .on("click", function(d){tableClickFunction(matrix.sig_labels[d.row], matrix.proj_labels[d.col])});
                }

                // Make signature names 
                content_row.filter(function(d,i) { return i == 0;})
                    .text(function(d){return d;});

            }));
    }

    return $.when(matrix_promise,
        api.signature.clusters(global_status.precomputed, global_status.filter_group))
        .then(function(matrix, cls){

            var clusarr = Object.keys( cls ).map(function(key) { return cls[key]; });
            var clusmax = Math.max.apply(null, clusarr);
            updateCurrentSelections(matrix);

            // Detach filter sig box for later
            var filterSig = $('#sig_filt_input');
            filterSig.detach();

            // (Re)Create filter signature box
            var th = $('#filter_sigs_container');
            $(th).append(filterSig);
            filterSig.show();
            filterSig.trigger('input'); 


            for (var curr_cl = 1; curr_cl <= clusmax; curr_cl++) {

                // Create the Header row
                if (curr_cl == 1) {
                    var header_row = d3.select("#table"+curr_cl).select("thead").select("#proj_row").selectAll("th")
                        .data([""].concat(matrix.proj_labels));
                    header_row.enter().append('th');
                    header_row.html(function(d){return "<div>" + d+"</div>";})
                        .filter(function(d,i) {return i > 0;})
                        .on("click", function(d) { sortByColumn(d);});
                    header_row.exit().remove();
                }


                if (typeof(matrix.sig_labels) == "string") {
                    matrix.sig_labels = [matrix.sig_labels];
                }

                // Format cell data for better d3 binding
                var sig_labels = matrix.sig_labels.filter(function(x) { return cls[x] == curr_cl; });
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

                var sort_col = matrix.proj_labels.indexOf(global_status.sorted_column);
                if(sort_col > -1){
                    var sortFun = function(a,b){
                        var a_precomp = global_data.sigIsPrecomputed[a[0]];
                        var b_precomp = global_data.sigIsPrecomputed[b[0]];
                        if(a_precomp && b_precomp || !a_precomp && !b_precomp){
                            if (global_status.main_vis == "sigvp") {
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

                if (global_status.main_vis == "pcannotator") {
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

                if (global_status.main_vis == "pcannotator") {
                    content_row
                        .filter(function(d,i) { return i > 0;})
                        .text(function(d){
                            if(d.val > .8) { return "> .8";}
                            else if(d.val < .8) { return d.val.toFixed(2);}
                            else {return d.val.toPrecision(2);}
                        })
                        .style('background-color', function(d){return colorScale(d.val);})
                        .on("click", function(d){tableClickFunction(matrix.sig_labels[d.row], matrix.proj_labels[d.col])});

                } else {
                    content_row
                        .filter(function(d,i) { return i > 0;})
                        .text(function(d){
                            if(d.val < -50) { return "< -50";}
                            else if(d.val > -1) { return d.val.toFixed(2);}
                            else {return d.val.toPrecision(2);}
                        })
                        .style('background-color', function(d){return colorScale(d.val);})
                        .on("click", function(d){tableClickFunction(matrix.sig_labels[d.row], matrix.proj_labels[d.col])});
                }

                // Make signature names click-able
                content_row.filter(function(d,i) { return i == 0;})
                    .text(function(d){return d;})

                $('#table'+curr_cl).children('tbody').children("tr:not(:first-child)").hide();
                $("#table"+curr_cl).children("tbody").children("tr:first-child").children("td:first-child").attr("id", "sigclust_" + curr_cl);
                $("#table"+curr_cl).children('tbody').children('tr:first-child').children("td:first-child")
                    .on("click", function(d) { clickSummaryRow(this); });

                // Add '>' sign to top sig name to indiciate expandability
                $("#table"+curr_cl).children("tbody").children("tr:first-child").children("td:first-child")
                    .text(function(i, origText) {
                        return origText + " \u25B6 ";
                    });

            }


            $(".sigclust").on("mouseover", function(d) {
                tooltip.showTooltip("Click To Toggle Cluster Display", d);
            })
                .on("mouseout", function(d) {
                    tooltip.hideTooltip();
                });


        });
}

function clickSummaryRow(d) {
    var clust = d.id.split("_")[1];

    var table_id = $("#table"+clust).children("tbody");
    if (table_id.children("tr:not(:first-child)").is(":visible")) {
        table_id.children("tr:not(:first-child)").hide();
        table_id.children("tr:first-child").children("td:first-child")
            .text(function(i, origText) {
                origText = origText.split(" ")[0];
                return origText + " \u25B6 ";
            });
    } else {
        table_id.children("tr:not(:first-child)").show();
        table_id.children("tr:first-child").children("td:first-child")
            .text(function(i, origText) {
                origText = origText.split(" ")[0];
                return origText + " \u25BC";
            });
    }
}

function computeDataAverages(data) {

    // We want one average per projection 
    var averages = Array.apply(null, Array(data[0].length)).map(Number.prototype.valueOf, 0);
    for (var i = 0; i < data.length; i++) {
        for (var j = 0; j < data[i].length; j++) {
            averages[j] += data[i][j].val;
        }
    }
    for (i = 0; i < averages.length; i++) {
        averages[i] /= data.length;
    }
    return averages;
}

function createSigModal(sig_key){

    return api.signature.info(sig_key).then(function(sig_info){

        var sig_data = [];
        for(var gene in sig_info['sigDict'])
        {
            sig_data.push({'Gene': gene, 'Sign': sig_info['sigDict'][gene]});
        }
        var sigModal = $('#signatureModal');
        sigModal.find('h4').text(sig_key);
        var tableRows = d3.select('#signatureModal').select('tbody').selectAll('tr')
            .data(sig_data);
        tableRows.enter().append('tr');
        tableRows.exit().remove();

        var tableCells = tableRows.selectAll('td').data(function(d){return [d.Gene, d.Sign];});
        tableCells.enter().append('td');
        tableCells.text(function(d, i){
            if(i == 0){return d;}
            else{
                if(d == 1)  {return "+";}
                if(d == -1) {return "-";}
                if(d == 0)  {return "Unsigned";}
                return "Unknown";
            }
        });

        tableCells.exit().remove();

        sigModal.modal();
    });
}

function createGeneModal()
{
    return api.filterGroup.genes(global_status.filter_group)
        .then(function(genes){

            //Calculate max width
            var width_and_index = genes.map(function(e,i){return [e.length, i]});
            width_and_index.sort(function(a,b){return Math.sign(b[0] - a[0]);});
            var top10 = width_and_index.slice(0,10).map(function(e){return genes[e[1]];});
            var widths = [];
            for(var i = 0; i < top10.length; i++)
            {
                var div = document.createElement("div");
                $(div).text(top10[i]).css("position","absolute").css("left", "-9999px");
                $('body').append(div);
                widths.push($(div).width());
            }

            var maxWidth = d3.max(widths);

            /* var geneDivs = d3.select('#geneModal').select('.modal-body').selectAll('di
                .data(genes.sort());
                */

            geneDivs.enter().append('div');
            geneDivs.exit().remove();

            geneDivs
                .text(function(d){return d;})
                .style("width", maxWidth + "px");

            $('#geneModal').modal();
        });
}

function goBack(d) {

    if (global_stack.length == 1) {
        return;
    }

    curr_status = global_stack.pop();
    global_status = global_stack.pop();
    createTableFromData();
    drawChart();
    drawHeat();
}

function changeTableView() {
    var precomp = document.getElementById("precomputed_button").innerHTML;
    if (precomp == "Show Precomputed") {
        document.getElementById("precomputed_button").innerHTML = "Show Computed";
    } else {
        document.getElementById("precomputed_button").innerHTML = "Show Precomputed";
    }

    global_status.precomputed = !(global_status.precomputed);
    return createTableFromData();
}

function exprotSigProj() {

    var sig_key = global_status.plotted_signature;
    var proj_key = global_status.plotted_projection;
    var filter_group = global_status.filter_group;

    if (sig_key.length == 0 && proj_key.length == 0) {
        $("#plot_title_div").children().eq(0).text("");
        $("#plot_title_div").children().eq(1).text("");
        global_scatter.setData([], false);
        return $().promise();
    }

    var proj_promise = api.projection.coordinates(filter_group, proj_key);

    var sig_promise;
    if (global_status.scatterColorOption == "value" || global_data.sigIsPrecomputed[sig_key]) {
        sig_promise = api.signature.scores(sig_key)
    } else if (global_status.scatterColorOption == "rank") {
        sig_promise = api.signature.ranks(sig_key)
    }

    return $.when(proj_promise, sig_promise)
        .then(function(projection, signature) {

            var points = [];
            for (var sample_label in signature) {
                var x = projection[sample_label][0];
                var y = projection[sample_label][1];
                var sig_score = signature[sample_label];
                points.push([x, y, sig_score, sample_label]);
            } 

            var lineArray = [];
            points.forEach(function(infoArray, index) {
                var line = inforArray.join(",");
                lineArray.push(index == 0 ? "data:text/csv;charset=utf-8," + line : line);
            });
            var csvContent = lineArray.join("\n");

            var encodeURI = encodeURI(csvContent);
            var downloadLink = document.createElement("a");
            downloadLink.setAttribute("href", encodedURI);
            downloadLink.setAttribute("download", proj_key + ".csv");
            downloadLink.onclick = destroyClickedElement;
            downloadLink.style.display = "none";
            document.body.appendChild(downloadLink);

            downloadLink.click();
        });

}

function destroyClickedElement(event) {
    document.body.removeChild(event.target);
}

function selectRange() {
    var sig_key = global_status.plotted_signature;
    var proj_key = global_status.plotted_projection;
    var filter_group = global_status.filter_group;

    var proj_promise = api.projection.coordinates(filter_group, proj_key);

    if (global_status.subset_criteria == "Rank") {
        var sig_promise = api.signature.ranks(sig_key);
    } else {
        sig_promise = api.signature.scores(sig_key);
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
    var filter_group = global_status.filter_group;
    var pc1 = global_status.pc1;
    var pc2 = global_status.pc2;

    if (global_status.subset_criteria == "Rank") {
        var sig_promise = api.signature.ranks(sig_key);
    } else {
        sig_promise = api.signature.scores(sig_key);
    }

    var sig_info_promise = api.signature.info(sig_key);
    var proj_promise = api.pc.versus(filter_group, pc1, pc2);

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

