/*
A library for general purpose utility functions
Requires:
   d3.js

 */



/*
Exports a zip with data in it
To be used with the FastProject results viewer
 */
function exportSigProj()
{
    var zip = new JSZip();

    var data = getDataContext();

    var choice = $('#cluster_select').val();
    var sig_key = global_status.plotted_signature;
    var proj_key = global_status.plotted_projection;

    //Convert the data that's in the scatter plot to a tab-delimited table

    var proj = data.Projections[proj_key];
    var sig = data.SigScores[sig_key];
    var cluster_assignments = data.Clusters[proj_key][choice];
    
    var table;
    if(sig.isFactor){ 
        table = [proj[0], proj[1], sig.scores, cluster_assignments];}
    else{
        table = [proj[0], proj[1], sig.scores, sig.ranks, cluster_assignments];}


    table = d3.transpose(table);
    if(sig.isFactor){
        table = [["X", "Y", "Signature Score", "Cluster: "+choice]].concat(table);}
    else{
        table = [["X", "Y", "Signature Score", "Signature Rank", "Cluster: "+choice]].concat(table);}

    table = table.map(function(x){ return x.join("\t");});
    var scatter_csv_str = table.join("\n");
    zip.file("Scatter.txt", scatter_csv_str);

    //Convert the heatmap into a tab-delimited table
    if(!sig.isPrecomputed)
    {
        var heat_data_plus = global_heatmap.data_plus;
        var heat_data_minus = global_heatmap.data_minus;
        var heat_data = [];

        //stitch together heat_data_plus and heat_data_minus to make heat_data
        for(var i=0; i<heat_data_plus.length; i++)
        {
            var new_clust = {};
            new_clust.data = heat_data_plus[i].data;
            new_clust.weight = heat_data_plus[i].weight;
            new_clust.index = heat_data_plus[i].index;
            heat_data.push(new_clust);
        }

        //stitch together heat_data_plus and heat_data_minus to make heat_data
        for(var i=0; i<heat_data_minus.length; i++)
        {
            new_clust = heat_data[i];
            new_clust.data = new_clust.data.concat(heat_data_minus[i].data);
        }
        

        var header_row = ["Gene"];
        var heat_table = [];
        for(var i = 0; i<heat_data.length; i++)
        {
            header_row.push("Cluster " + heat_data[i].index);
            heat_table.push(heat_data[i].data.map(function(x){return x.value;}));
        }

        var row_labels = heat_data[0].data.map(function(x){return x.gene});
        heat_table = [row_labels].concat(heat_table);

        heat_table = d3.transpose(heat_table);

        heat_table = [header_row].concat(heat_table);

        var heat_csv_str = heat_table.map(function(x){ return x.join("\t");});
        heat_csv_str = heat_csv_str.join("\n");
        zip.file("HeatMap.txt", heat_csv_str);
    }
   
    //Get the scatter plot and convert to a PNG
    var svg = d3.select("#scatter_div").select("svg");
    svg.attr("version", 1.1)
       .attr("xmlns", "http://www.w3.org/2000/svg");
    var svg2 = svgCopy(svg.node());

    var html_data = svg2.parentNode.innerHTML;
    zip.file("Scatter.svg", html_data);

    var imgsrc = "data:image/svg+xml;base64," + btoa(html_data);

    var image = new Image();
    image.onload = function()
    {
        var canvas = document.createElement("canvas");
        canvas.width = image.width;
        canvas.height = image.height;
        var context = canvas.getContext("2d");
        context.drawImage(image, 0,0);

        var canvasdata = canvas.toDataURL("image/png");
        //Strip off the data URI portion
        var scatter_png = canvasdata.substring(canvasdata.indexOf(",")+1);
        
        //Take the result and stick it into a zip
        
        zip.file("Scatter.png", scatter_png   , {base64: true});

        var zip_uri = "data:application/zip;base64," + zip.generate({type:"base64"});

        var a = document.createElement("a");
        a.download = sig_key+"_"+proj_key+".zip";
        a.href = zip_uri;
        a.click();

    };
    image.src = imgsrc;

}


/*
 * Creates a detached copy of an SVG image
 * Copies all applied css styles directly
 * Used for exporting an SVG
 */
function svgCopy(svg)
{
    var svg_copy = svg.cloneNode(true);

    //Need to stick it in the DOM so basic styles are computed
    var div = document.createElement("div");
    div.style.display = "none";
    document.body.appendChild(div);
    div.appendChild(svg_copy);
    
    //Apply the differences here
    apply_styles(svg, svg_copy);
    return svg_copy;
}


/*
 * Recursively copy computed styles from one node to another
 * Repeats for all children
 * Assumes nodes between original and unstyled correspond directly
 */
function apply_styles(original_node, unstyled_node)
{
    var cssStyles = getComputedStyle(original_node);
    var unstyledStyles = getComputedStyle(unstyled_node);
    var computedStyleStr = unstyled_node.getAttribute("style");
    if(computedStyleStr === null){ computedStyleStr = "";}
    var i;
    for(i = 0; i < cssStyles.length;i++)
    {
        var key = cssStyles[i];
        var value = cssStyles.getPropertyValue(key);

        //For some reason, if height or width is set explicitly as an attribute, the computed style is "auto"
        //If we then add width:auto to the elements style string, we override the attribute
        //This prevents that from happening
        if(key === "height" && original_node.hasAttribute(key)){continue;}
        if(key === "width" && original_node.hasAttribute(key)){continue;}

        if (value !== unstyledStyles.getPropertyValue(key))
        {
            computedStyleStr += key + ":" + value + ";";
        }
    }

    unstyled_node.setAttribute("style", computedStyleStr);

    var original_children = original_node.children;
    var unstyled_children = unstyled_node.children;

    for(i=0; i < original_children.length; i++)
    {
        apply_styles(original_children[i], unstyled_children[i]);
    }
}
