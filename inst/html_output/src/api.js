var api = (function(){

    /* Coordinates come in in the form:
     * [
     *  [x, y, sample_name],
     *  [x, y, sample_name2],
     *  ...
     * ]
     * Converting back to the format we want is faster in Javascript than in R
     */

    // location.pathname is either:
    // this is either "/" or "/Results.html" if running normally
    // if running from RStudio Server this will be:
    // "/p/<port>/Results.html" or "/p/<port>/"
    // Need to set this dynamically so that it works with RStudio server

    var prefix = location.pathname.replace("Results.html", "")

    var fix_coordinates = function(x){

        var ii = x[0].length - 1;

        var result =  _(x)
            .keyBy(x => x[ii])
            .mapValues(x => x.slice(0, ii))
            .value();

        return result;
    }

    var output = {};

    // Signature API

    output.signature = {}

    output.signature.info = function(sig_name){
        var query = prefix.concat("Signature/Info/")
        query = query.concat(encodeURI(sig_name));
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.scores = function(sig_name){
        var query = prefix.concat("Signature/Scores/")
        query = query.concat(encodeURI(sig_name))
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.meta = function(meta_name){
        var query = prefix.concat("Signature/Meta/")
        query = query.concat(encodeURI(meta_name))
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.expression = function(sig_name){
        var query = prefix.concat("Signature/Expression/")
        query = query.concat(encodeURI(sig_name))
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.clusters = function(meta) {
        var query
        if (meta) {
            query = prefix.concat("FilterGroup/SigClusters/Meta");
        } else {
            query = prefix.concat("FilterGroup/SigClusters/Normal");
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Pools API

    output.pool = {}

    output.pool.stat = function() {
        var query = prefix.concat("Pool/Status");

        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.pool.values = function(pool, data_type, key) {
        var query = prefix.concat("Pool");

        if (data_type == "meta") {
            query = query.concat(pool, "/Meta/", key);
        }  else if (data_type == "gene") {
            query = query.concat(pool, "/Gene/", key);
        }

        return $.ajax(query, {dataType: "json"}).then(x => x);
    }

    // Clusters API

    output.clusters = {}

    output.clusters.sigProjMatrix = function(cluster_variable, meta)
    {
        var query = prefix.concat("Clusters/")
        query = query.concat(cluster_variable)
        if (meta) {
            query = query.concat("/SigProjMatrix/Meta")
        } else {
            query = query.concat("/SigProjMatrix/Normal")
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.clusters.cells = function(cluster_variable) {
        var query = prefix.concat("Clusters/")
        query = query.concat(cluster_variable)
        query = query.concat("/Cells")
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.clusters.list = function() {
        var query = prefix.concat("Clusters/list")
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup = {}

    output.filterGroup.pCorr = function(meta) {
        var query = prefix.concat("FilterGroup");
        if (meta) {
            query = query.concat("/PearsonCorr/Meta");
        } else {
            query = query.concat("/PearsonCorr/Normal");
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.listPCs = function()
    {
        var query = prefix.concat("FilterGroup/PearsonCorr/list")
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.loadings_pos = function(pcnum) {
        var query = prefix.concat("FilterGroup/");
        query = query.concat(encodeURI(pcnum), "/Loadings");
        query = query.concat("/Positive");
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.loadings_neg= function(pcnum) {
        var query = prefix.concat("FilterGroup/");
        query = query.concat(encodeURI(pcnum), "/Loadings");
        query = query.concat("/Negative");
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Projections API

    output.projections = {}

    output.projections.coordinates = function(projection_name)
    {
        var query = prefix.concat("Projections/")
        query = query.concat(encodeURI(projection_name), "/coordinates")
        return $.ajax(query, {dataType: "json"}).then(x => fix_coordinates(x))
    }

    output.projections.list = function()
    {
        var query = prefix.concat("Projections/list")
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Tree API

    output.tree = {}

    // Get milestones and connectivity for this 2d projection
    output.tree.milestones = function(projection) {
        var query = prefix.concat("Tree/Projections/");
        query = query.concat(encodeURI(projection), "/milestones");
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Get cell coordinates for this 2d projection
    output.tree.coordinates = function(projection) {
        var query = prefix.concat("Tree/Projections/");
        query = query.concat(encodeURI(projection), "/coordinates")
        return $.ajax(query, {dataType: "json"}).then(x => fix_coordinates(x))
    }

    output.tree.sigProjMatrix = function(meta)
    {
        var query = prefix.concat("Tree")
        if (meta) {
            query = query.concat("/SigProjMatrix/Meta")
        } else {
            query = query.concat("/SigProjMatrix/Normal")
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }



    // PC API

    output.pc = {}

    output.pc.coordinates = function() {
        var query = prefix.concat("FilterGroup/PCA/Coordinates")
        return $.ajax(query, {dataType: "json"}).then(x => fix_coordinates(x))
    }	

    output.pc.versus = function(pc1, pc2) {
        var query = prefix.concat("FilterGroup");
        query = query.concat("/PCVersus/", encodeURI(pc1), "/", encodeURI(pc2));
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Expression API


    output.expression = {}

    output.expression.gene = function(gene_name) {
        var query = prefix.concat("Expression/Gene/");
        query = query.concat(encodeURI(gene_name));
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.expression.genes = {}
    output.expression.genes.list = function() {
        var query = prefix.concat("Expression/Genes/List");
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Analysis API

    output.analysis = {}

    output.analysis.run = function(subset) {
        var query = prefix.concat("Analysis/Run/");
        return $.ajax({
            type: "POST",
            url: query,
            data: JSON.stringify(subset),
        }).done(alert("Running Subset Analysis"));	
    }

    // Session Info Api

    output.sessionInfo = function() {
        var query = prefix.concat("SessionInfo")
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    return output;

})();
