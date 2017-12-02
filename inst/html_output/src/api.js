var api = (function(){

    /* Coordinates come in in the form:
     * [
     *  [x, y, sample_name],
     *  [x, y, sample_name2],
     *  ...
     * ]
     * Converting back to the format we want is faster in Javascript than in R
     */

    var fix_coordinates = function(x){
        var result =  _(x)
            .keyBy(x => x[2])
            .mapValues(x => x.slice(0, 2))
            .value();

        return result;
    }

    var output = {};

    // Signature API
    
    output.signature = {}

    output.signature.info = function(sig_name){
        var query = "/Signature/Info/"
        query = query.concat(encodeURI(sig_name));
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.listPrecomputed = function(){
        var query = "/Signature/ListPrecomputed"
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.scores = function(sig_name){
        var query = "/Signature/Scores/"
        query = query.concat(encodeURI(sig_name))
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.ranks = function(sig_name){
        var query = "/Signature/Ranks/"
        query = query.concat(encodeURI(sig_name))
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.expression = function(sig_name){
        var query = "/Signature/Expression/"
        query = query.concat(encodeURI(sig_name))
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.clusters = function(precomputed, filter) {
        var query = ""
        if (precomputed) {
            query = query.concat("/FilterGroup/", filter, "/SigClusters/Precomputed");

        } else {
            query = query.concat("/FilterGroup/", filter, "/SigClusters/Normal");
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // FilterGroup API

    output.filterGroup = {}

    output.filterGroup.listProjections = function(filter_group)
    {
        var query = "/FilterGroup/"
        query = query.concat(encodeURI(filter_group), "/projections/list")
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.sigProjMatrix = function(filter_group, precomputed)
    {
        var query = "/FilterGroup/"
        if (precomputed) {
            query = query.concat(encodeURI(filter_group), "/SigProjMatrix/Precomputed")
        } else {
            query = query.concat(encodeURI(filter_group), "/SigProjMatrix/Normal")
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.sigProjMatrixP = function(filter_group, precomputed)
    {
        var query = "/FilterGroup/"
        if (precomputed) {
            query = query.concat(encodeURI(filter_group), "/SigProjMatrix_P/Precomputed")
        } else {
            query = query.concat(encodeURI(filter_group), "/SigProjMatrix_P/Normal")
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.treeSigProjMatrixP = function(filter_group, precomputed)
    {
        var query = "/FilterGroup/"
        if (precomputed) {
            query = query.concat(encodeURI(filter_group), "/Tree/SigProjMatrix_P/Precomputed")
        } else {
            query = query.concat(encodeURI(filter_group), "/Tree/SigProjMatrix_P/Normal")
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.list = function()
    {
        var query = "/FilterGroup/list";
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.genes = function(filter_group)
    {
        var query = "/FilterGroup/"
        query = query.concat(encodeURI(filter_group), "/genes")
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.pCorr = function(filter_group, precomputed) {
        var query = "/FilterGroup/";
        if (precomputed) {
            query = query.concat(encodeURI(filter_group), "/PearsonCorr/Precomputed");
        } else {
            query = query.concat(encodeURI(filter_group), "/PearsonCorr/Normal");
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.loadings_pos = function(filter_group, pcnum) {
        var query = "/FilterGroup/";
        query = query.concat(encodeURI(filter_group), "/", encodeURI(pcnum), "/Loadings");
        query = query.concat("/Positive");
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.loadings_neg= function(filter_group, pcnum) {
        var query = "/FilterGroup/";
        query = query.concat(encodeURI(filter_group), "/", encodeURI(pcnum), "/Loadings");
        query = query.concat("/Negative");
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Projection API

    output.projection = {}
     
    output.projection.coordinates = function(filter_group, projection_name)
    {
        var query = "/FilterGroup/"
        query = query.concat(encodeURI(filter_group), "/",
            encodeURI(projection_name), "/coordinates")
        return $.ajax(query, {dataType: "json"}).then(x => fix_coordinates(x))
    }

    output.projection.clusters = function(filter_group, projection_name,
        cluster_method, parameter)
    {
        var query = "/FilterGroup/"
        query = query.concat(encodeURI(filter_group), "/", encodeURI(projection_name),
            "/clusters/", encodeURI(cluster_method), "/", encodeURI(parameter))
        return $.ajax(query, {dataType: "json"}).then(x => x)

    }

    // Tree API
    
    output.tree = {}

    output.tree.tree = function(filter_group) 
    {
        var query = "/FilterGroup/";
        query = query.concat(encodeURI(filter_group), "/Tree/List");
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.tree.tree_points = function(filter_group, projection) {
        var query = "/FilterGroup/";
        query = query.concat(encodeURI(filter_group), "/", encodeURI(projection), "/Tree/Points");
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.tree.coordinates = function(filter_group, projection) {
        var query = "/FilterGroup/";
        query = query.concat(encodeURI(filter_group), "/", encodeURI(projection), "/Tree/Projection")
        return $.ajax(query, {dataType: "json"}).then(x => fix_coordinates(x))
    }


    // PC API

    output.pc = {}

    output.pc.coordinates = function(filter_group, sig_name, pcnum) {
        var query = "/FilterGroup/";
        query = query.concat(encodeURI(filter_group), "/", encodeURI(sig_name), "/", encodeURI(pcnum), "/Coordinates");
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }	

    output.pc.versus = function(filter_group, pc1, pc2) {
        var query = "/FilterGroup/";
        query = query.concat(encodeURI(filter_group), "/PCVersus/", encodeURI(pc1), "/", encodeURI(pc2));
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Expression API


    output.expression = {}

    output.expression.gene = function(gene_name) {
        var query = "/Expression/Gene/";
        query = query.concat(encodeURI(gene_name));
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.expression.genes = {}
    output.expression.genes.list = function() {
        var query = "/Expression/Genes/List";
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Analysis API

    output.analysis = {}

    output.analysis.run = function(subset) {
        var query = "/Analysis/Run/";
        return $.ajax({
            type: "POST",
            url: query,
            data: JSON.stringify(subset),
        }).done(alert("Running Subset Analysis"));	
    }

    return output;

})();
