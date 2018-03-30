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
        var query = "/Signature/Info/"
        query = query.concat(encodeURI(sig_name));
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.listMeta = function(){
        var query = "/Signature/ListMeta"
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.scores = function(sig_name){
        var query = "/Signature/Scores/"
        query = query.concat(encodeURI(sig_name))
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.meta = function(meta_name){
        var query = "/Signature/Meta/"
        query = query.concat(encodeURI(meta_name))
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.expression = function(sig_name){
        var query = "/Signature/Expression/"
        query = query.concat(encodeURI(sig_name))
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.clusters = function(meta) {
        var query = ""
        if (meta) {
            query = query.concat("/FilterGroup/SigClusters/Meta");
        } else {
            query = query.concat("/FilterGroup/SigClusters/Normal");
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Clusters API

    output.clusters = {}

    output.clusters.sigProjMatrix = function(meta, pvalue)
    {
        var query = "/Clusters"
        if (meta) {
            query = query.concat("/SigProjMatrix/Meta")
        } else {
            query = query.concat("/SigProjMatrix/Normal")
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.clusters.cells = function() {
        var query = "/Clusters"
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup = {}

    output.filterGroup.pCorr = function(meta) {
        var query = "/FilterGroup";
        if (meta) {
            query = query.concat("/PearsonCorr/Meta");
        } else {
            query = query.concat("/PearsonCorr/Normal");
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.listPCs = function()
    {
        var query = "/FilterGroup/PearsonCorr/list"
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.loadings_pos = function(pcnum) {
        var query = "/FilterGroup/";
        query = query.concat(encodeURI(pcnum), "/Loadings");
        query = query.concat("/Positive");
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.loadings_neg= function(pcnum) {
        var query = "/FilterGroup/";
        query = query.concat(encodeURI(pcnum), "/Loadings");
        query = query.concat("/Negative");
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Projections API

    output.projections = {}

    output.projections.coordinates = function(projection_name)
    {
        var query = "/Projections/"
        query = query.concat(encodeURI(projection_name), "/coordinates")
        return $.ajax(query, {dataType: "json"}).then(x => fix_coordinates(x))
    }

    output.projections.list = function()
    {
        var query = "/Projections/list"
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.projections.sigProjMatrix = function(meta)
    {
        var query = "/Projections"
        if (meta) {
            query = query.concat("/SigProjMatrix/Meta")
        } else {
            query = query.concat("/SigProjMatrix/Normal")
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }


    // Tree API

    output.tree = {}

    output.tree.tree = function()
    {
        var query = "/Tree/List"
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.tree.tree_points = function(projection) {
        var query = "/Tree/";
        query = query.concat(encodeURI(projection), "/Points");
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.tree.coordinates = function(projection) {
        var query = "/Tree/";
        query = query.concat(encodeURI(projection), "/Projection")
        return $.ajax(query, {dataType: "json"}).then(x => fix_coordinates(x))
    }

    output.tree.sigProjMatrix = function(meta)
    {
        var query = "/Tree"
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
        var query = "/FilterGroup/PCA/Coordinates"
        return $.ajax(query, {dataType: "json"}).then(x => fix_coordinates(x))
    }	

    output.pc.versus = function(pc1, pc2) {
        var query = "/FilterGroup";
        query = query.concat("/PCVersus/", encodeURI(pc1), "/", encodeURI(pc2));
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

    // Session Info Api

    output.sessionInfo = function() {
        var query = "/SessionInfo"
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    return output;

})();
