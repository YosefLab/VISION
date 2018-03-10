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

    // FilterGroup API

    output.filterGroup = {}

    output.filterGroup.listProjections = function()
    {
        var query = "/FilterGroup/projections/list"
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.sigProjMatrix = function(meta)
    {
        var query = "/FilterGroup"
        if (meta) {
            query = query.concat("/SigProjMatrix/Meta")
        } else {
            query = query.concat("/SigProjMatrix/Normal")
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.sigProjMatrixP = function(meta, pvalue)
    {
        var query = "/FilterGroup"
        if (meta) {
            query = query.concat("/SigProjMatrix_P/Meta")
        } else {
            if (pvalue == "nominal") {
                query = query.concat("/SigProjMatrix_P/Normal")
            } else {
                query = query.concat("/SigProjMatrix_Pemp/Normal")
            }
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.sigProjMatrixPClusters = function(meta, pvalue)
    {
        var query = "/FilterGroup"
        if (meta) {
            query = query.concat("/SigProjMatrix_P_Clusters/Meta")
        } else {
            if (pvalue == "nominal") {
                query = query.concat("/SigProjMatrix_P_Clusters/Normal")
            } else {
                query = query.concat("/SigProjMatrix_Pemp_Clusters/Normal")
            }
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.treeSigProjMatrixP = function(meta)
    {
        var query = "/FilterGroup"
        if (meta) {
            query = query.concat("/Tree/SigProjMatrix_P/Meta")
        } else {
            query = query.concat("/Tree/SigProjMatrix_P/Normal")
        }
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.filterGroup.list = function()
    {
        var query = "/FilterGroup/list";
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

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

    // Projection API

    output.projection = {}
     
    output.projection.coordinates = function(projection_name)
    {
        var query = "/FilterGroup/"
        query = query.concat(encodeURI(projection_name), "/coordinates")
        return $.ajax(query, {dataType: "json"}).then(x => fix_coordinates(x))
    }

    output.projection.clusters = function(projection_name,
        cluster_method, parameter)
    {
        var query = "/FilterGroup/"
        query = query.concat(encodeURI(projection_name), "/clusters/",
            encodeURI(cluster_method), "/", encodeURI(parameter))
        return $.ajax(query, {dataType: "json"}).then(x => x)

    }

    // Tree API
    
    output.tree = {}

    output.tree.tree = function() 
    {
        var query = "/FilterGroup/Tree/List"
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.tree.tree_points = function(projection) {
        var query = "/FilterGroup/";
        query = query.concat(encodeURI(projection), "/Tree/Points");
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.tree.coordinates = function(projection) {
        var query = "/FilterGroup/";
        query = query.concat(encodeURI(projection), "/Tree/Projection")
        return $.ajax(query, {dataType: "json"}).then(x => fix_coordinates(x))
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

    // Misc

    output.cellClusters = function() {
        var query = "/Clusters"
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    return output;

})();
