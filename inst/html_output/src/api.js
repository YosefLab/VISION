var api = (function(){


    // location.pathname is either:
    // this is either "/" or "/Results.html" if running normally
    // if running from RStudio Server this will be:
    // "/p/<port>/Results.html" or "/p/<port>/"
    // Need to set this dynamically so that it works with RStudio server

    var prefix = location.pathname.replace("Results.html", "")

    // This is needed to clear the cache
    // Each time the page is loaded, a new random version is chosen so that
    // the first call to every api endpoint is not cached
    var version = Math.floor(Math.random() * parseInt("ffffff", 16)).toString(16)
    var postProcess = function(x){
        return prefix + x + "?v=" + version
    }

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

    var fix_coordinates_1d = function(x){

        var result =  _(x)
            .keyBy(x => x[1])
            .mapValues(x => x[0])
            .value();

        return result;
    }

    var output = {};

    // Signature API

    output.signature = {}

    output.signature.info = function(sig_name){
        var query = "Signature/Info/"
        query = query.concat(encodeURIComponent(sig_name));
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.scores = function(sig_name){
        var query = "Signature/Scores/"
        query = query.concat(encodeURIComponent(sig_name))
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => {
            return _.fromPairs(_.zip(x['cells'], x['values']))
        })
    }

    output.signature.meta = function(meta_name){
        var query = "Signature/Meta/"
        query = query.concat(encodeURIComponent(meta_name))
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => {
            return _.fromPairs(_.zip(x['cells'], x['values']))
        })
    }

    output.signature.expression = function(sig_name, cluster_var){
        var query = "Signature/Expression/"
        query = query.concat(encodeURIComponent(sig_name))
        query = query.concat("/" + encodeURIComponent(cluster_var))
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.signature.clusters = function(meta) {
        var query
        if (meta) {
            query = "FilterGroup/SigClusters/Meta"
        } else {
            query = "FilterGroup/SigClusters/Normal"
        }
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Protein API

    output.protein = {}

    output.protein.clusters = function() {
        var query= "FilterGroup/SigClusters/Proteins"
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.protein.values = function(protein_name){
        var query = "Proteins/"
        query = query.concat(encodeURIComponent(protein_name))
        query = query.concat("/Values")
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => {
            return _.fromPairs(_.zip(x['cells'], x['values']))
        })
    }

    // Clusters API

    output.clusters = {}

    output.clusters.sigProjMatrix = function(cluster_variable, meta)
    {
        var query = "Clusters/"
        query = query.concat(encodeURIComponent(cluster_variable))
        if (meta) {
            query = query.concat("/SigProjMatrix/Meta")
        } else {
            query = query.concat("/SigProjMatrix/Normal")
        }
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.clusters.proteinMatrix = function(cluster_variable)
    {
        var query = "Clusters/"
        query = query.concat(encodeURIComponent(cluster_variable))
        query = query.concat("/ProteinMatrix")
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.clusters.cells = function(cluster_variable) {
        var query = "Clusters/"
        query = query.concat(encodeURIComponent(cluster_variable))
        query = query.concat("/Cells")
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => {
            return _.fromPairs(_.zip(x['cells'], x['values']))
        })
    }

    output.clusters.meta_levels = function() {
        var query = "Clusters/MetaLevels"
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.clusters.list = function() {
        var query = "Clusters/list"
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    //Yanay
    output.de = function(type_n, subtype_n, group_num, type_d, subtype_d, group_denom, min_cells, subsample_groups, subsample_N) {
        var query = "DE"

        query = postProcess(query)
        return $.ajax(query, {
            type: "POST",
            data: JSON.stringify({"type_n":type_n, "type_d":type_d, "subtype_n":subtype_n, "subtype_d":subtype_d,  "group_num":group_num, "group_denom":group_denom, "min_cells":min_cells, "subsample_groups":subsample_groups, "subsample_N":subsample_N}),
            dataType: "json"
        }).then(x => x);
    }

    // Pearson Correlations (LCA) api

    output.pCorr = {}

    output.pCorr.signatures = function() {
        var query = "PearsonCorr/Normal"
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.pCorr.proteins = function() {
        var query = "PearsonCorr/Proteins"
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.pCorr.meta = function() {
        var query = "PearsonCorr/Meta"
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Projections API

    output.projections = {}

    output.projections.coordinates = function(projection_name, projection_column)
    {
        var query = "Projections/"
        query = query.concat(encodeURIComponent(projection_name), "/coordinates/", encodeURIComponent(projection_column))
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => fix_coordinates_1d(x))
    }

    output.projections.list = function()
    {
        var query = "Projections/list"
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Tree API

    output.tree = {}

    // Get cell coordinates for this 2d projection
    output.tree.coordinates = function(projection) {
        var query = "Tree/Projections/"
        query = query.concat(encodeURIComponent(projection), "/coordinates")
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => [fix_coordinates(x[0]), [x[1], x[2]]])
    }


    output.tree.list = function()
    {
        var query = "Tree/Projections/list"
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.tree.sigProjMatrix = function(meta)
    {
        var query = "Tree"
        if (meta) {
            query = query.concat("/SigProjMatrix/Meta")
        } else {
            query = query.concat("/SigProjMatrix/Normal")
        }
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    output.tree.proteinMatrix = function()
    {
        var query = "Tree/ProteinMatrix"
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }


    // Expression API


    output.expression = {}

    output.expression.gene = function(gene_name) {
        var query = "Expression/Gene/"
        query = query.concat(encodeURIComponent(gene_name));
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => {
            return _.fromPairs(_.zip(x['cells'], x['values']))
        })
    }

    output.expression.genes = {}
    output.expression.genes.list = function() {
        var query = "Expression/Genes/List"
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Analysis API

    output.analysis = {}

    output.analysis.run = function(subset) {
        var query = "Analysis/Run/"
        query = postProcess(query)
        return $.ajax({
            type: "POST",
            url: query,
            data: JSON.stringify(subset),
        }).done(alert("Running Subset Analysis"));
    }

    // Cell API

    output.cell = {}

    output.cell.meta = function(cellId){
        var query = "Cell/" + encodeURIComponent(cellId) + "/Meta"

        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    // Cells API

    output.cells = {}

    output.cells.meta = function(subset) {
        var query = "Cells/Meta"

        query = postProcess(query)
        return $.ajax(query, {
            type: "POST",
            data: JSON.stringify(subset),
            dataType: "json"
        }).then(x => x);

    }

    output.cells.saveSelection = function(selectionName, subset) {
        var query = "Cells/Selections/"
        query = query.concat(encodeURIComponent(selectionName))

        query = postProcess(query)
        return $.ajax(query, {
            type: "POST",
            data: JSON.stringify(subset),
            dataType: "json"
        }).then(x => x);

    }

    output.cells.getSelection = function(selectionName) {
        var query = "Cells/Selections/"
        query = query.concat(encodeURIComponent(selectionName))
        query = postProcess(query)
        return $.ajax(query, {dataType: "json" }).then(x => x);

    }

    output.cells.listSelections = function() {
        var query = "Cells/Selections"
        query = postProcess(query)
        return $.ajax(query, {dataType: "json" }).then(x => x);

    }
    
    // Modules Api
    // the modules api uses existing signature api endpoints with only some specific extra ones
    output.modules = {}
    
    output.modules.lcs = function() {
        var query = "Modules/LC"
        query = postProcess(query)
        return $.ajax(query, {dataType: "json" }).then(x => x);
    }

    // Session Info Api

    output.sessionInfo = function() {
        var query = "SessionInfo"
        query = postProcess(query)
        return $.ajax(query, {dataType: "json"}).then(x => x)
    }

    return output;

})();
